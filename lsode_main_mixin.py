from __future__ import annotations
from typing import TYPE_CHECKING, Iterable
import sys

# Satisfy Mypy Statically
if TYPE_CHECKING:
    # Use a simple import for mypy. This code block is IGNORED at runtime.
    from formatter_lsode import LsodeFormatter

class LsodeMainMixin:
    def make_main(self,
                  y_initial: Iterable[str],
                  t_initial: str,
                  tout_initial: str,
                  rtol: str | Iterable[str],
                  atol: str | Iterable[str],
                  num_steps: int = 12,
                  tout_multiplier: str = '10.D0',
                  mf: int = 21,
                  itask: int = 1,
                  iopt: int = 0,
                  ) -> str:
        """
        Generates a modern Fortran main program to run a DLSODE simulation,
        based on the dlsode_example.txt driver.
        """

        LsodeFormatter = sys.modules['formatter_lsode'].LsodeFormatter
        if not isinstance(self, LsodeFormatter):
            raise TypeError('This is a mixin for LsodeFormatter only.')
        del LsodeFormatter

        lines = []
        n = self.neq  # This is NEQ

        # --- Validate Inputs based on ITOL ---
        y_init_list = list(y_initial)
        if len(y_init_list) != n:
            raise ValueError(f"y_initial must have length {n} (NEQ), but has {len(y_init_list)}")

        # Determine if RTOL and ATOL are scalar (str) or vector (Iterable)
        is_rtol_scalar, rtol_list = (True, [rtol]) if isinstance(rtol, str) else (False, list(rtol))
        is_atol_scalar, atol_list = (True, [atol]) if isinstance(atol, str) else (False, list(atol))

        # --- 1. Derive ITOL based on argument types ---

        if is_rtol_scalar and is_atol_scalar:
            # Case 1: Scalar RTOL, Scalar ATOL
            itol = 1
        elif is_rtol_scalar and not is_atol_scalar:
            # Case 2: Scalar RTOL, Array ATOL
            itol = 2
        elif not is_rtol_scalar and is_atol_scalar:
            # Case 3: Array RTOL, Scalar ATOL
            itol = 3
        elif not is_rtol_scalar and not is_atol_scalar:
            # Case 4: Array RTOL, Array ATOL
            itol = 4

        # --- 2. Validate vector lengths (must equal NEQ) ---

        # Validate ATOL length if it is an array (ITOL=2 or 4)
        if not is_atol_scalar and len(atol_list) != n:
            raise ValueError(f"ATOL must be an array of length NEQ ({n}) for ITOL={itol}.")

        # Validate RTOL length if it is an array (ITOL=3 or 4)
        if not is_rtol_scalar and len(rtol_list) != n:
            raise ValueError(f"RTOL must be an array of length NEQ ({n}) for ITOL={itol}.")

        # --- Calculate LRW and LIW ---
        # Based on MF = 21 (BDF, user-supplied dense Jacobian) [cite: 3]
        # LRW = 22 + NEQ*NEQ + 9*NEQ
        # LIW = 20 + NEQ
        if mf in (21, 22): # 21=user jac, 22=internal jac
            lrw = 22 + n*n + 9*n
            liw = 20 + n
        # Add other MF calculations as needed
        else:
            raise NotImplementedError(f"LRW/LIW calculation for MF={mf} is not implemented.")

        # --- Start writing the program ---
        lines.append('program main_driver')
        lines.append('  use, intrinsic :: iso_c_binding, only: dp => c_double')
        lines.append('  use functions')
        lines.append('  implicit none')
        lines.append('')

        # --- Declarations ---
        lines.append('  ! Solver parameters')
        lines.append('  integer :: iopt, istate, itask, itol, liw, lrw, mf, neq, iout')
        lines.append(f'  integer :: iwork({liw})')
        lines.append('')
        lines.append('  ! State and time variables')
        lines.append('  double precision :: t, tout')
        lines.append(f'  double precision :: y({self.n_yca})')

        # RTOL declaration
        if is_rtol_scalar:
            lines.append('  double precision :: rtol')
        else:
            lines.append(f'  double precision :: rtol({n})')

        # ATOL declaration
        if is_atol_scalar:
            lines.append('  double precision :: atol')
        else:
            lines.append(f'  double precision :: atol({n})')

        lines.append(f'  double precision :: rwork({lrw})')
        lines.append('')

        # --- Initialization ---
        lines.append('  ! --- Initialization ---')
        lines.append(f'  neq = {n}')
        lines.append('')

        # Y initial conditions [cite: 1, 2]
        for i, val in enumerate(y_init_list):
            lines.append(f'  y({1+i}) = {val}')

        lines.append(f'  t = {t_initial}')
        lines.append(f'  tout = {tout_initial}')
        lines.append('')

        # --- Solver Configuration ---
        lines.append('  ! --- Solver Configuration ---')
        lines.append(f'  itol = {itol}')

        # RTOL assignment [cite: 2]
        if is_rtol_scalar:
            lines.append(f'  rtol = {rtol_list[0]}')
        else:
            for i, val in enumerate(rtol_list):
                lines.append(f'  rtol({1+i}) = {val}')

        # ATOL assignment [cite: 2]
        if is_atol_scalar:
            lines.append(f'  atol = {atol_list[0]}')
        else:
            for i, val in enumerate(atol_list):
                lines.append(f'  atol({1+i}) = {val}')

        lines.append(f'  itask = {itask}       ! 1 = Normal computation; 2 = one step only and return')
        lines.append('  istate = 1       ! 1 = First call [cite: 2]')
        lines.append(f'  iopt = {iopt}        ! 0 = No optional inputs [cite: 3]')
        lines.append(f'  lrw = {lrw}      ! [cite: 3]')
        lines.append(f'  liw = {liw}      ! [cite: 3]')
        lines.append(f'  mf = {mf}        ! 21 = BDF, user-supplied dense Jacobian [cite: 3]')
        lines.append('')

        # Dynamically create the format string
        y_format = f"{n}ES14.6"

        # Dynamically create the WRITE statement using an inline format string
        y_list = ", ".join([f"y({i+1})" for i in range(n)])
        inline_write_format = f"'(\" At t =\", ES12.4, \"   y =\", {y_format})'"

        lines.append('  ! --- Integration Loop ---')
        lines.append(f'  do iout = 1, {num_steps}')
        lines.append('    call dlsode(f, neq, y, t, tout, itol, rtol, atol, itask, &')
        lines.append(' &              istate, iopt, rwork, lrw, iwork, liw, jac, mf)')
        lines.append('')

        # --- MODERN WRITE: Replaces 'write(*, 20)' ---
        lines.append(f'    write(*, {inline_write_format}) t, {y_list}')
        lines.append('')

        # Error handling
        lines.append('    if (istate < 0) then')
        # --- MODERN WRITE: Replaces 'write(*, 90)' ---
        lines.append("      write(*, '(/\" Error halt.. ISTATE =\", I3)') istate")
        lines.append('      stop "Error halt"')
        lines.append('    end if')
        lines.append('')
        lines.append(f'    tout = tout * {tout_multiplier} !')
        lines.append('  end do')
        lines.append('')

        # --- Final Statistics ---
        lines.append('  ! --- Final Statistics ---')
        # --- MODERN WRITE: Replaces 'write(*, 60)' ---
        lines.append("  write(*, '(/\" No. steps =\", i4, \",  No. f-s =\", i4, \",  No. J-s =\", i4)') iwork(11), iwork(12), iwork(13)")
        lines.append('')

        lines.append('end program main_driver')

        return '\n'.join(lines)
