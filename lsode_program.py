from typing import Iterable

def make_program_file(*,
                      neq: int,
                      n_yca: int,
                      y0: Iterable[float],
                      t0: float,
                      consts: Iterable[tuple[int, str, float]],
                      tout: float,
                      rtol: float | Iterable[float],
                      atol: float | Iterable[float],
                      num_steps: int,
                      tout_multiplier: float,
                      mf: int,
                      itask: int,
                      iopt: int,
                      ) -> str:
    """
    Generates a modern Fortran main program to run a DLSODE simulation,
    based on the dlsode.f example driver.
    """

    # --- Validate Inputs based on ITOL ---
    y_init_list = list(y0)
    if len(y_init_list) != neq:
        raise ValueError(f"y0 must have length {neq} (NEQ), but has {len(y_init_list)}")

    # Determine if RTOL and ATOL are scalar (str) or vector (Iterable)
    is_rtol_scalar, rtol_list = (True, [rtol]) if isinstance(rtol, float) else (False, list(rtol))
    is_atol_scalar, atol_list = (True, [atol]) if isinstance(atol, float) else (False, list(atol))

    # --- 1. Derive ITOL based on argument types ---
    if is_rtol_scalar and is_atol_scalar:
        itol = 1
    elif is_rtol_scalar and not is_atol_scalar:
        itol = 2
    elif not is_rtol_scalar and is_atol_scalar:
        itol = 3
    elif not is_rtol_scalar and not is_atol_scalar:
        itol = 4

    rtol_decl, rtol_init = _process_xtol(neq, itol, "RTOL", is_rtol_scalar, rtol_list)
    atol_decl, atol_init = _process_xtol(neq, itol, "ATOL", is_atol_scalar, atol_list)

    # --- Calculate LRW and LIW ---
    # Based on MF = 21 (BDF, user-supplied dense Jacobian)
    # LRW = 22 + NEQ*NEQ + 9*NEQ
    # LIW = 20 + NEQ
    if mf in (21, 22): # 21=user jac, 22=internal jac
        lrw = 22 + neq*neq + 9*neq
        liw = 20 + neq
    # Add other MF calculations as needed
    else:
        raise NotImplementedError(f"LRW/LIW calculation for MF={mf} is not implemented.")

    y_init = [f'  y({1+i}) = {_dp(val)}' for i, val in enumerate(y_init_list)]
    consts_init = [f'  y({1+i}) = {_dp(val):30} ! {name}' for i, name, val in consts]

    # Dynamically create the WRITE statement using an inline format string
    y_list = ", ".join([f"y({i+1})" for i in range(neq)])
    inline_write_format = f"'(\" At t =\", ES12.4, \"   y =\", {neq}ES14.6)'"

    tout_expression = f'+ {_dp(tout)}' if tout_multiplier == 1.0 else f'* {_dp(tout_multiplier)}'

    return _TEMPLATE_PROGRAM.format(
            neq=neq,
            n_yca=n_yca,
            itol=itol,
            itask=itask,
            iopt=iopt,
            mf=mf,
            lrw=lrw,
            liw=liw,
            rtol_decl=rtol_decl,
            atol_decl=atol_decl,
            rtol_init='\n'.join(rtol_init),
            atol_init='\n'.join(atol_init),
            t0=_dp(t0),
            tout=_dp(tout),
            tout_expression=tout_expression,
            y_init='\n'.join(y_init),
            consts_init='\n'.join(consts_init),
            num_steps=num_steps,
            y_list=y_list,
            inline_write_format=inline_write_format,
            )

def _dp(f: float) -> str:
    return f'{f}_dp'

def _process_xtol(neq: int,
                  itol: int,
                  xtol_name: str,
                  is_xtol_scalar: bool,
                  xtol_list: list[float],
                  ) -> tuple[str, list[str]]:

    name_up = xtol_name.upper()
    name_lo = xtol_name.lower()

    if not is_xtol_scalar and len(xtol_list) != neq:
        raise ValueError(f"{name_up} must be an array of length {neq} (NEQ) for ITOL={itol}.")

    xtol_decl = f'  real(dp) :: {name_lo}' if is_xtol_scalar \
                    else f'  real(dp) :: {name_lo}(neq)'

    xtol_init = []
    if is_xtol_scalar:
        xtol_init.append(f'  {name_lo} = {_dp(xtol_list[0])}')
    else:
        for i, val in enumerate(xtol_list):
            xtol_init.append(f'  {name_lo}({1+i}) = {_dp(val)}')

    return (xtol_decl, xtol_init)


_TEMPLATE_PROGRAM = """\
program lsode_program
  use, intrinsic :: iso_c_binding, only: dp => c_double
  use functions
  implicit none

  ! Solver parameters
  integer, parameter :: neq = {neq}
  integer, parameter :: itol = {itol}
  integer, parameter :: itask = {itask}  ! 1 = Normal computation; 2 = one step only and return
  integer, parameter :: iopt = {iopt}   ! 0 = No optional inputs
  integer, parameter :: mf = {mf}    ! 21 = BDF, user-supplied dense Jacobian
  integer :: istate, iout

  ! State and time variables
  real(dp) :: t, tout
  real(dp) :: y({n_yca})  ! neq + consts + functions (n_yca)

  ! Relative and absolute tolerances
{rtol_decl}
{atol_decl}

  ! Work arrays (real and integer)
  integer, parameter :: lrw = {lrw}
  integer, parameter :: liw = {liw}
  integer :: iwork(liw)
  real(dp) :: rwork(lrw)

{rtol_init}
{atol_init}

  istate = 1  ! First call

  ! --- Initialization ---

  t = {t0}
  tout = {tout}

  ! y0
{y_init}

  ! Constants
{consts_init}

  ! --- Integration Loop ---
  do iout = 1, {num_steps}
    call dlsode(f, neq, y, t, tout, itol, rtol, atol, itask, &
 &              istate, iopt, rwork, lrw, iwork, liw, jac, mf)

    write(*, {inline_write_format}) t, {y_list}

    if (istate < 0) then
      write(*, '(/" Error halt.. ISTATE =", I3)') istate
      stop "Error halt"
    end if

    tout = tout {tout_expression}
  end do

  ! --- Final Statistics ---

  write(*, '(/" No. steps =", i4, ",  No. f-s =", i4, ",  No. J-s =", i4)') iwork(11), iwork(12), iwork(13)

end
"""
