program lsode_program
  use, intrinsic :: iso_c_binding, only: dp => c_double
  use functions
  implicit none

  ! Solver parameters
  integer, parameter :: neq = 3
  integer, parameter :: itol = 2
  integer, parameter :: itask = 1  ! 1 = Normal computation; 2 = one step only and return
  integer, parameter :: iopt = 0   ! 0 = No optional inputs
  integer, parameter :: mf = 21    ! 21 = BDF, user-supplied dense Jacobian
  integer :: istate, iout

  ! State and time variables
  real(dp) :: t, tout
  real(dp) :: y(15)  ! neq + consts + functions (n_yca)

  ! Relative and absolute tolerances
  real(dp) :: rtol
  real(dp) :: atol(neq)

  ! Work arrays (real and integer)
  integer, parameter :: lrw = 58
  integer, parameter :: liw = 23
  integer :: iwork(liw)
  real(dp) :: rwork(lrw)

  rtol = 0.0001_dp
  atol(1) = 1e-06_dp
  atol(2) = 1e-10_dp
  atol(3) = 1e-06_dp

  istate = 1  ! First call

  ! --- Initialization ---

  t = 0.0_dp
  tout = 0.4_dp

  ! y0
  y(1) = 1.0_dp
  y(2) = 0.0_dp
  y(3) = 0.0_dp

  ! Constants
  y(4) = 0.04_dp                        ! c1
  y(5) = 10000.0_dp                     ! c2
  y(6) = 30000000.0_dp                  ! c3

  ! --- Integration Loop ---
  do iout = 1, 12
    call dlsode(f, neq, y, t, tout, itol, rtol, atol, itask, &
 &              istate, iopt, rwork, lrw, iwork, liw, jac, mf)

    write(*, '(" At t =", ES12.4, "   y =", 3ES14.6)') t, y(1), y(2), y(3)

    if (istate < 0) then
      write(*, '(/" Error halt.. ISTATE =", I3)') istate
      stop "Error halt"
    end if

    tout = tout + 0.4_dp
  end do

  ! --- Final Statistics ---

  write(*, '(/" No. steps =", i4, ",  No. f-s =", i4, ",  No. J-s =", i4)') iwork(11), iwork(12), iwork(13)

end