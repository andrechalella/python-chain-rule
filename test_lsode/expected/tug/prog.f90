program lsode_program
  use, intrinsic :: iso_c_binding, only: dp => c_double
  use functions
  implicit none

  ! Solver parameters
  integer, parameter :: neq = 6
  integer, parameter :: itol = 2
  integer, parameter :: itask = 1  ! 1 = Normal computation; 2 = one step only and return
  integer, parameter :: iopt = 0   ! 0 = No optional inputs
  integer, parameter :: mf = 21    ! 21 = BDF, user-supplied dense Jacobian
  integer :: istate, iout

  ! State and time variables
  real(dp) :: t, tout
  real(dp) :: y(90)  ! neq + consts + functions (n_yca)

  ! Relative and absolute tolerances
  real(dp) :: rtol
  real(dp) :: atol(neq)

  ! Work arrays (real and integer)
  integer, parameter :: lrw = 112
  integer, parameter :: liw = 26
  integer :: iwork(liw)
  real(dp) :: rwork(lrw)

  rtol = 0.0001_dp
  atol(1) = 1.0_dp
  atol(2) = 1.0_dp
  atol(3) = 0.01_dp
  atol(4) = 1.0_dp
  atol(5) = 1.0_dp
  atol(6) = 0.01_dp

  istate = 1  ! First call

  ! --- Initialization ---

  t = 0.0_dp
  tout = 10.0_dp

  ! y0
  y(1) = 100.0_dp
  y(2) = 0.0_dp
  y(3) = 1.5707963267948966_dp
  y(4) = 0.0_dp
  y(5) = 0.0_dp
  y(6) = 0.0_dp

  ! Constants
  y(7) = 0.2_dp                         ! FPDaft
  y(8) = 0.1_dp                         ! FPDbow
  y(9) = 98066.5_dp                     ! FPDmax
  y(10) = 29419.949999999997_dp          ! FPLmax
  y(11) = 4500000.0_dp                   ! Iz
  y(12) = 100.0_dp                       ! LTL0
  y(13) = 3.141592653589793_dp           ! PI
  y(14) = 212322.33946969695_dp          ! kTL
  y(15) = 300000.0_dp                    ! m
  y(16) = 15.0_dp                        ! xCPmax
  y(17) = 9.0_dp                         ! xCPmed
  y(18) = 15.0_dp                        ! xG
  y(19) = 9.0_dp                         ! xT

  ! --- Integration Loop ---
  do iout = 1, 100
    call dlsode(f, neq, y, t, tout, itol, rtol, atol, itask, &
 &              istate, iopt, rwork, lrw, iwork, liw, jac, mf)

    write(*, '(" At t =", ES12.4, "   y =", 6ES14.6)') t, y(1), y(2), y(3), y(4), y(5), y(6)

    if (istate < 0) then
      write(*, '(/" Error halt.. ISTATE =", I3)') istate
      stop "Error halt"
    end if

    tout = tout + 10.0_dp
  end do

  ! --- Final Statistics ---

  write(*, '(/" No. steps =", i4, ",  No. f-s =", i4, ",  No. J-s =", i4)') iwork(11), iwork(12), iwork(13)

end