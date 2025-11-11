module functions
  use, intrinsic :: iso_c_binding, only: dp => c_double
  implicit none
contains

subroutine f(neq, t, y, ydot)
  integer, intent(in) :: neq
  real(dp), intent(in) :: t
  real(dp), intent(inout) :: y(*)
  real(dp), intent(out) :: ydot(*)
  call update_a(y)
  ydot(1) = y(7)
  ydot(2) = y(11)
  ydot(3) = y(14)
end subroutine

subroutine jac(neq, t, y, ml, mu, pd, nrpd)
  integer, intent(in) :: neq, ml, mu, nrpd
  real(dp), intent(in) :: t, y(*)
  real(dp), intent(out) :: pd(nrpd, *)
  pd(1,1) = y(8)
  pd(1,2) = y(9)
  pd(1,3) = y(10)
  pd(2,1) = y(4)
  pd(2,2) = y(12)
  pd(2,3) = y(13)
  pd(3,2) = y(15)
end subroutine

subroutine update_a(y)
  real(dp), intent(inout) :: y(*)
  y(7) = ef_n1__f1(y)
  y(8) = ef_n1_d0_f1(y)
  y(9) = ef_n1_d1_f1(y)
  y(10) = ef_n1_d2_f1(y)
  y(11) = ef_n1__f2(y)
  y(12) = ef_n1_d1_f2(y)
  y(13) = ef_n1_d2_f2(y)
  y(14) = ef_n1__f3(y)
  y(15) = ef_n1_d1_f3(y)
end subroutine

real(dp) function ef_n1__f1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-y(4)*y(1) + y(5)*y(2)*y(3))
end function

real(dp) function ef_n1__f2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-(y(5)*y(2)*y(3) + y(6)*(y(2))**2) + y(4)*y(1))
end function

real(dp) function ef_n1__f3(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(6)*(y(2))**2
end function

real(dp) function ef_n1_d0_f1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -y(4)
end function

real(dp) function ef_n1_d1_f1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(5)*y(3)
end function

real(dp) function ef_n1_d2_f1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(5)*y(2)
end function

real(dp) function ef_n1_d1_f2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-y(5)*y(3) - 2.0_dp*y(6)*y(2))
end function

real(dp) function ef_n1_d2_f2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -y(5)*y(2)
end function

real(dp) function ef_n1_d1_f3(y) result(res)
  real(dp), intent(in) :: y(*)
  res = 2.0_dp*y(6)*y(2)
end function

end