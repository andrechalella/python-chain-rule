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
  ydot(1) = y(4)
  ydot(2) = y(5)
  ydot(3) = y(6)
  ydot(4) = y(73)
  ydot(5) = y(76)
  ydot(6) = y(86)
end subroutine

subroutine jac(neq, t, y, ml, mu, pd, nrpd)
  integer, intent(in) :: neq, ml, mu, nrpd
  real(dp), intent(in) :: t, y(*)
  real(dp), intent(out) :: pd(nrpd, *)
  pd(1,4) = 1.0_dp
  pd(2,5) = 1.0_dp
  pd(3,6) = 1.0_dp
  pd(4,1) = y(74)
  pd(4,2) = y(75)
  pd(4,3) = y(87)
  pd(5,1) = y(77)
  pd(5,2) = y(78)
  pd(5,3) = y(88)
  pd(6,1) = y(81)
  pd(6,2) = y(82)
  pd(6,3) = y(90)
end subroutine

subroutine update_a(y)
  real(dp), intent(inout) :: y(*)
  y(20) = ef_n1__XTL(y)
  y(21) = ef_n1_d0_XTL(y)
  y(22) = ef_n1_d2_XTL(y)
  y(23) = ef_n1__YTL(y)
  y(24) = ef_n1_d1_YTL(y)
  y(25) = ef_n1_d2_YTL(y)
  y(26) = ef_n1__a2(y)
  y(27) = ef_n1_d2_a2(y)
  y(28) = ef_n1__xCP(y)
  y(29) = ef_n1_d2_xCP(y)
  y(30) = ef_n2__FPD1(y)
  y(31) = ef_n2__FPD2(y)
  y(32) = ef_n2__FPL1(y)
  y(33) = ef_n2__FPL2(y)
  y(34) = ef_n3_d2_FPD1(y)
  y(35) = ef_n3_d2_FPD2(y)
  y(36) = ef_n3_d2_FPL1(y)
  y(37) = ef_n3_d2_FPL2(y)
  y(38) = ef_n3__LTL(y)
  y(39) = ef_n3__tang(y)
  y(40) = ef_n3_d1_tang(y)
  y(41) = ef_n4__FPD(y)
  y(42) = ef_n4__FPL(y)
  y(43) = ef_n4_d0_LTL(y)
  y(44) = ef_n4_d1_LTL(y)
  y(45) = ef_n4__cosg(y)
  y(46) = ef_n4__dLTL(y)
  y(47) = ef_n4__gamma(y)
  y(48) = ef_n4__sing(y)
  y(49) = ef_n4_d0_tang(y)
  y(50) = ef_n5_d2_FPD(y)
  y(51) = ef_n5_d2_FPL(y)
  y(52) = ef_n5__FTL(y)
  y(53) = ef_n5_d2_LTL(y)
  y(54) = ef_n5_d2_tang(y)
  y(55) = ef_n6_d0_cosg(y)
  y(56) = ef_n6_d1_cosg(y)
  y(57) = ef_n6_d0_gamma(y)
  y(58) = ef_n6_d1_gamma(y)
  y(59) = ef_n6_d0_sing(y)
  y(60) = ef_n6_d1_sing(y)
  y(61) = ef_n7_d0_FTL(y)
  y(62) = ef_n7_d1_FTL(y)
  y(63) = ef_n7_d2_cosg(y)
  y(64) = ef_n7_d2_gamma(y)
  y(65) = ef_n7_d2_sing(y)
  y(66) = ef_n8_d2_FTL(y)
  y(67) = ef_n11__RX(y)
  y(68) = ef_n11_d0_RX(y)
  y(69) = ef_n11_d1_RX(y)
  y(70) = ef_n11__RY(y)
  y(71) = ef_n11_d0_RY(y)
  y(72) = ef_n11_d1_RY(y)
  y(73) = ef_n12__f4(y)
  y(74) = ef_n12_d0_f4(y)
  y(75) = ef_n12_d1_f4(y)
  y(76) = ef_n12__f5(y)
  y(77) = ef_n12_d0_f5(y)
  y(78) = ef_n12_d1_f5(y)
  y(79) = ef_n13_d0_M(y)
  y(80) = ef_n13_d1_M(y)
  y(81) = ef_n14_d0_f6(y)
  y(82) = ef_n14_d1_f6(y)
  y(83) = ef_n16__M(y)
  y(84) = ef_n17_d2_RX(y)
  y(85) = ef_n17_d2_RY(y)
  y(86) = ef_n17__f6(y)
  y(87) = ef_n18_d2_f4(y)
  y(88) = ef_n18_d2_f5(y)
  y(89) = ef_n30_d2_M(y)
  y(90) = ef_n31_d2_f6(y)
end subroutine

real(dp) function ef_n12__f4(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(67))/(y(15))
end function

real(dp) function ef_n12__f5(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(70))/(y(15))
end function

real(dp) function ef_n17__f6(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(83))/(y(11))
end function

real(dp) function ef_n4__FPD(y) result(res)
  real(dp), intent(in) :: y(*)
  if (abs(y(26)) < 0.5_dp*y(13)) then
    res = y(30)
  else
    res = y(31)
  end if
end function

real(dp) function ef_n5_d2_FPD(y) result(res)
  real(dp), intent(in) :: y(*)
  if (abs(y(26)) < 0.5_dp*y(13)) then
    res = y(34)
  else
    res = y(35)
  end if
end function

real(dp) function ef_n2__FPD1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(9)*(y(8) + (1.0_dp - y(8))*(sin(y(26)))**2)
end function

real(dp) function ef_n3_d2_FPD1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(9)*(1.0_dp - y(8))*sin(2.0_dp*y(26))*y(27)
end function

real(dp) function ef_n2__FPD2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(9)*(y(7) + (1.0_dp - y(7))*(sin(y(26)))**2)
end function

real(dp) function ef_n3_d2_FPD2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(9)*(1.0_dp - y(7))*sin(2.0_dp*y(26))*y(27)
end function

real(dp) function ef_n4__FPL(y) result(res)
  real(dp), intent(in) :: y(*)
  if (abs(y(3)) < 0.25_dp*y(13)) then
    res = y(32)
  else if (abs(y(3)) < 0.375_dp*y(13)) then
    res = y(33)
  else
    res = 0.0_dp
  end if
end function

real(dp) function ef_n5_d2_FPL(y) result(res)
  real(dp), intent(in) :: y(*)
  if (abs(y(3)) < 0.25_dp*y(13)) then
    res = y(36)
  else if (abs(y(3)) < 0.375_dp*y(13)) then
    res = y(37)
  else
    res = 0.0_dp
  end if
end function

real(dp) function ef_n2__FPL1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(10)*sin(2.0_dp*y(26))
end function

real(dp) function ef_n3_d2_FPL1(y) result(res)
  real(dp), intent(in) :: y(*)
  res = 2.0_dp*y(10)*cos(2.0_dp*y(26))*y(27)
end function

real(dp) function ef_n2__FPL2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = y(10)*(cos(4.0_dp*(-0.25_dp*y(13) + y(26))))**2*sign(1.0_dp, y(26))
end function

real(dp) function ef_n3_d2_FPL2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -4.0_dp*y(10)*sin(8.0_dp*(-0.25_dp*y(13) + y(26)))*sign(1.0_dp, y(26))*y(27)
end function

real(dp) function ef_n5__FTL(y) result(res)
  real(dp), intent(in) :: y(*)
  if (y(46) > 0.0_dp) then
    res = y(14)*(-y(12) + y(38))
  else
    res = 0.0_dp
  end if
end function

real(dp) function ef_n7_d0_FTL(y) result(res)
  real(dp), intent(in) :: y(*)
  if (y(46) > 0.0_dp) then
    res = y(14)*y(43)
  else
    res = 0.0_dp
  end if
end function

real(dp) function ef_n7_d1_FTL(y) result(res)
  real(dp), intent(in) :: y(*)
  if (y(46) > 0.0_dp) then
    res = y(14)*y(44)
  else
    res = 0.0_dp
  end if
end function

real(dp) function ef_n8_d2_FTL(y) result(res)
  real(dp), intent(in) :: y(*)
  if (y(46) > 0.0_dp) then
    res = y(14)*y(53)
  else
    res = 0.0_dp
  end if
end function

real(dp) function ef_n3__LTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = sqrt(((y(20))**2 + (y(23))**2))
end function

real(dp) function ef_n4_d0_LTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(21)*y(20))/(sqrt(((y(20))**2 + (y(23))**2)))
end function

real(dp) function ef_n4_d1_LTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(24)*y(23))/(sqrt(((y(20))**2 + (y(23))**2)))
end function

real(dp) function ef_n5_d2_LTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((y(20)*y(22) + y(23)*y(25)))/(sqrt(((y(20))**2 + (y(23))**2)))
end function

real(dp) function ef_n16__M(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-(y(18) - y(28))*sin(y(3))*y(41) + (y(18) - y(19))*sin((y(3) + y(47)))*y(52) + (y(18) - y(28))*cos(y(3))*y(42))
end function

real(dp) function ef_n13_d0_M(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((y(18) - y(19))*sin((y(3) + y(47)))*y(61) + (y(18) - y(19))*cos((y(3) + y(47)))*y(52)*y(57))
end function

real(dp) function ef_n13_d1_M(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((y(18) - y(19))*sin((y(3) + y(47)))*y(62) + (y(18) - y(19))*cos((y(3) + y(47)))*y(52)*y(58))
end function

real(dp) function ef_n30_d2_M(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-((y(18) - y(28))*sin(y(3))*y(50) + (y(18) - y(28))*sin(y(3))*y(42) + (y(18) - y(28))*cos(y(3))*y(41) + cos(y(3))*y(42)*y(29)) + (1.0_dp + y(64))*(y(18) - y(19))*cos((y(3) + y(47)))*y(52) + (y(18) - y(19))*sin((y(3) + y(47)))*y(66) + (y(18) - y(28))*cos(y(3))*y(51) + sin(y(3))*y(41)*y(29))
end function

real(dp) function ef_n11__RX(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-y(52)*y(45) + y(41))
end function

real(dp) function ef_n11_d0_RX(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(52)*y(55) + y(61)*y(45))
end function

real(dp) function ef_n11_d1_RX(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(52)*y(56) + y(62)*y(45))
end function

real(dp) function ef_n17_d2_RX(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-(y(52)*y(63) + y(66)*y(45)) + y(50))
end function

real(dp) function ef_n11__RY(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-y(52)*y(48) + y(42))
end function

real(dp) function ef_n11_d0_RY(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(52)*y(59) + y(61)*y(48))
end function

real(dp) function ef_n11_d1_RY(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(52)*y(60) + y(62)*y(48))
end function

real(dp) function ef_n17_d2_RY(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-(y(52)*y(65) + y(66)*y(48)) + y(51))
end function

real(dp) function ef_n1__XTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(1) + (y(18) - y(19))*cos(y(3)))
end function

real(dp) function ef_n1_d0_XTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = 1.0_dp
end function

real(dp) function ef_n1_d2_XTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(18) - y(19))*sin(y(3))
end function

real(dp) function ef_n1__YTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(2) - (y(18) - y(19))*sin(y(3)))
end function

real(dp) function ef_n1_d1_YTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = 1.0_dp
end function

real(dp) function ef_n1_d2_YTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(18) - y(19))*cos(y(3))
end function

real(dp) function ef_n1__a2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-3.141592653589793_dp + modulo((3.141592653589793_dp + y(3)), 6.283185307179586_dp))
end function

real(dp) function ef_n1_d2_a2(y) result(res)
  real(dp), intent(in) :: y(*)
  res = 1.0_dp
end function

real(dp) function ef_n4__cosg(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(20))/(y(38))
end function

real(dp) function ef_n6_d0_cosg(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((-y(43)*y(20) + y(38)*y(21)))/((y(38))**2)
end function

real(dp) function ef_n6_d1_cosg(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(44)*y(20))/((y(38))**2)
end function

real(dp) function ef_n7_d2_cosg(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((-y(53)*y(20) + y(38)*y(22)))/((y(38))**2)
end function

real(dp) function ef_n4__dLTL(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-y(12) + y(38))
end function

real(dp) function ef_n12_d0_f4(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(68))/(y(15))
end function

real(dp) function ef_n12_d1_f4(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(69))/(y(15))
end function

real(dp) function ef_n18_d2_f4(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(84))/(y(15))
end function

real(dp) function ef_n12_d0_f5(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(71))/(y(15))
end function

real(dp) function ef_n12_d1_f5(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(72))/(y(15))
end function

real(dp) function ef_n18_d2_f5(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(85))/(y(15))
end function

real(dp) function ef_n14_d0_f6(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(79))/(y(11))
end function

real(dp) function ef_n14_d1_f6(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(80))/(y(11))
end function

real(dp) function ef_n31_d2_f6(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(89))/(y(11))
end function

real(dp) function ef_n4__gamma(y) result(res)
  real(dp), intent(in) :: y(*)
  res = atan(y(39))
end function

real(dp) function ef_n6_d0_gamma(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(49))/((1.0_dp + (y(39))**2))
end function

real(dp) function ef_n6_d1_gamma(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(40))/((1.0_dp + (y(39))**2))
end function

real(dp) function ef_n7_d2_gamma(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(54))/((1.0_dp + (y(39))**2))
end function

real(dp) function ef_n4__sing(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(23))/(y(38))
end function

real(dp) function ef_n6_d0_sing(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(43)*y(23))/((y(38))**2)
end function

real(dp) function ef_n6_d1_sing(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((-y(44)*y(23) + y(38)*y(24)))/((y(38))**2)
end function

real(dp) function ef_n7_d2_sing(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((-y(53)*y(23) + y(38)*y(25)))/((y(38))**2)
end function

real(dp) function ef_n3__tang(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(23))/(y(20))
end function

real(dp) function ef_n4_d0_tang(y) result(res)
  real(dp), intent(in) :: y(*)
  res = -(y(21)*y(23))/((y(20))**2)
end function

real(dp) function ef_n3_d1_tang(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(24))/(y(20))
end function

real(dp) function ef_n5_d2_tang(y) result(res)
  real(dp), intent(in) :: y(*)
  res = ((-y(22)*y(23) + y(20)*y(25)))/((y(20))**2)
end function

real(dp) function ef_n1__xCP(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (y(16)*(cos(0.5_dp*y(3)))**2 + (y(17) - 0.5_dp*y(16))*(sin(y(3)))**2)
end function

real(dp) function ef_n1_d2_xCP(y) result(res)
  real(dp), intent(in) :: y(*)
  res = (-0.5_dp*y(16)*sin(y(3)) + (y(17) - 0.5_dp*y(16))*sin(2.0_dp*y(3)))
end function

end