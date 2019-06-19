module common_mod02

  use parameters_mod, only : MAX_SEGMENTS, MAX_ELEMS
  use precision_mod
  implicit none

!      common/xxx05/X0,Y0,T0,S0,dphidn0
!      common/xxx06/tnx0,tny0,vnx0,vny0,arel
!      common/xxx07/xcenter,ycenter
!      common/xxx08/xwmin,ywmin,xwmax,ywmax

  real(FP_KIND), dimension(MAX_SEGMENTS * MAX_ELEMS) :: x0, y0, t0, s0, dphidn0, tnx0, tny0, vnx0, vny0, arel
  real(FP_KIND)                                      :: xcenter, ycenter, xwmin, ywmin, xwmax, ywmax

end module common_mod02
