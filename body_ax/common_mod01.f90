module common_mod01

  use parameters_mod, only : MAX_SEGMENTS, MAX_ELEMS
  use precision_mod
  implicit none

!      common/xxx01/Iflow,NSG,NGL,NE,Itp
!      common/xxx02/xw,yw,tw
!      common/xxx03/actis,xcntr,ycntr
!      common/xxx04/Vx,cr,Xlvr,Ylvr

! iflow == 50 => sphere; iflow == 51 => thorus

! NSG: Number of segments defining the body contour
! NGL: Number of Gaussian points for integration
!      over each element
! NE(i): Number of elements on the ith segment
! Itp(i): Index for shape of the ith segment:
!         1 for a straight segment
!         2 for a circular arc
! (Xe, Ye):  end-nodes of elements on a segment
! (Xm, Ym):  mid-nodes of elements on a segment
! (Xw, Yw):  end-nodes of elements on all segments

  integer                                               :: iflow, nsg, ngl
  character(LEN=128)                                    :: ifname
  integer,       dimension(MAX_SEGMENTS)                :: ne, itp
  real(FP_KIND), dimension(MAX_SEGMENTS, MAX_ELEMS + 1) :: xw, yw, tw
  real(FP_KIND), dimension(MAX_SEGMENTS)                :: actis, xcntr, ycntr
  real(FP_KIND)                                         :: vx, cr, xlvr,ylvr

end module common_mod01
