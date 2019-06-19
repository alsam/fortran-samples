
module Gauss_Legendre

    use precision_mod
    use parameters_mod
    implicit none

    integer                                         :: num_points  
    real(FP_KIND), dimension(MAX_QUADRATURE_POINTS) :: zz, ww

contains

subroutine Gauss_Legendre_Init (n)
      
    implicit none
    integer, intent(inout) :: n

    call Gauss_Legendre_Table (n, zz, ww)
    num_points = n

end subroutine Gauss_Legendre_Init

subroutine Gauss_Legendre_Table (n,z,w)

!-----------------------------------------
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! fortran 90 version by Alexander Samoilov
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!----------------------------------------

!------------------------------------------------
! This program accompanies the book:
!
!   C. Pozrikidis
! ``Numerical Computation in Science and Engineering''
!   Oxford University Press, 1998
!------------------------------------------------

!------------------------------------------------
! Abscissas and weights for the Gauss-Legendre
! quadrature with N points
!
! This table contains values for
!
!  N = 1,2,3,4,5,6,8,12
!
!  Default value is 20
!------------------------------------------------

      use precision_mod
      use parameters_mod
      implicit none

      integer, intent(IN OUT) :: n
      real(FP_KIND), dimension(MAX_QUADRATURE_POINTS), intent(OUT) :: z,w


      select case (n)

!--------------------
      case (1)
!--------------------

        Z(1) = 0.0_FP_KIND

        W(1) = 2.0_FP_KIND

!-------------------------
      case (2)
!-------------------------

        Z(1) = -0.57735026918962576450
        Z(2) = -Z(1)

        W(1:2) = [ 1.0_FP_KIND, 1.0_FP_KIND ]

!-------------------------
      case (3)
!-------------------------

        Z(1:2) = [ -0.77459666924148337703_FP_KIND, 0.0_FP_KIND ]
        Z(3) =     -Z(1)

        W(1:2) = [ 0.555555555555555555555_FP_KIND, 0.888888888888888888888_FP_KIND ]
        W(3) =      W(1)

!-------------------------
      case (4)
!-------------------------

        Z(1:2) = [ -0.86113631159405257522_FP_KIND, -0.33998104358485626480_FP_KIND ]
        Z(3:4) =   -Z(2:1:-1)

        W(1:2) = [ 0.34785484513745385737_FP_KIND,   0.65214515486254614262_FP_KIND ]
        W(3:4) =    W(2:1:-1)

!-------------------------
      case (5)
!-------------------------

        Z(1:3) = [ -0.90617984593866399279_FP_KIND, -0.53846931010568309103_FP_KIND, 0.0_FP_KIND ]
        Z(4:5) =   -Z(2:1:-1)

        W(1:3) = [ 0.23692688505618908751_FP_KIND,   0.47862867049936646804_FP_KIND, 0.56888888888888888889_FP_KIND ]
        W(4:5) =    W(2:1:-1)

!-------------------------
      case (6)
!-------------------------

        Z(1:3) = [ -0.932469514203152_FP_KIND, -0.661209386466265_FP_KIND, -0.238619186083197_FP_KIND ]
        Z(4:6) =   -Z(3:1:-1)

        W(1:3) = [ 0.171324492379170_FP_KIND,   0.360761573048139_FP_KIND,  0.467913934572691_FP_KIND ]
        W(4:6) =    W(3:1:-1)

!-------------------------
      case (8)
!-------------------------

        Z(1:4) = [ -0.96028985649753623168_FP_KIND, -0.79666647741362673959_FP_KIND, &
                   -0.52553240991632898581_FP_KIND, -0.18343464249564980493_FP_KIND ]
        Z(5:8) =   -Z(4:1:-1)

        W(1:4) = [  0.10122853629037625915_FP_KIND,  0.22238103445337447054_FP_KIND, &
                    0.31370664587788728733_FP_KIND,  0.36268378337836198296_FP_KIND ]
        W(5:8) =    W(4:1:-1)

!--------------------------
      case (12)
!--------------------------

        Z(1:6) = [ -0.981560634246719_FP_KIND, -0.904117256370475_FP_KIND, -0.769902674194305_FP_KIND, &
                   -0.587317954286617_FP_KIND, -0.367831498998180_FP_KIND, -0.125233408511469_FP_KIND ]

        Z(7:12) =  -Z(6:1:-1)

        W(1:6) = [  0.047175336386511_FP_KIND,  0.106939325995318_FP_KIND,  0.160078328543346_FP_KIND, &
                    0.203167426723066_FP_KIND,  0.233492536538355_FP_KIND,  0.249147045813403_FP_KIND ]

        W(7:12) =   W(6:1:-1)

!---------------------------
!     case (20)
!---------------------------

!-----
! trap
!-----

      case default ! e.g. n .ne. 1..6, 8, 12, OR == 20

        if (n /= 20) then

          write (*,*)
          write (*,*) ' Gauss_Legendre:'
          write (*,*)
          write (*,*) '   Chosen number of Gaussian points'
          write (*,*) '   is not available; Will take N=20'

          n = 20
        endif



        Z(1:10) = [ -0.993128599185094924786_FP_KIND, -0.963971927277913791268_FP_KIND, &
                    -0.912234428251325905868_FP_KIND, -0.839116971822218823395_FP_KIND, &
                    -0.746331906460150792614_FP_KIND, -0.636053680726515025453_FP_KIND, &
                    -0.510867001950827098004_FP_KIND, -0.373706088715419560673_FP_KIND, &
                    -0.227785851141645078080_FP_KIND, -0.076526521133497333755_FP_KIND ]

        Z(11:20) =  -Z(10:1:-1)
       
        W(1:10) = [  0.017614007139152118312_FP_KIND,  0.040601429800386941331_FP_KIND, &
                     0.062672048334109063570_FP_KIND,  0.083276741576704748725_FP_KIND, &
                     0.101930119817240435037_FP_KIND,  0.118194531961518417312_FP_KIND, &
                     0.131688638449176626898_FP_KIND,  0.142096109318382051329_FP_KIND, &
                     0.149172986472603746788_FP_KIND,  0.152753387130725850698_FP_KIND ]

        W(11:20) =   W(10:1:-1)

!-----------
      end select
!-----------

      return
    end subroutine Gauss_Legendre_Table

end module Gauss_Legendre

