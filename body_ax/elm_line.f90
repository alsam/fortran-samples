      subroutine elm_line &
 
        (N &
        ,ratio &
        ,X1,Y1 &
        ,X2,Y2 &
        ,sinit &
        ,Isym &
        ,Xe,Ye,se &
        ,Xm,Ym,sm &
        )

!-----------------------------------------
! FDLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved.
!
! This program is to be used only under the
! stipulations of the licensing agreement.
!----------------------------------------

!-----------------------------------------------------
! Disretization of a line segment into N elements
!
!  X1,Y1: coordinates of the first point
!  X2,Y2: coordinates of the last point
!
!  ratio:
!
!     If Isym = 0, ratio of length of LAST to FIRST element
!     If Isym = 1, ratio of length of MID  to FIRST element
!
!  alpha: geometric factor ratio
!
!  sinit: specified arc length at (X1, Y1)
!
!  se: arc length at the element end-nodes
!  sm: arc length at the element mid-nodes
!
!  xe,ye: end nodes
!  xm,ym: mid nodes
!
!-----------------------------------------------------

      use precision_mod
      use parameters_mod
      implicit none

      ! parameters
      integer,                                 intent(in)  :: n, isym
      real(FP_KIND),                           intent(in)  :: ratio, x1, y1, x2, y2, sinit
      real(FP_KIND), dimension(MAX_ELEMS + 1), intent(out) :: xe, ye, se
      real(FP_KIND), dimension(MAX_ELEMS),     intent(out) :: xm, ym, sm

      ! internal variables
      integer       :: n1, nh, nh1, i
      real(FP_KIND) :: alpha, deltax, deltay, factor, tmp1, tmp2, texp, xh, yh

!------------
! one element
!------------

      if (N == 1) then
       xe(1) = X1
       ye(1) = Y1
       xe(2) = X2
       ye(2) = Y2
       se(1) = sinit
       se(2) = se(1)+sqrt( (X2-X1)**2+(Y2-Y1)**2)
       Go to 99
      end if

!--------
! prepare
!--------

      N1 = N+1

!--------------------
! biased distribution
!--------------------

      if (Isym == 0) then

      if (ratio ==  ONE) then
        alpha  = ONE
        factor = ONE/N
      else
        texp   = ONE/(N-ONE)
        alpha  = ratio**texp
        factor = (ONE-alpha)/(ONE-alpha**N)
      end if

      deltax = (x2-x1) * factor   ! x length of first element
      deltay = (y2-y1) * factor   ! y length of first element

      Xe(1) = X1     ! first point
      Ye(1) = Y1
      se(1) = sinit

      do i=2,N1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax*alpha
        deltaY = deltay*alpha
      end do

      Go to 99

      end if

!----------------------
! symmetric distribution
! even number of points
!----------------------

      If(mod(N,2) == 0) then

      Xh = HALF*(X1+X2)      ! mid-point
      Yh = HALF*(Y1+Y2)

      If(N == 2) then
       xe(1) = X1
       ye(1) = Y1
       xe(2) = Xh
       ye(2) = Yh
       xe(3) = X2
       ye(3) = Y2
       se(1) = sinit
       se(2) = se(1)+sqrt( (Xh-X1)**2+(Yh-Y1)**2)
       se(3) = se(2)+sqrt( (X2-Xh)**2+(Y2-Yh)**2)
       Go to 99
      End If

      Nh  = N/2
      Nh1 = Nh+1

      If(ratio == ONE) then
        alpha  = ONE
        factor = ONE/Nh
      Else
        texp   = ONE/(Nh-ONE)
        alpha  = ratio**texp
        factor = (ONE-alpha)/(ONE-alpha**Nh)
      End If

      deltax = (xh-x1) * factor   ! x length of first element
      deltay = (yh-y1) * factor   ! y length of first element

      Xe(1) = X1    ! first point
      Ye(1) = Y1
      se(1) = sinit

      Do i=2,Nh1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax*alpha
        deltaY = deltay*alpha
      End Do

      deltax = deltax/alpha
      deltaY = deltay/alpha

      Do i=Nh1+1,N1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax/alpha
        deltay = deltay/alpha
      End Do

      Go to 99

      End If

!-----------------------
! symmetric distribution
! odd number of points
!-----------------------

      If(ratio == ONE) then
        alpha  = ONE
        factor = ONE/N1
      Else
        texp   = TWO/(N-ONE)
        alpha  = ratio**texp
        tmp1   = HALF*(N+ONE)
        tmp2   = HALF*(N-ONE)
        factor = (ONE-alpha)/(TWO-alpha**tmp1-alpha**tmp2)
      End If

      deltax = (x2-x1) * factor   ! x length of first element
      deltay = (y2-y1) * factor   ! y length of first element

      Xe(1) = X1     ! first point
      Ye(1) = Y1
      se(1) = sinit

      Do i=2,(N+3)/2
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax*alpha
        deltaY = deltay*alpha
      End Do

      deltax = deltax/(alpha**2)
      deltay = deltay/(alpha**2)

      Do i=(N+5)/2,N1
        Xe(i)  = Xe(i-1)+deltax
        Ye(i)  = Ye(i-1)+deltaY
        se(i)  = se(i-1)+sqrt(deltax**2+deltay**2)
        deltax = deltax/alpha
        deltaY = deltay/alpha
      End Do

!-----
! Done
!-----

  99  Continue

!-------------------
! compute mid-points
!-------------------

      Do i=1,N
       xm(i) = HALF*(xe(i)+xe(i+1))
       ym(i) = HALF*(ye(i)+ye(i+1))
       sm(i) = HALF*(se(i)+se(i+1))
      End Do

!-----
! Done
!-----

      write (6,*) "elm_line: Geometric ratio: ",alpha

      Return
      End
