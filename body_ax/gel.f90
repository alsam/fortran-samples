      subroutine gel &
 
       (n  &    ! system size &
       ,A  &    ! coefficient matrix &
       ,rhs &   ! right-hand side &
       ,x &
       ,Isym &
       ,Iwlpvt & !    + ,l,u
       ,det &
       ,Istop &
       )

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
! C. Pozrikidis
! ``Numerical Computation in Science and Engineering'
! Oxford University Press, 1998
!------------------------------------------------

!-----------------------------------------------
!  Gauss elimination with option of row pivoting.
!
!  Algorithm 3.5.1
!
!  This subroutine returns:
!
!  (a) the solution vector,
!  (b) the and lower triangular matrix factors L and U
!  (c) a flag for completion
!  (d) the determinant
!
!  SYMBOLS:
!  --------
!
!   a ...... square matrix
!   n ...... size (rows/columns) of matrix a
!   rhs .... right hand side vector (e.g. b, as in Ax=b)
!   c....... extended matrix
!   x ...... solution vector
!
!   Isym ... 1 if a is symmetric; 0 if nonsymmetric
!   Iwlpvt.. 0 for no pivoting, 1 for pivoting
!
!   eps1..... tolerance to identify a singular matrix
!   tol..... tolerance for the residuals
!
!   l ...... lower triangular matrix
!   u ...... upper triangular matrix
!   det .... determinant (det(a) = +- det(l)*det(u))
!
!   Istop... flag: Istop = 1 if something is wrong
!
!   pivot .. absolute value of pivot candidates
!   Ipv..... location of pivotal element
!   Icount.. counter for number of row interchanges
!
!-------------------------------------------------------------

      use precision_mod
      use parameters_mod
      implicit none

      ! parameters
      integer,                                    intent(in)    :: n, isym
      integer,                                    intent(inout) :: iwlpvt
      integer,                                    intent(out)   :: istop
      real(FP_KIND), dimension(MAX_DIM, MAX_DIM), intent(in)    :: a
      real(FP_KIND), dimension(MAX_DIM),          intent(in)    :: rhs
      real(FP_KIND), dimension(MAX_DIM),          intent(out)   :: x
      real(FP_KIND),                              intent(out)   :: det

      ! internal variables
      real(FP_KIND), dimension(MAX_DIM, MAX_DIM)                :: u, l
      real(FP_KIND), dimension(MAX_DIM, MAX_DIM + 1)            :: c

      real(FP_KIND), parameter  :: eps1=1.0E-8_FP_KIND, tol=1.0E-8_FP_KIND
      integer                   :: i, j, m, ma, m1, na, n1, icount, ipv
      real(FP_KIND)             :: pivot, summ, relax, save1

!----------------------------------------
! Disable pivoting for symmetric systems
!---------------------------------------

      if (isym == 1) then
        write (6,*) " gel: system is symmetric"
        write (6,*) "      pivoting disabled"
        Iwlpvt = 0
      end if

!-----------
! initialize
!-----------

      istop  = 0
      icount = 0     ! counts row interchanges

!--------
! prepare
!--------

      na = n-1
      n1 = n+1

!-------------------
! Initialize l and c
!-------------------

      open (unit = 77, form = 'formatted', status = 'replace', file = 'linear_system.text')
      write (77,'(a,i0,a,i0,a,f18.14)'), (('  a (',i,',',j,') = ', a(i,j),j=1,n),i=1,n)
      write (77,'(a,i0,a,f18.14)'), ('  rhs (',i,') = ', rhs(i), i=1,n)
      close (77)

      ! save a(:,:) and rhs(:) in unformatted file
      !open(unit = 77, form = 'unformatted', status = 'new', file = 'linear_system.data')
      open (unit = 77, form = 'unformatted', status = 'replace', file = 'linear_system.data')
      write (77), n
      write (77), ((a(i,j),j=1,n),i=1,n)
      write (77), (rhs(i), i=1,n)
      close (77)

      l(1:n,1:n) = ZERO
      c(1:n,1:n) = a(1:n,1:n)
      c(1:n,n1)  = rhs(1:n)

!---------------------
! Begin row reductions
!---------------------

      outer_loop: do, m=1,na           ! outer loop for working row

        ma = m-1
        m1 = m+1
 
        !-----------------------------
        ! Pivoting module
        !
        ! begin by searching column i
        ! for largest element
        !----------------------------
  
        pivoting: if (iwlpvt == 1) then
  
            ipv = m
            pivot = abs(c(m,m))
  
            do j=m1,n
              if (abs(c(j,m)) > pivot) then
                ipv = j
                pivot = abs(c(j,m))
              end if
            end do
  
            !write (*,*) 'pivot = ', pivot
  
            if (pivot < eps1) then
              write (6,*)
              write (6,*) " gel: trouble at station 1"
              write (6,*)
              istop = 1
              return
            end if
  
            !--------------------------------------
            ! switch the working row with
            ! the row containing the pivot element
            !
            ! also switch rows in l
            !--------------------------------------
  
            swap_rows: if (ipv /= m) then
  
              do, j=m,n1
                save1    = c(m,j)
                c(m,j)   = c(ipv,j)
                c(ipv,j) = save1
              end do
  
              do, j=1,ma
                save1    = l(m,j)
                l(m,j)   = l(ipv,j)
                l(ipv,j) = save1
              end do
  
              icount = icount+1
  
            end if swap_rows
  
        end if pivoting
 
 !---------------------------------------
 ! reduce column i beneath element c(m,m)
 !---------------------------------------
 
       reduce: do, i=m1,n
 
        if (isym == 1) then        ! symmetric matrix
 
          relax  = c(m,i)/c(m,m)
          l(i,m) = relax
          c(i,m) = ZERO
 
          c(i,i:n1) = c(i,i:n1)-relax*c(m,i:n1)
 
         else                     ! non-symmetric matrix
 
          relax  = c(i,m)/c(m,m)
          l(i,m) = relax
          c(i,m) = ZERO
 
          c(i,m1:n1) = c(i,m1:n1)-l(i,m)*c(m,m1:n1)
 
        end if
 
       end do reduce

      end do outer_loop     ! for working row

!---------------------------------
! check the last diagonal element
! for singularity
!--------------------------------

      if (abs(c(n,n)) < eps1) then

        write (6,*)
        write (6,*) " gel: trouble at station 2"
        write (6,*)
        istop = 1
        return

      end if

!----------------------
! complete the matrix l
!----------------------
      
      do, i=1,n
        l(i,i)=ONE
      end do

!--------------------
! define the matrix u
!--------------------

      u(1:n,1:n) = c(1:n,1:n)

!-----------------------------------
! perform back-substitution to solve
! the reduced system
! using the upper triangular matrix c
!------------------------------------

      x(n) = c(n,n1)/c(n,n)

      do i=na,1,-1
        x(i) = (c(i,n1) - dot_product(c(i,i+1:n), x(i+1:n))) / c(i,i)
      end do

!-----------------------
! compute the determinant as
!
! det(a) = (+-) det(l)*det(u)
!
!-----------------------

      det = ONE

      do, i=1,n
        det=det*c(i,i)
      end do

      if (iwlpvt == 1) then

        write (6,*) " gel: number of row interchanges : ",Icount

        do, i=1,icount
          det = -det
        end do

      end If

!----------------------
! compute the residuals
!----------------------

      do i=1,n

        summ = rhs(i) - dot_product(a(i, 1:n), x(1:n))
        
        if (abs(summ) > tol) then
          Istop = 1
          write (6,*) " gel: problem in solving the linear system"
          write (6,100) i,summ
        end if

      end do

!----------------------
! write to file residual, det, l, u, solution x
!----------------------

      open (unit = 77, form = 'formatted', status = 'replace', file = 'linear_system.solution.text')
      write (77,*), 'residual = ', summ, ' det = ', det
      write (77,'(a,i0,a,i0,a,f18.14)'), (('  l (',i,',',j,') = ', l(i,j),j=1,n),i=1,n)
      write (77,'(a,i0,a,i0,a,f18.14)'), (('  u (',i,',',j,') = ', u(i,j),j=1,n),i=1,n)
      write (77,'(a,i0,a,f18.14)'), ('  x (',i,') = ', x(i), i=1,n)
      close (77)

!-----
! Done
!-----

  100 Format (1x,i4,f15.10)
  101 Format(16(16(1x,f5.3),/))

      Return
      end
