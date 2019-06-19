! This module provides a simple facility for changing between single
! and double precision

module precision_mod
  integer, parameter, public :: SINGLE_PRECISION = kind(0.0) ! Single precision
  integer, parameter, public :: DOUBLE_PRECISION = kind(0.0d0) ! Double precision
  
  ! Comment out one of the lines below
  !integer, parameter, public :: FP_KIND = SINGLE_PRECISION
  integer, parameter, public :: FP_KIND = DOUBLE_PRECISION
end module precision_mod
