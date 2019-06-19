module parameters_mod

  use precision_mod
  implicit none

  !integer,      parameter :: WPREC = selected_real_kind (15, 307)
  integer,       parameter :: MAX_SEGMENTS = 10, MAX_ELEMS = 512, MAX_QUADRATURE_POINTS = 20, MAX_DIM = MAX_SEGMENTS * MAX_ELEMS

  ! some math constants
  real(FP_KIND), parameter :: eps   = 1e-2_FP_KIND, acc = 1e-9_FP_KIND
  real(FP_KIND), parameter :: ZERO  = 0.0_FP_KIND, &
                              HALF  = 0.5_FP_KIND, &
                              ONE   = 1.0_FP_KIND, &
                              TWO   = 2.0_FP_KIND, &
                              THREE = 3.0_FP_KIND, &
                              FOUR  = 4.0_FP_KIND, &
                              FIVE  = 5.0_FP_KIND 
  real(FP_KIND), parameter :: PI    = 3.1415926535797_FP_KIND, &
                              PIH   = HALF  * PI,              &
                              PI2   = TWO   * PI,              &
                              PI4   = FOUR  * PI

end module parameters_mod

