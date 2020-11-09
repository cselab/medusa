! mpif77 main.f90 -lfftw3f
program main
  implicit none
  include "fftw3.f"
  integer :: Plan
  call sfftw_destroy_plan(Plan)
end program main
