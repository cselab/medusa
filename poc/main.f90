! mpif90 main.f90 -I$HOME/.local/include -L$HOME/.local/lib -lfftw3 -lfftw3f
program main
  USE, INTRINSIC :: ISO_C_BINDING
  INCLUDE 'fftw3.f03'
  type(C_PTR) :: Plan
  call sfftw_execute(Plan)
  call dfftw_execute(Plan)
  call fftw_destroy_plan(Plan)
end program main
