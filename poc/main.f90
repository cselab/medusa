program main
  USE, INTRINSIC :: ISO_C_BINDING
  INCLUDE 'fftw3.f03'
  type(C_PTR) :: Plan
  call sfftw_execute(Plan)
  call dfftw_execute(Plan)
  call fftw_destroy_plan(Plan)
end program main
