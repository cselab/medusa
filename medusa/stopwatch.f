subroutine stopwatch(oper)
! This subrotuine computes elapsed cpu time after a second call.
   
  character(len = 4), intent(in) :: oper
! real time
!  real :: dtime, times(2), cputime
!  external dtime
!  integer :: realtid, times, tms(4), tstart, tstop
! save tstart
!!character(len = 10) :: date, time
!    integer :: realtid, tms(4), tstart, tstop
        real :: tstart, tstop
   save tstart

   if (oper == 'STAR' .or. oper == 'star') then
!       realtid = times(tms)
!       tstart = tms(1)+tms(2)+tms(3)+tms(4)
        call cpu_time(tstart)
   elseif (oper == 'STOP' .or. oper == 'stop') then
!       realtid = times(tms)
        call cpu_time(tstop)
!       tstop = tms(1)+tms(2)+tms(3)+tms(4) - tstart
        tstop = tstop - tstart
        write (*, '("Elapsed cpu time: ", f16.2, " seconds")') tstop
   else
      write (*, *)
      write (*, '("Error: in stopwatch.")')
      stop
   endif
end subroutine stopwatch
