!-------------------------------------------------------------------------------
!* filename: wvic_dump_restart                                                 *!
!* project : ppm                                                              *!
!* purpose : dump restart                                                     *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Sep 22 10:04:21 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
! $Log: wvic_dump_restart.F,v $
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.1  2005/09/28 11:40:33  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------


SUBROUTINE wvic_dump_restart (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  INTEGER, INTENT(inout) :: info
  CHARACTER(len=256) :: msg
  INTEGER            :: irfile
  CHARACTER(len=256) :: rfile
  INTEGER            :: i,j
  
  !-----------------------------------------------------------------------------
  ! construct filename
  !-----------------------------------------------------------------------------
  if(verbose) CALL ppm_write(rank,'wvic_dump_restart','entered',info)

  WRITE(rfile,'(A,A,I5.5,A)') &
       &runtag(1:iruntag),'_R',rank,'.rst'

  irfile = 20 + rank
    
  OPEN(irfile,file=rfile,form='unformatted')

  WRITE(irfile) np
  WRITE(irfile) dt,time,itime, dt_max
  
  !-----------------------------------------------------------------------------
  ! write positions
  !-----------------------------------------------------------------------------
  WRITE(irfile) xp
  !-----------------------------------------------------------------------------
  ! write strengths
  !-----------------------------------------------------------------------------
  WRITE(irfile) wp
  CLOSE(irfile)

  if(verbose) CALL ppm_write(rank,'wvic_dump_restart','complete',info)
END SUBROUTINE wvic_dump_restart
