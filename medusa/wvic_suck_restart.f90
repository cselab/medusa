!-------------------------------------------------------------------------------
!* filename: wvic_suck_restart                                                 *!
!* project : ppm                                                              *!
!* purpose : read a restart file                                              *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Sep 22 10:15:17 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!-------------------------------------------------------------------------------
! $Log: wvic_suck_restart.F,v $
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.1  2005/09/28 11:40:45  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------

SUBROUTINE wvic_suck_restart (info)
  
  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE
  
  INTEGER, INTENT(inout) :: info
  CHARACTER(len=256) :: msg
  INTEGER            :: irfile
  CHARACTER(len=128) :: rfile
  INTEGER            :: i,j,ios,istat
  LOGICAL            :: lexist
  
  
  info = 0
  !-----------------------------------------------------------------------------
  ! construct filename
  !-----------------------------------------------------------------------------
  IF(verbose) CALL ppm_write(rank,'wvic_suck_restart','entered',info)
  
  WRITE(rfile,'(A,A,I5.5,A)') &
       &runtag(1:iruntag),'_R',rank,'.rst'
  WRITE(msg,*) 'trying to suck ',rfile
  IF(verbose) CALL ppm_write(rank,'wvic_suck_restart',msg,info)
  
  irfile = 20 + rank
  
  INQUIRE(file=rfile,exist=lexist,iostat=ios)
  IF(.NOT.lexist) THEN
     info = -1
     GOTO 9999
  END IF
  
  OPEN(irfile,file=rfile,form='unformatted')

  READ(irfile) np
  READ(irfile) dt,time,itime,dt_max

  ALLOCATE(xp(dime,np),wp(lda,np),stat=istat)
  IF(istat.NE.0) THEN
     IF(verbose) THEN
        CALL ppm_write(rank,'wvic_suck_restart','failed to allocate ps',info)
     END IF
     info = -1
     GOTO 9999
  END IF
  
  READ(irfile) xp
  READ(irfile) wp
  max_vorticity = SQRT(MAXVAL(wp(1,:)**2 + wp(2,:)**2 + wp(3,:)**2))
  CLOSE(irfile)
  
9999 CONTINUE
  IF(info.NE.0) THEN
     IF(verbose) THEN
      CALL ppm_write(rank,'wvic_suck_restart',&
        &             'failed to open restart file',istat)
     END IF
  END IF
  IF(verbose) CALL ppm_write(rank,'wvic_suck_restart','complete',istat)
END SUBROUTINE wvic_suck_restart
