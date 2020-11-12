!-------------------------------------------------------------------------------
!* filename: wvic_write_restart                                                *!
!* project : ppm                                                              *!
!* purpose : write a restart file                                             *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Mon Aug 23 12:37:15 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_write_restart.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.1  2005/09/28 11:40:47  michaebe
!  Fork from ppm_pvc
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  restart writing, based on Ivo Sbalzarini''s pse_write_restart
!-------------------------------------------------------------------------------
SUBROUTINE wvic_write_restart (info)

  USE module_wvic
  USE ppm_module_user_io
  USE ppm_module_data
  USE ppm_module_write

  !-----------------------------------------------------------------------------
  !  Arguments
  !-----------------------------------------------------------------------------
  INTEGER,  INTENT(inout)                :: info
  
  !-----------------------------------------------------------------------------
  !  Localities
  !-----------------------------------------------------------------------------
  INTEGER,  DIMENSION(9 )                :: ivec
  INTEGER                                :: iounit
  REAL(mk), DIMENSION(8 )                :: fvec
  CHARACTER(len=256)                     :: filename

  INCLUDE 'mpif.h'


  WRITE(filename,'(A,I4.4,A)') runtag(1:iruntag),itime,'.restart'

  !-----------------------------------------------------------------------------
  !  open io unit
  !-----------------------------------------------------------------------------
  iounit = 26
  CALL ppm_io_open(iounit,filename,ppm_param_io_write,ppm_param_io_replace,&
       &           ppm_param_io_binary,ppm_param_io_centralized,info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_write_restart','failed to open write unit',info)
     info = -1
     GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  !  write  total number of particles to file
  !-----------------------------------------------------------------------------
  CALL ppm_io(iounit,np,ppm_param_io_write,ppm_param_io_sum,STAT=info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_write_restart','failed to write ntot',info)
     info = -1
     GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  !  write mk,itime,dime to file
  !-----------------------------------------------------------------------------
  ivec(1)   = mk
  ivec(2)   = itime
  ivec(3)   = dime
  ivec(4:6) = nx(1:3)
  ivec(7:9) = ghostsize(1:3)

  CALL ppm_io(iounit,ivec,ppm_param_io_write,ppm_param_io_root,STAT=info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_write_restart','failed to write ivec',info)
     info = -1
     GOTO 9999
  END IF
  
  
  !-----------------------------------------------------------------------------
  !  write time, dt, minphys, maxphys to file
  !-----------------------------------------------------------------------------
  fvec(1)   = dt
  fvec(2)   = time
  fvec(3:5) = min_physg(1:3)
  fvec(6:8) = max_physg(1:3)

  CALL ppm_io(iounit,fvec,ppm_param_io_write,ppm_param_io_root,STAT=info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_write_restart','failed to write fvec',info)
     info = -1
     GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  !  write xp
  !-----------------------------------------------------------------------------
  CALL ppm_io(iounit,xp,ppm_param_io_write, ppm_param_io_concat,STAT=info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_write_restart','failed to write xp',info)
     info = -1
     GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  !  write wp
  !-----------------------------------------------------------------------------
  CALL ppm_io(iounit,wp,ppm_param_io_write, ppm_param_io_concat,STAT=info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_write_restart','failed to write wp',info)
     info = -1
     GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  !  close file
  !-----------------------------------------------------------------------------
  CALL ppm_io_close(iounit,info)

  !-----------------------------------------------------------------------------
  !  Confirm
  !-----------------------------------------------------------------------------
  IF (rank .EQ. 0) THEN
     if(verbose) CALL ppm_write(rank,'pse_write_restart',filename,info)
  ENDIF
  

9999 CONTINUE
  RETURN
END SUBROUTINE wvic_write_restart




  
