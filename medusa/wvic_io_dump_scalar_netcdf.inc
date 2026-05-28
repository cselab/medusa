!-------------------------------------------------------------------------------
!* filename: wvic_io_dump_scalar_netcdf                                       *!
!* project : ppm                                                              *!
!* purpose : dump a scalar field to a raw binary file                         *!
!*         :                                                                  *!
!* author  : Philippe Chatelain                                               *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!  netcdf dependency removed.  Layout (single precision):
!     INTEGER  :: nx, ny, nz
!     REAL(ms) :: grid(2,3)
!     REAL(ms) :: time, dt
!     REAL(ms) :: state(nx,ny,nz)
!------------------------------------------------------------------------------!


SUBROUTINE wvic_io_dump_scalar_netcdf(scalfield, n_data, scalname, iter,info)

  USE module_wvic
  USE ppm_module_write

  !-----------------------------------------------------------------------------
  ! arguments
  !-----------------------------------------------------------------------------
  REAL(mk), DIMENSION(:,:,:,:),    POINTER :: scalfield
  CHARACTER(len=256)                       :: scalname
  INTEGER, DIMENSION(:,:), POINTER         :: n_data
  INTEGER, INTENT(in)                      :: iter
  INTEGER, INTENT(inout)                   :: info

  !-----------------------------------------------------------------------------
  ! localities
  !-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: ms = kind(1.0e0)
  CHARACTER(len=256) :: ncfile
  INTEGER            :: ncid
  CHARACTER(len=256) :: msg
  INTEGER            :: stat, isub, isubl
  INTEGER            :: i,j,k, iscalname
  INTEGER            :: lnx, lny, lnz
  REAL(ms), DIMENSION(2,3) :: grid
  REAL(ms), DIMENSION(:,:,:), ALLOCATABLE :: state

  !-----------------------------------------------------------------------------
  ! construct filename
  !-----------------------------------------------------------------------------
  iscalname = LEN_TRIM(scalname)
  WRITE(ncfile,'(A,A,A,A,I5.5,A,I5.5,A)')&
       &scalname(1:iscalname),'_',runtag(1:iruntag),'R',rank,'I',iter,'.nc'

  isub = 1; isubl = isublist(isub)
  lnx = n_data(3,isubl); lny = n_data(2,isubl); lnz = n_data(1,isubl)

  !-----------------------------------------------------------------------------
  ! grid data
  !-----------------------------------------------------------------------------
  grid(1,1) = REAL(min_sub(3,isubl),MS)
  grid(1,2) = REAL(min_sub(2,isubl),MS)
  grid(1,3) = REAL(min_sub(1,isubl),MS)
  grid(2,1) = REAL(dZ);grid(2,2) = REAL(dy); grid(2,3) = REAL(dx)

  !-----------------------------------------------------------------------------
  ! state data
  !-----------------------------------------------------------------------------
  ALLOCATE(state(lnx,lny,lnz),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'scalar dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_dump_scalar_field',msg,info)
     GOTO 9999
  END IF
  DO k=1,n_data(3,isubl)
     DO j=1,n_data(2,isubl)
        DO i=1,n_data(1,isubl)
           state(k,j,i) = REAL(scalfield(i,j,k,isub),ms)
        END DO
     END DO
  END DO

  ncid = 72
  OPEN(unit=ncid, file=ncfile, form='unformatted', access='stream', &
       & status='replace', action='write', iostat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'could not open ', TRIM(ncfile)
     CALL ppm_write(rank,'wvic_io_dump_scalar_netcdf',msg,info)
     DEALLOCATE(state,stat=stat)
     GOTO 9999
  END IF
  WRITE(ncid) lnx, lny, lnz
  WRITE(ncid) grid
  WRITE(ncid) REAL(time,ms), REAL(dt,ms)
  WRITE(ncid) state
  CLOSE(ncid)
  DEALLOCATE(state,stat=stat)

  IF (verbose) THEN
      WRITE(msg,*) 'data written'
      CALL ppm_write(rank,'wvic_io_dump_scalar_netcdf',msg,info)
  END IF

9999 CONTINUE
  RETURN
END SUBROUTINE wvic_io_dump_scalar_netcdf
