!------------------------------------------------------------------------------!
!* filename: wvic_field2netcdf                                                 *!
!* project : ppm                                                              *!
!* purpose : write field into a raw binary file (one file per rank)            *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 19:42:31 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  netcdf dependency removed: each rank dumps its own subdomain to a simple
!  stream-access binary file.  Layout (single precision):
!     INTEGER  :: nx, ny, nz
!     REAL(ms) :: grid(2,3)
!     REAL(ms) :: time, dt
!     REAL(ms) :: state(3,nx,ny,nz)
!  The matching reader is wvic_netcdf2field.
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic FIELD2NETCDF =
! dump stuff to a raw binary file
!------------------------------------------------------------------------------!
SUBROUTINE wvic_field2netcdf (info)

  USE module_wvic
  USE ppm_module_write
  USE MPI

  !----------------------------------------------------------------------------!
  ! arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(inout)    :: info

  !----------------------------------------------------------------------------!
  ! localities
  !----------------------------------------------------------------------------!
  INTEGER, PARAMETER :: ms = kind(1.0e0)
  CHARACTER(len=256) :: ncfile
  INTEGER            :: ncid
  CHARACTER(len=256) :: msg
  INTEGER            :: stat, isub, isubl
  INTEGER            :: i,j,k
  INTEGER            :: lnx, lny, lnz
  INTEGER, PARAMETER :: vsize = 3
  REAL(ms), DIMENSION(2,3) :: grid
  REAL(ms), DIMENSION(:,:,:,:), ALLOCATABLE :: state
  !----------------------------------------------------------------------------!
  ! construct filename
  !----------------------------------------------------------------------------!
  WRITE(ncfile,'(A,A,I5.5,A,I5.5,A)') &
       &runtag(1:iruntag),'R',rank,'I',itime,'.nc'

  isub = 1; isubl = isublist(isub)
  lnx = ndata(1,isubl); lny = ndata(2,isubl); lnz = ndata(3,isubl)

  !----------------------------------------------------------------------------!
  ! grid data
  !----------------------------------------------------------------------------!
  grid(1,1) = REAL(min_sub(3,isubl),MS)
  grid(1,2) = REAL(min_sub(2,isubl),MS)
  grid(1,3) = REAL(min_sub(1,isubl),MS)
  grid(2,1) = REAL(dZ);grid(2,2) = REAL(dy); grid(2,3) = REAL(dx)

  !----------------------------------------------------------------------------!
  ! state vector data
  !----------------------------------------------------------------------------!
  ALLOCATE(state(vsize,lnx,lny,lnz),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'field dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_field2netcdf',msg,info)
     GOTO 9999
  END IF
  DO k=1,lnz
     DO j=1,lny
        DO i=1,lnx
           state(1,i,j,k) = REAL(field_wp(1,i,j,k,isub),ms)
           state(2,i,j,k) = REAL(field_wp(2,i,j,k,isub),ms)
           state(3,i,j,k) = REAL(field_wp(3,i,j,k,isub),ms)
        END DO
     END DO
  END DO

  !----------------------------------------------------------------------------!
  ! write everything in one shot
  !----------------------------------------------------------------------------!
  ncid = 71
  OPEN(unit=ncid, file=ncfile, form='unformatted', access='stream', &
       & status='replace', action='write', iostat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'could not open ', TRIM(ncfile)
     CALL ppm_write(rank,'wvic_field2netcdf',msg,info)
     DEALLOCATE(state,stat=stat)
     GOTO 9999
  END IF
  WRITE(ncid) lnx, lny, lnz
  WRITE(ncid) grid
  WRITE(ncid) REAL(time,ms), REAL(dt,ms)
  WRITE(ncid) state
  CLOSE(ncid)

  DEALLOCATE(state,stat=stat)

  IF(verbose) THEN
     WRITE(msg,*) 'data written'
     CALL ppm_write(rank,'wvic_field2netcdf',msg,info)
  END IF

9999 CONTINUE
  WRITE(msg,*) 'complete'

END SUBROUTINE wvic_field2netcdf
