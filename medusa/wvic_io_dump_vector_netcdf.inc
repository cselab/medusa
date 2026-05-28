!-------------------------------------------------------------------------------
!* filename: wvic_io_dump_vector_netcdf                                       *!
!* project : wvic                                                             *!
!* purpose : dump a vector field to a raw binary file                         *!
!*         :                                                                  *!
!* author  : Philippe Chatelain                                               *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!  netcdf dependency removed.  Layout (single precision):
!     INTEGER  :: nx, ny, nz
!     REAL(ms) :: grid(2,3)
!     REAL(ms) :: time, dt
!     REAL(ms) :: state(3,nx,ny,nz)
!  For complex fields the real part is written to <name>_rp_... and the
!  imaginary part to <name>_ip_... .
!------------------------------------------------------------------------------!

#if __WHAT == __REAL
SUBROUTINE wvic_io_dump_vector_netcdf(vectfield, n_data, vectname, iter, info)
#elif __WHAT == __COMPLEX
SUBROUTINE wvic_io_dump_vectorc_netcdf(vectfield, n_data, vectname, iter, info)
#endif

  USE module_wvic
  USE ppm_module_write

  !-----------------------------------------------------------------------------
  ! arguments
  !-----------------------------------------------------------------------------
#if __WHAT == __REAL
  REAL(mk), DIMENSION(:,:,:,:,:),  POINTER :: vectfield
#elif __WHAT == __COMPLEX
  COMPLEX(mk), DIMENSION(:,:,:,:,:),  POINTER :: vectfield
#endif
  INTEGER, DIMENSION(:,:), POINTER         :: n_data
  CHARACTER(len=256)                       :: vectname
  INTEGER, INTENT(IN)                      :: iter
  INTEGER, INTENT(inout)                   :: info

  !-----------------------------------------------------------------------------
  ! localities
  !-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: ms = kind(1.0e0)
  CHARACTER(len=256) :: ncfile
  INTEGER            :: ncid
  CHARACTER(len=256) :: msg
  INTEGER            :: stat, isub, isubl
  INTEGER            :: i,j,k,ivectname
  INTEGER            :: lnx, lny, lnz
  INTEGER, PARAMETER :: vsize = 3
  REAL(ms), DIMENSION(2,3) :: grid
  REAL(ms), DIMENSION(:,:,:,:), ALLOCATABLE :: state

  !-----------------------------------------------------------------------------
  ! construct filename
  !-----------------------------------------------------------------------------
  ivectname = LEN_TRIM(vectname)
#if __WHAT == __REAL
  WRITE(ncfile,'(A,A,A,A,I5.5,A,I5.5,A)') &
       &vectname(1:ivectname),'_',runtag(1:iruntag),'R',rank,'I',iter,'.nc'
#elif __WHAT == __COMPLEX
  WRITE(ncfile,'(A,A,A,A,I5.5,A,I5.5,A)') &
       &vectname(1:ivectname),'_rp_',runtag(1:iruntag),'R',rank,'I',iter,'.nc'
#endif

  isub = 1; isubl = isublist(isub)
  lnx = n_data(1,isubl); lny = n_data(2,isubl); lnz = n_data(3,isubl)

  !-----------------------------------------------------------------------------
  ! grid data
  !-----------------------------------------------------------------------------
  grid(1,1) = REAL(min_sub(3,isubl),MS)
  grid(1,2) = REAL(min_sub(2,isubl),MS)
  grid(1,3) = REAL(min_sub(1,isubl),MS)
  grid(2,1) = REAL(dZ);grid(2,2) = REAL(dy); grid(2,3) = REAL(dx)

  !-----------------------------------------------------------------------------
  ! state vector data
  !-----------------------------------------------------------------------------
  ALLOCATE(state(vsize,lnx,lny,lnz),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'vector dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_dump_vector_field',msg,info)
     GOTO 9999
  END IF
  DO k=1,lnz
     DO j=1,lny
        DO i=1,lnx
#if __WHAT == __REAL
           state(1,i,j,k) = REAL(vectfield(1,i,j,k,isub),ms)
           state(2,i,j,k) = REAL(vectfield(2,i,j,k,isub),ms)
           state(3,i,j,k) = REAL(vectfield(3,i,j,k,isub),ms)
#elif __WHAT == __COMPLEX
           state(1,i,j,k) = REAL(vectfield(1,k,j,i,isub),ms)
           state(2,i,j,k) = REAL(vectfield(2,k,j,i,isub),ms)
           state(3,i,j,k) = REAL(vectfield(3,k,j,i,isub),ms)
#endif
        END DO
     END DO
  END DO

  ncid = 72
  OPEN(unit=ncid, file=ncfile, form='unformatted', access='stream', &
       & status='replace', action='write', iostat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'could not open ', TRIM(ncfile)
     CALL ppm_write(rank,'wvic_io_dump_vector_netcdf',msg,info)
     DEALLOCATE(state,stat=stat)
     GOTO 9999
  END IF
  WRITE(ncid) lnx, lny, lnz
  WRITE(ncid) grid
  WRITE(ncid) REAL(time,ms), REAL(dt,ms)
  WRITE(ncid) state
  CLOSE(ncid)
  DEALLOCATE(state,stat=stat)

#if __WHAT == __COMPLEX
  !-----------------------------------------------------------------------------
  ! imaginary part
  !-----------------------------------------------------------------------------
  ivectname = LEN_TRIM(vectname)
  WRITE(ncfile,'(A,A,A,A,I5.5,A,I5.5,A)') &
       &vectname(1:ivectname),'_ip_',runtag(1:iruntag),'R',rank,'I',iter,'.nc'

  ALLOCATE(state(vsize,lnx,lny,lnz),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'vector dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_dump_vector_field',msg,info)
     GOTO 9999
  END IF
  DO k=1,lnz
     DO j=1,lny
        DO i=1,lnx
           state(1,i,j,k) = REAL(AIMAG(vectfield(1,k,j,i,isub)),ms)
           state(2,i,j,k) = REAL(AIMAG(vectfield(2,k,j,i,isub)),ms)
           state(3,i,j,k) = REAL(AIMAG(vectfield(3,k,j,i,isub)),ms)
        END DO
     END DO
  END DO
  WRITE(*,*) 'Maxval(state): ',MAXVAL(REAL(state))

  ncid = 72
  OPEN(unit=ncid, file=ncfile, form='unformatted', access='stream', &
       & status='replace', action='write', iostat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'could not open ', TRIM(ncfile)
     CALL ppm_write(rank,'wvic_io_dump_vector_netcdf',msg,info)
     DEALLOCATE(state,stat=stat)
     GOTO 9999
  END IF
  WRITE(ncid) lnx, lny, lnz
  WRITE(ncid) grid
  WRITE(ncid) REAL(time,ms), REAL(dt,ms)
  WRITE(ncid) state
  CLOSE(ncid)
  DEALLOCATE(state,stat=stat)
#endif

  IF (verbose) THEN
     WRITE(msg,*) 'data written'
     CALL ppm_write(rank,'wvic_io_dump_vector_netcdf',msg,info)
  END IF

9999 CONTINUE

#if __WHAT == __REAL
END SUBROUTINE wvic_io_dump_vector_netcdf
#elif __WHAT == __COMPLEX
END SUBROUTINE wvic_io_dump_vectorc_netcdf
#endif
