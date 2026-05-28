!------------------------------------------------------------------------------!
!* filename: wvic_netcdf2field                                                 *!
!* project : ppm                                                              *!
!* purpose : read field from a raw binary file (one file per rank)             *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!  netcdf dependency removed.  Reads the layout written by wvic_field2netcdf.
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic NETCDF2FIELD =
! read state back from a raw binary file
!------------------------------------------------------------------------------!
SUBROUTINE wvic_netcdf2field (info)

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
  REAL(ms)           :: rdt, rtime
  !----------------------------------------------------------------------------!
  ! construct filename
  !----------------------------------------------------------------------------!
  WRITE(ncfile,'(A,A,I5.5,A,I5.5,A)') &
       &runtag(1:iruntag),'R',rank,'I',netcdf_itime,'.nc'

  isub = 1; isubl = isublist(isub)

  !----------------------------------------------------------------------------!
  ! open the data set
  !----------------------------------------------------------------------------!
  ncid = 71
  OPEN(unit=ncid, file=ncfile, form='unformatted', access='stream', &
       & status='old', action='read', iostat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'could not open ', TRIM(ncfile)
     CALL ppm_write(rank,'wvic_netcdf2field',msg,info)
     RETURN
  END IF

  READ(ncid) lnx, lny, lnz

  !----------------------------------------------------------------------------!
  ! allocate temporary memory, ...
  !----------------------------------------------------------------------------!
  ALLOCATE(state(vsize,lnx,lny,lnz),stat=stat)

  !----------------------------------------------------------------------------!
  ! and read the variables
  !----------------------------------------------------------------------------!
  READ(ncid) grid
  READ(ncid) rtime, rdt
  READ(ncid) state
  CLOSE(ncid)

  !----------------------------------------------------------------------------!
  ! copy the contents to field_wp
  !----------------------------------------------------------------------------!
  DO k=1,lnz
     DO j=1,lny
        DO i=1,lnx
           field_wp(1,i,j,k,isub) = REAL(state(1,i,j,k),mk)
           field_wp(2,i,j,k,isub) = REAL(state(2,i,j,k),mk)
           field_wp(3,i,j,k,isub) = REAL(state(3,i,j,k),mk)
        END DO
     END DO
  END DO

  time = REAL(rtime,mk)
  dt   = REAL(rdt,  mk)

  !----------------------------------------------------------------------------!
  ! deallocate tempoarary memory
  !----------------------------------------------------------------------------!
  DEALLOCATE(state)

END SUBROUTINE wvic_netcdf2field
