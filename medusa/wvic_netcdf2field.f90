!------------------------------------------------------------------------------!
!* filename: wvic_field2netcdf                                                 *!
!* project : ppm                                                              *!
!* purpose : write field into a dx compatible netcdf file                     *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 19:42:31 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_netcdf2field.f,v $
!  Revision 1.3  2006/08/24 11:29:23  menahel
!  ibncluded a to-real cast for time and dt
!
!  Revision 1.2  2006/08/24 09:47:46  menahel
!  changed rank size from 4.4 to 5.5 plus writing/reading dt now.
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.2  2005/11/21 17:36:37  michaebe
!  major session. now let''s debug
!
!  Revision 1.1  2005/09/28 11:40:34  michaebe
!  Fork from ppm_pvc
!
!  Revision 1.4  2005/01/06 12:28:42  michaebe
!  all the changes needed for the tubes
!
!  Revision 1.3  2004/12/04 17:34:17  michaebe
!  some mods, added the dump of time
!
!  Revision 1.2  2004/12/04 17:19:05  michaebe
!  some typos
!
!  Revision 1.1  2004/12/04 17:11:11  michaebe
!  second implementation, as dalco deleted the first....
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic FIELD2NETCDF =
! dump stuff to netcdf format
!------------------------------------------------------------------------------!
SUBROUTINE wvic_netcdf2field (info)
  
  USE module_wvic
  USE ppm_module_write
  USE netcdf
  
  INCLUDE 'mpif.h'
  
  !----------------------------------------------------------------------------!
  ! arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(inout)    :: info
  
  !----------------------------------------------------------------------------!
  ! localities
  !----------------------------------------------------------------------------!
  INTEGER, PARAMETER :: ms = kind(1.0e0)
  CHARACTER(len=256) :: ncfile
  INTEGER            :: ncid, idnx, idny, idnz, idvsize
  INTEGER            :: idnaxes, idndeltas, idntsteps
  INTEGER            :: ivgrid, ivstate, ivtime, varid, ivdt
  CHARACTER(len=256) :: msg
  INTEGER            :: stat, isub, isubl,iistart,iifinish,ibb
  INTEGER            :: i,j,k,initialsize,chunksize
  INTEGER, PARAMETER :: naxes = 3, ndeltas = 2, ntsteps = 1, vsize = 3
  REAL(ms), DIMENSION(2,3) :: grid
  REAL(ms), DIMENSION(:,:,:,:), ALLOCATABLE :: state
  REAL(ms)           :: rdt, rtime
  INTEGER            :: NETCDFPREC
  !----------------------------------------------------------------------------!
  ! construct filename
  !----------------------------------------------------------------------------!
  WRITE(ncfile,'(A,A,I5.5,A,I5.5,A)') &
       &runtag(1:iruntag),'R',rank,'I',netcdf_itime,'.nc'
  IF(ms.EQ.KIND(1.0e0)) THEN
     NETCDFPREC = NF90_FLOAT
  ELSE
     NETCDFPREC = NF90_DOUBLE
  END IF
  
  !----------------------------------------------------------------------------!
  ! open the data set
  !----------------------------------------------------------------------------!
  stat = nf90_open(ncfile, NF90_NOWRITE,ncid)

  !----------------------------------------------------------------------------!
  ! get the variable id
  !----------------------------------------------------------------------------!
  stat = nf90_inq_varid(ncid,"state",varid)

  !----------------------------------------------------------------------------!
  ! allocate temporary memory, ...
  !----------------------------------------------------------------------------!
  isub = 1; isubl = isublist(isub)
  ALLOCATE(state(vsize,ndata(1,isubl),ndata(2,isubl),ndata(3,isubl)),stat=stat)
  
  !----------------------------------------------------------------------------!
  ! and read the variable
  !----------------------------------------------------------------------------!
  stat = nf90_get_var(ncid,varid,state)

  !----------------------------------------------------------------------------!
  ! copy the contents to field_wp
  !----------------------------------------------------------------------------!
  DO k=1,ndata(3,isubl)
     DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
           field_wp(1,i,j,k,isub) = REAL(state(1,i,j,k),mk)
           field_wp(2,i,j,k,isub) = REAL(state(2,i,j,k),mk)
           field_wp(3,i,j,k,isub) = REAL(state(3,i,j,k),mk)
        END DO
     END DO
  END DO

  !----------------------------------------------------------------------------!
  ! get time and timestep info
  !----------------------------------------------------------------------------!
  stat = nf90_inq_varid(ncid,"time",varid)
  stat = nf90_get_var(ncid,varid,rtime)
  stat = nf90_inq_varid(ncid,"dt"  ,varid)
  IF(stat.EQ.0) stat = nf90_get_var(ncid,varid,rdt)
  time = REAL(rtime,mk)
  dt   = REAL(rdt,  mk)
  
  !----------------------------------------------------------------------------!
  ! deallocate tempoarary memory
  !----------------------------------------------------------------------------!
  DEALLOCATE(state)

  !----------------------------------------------------------------------------!
  ! close the netcdf file
  !----------------------------------------------------------------------------!
  stat = nf90_close(ncid)
  
END SUBROUTINE wvic_netcdf2field
