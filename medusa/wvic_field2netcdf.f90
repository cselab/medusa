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
!  $Log: wvic_field2netcdf.f,v $
!  Revision 1.3  2006/08/24 11:29:23  menahel
!  ibncluded a to-real cast for time and dt
!
!  Revision 1.2  2006/08/24 09:47:45  menahel
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
SUBROUTINE wvic_field2netcdf (info)
  
  USE module_wvic
  USE ppm_module_write
  USE netcdf
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
  INTEGER            :: ncid, idnx, idny, idnz, idvsize
  INTEGER            :: idnaxes, idndeltas, idntsteps
  INTEGER            :: ivgrid, ivstate, ivtime, ivdt
  CHARACTER(len=256) :: msg
  INTEGER            :: stat, isub, isubl,iistart,iifinish,ibb
  INTEGER            :: i,j,k,initialsize,chunksize
  INTEGER, PARAMETER :: naxes = 3, ndeltas = 2, ntsteps = 1, vsize = 3
  REAL(ms), DIMENSION(2,3) :: grid
  REAL(ms), DIMENSION(:,:,:,:), ALLOCATABLE :: state
  INTEGER            :: NETCDFPREC
  !----------------------------------------------------------------------------!
  ! construct filename
  !----------------------------------------------------------------------------!
  WRITE(ncfile,'(A,A,I5.5,A,I5.5,A)') &
       &runtag(1:iruntag),'R',rank,'I',itime,'.nc'
  IF(ms.EQ.KIND(1.0e0)) THEN
     NETCDFPREC = NF90_FLOAT
  ELSE
     NETCDFPREC = NF90_DOUBLE
  END IF
  ! IIFINISH=1
  ! DO IBB=1,2048,128
  !    IISTART=IIFINISH
  !    IIFINISH=IBB+128
  !    IF(((rank+1).GE.IISTART).AND.((rank+1).LT.IIFINISH)) THEN
  
  !----------------------------------------------------------------------------!
  ! creeate the dataset
  !----------------------------------------------------------------------------!
  !initialsize = ndata(1,isublist(1))*ndata(2,isublist(1))*ndata(3,isublist(1))*3*4 
  chunksize = 65536
  stat = nf90_create(path=ncfile,cmode=NF90_CLOBBER, ncid=ncid, & 
       &  chunksize=chunksize)
  !----------------------------------------------------------------------------!
  ! put the dataset into definition mode
  !----------------------------------------------------------------------------!
  stat = nf90_redef(ncid)
  
  !----------------------------------------------------------------------------!
  ! define the dimensions
  isub = 1; isubl = isublist(isub)
  stat = nf90_def_dim(ncid, 'nx', ndata(1,isubl), idnx)
  stat = nf90_def_dim(ncid, 'ny', ndata(2,isubl), idny)
  stat = nf90_def_dim(ncid, 'nz', ndata(3,isubl), idnz)
  stat = nf90_def_dim(ncid, 'vsize', vsize, idvsize)
  stat = nf90_def_dim(ncid, 'naxes', naxes, idnaxes) 
  stat = nf90_def_dim(ncid, 'ndeltas', ndeltas, idndeltas)
  stat = nf90_def_dim(ncid, 'ntsteps', ntsteps, idntsteps)
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  ! define variables
  stat = nf90_def_var(ncid, 'grid',  NETCDFPREC, &
       &             (/idndeltas, idnaxes/), ivgrid)
  stat = nf90_def_var(ncid, 'state', NETCDFPREC, &
       &             (/idvsize,idnx,idny,idnz/), ivstate)
  stat = nf90_def_var(ncid, 'time',  NETCDFPREC, &
       &             idntsteps, ivtime)
  stat = nf90_def_var(ncid, 'dt',    NETCDFPREC, &
       &             idntsteps, ivdt)
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! put attributes
  stat = nf90_put_att(ncid, ivstate, 'field', 'state, vector')
  stat = nf90_put_att(ncid, ivstate, 'positions', 'grid, regular')
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! leave define mode
  !----------------------------------------------------------------------------!
  stat = nf90_enddef(ncid)

  !----------------------------------------------------------------------------!
  ! now were going to write the data
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! grid data
  isub = 1; isubl = isublist(isub)
  grid(1,1) = REAL(min_sub(3,isubl),MS)
  grid(1,2) = REAL(min_sub(2,isubl),MS)
  grid(1,3) = REAL(min_sub(1,isubl),MS)
  grid(2,1) = REAL(dZ);grid(2,2) = REAL(dy); grid(2,3) = REAL(dx)
  stat = nf90_put_var(ncid=ncid, varid=ivgrid, values=grid)
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! lets remember what time it is ..
  stat = nf90_put_var(ncid=ncid,varid=ivtime, values=REAL(time))
  stat = nf90_put_var(ncid=ncid,varid=ivdt,   values=REAL(dt))
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! state vector data
  ALLOCATE(state(vsize,ndata(1,isubl),ndata(2,isubl),ndata(3,isubl)),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'netcdf field dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_field2netcdf',msg,info)
     GOTO 9999
  END IF
  DO k=1,ndata(3,isubl)
     DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
           state(1,i,j,k) = REAL(field_wp(1,i,j,k,isub))
           state(2,i,j,k) = REAL(field_wp(2,i,j,k,isub))
           state(3,i,j,k) = REAL(field_wp(3,i,j,k,isub))
        END DO
     END DO
  END DO
  stat = nf90_put_var(ncid=ncid,varid=ivstate, values=state)
  DEALLOCATE(state,stat=stat)
  
  !----------------------------------------------------------------------------!
  
  IF(verbose) THEN
     WRITE(msg,*) 'netcdf data written'
     CALL ppm_write(rank,'wvic_field2netcdf',msg,info)
  END IF
  
9999 CONTINUE
  stat = nf90_close(ncid)
! END IF
! CALL MPI_BARRIER(COMM,INFO)
! END DO
  WRITE(msg,*) 'complete'
!  CALL ppm_Write(rank,'wvic_field2netcdf',msg,info)

END SUBROUTINE wvic_field2netcdf
