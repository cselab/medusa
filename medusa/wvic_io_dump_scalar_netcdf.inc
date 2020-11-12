!-------------------------------------------------------------------------------
!* filename: wvic_io_dump_scalar_netcdf                                       *!
!* project : ppm                                                              *!
!* purpose : dump a scalar field to a netcdf file                             *!
!*         :                                                                  *!
!* author  : Philippe Chatelain                                               *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 18:17:05 2004                                         *!
!* please return to <pchatela@inf.ethz.ch>                                    *!
!
! $Log: wvic_io_dump_scalar_netcdf.F,v $
! Revision 1.3  2006/09/11 16:52:46  pchatela
! Fixes
!
! Revision 1.2  2006/09/11 14:57:27  pchatela
! Fixed velocity computation with odd/even symmetries
! Added flag for adaptive time step
!
! Revision 1.1  2006/09/01 15:45:19  pchatela
! IO module stuff
!
!
! 
!------------------------------------------------------------------------------!


SUBROUTINE wvic_io_dump_scalar_netcdf(scalfield, n_data, scalname, iter,info)

  USE module_wvic
  USE ppm_module_write
  USE netcdf

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
  INTEGER            :: ncid, idnx, idny, idnz
  INTEGER            :: idnaxes, idndeltas, idntsteps
  INTEGER            :: ivgrid, ivstate, ivtime, ivdt
  CHARACTER(len=256) :: msg
  INTEGER            :: stat, isub, isubl
  INTEGER            :: i,j,k, iscalname,chunksize
  INTEGER, PARAMETER :: naxes = 3, ndeltas = 2, ntsteps = 1
  REAL(ms), DIMENSION(2,3) :: grid
  REAL(ms), DIMENSION(:,:,:), ALLOCATABLE :: state
  INTEGER            :: NETCDFPREC
  
  !-----------------------------------------------------------------------------
  ! construct filename
  !-----------------------------------------------------------------------------
  iscalname = LEN_TRIM(scalname)
  WRITE(ncfile,'(A,A,A,A,I5.5,A,I5.5,A)')&
       &scalname(1:iscalname),'_',runtag(1:iruntag),'R',rank,'I',iter,'.nc'
  IF(ms.EQ.KIND(1.0e0)) THEN
     NETCDFPREC = NF90_FLOAT
  ELSE
     NETCDFPREC = NF90_DOUBLE
  END IF
  !-----------------------------------------------------------------------------
  ! creeate the dataset
  !-----------------------------------------------------------------------------
  chunksize = 65536
  stat = nf90_create(path=ncfile,cmode=NF90_CLOBBER, ncid=ncid, & 
       &  chunksize=chunksize)
  !-----------------------------------------------------------------------------
  ! put the dataset into definition mode
  !-----------------------------------------------------------------------------
  stat = nf90_redef(ncid)

  !-----------------------------------------------------------------------------
  ! define the dimensions
  
  isub = 1; isubl = isublist(isub)
  stat = nf90_def_dim(ncid, 'nx', n_data(3,isubl), idnx)
  stat = nf90_def_dim(ncid, 'ny', n_data(2,isubl), idny)
  stat = nf90_def_dim(ncid, 'nz', n_data(1,isubl), idnz)
  stat = nf90_def_dim(ncid, 'naxes', naxes, idnaxes)
  stat = nf90_def_dim(ncid, 'ndeltas', ndeltas, idndeltas)
  stat = nf90_def_dim(ncid, 'ntsteps', ntsteps, idntsteps)
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! define variables
  stat = nf90_def_var(ncid, 'grid',  NETCDFPREC, &
       &             (/idndeltas, idnaxes/), ivgrid)
  stat = nf90_def_var(ncid, 'state', NETCDFPREC, &
       &             (/idnx,idny,idnz/), ivstate)
  stat = nf90_def_var(ncid, 'time',  NETCDFPREC, &
       &             idntsteps, ivtime)
  stat = nf90_def_var(ncid, 'dt',    NETCDFPREC, &
       &             idntsteps, ivdt)
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! put attributes
  stat = nf90_put_att(ncid, ivstate, 'field', 'state, scalar')
  stat = nf90_put_att(ncid, ivstate, 'positions', 'grid, regular')
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! leave define mode
  !-----------------------------------------------------------------------------
  stat = nf90_enddef(ncid)

  !-----------------------------------------------------------------------------
  ! now were going to write the data
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! grid data
  isub = 1; isubl = isublist(isub)
  grid(1,1) = REAL(min_sub(3,isubl),MS)
  grid(1,2) = REAL(min_sub(2,isubl),MS)
  grid(1,3) = REAL(min_sub(1,isubl),MS)
  grid(2,1) = REAL(dZ);grid(2,2) = REAL(dy); grid(2,3) = REAL(dx)
  stat = nf90_put_var(ncid=ncid, varid=ivgrid, values=grid)
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! lets remember what time it is ..
  stat = nf90_put_var(ncid=ncid,varid=ivtime, values=REAL(time))
  stat = nf90_put_var(ncid=ncid,varid=ivdt,   values=REAL(dt))
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! state vector data
  ALLOCATE(state(n_data(3,isubl),n_data(2,isubl),n_data(1,isubl)),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'scalar netcdf dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_dump_scalar_field',msg,info)
     GOTO 9999
  END IF
  DO k=1,n_data(3,isubl)
     DO j=1,n_data(2,isubl)
        DO i=1,n_data(1,isubl)
           state(k,j,i) = REAL(scalfield(i,j,k,isub))
        END DO
     END DO
  END DO
  stat = nf90_put_var(ncid=ncid,varid=ivstate, values=state)
  DEALLOCATE(state,stat=stat)
  
  !-----------------------------------------------------------------------------
  
  IF (verbose) THEN
      WRITE(msg,*) 'netcdf data written'
      CALL ppm_write(rank,'wvic_io_dump_scalar_netcdf',msg,info)
  END IF

9999 CONTINUE
  stat = nf90_close(ncid)
  RETURN
END SUBROUTINE wvic_io_dump_scalar_netcdf

