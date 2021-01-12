!-------------------------------------------------------------------------------
!* filename: wvic_io_dump_vector_netcdf                                       *!
!* project : wvic                                                             *!
!* purpose : dump a vector field to a netcdf file                             *!
!*         :                                                                  *!
!* author  : Philippe Chatelain                                               *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 18:17:05 2004                                         *!
!* please return to <pchatela@inf.ethz.ch>                                    *!
!
! $Log: wvic_io_dump_vector_netcdf.F,v $
! Revision 1.5  2006/09/16 00:22:03  pchatela
! Implemented the kinetic energy spectrum, dumped into an ascii file.
!
! Revision 1.4  2006/09/11 16:52:46  pchatela
! Fixes
!
! Revision 1.3  2006/09/11 14:57:27  pchatela
! Fixed velocity computation with odd/even symmetries
! Added flag for adaptive time step
!
! Revision 1.2  2006/09/01 17:33:29  pchatela
! small bugfix
!
! Revision 1.1  2006/09/01 15:45:19  pchatela
! IO module stuff
!
!
! 
!------------------------------------------------------------------------------!

#if __WHAT == __REAL
SUBROUTINE wvic_io_dump_vector_netcdf(vectfield, n_data, vectname, iter, info)
#elif __WHAT == __COMPLEX
SUBROUTINE wvic_io_dump_vectorc_netcdf(vectfield, n_data, vectname, iter, info)
#endif

  USE module_wvic
  USE ppm_module_write
  USE netcdf

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
  INTEGER            :: ncid, idnx, idny, idnz, idvsize
  INTEGER            :: idnaxes, idndeltas, idntsteps
  INTEGER            :: ivgrid, ivstate, ivtime, ivdt
  CHARACTER(len=256) :: msg
  INTEGER            :: stat, isub, isubl
  INTEGER            :: i,j,k,ivectname,chunksize
  INTEGER, PARAMETER :: naxes = 3, ndeltas = 2, ntsteps = 1, vsize = 3
  REAL(ms), DIMENSION(2,3) :: grid
  REAL(ms), DIMENSION(:,:,:,:), ALLOCATABLE :: state
  INTEGER            :: NETCDFPREC

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
  IF(mk.EQ.KIND(1.0e0)) THEN
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
  stat = nf90_def_dim(ncid, 'nx', n_data(1,isubl), idnx)
  stat = nf90_def_dim(ncid, 'ny', n_data(2,isubl), idny)
  stat = nf90_def_dim(ncid, 'nz', n_data(3,isubl), idnz)
  stat = nf90_def_dim(ncid, 'vsize', vsize, idvsize)
  stat = nf90_def_dim(ncid, 'naxes', naxes, idnaxes)
  stat = nf90_def_dim(ncid, 'ndeltas', ndeltas, idndeltas)
  stat = nf90_def_dim(ncid, 'ntsteps', ntsteps, idntsteps)
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! define variables
  stat = nf90_def_var(ncid, 'grid',  NETCDFPREC, &
       &             (/idndeltas, idnaxes/), ivgrid)
  stat = nf90_def_var(ncid, 'state', NETCDFPREC, &
       &             (/idvsize,idnx,idny,idnz/), ivstate)
  stat = nf90_def_var(ncid, 'time',  NETCDFPREC, &
       &             idntsteps, ivtime)
  stat = nf90_def_var(ncid, 'dt',    NETCDFPREC, &
       &             idntsteps, ivdt)
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! put attributes
  stat = nf90_put_att(ncid, ivstate, 'field', 'state, vector')
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
  stat = nf90_put_var(ncid=ncid,varid=ivtime, values=time)
  stat = nf90_put_var(ncid=ncid,varid=ivdt,   values=REAL(dt))
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! state vector data
  ALLOCATE(state(vsize,n_data(1,isubl),n_data(2,isubl),n_data(3,isubl)),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'scalar netcdf dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_dump_vector_field',msg,info)
     GOTO 9999
  END IF
  DO k=1,n_data(3,isubl)
     DO j=1,n_data(2,isubl)
        DO i=1,n_data(1,isubl)
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
  stat = nf90_put_var(ncid=ncid,varid=ivstate, values=state)
  DEALLOCATE(state,stat=stat)
  stat = nf90_close(ncid)
  !-----------------------------------------------------------------------------

#if __WHAT == __COMPLEX
  !-----------------------------------------------------------------------------
  ! construct filename
  !-----------------------------------------------------------------------------
  ivectname = LEN_TRIM(vectname)
  WRITE(ncfile,'(A,A,A,A,I5.5,A,I5.5,A)') &
       &vectname(1:ivectname),'_ip_',runtag(1:iruntag),'R',rank,'I',iter,'.nc'
  IF(mk.EQ.KIND(1.0e0)) THEN
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
  stat = nf90_def_dim(ncid, 'nx', n_data(1,isubl), idnx)
  stat = nf90_def_dim(ncid, 'ny', n_data(2,isubl), idny)
  stat = nf90_def_dim(ncid, 'nz', n_data(3,isubl), idnz)
  stat = nf90_def_dim(ncid, 'vsize', vsize, idvsize)
  stat = nf90_def_dim(ncid, 'naxes', naxes, idnaxes)
  stat = nf90_def_dim(ncid, 'ndeltas', ndeltas, idndeltas)
  stat = nf90_def_dim(ncid, 'ntsteps', ntsteps, idntsteps)
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! define variables
  stat = nf90_def_var(ncid, 'grid',  NETCDFPREC, &
       &             (/idndeltas, idnaxes/), ivgrid)
  stat = nf90_def_var(ncid, 'state', NETCDFPREC, &
       &             (/idvsize,idnx,idny,idnz/), ivstate)
  stat = nf90_def_var(ncid, 'time',  NETCDFPREC, &
       &             idntsteps, ivtime)
  stat = nf90_def_var(ncid, 'dt',    NETCDFPREC, &
       &             idntsteps, ivdt)
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! put attributes
  stat = nf90_put_att(ncid, ivstate, 'field', 'state, vector')
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
  stat = nf90_put_var(ncid=ncid,varid=ivtime, values=time)
  stat = nf90_put_var(ncid=ncid,varid=ivdt,   values=REAL(dt))
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! state vector data
  ALLOCATE(state(vsize,n_data(1,isubl),n_data(2,isubl),n_data(3,isubl)),stat=stat)
  IF(stat.NE.0) THEN
     WRITE(msg,*) 'vector netcdf dump failed, not enough memory'
     CALL ppm_write(rank,'wvic_dump_vector_field',msg,info)
     GOTO 9999
  END IF
  DO k=1,n_data(3,isubl)
     DO j=1,n_data(2,isubl)
        DO i=1,n_data(1,isubl)
           state(1,i,j,k) = REAL(AIMAG(vectfield(1,k,j,i,isub)),ms)
           state(2,i,j,k) = REAL(AIMAG(vectfield(2,k,j,i,isub)),ms)
           state(3,i,j,k) = REAL(AIMAG(vectfield(3,k,j,i,isub)),ms)
        END DO
     END DO
  END DO
  WRITE(*,*) 'Maxval(state): ',MAXVAL(REAL(state))
  stat = nf90_put_var(ncid=ncid,varid=ivstate, values=state)
  DEALLOCATE(state,stat=stat)
  stat = nf90_close(ncid)
  !-----------------------------------------------------------------------------
#endif

  IF (verbose) THEN
     WRITE(msg,*) 'netcdf data written'
     CALL ppm_write(rank,'wvic_io_dump_vector_netcdf',msg,info)
  END IF

9999 CONTINUE
  
#if __WHAT == __REAL
END SUBROUTINE wvic_io_dump_vector_netcdf
#elif __WHAT == __COMPLEX
END SUBROUTINE wvic_io_dump_vectorc_netcdf
#endif
