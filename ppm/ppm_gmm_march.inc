      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_gmm_march
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs the group marching to
      !                 extend the infprmation from the interface to a
      !                 narrow band. This is used for level set
      !                 reinitialization or function extension. The
      !                 interface is assumed to be at the zero level.
      !                 Use postprocess shifting if this is not the
      !                 case. 
      !
      !  Input        : width           (F) Width of the narrow band to
      !                                     be produced on each side of
      !                                     the interface.
      !                 order           (I) Order of the method to be
      !                                     used. One of 
      !                                        ppm_param_order_1
      !                                        ppm_param_order_2
      !                                        ppm_param_order_3
      !                 rhscst          (F) constant value for the right
      !                                     hand side of grad u * grad f
      !                                     = c. If speed is present,
      !                                     this argument will be
      !                                     ignored.
      !                 MaxIter         (I) argument specifying the
      !                                     maximum number of allowed
      !                                     iterations. This can be useful
      !                                     since a cyclic dependency in the
      !                                     GMM algorithms could cause 
      !                                     infinite loops.
      !                 speed(:,:,:,:)  (F) rank 4 (3d) or rank 3 (2d)
      !                                     field of front speeds.
      !                                     OPTIONAL to override rhscst.
      !                 chi([:],:,:,:,:)(F) rank 5 (3d) or rank 4 (2d)
      !                                     field specifying the positions
      !                                     of the grid nodes. 1st index:
      !                                     1..ppm_dim, then i,j,[k],isub.
      !                                     OPTIONAL. Uniform grid is
      !                                     assumed if absent. Ghostlayers
      !                                     of size >=1 must be pre-filled.
      !
      !  Input/output : fdata(...)      (F) Level data. Either rank 3
      !                                     (for 2D scalar fields), or rank
      !                                     4 (for 3D scalar fields).
      !                                     Indices: (i,k,[k],isub).
      !                                     On input: signed distance
      !                                     function after
      !                                     initialization as returned
      !                                     by ppm_gmm_kickoff.
      !                                     Uninitialized values
      !                                     must be set to HUGE.
      !                                     On output: 
      !                                     approximation of the
      !                                     extended signed distance
      !                                     function in the whole band.
      !                                     Points far
      !                                     from the interface will have
      !                                     the value HUGE.
      !                 udata(...)      (F) Function data. OPTIONAL. If
      !                                     present, udata will be
      !                                     extended in the whole band
      !                                     such that its gradient is
      !                                     orthogonal to the gradient
      !                                     of fdata. fdata in this case
      !                                     is a pure input and needs to
      !                                     be available in the whole
      !                                     band. Needs to be a scalar
      !                                     field.
      !
      !  Output       : info            (I) return status. 0 on success.
      !
      !  Remarks      : This routine uses the o(N) group marching method
      !                 of Kim, rather than o(NlogN) fast marching. If
      !                 combines this with the high order equations by
      !                 Chopp to get 2nd or 3rd order accurate
      !                 extrapolations in linear time.
      !
      !                 speed should maybe be passed as speed**2. This
      !                 would save the multiplication in gmm_slvupwind,
      !                 but would require SQRTs in the initial
      !                 conversion from distance to travel time.
      !
      !                 When recomputing the close neighbors one could
      !                 detect local minima and directly remove (i.e.
      !                 accept) them to avoid unnecessary double
      !                 computations.
      !
      !                 To be more memory efficient, we could maybe
      !                 allocate gmm_state only for the sub we are
      !                 currently working on instead of for the full local
      !                 field. This will however only give cashback if
      !                 there are lots of subs per proc, which is currently
      !                 not the case for mesh topologies.
      !
      !  References   : Chopp:2001, Kim:2001b
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_march.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.9  2006/04/06 14:40:29  ivos
      !  Added the MaxIter argument to specify the maximum number of allowed
      !  iterations for the GMM marching.
      !
      !  Revision 1.8  2005/08/17 13:55:41  ivos
      !  fix: ghost layers are now only grown (if needed), but never shrunk.
      !  The latter confused the clients...
      !
      !  Revision 1.7  2005/07/14 19:58:13  ivos
      !  Added OPTIONAL argument chi for mesh node positions in distorted
      !  (mapped) meshes. For use with AGM for example.
      !
      !  Revision 1.6  2005/05/24 23:20:11  ivos
      !  bugfix in parallel version: ppm_gmm_state is now properly mapped
      !  and the pointers are restored at the end to avoid loss of data
      !  if the mapping routines reallocate a field.
      !
      !  Revision 1.5  2005/05/10 04:48:45  ivos
      !  Split marching and extension routines for faster compilation,
      !  Sharked extension routines, moved all initialization to gmm_init, and
      !  code cosmetics.
      !
      !  Revision 1.4  2005/04/27 01:06:11  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.3  2005/03/16 06:20:08  ivos
      !  Several bugfixes. 1st order version is now tested. Moved all large
      !  data to the module.
      !
      !  Revision 1.2  2005/03/11 04:17:01  ivos
      !  Added possiblity to for velocity extension and to get back
      !  the closest point transform.
      !
      !  Revision 1.1  2005/03/10 01:37:20  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_march_2ds(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_2dd(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_march_3ds(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_3dd(width,order,fdata,rhscst,MaxIter,info, &
     &    speed,udata,chi)
#endif 
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_data_gmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_map_field_ghost
      USE ppm_module_gmm_march_fwd
      USE ppm_module_gmm_march_bkwd
      USE ppm_module_gmm_extend_fwd
      USE ppm_module_gmm_extend_bkwd
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      USE MPI
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , INTENT(IN), OPTIONAL :: speed
      REAL(MK), DIMENSION(:,:,:)     , POINTER              :: fdata
      REAL(MK), DIMENSION(:,:,:)     , POINTER, OPTIONAL    :: udata
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: chi
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: speed
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER              :: fdata
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER, OPTIONAL    :: udata
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN), OPTIONAL :: chi
#endif
      REAL(MK)                       , INTENT(IN   )        :: width,rhscst
      INTEGER                        , INTENT(IN   )        :: order,MaxIter
      INTEGER                        , INTENT(  OUT)        :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                          :: i,j,k,p,iopt,xhi,yhi,zhi
      INTEGER                          :: npos,npos0,jsub,isub
      INTEGER                          :: nord
      INTEGER                          :: Marchit
      INTEGER, DIMENSION(3)            :: ghostsize
      INTEGER, DIMENSION(4)            :: ldl,ldu
      REAL(MK),DIMENSION(2)            :: coefs
      REAL(MK)                         :: t0,deltaT,sqrtdiminv,smin
      REAL(MK)                         :: TM,dx,dy,dz,dxinv,dyinv,dzinv
      REAL(MK)                         :: hsave,mindx
      REAL(MK)                         :: ainv,big
#ifdef __MPI
      REAL(MK)                         :: allsmin
#endif
      LOGICAL                          :: lok
      CHARACTER(LEN=ppm_char)          :: cbuf
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:    ), POINTER :: fdta
      REAL(MK), DIMENSION(:,:,:    ), POINTER :: dta
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:  ), POINTER :: fdta
      REAL(MK), DIMENSION(:,:,:,:  ), POINTER :: dta
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_march',t0,info)
      big        = HUGE(big)
      hsave      = 0.9_MK*big
#if   __DIM == __2D
      sqrtdiminv = 1.0_MK/SQRT(2.0_MK)
#else
      sqrtdiminv = 1.0_MK/SQRT(3.0_MK)
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (width .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march',  &
     &            'width must be positive!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((order.NE.ppm_param_order_1).AND.(order.NE.ppm_param_order_2)  &
     &         .AND.(order.NE.ppm_param_order_3)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march',  &
     &            'order must be 1, 2, or 3!',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __DIM == __3D
          IF (SIZE(fdata,4) .LT. ppm_nsublist(gmm_topoid)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,3) .LT. maxzhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'z dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == __2D
          IF (SIZE(fdata,3) .LT. ppm_nsublist(gmm_topoid)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_march',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF          ! ppm_debug for argument check

      !-------------------------------------------------------------------------
      !  Determine ghostsize needed. Only change if larger than what it
      !  already is. Do not shrink as this would freak out the client...
      !-------------------------------------------------------------------------
#if   __DIM == __3D
      i = MIN(1-LBOUND(fdata,1),1-LBOUND(fdata,2),1-LBOUND(fdata,3))
#else
      i = MIN(1-LBOUND(fdata,1),1-LBOUND(fdata,2))
#endif
      IF (order .EQ. ppm_param_order_1) THEN
          ghostsize = MAX(1,i)
          nord = 1
       !  IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A)') 'Using first order upwinding'
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
       !  ENDIF
      ELSEIF (order .EQ. ppm_param_order_2) THEN
          ghostsize = MAX(2,i)
          nord = 2
       !  IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A)') 'Using second order upwinding'
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
       !  ENDIF
      ELSE
          ghostsize = MAX(3,i)
          nord = 3
       !  IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A)') 'Using third order upwinding'
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
       !  ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Set the pointers
      !-------------------------------------------------------------------------
      IF (PRESENT(udata)) THEN
          !---------------------------------------------------------------------
          !  March udata and use fdata as level set 
          !---------------------------------------------------------------------
          dta => udata
          fdta => fdata
      ELSE
          !---------------------------------------------------------------------
          !  March the level set itself
          !---------------------------------------------------------------------
          dta => fdata
          fdta => fdata
      ENDIF

      !-------------------------------------------------------------------------
      !  Find mesh spacing
      !-------------------------------------------------------------------------
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          dx = (ppm_max_physs(1,gmm_topoid)-ppm_min_physs(1,gmm_topoid))/   &
     &        REAL(ppm_cart_mesh(gmm_meshid,gmm_topoid)%Nm(1)-1,ppm_kind_single)
          dy = (ppm_max_physs(2,gmm_topoid)-ppm_min_physs(2,gmm_topoid))/   &
     &        REAL(ppm_cart_mesh(gmm_meshid,gmm_topoid)%Nm(2)-1,ppm_kind_single)
          IF (ppm_dim .GT. 2) THEN
              dz = (ppm_max_physs(3,gmm_topoid)-ppm_min_physs(3,gmm_topoid))/ &
     &            REAL(ppm_cart_mesh(gmm_meshid,gmm_topoid)%Nm(3)-1,     &
     &            ppm_kind_single)
          ENDIF
      ELSE
          dx = (ppm_max_physd(1,gmm_topoid)-ppm_min_physd(1,gmm_topoid))/   &
     &        REAL(ppm_cart_mesh(gmm_meshid,gmm_topoid)%Nm(1)-1,ppm_kind_double)
          dy = (ppm_max_physd(2,gmm_topoid)-ppm_min_physd(2,gmm_topoid))/   &
     &        REAL(ppm_cart_mesh(gmm_meshid,gmm_topoid)%Nm(2)-1,ppm_kind_double)
          IF (ppm_dim .GT. 2) THEN
              dz = (ppm_max_physd(3,gmm_topoid)-ppm_min_physd(3,gmm_topoid))/ &
     &            REAL(ppm_cart_mesh(gmm_meshid,gmm_topoid)%Nm(3)-1,     &
     &            ppm_kind_double)
          ENDIF
      ENDIF
      dxinv = 1.0_MK/dx
      dyinv = 1.0_MK/dy
#if   __DIM == __3D
      dzinv = 1.0_MK/dz
      mindx = MIN(dx,dy,dz)
#else
      dzinv = 0.0_MK
      mindx = MIN(dx,dy)
#endif

      !-------------------------------------------------------------------------
      !  Initialize the ghost layers.
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_init,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'ghost mapping init failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_ghost_get,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'ghost get mapping failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_push,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing level data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (PRESENT(udata)) THEN
          CALL ppm_map_field_ghost(dta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_push,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing field data failed',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_send,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (PRESENT(udata)) THEN
          CALL ppm_map_field_ghost(dta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_pop,info)
          IF (info .NE. ppm_param_success) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping field data failed',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_pop,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping level data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      ! because it has not been mapped
      IF (.NOT.(PRESENT(udata))) dta => fdta

      !-------------------------------------------------------------------------
      !  Allocate and initialize status array
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow
#if   __DIM == __3D
      ldl(1) = 1-ghostsize(1)
      ldl(2) = 1-ghostsize(2)
      ldl(3) = 1-ghostsize(3)
      ldl(4) = 1
      ldu(1) = maxxhi+ghostsize(1)
      ldu(2) = maxyhi+ghostsize(2)
      ldu(3) = maxzhi+ghostsize(3)
      ldu(4) = ppm_nsublist(gmm_topoid)
      CALL ppm_alloc(gmm_state3d,ldl,ldu,iopt,info)
#else
      ldl(1) = 1-ghostsize(1)
      ldl(2) = 1-ghostsize(2)
      ldl(3) = 1
      ldu(1) = maxxhi+ghostsize(1)
      ldu(2) = maxyhi+ghostsize(2)
      ldu(3) = ppm_nsublist(gmm_topoid)
      CALL ppm_alloc(gmm_state2d,ldl,ldu,iopt,info)
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_march',      &
     &        'point status array GMM_STATE',__LINE__,info)
          GOTO 9999
      ENDIF

#if   __DIM == __3D
      gmm_state3d = ppm_gmm_param_far

      !-------------------------------------------------------------------------
      !  Convert distance to travel time and find minimum slowness (=
      !  max speed)
      !-------------------------------------------------------------------------
      DO jsub=1,ppm_nsublist(gmm_topoid)
          IF (PRESENT(speed)) THEN
              isub = ppm_isublist(jsub,gmm_topoid)
              xhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub)
              yhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub)
              zhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(3,isub)
              smin  = -big
              DO k=1,zhi
                  DO j=1,yhi
                      DO i=1,xhi
                          IF (speed(i,j,k,jsub) .GT. smin)    &
     &                        smin = speed(i,j,k,jsub)
                          IF (dta(i,j,k,jsub) .LT. hsave) THEN
                              dta(i,j,k,jsub)=dta(i,j,k,jsub)/speed(i,j,k,jsub)
                          ENDIF
                      ENDDO
                  ENDDO
              ENDDO
              smin = 1.0_MK/smin
#ifdef __MPI
              !-----------------------------------------------------------------
              !  Ensure consistent computations on all processors
              !-----------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
              CALL MPI_Allreduce(smin,allsmin,1,MPI_REAL,MPI_MIN,   &
     &            ppm_comm,info)
#else
              CALL MPI_Allreduce(smin,allsmin,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
     &            ppm_comm,info)
#endif
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_gmm_march',      &
     &                'MPI_ALLREDUCE of smin',__LINE__,info)
                  GOTO 9999
              ENDIF
              smin = allsmin
#endif
          ELSE
              smin = 1.0_MK
          ENDIF
      ENDDO        ! jsub

      !-------------------------------------------------------------------------
      !  Initialize travel limit (see Kim:2001a)
      !-------------------------------------------------------------------------
      deltaT = 0.99_MK*sqrtdiminv*mindx*smin
      !IF (.NOT.PRESENT(udata)) deltaT = (10.0_MK**REAL(-nord+1,MK))*deltaT

      !-------------------------------------------------------------------------
      !  Make a list of all initialized points and find the minimum
      !  travel time on the boundary of the band.
      !-------------------------------------------------------------------------
      npos0 = 0
      DO jsub=1,ppm_nsublist(gmm_topoid)
          isub = ppm_isublist(jsub,gmm_topoid)
          xhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub)
          yhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub)
          zhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(3,isub)
          TM   = big
          DO k=1,zhi
              DO j=1,yhi
                  DO i=1,xhi
                      coefs(1) = dta(i,j,k,jsub)
                      coefs(2) = fdta(i,j,k,jsub)
                      IF ((coefs(1).LT.hsave).AND.(gmm_state3d(i,j,k,jsub).EQ.&
     &                    ppm_gmm_param_far)) THEN
                          IF ((dta(i-1,j,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i+1,j,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j-1,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j+1,k,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j,k-1,jsub) .GT. hsave) .OR.   &
     &                        (dta(i,j,k+1,jsub) .GT. hsave)) THEN
                              !-------------------------------------------------
                              !  Add points on the surface to close set and
                              !  set TM to min distance to interface
                              !-------------------------------------------------
                              IF (ABS(coefs(2)) .LT. TM) TM = ABS(coefs(2))
                              gmm_state3d(i,j,k,jsub) = ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                          ELSE
                              gmm_state3d(i,j,k,jsub) = ppm_gmm_param_accepted
                          ENDIF
                      ENDIF
                  ENDDO       ! mesh cell i
              ENDDO       ! mesh cell j
          ENDDO       ! mesh cell k
      ENDDO       ! jsub
      npos = npos0

      !-------------------------------------------------------------------------
      !  Initialize ghost layers for state array
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost(gmm_state3d,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_push,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(gmm_state3d,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_send,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending status ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(gmm_state3d,gmm_topoid,gmm_meshid,ghostsize,&
     &                         ppm_param_map_pop,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping status data failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  March 
      !-------------------------------------------------------------------------
      Marchit = 0
      DO WHILE ((npos .GT. 0) .AND. (Marchit .LT. MaxIter))
          Marchit = Marchit + 1
          TM = TM + deltaT
          IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A,I4,A,I9)') 'Iteration ',Marchit,   &
     &            ': Number of points in list: ',npos
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Hard-core debuging
          !---------------------------------------------------------------------
!         WRITE(cbuf,'(A,I4.4,A)') 'state_',Marchit,'.out'
!         OPEN(50,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
!         WRITE(cbuf,'(A,I4.4,A)') 'level_',Marchit,'.out'
!         OPEN(40,FILE=cbuf,STATUS='REPLACE',ACTION='WRITE')
!         DO k=1,zhi
!             DO j=1,yhi
!                 DO i=1,xhi
!                     IF (ABS(fdta(i,j,k,1)) .LT. hsave) THEN
!                         WRITE(40,*) fdta(i,j,k,1)
!                     ELSE
!                         WRITE(40,*) 0.0_MK
!                     ENDIF
!                     WRITE(50,*) gmm_state3d(i,j,k,1)
!                 ENDDO
!             ENDDO
!         ENDDO
!         CLOSE(40)
!         CLOSE(50)

          !---------------------------------------------------------------------
          !  Reverse order: recompute neighbors of points in ggm_ipos to
          !  fix stability.
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              DO i=1,nord
                  IF (PRESENT(speed)) THEN
                      CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ELSE
                      CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &                    'Backward extension failed.',__LINE__,info)
                      GOTO 8888
                  ENDIF
              ENDDO
          ELSE
              IF (PRESENT(speed)) THEN
                  CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
              ELSE
                  CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info)
              ENDIF
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &                'Backward marching failed.',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDIF

          !---------------------------------------------------------------------
          !  Forward order: update neighbors of points in ggm_ipos and
          !  advance front
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              IF (PRESENT(speed)) THEN
                  CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
              ELSE
                  CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info)
              ENDIF
          ELSE
              IF (PRESENT(speed)) THEN
                  CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
              ELSE
                  CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                rhscst,dxinv,dyinv,dzinv,ghostsize,info)
              ENDIF
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &            'Forward marching failed.',__LINE__,info)
              GOTO 8888
          ENDIF
      ENDDO          ! while npos.GT.0

#elif __DIM == __2D
      gmm_state2d = ppm_gmm_param_far

      !-------------------------------------------------------------------------
      !  Convert distance to travel time and find minimum slowness (=
      !  max speed)
      !-------------------------------------------------------------------------
      DO jsub=1,ppm_nsublist(gmm_topoid)
          IF (PRESENT(speed)) THEN
              isub = ppm_isublist(jsub,gmm_topoid)
              xhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub)
              yhi = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub)
              smin  = -big
              DO j=1,yhi
                  DO i=1,xhi
                      IF (speed(i,j,jsub) .GT. smin)    &
     &                    smin = speed(i,j,jsub)
                      IF (dta(i,j,jsub) .LT. hsave) THEN
                          dta(i,j,jsub) =               &
     &                        dta(i,j,jsub)/speed(i,j,jsub)
                      ENDIF
                  ENDDO
              ENDDO
              smin = 1.0_MK/smin
#ifdef __MPI
              !-----------------------------------------------------------------
              !  Ensure consistent computations on all processors
              !-----------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
              CALL MPI_Allreduce(smin,allsmin,1,MPI_REAL,MPI_MIN,   &
     &            ppm_comm,info)
#else
              CALL MPI_Allreduce(smin,allsmin,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
     &            ppm_comm,info)
#endif
              IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_mpi_fail,'ppm_gmm_march',      &
     &                'MPI_ALLREDUCE of smin',__LINE__,info)
                  GOTO 9999
              ENDIF
              smin = allsmin
#endif
          ELSE
              smin = 1.0_MK
          ENDIF
      ENDDO        ! jsub

      !-------------------------------------------------------------------------
      !  Initialize travel limit (see Kim:2001a)
      !-------------------------------------------------------------------------
      deltaT = 0.99_MK*sqrtdiminv*mindx*smin
      !IF (.NOT.PRESENT(udata)) deltaT = (10.0_MK**REAL(-nord+1,MK))*deltaT

      !-------------------------------------------------------------------------
      !  Make a list of all initialized points and find the minimum
      !  travel time on the boundary of the band.
      !-------------------------------------------------------------------------
      npos0 = 0
      DO jsub=1,ppm_nsublist(gmm_topoid)
          isub = ppm_isublist(jsub,gmm_topoid)
          xhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub)
          yhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub)
          TM   = big
          DO j=1,yhi
              DO i=1,xhi
                  coefs(1) = dta(i,j,jsub)
                  coefs(2) = fdta(i,j,jsub)
                  IF ((coefs(1).LT.hsave).AND.(gmm_state2d(i,j,jsub).EQ.&
     &                ppm_gmm_param_far)) THEN
                      IF ((dta(i-1,j,jsub) .GT. hsave) .OR.   &
     &                    (dta(i+1,j,jsub) .GT. hsave) .OR.   &
     &                    (dta(i,j-1,jsub) .GT. hsave) .OR.   &
     &                    (dta(i,j+1,jsub) .GT. hsave)) THEN
                          !-----------------------------------------------------
                          !  Add points on the surface to close set and
                          !  set TM to min distance to interface
                          !-----------------------------------------------------
                          IF (ABS(coefs(2)) .LT. TM) TM = ABS(coefs(2))
                          gmm_state2d(i,j,jsub) = ppm_gmm_param_close
#include "ppm_gmm_add_to_list.inc"
                      ELSE
                          gmm_state2d(i,j,jsub) = ppm_gmm_param_accepted
                      ENDIF
                  ENDIF
              ENDDO       ! mesh cell i
          ENDDO       ! mesh cell j
      ENDDO       ! jsub
      npos = npos0

      !-------------------------------------------------------------------------
      !  Initialize ghost layers for state array
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost(gmm_state2d,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_push,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing status data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(gmm_state2d,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_send,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending status ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(gmm_state2d,gmm_topoid,gmm_meshid,ghostsize,&
     &                         ppm_param_map_pop,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping status data failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  March 
      !-------------------------------------------------------------------------
      Marchit = 0
      DO WHILE ((npos .GT. 0) .AND. (Marchit .LT. MaxIter))
          Marchit = Marchit + 1
          TM = TM + deltaT
          IF (ppm_debug .GT. 0) THEN
              WRITE(cbuf,'(A,I4,A,I9)') 'Iteration ',Marchit,   &
     &            ': Number of points in list: ',npos
              CALL ppm_write(ppm_rank,'ppm_gmm_march',cbuf,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Reverse order: recompute neighbors of points in ggm_ipos to
          !  fix stability
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              DO i=1,nord
                  IF (PRESENT(speed)) THEN
                      IF (PRESENT(chi)) THEN
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,&
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info,    &
     &                        speed,chi)
                      ELSE
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos,&
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                      ENDIF
                  ELSE
                      IF (PRESENT(chi)) THEN
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos, &
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info,     &
     &                        chi=chi)
                      ELSE
                          CALL ppm_gmm_extend_bkwd(fdta,dta,width,order,npos, &
     &                        TM,rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                      ENDIF
                  ENDIF
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &                    'Backward extension failed.',__LINE__,info)
                      GOTO 8888
                  ENDIF
              ENDDO
          ELSE
              IF (PRESENT(speed)) THEN
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
                  ELSE
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ENDIF
              ELSE
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,chi=chi)
                  ELSE
                      CALL ppm_gmm_march_bkwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
              ENDIF
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &            'Backward marching failed.',__LINE__,info)
              GOTO 8888
          ENDIF

          !---------------------------------------------------------------------
          !  Forward order: update neighbors of points in ggm_ipos and
          !  advance front
          !---------------------------------------------------------------------
          IF (PRESENT(udata)) THEN
              IF (PRESENT(speed)) THEN
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
                  ELSE
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ENDIF
              ELSE
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,chi=chi)
                  ELSE
                      CALL ppm_gmm_extend_fwd(fdta,dta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
              ENDIF
          ELSE
              IF (PRESENT(speed)) THEN
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
                  ELSE
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed)
                  ENDIF
              ELSE
                  IF (PRESENT(chi)) THEN
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info,chi=chi)
                  ELSE
                      CALL ppm_gmm_march_fwd(fdta,width,order,npos,TM,  &
     &                    rhscst,dxinv,dyinv,dzinv,ghostsize,info)
                  ENDIF
              ENDIF
          ENDIF
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',      &
     &            'Forward marching failed.',__LINE__,info)
              GOTO 8888
          ENDIF
      ENDDO          ! while npos.GT.0
#endif 

      !-------------------------------------------------------------------------
      !  Deallocate status array
      !-------------------------------------------------------------------------
 8888 iopt = ppm_param_dealloc
#if   __DIM == __3D
      CALL ppm_alloc(gmm_state3d,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_march',      &
     &        'point status array GMM_STATE3D',__LINE__,info)
      ENDIF
#elif __DIM == __2D
      CALL ppm_alloc(gmm_state2d,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_march',      &
     &        'point status array GMM_STATE2D',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  In case one of the mappings reallocated the data, we might want to
      !  pass back out the new data and save it from destruction upon
      !  subroutine termination.
      !-------------------------------------------------------------------------
      fdata => fdta
      IF (PRESENT(udata)) udata => dta
      NULLIFY(dta)
      NULLIFY(fdta)

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_march',t0,info)
      RETURN
#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_2dd
#endif 

#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_3dd
#endif 
#endif
