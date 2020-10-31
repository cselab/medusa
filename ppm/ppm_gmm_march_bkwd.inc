      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_gmm_march_bkwd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs the backward marching 
      !                 step of the GMM. See ppm_gmm_march for details.
      !
      !  Input        : width           (F) Width of the narrow band to
      !                                     be produced on each side of
      !                                     the interface.
      !                 order           (I) Order of the method to be
      !                                     used. One of 
      !                                        ppm_param_order_1
      !                                        ppm_param_order_2
      !                                        ppm_param_order_3
      !                 npos            (I) Current number of points in the 
      !                                     close set.
      !                 TM              (F) Current threshold for wave
      !                                     front location.
      !                 rhscst          (F) constant value for the right
      !                                     hand side of grad u * grad f
      !                                     = c. If speed is present,
      !                                     this argument will be
      !                                     ignored.
      !                 dxinv           (F) inverse of the x grid spacing.
      !                 dyinv           (F) inverse of the y grid spacing.
      !                 dzinv           (F) inverse of the z grid spacing
      !                                     (Not used in 2D version).
      !                 ghostsize(3)    (I) Size of the ghostlayer on all
      !                                     sides. 
      !                 speed(:,:,:,:)  (F) rank 4 (3d) or rank 3 (2d)
      !                                     field of front speeds.
      !                                     OPTIONAL to override rhscst.
      !                 chi([:],:,:,:,:)(F) rank 5 (3d) or rank 4 (2d)
      !                                     field specifying the positions
      !                                     of the grid nodes. 1st index:
      !                                     1..ppm_dim, then i,j,[k],isub.
      !                                     OPTIONAL. Uniform grid is
      !                                     assumed if absent.
      !
      !  Input/output : fdta(:,:,:,[:]) (F) pointer to level function.
      !
      !  Output       : info            (I) return status. 0 on success.
      !
      !  Remarks      : 
      !
      !  References   : Chopp:2001, Kim:2001b
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_march_bkwd.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.6  2005/07/14 19:58:14  ivos
      !  Added OPTIONAL argument chi for mesh node positions in distorted
      !  (mapped) meshes. For use with AGM for example.
      !
      !  Revision 1.5  2005/06/17 17:34:57  ivos
      !  Removed unused local variables.
      !
      !  Revision 1.4  2005/05/29 02:32:33  ivos
      !  bugfix: negative points were neglected due to incomplete check.
      !
      !  Revision 1.3  2005/05/22 03:41:26  ivos
      !  bugfix: values are now only replaced if they are of the same sign.
      !  The condition from Kim that the ABS is smaller is not sufficient
      !  for concave level sets.
      !
      !  Revision 1.2  2005/05/10 04:48:46  ivos
      !  Split marching and extension routines for faster compilation,
      !  Sharked extension routines, moved all initialization to gmm_init, and
      !  code cosmetics.
      !
      !  Revision 1.1  2005/04/27 01:08:38  ivos
      !  Initial commit, but tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_march_bkwd_2ds(fdta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_bkwd_2dd(fdta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#endif 
#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_march_bkwd_3ds(fdta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_march_bkwd_3dd(fdta,width,order,npos,TM,   &
     &    rhscst,dxinv,dyinv,dzinv,ghostsize,info,speed,chi)
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
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , POINTER          :: fdta
      REAL(MK), DIMENSION(:,:,:)     , INTENT(IN), OPTIONAL :: speed
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: chi
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER          :: fdta
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL :: speed
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN), OPTIONAL :: chi
#endif
      REAL(MK)                       , INTENT(IN   )    :: width,rhscst
      REAL(MK)                       , INTENT(IN   )    :: TM,dxinv,dyinv
      REAL(MK)                       , INTENT(IN   )    :: dzinv
      INTEGER                        , INTENT(IN   )    :: order
      INTEGER, DIMENSION(3)          , INTENT(IN   )    :: ghostsize
      INTEGER                        , INTENT(INOUT)    :: npos
      INTEGER                        , INTENT(  OUT)    :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                          :: i,j,k,p,xhi,yhi,zhi,ii,jj,kk
      INTEGER                          :: jsub,isub
      INTEGER                          :: i1,i2,i3
      INTEGER, DIMENSION(-3:3)         :: sx,sy,sz
      REAL(MK)                         :: t0,onethird,onetwelfth
      REAL(MK)                         :: dxihalf,dyihalf,dzihalf
      REAL(MK)                         :: dxitwelve,dyitwelve,dzitwelve
      REAL(MK)                         :: valijk,det,hsave,fdta0
      REAL(MK)                         :: lmyeps,ainv,big,absfdta0
      REAL(MK), DIMENSION(3)           :: coefs,gphi
      REAL(MK), DIMENSION(3,3)         :: jac,ji
      REAL(MK), DIMENSION(-3:3,ppm_dim):: phi,psi
      REAL(MK), DIMENSION(ppm_dim)     :: alpha,beta
      REAL(MK), DIMENSION(2)           :: roots
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_march_bkwd',t0,info)
      phi      = 0.0_MK
      psi      = 0.0_MK
      big      = HUGE(big)
      hsave    = 0.9_MK*big
      onethird = 1.0_MK/3.0_MK
#if   __KIND == __SINGLE_PRECISION
      lmyeps   = ppm_myepss
#else
      lmyeps   = ppm_myepsd
#endif
      dxihalf  = 0.5_MK*dxinv
      dyihalf  = 0.5_MK*dyinv
      dxitwelve  = onetwelfth*dxinv
      dyitwelve  = onetwelfth*dyinv
#if   __DIM == __3D
      dzihalf  = 0.5_MK*dzinv
      dzitwelve  = onetwelfth*dzinv
#endif 

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (width .LT. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march_bkwd',  &
     &            'width must be positive!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF ((order.NE.ppm_param_order_1).AND.(order.NE.ppm_param_order_2)  &
     &         .AND.(order.NE.ppm_param_order_3)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_march_bkwd',  &
     &            'order must be 1, 2, or 3!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF          ! ppm_debug for argument check

#if   __DIM == __3D
      !-------------------------------------------------------------------------
      !  Reverse order: recompute neighbors of points in ggm_ipos
      !-------------------------------------------------------------------------
      DO p=npos,1,-1
          ii   = gmm_ipos(1,p)
          jj   = gmm_ipos(2,p)
          kk   = gmm_ipos(3,p)
          jsub = gmm_ipos(4,p)
          isub = ppm_isublist(jsub,gmm_topoid)
          xhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub)
          yhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub)
          zhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(3,isub)
          fdta0= fdta(ii,jj,kk,jsub)
          absfdta0 = fdta0
          IF (absfdta0 .LT. 0.0_MK) absfdta0 = -absfdta0
          !---------------------------------------------------------------------
          !  GMM update condition (see Kim:2001a)
          !---------------------------------------------------------------------
          IF (.NOT.(absfdta0.GT.TM)) THEN
              !-----------------------------------------------------------------
              !  Recompute non-accepted close neighbors
              !-----------------------------------------------------------------
              i = ii - 1
              j = jj
              k = kk
              IF (i.GT.0) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                          IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,k,jsub).GT.hsave)) THEN 
                              fdta(i,j,k,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii + 1
              j = jj
              k = kk
              IF (i.LE.xhi) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                          IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,k,jsub).GT.hsave)) THEN 
                              fdta(i,j,k,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii
              j = jj - 1 
              k = kk
              IF (j.GT.0) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                          IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,k,jsub).GT.hsave)) THEN 
                              fdta(i,j,k,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii
              j = jj + 1 
              k = kk
              IF (j.LE.yhi) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                          IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,k,jsub).GT.hsave)) THEN 
                              fdta(i,j,k,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii
              j = jj 
              k = kk - 1
              IF (k.GT.0) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                          IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,k,jsub).GT.hsave)) THEN 
                              fdta(i,j,k,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii
              j = jj 
              k = kk + 1
              IF (k.LE.zhi) THEN
                  IF ((gmm_state3d(i,j,k,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,k,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j,k
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,k,jsub))) THEN
                          IF ((valijk*fdta(i,j,k,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,k,jsub).GT.hsave)) THEN 
                              fdta(i,j,k,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF           ! TT .LE. TM
      ENDDO         ! p=npos,1,-1

      !-------------------------------------------------------------------------
      !  Update ghost layers for fdta
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_push,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_send,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                     ppm_param_map_pop,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping field data failed',__LINE__,info)
          GOTO 9999
      ENDIF

#elif __DIM == __2D
      !-------------------------------------------------------------------------
      !  Reverse order: recompute neighbors of points in ggm_ipos
      !-------------------------------------------------------------------------
      DO p=npos,1,-1
          ii   = gmm_ipos(1,p)
          jj   = gmm_ipos(2,p)
          jsub = gmm_ipos(3,p)
          isub = ppm_isublist(jsub,gmm_topoid)
          xhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(1,isub)
          yhi  = ppm_cart_mesh(gmm_meshid,gmm_topoid)%nnodes(2,isub)
          fdta0= fdta(ii,jj,jsub)
          absfdta0 = fdta0
          IF (absfdta0 .LT. 0.0_MK) absfdta0 = -absfdta0
          !---------------------------------------------------------------------
          !  GMM update condition (see Kim:2001a)
          !---------------------------------------------------------------------
          IF (.NOT.(absfdta0.GT.TM)) THEN
              !-----------------------------------------------------------------
              !  Recompute non-accepted close neighbors
              !-----------------------------------------------------------------
              i = ii - 1
              j = jj
              IF (i.GT.0) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                          IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,jsub).GT.hsave)) THEN 
                              fdta(i,j,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii + 1
              j = jj
              IF (i.LE.xhi) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                          IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,jsub).GT.hsave)) THEN 
                              fdta(i,j,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii
              j = jj - 1 
              IF (j.GT.0) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                          IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,jsub).GT.hsave)) THEN 
                              fdta(i,j,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
              i = ii
              j = jj + 1 
              IF (j.LE.yhi) THEN
                  IF ((gmm_state2d(i,j,jsub) .NE.      &
     &                ppm_gmm_param_accepted) .AND.      &
     &                (ABS(fdta(i,j,jsub)).GT.absfdta0)) THEN
                      !---------------------------------------------------------
                      !  Update point i,j
                      !---------------------------------------------------------
#include "ppm_gmm_slvupwnd.inc"
                      IF (ABS(valijk).LT.ABS(fdta(i,j,jsub))) THEN
                          IF ((valijk*fdta(i,j,jsub).GE.0.0_MK) .OR.   &
     &                        (fdta(i,j,jsub).GT.hsave)) THEN 
                              fdta(i,j,jsub) = valijk 
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF           ! TT .LE. TM
      ENDDO         ! p=npos,1,-1

      !-------------------------------------------------------------------------
      !  Update ghost layers for fdta
      !-------------------------------------------------------------------------
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_push,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'pushing field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                         ppm_param_map_send,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'sending ghosts failed',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_map_field_ghost(fdta,gmm_topoid,gmm_meshid,ghostsize, &
     &                     ppm_param_map_pop,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_march',  &
     &        'popping field data failed',__LINE__,info)
          GOTO 9999
      ENDIF
#endif 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_march_bkwd',t0,info)
      RETURN
#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_bkwd_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_bkwd_2dd
#endif 

#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_march_bkwd_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_march_bkwd_3dd
#endif 
#endif
