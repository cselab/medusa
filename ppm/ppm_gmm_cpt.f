      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_gmm_cpt
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs a closest point transform.
      !                 For each grid point adjacent to the interface,
      !                 the closest point ON the interface is returned.
      !                 ppm_gmm_init must be called BEFORE this routine is
      !                 invoked.
      !
      !  Input        : fdata(...)      (F) Level data. Either rank 3
      !                                     (for 2D scalar fields), or rank
      !                                     4 (for 3D scalar fields).
      !                                     Indices: (i,j,[k],isub).
      !                                     Every zero-crossing is
      !                                     interpreted as an interface.
      !                                     A ghostsize of 1 is needed
      !                                     on all sides which must be
      !                                     filled with the old level
      !                                     function value on input!! Only
      !                                     scalar fdata supported.
      !                 tol             (F) Relative tolerance for the
      !                                     determined distance to the
      !                                     interface. 1E-3 is a good
      !                                     choice. The tolerance is in
      !                                     multiples of grid spacings.
      !                 chi([:],:,:,:,:)(F) rank 5 (3d) or rank 4 (2d)
      !                                     field specifying the positions
      !                                     of the grid nodes. 1st index:
      !                                     1..ppm_dim, then i,j,[k],isub.
      !                                     OPTIONAL. Uniform grid is
      !                                     assumed if absent. Ghostlayers
      !                                     of size >=1 must be pre-filled.
      !
      !  Input/output : 
      !
      !  Output       : npts            (I) Number of unique grid points
      !                                     adjacent to the interface.
      !                 ipts(:,:)       (I) Mesh indices of these
      !                                     points. 1st index:
      !                                     i,j,(k),isub (local sub ID);
      !                                     2nd: 1...npts. Will be
      !                                     allocated by this routine.
      !                 closest(:,:)    (F) Locations of the closest
      !                                     points ON the interface. 1st
      !                                     index: x,y,(z), 2nd:
      !                                     1...npts. Will be allocated
      !                                     by this routine.
      !                 info            (I) return status. 0 on success.
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_cpt.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.8  2005/07/14 19:58:11  ivos
      !  Added OPTIONAL argument chi for mesh node positions in distorted
      !  (mapped) meshes. For use with AGM for example.
      !
      !  Revision 1.7  2005/06/06 20:33:29  ivos
      !  Added pointer NULLIFication at the end to restore alloccation status.
      !
      !  Revision 1.6  2005/06/01 05:19:05  ivos
      !  Updated to new interface of gmm_kickoff.
      !
      !  Revision 1.5  2005/04/27 01:06:10  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.4  2005/04/21 04:49:56  ivos
      !  bugfix: pointers to array slabs were incorrectly used.
      !
      !  Revision 1.3  2005/03/16 06:20:07  ivos
      !  Several bugfixes. 1st order version is now tested. Moved all large
      !  data to the module.
      !
      !  Revision 1.2  2005/03/12 04:08:34  ivos
      !  Misc bug fixes.
      !
      !  Revision 1.1  2005/03/11 21:09:10  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_cpt_2ds(fdata,tol,npts,ipts,closest,   &
     &    info,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_cpt_2dd(fdata,tol,npts,ipts,closest,   &
     &    info,chi)
#endif 

#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_gmm_cpt_3ds(fdata,tol,npts,ipts,closest,   &
     &    info,chi)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_gmm_cpt_3dd(fdata,tol,npts,ipts,closest,   &
     &    info,chi)
#endif 
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

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
      USE ppm_module_util_qsort
      USE ppm_module_gmm_kickoff
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __2D
      REAL(MK), DIMENSION(:,:,:)     , POINTER             :: fdata
      REAL(MK), DIMENSION(:,:,:,:)   , INTENT(IN), OPTIONAL:: chi
#elif __DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:)   , POINTER             :: fdata
      REAL(MK), DIMENSION(:,:,:,:,:) , INTENT(IN), OPTIONAL:: chi
#endif
      INTEGER , DIMENSION(:,:)       , POINTER             :: ipts
      REAL(MK), DIMENSION(:,:)       , POINTER             :: closest
      REAL(MK)                       , INTENT(IN   )       :: tol
      INTEGER                        , INTENT(  OUT)       :: info,npts
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                          :: i,iopt,Nt,isub
      INTEGER                          :: n1,n2,n3,jsub,prev
      INTEGER, DIMENSION(2)            :: ldu
      REAL(MK),DIMENSION(:,:), POINTER :: clotmp
      REAL(MK)                         :: t0,x,y,z,xx,yy,zz,dx,dy,dz
      REAL(MK)                         :: s,sprev,thresh
      LOGICAL                          :: lok
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_cpt',t0,info)
#if   __KIND == __SINGLE_PRECISION
      clotmp => gmm_clos
#elif __KIND == __DOUBLE_PRECISION
      clotmp => gmm_clod
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_cpt',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (tol .LE. 0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'tolerance must be >0!',__LINE__,info)
              GOTO 9999
          ENDIF
#if   __DIM == __3D
          IF (SIZE(fdata,4) .LT. ppm_nsublist(gmm_topoid)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,3) .LT. maxzhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'z dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#elif __DIM == __2D
          IF (SIZE(fdata,3) .LT. ppm_nsublist(gmm_topoid)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'field data for some subs is missing',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,1) .LT. maxxhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'x dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (UBOUND(fdata,2) .LT. maxyhi) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_gmm_cpt',  &
     &            'y dimension of field data does not match mesh',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
      ENDIF          ! ppm_debug for argument check

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

      !-------------------------------------------------------------------------
      !  Get non-unique closest points. We assume that the largest occuring
      !  values in fdata are outside of the narrow band.
      !-------------------------------------------------------------------------
      Nt     = -1        ! leave fdata untouched !
      thresh = 0.99_MK*MAXVAL(ABS(fdata))
      IF (PRESENT(chi)) THEN
          CALL ppm_gmm_kickoff(fdata,tol,thresh,info,Nt,iptstmp,clotmp,chi=chi)
      ELSE
          CALL ppm_gmm_kickoff(fdata,tol,thresh,info,Nt,iptstmp,clotmp)
      ENDIF
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_cpt',  &
     &        'Starting GMM failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Hash mesh coordinates to key
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = Nt
      CALL ppm_alloc(key,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_cpt',        &
     &        'sort key KEY',__LINE__,info)
          GOTO 9999
      ENDIF
      n1 = maxxhi
      n2 = maxxhi*maxyhi
      n3 = n2*maxzhi
      DO i=1,Nt
          key(i) = iptstmp(1,i) + (iptstmp(2,i)-1)*n1 + (iptstmp(3,i)-1)*n2
#if   __DIM == __3D
          key(i) = key(i) + (iptstmp(4,i)-1)*n3
#endif
      ENDDO

      !-------------------------------------------------------------------------
      !  Sort 
      !-------------------------------------------------------------------------
      CALL ppm_util_qsort(key,idx,info)
      IF (info .NE. ppm_param_success) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_gmm_cpt',  &
     &        'Sorting failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the number of unique points 
      !-------------------------------------------------------------------------
      npts = 0
      prev = -1
      DO i=1,Nt
          IF (key(idx(i)) .NE. prev) THEN
              npts = npts + 1
              prev = key(idx(i))
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for unique CPT
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_dim
      ldu(2) = npts
      CALL ppm_alloc(closest,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_cpt',        &
     &        'closest point transform CLOSEST',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_dim+1
      ldu(2) = npts
      CALL ppm_alloc(ipts,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_cpt',        &
     &        'closest point transform CLOSEST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store unique CPT
      !-------------------------------------------------------------------------
      npts = 0
      prev = -1
#if   __DIM == __3D
      DO i=1,Nt
          IF (key(idx(i)) .NE. prev) THEN
              npts = npts + 1
              ipts(1,npts) = iptstmp(1,idx(i))
              ipts(2,npts) = iptstmp(2,idx(i))
              ipts(3,npts) = iptstmp(3,idx(i))
              ipts(4,npts) = iptstmp(4,idx(i))
              prev  = key(idx(i))
              sprev = HUGE(sprev)
          ENDIF
          isub = ipts(4,npts)
          jsub = ppm_isublist(isub,gmm_topoid)
          IF (PRESENT(chi)) THEN
              x = chi(1,ipts(1,npts),ipts(2,npts),ipts(3,npts),isub)
              y = chi(2,ipts(1,npts),ipts(2,npts),ipts(3,npts),isub)
              z = chi(3,ipts(1,npts),ipts(2,npts),ipts(3,npts),isub)
          ELSE
              IF (ppm_kind .EQ. ppm_kind_single) THEN
                  x = ppm_min_subs(1,jsub,gmm_topoid) + (ipts(1,npts)-1)*dx
                  y = ppm_min_subs(2,jsub,gmm_topoid) + (ipts(2,npts)-1)*dy
                  z = ppm_min_subs(3,jsub,gmm_topoid) + (ipts(3,npts)-1)*dz
              ELSE
                  x = ppm_min_subd(1,jsub,gmm_topoid) + (ipts(1,npts)-1)*dx
                  y = ppm_min_subd(2,jsub,gmm_topoid) + (ipts(2,npts)-1)*dy
                  z = ppm_min_subd(3,jsub,gmm_topoid) + (ipts(3,npts)-1)*dz
              ENDIF
          ENDIF
          xx = clotmp(1,idx(i))
          yy = clotmp(2,idx(i))
          zz = clotmp(3,idx(i))
          s  = (xx-x)*(xx-x) + (yy-y)*(yy-y) + (zz-z)*(zz-z)
          IF (s .LT. sprev) THEN
              closest(1,npts) = xx
              closest(2,npts) = yy
              closest(3,npts) = zz
              sprev = s
          ENDIF
      ENDDO
#elif __DIM == __2D
      DO i=1,Nt
          IF (key(idx(i)) .NE. prev) THEN
              npts = npts + 1
              ipts(1,npts) = iptstmp(1,idx(i))
              ipts(2,npts) = iptstmp(2,idx(i))
              ipts(3,npts) = iptstmp(3,idx(i))
              prev  = key(idx(i))
              sprev = HUGE(sprev)
          ENDIF
          isub = ipts(3,npts)
          jsub = ppm_isublist(isub,gmm_topoid)
          IF (PRESENT(chi)) THEN
              x = chi(1,ipts(1,npts),ipts(2,npts),isub)
              y = chi(2,ipts(1,npts),ipts(2,npts),isub)
          ELSE
              IF (ppm_kind .EQ. ppm_kind_single) THEN
                  x = ppm_min_subs(1,jsub,gmm_topoid) + (ipts(1,npts)-1)*dx
                  y = ppm_min_subs(2,jsub,gmm_topoid) + (ipts(2,npts)-1)*dy
              ELSE
                  x = ppm_min_subd(1,jsub,gmm_topoid) + (ipts(1,npts)-1)*dx
                  y = ppm_min_subd(2,jsub,gmm_topoid) + (ipts(2,npts)-1)*dy
              ENDIF
          ENDIF
          xx = clotmp(1,idx(i))
          yy = clotmp(2,idx(i))
          s  = (xx-x)*(xx-x) + (yy-y)*(yy-y)
          IF (s .LT. sprev) THEN
              closest(1,npts) = xx
              closest(2,npts) = yy
              sprev = s
          ENDIF
      ENDDO
#endif

      !-------------------------------------------------------------------------
      !  Free temp memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(iptstmp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_cpt',        &
     &        'temporary point list IPTSTMP',__LINE__,info)
      ENDIF
      CALL ppm_alloc(clotmp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_cpt',        &
     &        'temporary closest points CLOTMP',__LINE__,info)
      ENDIF
#if   __KIND == __SINGLE_PRECISION
      NULLIFY(gmm_clos)
#elif __KIND == __DOUBLE_PRECISION
      NULLIFY(gmm_clod)
#endif
      CALL ppm_alloc(idx,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_cpt',        &
     &        'Permutation index IDX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(key,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_gmm_cpt',        &
     &        'sort key KEY',__LINE__,info)
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_cpt',t0,info)
      RETURN
#if    __DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_cpt_2ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_cpt_2dd
#endif 

#elif  __DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_gmm_cpt_3ds
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_gmm_cpt_3dd
#endif 
#endif
