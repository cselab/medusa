      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_gmm_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine finalizes the ppm_gmm module and 
      !                 deallocates all data structures.
      !
      !  Input        : 
      !
      !  Input/output : 
      !
      !  Output       : info            (I) return status. 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_gmm_finalize.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.4  2005/05/10 04:48:44  ivos
      !  Split marching and extension routines for faster compilation,
      !  Sharked extension routines, moved all initialization to gmm_init, and
      !  code cosmetics.
      !
      !  Revision 1.3  2005/04/21 04:49:57  ivos
      !  bugfix: pointers to array slabs were incorrectly used.
      !
      !  Revision 1.2  2005/03/16 06:20:07  ivos
      !  Several bugfixes. 1st order version is now tested. Moved all large
      !  data to the module.
      !
      !  Revision 1.1  2005/03/10 01:38:43  ivos
      !  Initial check-in.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_gmm_finalize(info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_gmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)  :: ldu
      INTEGER                :: iopt
      REAL(ppm_kind_double)  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_gmm_finalize',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_gmm_finalize',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate work space structures
      !-------------------------------------------------------------------------
      gmm_lsiz = -1
      iopt = ppm_param_dealloc
      CALL ppm_alloc(gmm_phis,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'sparse data values GMM_PHIS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_phid,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'sparse data values GMM_PHID',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_ipos,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'sparse data locations GMM_IPOS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clod,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clos,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clod2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOD2',__LINE__,info)
      ENDIF
      CALL ppm_alloc(gmm_clos2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_gmm_finalize',     &
     &        'close point locations GMM_CLOS2',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_gmm_finalize',t0,info)
      RETURN

      END SUBROUTINE ppm_gmm_finalize
