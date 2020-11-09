      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_fmm_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine deallocates all the arrays 
      !                 from ppm_module_data_fmm
      !                 
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : info    (I) 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_finalize.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.9  2005/09/19 13:03:28  polasekb
      !  code cosmetics
      !
      !  Revision 1.8  2005/09/12 13:31:16  polasekb
      !  added ppm_subid
      !
      !  Revision 1.7  2005/08/04 16:03:58  polasekb
      !  some new to deallocate
      !
      !  Revision 1.6  2005/07/29 12:36:28  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.5  2005/07/25 14:43:01  polasekb
      !  typo
      !
      !  Revision 1.4  2005/07/25 14:42:24  polasekb
      !  added new variables for dealloc
      !
      !  Revision 1.3  2005/06/02 14:46:39  polasekb
      !  removed variable totalmass
      !
      !  Revision 1.2  2005/05/27 08:42:01  polasekb
      !  removed dummy argument
      !
      !  Revision 1.1  2005/05/27 08:02:42  polasekb
      !  initial implementation
      !
      !  Revision 0  2004/11/16 14:36:49  polasekb
      !  start
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_fmm_finalize(info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
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
      INTEGER, INTENT(INOUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1)             :: lda1 
      INTEGER, DIMENSION(2)             :: lda2 
      INTEGER, DIMENSION(3)             :: lda3 
      INTEGER                           :: iopt,istat
      INTEGER                           :: i,j
      REAL(8)                           :: t0
      CHARACTER(LEN=ppm_char)           :: mesg

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_fmm_finalize',t0,info)

      lda1(1)=0
      lda2(1:2) = 0
      lda3(1:3) = 0

      !-------------------------------------------------------------------------
      !  Deallocate global arrays (from ppm_module_data_fmm)
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Deallocate single/double variables
      !-------------------------------------------------------------------------
      istat = 0
      iopt = ppm_param_dealloc
     
!not used, too complicated
!#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(radius_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(expansion_s,lda3,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(centerofbox_s,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(totalmass_s,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(min_box_s,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(max_box_s,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(boxcost_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Anm_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(sqrtfac_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(fracfac_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(fac_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(rho_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(theta_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(phi_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Cnm_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Pnm_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Ynm_s,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Inner_s,lda1,iopt,info)
      istat=istat+info  
!#else
      CALL ppm_alloc(radius_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(expansion_d,lda3,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(centerofbox_d,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(min_box_d,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(max_box_d,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(boxcost_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Anm_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(sqrtfac_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(fracfac_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(fac_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(rho_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(theta_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(phi_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Cnm_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Pnm_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Ynm_d,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(Inner_d,lda1,iopt,info)
      istat=istat+info  
!#endif

      !-------------------------------------------------------------------------
      !  Deallocate integer variables
      !-------------------------------------------------------------------------
      CALL ppm_alloc(lhbx,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(lpdx,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(nchld,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(parent,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(child,lda2,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(blevel,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(nbpl,lda1,iopt,info)
      istat=istat+info  
      CALL ppm_alloc(ppm_boxid,lda2,iopt,info)
      istat=istat+info        
      CALL ppm_alloc(ppm_subid,lda2,iopt,info)
      istat=istat+info        
      CALL ppm_alloc(boxpart,lda1,iopt,info)
      istat=istat+info        

      IF (istat .NE. 0) THEN
          WRITE(mesg,'(A,I3,A)') 'for ',istat,'error while deallc' 
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fmm_finalize',mesg,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fmm_finalize',t0,info)
      RETURN
      END SUBROUTINE ppm_fmm_finalize

