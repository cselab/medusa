#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_alloc_1d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine (re)allocates the memory of
      !                 one-dimensional arrays (pointers) based on the
      !                 number of elements.
      !
      !  Input        : lda(1)   (I) Number of desired elements in leading 
      !                              dimension of array. (>0)
      !                 iopt     (I) allocation mode. One of:
      !                                  ppm_param_alloc_fit
      !                                  ppm_param_alloc_fit_preserve
      !                                  ppm_param_alloc_grow
      !                                  ppm_param_alloc_grow_preserve
      !                                  ppm_param_dealloc
      !
      !  Input/output : adata    (P) Pointer to array which is to be
      !                              (re)allocated.
      !
      !  Output       : info     (I) 0 upon success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_alloc_1d.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.14  2004/11/11 15:22:59  ivos
      !  Moved allocatable work data to the module.
      !
      !  Revision 1.13  2004/10/01 16:33:29  ivos
      !  cosmetics.
      !
      !  Revision 1.12  2004/10/01 16:08:53  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.11  2004/07/26 07:45:20  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.10  2004/07/19 15:47:01  ivos
      !  bugfix: grow and grow_preserve could shrink non-growing dimensions
      !  since they would only pay attention to the growing one.
      !  Fixed by introducing lda_new,ldl_new,ldu_new.
      !
      !  Revision 1.9  2004/02/20 14:23:19  ivos
      !  optimized: reallocation is now done with pointers. this saves half of
      !  the work (1 alloc, 1 dealloc and the copying back).
      !
      !  Revision 1.8  2004/02/19 10:11:28  walther
      !  Cosmetics: changed the order of the defines.
      !
      !  Revision 1.7  2004/02/12 17:35:35  ivos
      !  Bugfix: switched argument checking off if iopt=ppm_param_dealloc.
      !
      !  Revision 1.6  2004/02/06 13:48:53  walther
      !  Added complex variables.
      !
      !  Revision 1.5  2004/01/27 12:52:48  ivos
      !  Relaxed argument checks from ldl.GE.ldu to ldl.GT.ldu and
      !  lda.LE.0 to lda.LT.0.
      !
      !  Revision 1.4  2004/01/26 12:32:29  ivos
      !  Now returns the correct status in info and each allocate or deallocate
      !  is checked for proper execution. renamed data to adata to avoid
      !  name conflicts with Fortran keyword DATA.
      !
      !  Revision 1.3  2004/01/23 17:24:12  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.2  2004/01/22 14:33:37  ivos
      !  Did (1) insert comment headers, (2) update address in header, (3) added
      !  validity check for arguments where needed, (4) uses ppm_error now.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_alloc_1ds(adata,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_alloc_1dd(adata,lda,iopt,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_1dsc(adata,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_alloc_1ddc(adata,lda,iopt,info)
#elif __KIND == __INTEGER
      SUBROUTINE ppm_alloc_1di(adata,lda,iopt,info)
#elif __KIND == __LOGICAL
      SUBROUTINE ppm_alloc_1dl(adata,lda,iopt,info)
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      REAL(ppm_kind_single)     , DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)     , DIMENSION(:), POINTER :: adata
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single)  , DIMENSION(:), POINTER :: adata
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double)  , DIMENSION(:), POINTER :: adata
#elif __KIND == __INTEGER
      INTEGER                   , DIMENSION(:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                   , DIMENSION(:), POINTER :: adata
#elif __KIND == __LOGICAL
      LOGICAL                   , DIMENSION(:), POINTER :: adata
#endif
      INTEGER, DIMENSION(:)     , INTENT(IN)    :: lda
      INTEGER                   , INTENT(IN)    :: iopt
      INTEGER                   , INTENT(OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      REAL(ppm_kind_single)   , DIMENSION(:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double)   , DIMENSION(:), POINTER :: work
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_single), DIMENSION(:), POINTER :: work
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(ppm_kind_double), DIMENSION(:), POINTER :: work
#elif __KIND == __INTEGER
      INTEGER                 , DIMENSION(:), POINTER :: work
#elif __KIND == __LOGICAL
      LOGICAL                 , DIMENSION(:), POINTER :: work
#endif
      INTEGER               :: i,j,k,ldb,ldc
      LOGICAL               :: lcopy,lalloc,lrealloc
      REAL(ppm_kind_double) :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_alloc_1d',t0,info)

      !-------------------------------------------------------------------------
      !  Point to proper work array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      work => work_1ds
#elif __KIND == __DOUBLE_PRECISION
      work => work_1dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      work => work_1dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      work => work_1ddc
#elif __KIND == __INTEGER
      work => work_1di
#elif __KIND == __LOGICAL
      work => work_1dl
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0 .AND. iopt .NE. ppm_param_dealloc) THEN
          IF (lda(1) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_alloc_1d',  &
     &            'lda(1) must be >= 0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check the allocation type
      !-------------------------------------------------------------------------
      lcopy    = .FALSE.
      lalloc   = .FALSE.
      lrealloc = .FALSE.
      IF     (iopt.EQ.ppm_param_alloc_fit_preserve) THEN
         !----------------------------------------------------------------------
         !  fit memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            ldb = SIZE(adata)
            IF (ldb.NE.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF            
         ELSE
            lalloc = .TRUE.
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_fit) THEN
         !----------------------------------------------------------------------
         !  fit memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            ldb = SIZE(adata)
            IF (ldb.NE.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE. 
            ENDIF
         ELSE
            lalloc   = .TRUE. 
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow_preserve) THEN
         !----------------------------------------------------------------------
         !  grow memory and preserve the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            ldb = SIZE(adata)
            IF (ldb.LT.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
               lcopy    = .TRUE.
            ENDIF            
         ELSE
            lalloc = .TRUE.
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_alloc_grow) THEN
         !----------------------------------------------------------------------
         !  grow memory but skip the present contents
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) THEN
            ldb = SIZE(adata)
            IF (ldb.LT.lda(1)) THEN
               lrealloc = .TRUE.
               lalloc   = .TRUE.
            ENDIF            
         ELSE
            lalloc = .TRUE.
         ENDIF
      ELSEIF (iopt.EQ.ppm_param_dealloc) THEN
         !----------------------------------------------------------------------
         !  deallocate
         !----------------------------------------------------------------------
         IF (ASSOCIATED(adata)) DEALLOCATE(adata,STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_dealloc,'ppm_alloc_1d',   &
     &           'DATA',__LINE__,info)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate new memory
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         ALLOCATE(work(lda(1)),STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_alloc_1d',   &
     &           'WORK',__LINE__,info)
             GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Make a copy of the present contents
      !-------------------------------------------------------------------------
      IF (lcopy) THEN
         ldc = MIN(ldb,lda(1))
         DO i=1,ldc
            work(i) = adata(i)
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  If reallocating, deallocate old data
      !-------------------------------------------------------------------------
      IF (lrealloc) THEN
         DEALLOCATE(adata,STAT=info)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_dealloc,'ppm_alloc_1d',   &
     &           'DATA',__LINE__,info)
         ENDIF
      ENDIF 

      !-------------------------------------------------------------------------
      !  Set the pointer to the new array
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
         adata => work
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_alloc_1d',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_alloc_1ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_alloc_1dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_1dsc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_alloc_1ddc
#elif __KIND == __INTEGER
      END SUBROUTINE ppm_alloc_1di
#elif __KIND == __LOGICAL
      END SUBROUTINE ppm_alloc_1dl
#endif
