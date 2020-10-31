      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_fmm_pretraverse 
      !-------------------------------------------------------------------------
      !
      !  Purpose      :    Do a pretraversal for the target points, to determine
      !                    which particles are needed for the computation. 
      !
      !  Input        :    tp(:,:)      (F)  :  target points
      !                    Ntp          (I)  :  nr of target points
      !                    topoid       (I)  :  topolgy id onto which particles
      !                                         are mapped
      !                    
      !
      !  Input/Output :  
      !                     
      !  Output       :    coeff_subtop(:i,:) (I) : array containing needed
      !                                             coefficients from other
      !                                             procs
      !                                             1st index: subid
      !                                             2nd index: topoid
      !                    ccount             (I) : length of coeff_subtop
      !                    part_subtop(:)     (I) : array containing needed
      !                                             particles from other procs
      !                                             1st index: subid
      !                    pcount             (I) : length of part_subtop      
      !                    info               (I) : return status
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_pretraverse.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.14  2005/09/19 13:03:29  polasekb
      !  code cosmetics
      !
      !  Revision 1.13  2005/09/12 13:30:54  polasekb
      !  added ppm_subid
      !
      !  Revision 1.12  2005/09/12 11:38:03  polasekb
      !  changed duplex check to flag-arrays
      !
      !  Revision 1.11  2005/09/11 18:05:31  polasekb
      !  (final?) corrected version
      !  (also works parallel :-)
      !
      !  Revision 1.10  2005/09/11 11:44:46  polasekb
      !  now using correct topoid for leaf boxes
      !
      !  Revision 1.9  2005/09/10 07:50:22  polasekb
      !  changed init of stack for parallel version
      !
      !  Revision 1.8  2005/09/05 06:24:42  polasekb
      !  deallocation of local variables
      !  checking topology
      !
      !  Revision 1.7  2005/08/23 14:35:04  polasekb
      !  added parameter theta (as in potential)
      !
      !  Revision 1.6  2005/08/23 14:32:14  polasekb
      !  corrected acceptance criterion
      !
      !  Revision 1.5  2005/07/29 12:37:32  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.4  2005/07/21 12:42:28  polasekb
      !  adapted to target points
      !
      !  Revision 1.3  2005/06/04 00:28:29  ivos
      !  Fixed syntax error (XLF) in logical comparisons.
      !
      !  Revision 1.2  2005/06/02 14:35:01  polasekb
      !  changed allocation of variables
      !
      !  Revision 1.1  2005/05/27 08:02:14  polasekb
      !  initial implementation
      !
      !  Revision 0  2004/11/17 16:02:03  polasekb
      !  start.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fmm_pretraverse_s(tp,Ntp,topoid,theta,coeff_subtop, &
      &          ccount,part_subtop,pcount,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fmm_pretraverse_d(tp,Ntp,topoid,theta,coeff_subtop, &
      &          ccount,part_subtop,pcount,info)
#endif 


      !------------------------------------------------------------------------- 
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_alloc
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
        
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif


#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:),     POINTER        :: tp
      INTEGER,                      INTENT(IN   )  :: Ntp
      INTEGER,                      INTENT(IN   )  :: topoid
      REAL(MK),                     INTENT(IN   )  :: theta
      INTEGER, DIMENSION(:,:),      POINTER        :: coeff_subtop
      INTEGER,                      INTENT(  OUT)  :: ccount      
      INTEGER, DIMENSION(:  ),      POINTER        :: part_subtop
      INTEGER,                      INTENT(  OUT)  :: pcount    
      INTEGER,                      INTENT(  OUT)  :: info      
      
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      LOGICAL                              :: drct
      LOGICAL,DIMENSION(:),   POINTER      :: flagcoeff,flagpart
      INTEGER                              :: i,j,cnt,level
      INTEGER                              :: root,iopt,istat
      INTEGER                              :: curtopoid,in_topoid
      INTEGER                              :: stackpointer,curbox,cursub
      INTEGER,DIMENSION(1)                 :: ldu1
      INTEGER,DIMENSION(2)                 :: ldu2
      INTEGER,DIMENSION(:),   POINTER      :: lpart_subtop,stack
      INTEGER,DIMENSION(:,:), POINTER      :: lcoeff_subtop
      REAL(MK)                             :: dx,dy,dz,dist,rad,t0 
      REAL(MK),DIMENSION(:),  POINTER      :: radius      
      REAL(MK),DIMENSION(:,:),POINTER      :: centerofbox
      CHARACTER(LEN=256)                   :: cbuf
        

      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
        
      CALL substart('ppm_fmm_pretraverse',t0,info)
        
      !-------------------------------------------------------------------------
      !  Checking precision and pointing to the correct variables
      !-------------------------------------------------------------------------

#if     __KIND == __SINGLE_PRECISION
      centerofbox => centerofbox_s
      radius      => radius_s
#else
      centerofbox => centerofbox_d
      radius      => radius_d
#endif

      !-------------------------------------------------------------------------
      ! Finding the root of the tree (serial)
      ! Finding top level of tree (parallel)
      !-------------------------------------------------------------------------

      IF (ppm_nproc .GT. 1) THEN
         ! finding top level
         DO i=1,nlevel
           IF (nbpl(i) .GE. ppm_nproc) THEN
              level = i
              EXIT
           ENDIF
         ENDDO
      ELSE
        IF (parent(1) .EQ. ppm_param_undefined) THEN
           root = 1
        ELSE
           DO i=1,nbox
              IF (parent(i) .EQ. ppm_param_undefined) THEN
                 root = i
                 EXIT
              ENDIF
           ENDDO
        ENDIF
      ENDIF      

      ccount = 0
      pcount = 0

      !-------------------------------------------------------------------------
      !  Allocating and initializing local variables
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit
      istat = 0
      ldu1(1) = nbox
      CALL ppm_alloc(stack,ldu1,iopt,info)
      istat = istat + info
      CALL ppm_alloc(lpart_subtop,ldu1,iopt,info)
      istat = istat + info
      CALL ppm_alloc(flagcoeff,ldu1,iopt,info)
      istat = istat + info
      CALL ppm_alloc(flagpart,ldu1,iopt,info)
      istat = istat + info

      ldu2(1) = 2
      ldu2(2) = nbox
      CALL ppm_alloc(lcoeff_subtop,ldu2,iopt,info)
      istat = istat + info

      IF (istat .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_pretraverse', &
      &       'error allocating variables',__LINE__,info)
      GOTO 9999
      ENDIF

      DO i=1,nbox
        stack(i)           = 0
        lpart_subtop(i)    = 0
        lcoeff_subtop(1,i) = 0
        lcoeff_subtop(2,i) = 0
        flagcoeff(i)       = .FALSE.
        flagpart(i)        = .FALSE.
      ENDDO

      !-------------------------------------------------------------------------
      ! Traversing tree and building list part_subtop and coeff_subtop
      !-------------------------------------------------------------------------
      DO i=1,Ntp      
         IF (ppm_nproc .GT. 1) THEN
!          init stack parallel
           stackpointer = 1
           cnt = 0
           DO j=1,nbox
             IF (blevel(j) .EQ. level) THEN
               stack(stackpointer) = j
               stackpointer = stackpointer + 1
               cnt = cnt +1
               IF (cnt .EQ. nbpl(level)) THEN
                  EXIT
               ENDIF
             ENDIF
           ENDDO
         ELSE
!         init stack serial
          stackpointer = 1
          stack(stackpointer) = root
          stackpointer = stackpointer + 1
         ENDIF
         DO WHILE (stackpointer .GT. 1)
           curbox = stack(stackpointer-1)
           dx = tp(1,i) - centerofbox(1,curbox)
           dy = tp(2,i) - centerofbox(2,curbox)
           dz = tp(3,i) - centerofbox(3,curbox)
           dist = SQRT(dx*dx + dy*dy + dz*dz)
           !pop top box
           stackpointer = stackpointer -1
           !--------------------------------------------------------------------
           ! Checking Barnes-Hut Criterium
           !--------------------------------------------------------------------
           IF (radius(curbox) .EQ. 0) THEN
              !only one particle in box, do direct computation
              drct = .TRUE.
           ENDIF
           IF ((dist/(2*radius(curbox)) .GT. theta) .AND. (.NOT. drct)) THEN
              curtopoid = ppm_internal_topoid(blevel(curbox))
              cursub = ppm_subid(curbox,blevel(curbox))
              !far enough, compute part-box interaction
              IF (ppm_rank .NE. ppm_subs2proc(cursub,curtopoid)) THEN
                 ! check if its a duplicate
                 IF (.NOT. flagcoeff(curbox)) THEN
                   ccount                  = ccount+1
                   lcoeff_subtop(1,ccount) = cursub
                   lcoeff_subtop(2,ccount) = blevel(curbox)
                   flagcoeff(curbox)       = .TRUE.
                 ENDIF
              ENDIF
              !ELSE : ok, on same processor
           ELSE
             !not far enough, push childern
             IF (nchld(curbox) .GT. 0) THEN
               DO j=1,nchld(curbox)
                  IF (stackpointer .GT. nbox) THEN
                     CALL ppm_write(ppm_rank,'ppm_fmm_pretraverse', &
                        &         'stack overflow',info)
                  ENDIF
                  stack(stackpointer) = child(j,curbox)
                  stackpointer = stackpointer + 1
               ENDDO
             ELSE
               cursub = ppm_subid(curbox,topoid)
               in_topoid = ppm_internal_topoid(topoid)
               !no children, direct computation
               !particles only needed on the finest topolgy
               IF (ppm_rank .NE. ppm_subs2proc(cursub,in_topoid)) THEN
                 IF (.NOT. flagpart(curbox)) THEN
                   pcount               = pcount +1               
                   lpart_subtop(pcount) = cursub
                   flagpart(curbox)     = .TRUE.
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDDO 
      ENDDO      
      !-------------------------------------------------------------------------
      !  Allocating the correct size of return variables 
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu2(1) = 2
      ldu2(2) = ccount
      CALL ppm_alloc(coeff_subtop,ldu2,iopt,info)
      IF (info .NE. 0) THEN 
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_pretraverse', &
      &        'error allocating coeff_subtop',__LINE__,info)
      GOTO 9999
      ENDIF
      ldu1(1) = pcount
      CALL ppm_alloc(part_subtop,ldu1,iopt,info)
      IF (info .NE. 0) THEN 
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_pretraverse', &
      &        'error allocating part_subtop',__LINE__,info)
      GOTO 9999
      ENDIF

      coeff_subtop(:,:) = lcoeff_subtop(:,1:ccount)
      part_subtop(:)  = lpart_subtop(1:pcount)

      !-------------------------------------------------------------------------
      !  Deallocating local variables
      !-------------------------------------------------------------------------

      istat     = 0
      ldu1(1)   = 0
      ldu2(1:2) = 0
      CALL ppm_alloc(stack,ldu1,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(lpart_subtop,ldu1,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(lcoeff_subtop,ldu2,ppm_param_dealloc,info)
      istat = istat + info
      IF (istat .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expansion', &
      &       'error deallocating newlpdx',__LINE__,info)
      GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_fmm_pretraverse',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fmm_pretraverse_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fmm_pretraverse_d
#endif
