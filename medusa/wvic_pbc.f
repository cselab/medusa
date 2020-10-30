!------------------------------------------------------------------------------!
!* filename: wvic_pbc                                                        *!
!* project : ppm                                                              *!
!* purpose : enforce periodic domain (reinsert particles)                     *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Aug 11 08:14:23 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!------------------------------------------------------------------------------!
!  $Log: wvic_pbc.F,v $
!  Revision 1.4  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.3  2006/08/11 12:32:14  pchatela
!  Now compiles
!
!  Revision 1.2  2006/07/26 07:51:26  pchatela
!  Added boolean trailvortex
!  Added periodic bcs with reset of particle strengths
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.3  2005/11/21 17:36:38  michaebe
!  major session. now let''s debug
!
!  Revision 1.2  2005/11/11 20:07:15  michaebe
!  added wall conditions
!
!  Revision 1.1  2005/09/28 11:40:40  michaebe
!  Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = VOCLE PBC =
! Reinsert particles into the computational domain
!------------------------------------------------------------------------------!
SUBROUTINE wvic_pbc (info)
  
  USE module_wvic
  USE ppm_module_data
  USE ppm_module_error
  USE ppm_module_write
  USE ppm_module_impose_part_bc
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  !  arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out) :: info
  
  !----------------------------------------------------------------------------!
  !  locals
  !----------------------------------------------------------------------------!
  INTEGER              :: kp,dispair
  CHARACTER(len=256)   :: msg,pnm
  REAL(mk), DIMENSION(3) :: len_physg
  REAL(mk), DIMENSION(3) :: xptmp, wptmp, xp0tmp, wp0tmp, uptmp, dwptmp
  REAL(mk)             :: xdel,ydel,zdel
  INTEGER              :: iswap
  INTEGER              :: lp, npold
  
  
  npold = np
  WRITE(pnm,*) 'wvic_pbc'
  
  !----------------------------------------------------------------------------!
  ! obviously (hreidar) my stuff doesnt work anymore, so Ill now use jenss
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!
  ! output bounding box of particles
  !----------------------------------------------------------------------------!
  IF(verbose) WRITE(msg,*) 'bounding box low ',&
       & MINVAL(xp(1,1:np)),MINVAL(xp(2,1:np)),MINVAL(xp(3,1:np))
  IF(verbose) CALL ppm_write(rank,'wvic_pbc',msg,info)
  IF(verbose) WRITE(msg,*) 'bounding box hai ',&
       & MAXVAL(xp(1,1:np)),MAXVAL(xp(2,1:np)),MAXVAL(xp(3,1:np))
  IF(verbose) CALL ppm_write(rank,'wvic_pbc',msg,info)
#ifdef __WRONG
  CALL ppm_impose_part_bc(xp,np,topo_id,info)
#else
  len_physg = max_physg - min_physg
  
  DO dispair=1,2
     !----------------------------------------------------------------------!
     !  x direction
     !----------------------------------------------------------------------!
     DO kp=1,np
        IF(xp(1,kp).GE.max_physg(1).AND..NOT.wbcdef(2)) THEN
           xp(1,kp) = xp(1,kp) - len_physg(1)
        ELSE
           IF(xp(1,kp).LT.min_physg(1).AND..NOT.wbcdef(1)) THEN
              xp(1,kp) = xp(1,kp) + len_physg(1)
           END IF
        END IF
     END DO

     !----------------------------------------------------------------------!
     !  y direction
     !----------------------------------------------------------------------!
     DO kp=1,np
        IF(xp(2,kp).GE.max_physg(2).AND..NOT.wbcdef(4)) THEN
           xp(2,kp) = xp(2,kp) - len_physg(2)
        ELSE
           IF(xp(2,kp).LT.min_physg(2).AND..NOT.wbcdef(3)) THEN
              xp(2,kp) = xp(2,kp) + len_physg(2)
           END IF
        END IF
     END DO
     !----------------------------------------------------------------------!
     !  z direction
     !----------------------------------------------------------------------!
     IF (trailvortex) THEN
        kp = 1     ! Scan index
        lp = np+1  ! Index of first particle to be swapped
        DO 
           IF ((xp(3,kp).GT.max_physg(3)).OR.(xp(3,kp).LT.min_physg(3))) THEN
               ! Needs to be killed
               IF (xp(3,kp).LT.min_physg(3)) THEN
                   IF (SUM(ABS(wp(:,kp))).GT.0.0_mk) THEN
				      IF(verbose) WRITE(msg,*) 'Removing an non-zero upstream particle ',&
       & xp(1,kp), xp(2,kp), xp(3,kp), 'with weight ', wp(1,kp), wp(2,kp), wp(3,kp)
                      IF(verbose) CALL ppm_write(rank,'wvic_pbc',msg,info)
                   END IF
               END IF
               ! We swap, well, more like overwrite this guy with the last one
               lp = lp - 1
               IF (lp.LT.1) THEN
                   info = ppm_error_warning
                   CALL ppm_error(ppm_err_part_lost,pnm,'All particles flew out',__LINE__,info)
               END IF
               xp(:,kp) = xp(:,lp)
               wp(:,kp) = wp(:,lp)
               up(:,kp) = up(:,lp)
               dwp(:,kp) = dwp(:,lp)
               xp0(:,kp) = xp0(:,lp)
               wp0(:,kp) = wp0(:,lp)
           ELSE
               kp = kp + 1
           END IF
           IF (kp.EQ.lp) EXIT
        END DO
        np = kp - 1
     ELSE
        DO kp=1,np
           IF(xp(3,kp).GE.max_physg(3).AND..NOT.wbcdef(6)) THEN
              xp(3,kp) = xp(3,kp) - len_physg(3)
           ELSE
               IF(xp(3,kp).LT.min_physg(3).AND..NOT.wbcdef(5)) THEN
                   xp(3,kp) = xp(3,kp) + len_physg(3)
               END IF
           END IF
        END DO
     END IF
  END DO
#endif

  !----------------------------------------------------------------------------!
  ! output bounding box of particles
  !----------------------------------------------------------------------------!
  IF(verbose) WRITE(msg,*) 'bounding box low ',&
       & MINVAL(xp(1,1:np)),MINVAL(xp(2,1:np)),MINVAL(xp(3,1:np))
  IF(verbose) CALL ppm_write(rank,'wvic_pbc',msg,info)
  IF(verbose) WRITE(msg,*) 'bounding box hai ',&
       & MAXVAL(xp(1,1:np)),MAXVAL(xp(2,1:np)),MAXVAL(xp(3,1:np))
  IF(verbose) CALL ppm_write(rank,'wvic_pbc',msg,info)     
  IF(trailvortex.AND.verbose) WRITE(msg,*) 'np went from ', &
       & npold, ' to ', np
  IF(trailvortex.AND.verbose) CALL ppm_write(rank,'wvic_pbc',msg,info)

9999 CONTINUE  
  
END SUBROUTINE wvic_pbc
