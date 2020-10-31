!------------------------------------------------------------------------------!
!* filename: wvic_rhs_loca                                                   *!
!* project : ppm                                                              *!
!* purpose : right hand side for positions                                    *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 14:29:14 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!------------------------------------------------------------------------------!
! $Log: wvic_rhs_loca.F,v $
! Revision 1.3  2006/10/24 16:52:53  pchatela
! Added KEEPWCENTERED option if one wants to keep vorticity centered in domain
!
! Revision 1.2  2006/09/16 00:22:03  pchatela
! Implemented the kinetic energy spectrum, dumped into an ascii file.
!
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.3  2005/11/11 17:16:01  michaebe
! removed m2p p2m interfaces
!
! Revision 1.2  2005/11/11 14:04:25  michaebe
! clean up, additions, comments
!
! Revision 1.1  2005/09/28 11:40:42  michaebe
! Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!  right hand side function conforming with ppm ode solver standard
!------------------------------------------------------------------------------!
FUNCTION wvic_rhs_loca (vdummy, vxp, vup, vdime, vnp, &
     & ipack, lpack, rpack, info)
  
  USE module_wvic
  USE ppm_module_write
  USE ppm_module_data
  USE ppm_module_interp_m2p
  IMPLICIT NONE
  ! interface needed because pointer
  
  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(in)                        :: vdime, vnp
  REAL(mk),DIMENSION(:,:), POINTER           :: vdummy, vxp, vup
  REAL(MK),DIMENSION(:,:), POINTER,OPTIONAL  :: rpack
  INTEGER, DIMENSION(:,:), POINTER,OPTIONAL  :: ipack
  LOGICAL, DIMENSION(:,:), POINTER,OPTIONAL  :: lpack  
  INTEGER, INTENT(inout)            :: info
  
  !----------------------------------------------------------------------------!
  !  Return values
  !----------------------------------------------------------------------------!
  INTEGER                                    :: wvic_rhs_loca,kp
  CHARACTER(len=256)                         :: msg
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  !  none.

  
  IF(verbose) THEN
     WRITE(msg,*) 'entered'
     CALL ppm_write(rank,'wvic_rhs_loca',msg,info)
  END IF

  !----------------------------------------------------------------------------!
  !  Interpolate the velocity onto the particles (velocity field ghosts have
  !  already been sync'd in wvic_rhs_vort()
  !----------------------------------------------------------------------------!
  CALL ppm_interp_m2p(vxp,np,vup,dime,topo_id,mesh_id,&
       &              ppm_param_rmsh_kernel_mp4,&
       &              ghostsize, field_up,info)
  
  WRITE(msg,*) 'particle max vel: ',SQRT(MAXVAL((vup(1,:)**2 + vup(2,:)**2 + vup(3,:)**2))), &
           & ' at index ', MAXLOC((vup(1,:)**2 + vup(2,:)**2 + vup(3,:)**2)), & 
           & ' : xp = ', vxp(:,MAXLOC((vup(1,:)**2 + vup(2,:)**2 + vup(3,:)**2))), &
           & ', up = ',  vup(:,MAXLOC((vup(1,:)**2 + vup(2,:)**2 + vup(3,:)**2)))
  IF(verbose) CALL ppm_write(rank,'wvic_rhs_loca',msg,info)
  
  WRITE(msg,*) 'particle max velx: ',MAXVAL(ABS(vup(1,:))), ' at xp = ', vxp(:,(MAXLOC(ABS(vup(1,:)))))
  IF(verbose) CALL ppm_write(rank,'wvic_rhs_loca',msg,info)
  WRITE(msg,*) 'particle max vely: ',MAXVAL(ABS(vup(2,:))), ' at xp = ', vxp(:,(MAXLOC(ABS(vup(2,:)))))
  IF(verbose) CALL ppm_write(rank,'wvic_rhs_loca',msg,info)
  WRITE(msg,*) 'particle max velz: ',MAXVAL(ABS(vup(3,:))), ' at xp = ', vxp(:,(MAXLOC(ABS(vup(3,:)))))
  IF(verbose) CALL ppm_write(rank,'wvic_rhs_loca',msg,info)
  

  !-----------------------------------------------------
  ! correct for center of mass if requested
  !-----------------------------------------------------
  IF (keepwcentered) THEN
     DO kp=1,np
        vup(1,kp) = vup(1,kp) - u_cmass(1)
        vup(2,kp) = vup(2,kp) - u_cmass(2)
        vup(3,kp) = vup(3,kp) - u_cmass(3)
     END DO
  END IF
  wvic_rhs_loca = 12
  
END FUNCTION wvic_rhs_loca
  
  
  

  
  
