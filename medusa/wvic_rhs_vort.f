!------------------------------------------------------------------------------!
!* filename: wvic_rhs_vort                                                   *!
!* project : ppm                                                              *!
!* purpose : right hand side for the vorticity                                *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 14:29:14 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!------------------------------------------------------------------------------!
! $Log: wvic_rhs_vort.F,v $
! Revision 1.10  2006/10/23 08:19:11  pchatela
! LES models
! Bugfixes in KE spectra and factor for Parseval identity
! Removed the reset of noise in init_tvphysics_0 and_1
!
! Revision 1.9  2006/10/17 14:17:40  pchatela
! Added Ctrl parameters to handle LES models
!
! Revision 1.8  2006/10/16 15:38:46  pchatela
! Added divergence free (analytical) wvic_init_tvphysics_1.F
!
! Revision 1.7  2006/09/27 09:30:21  pchatela
! Fixes, spectra calculation,
! most importantly: moved the u_infty out, so it does not kill the dgammadt
!
! Revision 1.6  2006/09/16 00:22:03  pchatela
! Implemented the kinetic energy spectrum, dumped into an ascii file.
!
! Revision 1.5  2006/09/03 14:34:15  pchatela
! Added trivial routine to add U_infty in case of spectral velocities
! Fixed inverse transforms in the fft_velocities_bgw, the inverse transform should depend on the qty considered
!
! Revision 1.4  2006/08/29 13:40:29  menahel
! uncommented the enforcing and added a conditional call to
! ppm_interp_p2m_reno
!
! Revision 1.3  2006/08/24 11:25:07  menahel
! commented the enforcetv
!
! Revision 1.2  2006/08/23 09:54:50  pchatela
! Added calls to wvic_enforcetv0 in the case of the trailing vortex
! formation
!
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.6  2005/12/10 20:13:02  michaebe
! streamfunction mapping bugfix and clean up
!
! Revision 1.5  2005/11/21 17:36:39  michaebe
! major session. now let''s debug
!
! Revision 1.4  2005/11/11 20:00:50  michaebe
! added call to dgammadt extrapolation routine
!
! Revision 1.3  2005/11/11 17:14:38  michaebe
! added interface for allocation call
!
! Revision 1.2  2005/11/11 14:04:25  michaebe
! clean up, additions, comments
!
! Revision 1.1  2005/09/28 11:40:43  michaebe
! Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!  right hand side function conforming with ppm ode solver standard
!------------------------------------------------------------------------------!
FUNCTION wvic_rhs_vort ( vxp, vwp, vdwp, vdime, vnp, ipack, lpack, rpack, &
                       & rkstep, info)

  USE module_wvic
  USE ppm_module_write
  USE ppm_module_interp_p2m
  USE ppm_module_interp_m2p
  USE ppm_module_data
  IMPLICIT NONE
  !----------------------------------------------------------------------------!
  !  Interfaces
  !----------------------------------------------------------------------------!
  INTERFACE
     SUBROUTINE wvic_alloc_field (vfield_up, ilda, info)
       USE module_wvic
       IMPLICIT NONE
       REAL (mk), DIMENSION (:, :, :, :, :), POINTER :: vfield_up
       INTEGER                          , INTENT(in) :: ilda
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field
  END INTERFACE
  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(in)               :: vdime, vnp
  REAL(mk),DIMENSION(:,:),POINTER     :: vxp, vwp, vdwp
  REAL(MK),DIMENSION(:,:), POINTER,OPTIONAL  :: rpack
  INTEGER, DIMENSION(:,:), POINTER,OPTIONAL  :: ipack
  LOGICAL, DIMENSION(:,:), POINTER,OPTIONAL  :: lpack  
  INTEGER, INTENT(inout)            :: info
  INTEGER, INTENT(in)               :: rkstep
  INTEGER :: imax(2)
  CHARACTER(len=8) :: file
  !----------------------------------------------------------------------------!
  !  Return values
  !----------------------------------------------------------------------------!
  INTEGER :: wvic_rhs_vort
  INTEGER :: i,j,isub,isubl
  CHARACTER(len=256) :: msg
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  !  none.
  CALL wvic_pbc(info)
  IF(verbose) THEN
     WRITE(msg,*) 'endered, now p2m-ing'
     CALL ppm_write(rank,'wvic_rhs_vort',msg,info)
  END IF

  !----------------------------------------------------------------------------!
  !  Interpolate the vorticity onto the mesh and put it in field_wp
  !  JTR!: This needs not be done the first RK step
  !----------------------------------------------------------------------------!
  !IF (rkstep .EQ. 0) THEN
    IF (trailvortex) THEN
       CALL ppm_interp_p2m_reno(xp,np,wp,lda,topo_id,mesh_id,&
            &              ppm_param_rmsh_kernel_mp4,&
            &              ghostsize, field_wp,info)
    ELSE
       CALL ppm_interp_p2m(xp,np,wp,lda,topo_id,mesh_id,&
            &              ppm_param_rmsh_kernel_mp4,&
            &              ghostsize, field_wp,info)
    END IF
  
    IF (trailvortex) THEN
       SELECT CASE(flow_case)
       CASE(0)
          !-----------------------------------------------------
          !  4 straight tubes
          !-----------------------------------------------------
          CALL wvic_enforcetv0
       CASE(1)
          !-----------------------------------------------------
          !  4 tubes with divergence-free perturbation 
          !-----------------------------------------------------
          CALL wvic_enforcetv1
       CASE(2)
          !-----------------------------------------------------
          !  4 tubes with divergence free perturbation at x=0
          !-----------------------------------------------------
          CALL wvic_enforcetv2
       END SELECT
    END IF
  
    isub = 1
    CALL wvic_ghost(wvic_prm_vorticity,info)
    
    IF(verbose) THEN
       WRITE(msg,*) 'p2m - complete, now poisson...'
       CALL ppm_write(rank,'wvic_rhs_vort',msg,info)
    END IF
  !END IF! rkstep

  !----------------------------------------------------------------------------!
  !  Solve the poisson equation for the stream function... 
  !----------------------------------------------------------------------------!
  !IF (rkstep .LT. 2) THEN
    IF(wvic_multigrid) THEN
#ifdef __WITH_FISHPACK__
       CALL wvic_poisson_fishpack(info)
#else
  !#error you dont want to do this
       CALL wvic_poisson_mg(info)                 
#endif     
    ELSE
       CALL wvic_poisson_fft(info)
    END IF
  
  !----------------------------------------------------------------------------!
  !  ...and update the ghosts
  !----------------------------------------------------------------------------!
    IF(verbose) THEN
       WRITE(msg,*) 'poisson complete, now ghosting..'
       CALL ppm_write(rank,'wvic_rhs_vort',msg,info)
    END IF

    IF(wvic_multigrid) THEN
       CALL wvic_ghost(wvic_prm_stream,info)        
    END IF

  !----------------------------------------------------------------------------!
  !  Apply wall boundary conditions
  !  [TODO] Account for velocity potential contribution
  !  Fri Nov 11 14:38:07 CET 2005
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!
  !  Compute the velocity on the mesh...
  !  [TODO] Account for contribs of velocity potential
  !  Fri Nov 11 14:38:50 CET 2005
  !----------------------------------------------------------------------------!
  SELECT CASE(wvic_compvel_scheme)
  CASE(0)
     CALL wvic_compvel(info)
  CASE(1)
     CALL wvic_compvel_4th(info)
  CASE(2)
     CALL wvic_compvel_spec(info)
  END SELECT

  if(verbose) WRITE(msg,*) 'velx max before ghosts', MAXVAL(field_up(1,:,:,:,:))
  if(verbose) CALL ppm_write(rank,'wvic_rhs_vort',msg,info)
  if(verbose) WRITE(msg,*) 'vely max before ghosts', MAXVAL(field_up(2,:,:,:,:))
  if(verbose) CALL ppm_write(rank,'wvic_rhs_vort',msg,info)
  if(verbose) WRITE(msg,*) 'velz max before ghosts', MAXVAL(field_up(3,:,:,:,:))
  if(verbose) CALL ppm_write(rank,'wvic_rhs_vort',msg,info)
  !----------------------------------------------------------------------------!
  !  ...and update the ghosts
  !----------------------------------------------------------------------------!
  !JTR? why ghost vorticity?
  CALL wvic_ghost(wvic_prm_vorticity,info)
  !----------------------------------------------------------------------------!
  !  Correct the wall velocity
  !  [TODO] including possibly linear extrapolation for u(x=-h)
  !  OR!! use the vorticity to get it!
  !  Fri Nov 11 14:59:10 CET 2005
  !----------------------------------------------------------------------------!
  CALL wvic_ghost(wvic_prm_velocity,info)

  !END IF! rkstep
  !----------------------------------------------------------------------------!
  !  Allocate field dgammadt
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  !  Compute the stretching term on the mesh and the laplacian...
  !----------------------------------------------------------------------------!
  SELECT CASE(wvic_dgamma_scheme)
  CASE(0)
     CALL wvic_dgammadt(info)                             ! i
  CASE(1)
     CALL wvic_dgammadt_muscle (info)
  CASE(2)
     CALL wvic_dgammadt_pcons (info)
  CASE(3)
     CALL wvic_dgammadt_ncons (info)
  CASE(4)
     CALL wvic_dgammadt_4th (info)
  END SELECT
  
  !----------------------------------------------------------------------------!
  !  LES model
  !----------------------------------------------------------------------------!
  !JTR? has a turbulence model already been implemented?
  SELECT CASE(wvic_les_model)
  CASE(1)
     CALL wvic_les_tdm(info)
  CASE(2)
     CALL wvic_les_tdmclipped(info)
  CASE(3)
     CALL wvic_les_hypervisc(info)
  END SELECT
  
  !----------------------------------------------------------------------------!
  !  Add the u_infty
  !----------------------------------------------------------------------------!
  !IF (rkstep .LT. 2) THEN
    CALL wvic_compvel_adduinfty(info)
  !END IF
  !----------------------------------------------------------------------------!
  !  Impose explicit penalization
  !----------------------------------------------------------------------------!
  !JTR make logical that determines if explicit penalization should be used instead of testing the flow_case
  IF ((flow_case .gt. 8) .AND. (penalization_implicit .EQ. .FALSE.)) THEN
    CALL wvic_penalization_explicit(info) 
  END IF

  !----------------------------------------------------------------------------!
  !  ...and interpolate that guy back on the particles
  !----------------------------------------------------------------------------!
  CALL wvic_ghost(wvic_prm_dgammadt,info)
  CALL ppm_interp_m2p(xp,np,dwp,lda,topo_id,mesh_id,&
       &              ppm_param_rmsh_kernel_mp4,&
       &              ghostsize, field_dwp,info)

  wvic_rhs_vort = 11
  
END FUNCTION wvic_rhs_vort

  
  

  
  
