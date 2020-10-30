!------------------------------------------------------------------------------!
!* filename: wvic_init_physics_0                                              *!
!* project : ppm                                                              *!
!* purpose : initial condition for crow instability                           *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 09:45:17 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_init_physics_0.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.3  2005/11/21 17:36:38  michaebe
!  major session. now let''s debug
!
!  Revision 1.2  2005/11/11 14:04:24  michaebe
!  clean up, additions, comments
!
!  Revision 1.1  2005/11/11 10:01:06  michaebe
!  refurbished
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic INIT PHYSICS 0 =
!  initialize vorticity here
!------------------------------------------------------------------------------!
SUBROUTINE wvic_init_physics_0

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  ! interfaces
  INTERFACE
     SUBROUTINE wvic_alloc_field_s (vfield_up, info)
       USE module_wvic
       IMPLICIT NONE
       REAL (mk), DIMENSION (:, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_s
     SUBROUTINE wvic_alloc_field(vfield_up, ilda, info)
       USE module_wvic
       IMPLICIT NONE
       REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: vfield_up
       INTEGER                   , INTENT(out) :: info
       INTEGER                   , INTENT(in ) :: ilda
     END SUBROUTINE wvic_alloc_field
  END INTERFACE
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  ! localities: noise stuff
  INTEGER, PARAMETER       :: mkd = KIND(1.0d0)
  INTEGER(mkd)             :: an, a32
  REAL(mk), DIMENSION(500) :: aky,bky,akz,bkz
  INTEGER                  :: nk,kk
  REAL(mk)                 :: rnk, amp
  INTEGER, DIMENSION(3)    :: ilo, ihi
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  REAL(mk)                 :: a0, b0, gamma
  REAL(mk)                 :: ty_center, ty_vtx1, ty_vtx2, r_vtx1, r_vtx2
  REAL(mk)                 :: tz_center, tz_vtx1, tz_vtx2
  REAL(mk)                 :: fac1, ra0, tx, ty, tz, noise, typ, tzp, omega_x
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
  REAL(mk)                 :: ink_a0, ink_b0
  REAL(mk)                 :: tz_ink1, ty_ink1
  REAL(mk)                 :: tz_ink2, ty_ink2, pi
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: ig,jg,kg,i,j,k,isub,isubl,kmax
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: mg stuff
  INTEGER,  DIMENSION(6)     :: ibcdef
  INTEGER,  DIMENSION(1,6)   :: ibcdef_v
  REAL(mk), DIMENSION(1,1,1) :: ibcvalue
  REAL(mk), DIMENSION(1,1,1,1) :: ibcvalue_v
  REAL(mk)                     :: cresid
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  REAL(mk)                     :: fac2,fac3,fac4,fac5,fac6,ldiv
  REAL(mk), DIMENSION(:,:,:,:), POINTER :: divergence, potential
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: divergence_v, potential_v
  
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! - generate noise (angle of dislocation and amplitute of dislocation)
  ! - calculate the position of the vortices using b0 (intervortex distance)
  ! - initialize the vortices using a0 (size of vortices)
  ! - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  pi = ACOS(-1.0_mk)
  !###########################################################################
  !----------------------------------------------------------------------------!
  ! =  fixed parameter settings from Leveke:1997
  !----------------------------------------------------------------------------!
  Gamma = 20.2488_mk !/ 8.5_mk
  b0    = 2.1265_mk !/ 8.5_mk
  a0    = 0.5035_mk !/ 8.5_mk
  !###########################################################################


  !###########################################################################
  !----------------------------------------------------------------------------!
  ! =  ink injection parameters
  !----------------------------------------------------------------------------!
  ink_a0 = 0.5035_mk !/ 8.5_mk
  ink_b0 = b0
  !###########################################################################

  

  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_init_physics_0',msg,info)
     STOP ! brutal
  END IF

  isubl = isublist(1)
  ilo = istart(:,isubl)
  ihi = ilo + ndata(:,isubl)-1

  !----------------------------------------------------------------------------!
  ! jump start the random number generator
  !----------------------------------------------------------------------------!
  a32 = 0_mkd
  a32 = IBSET(a32,32)
  an  = 1_mkd
  an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
  !----------------------------------------------------------------------------!
  ! print 1st seed as debug
  !----------------------------------------------------------------------------!
  WRITE(msg,*) 'seed = ',an
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_physics_0',msg,info)

  !----------------------------------------------------------------------------!
  ! compute the random amplitudes
  !----------------------------------------------------------------------------!
  nk = wvic_noise_nmodes
  rnk = 1.0_mk/(2.0_mk*REAL(nk,mk))
  amp = wvic_noise_amp * rnk * b0
  DO i=1,nk

     aky(i) = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     bky(i) = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     akz(i)   = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     bkz(i)   = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)

  END DO
  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  !----------------------------------------------------------------------------!
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))
  ty_vtx1   = ty_center - 0.5_mk * b0
  ty_vtx2   = ty_center + 0.5_mk * b0
  tz_center = 0.5_mk*(min_physg(3) + max_physg(3))
  tz_vtx1   = tz_center
  tz_vtx2   = tz_center
  fac1      = Gamma/(a0**2 * M_PI)
  ra0       = 1.0_mk/a0**2
  length    = max_physg - min_physg
  !----------------------------------------------------------------------------!


  !----------------------------------------------------------------------------!
  ! derived parameters - inkicity 
  !----------------------------------------------------------------------------!
  ty_ink1   = ty_center - 0.5_mk * ink_b0 - 0.5_mk * ink_a0
  ty_ink2   = ty_center + 0.5_mk * ink_b0 - 0.5_mk * ink_a0
  tz_ink1   = tz_center
  tz_ink2   = tz_center
  !----------------------------------------------------------------------------!

  
  !----------------------------------------------------------------------------!
  ! create vortices
  !----------------------------------------------------------------------------!
  max_vorticity = fac1

  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)
        
        DO j=1,ndata(2,isubl)
           
           DO i=1,ndata(1,isubl)
              
              tx = REAL(i+istart(1,isubl)-2,mk)*dx + min_physg(1) 
              ty = REAL(j+istart(2,isubl)-2,mk)*dy + min_physg(2)
              tz = REAL(k+istart(3,isubl)-2,mk)*dz + min_physg(3)
              !----------------------------------------------------------------!
              ! or better
              !----------------------------------------------------------------!
              tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
              ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
              tz = min_sub(3,isubl) + REAL(k-1,mk)*dz
              
              
              omega_x = 2.0_mk * M_PI * tx / length(1)

              !----------------------------------------------------------------!
              ! compute perturbations
              noise = 0.0_mk
              IF(nk.EQ.1) THEN
                 noise = SIN(omega_x)
              ELSE
                 DO kk=1,nk
                    noise = noise &
                         & + SIN(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*aky(kk))) &
                         & + COS(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*bky(kk)))
                 END DO
              END IF
              typ = noise * amp
              noise = 0.0_mk
              IF(nk.EQ.1) THEN
                 noise = 0.0_mk
              ELSE
                 DO kk=1,nk
                    noise = noise &
                         & + SIN(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*akz(kk))) &
                         & + COS(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*bkz(kk)))
                 END DO
              END IF
              tzp = noise * amp
              !----------------------------------------------------------------!

              twp = 0.0_mk

              !----------------------------------------------------------------!
              ! compute vorticity - vtx1 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tz+tzp-tz_vtx1)**2
              twp(1) =  twp(1) + fac1 * EXP(-r_vtx1*ra0)
              ! left neighbor
              r_vtx1 = (ty+typ-ty_vtx1-length(2))**2 + (tz+tzp-tz_vtx1)**2
              twp(1) =  twp(1) + fac1 * EXP(-r_vtx1*ra0)
              ! right neighbor
              r_vtx1 = (ty+typ-ty_vtx1+length(2))**2 + (tz+tzp-tz_vtx1)**2
              twp(1) =  twp(1) + fac1 * EXP(-r_vtx1*ra0)
              ! bottom neighbor
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tz+tzp-tz_vtx1-length(3))**2
              twp(1) =  twp(1) + fac1 * EXP(-r_vtx1*ra0)
              ! top neighbor
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tz+tzp-tz_vtx1+length(3))**2
              twp(1) =  twp(1) + fac1 * EXP(-r_vtx1*ra0)
              !----------------------------------------------------------------!
              ! compute vorticity - vtx2 in ...
              !----------------------------------------------------------------!
              ! center box
              IF(nk.EQ.1) typ = 0.0_mk
              r_vtx2 = (ty+tzp-ty_vtx2)**2 + (tz+typ-tz_vtx2)**2
              twp(1) =  twp(1) - fac1 * EXP(-r_vtx2*ra0)
              ! left neighbor
              r_vtx2 = (ty+tzp-ty_vtx2-length(2))**2 + (tz+typ-tz_vtx2)**2
              twp(1) =  twp(1) - fac1 * EXP(-r_vtx2*ra0)
              ! right neighbor
              r_vtx2 = (ty+tzp-ty_vtx2+length(2))**2 + (tz+typ-tz_vtx2)**2
              twp(1) =  twp(1) - fac1 * EXP(-r_vtx2*ra0)
              ! bottom neighbor
              r_vtx2 = (ty+tzp-ty_vtx2)**2 + (tz+typ-tz_vtx2-length(3))**2
              twp(1) =  twp(1) - fac1 * EXP(-r_vtx2*ra0)
              ! top neighbor
              r_vtx2 = (ty+tzp-ty_vtx2)**2 + (tz+typ-tz_vtx2+length(3))**2
              twp(1) =  twp(1) - fac1 * EXP(-r_vtx2*ra0)

              !-----------------------------------------------------
              !  SET FIELD VORTICITY VALUES
              !-----------------------------------------------------
              field_wp(1,i,j,k,isub) = twp(1)
              field_wp(2,i,j,k,isub) = twp(2)
              field_wp(3,i,j,k,isub) = twp(3)
              !-----------------------------------------------------
              
           END DO

        END DO

     END DO

  END DO
  !---- Vorticity has been initialized on the field.

  !----------------------------------------------------------------------------!
  ! allocate auxilliary fields
  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
  CALL ppm_write(rank,'wvic_init_physics_0','ghosting',info)
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  CALL ppm_write(rank,'wvic_init_physics_0','ghst complete',info)


  
  !----------------------------------------------------------------------------!
  ! create particles with style
  !----------------------------------------------------------------------------!
  CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,topo_id,mesh_id,&
       & (/cutoff,HUGE(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
  WRITE(msg,*) ' created ',np,' particles'
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_physics_0',msg,info)
  !----------------------------------------------------------------------------!
  ! all set
  !----------------------------------------------------------------------------!
  ! we dont like plot3d no more
  ! CALL wvic_dumpfield_plot3d(info)


END SUBROUTINE wvic_init_physics_0
