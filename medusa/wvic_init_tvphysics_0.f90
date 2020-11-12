!------------------------------------------------------------------------------!
!* filename: wvic_init_tvphysics_0                                            *!
!* project : ppm                                                              *!
!* purpose : initial condition for crow instability                           *!
!*         :                                                                  *!
!* author  : Michael Bergdorf  Philippe Chatelain                             *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 09:45:17 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_init_tvphysics_0.F,v $
!  Revision 1.14  2006/10/23 08:19:11  pchatela
!  LES models
!  Bugfixes in KE spectra and factor for Parseval identity
!  Removed the reset of noise in init_tvphysics_0 and_1
!
!  Revision 1.13  2006/10/16 17:24:21  pchatela
!  Projection to solenoidal vorticity field for the periodic case of tbphysics_0
!  spatial diags have to be computed for Flow_case 5 and 6 in the periodic case
!
!  Revision 1.12  2006/10/04 11:16:38  pchatela
!  Added a fundamental mode to the noise parameters.
!  If not given, the length of the domain is used.
!
!  Revision 1.11  2006/09/28 15:35:52  pchatela
!  Noise fixes...
!
!  Revision 1.10  2006/09/27 09:30:21  pchatela
!  Fixes, spectra calculation,
!  most importantly: moved the u_infty out, so it does not kill the dgammadt
!
!  Revision 1.9  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.8  2006/08/30 17:42:20  pchatela
!  Bugfixes in initial conditions
!  Added dump at t=0
!
!  Revision 1.7  2006/08/30 08:48:00  pchatela
!  Removed a useless comment
!
!  Revision 1.6  2006/08/24 11:28:22  menahel
!  added noise and netcdf restart capabilites
!
!  Revision 1.5  2006/08/11 15:32:26  pchatela
!  Fixed msgs
!
!  Revision 1.4  2006/08/11 12:43:51  pchatela
!  now compiles
!
!  Revision 1.3  2006/07/26 14:14:58  pchatela
!  Fixes
!
!  Revision 1.2  2006/07/26 14:09:32  pchatela
!  Clean ups and fixes
!
!  Revision 1.1  2006/07/26 13:52:48  pchatela
!  Initial insertion
!
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic INIT TRAILING VORTICES PHYSICS 0 =
!  initialize vorticity here
!------------------------------------------------------------------------------!
SUBROUTINE wvic_init_tvphysics_0

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft

  !----------------------------------------------------------------------------!
  ! interfaces
  INTERFACE
     SUBROUTINE wvic_alloc_field_s (vfield_up, info)
       USE module_wvic
       REAL (mk), DIMENSION (:, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_s
     SUBROUTINE wvic_alloc_field(vfield_up, ilda, info)
       USE module_wvic
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
  REAL(mk)                 :: tx_center, tz_vtx1, tz_vtx2
  REAL(mk)                 :: ty_vtx3, ty_vtx4, tx_vtx3, tx_vtx4
  REAL(mk)                 :: tx_vtx1, tx_vtx2, r_vtx3, r_vtx4
  REAL(mk)                 :: fac1, ra1, tx, ty, tz, noise, txp, typ, tzp, omega_x
  REAL(mk)                 :: fac2, ra2
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
  REAL(mk)                 :: ink_a0, ink_b0
  REAL(mk)                 :: tz_ink1, ty_ink1
  REAL(mk)                 :: tz_ink2, ty_ink2, pi
  REAL(mk)                 :: phixy, phixz, phiyx, phiyz, phizx, phizy
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
  REAL(mk)                     :: fac3,fac4,fac5,fac6,ldiv
  REAL(mk), DIMENSION(:,:,:,:), POINTER :: divergence, potential
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: divergence_v, potential_v
  
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid
  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! - generate noise (angle of dislocation and amplitute of dislocation)
  ! - calculate the position of the vortices using b0 (intervortex distance)
  ! - initialize the vortices using a0 (size of vortices)
  ! - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
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
     CALL ppm_write(rank,'wvic_init_tvphysics_0',msg,info)
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


  !-----------------------------------------------------
  !                          Restart from netcdf file
  !-----------------------------------------------------
  IF(netcdf_restart) THEN
     time  = netcdf_time
     itime = netcdf_itime
     dt    = netcdf_dt
     CALL wvic_netcdf2field(info)
     IF(rank.EQ.0) THEN
        WRITE(msg,*) 'restarted from netcdf, itime = ',netcdf_itime
        CALL ppm_write(rank,'wvic_init_physics_5',msg,info)
     END IF
     GOTO 1122
  END IF
  
  
  !----------------------------------------------------------------------------!
  ! compute the random amplitudes
  !----------------------------------------------------------------------------!
  IF (wvic_noise_basemode.EQ.0.0_mk) THEN
     wvic_noise_basemode = max_physg(3) - min_physg(3)
  END IF
  
  nk = wvic_noise_nmodes
  rnk = 1.0_mk/(2.0_mk*REAL(nk,mk))
  amp = wvic_noise_amp * rnk * trailvortex_b1
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
  ty_vtx1   = ty_center - 0.5_mk * trailvortex_b1
  ty_vtx2   = ty_center + 0.5_mk * trailvortex_b1
  ty_vtx3   = ty_center - 0.5_mk * trailvortex_b2
  ty_vtx4   = ty_center + 0.5_mk * trailvortex_b2
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))
  tx_vtx1   = tx_center
  tx_vtx2   = tx_center
  tx_vtx3   = tx_center + trailvortex_z12
  tx_vtx4   = tx_center + trailvortex_z12
  
  fac1      = trailvortex_gamma/(trailvortex_a1**2 * M_PI)
  ra1       = 1.0_mk/trailvortex_a1**2
  fac2      = trailvortex_r * trailvortex_gamma/(trailvortex_a2**2 * M_PI)
  ra2       = 1.0_mk/trailvortex_a2**2
  length    = max_physg - min_physg
  !----------------------------------------------------------------------------!

  
  !----------------------------------------------------------------------------!
  ! create vortices
  !----------------------------------------------------------------------------!
  max_vorticity = MAX(fac1,fac2)

  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)
        
        DO j=1,ndata(2,isubl)
           
           DO i=1,ndata(1,isubl)
              
              !tx = REAL(i+istart(1,isubl)-2,mk)*dx + min_physg(1) 
              !ty = REAL(j+istart(2,isubl)-2,mk)*dy + min_physg(2)
              !tz = REAL(k+istart(3,isubl)-2,mk)*dz + min_physg(3)
              !----------------------------------------------------------------!
              ! or better
              !----------------------------------------------------------------!
              tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
              ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
              tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

              !-----------------------------------------------------
              !  noise
              !-----------------------------------------------------
              omega_x = 2.0_mk*M_PI*tz/wvic_noise_basemode
              noise   = 0.0_mk
              DO kk=1,nk
                 noise = noise &
                      & + SIN(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*akz(kk))) &
                      & + COS(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*bkz(kk)))
              END DO
              txp = noise * amp
              noise   = 0.0_mk
              DO kk=1,nk
                 noise = noise &
                      & + SIN(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*aky(kk))) &
                      & + COS(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*bky(kk)))
              END DO
              typ = noise * amp
              tzp = 0.0_mk
              
              twp = 0.0_mk
              
              !----------------------------------------------------------------!
              ! compute vorticity - vtx1 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tx+txp-tx_vtx1)**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! left neighbor
              r_vtx1 = (ty+typ-ty_vtx1-length(2))**2 + (tx+txp-tx_vtx1)**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! right neighbor
              r_vtx1 = (ty+typ-ty_vtx1+length(2))**2 + (tx+txp-tx_vtx1)**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! bottom neighbor
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tx+txp-tx_vtx1-length(1))**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! top neighbor
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tx+txp-tx_vtx1+length(1))**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              !----------------------------------------------------------------!
              ! compute vorticity - vtx2 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx2 = (ty+txp-ty_vtx2)**2 + (tx+typ-tx_vtx2)**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! left neighbor
              r_vtx2 = (ty+txp-ty_vtx2-length(2))**2 + (tx+typ-tx_vtx2)**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! right neighbor
              r_vtx2 = (ty+txp-ty_vtx2+length(2))**2 + (tx+typ-tx_vtx2)**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! bottom neighbor
              r_vtx2 = (ty+txp-ty_vtx2)**2 + (tx+typ-tx_vtx2-length(1))**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! top neighbor
              r_vtx2 = (ty+txp-ty_vtx2)**2 + (tx+typ-tx_vtx2+length(1))**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              !----------------------------------------------------------------!
              ! compute vorticity - vtx3 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx3 = (ty-txp-ty_vtx3)**2 + (tx-typ-tx_vtx3)**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! left neighbor
              r_vtx3 = (ty-txp-ty_vtx3-length(2))**2 + (tx-typ-tx_vtx3)**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! right neighbor
              r_vtx3 = (ty-txp-ty_vtx3+length(2))**2 + (tx-typ-tx_vtx3)**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! bottom neighbor
              r_vtx3 = (ty-txp-ty_vtx3)**2 + (tx-typ-tx_vtx3-length(1))**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! top neighbor
              r_vtx3 = (ty-txp-ty_vtx3)**2 + (tx-typ-tx_vtx3+length(1))**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              !----------------------------------------------------------------!
              ! compute vorticity - vtx4 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx4 = (ty-typ-ty_vtx4)**2 + (tx-txp-tx_vtx4)**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! left neighbor
              r_vtx4 = (ty-typ-ty_vtx4-length(2))**2 + (tx-txp-tx_vtx4)**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! right neighbor
              r_vtx4 = (ty-typ-ty_vtx4+length(2))**2 + (tx-txp-tx_vtx4)**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! bottom neighbor
              r_vtx4 = (ty-typ-ty_vtx4)**2 + (tx-txp-tx_vtx4-length(1))**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! top neighbor
              r_vtx4 = (ty-typ-ty_vtx4)**2 + (tx-txp-tx_vtx4+length(1))**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)

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

  IF (trailvortex.EQV..FALSE.) THEN
  !-----------------------------------------------------
  !  compute velocity from this
  !-----------------------------------------------------
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_up(1,i,j,k,isub) = field_wp(1,i,j,k,isub) 
              field_up(2,i,j,k,isub) = field_wp(2,i,j,k,isub) 
              field_up(3,i,j,k,isub) = field_wp(3,i,j,k,isub) 
           END DO
        END DO
     END DO
  END DO
  
  CALL ppm_fft_velocities(field_up,&
       & mesh_id,topo_id,t_topoid,fmesh_id,ghostsize,info)
       
  !-----------------------------------------------------
  !  take curl of velocity
  !-----------------------------------------------------
  CALL wvic_ghost(wvic_prm_velocity,info)
  fac1 = 8.0_mk / dx / 12.0_mk
  fac2 = 8.0_mk / dy / 12.0_mk
  fac3 = 8.0_mk / dz / 12.0_mk
  fac4 = 1.0_mk / dx / 12.0_mk
  fac5 = 1.0_mk / dy / 12.0_mk
  fac6 = 1.0_mk / dz / 12.0_mk
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              ! d u / d y
              phixy= fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))&
                   &-fac5*(field_up(1,i,j+2,k,isub)-field_up(1,i,j-2,k,isub))
              ! d u / d z
              phixz= fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))&
                   &-fac6*(field_up(1,i,j,k+2,isub)-field_up(1,i,j,k-2,isub))
              ! d v / d x
              phiyx= fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))&
                   &-fac4*(field_up(2,i+2,j,k,isub)-field_up(2,i-2,j,k,isub))
              ! d v / d z
              phiyz= fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))&
                   &-fac6*(field_up(2,i,j,k+2,isub)-field_up(2,i,j,k-2,isub))
              ! d w / d x
              phizx= fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))&
                   &-fac4*(field_up(3,i+2,j,k,isub)-field_up(3,i-2,j,k,isub))
              ! d w / d y
              phizy= fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))&
                   &-fac5*(field_up(3,i,j+2,k,isub)-field_up(3,i,j-2,k,isub))
              
              field_wp(1,i,j,k,isub) = phizy - phiyz
              field_wp(2,i,j,k,isub) = phixz - phizx
              field_wp(3,i,j,k,isub) = phiyx - phixy
           END DO
        END DO
     END DO
  END DO
  
  END IF
  
  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
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
1122 CONTINUE

  !----------------------------------------------------------------------------!
  ! create particles with style
  !----------------------------------------------------------------------------!
  CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,topo_id,mesh_id,&
       & (/cutoff,HUGE(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
  WRITE(msg,*) ' created ',np,' particles'
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_tvphysics_0',msg,info)
END SUBROUTINE wvic_init_tvphysics_0
