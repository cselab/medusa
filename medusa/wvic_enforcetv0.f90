!------------------------------------------------------------------------------!
!* filename: wvic_enforcetv0                                                  *!
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
!  $Log: wvic_enforcetv0.F,v $
!  Revision 1.7  2006/10/16 15:38:46  pchatela
!  Added divergence free (analytical) wvic_init_tvphysics_1.F
!
!  Revision 1.6  2006/10/04 11:16:38  pchatela
!  Added a fundamental mode to the noise parameters.
!  If not given, the length of the domain is used.
!
!  Revision 1.5  2006/09/28 15:35:52  pchatela
!  Noise fixes...
!
!  Revision 1.4  2006/09/27 09:30:21  pchatela
!  Fixes, spectra calculation,
!  most importantly: moved the u_infty out, so it does not kill the dgammadt
!
!  Revision 1.3  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.2  2006/08/31 08:00:50  pchatela
!  Fixed bugs similar to wvic_init_tvphysics_0.F
!  Non-empty message now
!
!  Revision 1.1  2006/08/11 13:36:04  pchatela
!  Initial insertion
!
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic ENFORCE TRAILING VORTICES CONDITION 0
!  initialize vorticity here
!------------------------------------------------------------------------------!
SUBROUTINE wvic_enforcetv0

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map
  
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
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: ig,jg,kg,i,j,k,isub,isubl,kmax,kl
  REAL(mk), PARAMETER :: onethird         = 0.3333333333333332_mk
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
  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! - generate noise (angle of dislocation and amplitute of dislocation)
  ! - calculate the position of the vortices using b0 (intervortex distance)
  ! - initialize the vortices using a0 (size of vortices)
  ! - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!


  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_enforcetv0',msg,info)
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
  max_vorticity = fac1

  DO isub=1,nsublist
     isubl = isublist(isub)
	 IF (min_sub(3,isubl).LT.(min_physg(3)+0.5*dz)) THEN
        DO k=1-ghostsize(3),1
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)
			  tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
              ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
			  tz = min_sub(3,isubl) + REAL(k-1,mk)*dz
              
			  !-----------------------------------------------------
              !  noise
              !-----------------------------------------------------
              omega_x = 2.0_mk*M_PI*(tz-u_infty(3)*time)/wvic_noise_basemode
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
			  IF(nk.EQ.1) typ = 0.0_mk
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
     END IF
     IF (max_sub(3,isubl).GT.(max_physg(3)-0.5*dz)) THEN
        k = ndata(3,isubl)
        DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
           DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
              field_wp(1,i,j,k,isub) = onethird*(4.0_mk*field_wp(1,i,j,k-1,isub) &
    &                                               - field_wp(1,i,j,k-2,isub))
              field_wp(2,i,j,k,isub) = onethird*(4.0_mk*field_wp(2,i,j,k-1,isub) &
    &                                               - field_wp(2,i,j,k-2,isub))
              field_wp(3,i,j,k,isub) = onethird*(4.0_mk*field_wp(3,i,j,k-1,isub) &
    &                                               - field_wp(3,i,j,k-2,isub))
           END DO
        END DO
     END IF
  END DO
  !---- Vorticity has been initialized on the upstream boundary of the field.
  WRITE(msg,*) ' enforced the upstream dirichlet condition '
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_enforcetv0',msg,info)


END SUBROUTINE wvic_enforcetv0
