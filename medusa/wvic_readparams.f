!                                                                              !
!                               %m           '@p                               !
!                               #@F           #@                    9@_        !
!           ___g_              #@*           q@*      ___       _    #@        !
!  @__gggE###""@@p           _@#             @L___gg@####       '@_  #|        !
!  '##""      g@#"          g@"             #F """d@             ##  #|        !
!            ##"          _#F @            g*      @             @   #      _F !
!      #_  _@*          _gW   #L          #       #@            q@   #    _gF  !
!       7@gW          _#*     @L                  @*            @*   #   _@^   !
!         @@                  @L                 #F            gF    #_g@#     !
!         '@|                 @L                g#            _F     #@#"      !
!          7                  @k               gF            g"       ^        !
!                             #"              #`                               !
!                                                                              !
!-------------------------------------------------------------------------------
!* filename: pmira_readparams                                                  *!
!* project : pmira                                                             *!
!* purpose : read parameters                                                  *!
!*         :                                                                  *!
!* author  : Michael Bergdorf, modified Ivo Sbalzarinis code                  *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Thu Jun  3 13:23:36 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!-------------------------------------------------------------------------------
! $Log: wvic_readparams.F,v $
! Revision 1.15  2006/10/25 11:35:50  menahel
! re-enabled .rst dumping and added NRESTART parameter
!
! Revision 1.14  2006/10/24 16:52:53  pchatela
! Added KEEPWCENTERED option if one wants to keep vorticity centered in domain
!
! Revision 1.13  2006/10/23 14:28:03  pchatela
! Added parameters for hyperviscosity model
! Fixed double Laplacian
!
! Revision 1.12  2006/10/23 08:19:11  pchatela
! LES models
! Bugfixes in KE spectra and factor for Parseval identity
! Removed the reset of noise in init_tvphysics_0 and_1
!
! Revision 1.11  2006/10/17 14:17:40  pchatela
! Added Ctrl parameters to handle LES models
!
! Revision 1.10  2006/10/04 11:16:38  pchatela
! Added a fundamental mode to the noise parameters.
! If not given, the length of the domain is used.
!
! Revision 1.9  2006/10/03 16:11:58  pchatela
! Added spectra in the z direction (misnamed kx spectrum...)
! Added spatial diagnostics, like kinetic energy, enstrophy, circulation
! as functions of z, dumped at the frequency ndump
!
! Revision 1.8  2006/09/27 09:30:21  pchatela
! Fixes, spectra calculation,
! most importantly: moved the u_infty out, so it does not kill the dgammadt
!
! Revision 1.7  2006/09/11 14:57:27  pchatela
! Fixed velocity computation with odd/even symmetries
! Added flag for adaptive time step
!
! Revision 1.6  2006/09/01 15:47:03  pchatela
! Added option to dump velocity and streamfunction fields
!
! Revision 1.5  2006/08/30 17:42:20  pchatela
! Bugfixes in initial conditions
! Added dump at t=0
!
! Revision 1.4  2006/08/22 19:41:31  menahel
! corrected ...VORTICES_GAMMA to ...VORTICES_PRIMARY_GAMMA
!
! Revision 1.3  2006/07/26 14:09:32  pchatela
! Clean ups and fixes
!
! Revision 1.2  2006/07/26 07:51:26  pchatela
! Added boolean trailvortex
! Added periodic bcs with reset of particle strengths
!
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.3  2005/12/10 20:28:39  michaebe
! added ndump to ctrl parameters
!
! Revision 1.2  2005/11/21 17:36:39  michaebe
! major session. now let''s debug
!
! Revision 1.1  2005/09/28 11:40:41  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------

SUBROUTINE wvic_readparams(ctrlfile, info)

  USE  module_wvic
  IMPLICIT NONE
  

  INTEGER, INTENT(inout)                :: info
  CHARACTER(len=256), INTENT(in)        :: ctrlfile

  !-----------------------------------------------------------------------------
  ! Declaration of local variables
  !-----------------------------------------------------------------------------
  
  INTEGER                               :: i,j
  INTEGER                               :: idx,i1,i2
  CHARACTER(LEN=256)                    :: cbuf
  INTEGER                               :: iUnit
  INTEGER                               :: ILEN,ios
  INTEGER                               :: ilenctrl
  INTEGER                               :: ibc,iline
  CHARACTER(LEN=256)                    :: cvalue,carg
  CHARACTER(LEN=256)                    :: bcloc,tvalue
  LOGICAL                               :: lExist
  
  !-----------------------------------------------------------------------------
  ! Definition of file unit
  !-----------------------------------------------------------------------------
  
  iUnit = 20
  info  = 0
  !-----------------------------------------------------------------------------
  ! Check that the parameter file exists
  !-----------------------------------------------------------------------------
  
  ilenctrl = LEN_TRIM(ctrlfile)         ! length of filename string
  INQUIRE(FILE=ctrlfile(1:ilenctrl), EXIST=lExist)
  IF(.NOT. lExist) THEN
     WRITE(*,'(2A)')'No such file: ',ctrlfile(1:ilenctrl)
     Info = 1
     GOTO 9999
  END IF
  
  !-----------------------------------------------------------------------------
  ! open the file
  !-----------------------------------------------------------------------------
  
  OPEN(iUnit, FILE=ctrlfile(1:ilenctrl), IOSTAT=ios, ACTION='READ')
  IF(ios .NE. 0) THEN
     WRITE(*,'(2A)')'Failed to open file: ',ctrlfile(1:ilenctrl)
     Info = 1
     GOTO 9999
  END IF
  
  !-----------------------------------------------------------------------------
  ! scan file
  !-----------------------------------------------------------------------------
  
  iline = 0
  DO 
     iline = iline + 1        ! increment line 
     READ(iUnit,'(A)',END=100,ERR=200) cbuf 
     ILEN = LEN_TRIM(cbuf)
     
     !--------------------------------------------------------------------------
     !  Skip comment or empty lines 
     !--------------------------------------------------------------------------
     IF(ILEN .GT. 0 .AND. cbuf(1:1) .NE. '#') THEN
        !-----------------------------------------------------------------------
        !  Remove all spaces
        !-----------------------------------------------------------------------
        j = 0
        DO i=1,ILEN
           IF(cbuf(i:i) .NE. ' ') THEN
              j = j + 1
              cbuf(j:j) = cbuf(i:i)
           END IF
        END DO
        ILEN = j     ! update length of string
        
        !-----------------------------------------------------------------------
        !  Find position of =
        !-----------------------------------------------------------------------
        idx = INDEX(cbuf,'=')
        
        !-----------------------------------------------------------------------
        !  Exit if = is missing
        !-----------------------------------------------------------------------
        IF(idx .LT. 0) THEN
           WRITE(*,'(A,I5)')'Incorrect line: ',iline
           Info = 1
           GOTO 9999
        END IF
        
        !-----------------------------------------------------------------------
        !  Get argument and value
        !-----------------------------------------------------------------------
        carg   = ADJUSTL(cbuf(1:idx-1))
        cvalue = ADJUSTL(cbuf(idx+1:ILEN))
        
        !-----------------------------------------------------------------------
        !  Convert to upper case
        !-----------------------------------------------------------------------
        CALL UpperCase(carg,idx-1)

        !-----------------------------------------------------------------------
        !  Parse and read input data
        !-----------------------------------------------------------------------
        IF(carg(1:(idx-1)) .EQ. 'RUNTAG') THEN
           READ(cvalue,'(A)',IOSTAT=ios,ERR=200) runtag
           iruntag = LEN_TRIM(runtag)
        ELSEIF (carg(1:(idx-1)).EQ.'NX') THEN
           READ(cvalue,*,iostat=ios,err=200) nx
        ELSEIF (carg(1:(idx-1)).EQ.'MIN_PHYS') THEN
           READ(cvalue,*,iostat=ios,err=200) min_physg
         ELSEIF (carg(1:(idx-1)).EQ.'MAX_PHYS') THEN
           READ(cvalue,*,iostat=ios,err=200) max_physg
         ELSEIF (carg(1:(idx-1)).EQ.'GHOSTSIZE') THEN
            READ(cvalue,*,iostat=ios,err=200) ghostsize
         ELSEIF (carg(1:(idx-1)).EQ.'TOPOLOGY') THEN
            READ(cvalue,*,iostat=ios,err=200) ndims
         ELSEIF (carg(1:(idx-1)).EQ.'DT') THEN
            READ(cvalue,*,iostat=ios,err=200) dt
         ELSEIF (carg(1:(idx-1)).EQ.'DT_ADAPT') THEN
            READ(cvalue,*,iostat=ios,err=200) dt_adapt
         ELSEIF (carg(1:(idx-1)).EQ.'DTMAX') THEN
            READ(cvalue,*,iostat=ios,err=200) dt_max
         ELSEIF (carg(1:(idx-1)).EQ.'NDUMP') THEN
            READ(cvalue,*,iostat=ios,err=200) ndump
         ELSEIF (carg(1:(idx-1)).EQ.'NRESTART') THEN
            READ(cvalue,*,iostat=ios,err=200) nrestart
         ELSEIF (carg(1:(idx-1)).EQ.'NSOLENOID') THEN
            READ(cvalue,*,iostat=ios,err=200) nkilldiv
         ELSEIF (carg(1:(idx-1)).EQ.'DUMPVEL') THEN
            READ(cvalue,*,iostat=ios,err=200) dumpvel
         ELSEIF (carg(1:(idx-1)).EQ.'DUMPPSI') THEN
            READ(cvalue,*,iostat=ios,err=200) dumppsi
         ELSEIF (carg(1:(idx-1)).EQ.'DUMPKESPEC') THEN
            READ(cvalue,*,iostat=ios,err=200) dumpkespec
         ELSEIF (carg(1:(idx-1)).EQ.'DUMPKESPECKX') THEN
            READ(cvalue,*,iostat=ios,err=200) dumpkespeckx
         ELSEIF (carg(1:(idx-1)).EQ.'KEEPWCENTERED') THEN
            READ(cvalue,*,iostat=ios,err=200) keepwcentered
         ELSEIF (carg(1:(idx-1)).EQ.'TEND') THEN
            READ(cvalue,*,iostat=ios,err=200) tend
         ELSEIF (carg(1:(idx-1)).EQ.'NU') THEN
            READ(cvalue,*,iostat=ios,err=200) nu
         ELSEIF (carg(1:(idx-1)).EQ.'TARGET_RE') THEN
            READ(cvalue,*,iostat=ios,err=200) target_re
         ELSEIF (carg(1:(idx-1)).EQ.'SC') THEN
            READ(cvalue,*,iostat=ios,err=200) sc
         ELSEIF (carg(1:(idx-1)).EQ.'MULTIGRID') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_multigrid
         ELSEIF (carg(1:(idx-1)).EQ.'LES') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_les
            IF (wvic_les) THEN 
                wvic_les_model = 1
                les_tdm_C_comp = 0.5_mk
                les_tdm_C_dil  = 0.0_mk
            END IF
         ELSEIF (carg(1:(idx-1)).EQ.'LES_MODEL') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_les_model
            IF (wvic_les_model.GT.0) wvic_les = .TRUE.
         ELSEIF (carg(1:(idx-1)).EQ.'LES_TDM_C_COMP') THEN
            READ(cvalue,*,iostat=ios,err=200) les_tdm_C_comp
            wvic_les_model = 1
            wvic_les = .TRUE.
         ELSEIF (carg(1:(idx-1)).EQ.'LES_TDM_C_DIL') THEN
            READ(cvalue,*,iostat=ios,err=200) les_tdm_C_dil
            wvic_les_model = 1
            wvic_les = .TRUE.
         ELSEIF (carg(1:(idx-1)).EQ.'LES_TDMCLIPPED_C') THEN
            READ(cvalue,*,iostat=ios,err=200) les_tdmclipped_C
            wvic_les_model = 2
            wvic_les = .TRUE.
         ELSEIF (carg(1:(idx-1)).EQ.'LES_HYPERVISC_C') THEN
            READ(cvalue,*,iostat=ios,err=200) les_hypervisc_C
            wvic_les_model = 3
            wvic_les = .TRUE.
         ELSEIF (carg(1:(idx-1)).EQ.'LES_HYPERVISC_TZERO') THEN
            READ(cvalue,*,iostat=ios,err=200) les_hypervisc_tzero
            wvic_les_model = 3
            wvic_les = .TRUE.
         ELSEIF (carg(1:(idx-1)).EQ.'RENORM') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_renormalize
         ELSEIF (carg(1:(idx-1)).EQ.'MUSCLE') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_muscle
         ELSEIF (carg(1:(idx-1)).EQ.'SCHEME') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_dgamma_scheme
         ELSEIF (carg(1:(idx-1)).EQ.'W_SCHEME') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_dgamma_scheme
         ELSEIF (carg(1:(idx-1)).EQ.'U_SCHEME') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_compvel_scheme
         ELSEIF (carg(1:(idx-1)).EQ.'UNBOUNDEDXY') THEN
            READ(cvalue,*,iostat=ios,err=200) unboundedxy
         ELSEIF (carg(1:(idx-1)).EQ.'TALK') THEN
            READ(cvalue,*,iostat=ios,err=200) verbose
         ELSEIF (carg(1:(idx-1)).EQ.'MAXLEV') THEN
            READ(cvalue,*,iostat=ios,err=200) maxlev
         ELSEIF (carg(1:(idx-1)).EQ.'NOISE_NMODES') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_noise_nmodes
         ELSEIF (carg(1:(idx-1)).EQ.'TUBE_RADIUS') THEN
            READ(cvalue,*,iostat=ios,err=200) tube_radius
         ELSEIF (carg(1:(idx-1)).EQ.'NOISE_AMP') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_noise_amp
         ELSEIF (carg(1:(idx-1)).EQ.'NOISE_FUNDMODE') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_noise_basemode
         ELSEIF (carg(1:(idx-1)).EQ.'RE') THEN
            READ(cvalue,*,iostat=ios,err=200) wvic_re
         ELSEIF (carg(1:(idx-1)).EQ.'U_INFTY') THEN
            READ(cvalue,*,iostat=ios,err=200) u_infty
         ELSEIF (carg(1:(idx-1)).EQ.'WALLS') THEN
            READ(cvalue,*,iostat=ios,err=200) wbcdef
		 ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex
         ELSEIF (carg(1:(idx-1)).EQ.'FLOW_CASE') THEN
            READ(cvalue,*,iostat=ios,err=200) flow_case
         ELSEIF (carg(1:(idx-1)).EQ.'DIE_AFTER') THEN
            READ(cvalue,*,iostat=ios,err=200) tot
            !-----------------------------------------------------
            !  netcdf restart
            !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'NETCDF_RESTART') THEN
            READ(cvalue,*,iostat=ios,err=200) netcdf_restart
         ELSEIF (carg(1:(idx-1)).EQ.'NETCDF_ITIME') THEN
            READ(cvalue,*,iostat=ios,err=200) netcdf_itime
         ELSEIF (carg(1:(idx-1)).EQ.'NETCDF_TIME') THEN
            READ(cvalue,*,iostat=ios,err=200) netcdf_time
         ELSEIF (carg(1:(idx-1)).EQ.'NETCDF_DT') THEN
            READ(cvalue,*,iostat=ios,err=200) netcdf_dt
            !-----------------------------------------------------
            !  piston parameters
            !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'PISTON_RADIUS') THEN
            READ(cvalue,*,iostat=ios,err=200) piston_radius
         ELSEIF (carg(1:(idx-1)).EQ.'PISTON_RE_RATIO') THEN
            READ(cvalue,*,iostat=ios,err=200) piston_re_ratio
         ELSEIF (carg(1:(idx-1)).EQ.'PISTON_VORT_RATIO') THEN
            READ(cvalue,*,iostat=ios,err=200) piston_vort_ratio
            !-----------------------------------------------------
            !  trailing vortices parameters
            !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_PRIMARY_SPAN') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_b1
         ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_SECONDARY_SPAN') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_b2
		 ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_PRIMARY_SIGMA') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_a1
         ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_SECONDARY_SIGMA') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_a2
         ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_PRIMARY_VERTDIFF') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_z12
		 ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_PRIMARY_GAMMA') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_gamma	
		 ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_GAMMA_RATIO') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_r
		 ELSEIF (carg(1:(idx-1)).EQ.'TRAILING_VORTICES_PRIMARY_GAMMA_RATIO') THEN
            READ(cvalue,*,iostat=ios,err=200) trailvortex_r
            !-----------------------------------------------------
            !  helical vortices parameters
            !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_N') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_n
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_R1') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_r1
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_R2') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_r2
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_A1') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_a1
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_A2') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_a2
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_PERIOD') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_period
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_GAMMA') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_gamma
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_GAMMA_RATIO') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_ratio !JTR
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_VORTICES_NVORTEX') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_nvortex 
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_PREC_ISEC1') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_precision_isec1
         ELSEIF (carg(1:(idx-1)).EQ.'HELICAL_PREC_ISEC2') THEN
            READ(cvalue,*,iostat=ios,err=200) helvortex_precision_isec2
         !-----------------------------------------------------
         !  Geometry - plate / block (Poiseuille)
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'BLOCK_WX') THEN
            READ(cvalue,*,iostat=ios,err=200) block_wx
         ELSEIF (carg(1:(idx-1)).EQ.'BLOCK_WY') THEN
            READ(cvalue,*,iostat=ios,err=200) block_wy
         ELSEIF (carg(1:(idx-1)).EQ.'BLOCK_WZ') THEN
            READ(cvalue,*,iostat=ios,err=200) block_wz
         ELSEIF (carg(1:(idx-1)).EQ.'PLATE_THICKNESS05INIT') THEN
            READ(cvalue,*,iostat=ios,err=200) thickness05init
         !-----------------------------------------------------
         !  Geometry - sphere
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'SPHERE_RADIUS') THEN
            READ(cvalue,*,iostat=ios,err=200) sphere_radius
         ELSEIF (carg(1:(idx-1)).EQ.'SPHERE_CENTER') THEN !JTR temp
            READ(cvalue,*,iostat=ios,err=200) object_offset !JTR temp
         ELSEIF (carg(1:(idx-1)).EQ.'OBJECT_OFFSET') THEN
            READ(cvalue,*,iostat=ios,err=200) object_offset
         ELSEIF (carg(1:(idx-1)).EQ.'LEVELSET_BAND') THEN  !JTR temp
            READ(cvalue,*,iostat=ios,err=200) stepfunction_band  !JTR temp
         !-----------------------------------------------------
         !  Geometry - cylinder(s)
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'CYLINDER_RADIUS') THEN
            READ(cvalue,*,iostat=ios,err=200) cylinder_radius
         ELSEIF (carg(1:(idx-1)).EQ.'CYLINDER_NO') THEN
            READ(cvalue,*,iostat=ios,err=200) cylinder_no
         !-----------------------------------------------------
         !  Step function specific
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'STEPFUNCTION_BAND') THEN
            READ(cvalue,*,iostat=ios,err=200) stepfunction_band
         ELSEIF (carg(1:(idx-1)).EQ.'LEVELSET_LAMBDA') THEN !JTR temp
            READ(cvalue,*,iostat=ios,err=200) penalization_lambda !JTR temp
         ELSEIF (carg(1:(idx-1)).EQ.'STEP_FUNCTION') THEN
            READ(cvalue,*,iostat=ios,err=200) step_function  
         ELSEIF (carg(1:(idx-1)).EQ.'STEP1_INTERVAL') THEN
            READ(cvalue,*,iostat=ios,err=200) step1_interval 
         ELSEIF (carg(1:(idx-1)).EQ.'STEP1_LINEARFRACTION') THEN
            READ(cvalue,*,iostat=ios,err=200) step1_linearfraction
         ELSEIF (carg(1:(idx-1)).EQ.'STEP1_OFFSET') THEN
            READ(cvalue,*,iostat=ios,err=200) step1_offset
         !-----------------------------------------------------
         !  Penalization specific
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'PENALIZATION_LAMBDA') THEN
            READ(cvalue,*,iostat=ios,err=200) penalization_lambda
         ELSEIF (carg(1:(idx-1)).EQ.'PENALIZATION_IMPLICIT') THEN
            READ(cvalue,*,iostat=ios,err=200) penalization_implicit
         ELSEIF (carg(1:(idx-1)).EQ.'PENALIZATION_CLEARRHS') THEN
            READ(cvalue,*,iostat=ios,err=200) penalization_clearrhs
         ELSEIF (carg(1:(idx-1)).EQ.'PENALIZATION_ONETERM') THEN
            READ(cvalue,*,iostat=ios,err=200) penalization_oneterm
         !-----------------------------------------------------
         !  Ramping, force calculation, misc
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'DT_ADAPT_FRACTION') THEN
            READ(cvalue,*,iostat=ios,err=200) dt_adapt_fraction
         ELSEIF (carg(1:(idx-1)).EQ.'DT_RAMPSTEPS') THEN
            READ(cvalue,*,iostat=ios,err=200) dt_rampsteps
         ELSEIF (carg(1:(idx-1)).EQ.'DT_RAMPSTART') THEN
            READ(cvalue,*,iostat=ios,err=200) dt_rampstart
         ELSEIF (carg(1:(idx-1)).EQ.'ONSET_RAMPTIME') THEN
            READ(cvalue,*,iostat=ios,err=200) onset_ramptime
         ELSEIF (carg(1:(idx-1)).EQ.'NFORCES') THEN
            READ(cvalue,*,iostat=ios,err=200) nforces
            ALLOCATE(forcecv(6*nforces))
         ELSEIF (carg(1:(idx-1)).EQ.'FORCECV') THEN
            READ(cvalue,*,iostat=ios,err=200) forcecv
         !-----------------------------------------------------
         !  Movement specific
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'U_SOLID') THEN
            READ(cvalue,*,iostat=ios,err=200) u_solid
         ELSEIF (carg(1:(idx-1)).EQ.'HARMONIC_AMPLITUDE') THEN
            READ(cvalue,*,iostat=ios,err=200) harmonic_amplitude
         ELSEIF (carg(1:(idx-1)).EQ.'HARMONIC_PERIOD') THEN
            READ(cvalue,*,iostat=ios,err=200) harmonic_period
         ELSEIF (carg(1:(idx-1)).EQ.'HARMONIC_PHASE') THEN
            READ(cvalue,*,iostat=ios,err=200) harmonic_phase
         !-----------------------------------------------------
         !  Medusa               
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'MEDUSA_FILE') THEN
            READ(cvalue,*,iostat=ios,err=200) medusa_file
         ELSEIF (carg(1:(idx-1)).EQ.'MEDUSA_SCALE') THEN
            READ(cvalue,*,iostat=ios,err=200) medusa_scale
         ELSEIF (carg(1:(idx-1)).EQ.'MEDUSA_MASS_ORIGINAL') THEN
            READ(cvalue,*,iostat=ios,err=200) medu_mass_orig
         ELSEIF (carg(1:(idx-1)).EQ.'MEDUSA_PERIOD') THEN
            READ(cvalue,*,iostat=ios,err=200) medusa_period
         ELSEIF (carg(1:(idx-1)).EQ.'MEDUSA_FRAMES') THEN
            READ(cvalue,*,iostat=ios,err=200) medusa_frames
         !-----------------------------------------------------
         !  STL input            
         !-----------------------------------------------------
         ELSEIF (carg(1:(idx-1)).EQ.'STL_FILE') THEN
            READ(cvalue,*,iostat=ios,err=200) stlfile
         ELSEIF (carg(1:(idx-1)).EQ.'STL_TRANSLATE') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_translate
         ELSEIF (carg(1:(idx-1)).EQ.'STL_SCALE') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_scale
         ELSEIF (carg(1:(idx-1)).EQ.'STL_CHECK_BOUNDING') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_check_bounding
         ELSEIF (carg(1:(idx-1)).EQ.'STL_CHECK_INTERSECTIONS') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_check_intersections
         ELSEIF (carg(1:(idx-1)).EQ.'STL_DOUBLE_CHECK') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_double_check
         ELSEIF (carg(1:(idx-1)).EQ.'STL_NONVERBOSE') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_nonverbose
         ELSEIF (carg(1:(idx-1)).EQ.'STL_MIN_PHYS_INSIDE') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_min_phys_inside
         ELSEIF (carg(1:(idx-1)).EQ.'STL_MOLLIFY') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_mollify
         ELSEIF (carg(1:(idx-1)).EQ.'STL_STENCIL') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_stencil
         ELSEIF (carg(1:(idx-1)).EQ.'STL_SUBDIVIDE') THEN
            READ(cvalue,*,iostat=ios,err=200) stl_subdivide
         END IF
        
     END IF
  END DO
  
  !-----------------------------------------------------------------------------
  !  Something went wrong if we got here
  !-----------------------------------------------------------------------------
200 CONTINUE
  WRITE(*,'(A,I5,2A)') 'Error reading line: ',iline,                      &
       &                 ' of file: ',ctrlfile(1:ilenctrl)
  ILEN = LEN_TRIM(cbuf)
  WRITE(*,'(A)') cbuf(1:ILEN)
  Info = 1
  GOTO 9999
  
  !-----------------------------------------------------------------------------
  !  End of file
  !-----------------------------------------------------------------------------
100 Info = 0
  
  !-----------------------------------------------------------------------------
  !  Close file
  !-----------------------------------------------------------------------------
  CLOSE(iUnit)
  
  !-----------------------------------------------------------------------------
  !  Return 
  !-----------------------------------------------------------------------------
9999 CONTINUE
  
END SUBROUTINE wvic_readparams




