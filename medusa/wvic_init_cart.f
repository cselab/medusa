!------------------------------------------------------------------------------!
!* filename: wvic_init_cart                                                        *!
!* project : ppm                                                              *!
!* purpose : initialize the wvic tailored for                                 *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 13:14:26 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_init_cart.F,v $
!  Revision 1.18  2006/10/25 17:55:47  pchatela
!  Bugfix: divergence-free shaking of tubes
!
!  Revision 1.17  2006/10/25 11:35:50  menahel
!  re-enabled .rst dumping and added NRESTART parameter
!
!  Revision 1.16  2006/10/24 16:52:53  pchatela
!  Added KEEPWCENTERED option if one wants to keep vorticity centered in domain
!
!  Revision 1.15  2006/10/23 14:28:03  pchatela
!  Added parameters for hyperviscosity model
!  Fixed double Laplacian
!
!  Revision 1.14  2006/10/23 08:19:11  pchatela
!  LES models
!  Bugfixes in KE spectra and factor for Parseval identity
!  Removed the reset of noise in init_tvphysics_0 and_1
!
!  Revision 1.13  2006/10/17 14:17:40  pchatela
!  Added Ctrl parameters to handle LES models
!
!  Revision 1.12  2006/10/16 15:38:46  pchatela
!  Added divergence free (analytical) wvic_init_tvphysics_1.F
!
!  Revision 1.11  2006/10/04 11:16:38  pchatela
!  Added a fundamental mode to the noise parameters.
!  If not given, the length of the domain is used.
!
!  Revision 1.10  2006/10/03 16:11:58  pchatela
!  Added spectra in the z direction (misnamed kx spectrum...)
!  Added spatial diagnostics, like kinetic energy, enstrophy, circulation
!  as functions of z, dumped at the frequency ndump
!
!  Revision 1.9  2006/09/27 09:30:21  pchatela
!  Fixes, spectra calculation,
!  most importantly: moved the u_infty out, so it does not kill the dgammadt
!
!  Revision 1.8  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.7  2006/09/01 15:47:03  pchatela
!  Added option to dump velocity and streamfunction fields
!
!  Revision 1.6  2006/08/24 12:04:39  pchatela
!  Nuked the mg_init stuff
!
!  Revision 1.5  2006/08/24 11:51:59  pchatela
!  Added simple test for single prec domain decomposition problem
!
!  Revision 1.4  2006/08/24 09:26:01  menahel
!  removed logfile stuff, cleanup.
!
!  Revision 1.3  2006/07/26 14:09:32  pchatela
!  Clean ups and fixes
!
!  Revision 1.2  2006/07/26 07:51:26  pchatela
!  Added boolean trailvortex
!  Added periodic bcs with reset of particle strengths
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.7  2005/12/10 20:28:56  michaebe
!  broadcast of ndump
!
!  Revision 1.6  2005/11/21 17:36:38  michaebe
!  major session. now let''s debug
!
!  Revision 1.5  2005/11/11 17:15:23  michaebe
!  adapted field alloc calls
!
!  Revision 1.4  2005/11/11 14:04:24  michaebe
!  clean up, additions, comments
!
!  Revision 1.3  2005/11/11 09:57:18  michaebe
!  new
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic INIT CROW =
!------------------------------------------------------------------------------!
SUBROUTINE wvic_init_cart(ctrlfile, info)

  USE ppm_module_data
  USE ppm_module_init
  USE ppm_module_mktopo
  USE ppm_module_set_proc_speed
  USE ppm_module_estimate_proc_speed
  USE ppm_module_write
  USE ppm_module_mg
  USE module_wvic
!JTMP
!  USE wvic_module_io
!  USE netcdf
!  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE wvic_alloc_field (vfield_up, ilda, info)
       USE module_wvic
       IMPLICIT NONE
       REAL (mk), DIMENSION (:, :, :, :, :), POINTER :: vfield_up
       INTEGER                          , INTENT(in) :: ilda
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field

     SUBROUTINE wvic_alloc_field_s (vfield_up, info)
       USE module_wvic
       IMPLICIT NONE
       REAL (mk), DIMENSION (:,  :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_s
  END INTERFACE
  
  !----------------------------------------------------------------------------!
  !  arguments
  !----------------------------------------------------------------------------!
  CHARACTER(len=256), INTENT(in) :: ctrlfile      ! name of the control file
  INTEGER, INTENT(inout)         :: info
  !----------------------------------------------------------------------------!
  !  localities
  !----------------------------------------------------------------------------!
  CHARACTER(len=256) :: logfile, dumpfile, msg
  INTEGER            :: assigning, decomposition
  INTEGER            :: prec    !
  INTEGER            :: nsubs   !
  INTEGER            :: maptype
  INTEGER            :: tol, istat, isub
  INTEGER            :: i,j,kp,k,nxl,nxu
  INTEGER            :: iwvic_multigrid, iverbose,isubl
  REAL(mk)           :: r,tx,ty,tz, SIGMA, SNORM, PERTU, PERTUF, rpos
  LOGICAL            :: restarted
  REAL(mk), dimension(3) :: rbuf
  integer,  dimension(3) :: ibuf
  CHARACTER(len=256)     :: cbuf
  INTEGER                :: icoords(3), icpu,jcpu,kcpu, irank
  REAL(mk)               :: dxsub(3)
  REAL(mk), DIMENSION(4)       :: tsp
  INTEGER,  DIMENSION(3,6)     :: ibcdef
  REAL(mk), DIMENSION(3,1,1,1) :: ibcvalue
  REAL(mk), DIMENSION(:,:),POINTER     :: XXP => NULL()
  LOGICAL :: periods(3) = .TRUE.
  INCLUDE 'mpif.h'
  !-----------------------------------------------------
  !  memory measurements
  !-----------------------------------------------------
  INTEGER    :: eye
  INTEGER*4  :: fragments
  INTEGER*8  :: total_free, largest_free, total_used
  INTEGER    :: heap_info

  REAL(mk) :: rad1t, rad1r, theta1, radstr, rad1sq, radstrength, vrrad, rletrad
  REAL(mk) :: gamma, rad1sqTILDA


  tot = 200000
  dt_adapt = .FALSE.
  wbcdef = .FALSE.
  wvic_muscle = .FALSE.
  wvic_dgamma_scheme = 0
  wvic_compvel_scheme = 0
  
  wvic_les = .FALSE.
  wvic_les_model = 0
  
  les_tdm_C_comp = 0.0_mk
  les_tdm_C_dil  = 0.0_mk
  les_tdmclipped_C = 0.0_mk
  
  les_hypervisc_C = 0.0_mk
  les_hypervisc_tzero = 1.0_mk
  
  unboundedxy = .FALSE.
  trailvortex = .FALSE.
  wvic_noise_amp = 0.0_mk
  wvic_noise_basemode = 0.0_mk
  
  helvortex_period = 0.0_mk
  helvortex_ratio = 0.0_mk
  
  target_re = -1.0_mk
  ALLOCATE(cow(4))
  cow(1) = 1.0_mk
  cow(2) = 1.0_mk
  cow(3) = 1.0_mk
  cow(4) = 0.0_mk
  piston_RE_ratio = -1.0_mk
  
  u_infty = 0.0_mk
  
  nkilldiv = 0
  
  ndump = 10
  dumpvel = .FALSE.
  dumppsi = .FALSE.
  dumpkespec = .FALSE.
  dumpkespeckx = .FALSE.
  
  keepwcentered = .FALSE.

  !---------------------------------------
  ! JTR and the julenisser were here
  !---------------------------------------
  stepfunction_band = 1.0_mk
  penalization_lambda = 1.0_mk
  lambda = 0.0_mk
  penalization_clearrhs = .true.
  penalization_implicit = .false.
  penalization_oneterm = .true.
  dt_adapt_fraction = 1.0_mk
  onset_ramptime = 0.0_mk
  dt_rampsteps = 0
  dt_rampstart = 1.0E-3_mk

  object_offset = 0.0_mk

  cylinder_no = 2

  thickness05init = 0.0_mk

  step_function = 0
  step1_interval = 6.0_mk
  step1_linearfraction = 0.5_mk
  step1_offset = 0.0_mk

  u_solid = 0.0_mk
  harmonic_amplitude = 0.0_mk
  harmonic_period = 0.0_mk
  harmonic_phase = 0.0_mk

  stl_check_bounding       = .true.
  stl_double_check         = .false.
  stl_check_intersections  = .true.
  stl_min_phys_inside      = .false.
  stl_nonverbose           = .true.
  stl_translate = 0.0_mk
  stl_scale = 0.0_mk
  stl_mollify = 9
  stl_stencil = 1
  stl_subdivide = 1

  medusa_scale   = 1.0_mk
  medu_move_pos  = 0.0_mk
  medu_move_vel  = 0.0_mk

  !---------------------------------------
  ! JTR out
  !---------------------------------------
  !----------------------------------------------------------------------------!
  !  dimension 3, lda = 3
  !----------------------------------------------------------------------------!
  dime  = 3; lda = 3
  time  = 0.0_mk
  itime = 0
  
  !----------------------------------------------------------------------------!
  !  first of all initialize mpi
  !----------------------------------------------------------------------------!
  ALLOCATE(min_physg(dime),max_physg(dime),stat=istat)
  IF(istat.NE.0) CALL wvic_died
  ALLOCATE(bcdef(2*dime),nx(dime),stat=istat)
  IF(istat.NE.0) CALL wvic_died
  ALLOCATE(ghostsize(dime),stat=istat)
  IF(istat.NE.0) CALL wvic_died


  CALL MPI_Init(info)
  CALL MPI_Comm_Size(MPI_COMM_WORLD,nproc,info)
  CALL MPI_Comm_Rank(MPI_COMM_WORLD,rank, info)
  IF(rank.EQ.0) CALL wvic_readparams (ctrlfile,info)
  CALL MPI_BCast(ndims,3,MPI_INTEGER,0,MPI_COMM_WORLD,info)
! CALL MPI_Get_Processor_Name(procname,iprocname,info)
!call mpi_barrier(MPI_COMM_WORLD,info)
!print*,rank,' procname iprocname info = ',info,iprocname
  CALL MPI_CART_CREATE(MPI_COMM_WORLD,lda,ndims,periods,.TRUE.,comm,info)
  PRINT*,'rank ',rank,' communicator: ',comm,info
  CALL MPI_CART_COORDS(comm,rank,lda,coord,info)
  !  lets say hello
  IF(rank.EQ.0) WRITE(*,'(I5,2A)') rank,' started on ',procname(1:iprocname)
  maxsubs = 10000
  !----------------------------------------------------------------------------!
  !  prepare the globals
  !----------------------------------------------------------------------------!
  IF(mk.EQ.KIND(1.0d0)) THEN
     mpi_prec = MPI_DOUBLE_PRECISION
     prec     = ppm_kind_double
  ELSE
     mpi_prec = MPI_REAL
     prec     = ppm_kind_single
  END IF

  iverbose = 0

  !----------------------------------------------------------------------------!
  !  thanks to ivo: read the parameters
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!
  !  distribute the parameters to the other nodes
  !----------------------------------------------------------------------------!

  
  !----------------------------------------------------------------------------!
  ! send runtag
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(runtag,256,MPI_CHARACTER,0,comm,info)
  CALL MPI_BCast(iruntag,1,MPI_INTEGER,0,comm,info)

  !----------------------------------------------------------------------------!
  ! send nx
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(nx,dime,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send min_physg/max_physg
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(min_physg,dime,mpi_prec,0,comm,info)
  CALL MPI_BCast(max_physg,dime,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send ghostsize
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(ghostsize,dime,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send timestep
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(dt,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send timestep adaptivity flag
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(dt_adapt,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send target re
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(target_re,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send max timestep
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(dt_max,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send end time
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(tend,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send nu
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(nu,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send schmidt number
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(sc,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send number of levels
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(maxlev,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send renormalize?
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_renormalize,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send number of modes
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_noise_nmodes,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  ! noise amplitude
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_noise_amp,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! noise fundamental wavelength
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_noise_basemode,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! muscle scheme
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_muscle,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! which flow case?
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(flow_case,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  ! les model?
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_les,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! which les model?
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_les_model,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  ! les parameters
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(les_tdm_C_comp,  1,mpi_prec,0,comm,info)
  CALL MPI_BCast(les_tdm_C_dil,   1,mpi_prec,0,comm,info)
  
  CALL MPI_BCast(les_tdmclipped_C,1,mpi_prec,0,comm,info)
  
  CALL MPI_BCast(les_hypervisc_C, 1,mpi_prec,0,comm,info)
  CALL MPI_BCast(les_hypervisc_tzero, 1,mpi_prec,0,comm,info)
  
  !----------------------------------------------------------------------------!
  ! wall configuration
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wbcdef,6,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! Are we dealing with trailing vortices formation
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(trailvortex,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! Are we dealing with unbounded in ZY
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(unboundedxy,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  !  w scheme
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_dgamma_scheme,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  !  u scheme
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(wvic_compvel_scheme,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  !  n Divergence free projection
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(nkilldiv,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  !  n dump
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(ndump,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(nrestart,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  !  dump velocity
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(dumpvel,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  !  dump streamfunction
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(dumppsi,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  !  dump kinetic energy spectrum
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(dumpkespec,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  !  dump kinetic energy spectrum as a function of kz
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(dumpkespeckx,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  !  Keep omega centered in the domain???
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(keepwcentered,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send u_infty
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(u_infty,dime,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send neetcdf_restart
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(netcdf_restart,1,MPI_LOGICAL,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send netcdf_time
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(netcdf_time,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send netcdf_itime
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(netcdf_itime,1,MPI_INTEGER,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send netcdf_dt
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(netcdf_dt,1,mpi_prec,0,comm,info)

  !----------------------------------------------------------------------------!
  ! send piston_radius
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(piston_radius,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send piston_re_ratio
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(piston_re_ratio,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send piston_radius
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(piston_vort_ratio,1,mpi_prec,0,comm,info)

  !----------------------------------------------------------------------------!
  ! send trailing vortices spans
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(trailvortex_b1,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(trailvortex_b2,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send trailing vortices radii
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(trailvortex_a1,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(trailvortex_a2,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send trailing vortices vertical offset
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(trailvortex_z12,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send trailing vortices gamma
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(trailvortex_gamma,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send trailing vortices gamma ratio
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(trailvortex_r,1,mpi_prec,0,comm,info)
  
  !----------------------------------------------------------------------------!
  ! send helical vortices radii and number
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(helvortex_n,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(helvortex_nvortex,1,MPI_INTEGER,0,comm,info) ! JTR 200703
  CALL MPI_BCast(helvortex_ratio,1,mpi_prec,0,comm,info) ! added by JHW 20061221
  CALL MPI_BCast(helvortex_r1,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(helvortex_r2,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(helvortex_a1,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(helvortex_a2,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(helvortex_precision_isec1,1,mpi_prec,0,comm,info) ! JTR 200705
  CALL MPI_BCast(helvortex_precision_isec2,1,mpi_prec,0,comm,info) ! JTR 200705
  CALL MPI_BCast(helvortex_period,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(helvortex_gamma,1,mpi_prec,0,comm,info)

  !---------------------------------------
  ! JTR and the julenisser were here
  !---------------------------------------
  !----------------------------------------------------------------------------!
  ! send penalization arguments and sphere radius !JTR
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(penalization_lambda,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(stepfunction_band,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(penalization_clearrhs,1,MPI_LOGICAL,0,comm,info)
  CALL MPI_BCast(penalization_implicit,1,MPI_LOGICAL,0,comm,info)
  CALL MPI_BCast(penalization_oneterm,1,MPI_LOGICAL,0,comm,info)
  CALL MPI_BCast(sphere_radius,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(object_offset,dime,mpi_prec,0,comm,info)
  CALL MPI_BCast(cylinder_radius,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(cylinder_no,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(onset_ramptime,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(dt_adapt_fraction,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(nforces,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(block_wx,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(block_wy,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(block_wz,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(thickness05init,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send step function arguments
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(step_function,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(step1_interval,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(step1_linearfraction,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(step1_offset,1,mpi_prec,0,comm,info)
  !----------------------------------------------------------------------------!
  ! send solid movement arguments
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(u_solid,3,mpi_prec,0,comm,info)
  CALL MPI_BCast(harmonic_amplitude,3,mpi_prec,0,comm,info)
  CALL MPI_BCast(harmonic_period,3,mpi_prec,0,comm,info)
  CALL MPI_BCast(harmonic_phase,3,mpi_prec,0,comm,info)

  ALLOCATE(forcecv(6*nforces))
  CALL MPI_BCast(forcecv,6*nforces,mpi_prec,0,comm,info)
  ALLOCATE(vforcenocaold(3*nforces))
  vforcenocaold=0.0_mk  

  !----------------------------------------------------------------------------!
  ! send STL arguments
  !----------------------------------------------------------------------------!
  CALL MPI_BCast(stl_translate,3,mpi_prec,0,comm,info)
  CALL MPI_BCast(stl_scale,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(stl_check_bounding,1,MPI_LOGICAL,0,comm,info)
  CALL MPI_BCast(stl_double_check,1,MPI_LOGICAL,0,comm,info)
  CALL MPI_BCast(stl_check_intersections,1,MPI_LOGICAL,0,comm,info)
  CALL MPI_BCast(stl_min_phys_inside,1,MPI_LOGICAL,0,comm,info)
  CALL MPI_BCast(stl_mollify,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(stl_stencil,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(stl_subdivide,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(stlfile,256,MPI_CHARACTER,0,comm,info)
  CALL MPI_BCast(medusa_file,256,MPI_CHARACTER,0,comm,info)
  CALL MPI_BCast(medusa_scale,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(medu_mass_orig,1,mpi_prec,0,comm,info)
  CALL MPI_BCast(medusa_frames,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(medusa_period,1,mpi_prec,0,comm,info)

  IF(MINVAL(ghostsize) .LT. (stl_stencil+1) .AND. (flow_case .EQ. 14)) THEN
    WRITE(0,*) 'Ghostlayer insufficient for stl stencil, exiting\n'
    call mpi_finalize(info)
    stop
    call wvic_died
  END IF


  ALLOCATE(subcell(1:stl_subdivide,1:stl_subdivide,1:stl_subdivide))
  !---------------------------------------
  ! JTR out
  !---------------------------------------

  !----------------------------------------------------------------------------!
  ! send wvic_multigrid
  !----------------------------------------------------------------------------!
  IF(rank.EQ.0) THEN
     IF(wvic_multigrid) THEN
        iwvic_multigrid = 1
     ELSE
        iwvic_multigrid = 0
     END IF
  END IF
  ibuf(1) = iwvic_multigrid
  CALL MPI_BCast(iwvic_multigrid,1,MPI_INTEGER,0,comm,info)
  IF(rank.NE.0) THEN
     IF(iwvic_multigrid.EQ.1) THEN
        wvic_multigrid = .TRUE.
     ELSE
        wvic_multigrid = .FALSE.
     END IF
  END IF
  !----------------------------------------------------------------------------!
  ! send wvic_verbose 
  !----------------------------------------------------------------------------!
  IF(rank.EQ.0) THEN
     IF(verbose) THEN
        iverbose = 1
     ELSE
        iverbose = 0
     END IF
  END IF
  ibuf(1) = iverbose
  CALL MPI_BCast(iverbose,1,MPI_INTEGER,0,comm,info)
  IF(rank.NE.0) THEN
     IF(iverbose.EQ.1) THEN
        verbose = .TRUE.
     ELSE
        verbose = .FALSE.
     END IF
  END IF
  !----------------------------------------------------------------------------!
  ! done sending, now stuff
  !----------------------------------------------------------------------------!
  fast_rmsh = .TRUE.
  debug =  0
  tol   = INT(LOG10(EPSILON(1.0e-11)))+1
  if(mk.eq.kind(1.0d0)) then
    tol   = -10
  else
    tol   = -3
  end if
  cutoff=-TINY(cutoff)
  krnl = ppm_param_rmsh_kernel_mp4
  CALL ppm_init(dime,prec,tol,comm,debug,info, ilogfile)
  IF (info.NE.0) CALL wvic_died
  !  say that everythings fine
  if(verbose) CALL ppm_write(rank,'wvic_init_cart_crow','parallel things init''d..',info)
  IF(MAXVAL(ghostsize).GT.1) THEN
     mg_order = ppm_param_order_4
  ELSE
     mg_order = ppm_param_order_2
  END IF
  
  !----------------------------------------------------------------------------!
  ! load restart file if exists
  !----------------------------------------------------------------------------!
  CALL wvic_suck_restart(info)
  IF(info.EQ.0) THEN
     restarted = .TRUE.
  ELSE
     restarted = .FALSE.
  END IF
  info = 0 

  IF(dt.LE.0.0_mk)   dt    = 1.0e-9_mk
  Sci = 1.0_mk / Sc
  info = 0
  
  !----------------------------------------------------------------------------!
  !  THEN INITIALIZE THE PPM
  !----------------------------------------------------------------------------!
  

  !----------------------------------------------------------------------------!
  !  mesh spacing
  !----------------------------------------------------------------------------!
  dx = (max_physg(1)-min_physg(1))/REAL(nx(1)-1,mk)
  dy = (max_physg(2)-min_physg(2))/REAL(nx(2)-1,mk)
  dz = (max_physg(3)-min_physg(3))/REAL(nx(3)-1,mk)

  !----------------------------------------------------------------------------!
  !  Create the topology
  !----------------------------------------------------------------------------!
  bcdef         = ppm_param_bcdef_periodic
  !----------------------------------------------------------------------------!
  !  Patch for walls
  !----------------------------------------------------------------------------!
  !DO i=1,6
  !   IF(wbcdef(i)) bcdef(i) = ppm_param_bcdef_dirichlet
  !END DO
  !bcdef = ppm_param_bcdef_periodic

  !----------------------------------------------------------------------------!
  !  Patch for trailing vortices
  !----------------------------------------------------------------------------!
  IF (trailvortex) THEN
     bcdef(5) = ppm_param_bcdef_dirichlet
     bcdef(6) = ppm_param_bcdef_neumann
  END IF
  
  !----------------------------------------------------------------------------!
  !  [TODO] Adapt for eventual use of ppm Multigrid
  !  Fri Nov 11 13:49:00 CET 2005
  !----------------------------------------------------------------------------!
  ibcvalue(:,:,:,:)=0.0_MK
  ibcdef(1,:)=bcdef(:)
  ibcdef(2,:)=bcdef(:)
  ibcdef(3,:)=bcdef(:)
  
  !  let the ppm assign ids
  topo_id = -1
  mesh_id = -1

  NULLIFY(ndata,isublist,sub2proc,sub_cost,proc_speed,sub_cost)
  ALLOCATE(proc_speed(nproc))
  
#ifdef __TESTFLOATDECOMP !??? To delete?
  WRITE(*,*) 'nsubs requested', nsubs
  nsublist = flow_case
  ALLOCATE(isublist(nsublist))
  ALLOCATE(min_sub(3,nsubs))
  ALLOCATE(max_sub(3,nsubs))
  ALLOCATE(ndata(3,nsubs))
  ALLOCATE(istart(3,nsubs))
  ALLOCATE(sub2proc(nsubs))
  !-----------------------------------------------------
  decomposition = ppm_param_decomp_cartesian
  assigning     = ppm_param_assign_internal
  !----------------------------------------------------------------------------!
  ! create the topology
  CALL ppm_mktopo(xxp, -1, nx, decomposition, assigning, min_physg, max_physg,&
       &          bcdef, ghostsize, topo_id, mesh_id, min_sub, max_sub,      &
       &          sub_cost,  sub2proc, nsubs, isublist, nsublist, istart,    &
       &          ndata, info,ndom=flow_case)
  IF (INFO.EQ.0) WRITE(*,*) 'Created nsubs ', nsubs
  IF (INFO.NE.0) WRITE(*,*) 'Could not create nsubs ', flow_case
  CALL wvic_died
  CALL MPI_Abort(comm,info)
#endif

  !-----------------------------------------------------
  !  create user defined topology
  nsubs    = ndims(1)*ndims(2)*ndims(3)
  nsublist = 1
  ALLOCATE(isublist(nsublist))
  ALLOCATE(min_sub(3,nsubs))
  ALLOCATE(max_sub(3,nsubs))
  ALLOCATE(ndata(3,nsubs))
  ALLOCATE(istart(3,nsubs))
  ALLOCATE(sub2proc(nsubs))
  ndata(1,:) = (nx(1)-1)/ndims(1)+1
  ndata(2,:) = (nx(2)-1)/ndims(2)+1
  ndata(3,:) = (nx(3)-1)/ndims(3)+1
  dxsub(1)   = (max_physg(1)-min_physg(1))/REAL(ndims(1),mk)
  dxsub(2)   = (max_physg(2)-min_physg(2))/REAL(ndims(2),mk)
  dxsub(3)   = (max_physg(3)-min_physg(3))/REAL(ndims(3),mk)
  isubl = 1
  DO kcpu=1,ndims(3)
     DO jcpu=1,ndims(2)
        DO icpu=1,ndims(1)
           istart(1,isubl)  = (icpu-1)*(ndata(1,isubl)-1)+1
           istart(2,isubl)  = (jcpu-1)*(ndata(2,isubl)-1)+1
           istart(3,isubl)  = (kcpu-1)*(ndata(3,isubl)-1)+1
           icoords = (/icpu-1,jcpu-1,kcpu-1/)
           CALL MPI_CART_RANK(comm,icoords,irank,info)
           sub2proc(isubl)  = irank
           IF(rank.EQ.irank) isublist(1) = isubl
           min_sub(1,isubl) = min_physg(1) + REAL(icpu-1,mk)*dxsub(1)
           min_sub(2,isubl) = min_physg(2) + REAL(jcpu-1,mk)*dxsub(2)
           min_sub(3,isubl) = min_physg(3) + REAL(kcpu-1,mk)*dxsub(3)
           max_sub(1,isubl) = min_physg(1) + REAL(icpu  ,mk)*dxsub(1)
           max_sub(2,isubl) = min_physg(2) + REAL(jcpu  ,mk)*dxsub(2)
           max_sub(3,isubl) = min_physg(3) + REAL(kcpu  ,mk)*dxsub(3)
           isubl = isubl + 1
        END DO
     END DO
  END DO
  !-----------------------------------------------------
  decomposition = ppm_param_decomp_user_defined  
  assigning     = ppm_param_assign_user_defined
  !----------------------------------------------------------------------------!
  ! create the topology
  CALL ppm_mktopo(xxp, -1, nx, decomposition, assigning, min_physg, max_physg,&
       &          bcdef, ghostsize, topo_id, mesh_id, min_sub, max_sub,      &
       &          sub_cost,  sub2proc, nsubs, isublist, nsublist, istart,    &
       &          ndata, info)

  IF(info.NE.0) CALL wvic_died
  !----------------------------------------------------------------------------!

  IF (unboundedxy) THEN
      !  let the ppm assign ids
      dblxy_topo_id = -1
      dblxy_mesh_id = -1
      ALLOCATE(dblxy_min_physg(dime),dblxy_max_physg(dime),stat=istat)
      IF(istat.NE.0) CALL wvic_died
      dblxy_min_physg(:) = min_physg(:)
      dblxy_max_physg(1) = 2.0_mk*max_physg(1) - min_physg(1)
      dblxy_max_physg(2) = 2.0_mk*max_physg(2) - min_physg(2)
      dblxy_max_physg(3) = 2.0_mk*max_physg(3) - min_physg(3)
      
      nsubs    = ndims(1)*ndims(2)*ndims(3)
      nsublist = 1
      
      ALLOCATE(dblxy_isublist(nsublist))
      ALLOCATE(dblxy_min_sub(3,nsubs))
      ALLOCATE(dblxy_max_sub(3,nsubs))
      ALLOCATE(dblxy_ndata(3,nsubs))
      ALLOCATE(dblxy_istart(3,nsubs))
      ALLOCATE(dblxy_sub2proc(nsubs))
      dblxy_ndata(1,:) = 2*(nx(1)-1)/ndims(1)+1
      dblxy_ndata(2,:) = 2*(nx(2)-1)/ndims(2)+1
      dblxy_ndata(3,:) =   (nx(3)-1)/ndims(3)+1
      dxsub(1)   = (dblxy_max_physg(1)-dblxy_min_physg(1))/REAL(ndims(1),mk)
      dxsub(2)   = (dblxy_max_physg(2)-dblxy_min_physg(2))/REAL(ndims(2),mk)
      dxsub(3)   = (dblxy_max_physg(3)-dblxy_min_physg(3))/REAL(ndims(3),mk)
      isubl = 1
      DO kcpu=1,ndims(3)
         DO jcpu=1,ndims(2)
            DO icpu=1,ndims(1)
               dblxy_istart(1,isubl)  = (icpu-1)*(dblxy_ndata(1,isubl)-1)+1
               dblxy_istart(2,isubl)  = (jcpu-1)*(dblxy_ndata(2,isubl)-1)+1
               dblxy_istart(3,isubl)  = (kcpu-1)*(dblxy_ndata(3,isubl)-1)+1
               icoords = (/icpu-1,jcpu-1,kcpu-1/)
               CALL MPI_CART_RANK(comm,icoords,irank,info)
               dblxy_sub2proc(isubl)  = irank
               IF(rank.EQ.irank) dblxy_isublist(1) = isubl
               dblxy_min_sub(1,isubl) = dblxy_min_physg(1) + REAL(icpu-1,mk)*dxsub(1)
               dblxy_min_sub(2,isubl) = dblxy_min_physg(2) + REAL(jcpu-1,mk)*dxsub(2)
               dblxy_min_sub(3,isubl) = dblxy_min_physg(3) + REAL(kcpu-1,mk)*dxsub(3)
               dblxy_max_sub(1,isubl) = dblxy_min_physg(1) + REAL(icpu  ,mk)*dxsub(1)
               dblxy_max_sub(2,isubl) = dblxy_min_physg(2) + REAL(jcpu  ,mk)*dxsub(2)
               dblxy_max_sub(3,isubl) = dblxy_min_physg(3) + REAL(kcpu  ,mk)*dxsub(3)
               isubl = isubl + 1
            END DO
         END DO
      END DO
      !-----------------------------------------------------
      decomposition = ppm_param_decomp_user_defined  
      assigning     = ppm_param_assign_user_defined
      !----------------------------------------------------------------------------!
      ! create the topology
      CALL ppm_mktopo(xxp, -1, nx, decomposition, assigning, dblxy_min_physg, dblxy_max_physg,&
       &          bcdef, ghostsize, dblxy_topo_id, dblxy_mesh_id,  dblxy_min_sub,  dblxy_max_sub,      &
       &          dblxy_sub_cost,  dblxy_sub2proc, nsubs, dblxy_isublist, nsublist, dblxy_istart,    &
       &          dblxy_ndata, info)
      IF(info.NE.0) CALL wvic_died
      !----------------------------------------------------------------------------!
  END IF


  !-----------------------------------------------------
  !  measure memory here
  !-----------------------------------------------------
  !eye = heap_info(fragments, total_free, largest_free, total_used)
  !WRITE(msg,*) 'heap_info fragments =',fragments,' total_free =', &
  !     & total_Free,' largest_free = ',largest_free/1024/1024,' TOTAL USED = '&
  !     &,total_used/1024/1024, 'MB i=',eye
  !IF(rank.EQ.0.OR.rank.EQ.1) CALL ppm_Write(rank,'wvic_tvdrk3',msg,info)
  
  !----------------------------------------------------------------------------!
  !  Set the iswall arrays
  !----------------------------------------------------------------------------!
    CALL wvic_iswall(info)

    IF(verbose) THEN
     IF(rank.EQ.0) THEN
        DO isubl=1,nsubs
           WRITE(msg,*) 'S(',isubl,') [',min_sub(1,isubl),',',max_sub(1,isubl),&
                &']x[',min_sub(2,isubl),',',max_sub(2,isubl),&
                &']x[',min_sub(3,isubl),',',max_sub(3,isubl),']'
           CALL ppm_write(rank,'wvic_init_cart',msg,info)
        END DO
     END IF
  END IF
  
  !----------------------------------------------------------------------------!
  !  Allocate the armada of fields
  !----------------------------------------------------------------------------!
  CALL wvic_alloc_field(field_up,dime,info)
  IF(info.NE.0) CALL wvic_died
  CALL wvic_alloc_field(field_wp,lda,info)
  IF(info.NE.0) CALL wvic_died
  CALL wvic_alloc_field(field_dwp,lda,info)
  IF(info.NE.0) CALL wvic_died

  !---------------------------------------
  ! JTR and the julenisser were here
  !---------------------------------------
  CALL wvic_alloc_field_s(field_H,info)
  IF(info.NE.0) CALL wvic_died
  CALL wvic_alloc_field(field_ubar,lda,info)
  IF(info.NE.0) CALL wvic_died

  !-----------------------------------------------------------------------------
  ! Allocate array for storing z-velocity on the line of pressure integration.
  ! to use for time derivative of z-velocity. Is initialized as 0
  !-----------------------------------------------------------------------------
  ALLOCATE(u1old(1:MAXVAL(ndata(1,isublist(1:nsublist))),nsublist))
  ALLOCATE(u2old(1:MAXVAL(ndata(2,isublist(1:nsublist))),nsublist))
  ALLOCATE(u3old(1:MAXVAL(ndata(3,isublist(1:nsublist))),nsublist))
  u1old(:,:)=0.0_mk
  u2old(:,:)=0.0_mk
  u3old(:,:)=0.0_mk
  !---------------------------------------
  ! JTR out
  !---------------------------------------

  !----------------------------------------------------------------------------!
  !  here the INITIAL CONDITION IS SET UP
  !  JTR setting up step function
  !----------------------------------------------------------------------------!
  IF (flow_case .GT. 8) THEN
    IF (step_function .EQ. 1) THEN
      CALL init_stepfunction1(info)
    END IF
  END IF
  IF(.NOT.restarted) THEN
     IF (.NOT.trailvortex) THEN
        SELECT CASE(flow_case)
        CASE(0)
           !-----------------------------------------------------
           !  Tubes
           !-----------------------------------------------------
           CALL wvic_init_physics_0
        CASE(1)
           !-----------------------------------------------------
           !  Poiseuille
           !-----------------------------------------------------
           PRINT*,'CASE 1 HAS BEEN UNCONFIGURED, USE FLOWCASE=0'
        CASE(2)
           !-----------------------------------------------------
           !  Ring ring.
           !-----------------------------------------------------
           CALL wvic_init_physics_4
        CASE(3)
           !-----------------------------------------------------
           !  Ring ring.
           !-----------------------------------------------------
           CALL wvic_init_physics_5
        CASE(4)
           !-----------------------------------------------------
           !  Random vorticity
           !-----------------------------------------------------
           CALL wvic_init_physics_6
        CASE(5)
           !-----------------------------------------------------
           !  Four tubes
           !-----------------------------------------------------
           CALL wvic_init_tvphysics_0
        CASE(6)
           !-----------------------------------------------------
           !  Four tubes with analytic divergence free perturb.
           !-----------------------------------------------------
           CALL wvic_init_tvphysics_1
         CASE(7)
           !-----------------------------------------------------
           ! Four tubes with sine pertubations
           !-----------------------------------------------------
           CALL wvic_init_tvphysics_2
         CASE(8)
           !-----------------------------------------------------
           ! Helical vortices
           !-----------------------------------------------------
           CALL wvic_init_helphysics_0
         CASE(9)
       !---------------------------------------------------------
       ! Maybe remove all the wvic_init_stepfunc with the actual functions
       !---------------------------------------------------------
           !-----------------------------------------------------
           ! A block or a plate (equivalent of two parallel plates)
           !-----------------------------------------------------
           CALL wvic_init_stepfunc
	   CALL wvic_stepfunc_plate
	   CALL wvic_barycentric
           IF (thickness05init .EQ. 0.0_mk) THEN
             CALL wvic_init_onset
           ELSE
             WRITE(msg,*) 'Initialising analytic Poiseuille flow (z) profile\n'
             WRITE(0,*) msg
             CALL wvic_init_poiseulle
           ENDIF 
         CASE(10)
           !-----------------------------------------------------
           ! Onset flow past sphere
           !-----------------------------------------------------
           CALL wvic_init_stepfunc
	   CALL wvic_stepfunc_sphere
!           CALL wvic_calculate_mass(info)
           CALL wvic_init_onset
         CASE(11)
           !-----------------------------------------------------
           ! Medusa
           !-----------------------------------------------------
           CALL wvic_init_stepfunc
	   CALL wvic_readmedusa
           CALL wvic_init_onset
         CASE(12)
           !-----------------------------------------------------
           ! Cylinder array - Dong et al.
           !-----------------------------------------------------
           CALL wvic_init_stepfunc
           CALL wvic_stepfunc_cylinderarray
           CALL wvic_init_onset
         CASE(13)
           !-----------------------------------------------------
           ! Harmonic oscilations of plate
           !-----------------------------------------------------
           CALL wvic_init_stepfunc
	   CALL wvic_stepfunc_plate
	   CALL wvic_barycentric
           CALL wvic_init_onset
         CASE(14)
           !-----------------------------------------------------
           ! STL
	   ! wvic_readstl must be called in wvic_init_cart
           !-----------------------------------------------------
           CALL wvic_readstl (info)
           CALL wvic_init_stepfunc
	   CALL wvic_barycentric
           CALL wvic_init_onset

        END SELECT
     ELSE
        SELECT CASE(flow_case)
        CASE(0)
           !-----------------------------------------------------
           !  4 straight tubes
           !-----------------------------------------------------
           CALL wvic_init_tvphysics_0
        CASE(1)
           !-----------------------------------------------------
           !  4 tubes with divergence-free perturbation 
           !-----------------------------------------------------
           CALL wvic_init_tvphysics_1
        CASE(2)
           !-----------------------------------------------------
           !  4 tubes with divergence free perturbation at x=0
           !-----------------------------------------------------
           CALL wvic_init_tvphysics_0
           ! Create upstream group
           CALL MPI_COMM_GROUP(MPI_COMM_WORLD, whole_group,info)
           ALLOCATE(upstream_ranks(ndims(1)*ndims(2)),stat=istat)
           IF(istat.NE.0) CALL wvic_died
           n_upstream = 0
           kcpu = 1
           DO jcpu=1,ndims(2)
             DO icpu=1,ndims(1)
                icoords = (/icpu-1,jcpu-1,kcpu-1/)
                CALL MPI_CART_RANK(comm,icoords,irank,info)
                n_upstream = n_upstream + 1
                upstream_ranks(n_upstream) = irank
             END DO
           END DO
           CALL MPI_GROUP_INCL(whole_group, n_upstream, upstream_ranks, upstream_group,info)
           CALL MPI_COMM_CREATE(MPI_COMM_WORLD, upstream_group, upstream_comm,info)
           CALL MPI_GROUP_RANK(upstream_group, upstream_rank,info)
           WRITE(msg,*) '[upstream group creation] proc ranks: ', upstream_ranks(1:n_upstream)
           IF (verbose) CALL ppm_write(rank,'wvic_init_cart',msg,info)
        END SELECT
     END IF
  !--------------------------------------------------------------------------
  ! Interpolate step function to particles
  ! (requires ppm_module_rmsh_create_part)
  !--------------------------------------------------------------------------
  IF (object_move .EQV. .true.) THEN
    CALL ppm_rmsh_create_part(xpc,npc,hpc,field_H,topo_id,mesh_id,&
       & (/0.0_mk,1.01_mk/),info,resetpos=.TRUE.,field_wp = field_up,wp&
          & = up, lda2 = dime)
  ENDIF

  END IF
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  !  [TODO] Do we need dwp?
  !  Fri Nov 11 13:55:34 CET 2005
  !----------------------------------------------------------------------------!
  ALLOCATE(up(dime,np),dwp(lda,np),stat=istat)
  up = 0.0_mk
  dwp = 0.0_mk
  IF(istat.NE.0) THEN
     PRINT*,'died'
     CALL wvic_died
     PRINT*,'died done'
     STOP
  END IF
  !----------------------------------------------------------------------------!
  !  map the stuff correctly
  !----------------------------------------------------------------------------!
  CALL wvic_map_particles (ppm_param_map_global,info)

  !-----------------------------------------------------
  !  measure memory here
  !-----------------------------------------------------
  !eye = heap_info(fragments, total_free, largest_free, total_used)
  !WRITE(msg,*) '[after map_global] heap_info fragments =',fragments,' total_free =', &
  !     & total_Free,' largest_free = ',largest_free/1024/1024,' TOTAL USED = '&
  !     &,total_used/1024/1024, 'MB i=',eye
  !IF(rank.EQ.0.OR.rank.EQ.1) CALL ppm_Write(rank,'wvic_tvdrk3',msg,info)
  

  IF(verbose) CALL ppm_write(rank,'wvic_init_cart','done with initialization',info)
  
  CALL MPI_Barrier(comm,info)

END SUBROUTINE wvic_init_cart
