!------------------------------------------------------------------------------!
!* filename: module_wvic                                                     *!
!* project : ppm                                                              *!
!* purpose : module for the vortex client for ppm                             *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 10 14:12:59 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: module_wvic.F,v $
!  Revision 1.13  2006/10/25 11:35:50  menahel
!  re-enabled .rst dumping and added NRESTART parameter
!
!  Revision 1.12  2006/10/24 16:52:53  pchatela
!  Added KEEPWCENTERED option if one wants to keep vorticity centered in domain
!
!  Revision 1.11  2006/10/23 14:28:03  pchatela
!  Added parameters for hyperviscosity model
!  Fixed double Laplacian
!
!  Revision 1.10  2006/10/23 08:19:10  pchatela
!  LES models
!  Bugfixes in KE spectra and factor for Parseval identity
!  Removed the reset of noise in init_tvphysics_0 and_1
!
!  Revision 1.9  2006/10/17 14:17:40  pchatela
!  Added Ctrl parameters to handle LES models
!
!  Revision 1.8  2006/10/04 11:16:38  pchatela
!  Added a fundamental mode to the noise parameters.
!  If not given, the length of the domain is used.
!
!  Revision 1.7  2006/10/03 16:11:58  pchatela
!  Added spectra in the z direction (misnamed kx spectrum...)
!  Added spatial diagnostics, like kinetic energy, enstrophy, circulation
!  as functions of z, dumped at the frequency ndump
!
!  Revision 1.6  2006/09/27 09:30:21  pchatela
!  Fixes, spectra calculation,
!  most importantly: moved the u_infty out, so it does not kill the dgammadt
!
!  Revision 1.5  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.4  2006/09/01 15:47:03  pchatela
!  Added option to dump velocity and streamfunction fields
!
!  Revision 1.3  2006/07/26 14:09:32  pchatela
!  Clean ups and fixes
!
!  Revision 1.2  2006/07/26 07:51:25  pchatela
!  Added boolean trailvortex
!  Added periodic bcs with reset of particle strengths
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.6  2005/12/10 21:10:51  michaebe
!  added constant for TDM les model
!
!  Revision 1.5  2005/12/10 20:36:35  michaebe
!  added hyperviscosity parameter
!
!  Revision 1.4  2005/11/21 17:36:36  michaebe
!  major session. now let''s debug
!
!  Revision 1.3  2005/10/14 11:44:53  michaebe
!  added liswall variable
!
!  Revision 1.2  2005/10/07 14:33:53  pchatela
!  Generation of "Wall boundary conditions vs subdomain" flags
!  Enforcement of wall boundary conditions a la Thom
!  Added entries in Makefile
!
!  Revision 1.1  2005/09/28 11:40:23  michaebe
!  Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


MODULE module_wvic

  INTEGER, PARAMETER :: mk = KIND(1.0e0)           ! double precision
  CHARACTER(len=256)                :: runtag ! JHW moved here

  !----------------------------------------------------------------------------!
  !  data
  !----------------------------------------------------------------------------!
#ifdef __AGM_TESTRUN__
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_xp, field_xpn
  REAL(mk), DIMENSION(:,:,:,:  ), POINTER :: field_mon
  INTEGER                        :: agmiter
#endif
  REAL(mk), DIMENSION(:,:  )    , POINTER :: xp        ! positions
  REAL(mk), DIMENSION(:,:  )    , POINTER :: xp0       ! some odda positions
  REAL(mk), DIMENSION(:,:  )    , POINTER :: wp        ! vorticity
  REAL(mk), DIMENSION(:,:  )    , POINTER :: wp0       ! some odda vorticity
  REAL(mk), DIMENSION(:,:  )    , POINTER :: dwp       ! vorticity
  REAL(mk), DIMENSION(:,:  )    , POINTER :: dp        ! density gradient
  REAL(mk), DIMENSION(:,:  )    , POINTER :: up        ! velocities
  REAL(mk), DIMENSION(:,:  )    , POINTER :: rn        ! renormalization
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_up    ! up on the mesh
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_rhs   ! up on the mesh
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_wp    ! wp on the mesh 
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_wps   ! streamfunc
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_dwp   ! wp righthandside
  REAL(mk), DIMENSION(:,:,:,:  ), POINTER :: field_sp    ! sp on the mesh
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_dp    ! dp on the mesh
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_rn, field_l2q ! renormalization

  REAL(mk), DIMENSION(:,:,:),     POINTER :: wx, wy, wz  ! rmsh weights
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: field_ubar  ! Dirac on the mesh
  REAL(mk), DIMENSION(:,:,:,:  ), POINTER :: field_H     ! Heaviside on the mesh
  REAL(mk), DIMENSION(:,:  )    , POINTER :: xpc       ! levelset positions !JTR
  REAL(mk), DIMENSION(:    )    , POINTER :: hpc       ! levelset values
  REAL(mk), DIMENSION(:,:  )    , POINTER :: ubarp       ! solid velocities
  
  !----------------------------------------------------------------------------!
  !  resolution
  !----------------------------------------------------------------------------!
  INTEGER                           :: np          ! number of particles
  INTEGER                           :: npc         ! no. levelset particles
  INTEGER,  DIMENSION(:  ), POINTER :: nx          ! grid resolution
  INTEGER                           :: dime, lda
  REAL(mk), DIMENSION(:), POINTER            :: cow
  !----------------------------------------------------------------------------!
  !  parallel stuff
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)                :: procname    ! processor name
  INTEGER                           :: iprocname   ! length of processor name
  INTEGER                           :: rank, nproc        ! rank of this machine
  INTEGER                           :: ilogfile    ! log file unit
  INTEGER                           :: topo_id
  INTEGER                           :: mesh_id
  INTEGER                           :: debug
  INTEGER                           :: maxsubs, nsublist, maxlev
  INTEGER,  DIMENSION(:  ), POINTER :: sub2proc, isublist
  INTEGER,  DIMENSION(:,:), POINTER :: istart, ndata
  REAL(mk), DIMENSION(:  ), POINTER :: proc_speed, sub_cost
! CHARACTER(len=256)                :: runtag
  INTEGER                           :: iruntag, comm, mpi_prec
! INTEGER                           :: coord(3),ndims(3) ! JHW
  INTEGER, DIMENSION(3)             :: coord,ndims
  
  INTEGER                           :: whole_group, upstream_group, upstream_comm
  INTEGER                           :: upstream_rank
  INTEGER, DIMENSION(:), POINTER    :: upstream_ranks
  INTEGER                           :: n_upstream

  ! DOUBLE X-Y SIZE TOPO STUFF
  INTEGER                           :: dblxy_topo_id
  INTEGER                           :: dblxy_mesh_id
  INTEGER,  DIMENSION(:  ), POINTER :: dblxy_sub2proc, dblxy_isublist
  INTEGER,  DIMENSION(:,:), POINTER :: dblxy_istart, dblxy_ndata
  REAL(mk), DIMENSION(:  ), POINTER :: dblxy_proc_speed, dblxy_sub_cost

  !----------------------------------------------------------------------------!
  !  geometry
  !----------------------------------------------------------------------------!
  REAL(mk), DIMENSION(:  ), POINTER :: min_physg, max_physg
  REAL(mk), DIMENSION(:,:), POINTER :: min_sub, max_sub
  
  ! DOUBLE X-Y SIZE GEOMETRY STUFF
  REAL(mk), DIMENSION(:  ), POINTER :: dblxy_min_physg, dblxy_max_physg
  REAL(mk), DIMENSION(:,:), POINTER :: dblxy_min_sub, dblxy_max_sub
  
  INTEGER,  DIMENSION(:  ), POINTER :: bcdef
  LOGICAL,  DIMENSION(6)            :: wbcdef
  INTEGER,  DIMENSION(:,:,:), POINTER :: iswall
  LOGICAL,  DIMENSION(:,:  ), POINTER :: liswall
  LOGICAL                           :: trailvortex
  LOGICAL                           :: unboundedxy
  INTEGER,  DIMENSION(:  ), POINTER :: ghostsize
  REAL(mk)                          :: dx, dy, dz
  INTEGER                           :: tot ! time on target
  !----------------------------------------------------------------------------!
  !  other things
  !----------------------------------------------------------------------------!
  REAL(mk), PARAMETER               :: M_PI = 3.14159265358979
  REAL(mk)                          :: time, dt, mmt1, mmt2, cutoff, tend
  INTEGER                           :: itime, krnl, ndump, nrestart, nkilldiv
  LOGICAL                           :: dumpvel, dumppsi, dumpkespec,dumpkespeckx
  LOGICAL                           :: keepwcentered

  !----------------------------------------------------------------------------!
  !  physics
  !----------------------------------------------------------------------------!
  REAL(mk)                          :: nu, Sc, Sci ! viscosity, diffusivity
  REAL(mk)                          :: max_vorticity, target_re
  REAL(mk), DIMENSION(3)            :: u_cmass, u_infty
  INTEGER                           :: g_istage
  
  !----------------------------------------------------------------------------!
  !  poisson stuff
  !----------------------------------------------------------------------------!
  INTEGER, DIMENSION(3)             :: ffttopo_id, fftmesh_id
  

  !----------------------------------------------------------------------------!
  !  parameters
  !----------------------------------------------------------------------------!
  INTEGER, PARAMETER :: wvic_prm_vorticity    = 1
  INTEGER, PARAMETER :: wvic_prm_stream       = 2
  INTEGER, PARAMETER :: wvic_prm_velocity     = 3
  INTEGER, PARAMETER :: wvic_prm_dgammadt     = 4

  !----------------------------------------------------------------------------!
  !  some other stuff
  !----------------------------------------------------------------------------!
  CHARACTER(len=32)  :: version = '0.1a'
  LOGICAL            :: verbose
  LOGICAL            :: wvic_les
  LOGICAL            :: fast_rmsh
  LOGICAL            :: wvic_multigrid
  LOGICAL            :: wvic_renormalize
  REAL(mk)           :: wvic_noise_amp
  REAL(mk)           :: wvic_noise_basemode
  REAL(mk)           :: tube_radius
  INTEGER            :: wvic_noise_nmodes
  REAL(mk)           :: wvic_re
  LOGICAL            :: wvic_muscle
  INTEGER            :: wvic_dgamma_scheme
  INTEGER            :: wvic_compvel_scheme
  INTEGER            :: mg_order
  REAL(mk)           :: dt_max
  LOGICAL            :: dt_adapt
  INTEGER            :: flow_case

  !----------------------------------------------------------------------------!
  !  LES Thingies
  !----------------------------------------------------------------------------!
  INTEGER            :: wvic_les_model
  REAL(mk)           :: les_tdm_C_comp, les_tdm_C_dil
  REAL(mk)           :: les_hypervisc_C, les_hypervisc_tzero, les_tdmclipped_C

  !-----------------------------------------------------
  !  Restarting
  !-----------------------------------------------------
  LOGICAL            :: netcdf_restart
  INTEGER            :: netcdf_itime
  REAL(mk)           :: netcdf_time
  REAL(mk)           :: netcdf_dt


  !-----------------------------------------------------
  !  Piston
  !-----------------------------------------------------
  REAL(mk)           :: piston_RE_ratio
  REAL(mk)           :: piston_vort_ratio
  REAL(mk)           :: piston_radius
  
  !-----------------------------------------------------
  !  Trailing vortices
  !-----------------------------------------------------
  REAL(mk)           :: trailvortex_b1, trailvortex_b2
  REAL(mk)           :: trailvortex_a1, trailvortex_a2
  REAL(mk)           :: trailvortex_z12
  REAL(mk)           :: trailvortex_r, trailvortex_gamma
  
  !-----------------------------------------------------
  !  Helical vortices
  !-----------------------------------------------------
  REAL(mk)           :: helvortex_r1, helvortex_r2
  REAL(mk)           :: helvortex_a1, helvortex_a2
!  REAL(mk)           :: helvortex_a3
!  REAL(mk), DIMENSION(:), ALLOCATABLE :: gammadistribution
  REAL(mk)           :: helvortex_precision_isec1, helvortex_precision_isec2
  REAL(mk)           :: helvortex_ratio, helvortex_gamma
  REAL(mk)           :: helvortex_period
  INTEGER            :: helvortex_n, helvortex_nvortex
  !-----------------------------------------------------
  ! Geometry - general
  !-----------------------------------------------------
  REAL(mk), DIMENSION(3) :: object_offset
  !-----------------------------------------------------
  ! Spherical flow
  !-----------------------------------------------------
  REAL(mk)               :: sphere_radius
  !-----------------------------------------------------
  ! Cylinder flow
  !-----------------------------------------------------
  REAL(mk)               :: cylinder_radius
  INTEGER                :: cylinder_no
  !-----------------------------------------------------
  ! Plate / block   
  !-----------------------------------------------------
  REAL(mk)           :: block_wx,block_wy,block_wz
  REAL(mk)           :: thickness05init
  !-----------------------------------------------------
  ! Penalization factors
  !-----------------------------------------------------
  REAL(mk)           :: penalization_lambda           
  REAL(mk)           :: lambda
  REAL(mk)           :: dt_adapt_fraction
  REAL(mk)           :: stepfunction_band           
  LOGICAL            :: penalization_clearrhs
  LOGICAL            :: penalization_implicit
  LOGICAL            :: penalization_oneterm
  INTEGER            :: dt_rampsteps
  REAL(mk)           :: dt_rampstart
  REAL(mk)           :: onset_ramptime
  INTEGER            :: nforces
  REAL(mk), DIMENSION(3)                :: gforcewmomentold
  REAL(mk), DIMENSION(:), ALLOCATABLE   :: forcecv
!  REAL(mk), DIMENSION(3)               :: alphausum,alphawsum,alphaRHSsum ! for penalization force JTR slet
  REAL(mk), DIMENSION(:), ALLOCATABLE   :: vforcenocaold
  REAL(mk), DIMENSION(:,:), ALLOCATABLE :: u1old, u2old, u3old
!  REAL(mk), DIMENSION(:), ALLOCATABLE  :: u3old
  !-----------------------------------------------------
  ! Step function specific
  !-----------------------------------------------------
  INTEGER            :: step_function
  REAL(mk)           :: step1_interval, step1_linearfraction, step1_offset
  REAL(mk)           :: step1A,step1B,step1C,step1E
  REAL(mk)           :: step1varA,step1varB,step1varC,step1varD,step1varE, &
                      & step1varJ,step1varK,step1varL
  !-----------------------------------------------------
  ! Solid movement       
  !-----------------------------------------------------
  REAL(mk), DIMENSION(3) :: u_solid
  REAL(mk), DIMENSION(3) :: harmonic_amplitude
  REAL(mk), DIMENSION(3) :: harmonic_period
  REAL(mk), DIMENSION(3) :: harmonic_phase
  REAL(mk), DIMENSION(3) :: object_cmass
  REAL(mk)               :: object_mass
  LOGICAL                :: object_move
  !-----------------------------------------------------
  ! STL import
  !-----------------------------------------------------
  CHARACTER(len=256)                      :: medusa_file
  INTEGER                                 :: medusa_frames
  REAL(mk)                                :: medusa_scale
  REAL(mk)                                :: medu_mass_orig
  REAL(mk)                                :: medu_thrusty
  REAL(mk)                                :: medu_move_vel
  REAL(mk)                                :: medu_move_pos
  INTEGER                                 :: medusa_count
  REAL(mk), DIMENSION(:,:,:), ALLOCATABLE :: medusa_points_array
  REAL(mk), DIMENSION(:,:,:), ALLOCATABLE :: medusa_vel_array
  REAL(mk), DIMENSION(:,:), ALLOCATABLE   :: medusa_points
  REAL(mk), DIMENSION(:,:), ALLOCATABLE   :: medusa_vel
  REAL(mk), DIMENSION(:,:), ALLOCATABLE   :: medusa_panel
  REAL(mk), DIMENSION(:)  , ALLOCATABLE   :: medusa_panel_squared
  REAL(mk), DIMENSION(:,:), ALLOCATABLE   :: medusa_normal
  REAL(mk)                                :: medusa_dt
  REAL(mk)                                :: medusa_period
  REAL(mk), DIMENSION(:)  , ALLOCATABLE   :: medusa_maxx,medusa_maxy,medusa_miny
  REAL(mk)                                :: mmaxx,mmaxy,mminy
  !-----------------------------------------------------
  ! STL import
  !-----------------------------------------------------
  REAL(mk), DIMENSION(:), ALLOCATABLE     :: tri_denom,tri_udotv
  REAL(mk), DIMENSION(:,:), ALLOCATABLE   :: tri_norm,tri_base
  REAL(mk), DIMENSION(:), ALLOCATABLE     :: tri_udotu,tri_vdotv,tri_wdotw
  REAL(mk), DIMENSION(:,:), ALLOCATABLE   :: tri_vecu,tri_vecv,tri_vecw
  INTEGER                                 :: tri_count
  INTEGER                                 :: stl_subdivide
  REAL(mk), DIMENSION(:,:,:), ALLOCATABLE :: subcell
  REAL(mk)                                :: stl_scale
  LOGICAL                                 :: stl_check_bounding
  LOGICAL                                 :: stl_double_check
  LOGICAL                                 :: stl_check_intersections
  LOGICAL                                 :: stl_min_phys_inside
  LOGICAL                                 :: stl_nonverbose    
  INTEGER                                 :: stl_mollify
  INTEGER                                 :: stl_stencil
  REAL(mk), DIMENSION(3)                  :: stl_translate
  CHARACTER(len=256)                      :: stlfile 
  REAL(mk)                 :: bndminx,bndminy,bndminz,bndmaxx,bndmaxy,bndmaxz 
  
END MODULE module_wvic
