      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_data_fmm
      !-------------------------------------------------------------------------
      !
      ! Purpose       : fast mulipole method, data
      !               
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_fmm.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.10  2005/09/19 13:03:30  polasekb
      !  code cosmetics
      !
      !  Revision 1.9  2005/09/12 13:30:19  polasekb
      !  added ppm_subid
      !
      !  Revision 1.8  2005/08/11 15:13:33  polasekb
      !  added maxboxcost
      !
      !  Revision 1.7  2005/08/08 13:34:44  polasekb
      !  removed fmm_prec
      !
      !  Revision 1.6  2005/08/04 16:02:46  polasekb
      !  addes some new data
      !
      !  Revision 1.5  2005/07/29 12:36:51  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.4  2005/07/27 21:11:54  polasekb
      !  added totalmass (again)
      !
      !  Revision 1.3  2005/07/25 14:28:32  polasekb
      !  added some constants for the spherical harmonics
      !
      !  Revision 1.2  2005/06/02 13:55:16  polasekb
      !  removed totalmass
      !
      !  Revision 1.1  2005/05/27 08:04:09  polasekb
      !  initial implementation
      !
      !  
      !	 Revision 0 2004/11/11 4:04:15 polasekb
      !  Start.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

MODULE ppm_module_data_fmm   
  !--------------------------------------------------------------------------
  !Modules
  !-----------------------------------------------------------------------------
   USE ppm_module_data,ONLY:ppm_kind_single,ppm_kind_double
   PRIVATE :: ppm_kind_single,ppm_kind_double

  !-----------------------------------------------------------------------------
  ! Defining radius of tree boxes
  !-----------------------------------------------------------------------------

  REAL(ppm_kind_single),DIMENSION(:),POINTER     :: radius_s
  REAL(ppm_kind_double),DIMENSION(:),POINTER     :: radius_d

  !-----------------------------------------------------------------------------
  ! Defining expansions of all tree boxes
  ! 1st index: boxid
  ! 2nd/3rd index: expansion
  !-----------------------------------------------------------------------------

  COMPLEX(ppm_kind_single),DIMENSION(:,:,:),POINTER :: expansion_s
  COMPLEX(ppm_kind_double),DIMENSION(:,:,:),POINTER :: expansion_d

  !-----------------------------------------------------------------------------
  ! Defining center of mass of tree boxes
  !-----------------------------------------------------------------------------

  REAL(ppm_kind_single),DIMENSION(:,:),POINTER   :: centerofbox_s
  REAL(ppm_kind_double),DIMENSION(:,:),POINTER   :: centerofbox_d

  !-----------------------------------------------------------------------------
  ! Defining totalmass of tree boxes
  !-----------------------------------------------------------------------------

  REAL(ppm_kind_single),DIMENSION(:),POINTER     :: totalmass_s
  REAL(ppm_kind_double),DIMENSION(:),POINTER     :: totalmass_d

  !---------------------------------------------------------------------------
  ! Store tree output in data file
  !--------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Defining min_box, minimum extent of tree boxes
  !----------------------------------------------------------------------------

  REAL(ppm_kind_single),DIMENSION(:,:),POINTER   :: min_box_s
  REAL(ppm_kind_double),DIMENSION(:,:),POINTER   :: min_box_d


  !----------------------------------------------------------------------------
  ! Defining max_box, maximum extent of tree boxes
  !----------------------------------------------------------------------------

  REAL(ppm_kind_single),DIMENSION(:,:),POINTER   :: max_box_s
  REAL(ppm_kind_double),DIMENSION(:,:),POINTER   :: max_box_d


  !----------------------------------------------------------------------------
  ! Defining nbox, total number of boxes
  !----------------------------------------------------------------------------

  INTEGER                                        :: nbox 


  !----------------------------------------------------------------------------
  ! Defining nchld, number of children of the box
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:),POINTER                  :: nchld

  !----------------------------------------------------------------------------
  ! Defining lhbx, pointer to first and last point in box (in lpdx)
  ! 1st index: 1 or 2, first and last
  ! 2nd index: box id
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:,:),POINTER               :: lhbx

  !----------------------------------------------------------------------------
  ! Defining lpdx, permutation of xp, particles ordered according to tree
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:),POINTER                :: lpdx

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Defining boxcost
  !----------------------------------------------------------------------------

  REAL(ppm_kind_single),DIMENSION(:),POINTER   :: boxcost_s
  REAL(ppm_kind_double),DIMENSION(:),POINTER   :: boxcost_d


  !----------------------------------------------------------------------------
  ! Defining parent, the partent of the box
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:),POINTER                 :: parent


  !----------------------------------------------------------------------------
  ! Defining child, the child ids of the box
  ! 1st index: child number (1-8 in octtree)
  ! 2nd index: box id
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:,:),POINTER               :: child


  !----------------------------------------------------------------------------
  ! Defining blevel, level of box
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:),POINTER                 :: blevel


  !----------------------------------------------------------------------------
  ! Defining nbpl, number of boxes per level
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:),POINTER                 :: nbpl


  !----------------------------------------------------------------------------
  ! Defining nlevel, total number of levels
  !----------------------------------------------------------------------------

  INTEGER                                       :: nlevel

  !---------------------------------------------------------------------------
  ! End Tree data
  !--------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Defining boxid, box id of sub id
  ! 1st index: sub id
  ! 2nd index: user topology id
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:,:),POINTER               :: ppm_boxid

  !----------------------------------------------------------------------------
  ! Defining subid. sub id of box id
  ! 1st index: box id
  ! 2nd index: user topology id
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:,:),POINTER               :: ppm_subid

  !----------------------------------------------------------------------------
  ! Defining boxpart, which particle is in which box
  !----------------------------------------------------------------------------

  INTEGER,DIMENSION(:),POINTER               :: boxpart

  !----------------------------------------------------------------------------
  ! Defining maxboxcost, the maximum nr of particles per box
  !----------------------------------------------------------------------------

  REAL(ppm_kind_single)                      :: maxboxcost_s
  REAL(ppm_kind_double)                      :: maxboxcost_d

  !----------------------------------------------------------------------------
  ! Storing Info for spherical harmonics 
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Defining Anm
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:,:),POINTER :: Anm_s
  REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: Anm_d

  !----------------------------------------------------------------------------
  ! Defining sqrtfac
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:,:),POINTER :: sqrtfac_s
  REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: sqrtfac_d

  !----------------------------------------------------------------------------
  ! Defining Cnm
  !----------------------------------------------------------------------------
  COMPLEX(ppm_kind_single),DIMENSION(:,:),POINTER :: Cnm_s
  COMPLEX(ppm_kind_double),DIMENSION(:,:),POINTER :: Cnm_d

  !----------------------------------------------------------------------------
  ! Defining Inner
  !----------------------------------------------------------------------------
  COMPLEX(ppm_kind_single),DIMENSION(:,:),POINTER :: Inner_s
  COMPLEX(ppm_kind_double),DIMENSION(:,:),POINTER :: Inner_d
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Defining Ynm
  !----------------------------------------------------------------------------
  COMPLEX(ppm_kind_single),DIMENSION(:,:),POINTER :: Ynm_s
  COMPLEX(ppm_kind_double),DIMENSION(:,:),POINTER :: Ynm_d

  !----------------------------------------------------------------------------
  ! Defining Pnm
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:,:),POINTER :: Pnm_s
  REAL(ppm_kind_double),DIMENSION(:,:),POINTER :: Pnm_d

  !----------------------------------------------------------------------------
  ! Defining fracfac
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: fracfac_s
  REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: fracfac_d

  !----------------------------------------------------------------------------
  ! Defining rho
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: rho_s
  REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: rho_d

  !----------------------------------------------------------------------------
  ! Defining theta
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: theta_s
  REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: theta_d

  !----------------------------------------------------------------------------
  ! Defining phi
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: phi_s
  REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: phi_d

  !----------------------------------------------------------------------------
  ! Defining fac
  !----------------------------------------------------------------------------
  REAL(ppm_kind_single),DIMENSION(:  ),POINTER :: fac_s
  REAL(ppm_kind_double),DIMENSION(:  ),POINTER :: fac_d


END MODULE ppm_module_data_fmm
