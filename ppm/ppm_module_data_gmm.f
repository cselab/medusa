      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_data_gmm
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the group 
      !                 marching method routines.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_gmm.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.4  2005/04/27 01:06:14  ivos
      !  Convergence tests completed, cleaned up code, optmized code (Shark),
      !  and changed structure to allow faster compilation.
      !
      !  Revision 1.3  2005/04/21 04:48:25  ivos
      !  Cleaned interfaces and removed unnecessary overloaded versions.
      !
      !  Revision 1.2  2005/03/16 06:20:10  ivos
      !  Several bugfixes. 1st order version is now tested. Moved all large
      !  data to the module.
      !
      !  Revision 1.1  2005/03/10 01:37:17  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_data_gmm

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_single, ppm_kind_double 
         PRIVATE :: ppm_kind_single, ppm_kind_double
         
         !----------------------------------------------------------------------
         !  Parameters
         !----------------------------------------------------------------------
         ! It is important that far is 0 and far<close<accepted
         ! Otherwise the computation of the switches in gmm_slvupwnd needs to
         ! be adjusted !
         INTEGER, PARAMETER      :: ppm_gmm_param_far      = 0
         INTEGER, PARAMETER      :: ppm_gmm_param_close    = 1
         INTEGER, PARAMETER      :: ppm_gmm_param_accepted = 2

         !----------------------------------------------------------------------
         !  Global values
         !----------------------------------------------------------------------
         INTEGER                  :: incr,maxxhi,maxyhi,maxzhi
         INTEGER                  :: gmm_topoid,gmm_meshid
         
         !----------------------------------------------------------------------
         !  Sparse matrix structure as workspace
         !----------------------------------------------------------------------
         INTEGER              , DIMENSION(:,:)  , POINTER   :: gmm_ipos
         REAL(ppm_kind_double), DIMENSION(:  )  , POINTER   :: gmm_phid
         REAL(ppm_kind_single), DIMENSION(:  )  , POINTER   :: gmm_phis
         REAL(ppm_kind_double), DIMENSION(:,:)  , POINTER   :: gmm_clod
         REAL(ppm_kind_single), DIMENSION(:,:)  , POINTER   :: gmm_clos
         REAL(ppm_kind_double), DIMENSION(:,:)  , POINTER   :: gmm_clod2
         REAL(ppm_kind_single), DIMENSION(:,:)  , POINTER   :: gmm_clos2
         INTEGER                                            :: gmm_lsiz = -1
         INTEGER              , DIMENSION(:,:,:,:), POINTER :: gmm_state3d
         INTEGER              , DIMENSION(:,:,:), POINTER   :: gmm_state2d
         INTEGER              , DIMENSION(:,:)  , POINTER   :: iptstmp
         INTEGER              , DIMENSION(:,:)  , POINTER   :: iptstmp2
         INTEGER              , DIMENSION(:  )  , POINTER   :: idx,key

      END MODULE ppm_module_data_gmm
