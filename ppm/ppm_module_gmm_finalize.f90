      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_gmm_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the 
      !                 finalization routine of the marching method.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_gmm_finalize.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2005/03/10 01:37:16  ivos
      !  Initial check-in. BEWARE: Not tested in parallel yet!
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_gmm_finalize

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_gmm_finalize
         !----------------------------------------------------------------------
         INTERFACE ppm_gmm_finalize
            MODULE PROCEDURE ppm_gmm_finalize
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
        CONTAINS

#include "ppm_gmm_finalize.inc"

      END MODULE ppm_module_gmm_finalize
