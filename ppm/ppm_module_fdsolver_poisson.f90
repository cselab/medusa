#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_fdsolver_poisson
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for field solver 
      !                        routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __COMPLEX          3
#define __DOUBLE_COMPLEX   4



   MODULE ppm_module_fdsolver_poisson

         INTERFACE ppm_fdsolver_poisson
            MODULE PROCEDURE ppm_fdsolver_poisson_2dc
            MODULE PROCEDURE ppm_fdsolver_poisson_2dcc
            MODULE PROCEDURE ppm_fdsolver_poisson_3dc
            MODULE PROCEDURE ppm_fdsolver_poisson_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __COMPLEX
#include "ppm_fdsolver_poisson_2d.f"
#include "ppm_fdsolver_poisson_3d.f"
#undef  __KIND

#define __KIND __DOUBLE_COMPLEX
#include "ppm_fdsolver_poisson_2d.f"
#include "ppm_fdsolver_poisson_3d.f"
#undef  __KIND



      END MODULE ppm_module_fdsolver_poisson
