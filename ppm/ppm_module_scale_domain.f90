      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_scale_domain
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_scale_domain.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 14:26:48  ivos
      !  Initial implementation.
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

      MODULE ppm_module_scale_domain

         !----------------------------------------------------------------------
         !  Define interfaces to the domain scaling routine
         !----------------------------------------------------------------------
         INTERFACE ppm_scale_domain
            MODULE PROCEDURE ppm_scale_domain_s
            MODULE PROCEDURE ppm_scale_domain_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_scale_domain.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_scale_domain.f"
#undef __KIND

      END MODULE ppm_module_scale_domain
