      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_util_rank
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_rank.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 07:30:13  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
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

      MODULE ppm_module_util_rank

         !----------------------------------------------------------------------
         !  Define interfaces to the cell list ranking routines
         !----------------------------------------------------------------------
         INTERFACE ppm_util_rank3d
            MODULE PROCEDURE ppm_util_rank3d_s
            MODULE PROCEDURE ppm_util_rank3d_d
         END INTERFACE
         INTERFACE ppm_util_rank2d
            MODULE PROCEDURE ppm_util_rank2d_s
            MODULE PROCEDURE ppm_util_rank2d_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __KIND __SINGLE_PRECISION
#include "ppm_util_rank3d.inc"
#include "ppm_util_rank2d.inc"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_rank3d.inc"
#include "ppm_util_rank2d.inc"
#undef __KIND

      END MODULE ppm_module_util_rank
