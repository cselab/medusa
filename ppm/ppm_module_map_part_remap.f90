      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_map_part_remap
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the mapping
      !                 routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_map_part_remap.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 07:29:55  ivos
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
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2

      MODULE ppm_module_map_part_remap

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_remap
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_remap
            MODULE PROCEDURE ppm_map_part_remap_d
            MODULE PROCEDURE ppm_map_part_remap_s
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_remap.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_remap.f"
#undef __KIND

      END MODULE ppm_module_map_part_remap
