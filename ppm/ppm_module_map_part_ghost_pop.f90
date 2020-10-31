      !-------------------------------------------------------------------------
      !  Module       :          ppm_module_map_part_ghost_pop
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
      !  $Log: ppm_module_map_part_ghost_pop.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2006/02/03 09:38:23  ivos
      !  Newly added the PRELIMINARY version of the ghost_pop routines. They
      !  still need major clean-up.
      !
      !  Revision 1.1  2004/07/26 07:29:54  ivos
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
#define __INTEGER                  3
#define __LOGICAL                  4
#define __SINGLE_PRECISION_COMPLEX 5
#define __DOUBLE_PRECISION_COMPLEX 6

      MODULE ppm_module_map_part_ghost_pop

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_map_part_ghost_pop
         !----------------------------------------------------------------------
         INTERFACE ppm_map_part_ghost_pop
            ! scalar (1d) particle data
            MODULE PROCEDURE ppm_map_part_ghost_pop_1dd
            MODULE PROCEDURE ppm_map_part_ghost_pop_1ds
            MODULE PROCEDURE ppm_map_part_ghost_pop_1di
            MODULE PROCEDURE ppm_map_part_ghost_pop_1dl
            MODULE PROCEDURE ppm_map_part_ghost_pop_1ddc
            MODULE PROCEDURE ppm_map_part_ghost_pop_1dsc

            ! vector (2d) particle data
            MODULE PROCEDURE ppm_map_part_ghost_pop_2dd
            MODULE PROCEDURE ppm_map_part_ghost_pop_2ds
            MODULE PROCEDURE ppm_map_part_ghost_pop_2di
            MODULE PROCEDURE ppm_map_part_ghost_pop_2dl
            MODULE PROCEDURE ppm_map_part_ghost_pop_2ddc
            MODULE PROCEDURE ppm_map_part_ghost_pop_2dsc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM 1
#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#undef __DIM

#define __DIM 2
#define __KIND __SINGLE_PRECISION
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_map_part_ghost_pop.f"
#undef __KIND
#undef __DIM

      END MODULE ppm_module_map_part_ghost_pop
