      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_mesh_alloc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the mesh routines.
      !                 It also included those routines and provides
      !                 INTERFACEs.
      !                
      !  Remarks      : The terminology distinguishes between meshes and
      !                 fields (the data living on the meshes). Several
      !                 fields can use the same mesh. Meshes are defined as
      !                 ppm-internal TYPES, whereas fields are
      !                 user-provided arrays.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_mesh_alloc.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 07:29:57  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_mesh_alloc

         !----------------------------------------------------------------------
         !  Define interface to ppm_mesh_alloc
         !----------------------------------------------------------------------
         INTERFACE ppm_mesh_alloc
             MODULE PROCEDURE ppm_mesh_alloc_list
             MODULE PROCEDURE ppm_mesh_alloc_equi
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_mesh_alloc_list.f"
#include "ppm_mesh_alloc_equi.f"

      END MODULE ppm_module_mesh_alloc
