      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_tree_cutpos
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 routine that determines the optimal position where
      !                 to cut a box.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_tree_cutpos.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/09/22 10:32:09  ivos
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

      MODULE ppm_module_tree_cutpos

         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_cutpos
            MODULE PROCEDURE ppm_tree_cutpos_s
            MODULE PROCEDURE ppm_tree_cutpos_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_tree_cutpos.inc"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_tree_cutpos.inc"
#undef __KIND

      END MODULE ppm_module_tree_cutpos
