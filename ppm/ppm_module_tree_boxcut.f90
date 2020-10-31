      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_tree_boxcut
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 routine that subdivides tree boxes.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_tree_boxcut.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/09/22 10:32:08  ivos
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

      MODULE ppm_module_tree_boxcut

         !----------------------------------------------------------------------
         !  Define interfaces to the routine(s)
         !----------------------------------------------------------------------
         INTERFACE ppm_tree_boxcut
            MODULE PROCEDURE ppm_tree_boxcut_s
            MODULE PROCEDURE ppm_tree_boxcut_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_tree_boxcut.inc"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_tree_boxcut.inc"
#undef __KIND

      END MODULE ppm_module_tree_boxcut
