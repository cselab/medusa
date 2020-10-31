      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_gmm_extend_fwd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the 
      !                 forward extending routine of the marching method.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_gmm_extend_fwd.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2005/05/10 04:41:17  ivos
      !  Newly created during modularization of ppm_gmm_march. Marching
      !  and orthogonal extendion are now in separate routines for faster
      !  compilation.
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
#define __2D                       3
#define __3D                       4

      MODULE ppm_module_gmm_extend_fwd

         !----------------------------------------------------------------------
         !  Define interfaces to ppm_gmm_extend_fwd
         !----------------------------------------------------------------------
         INTERFACE ppm_gmm_extend_fwd
            MODULE PROCEDURE ppm_gmm_extend_fwd_2ds
            MODULE PROCEDURE ppm_gmm_extend_fwd_2dd
            MODULE PROCEDURE ppm_gmm_extend_fwd_3ds
            MODULE PROCEDURE ppm_gmm_extend_fwd_3dd
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM __2D
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend_fwd.inc"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend_fwd.inc"
#undef __KIND
#undef __DIM

#define __DIM __3D
#define __KIND __SINGLE_PRECISION
#include "ppm_gmm_extend_fwd.inc"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_gmm_extend_fwd.inc"
#undef __KIND
#undef __DIM

      END MODULE ppm_module_gmm_extend_fwd
