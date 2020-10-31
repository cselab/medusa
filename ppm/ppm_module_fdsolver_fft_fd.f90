#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_fdsolver_fft_fd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the fft 
      !                  routines
      !                 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fdsolver_fft_fd.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2005/02/16 12:07:09  hiebers
      !  initial implementation
      !
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
#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __SINGLE_PRECISION_COMPLEX   5
#define __DOUBLE_PRECISION_COMPLEX   6
#define __SINGLE_PRECISION_COMPLEX_Z 7
#define __DOUBLE_PRECISION_COMPLEX_Z 8

#define __SLAB                      10

      MODULE ppm_module_fdsolver_fft_fd

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_fft_fd
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_fft_fd_slab
            MODULE PROCEDURE ppm_fdsolver_fft_fd_slab_3ds
            MODULE PROCEDURE ppm_fdsolver_fft_fd_slab_3dd
         END INTERFACE

         INTERFACE ppm_fdsolver_fft_fd
            MODULE PROCEDURE ppm_fdsolver_fft_fd_2ds
            MODULE PROCEDURE ppm_fdsolver_fft_fd_2dd
            MODULE PROCEDURE ppm_fdsolver_fft_fd_2dc
            MODULE PROCEDURE ppm_fdsolver_fft_fd_2dcc

            MODULE PROCEDURE ppm_fdsolver_fft_fd_3ds
            MODULE PROCEDURE ppm_fdsolver_fft_fd_3dd
            MODULE PROCEDURE ppm_fdsolver_fft_fd_3dc
            MODULE PROCEDURE ppm_fdsolver_fft_fd_3dcc
         END INTERFACE

         INTERFACE ppm_fdsolver_fft_fd_z
            MODULE PROCEDURE ppm_fdsolver_fft_fd_z_3dc
            MODULE PROCEDURE ppm_fdsolver_fft_fd_z_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __CASE __SLAB
 
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND

#undef __CASE

 
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_fft_fd_2d.inc"
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_fft_fd_2d.inc"
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_fdsolver_fft_fd_2d.inc"
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_fdsolver_fft_fd_2d.inc"
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND



#define __KIND __SINGLE_PRECISION_COMPLEX_Z
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX_Z
#include "ppm_fdsolver_fft_fd_3d.inc"
#undef __KIND

      END MODULE ppm_module_fdsolver_fft_fd
