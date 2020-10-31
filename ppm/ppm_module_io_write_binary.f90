      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_io_write_binary
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the IO routines.
      !                 It also included those routines and provides
      !                 INTERFACEs.
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_io_write_binary.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 07:29:42  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
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

      MODULE ppm_module_io_write_binary

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_write_binary
         !----------------------------------------------------------------------
         INTERFACE ppm_io_write_binary
             MODULE PROCEDURE ppm_io_write_binarys
             MODULE PROCEDURE ppm_io_write_binaryd
             MODULE PROCEDURE ppm_io_write_binaryi
             MODULE PROCEDURE ppm_io_write_binaryl
             MODULE PROCEDURE ppm_io_write_binarysc
             MODULE PROCEDURE ppm_io_write_binarydc
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_io_write_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_io_write_binary.f"
#undef __KIND
#define __KIND __INTEGER
#include "ppm_io_write_binary.f"
#undef __KIND
#define __KIND __LOGICAL
#include "ppm_io_write_binary.f"
#undef __KIND
#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_io_write_binary.f"
#undef __KIND
#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_io_write_binary.f"
#undef __KIND

      END MODULE ppm_module_io_write_binary
