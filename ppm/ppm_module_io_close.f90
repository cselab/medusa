      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_io_close
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
      !  $Log: ppm_module_io_close.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 07:29:38  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_io_close

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_close
         !----------------------------------------------------------------------
         INTERFACE ppm_io_close
             MODULE PROCEDURE ppm_io_close
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_io_close.f"

      END MODULE ppm_module_io_close
