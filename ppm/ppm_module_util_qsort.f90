      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_util_qsort
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
      !  $Log: ppm_module_util_qsort.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2005/02/09 16:28:57  pchatela
      !  Initial insertion
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __INTEGER                    3
      
      MODULE ppm_module_util_qsort

         !----------------------------------------------------------------------
         !  Define interface to the list inversion routine
         !----------------------------------------------------------------------
         INTERFACE ppm_util_qsort
            MODULE PROCEDURE ppm_util_qsort_s
            MODULE PROCEDURE ppm_util_qsort_d
            MODULE PROCEDURE ppm_util_qsort_i
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
#define __KIND __SINGLE_PRECISION
#include "ppm_util_qsort.f"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_qsort.f"
#undef __KIND

#define __KIND __INTEGER
#include "ppm_util_qsort.f"
#undef __KIND

      END MODULE ppm_module_util_qsort
