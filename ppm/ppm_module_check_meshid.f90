      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_check_meshid
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routine which checks mesh IDs for their validity.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_check_meshid.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/08/31 12:13:46  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_check_meshid

         !----------------------------------------------------------------------
         !  Define interface to the meshid check routine
         !----------------------------------------------------------------------
         INTERFACE ppm_check_meshid
            MODULE PROCEDURE ppm_check_meshid
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_check_meshid.inc"

      END MODULE ppm_module_check_meshid
