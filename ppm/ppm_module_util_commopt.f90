      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_util_commopt
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
      !  $Log: ppm_module_util_commopt.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 07:30:11  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_util_commopt

         !----------------------------------------------------------------------
         !  Define interface to the communication optimization routine
         !----------------------------------------------------------------------
        ! here I do some hard core hacking
        ! Thu Apr 27 15:03:38 CEST 2006

        INTERFACE ppm_util_commopt
           MODULE PROCEDURE ppm_util_commopt_cart
        END INTERFACE
         !         INTERFACE ppm_util_commopt_cart_somethingelse
!            MODULE PROCEDURE ppm_util_commopt_somethingelse
!         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_util_commopt_cart.inc"
!#include "ppm_util_commopt.f"

      END MODULE ppm_module_util_commopt
