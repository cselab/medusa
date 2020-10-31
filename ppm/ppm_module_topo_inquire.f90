      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_topo_inquire
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 callable from the outside. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_topo_inquire.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2004/07/26 07:30:10  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_topo_inquire

         !----------------------------------------------------------------------
         !  Define interface to topology inquire routine
         !----------------------------------------------------------------------
         INTERFACE ppm_topo_inquire
            MODULE PROCEDURE ppm_topo_inquire
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_topo_inquire.inc"

      END MODULE ppm_module_topo_inquire
