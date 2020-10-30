      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_create_ode
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_create_ode
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_create_ode.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.2  2004/07/26 11:20:29  michaebe
      !  added dummy interface
      !
      !  Revision 1.1  2004/07/26 07:45:47  michaebe
      !  Procedure modules created in the course of atomization.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2


      MODULE ppm_module_ode_create_ode

        !-----------------------------------------------------
        !  Dummy interface
        !-----------------------------------------------------
        INTERFACE ppm_ode_create_ode
           MODULE PROCEDURE ppm_ode_create_ode
        END INTERFACE
        
      CONTAINS

#include "ppm_ode_create_ode.f"

      END MODULE ppm_module_ode_create_ode
