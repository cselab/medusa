      !-------------------------------------------------------------------------
      !  Module       :                  ppm_module_ode
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the user-callable functions of
      !                 the ode time integrator.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.3  2004/07/26 13:41:34  ivos
      !  Initial implementation. These are meta-modules for the user-
      !  callable functions. Only these modules will be given away to
      !  the user.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_ode

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_ode_alldone
         USE ppm_module_ode_create_ode
         USE ppm_module_ode_finalize
         USE ppm_module_ode_init
         USE ppm_module_ode_step
         USE ppm_module_ode_start
         
      END MODULE ppm_module_ode
