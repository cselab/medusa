      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_topo
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all user-callable routines
      !                 needed to create and manage ppm topologies.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_topo.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.2  2004/07/26 14:26:30  ivos
      !  Added ppm_module_scale_domain.
      !
      !  Revision 1.1  2004/07/26 13:40:32  ivos
      !  Initial implementation. These are meta-modules for the user-
      !  callable functions. Only these mod files will be given away
      !  to the user.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_topo

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_mktopo
         USE ppm_module_topo_check
         USE ppm_module_topo_inquire
         USE ppm_module_mesh_define
         USE ppm_module_scale_domain
         
      END MODULE ppm_module_topo
