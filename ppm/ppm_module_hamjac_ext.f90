      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_hamjac_ext
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_hamjac_ext
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_hamjac_ext.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2005/07/25 00:34:07  ivos
      !  Initial check-in.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __2D               3
#define __3D               4
#define __VEC              5
#define __SCA              6

      MODULE ppm_module_hamjac_ext

        !-----------------------------------------------------
        !  Interface
        !-----------------------------------------------------
        INTERFACE ppm_hamjac_ext_step

           MODULE PROCEDURE ppm_hamjac_ext_step_3ds
           MODULE PROCEDURE ppm_hamjac_ext_step_3dd
           MODULE PROCEDURE ppm_hamjac_ext_step_3dsv
           MODULE PROCEDURE ppm_hamjac_ext_step_3ddv

        END INTERFACE

        INTERFACE ppm_hamjac_ext

           MODULE PROCEDURE ppm_hamjac_ext_3ds
           MODULE PROCEDURE ppm_hamjac_ext_3dd
           MODULE PROCEDURE ppm_hamjac_ext_3dsv
           MODULE PROCEDURE ppm_hamjac_ext_3ddv

        END INTERFACE


      CONTAINS
#define __DIME  __3D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_step_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_step_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME

#define __DIME  __3D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_ext_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME


#define __DIME  __3D
#define __MODE  __VEC
#define __KIND  __SINGLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_step_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_step_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME

#define __DIME  __3D
#define __MODE  __VEC
#define __KIND  __SINGLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_ext_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME




      END MODULE ppm_module_hamjac_ext
        

        

