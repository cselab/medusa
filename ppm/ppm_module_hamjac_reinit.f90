      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_hamjac_reinit
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_hamjac_reinit
      !                 reinitialization routines.
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_hamjac_reinit.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.1  2005/07/25 00:34:08  ivos
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

      MODULE ppm_module_hamjac_reinit

        !-----------------------------------------------------
        !  Interface
        !-----------------------------------------------------
        INTERFACE ppm_hamjac_reinit_step
           MODULE PROCEDURE ppm_hamjac_reinit_step_3ds
           MODULE PROCEDURE ppm_hamjac_reinit_step_3dd
           MODULE PROCEDURE ppm_hamjac_reinit_step_3dsV
           MODULE PROCEDURE ppm_hamjac_reinit_step_3ddV
           MODULE PROCEDURE ppm_hamjac_reinit_step_2ds
           MODULE PROCEDURE ppm_hamjac_reinit_step_2dd
        END INTERFACE

        INTERFACE ppm_hamjac_reinit_loc_step
           MODULE PROCEDURE ppm_hamjac_reinit_loc_step_3ds
           MODULE PROCEDURE ppm_hamjac_reinit_loc_step_3dd
           MODULE PROCEDURE ppm_hamjac_reinit_loc_step_3dsV
           MODULE PROCEDURE ppm_hamjac_reinit_loc_step_3ddV
        END INTERFACE

        INTERFACE ppm_hamjac_reinit
           MODULE PROCEDURE ppm_hamjac_reinit_3ds
           MODULE PROCEDURE ppm_hamjac_reinit_3dd
           MODULE PROCEDURE ppm_hamjac_reinit_3dsV
           MODULE PROCEDURE ppm_hamjac_reinit_3ddV
           MODULE PROCEDURE ppm_hamjac_reinit_2ds
           MODULE PROCEDURE ppm_hamjac_reinit_2dd
        END INTERFACE

        INTERFACE ppm_hamjac_reinit_loc
           MODULE PROCEDURE ppm_hamjac_reinit_loc_3ds
           MODULE PROCEDURE ppm_hamjac_reinit_loc_3dd
           MODULE PROCEDURE ppm_hamjac_reinit_loc_3dsV
           MODULE PROCEDURE ppm_hamjac_reinit_loc_3ddV
        END INTERFACE

        INTERFACE ppm_hamjac_reinit_step_ref
           MODULE PROCEDURE ppm_hamjac_reinit_step_ref_3ds
           MODULE PROCEDURE ppm_hamjac_reinit_step_ref_3dd
        END INTERFACE

        INTERFACE ppm_hamjac_reinit_ref
           MODULE PROCEDURE ppm_hamjac_reinit_ref_3ds
           MODULE PROCEDURE ppm_hamjac_reinit_ref_3dd
        END INTERFACE

      CONTAINS
#define __DIME  __3D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_reinit_step_3d.inc"
#include "ppm_hamjac_reinit_loc_step_3d.inc"
#include "ppm_hamjac_reinit_step_ref_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_reinit_step_3d.inc"
#include "ppm_hamjac_reinit_loc_step_3d.inc"
#include "ppm_hamjac_reinit_step_ref_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME
#define __DIME  __3D
#define __MODE  __VEC
#define __KIND  __SINGLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_reinit_step_3d.inc"
#include "ppm_hamjac_reinit_loc_step_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_reinit_step_3d.inc"
#include "ppm_hamjac_reinit_loc_step_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME

        
        
#define __DIME  __2D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 2D SCA SINGLE
#include "ppm_hamjac_reinit_step_2d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 2D SCA DOUBLE
#include "ppm_hamjac_reinit_step_2d.inc"
#undef __KIND
#undef __MODE
#undef __DIME

#define __DIME  __3D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 3D SCA SINGLE
#include "ppm_hamjac_reinit_3d.inc"
#include "ppm_hamjac_reinit_loc_3d.inc"
#include "ppm_hamjac_reinit_ref_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D SCA DOUBLE
#include "ppm_hamjac_reinit_3d.inc"
#include "ppm_hamjac_reinit_loc_3d.inc"
#include "ppm_hamjac_reinit_ref_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME

#define __DIME  __3D
#define __MODE  __VEC
#define __KIND  __SINGLE_PRECISION
        ! 3D VEC SINGLE
#include "ppm_hamjac_reinit_3d.inc"
#include "ppm_hamjac_reinit_loc_3d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 3D VEC DOUBLE
#include "ppm_hamjac_reinit_3d.inc"
#include "ppm_hamjac_reinit_loc_3d.inc"
#undef __KIND
#undef __MODE
#undef __DIME


        
#define __DIME  __2D
#define __MODE  __SCA
#define __KIND  __SINGLE_PRECISION
        ! 2D SCA SINGLE
#include "ppm_hamjac_reinit_2d.inc"
#undef __KIND
#define __KIND  __DOUBLE_PRECISION
        ! 2D SCA DOUBLE
#include "ppm_hamjac_reinit_2d.inc"
#undef __KIND
#undef __MODE
#undef __DIME

      END MODULE ppm_module_hamjac_reinit
        

        

