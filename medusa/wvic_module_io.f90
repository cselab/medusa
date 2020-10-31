!-------------------------------------------------------------------------------
!* filename: wvic_module_io                                                   *!
!* project : ppm                                                              *!
!* purpose : module for the IO of WVIC client for ppm                         *!
!*         :                                                                  *!
!* author  : Philippe Chatelain                                               *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 10 14:12:59 2004                                         *!
!
!  $Log: wvic_module_io.F,v $
!  Revision 1.2  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.1  2006/09/01 15:45:19  pchatela
!  IO module stuff
!
!
!
!-------------------------------------------------------------------------------

MODULE wvic_module_io
	
    !---------------------------------------------------------------------------
    ! precision_tvk set the precision for vtk output
    ! my experience is that this can either be 4 or 8
    !---------------------------------------------------------------------------
    INTEGER, PARAMETER       :: precision_vtk = 4

    INTERFACE wvic_io_dump_field_netcdf
        MODULE PROCEDURE wvic_io_dump_vector_netcdf
        MODULE PROCEDURE wvic_io_dump_vectorc_netcdf
        MODULE PROCEDURE wvic_io_dump_scalar_netcdf
    END INTERFACE
    INTERFACE wvic_io_dump_field_vtk
        MODULE PROCEDURE wvic_io_dump_vector_vtk
        MODULE PROCEDURE wvic_io_dump_scalar_vtk
    END INTERFACE
    
CONTAINS

#define __REAL 1
#define __COMPLEX 2

#include "wvic_io_dump_scalar_netcdf.inc"

#define __WHAT __COMPLEX
#include "wvic_io_dump_vector_netcdf.inc"
#undef __WHAT
#define __WHAT __REAL
#include "wvic_io_dump_vector_netcdf.inc"
#undef __WHAT

#include "wvic_io_dump_scalar_vtk.inc"
#include "wvic_io_dump_vector_vtk.inc"

END MODULE wvic_module_io
