!-------------------------------------------------------------------------------
! Dummy module for ppm_fft_velocities
!-------------------------------------------------------------------------------
! $Log: ppm_module_fft.F,v $
! Revision 1.5  2006/10/03 16:11:58  pchatela
! Added spectra in the z direction (misnamed kx spectrum...)
! Added spatial diagnostics, like kinetic energy, enstrophy, circulation
! as functions of z, dumped at the frequency ndump
!
! Revision 1.4  2006/09/16 00:22:02  pchatela
! Implemented the kinetic energy spectrum, dumped into an ascii file.
!
! Revision 1.3  2006/09/03 14:34:15  pchatela
! Added trivial routine to add U_infty in case of spectral velocities
! Fixed inverse transforms in the fft_velocities_bgw, the inverse transform should depend on the qty considered
!
! Revision 1.2  2006/08/11 12:14:17  menahel
! *** empty log message ***
!
! Revision 1.1.1.1  2006/07/25 15:13:46  menahel
! initial import
!
! Revision 1.1  2005/09/28 11:40:24  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------

MODULE ppm_module_fft

  INTEGER :: sx234992_nothing

CONTAINS

#define __SINGLE_PRECISION 1
#define __VELOCITIES 2
#define __POTENTIAL 3
#define __KESPECTRUM 4
#define __KESPECTRUM_KX 5
#define __SOLENOIDAL 6


#define __KIND __SINGLE_PRECISION


#define __WHAT __VELOCITIES
#include "ppm_fft_velocities.F"
#undef __WHAT
#define __WHAT __POTENTIAL
#include "ppm_fft_velocities.F"
#undef __WHAT
#define __WHAT __VELOCITIES
#include "ppm_fft_velocities_bgw.F"
#undef __WHAT

#define __WHAT __KESPECTRUM
#include "ppm_fft_diagnostics.F"
#undef __WHAT

#define __WHAT __KESPECTRUM_KX
#include "ppm_fft_diagnostics.F"
#undef __WHAT

#define __WHAT __SOLENOIDAL
#include "ppm_fft_velocities.F"
#undef __WHAT

#include "ppm_fft_killdiv.F"

END MODULE ppm_module_fft
