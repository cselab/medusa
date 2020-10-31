!------------------------------------------------------------------------------!
!
!
!
!
!
!
!
!
!------------------------------------------------------------------------------!
!  project:  WVIC
!  purpose:  Hyperviscosity model for wvic
!   author:  Michael Bergdorf,
!         :  ETHZ Computational Science & Engineering Laboratory
!    email:  bergdorf@inf.ethz.ch
!------------------------------------------------------------------------------!
!    $Log: wvic_les_hypervisc.F,v $
!    Revision 1.3  2006/10/23 14:38:27  pchatela
!    Corrected preprocessor endif
!
!    Revision 1.2  2006/10/23 14:28:03  pchatela
!    Added parameters for hyperviscosity model
!    Fixed double Laplacian
!
!    Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!    initial import
!
!    Revision 1.1  2005/12/10 21:05:59  michaebe
!    implementation initiale
!
!------------------------------------------------------------------------------!

SUBROUTINE wvic_les_hypervisc

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  REAL(mk) :: fac0
  INTEGER  :: i,j,k,isub,isubl,info
  CHARACTER(len=256) :: msg

  info  = 0
  !les_hypervisc_tzero = 16.0_mk
  fac0  = -les_hypervisc_C/les_hypervisc_tzero

  IF(MINVAL(ghostsize).LT.2) THEN
     WRITE(msg,*) 'Hyperviscous scheme needs ghostsize >= 2'
     CALL ppm_write(rank,'wvic_les_hypervisc',msg,info)
     GOTO 9999
  END IF

  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
#ifndef __OLD_NABLA4
              field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub) + fac0 *( &
                   & field_wp(1,i - 2, j, k, isub) + field_wp(1,i + 2, j, k, isub) + &
                   & field_wp(1,i, j - 2, k, isub) + field_wp(1,i, j + 2, k, isub) + &
                   & field_wp(1,i, j, k - 2, isub) + field_wp(1,i, j, k + 2, isub) + &
                   & 2.0_mk*( field_wp(1,i - 1, j, k + 1, isub) + field_wp(1,i + 1, j, k + 1, isub) + &
                   &          field_wp(1,i, j - 1, k + 1, isub) + field_wp(1,i, j + 1, k + 1, isub) + &
                   &          field_wp(1,i - 1, j, k - 1, isub) + field_wp(1,i + 1, j, k - 1, isub) + &
                   &          field_wp(1,i, j - 1, k - 1, isub) + field_wp(1,i, j + 1, k - 1, isub) + &
                   &          field_wp(1,i - 1, j - 1, k, isub) + field_wp(1,i + 1, j - 1, k, isub) + &
                   &          field_wp(1,i - 1, j + 1, k, isub) + field_wp(1,i + 1, j + 1, k, isub))- &
                   &12.0_mk*( field_wp(1,i - 1, j, k, isub) + field_wp(1,i + 1, j, k, isub) + &
                   &          field_wp(1,i, j - 1, k, isub) + field_wp(1,i, j + 1, k, isub) + &
                   &          field_wp(1,i, j, k - 1, isub) + field_wp(1,i, j, k + 1, isub))+ &
                   &42.0_mk* field_wp(1, i, j, k, isub) )
                   
              field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub) + fac0 *( &
                   & field_wp(2,i - 2, j, k, isub) + field_wp(2,i + 2, j, k, isub) + &
                   & field_wp(2,i, j - 2, k, isub) + field_wp(2,i, j + 2, k, isub) + &
                   & field_wp(2,i, j, k - 2, isub) + field_wp(2,i, j, k + 2, isub) + &
                   & 2.0_mk*( field_wp(2,i - 1, j, k + 1, isub) + field_wp(2,i + 1, j, k + 1, isub) + &
                   &          field_wp(2,i, j - 1, k + 1, isub) + field_wp(2,i, j + 1, k + 1, isub) + &
                   &          field_wp(2,i - 1, j, k - 1, isub) + field_wp(2,i + 1, j, k - 1, isub) + &
                   &          field_wp(2,i, j - 1, k - 1, isub) + field_wp(2,i, j + 1, k - 1, isub) + &
                   &          field_wp(2,i - 1, j - 1, k, isub) + field_wp(2,i + 1, j - 1, k, isub) + &
                   &          field_wp(2,i - 1, j + 1, k, isub) + field_wp(2,i + 1, j + 1, k, isub))- &
                   &12.0_mk*( field_wp(2,i - 1, j, k, isub) + field_wp(2,i + 1, j, k, isub) + &
                   &          field_wp(2,i, j - 1, k, isub) + field_wp(2,i, j + 1, k, isub) + &
                   &          field_wp(2,i, j, k - 1, isub) + field_wp(2,i, j, k + 1, isub))+ &
                   &42.0_mk* field_wp(2, i, j, k, isub) )
                   
              field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub) + fac0 *( &
                   & field_wp(3,i - 2, j, k, isub) + field_wp(3,i + 2, j, k, isub) + &
                   & field_wp(3,i, j - 2, k, isub) + field_wp(3,i, j + 2, k, isub) + &
                   & field_wp(3,i, j, k - 2, isub) + field_wp(3,i, j, k + 2, isub) + &
                   & 2.0_mk*( field_wp(3,i - 1, j, k + 1, isub) + field_wp(3,i + 1, j, k + 1, isub) + &
                   &          field_wp(3,i, j - 1, k + 1, isub) + field_wp(3,i, j + 1, k + 1, isub) + &
                   &          field_wp(3,i - 1, j, k - 1, isub) + field_wp(3,i + 1, j, k - 1, isub) + &
                   &          field_wp(3,i, j - 1, k - 1, isub) + field_wp(3,i, j + 1, k - 1, isub) + &
                   &          field_wp(3,i - 1, j - 1, k, isub) + field_wp(3,i + 1, j - 1, k, isub) + &
                   &          field_wp(3,i - 1, j + 1, k, isub) + field_wp(3,i + 1, j + 1, k, isub))- &
                   &12.0_mk*( field_wp(3,i - 1, j, k, isub) + field_wp(3,i + 1, j, k, isub) + &
                   &          field_wp(3,i, j - 1, k, isub) + field_wp(3,i, j + 1, k, isub) + &
                   &          field_wp(3,i, j, k - 1, isub) + field_wp(3,i, j, k + 1, isub))+ &
                   &42.0_mk* field_wp(3, i, j, k, isub) )
#else
              field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub) + fac0&
                   &*(field_wp(1,i - 2, j, k, isub) - 8.0_mk * field_wp(1,i&
                   & - 1, j, k, isub) + 14.0_mk * field_wp(1,i, j, k, isub) +&
                   & field_wp(1,i - 1, j + 1, k, isub) + field_wp(1,i - 1,&
                   & j - 1, k, isub) + field_wp(1,i - 1, j, k + 1, isub) +&
                   & field_wp(1,i - 1, j, k - 1, isub) - 8.0_mk *&
                   & field_wp(1,i + 1, j, k, isub) - 2.0_mk * field_wp(1,i,&
                   & j + 1, k, isub) - 2.0_mk * field_wp(1,i, j - 1, k, isub)&
                   & - 2.0_mk * field_wp(1,i, j, k + 1,isub) - 2.0_mk *&
                   & field_wp(1,i, j, k - 1, isub) + field_wp(1,i + 2, j, k&
                   &, isub) + field_wp(1,i + 1, j + 1, k, isub) +&
                   & field_wp(1,i + 1, j - 1, k, isub) + field_wp(1,i + 1,&
                   & j, k + 1, isub) + field_wp(1,i + 1, j, k - 1, isub))
              
              field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub) + fac0&
                   &*(field_wp(2,i - 2, j, k, isub) - 8.0_mk * field_wp(2,i&
                   & - 1, j, k, isub) + 14.0_mk * field_wp(2,i, j, k, isub) +&
                   & field_wp(2,i - 1, j + 1, k, isub) + field_wp(2,i - 1,&
                   & j - 1, k, isub) + field_wp(2,i - 1, j, k + 1, isub) +&
                   & field_wp(2,i - 1, j, k - 1, isub) - 8.0_mk *&
                   & field_wp(2,i + 1, j, k, isub) - 2.0_mk * field_wp(2,i,&
                   & j + 1, k, isub) - 2.0_mk * field_wp(2,i, j - 1, k, isub)&
                   & - 2.0_mk * field_wp(2,i, j, k + 1,isub) - 2.0_mk *&
                   & field_wp(2,i, j, k - 1, isub) + field_wp(2,i + 2, j, k&
                   &, isub) + field_wp(2,i + 1, j + 1, k, isub) +&
                   & field_wp(2,i + 1, j - 1, k, isub) + field_wp(2,i + 1,&
                   & j, k + 1, isub) + field_wp(2,i + 1, j, k - 1, isub))
              
              field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub) + fac0&
                   &*(field_wp(3,i - 2, j, k, isub) - 8.0_mk * field_wp(3,i&
                   & - 1, j, k, isub) + 14.0_mk * field_wp(3,i, j, k, isub) +&
                   & field_wp(3,i - 1, j + 1, k, isub) + field_wp(3,i - 1,&
                   & j - 1, k, isub) + field_wp(3,i - 1, j, k + 1, isub) +&
                   & field_wp(3,i - 1, j, k - 1, isub) - 8.0_mk *&
                   & field_wp(3,i + 1, j, k, isub) - 2.0_mk * field_wp(3,i,&
                   & j + 1, k, isub) - 2.0_mk * field_wp(3,i, j - 1, k, isub)&
                   & - 2.0_mk * field_wp(3,i, j, k + 1,isub) - 2.0_mk *&
                   & field_wp(3,i, j, k - 1, isub) + field_wp(3,i + 2, j, k&
                   &, isub) + field_wp(3,i + 1, j + 1, k, isub) +&
                   & field_wp(3,i + 1, j - 1, k, isub) + field_wp(3,i + 1,&
                   & j, k + 1, isub) + field_wp(3,i + 1, j, k - 1, isub))
#endif
           END DO
        END DO
     END DO
  END DO
  
9999 CONTINUE
END SUBROUTINE wvic_les_hypervisc
