!------------------------------------------------------------------------------!
!* filename: wvic_dgammadt                                                     *!
!* project : ppm                                                              *!
!* purpose : compute right hand side for the vorticity equation               *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 17:23:11 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!------------------------------------------------------------------------------!
! $Log: wvic_dgammadt.F,v $
! Revision 1.1.1.1  2006/07/25 15:13:46  menahel
! initial import
!
! Revision 1.4  2005/11/21 17:36:37  michaebe
! major session. now let''s debug
!
! Revision 1.3  2005/11/11 17:13:11  michaebe
! removed allcoation calls
!
! Revision 1.2  2005/11/11 14:04:23  michaebe
! clean up, additions, comments
!
! Revision 1.1  2005/09/28 11:40:30  michaebe
! Fork from ppm_pvc
!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  compute laplacian and stretching term for the vorticity balance
!------------------------------------------------------------------------------!
SUBROUTINE wvic_dgammadt (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out)              :: info

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER :: isub, i , j , k, isubl
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6
  REAL(mk) :: tim1s, tim1e
  CHARACTER(len=256) :: msg

  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()

  fac1 = 1.0_mk/dx**2*nu
  fac2 = 1.0_mk/dy**2*nu
  fac3 = 1.0_mk/dz**2*nu

  fac4 = 0.5_mk/dx
  fac5 = 0.5_mk/dy
  fac6 = 0.5_mk/dz

  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)
#ifndef __CONSERVATIVE
              field_dwp(1,i,j,k,isub) = &
                   &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                   &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                   &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(1,i,j,k,isub)          +&
#ifndef __WITHOUT_STRETCHING                  
                   &fac4*(field_wp(1,i,j,k,isub)*( &
                   & field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub)))+    &
                   &fac5*(field_wp(2,i,j,k,isub)*( &
                   & field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub)))+    &
                   &fac6*(field_wp(3,i,j,k,isub)*( &
                   & field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub)))
#else
                   &0.0_mk
#endif                   

              field_dwp(2,i,j,k,isub) = &
                   &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                   &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                   &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(2,i,j,k,isub)          +&
#ifndef __WITHOUT_STRETCHING                  
                   &fac4*(field_wp(1,i,j,k,isub)*( &
                   & field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub)))+    &
                   &fac5*(field_wp(2,i,j,k,isub)*( &
                   & field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub)))+    &
                   &fac6*(field_wp(3,i,j,k,isub)*( &
                   & field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub)))
#else
                   &0.0_mk
#endif                   

              field_dwp(3,i,j,k,isub) = &
                   &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                   &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                   &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(3,i,j,k,isub)          +&
#ifndef __WITHOUT_STRETCHING                  
                   &fac4*(field_wp(1,i,j,k,isub)*( &
                   & field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub)))+    &
                   &fac5*(field_wp(2,i,j,k,isub)*( &
                   & field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub)))+    &
                   &fac6*(field_wp(3,i,j,k,isub)*( &
                   & field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub)))
#else
                   &0.0_mk
#endif                   
#else
#ifndef __REALLYCONSERVATIVE
                   
              field_dwp(1,i,j,k,isub) = &
                   &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                   &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                   &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(1,i,j,k,isub)          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k-1,isub))

              field_dwp(2,i,j,k,isub) = &
                   &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                   &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                   &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(2,i,j,k,isub)          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k-1,isub))

              field_dwp(3,i,j,k,isub) = &
                   &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                   &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                   &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(3,i,j,k,isub)          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k-1,isub))
#else
              field_dwp(1,i,j,k,isub) = &
                   &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                   &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                   &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(1,i,j,k,isub)          +&
                   & fac4*((field_wp(1,i  ,j,k,isub)*field_up(1,i+1,j,k,isub) + &
                   &        field_wp(1,i+1,j,k,isub)*field_up(1,i  ,j,k,isub))- &
                   &       (field_wp(1,i  ,j,k,isub)*field_up(1,i-1,j,k,isub) + &
                   &        field_wp(1,i-1,j,k,isub)*field_up(1,i  ,j,k,isub)))&
                   & + &
                   & fac5*((field_wp(2,i,j  ,k,isub)*field_up(1,i,j+1,k,isub) + &
                   &        field_wp(2,i,j+1,k,isub)*field_up(1,i,j  ,k,isub))- &
                   &       (field_wp(2,i  ,j,k,isub)*field_up(1,i,j-1,k,isub) + &
                   &        field_wp(2,i,j-1,k,isub)*field_up(1,i,j  ,k,isub)))&
                   & + &
                   & fac6*((field_wp(3,i  ,j,k,isub)*field_up(1,i,j,k+1,isub) + &
                   &        field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k  ,isub))- &
                   &       (field_wp(3,i,j,k  ,isub)*field_up(1,i,j,k-1,isub) + &
                   &        field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k  ,isub)))
              field_dwp(2,i,j,k,isub) = &
                   &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                   &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                   &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(2,i,j,k,isub)          +&
                   & fac4*((field_wp(1,i  ,j,k,isub)*field_up(2,i+1,j,k,isub) + &
                   &        field_wp(1,i+1,j,k,isub)*field_up(2,i  ,j,k,isub))- &
                   &       (field_wp(1,i  ,j,k,isub)*field_up(2,i-1,j,k,isub) + &
                   &        field_wp(1,i-1,j,k,isub)*field_up(2,i  ,j,k,isub)))&
                   & + &
                   & fac5*((field_wp(2,i,j  ,k,isub)*field_up(2,i,j+1,k,isub) + &
                   &        field_wp(2,i,j+1,k,isub)*field_up(2,i,j  ,k,isub))- &
                   &       (field_wp(2,i  ,j,k,isub)*field_up(2,i,j-1,k,isub) + &
                   &        field_wp(2,i,j-1,k,isub)*field_up(2,i,j  ,k,isub)))&
                   & + &
                   & fac6*((field_wp(3,i  ,j,k,isub)*field_up(2,i,j,k+1,isub) + &
                   &        field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k  ,isub))- &
                   &       (field_wp(3,i,j,k  ,isub)*field_up(2,i,j,k-1,isub) + &
                   &        field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k  ,isub)))
              field_dwp(3,i,j,k,isub) = &
                   &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                   &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                   &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(3,i,j,k,isub)          +&
                   & fac4*((field_wp(1,i  ,j,k,isub)*field_up(3,i+1,j,k,isub) + &
                   &        field_wp(1,i+1,j,k,isub)*field_up(3,i  ,j,k,isub))- &
                   &       (field_wp(1,i  ,j,k,isub)*field_up(3,i-1,j,k,isub) + &
                   &        field_wp(1,i-1,j,k,isub)*field_up(3,i  ,j,k,isub)))&
                   & + &
                   & fac5*((field_wp(2,i,j  ,k,isub)*field_up(3,i,j+1,k,isub) + &
                   &        field_wp(2,i,j+1,k,isub)*field_up(3,i,j  ,k,isub))- &
                   &       (field_wp(2,i  ,j,k,isub)*field_up(3,i,j-1,k,isub) + &
                   &        field_wp(2,i,j-1,k,isub)*field_up(3,i,j  ,k,isub)))&
                   & + &
                   & fac6*((field_wp(3,i  ,j,k,isub)*field_up(3,i,j,k+1,isub) + &
                   &        field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k  ,isub))- &
                   &       (field_wp(3,i,j,k  ,isub)*field_up(3,i,j,k-1,isub) + &
                   &        field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k  ,isub)))


              
#endif              
#endif                   
                   
           END DO

        END DO

     END DO

  END DO

  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_dgammadt',msg,info)
  
END SUBROUTINE wvic_dgammadt







!------------------------------------------------------------------------------!
!  compute laplacian and stretching term for the vorticity balance
!  use half a muscle scheme for the laplacian and the pcon scheme for the
!  stretching 
!------------------------------------------------------------------------------!
SUBROUTINE wvic_dgammadt_MUSCLE (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE


  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out)              :: info

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER :: isub, i , j , k, isubl
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9
  REAL(mk) :: tim1s, tim1e
  CHARACTER(len=256) :: msg

  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()

  fac1 = 1.0_mk/dx**2*nu
  fac2 = 1.0_mk/dy**2*nu
  fac3 = 1.0_mk/dz**2*nu

  fac4 = 0.5_mk/dx
  fac5 = 0.5_mk/dy
  fac6 = 0.5_mk/dz

  fac7 = 1.0_mk/(dx**2)*nu
  fac8 = 1.0_mk/(dy**2+dz**2)*nu
  fac9 = 1.0_mk/(dy**2+dz**2)*nu

  
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)

                   
              field_dwp(1,i,j,k,isub) = (&
                   &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                   &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                   &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(1,i,j,k,isub)          +&
                   &fac7*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                   &fac8*(field_wp(1,i,j+1,k+1,isub)+field_wp(1,i,j-1,k-1,isub))+&
                   &fac9*(field_wp(1,i,j+1,k-1,isub)+field_wp(1,i,j-1,k+1,isub))-&
                   &2.0_mk*(fac7+fac8+fac9)*field_wp(1,i,j,k,isub) )*0.5_mk +&
                   & fac4*((field_wp(1,i  ,j,k,isub)*field_up(1,i+1,j,k,isub) + &
                   &        field_wp(1,i+1,j,k,isub)*field_up(1,i  ,j,k,isub))- &
                   &       (field_wp(1,i  ,j,k,isub)*field_up(1,i-1,j,k,isub) + &
                   &        field_wp(1,i-1,j,k,isub)*field_up(1,i  ,j,k,isub)))&
                   & + &
                   & fac5*((field_wp(2,i,j  ,k,isub)*field_up(1,i,j+1,k,isub) + &
                   &        field_wp(2,i,j+1,k,isub)*field_up(1,i,j  ,k,isub))- &
                   &       (field_wp(2,i  ,j,k,isub)*field_up(1,i,j-1,k,isub) + &
                   &        field_wp(2,i,j-1,k,isub)*field_up(1,i,j  ,k,isub)))&
                   & + &
                   & fac6*((field_wp(3,i  ,j,k,isub)*field_up(1,i,j,k+1,isub) + &
                   &        field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k  ,isub))- &
                   &       (field_wp(3,i,j,k  ,isub)*field_up(1,i,j,k-1,isub) + &
                   &        field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k  ,isub)))

              field_dwp(2,i,j,k,isub) = (&
                   &fac7*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                   &fac8*(field_wp(2,i,j+1,k+1,isub)+field_wp(2,i,j-1,k-1,isub))+&
                   &fac9*(field_wp(2,i,j+1,k-1,isub)+field_wp(2,i,j-1,k+1,isub))-&
                   &2.0_mk*(fac7+fac8+fac9)*field_wp(2,i,j,k,isub)          +&
                   &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                   &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                   &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(2,i,j,k,isub) )*0.5_mk +&
                   & fac4*((field_wp(1,i  ,j,k,isub)*field_up(2,i+1,j,k,isub) + &
                   &        field_wp(1,i+1,j,k,isub)*field_up(2,i  ,j,k,isub))- &
                   &       (field_wp(1,i  ,j,k,isub)*field_up(2,i-1,j,k,isub) + &
                   &        field_wp(1,i-1,j,k,isub)*field_up(2,i  ,j,k,isub)))&
                   & + &
                   & fac5*((field_wp(2,i,j  ,k,isub)*field_up(2,i,j+1,k,isub) + &
                   &        field_wp(2,i,j+1,k,isub)*field_up(2,i,j  ,k,isub))- &
                   &       (field_wp(2,i  ,j,k,isub)*field_up(2,i,j-1,k,isub) + &
                   &        field_wp(2,i,j-1,k,isub)*field_up(2,i,j  ,k,isub)))&
                   & + &
                   & fac6*((field_wp(3,i  ,j,k,isub)*field_up(2,i,j,k+1,isub) + &
                   &        field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k  ,isub))- &
                   &       (field_wp(3,i,j,k  ,isub)*field_up(2,i,j,k-1,isub) + &
                   &        field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k  ,isub)))

              field_dwp(3,i,j,k,isub) = (&
                   &fac7*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                   &fac8*(field_wp(3,i,j+1,k+1,isub)+field_wp(3,i,j-1,k-1,isub))+&
                   &fac9*(field_wp(3,i,j+1,k-1,isub)+field_wp(3,i,j-1,k+1,isub))-&
                   &2.0_mk*(fac7+fac8+fac9)*field_wp(3,i,j,k,isub)          +&
                   &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                   &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                   &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(3,i,j,k,isub) )*0.5_mk +&
                   & fac4*((field_wp(1,i  ,j,k,isub)*field_up(3,i+1,j,k,isub) + &
                   &        field_wp(1,i+1,j,k,isub)*field_up(3,i  ,j,k,isub))- &
                   &       (field_wp(1,i  ,j,k,isub)*field_up(3,i-1,j,k,isub) + &
                   &        field_wp(1,i-1,j,k,isub)*field_up(3,i  ,j,k,isub)))&
                   & + &
                   & fac5*((field_wp(2,i,j  ,k,isub)*field_up(3,i,j+1,k,isub) + &
                   &        field_wp(2,i,j+1,k,isub)*field_up(3,i,j  ,k,isub))- &
                   &       (field_wp(2,i  ,j,k,isub)*field_up(3,i,j-1,k,isub) + &
                   &        field_wp(2,i,j-1,k,isub)*field_up(3,i,j  ,k,isub)))&
                   & + &
                   & fac6*((field_wp(3,i  ,j,k,isub)*field_up(3,i,j,k+1,isub) + &
                   &        field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k  ,isub))- &
                   &       (field_wp(3,i,j,k  ,isub)*field_up(3,i,j,k-1,isub) + &
                   &        field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k  ,isub)))

           END DO

        END DO

     END DO

  END DO

  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_dgammadt',msg,info)
  
END SUBROUTINE wvic_dgammadt_MUSCLE


!------------------------------------------------------------------------------!
!  compute laplacian and stretching term for the vorticity balance
!  use a fourth order laplacian
!------------------------------------------------------------------------------!
SUBROUTINE wvic_dgammadt_4TH (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE


  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out)              :: info

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER :: isub, i , j , k, isubl
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9
  REAL(mk) :: tim1s, tim1e
  CHARACTER(len=256) :: msg

  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()

  fac1 = 1.0_mk/dx**2*nu/12.0_mk
  fac2 = 1.0_mk/dy**2*nu/12.0_mk
  fac3 = 1.0_mk/dz**2*nu/12.0_mk

  fac4 = 8.0_mk/dx/12.0_mk
  fac5 = 8.0_mk/dy/12.0_mk
  fac6 = 8.0_mk/dz/12.0_mk

  fac7 = 1.0_mk/dx/12.0_mk
  fac8 = 1.0_mk/dy/12.0_mk
  fac9 = 1.0_mk/dz/12.0_mk

  
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)

                   
              field_dwp(1,i,j,k,isub) = (&
                   &16.0_mk*fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))-&
                   &        fac1*(field_wp(1,i+2,j,k,isub)+field_wp(1,i-2,j,k,isub))+&
                   &16.0_mk*fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))-&
                   &        fac2*(field_wp(1,i,j+2,k,isub)+field_wp(1,i,j-2,k,isub))+&
                   &16.0_mk*fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                   &        fac3*(field_wp(1,i,j,k+2,isub)+field_wp(1,i,j,k-2,isub))-&
                   &30.0_mk*(fac1+fac2+fac3)*field_wp(1,i,j,k,isub))          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k-1,isub))-&
                   & fac7*(field_wp(1,i+2,j,k,isub)*field_up(1,i+2,j,k,isub)-&
                   &       field_wp(1,i-2,j,k,isub)*field_up(1,i-2,j,k,isub))-&
                   & fac8*(field_wp(2,i,j+2,k,isub)*field_up(1,i,j+2,k,isub)-&
                   &       field_wp(2,i,j-2,k,isub)*field_up(1,i,j-2,k,isub))-&
                   & fac9*(field_wp(3,i,j,k+2,isub)*field_up(1,i,j,k+2,isub)-&
                   &       field_wp(3,i,j,k-2,isub)*field_up(1,i,j,k-2,isub))

              field_dwp(2,i,j,k,isub) = (&
                   &16.0_mk*fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))-&
                   &        fac1*(field_wp(2,i+2,j,k,isub)+field_wp(2,i-2,j,k,isub))+&
                   &16.0_mk*fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))-&
                   &        fac2*(field_wp(2,i,j+2,k,isub)+field_wp(2,i,j-2,k,isub))+&
                   &16.0_mk*fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                   &        fac3*(field_wp(2,i,j,k+2,isub)+field_wp(2,i,j,k-2,isub))-&
                   &30.0_mk*(fac1+fac2+fac3)*field_wp(2,i,j,k,isub))          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k-1,isub))-&
                   & fac7*(field_wp(1,i+2,j,k,isub)*field_up(2,i+2,j,k,isub)-&
                   &       field_wp(1,i-2,j,k,isub)*field_up(2,i-2,j,k,isub))-&
                   & fac8*(field_wp(2,i,j+2,k,isub)*field_up(2,i,j+2,k,isub)-&
                   &       field_wp(2,i,j-2,k,isub)*field_up(2,i,j-2,k,isub))-&
                   & fac9*(field_wp(3,i,j,k+2,isub)*field_up(2,i,j,k+2,isub)-&
                   &       field_wp(3,i,j,k-2,isub)*field_up(2,i,j,k-2,isub))

              field_dwp(3,i,j,k,isub) = (&
                   &16.0_mk*fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))-&
                   &        fac1*(field_wp(3,i+2,j,k,isub)+field_wp(3,i-2,j,k,isub))+&
                   &16.0_mk*fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))-&
                   &        fac2*(field_wp(3,i,j+2,k,isub)+field_wp(3,i,j-2,k,isub))+&
                   &16.0_mk*fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                   &        fac3*(field_wp(3,i,j,k+2,isub)+field_wp(3,i,j,k-2,isub))-&
                   &30.0_mk*(fac1+fac2+fac3)*field_wp(3,i,j,k,isub))          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k-1,isub))-&
                   & fac7*(field_wp(1,i+2,j,k,isub)*field_up(3,i+2,j,k,isub)-&
                   &       field_wp(1,i-2,j,k,isub)*field_up(3,i-2,j,k,isub))-&
                   & fac8*(field_wp(2,i,j+2,k,isub)*field_up(3,i,j+2,k,isub)-&
                   &       field_wp(2,i,j-2,k,isub)*field_up(3,i,j-2,k,isub))-&
                   & fac9*(field_wp(3,i,j,k+2,isub)*field_up(3,i,j,k+2,isub)-&
                   &       field_wp(3,i,j,k-2,isub)*field_up(3,i,j,k-2,isub))

           END DO

        END DO

     END DO

  END DO

  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_dgammadt',msg,info)
  
END SUBROUTINE wvic_dgammadt_4TH



!------------------------------------------------------------------------------!
!  compute laplacian and stretching term for the vorticity balance
!  use the pcon (pseudo-conservative) scheme
!------------------------------------------------------------------------------!
SUBROUTINE wvic_dgammadt_pcons (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE


  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out)              :: info

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER :: isub, i , j , k, isubl
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6, tripme(2)
  REAL(mk) :: tim1s, tim1e
  CHARACTER(len=256) :: msg
  LOGICAL :: stopme = .FALSE.
  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()

  fac1 = 1.0_mk/dx**2*nu
  fac2 = 1.0_mk/dy**2*nu
  fac3 = 1.0_mk/dz**2*nu

  fac4 = 0.5_mk/dx
  fac5 = 0.5_mk/dy
  fac6 = 0.5_mk/dz

  
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)
                   
              field_dwp(1,i,j,k,isub) = &
                   &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                   &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                   &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(1,i,j,k,isub)          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k-1,isub))

              field_dwp(2,i,j,k,isub) = &
                   &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                   &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                   &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(2,i,j,k,isub)          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k-1,isub))
              
              field_dwp(3,i,j,k,isub) = &
                   &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                   &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                   &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(3,i,j,k,isub)          +&
                   & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                   &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                   & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                   &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))+&
                   & fac6*(field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k+1,isub)-&
                   &       field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k-1,isub))
                   
           END DO

        END DO

     END DO

  END DO

  IF(itime.EQ.1.AND.g_istage.EQ.2) THEN
     
     IF(u_infty(1).GT.0.0_mk) THEN
        WRITE(msg,'(A,I3.3,A)') runtag(1:iruntag),itime,'.badie.o'
        OPEN(18,file=msg)
        DO i=1,ndata(1,1)
           WRITE(18,*) i,field_dwp(1:3,i,16,2,1)
        END DO
        CLOSE(18)
     END IF
     
     IF(u_infty(2).GT.0.0_mk) THEN
        WRITE(msg,'(A,I3.3,A)') runtag(1:iruntag),itime,'.badie.o'
        OPEN(18,file=msg)
        DO i=1,ndata(1,1)
           WRITE(18,*) i,field_dwp(1:3,16,i,2,1)
        END DO
        CLOSE(18)
     END IF
     
  END IF

  
  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_dgammadt',msg,info)
  
END SUBROUTINE wvic_dgammadt_pcons

!------------------------------------------------------------------------------!
!  compute laplacian and stretching term for the vorticity balance
!  use the nonconservative scheme
!------------------------------------------------------------------------------!
SUBROUTINE wvic_dgammadt_ncons (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out)              :: info

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER :: isub, i , j , k, isubl
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6
  REAL(mk) :: tim1s, tim1e
  CHARACTER(len=256) :: msg

  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()

  fac1 = 1.0_mk/dx**2*nu
  fac2 = 1.0_mk/dy**2*nu
  fac3 = 1.0_mk/dz**2*nu

  fac4 = 0.5_mk/dx
  fac5 = 0.5_mk/dy
  fac6 = 0.5_mk/dz

  
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)
              field_dwp(1,i,j,k,isub) = &
                   &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                   &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                   &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(1,i,j,k,isub)          +&
#ifndef __WITHOUT_STRETCHING                  
                   &fac4*(field_wp(1,i,j,k,isub)*( &
                   & field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub)))+    &
                   &fac5*(field_wp(2,i,j,k,isub)*( &
                   & field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub)))+    &
                   &fac6*(field_wp(3,i,j,k,isub)*( &
                   & field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub)))
#else
                   &0.0_mk
#endif                   

              field_dwp(2,i,j,k,isub) = &
                   &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                   &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                   &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(2,i,j,k,isub)          +&
#ifndef __WITHOUT_STRETCHING                  
                   &fac4*(field_wp(1,i,j,k,isub)*( &
                   & field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub)))+    &
                   &fac5*(field_wp(2,i,j,k,isub)*( &
                   & field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub)))+    &
                   &fac6*(field_wp(3,i,j,k,isub)*( &
                   & field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub)))
#else
                   &0.0_mk
#endif                   

              field_dwp(3,i,j,k,isub) = &
                   &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                   &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                   &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                   &2.0_mk*(fac1+fac2+fac3)*field_wp(3,i,j,k,isub)          +&
#ifndef __WITHOUT_STRETCHING                  
                   &fac4*(field_wp(1,i,j,k,isub)*( &
                   & field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub)))+    &
                   &fac5*(field_wp(2,i,j,k,isub)*( &
                   & field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub)))+    &
                   &fac6*(field_wp(3,i,j,k,isub)*( &
                   & field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub)))
#else
                   &0.0_mk
#endif                   

           END DO

        END DO

     END DO

  END DO

  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_dgammadt',msg,info)
  
END SUBROUTINE wvic_dgammadt_ncons



!------------------------------------------------------------------------------!
!  Artificial Dissipation Model
!------------------------------------------------------------------------------!
SUBROUTINE wvic_les_blue (info)

  USE module_wvic
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  !  Argumentation
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(inout) :: info
  
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER                :: i,j,k,isub,isubl
  REAL(mk)               :: un,us,ve,vw,wt,wb
  REAL(mk)               :: uc,vc,wc
  REAL(mk), DIMENSION(3) :: on,os,oe,ow,ot,ob,oc
  REAL(mk), DIMENSION(3) :: les
  REAL(mk)               :: dxi

  dxi = 1.0_mk/dx
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)

              !----------------------------------------------------------------!
              ! compute the velocity differeces
              !----------------------------------------------------------------!
              uc = field_up(1,i,j,k,isub)
              vc = field_up(2,i,j,k,isub)
              wc = field_up(3,i,j,k,isub)
              oc = field_wp(1:3,i,j,k,isub)
              
              un = uc - field_up(1,i+1,j,k,isub)
              us = uc - field_up(1,i-1,j,k,isub)

              ve = vc - field_up(2,i,j+1,k,isub)
              vw = vc - field_up(2,i,j-1,k,isub)

              wt = wc - field_up(3,i,j,k+1,isub)
              wb = wc - field_up(3,i,j,k-1,isub)

              !----------------------------------------------------------------!
              ! compute fluxes
              !----------------------------------------------------------------!
              on(1) = oc(1) - field_wp(1,i+1,j,k,isub)
              os(1) = oc(1) - field_wp(1,i-1,j,k,isub)
              oe(1) = oc(1) - field_wp(1,i,j+1,k,isub)
              ow(1) = oc(1) - field_wp(1,i,j-1,k,isub)
              ot(1) = oc(1) - field_wp(1,i,j,k+1,isub)
              ob(1) = oc(1) - field_wp(1,i,j,k-1,isub)

              on(2) = oc(2) - field_wp(2,i+1,j,k,isub)
              os(2) = oc(2) - field_wp(2,i-1,j,k,isub)
              oe(2) = oc(2) - field_wp(2,i,j+1,k,isub)
              ow(2) = oc(2) - field_wp(2,i,j-1,k,isub)
              ot(2) = oc(2) - field_wp(2,i,j,k+1,isub)
              ob(2) = oc(2) - field_wp(2,i,j,k-1,isub)

              on(3) = oc(3) - field_wp(3,i+1,j,k,isub)
              os(3) = oc(3) - field_wp(3,i-1,j,k,isub)
              oe(3) = oc(3) - field_wp(3,i,j+1,k,isub)
              ow(3) = oc(3) - field_wp(3,i,j-1,k,isub)
              ot(3) = oc(3) - field_wp(3,i,j,k+1,isub)
              ob(3) = oc(3) - field_wp(3,i,j,k-1,isub)

              !----------------------------------------------------------------!
              ! compute les term
              !----------------------------------------------------------------!
              les = 0.0_mk
              les(1) = les(1) + MIN( un,0.0_mk) * on(1)
              les(1) = les(1) + MIN(-us,0.0_mk) * os(1)
              les(1) = les(1) + MIN( ve,0.0_mk) * oe(1)
              les(1) = les(1) + MIN(-vw,0.0_mk) * ow(1)
              les(1) = les(1) + MIN( wt,0.0_mk) * ot(1)
              les(1) = les(1) + MIN(-wb,0.0_mk) * ob(1)

              les(2) = les(2) + MIN( un,0.0_mk) * on(2)
              les(2) = les(2) + MIN(-us,0.0_mk) * os(2)
              les(2) = les(2) + MIN( ve,0.0_mk) * oe(2)
              les(2) = les(2) + MIN(-vw,0.0_mk) * ow(2)
              les(2) = les(2) + MIN( wt,0.0_mk) * ot(2)
              les(2) = les(2) + MIN(-wb,0.0_mk) * ob(2)

              les(3) = les(3) + MIN( un,0.0_mk) * on(3)
              les(3) = les(3) + MIN(-us,0.0_mk) * os(3)
              les(3) = les(3) + MIN( ve,0.0_mk) * oe(3)
              les(3) = les(3) + MIN(-vw,0.0_mk) * ow(3)
              les(3) = les(3) + MIN( wt,0.0_mk) * ot(3)
              les(3) = les(3) + MIN(-wb,0.0_mk) * ob(3)

              field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub) &
                   & + 0.5_mk * les(1) * dxi
              field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub) &
                   & + 0.5_mk * les(2) * dxi
              field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub) &
                   & + 0.5_mk * les(3) * dxi
              
           END DO

        END DO

     END DO

  END DO

END SUBROUTINE wvic_les_blue
