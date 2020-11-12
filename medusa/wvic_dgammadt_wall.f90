!------------------------------------------------------------------------------!
!  filename: wvic_dgammadt_wall                                                !
!  project : WVIC                                                              !
!  purpose : compute dgammadt on the wall (one-sided)                          !
!          :                                                                   !
!  author  : Michael Bergdorf                                                  !
!          : Computational Science and Engineering Lab (CSE-Lab)               !
!          : ICOS, ETH Zurich                                                  !
!          :                                                                   !
!  date    : Wed Oct  5 09:21:43 2005                                          !
!  please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]      !
!------------------------------------------------------------------------------!
!  $Log: wvic_dgammadt_wall.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:46  menahel
!  initial import
!
!  Revision 1.6  2005/12/05 08:37:07  michaebe
!  added higher order Laplacian for one space direction
!
!  Revision 1.5  2005/11/21 17:48:47  michaebe
!  was turned off
!
!  Revision 1.4  2005/11/21 17:36:37  michaebe
!  major session. now let''s debug
!
!  Revision 1.3  2005/11/11 19:59:05  michaebe
!  added extrapolation of dgammadt to the ghost layer
!
!  Revision 1.2  2005/11/11 14:04:24  michaebe
!  clean up, additions, comments
!
!  Revision 1.1  2005/10/09 13:22:01  michaebe
!  initial implementation
!
!------------------------------------------------------------------------------!


#define __LOW_ORDER__
  SUBROUTINE wvic_dgammadt_wall

    USE module_wvic
    USE ppm_module_write

    INTEGER :: isub, i, j, k, isubl, iface, info
    REAL(mk) :: fac11(0:4), fac22(0:4), facx, facy, facz
    REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6
    REAL(mk) :: fac1x, fac1y, fac1z, fac2x, fac2y, fac2z
    REAL(mk) :: tim1s, tim1e, dxi, dyi, dzi
    CHARACTER(len=256) :: msg

    INCLUDE 'mpif.h'

    tim1s = MPI_WTIME()

    IF(wvic_dgamma_scheme.NE.4) THEN ! NOT fourth order

       ! gradient
       fac11(0) = -1.5_mk
       fac11(1) = 2.0_mk
       fac11(2) = -0.5_mk
       fac11(3) = 0.0_mk
       fac11(4) = 0.0_mk
       
       fac1x = 1.0_mk/dx
       fac1y = 1.0_mk/dy
       fac1z = 1.0_mk/dz

#ifdef __LOW_ORDER__
       ! laplacian
       fac22(0) = 1.0_mk
       fac22(1) = -2.0_mk
       fac22(2) = 1.0_mk
       fac22(3) = 0.0_mk
       fac22(4) = 0.0_mk
#else
       ! laplacian
       fac22(0) = 35.0_mk/12.0_mk
       fac22(1) = -26.0_mk/3.0_mk
       fac22(2) = 19.0_mk/2.0_mk
       fac22(3) = -14.0_mk/3.0_mk
       fac22(4) = 11.0_mk/12.0_mk
#endif       
       fac2x = nu/dx**2
       fac2y = nu/dy**2
       fac2z = nu/dz**2

       fac1 = 1.0_mk/dx**2*nu
       fac2 = 1.0_mk/dy**2*nu
       fac3 = 1.0_mk/dz**2*nu
       
       fac4 = 0.5_mk/dx
       fac5 = 0.5_mk/dy
       fac6 = 0.5_mk/dz

       
    ELSE

       WRITE(msg,*) '4-th order boundary handling not yet implemented'
       CALL ppm_write(rank,'wvic_dgamma_wall',msg,info)
       STOP

       ! gradient
       fac11(0) = -25.0_mk/12.0_mk
       fac11(1) = 4.0_mk
       fac11(2) = -3.0_mk
       fac11(3) = 4.0_mk/3.0_mk
       fac11(4) = -0.25_mk

       fac1x = 1.0_mk/dx
       fac1y = 1.0_mk/dy
       fac1z = 1.0_mk/dz
       
       ! laplacian
       fac22(0) = 35.0_mk/12.0_mk
       fac22(1) = -26.0_mk/3.0_mk
       fac22(2) = 19.0_mk/2.0_mk
       fac22(3) = -14.0_mk/3.0_mk
       fac22(4) = 11.0_mk/12.0_mk

       fac2x = nu/dx**2
       fac2y = nu/dy**2
       fac2z = nu/dz**2

    END IF

    !--------------------------------------------------------------------------!
    !  Go and correct the faces
    !--------------------------------------------------------------------------!
    DO isub = 1,nsublist

       DO iface = 1,6
          
          IF(.NOT. liswall( iface, isub ) ) CYCLE

          isubl = isublist(isub)
          
          SELECT CASE ( iface )
             
          CASE(1) ! lower YZ plane
             
             i=1
             DO k=1,ndata(3,isubl)
                DO j=1,ndata(2,isubl)

                   field_dwp(1,i,j,k,isub) = &
                        & fac2x*(fac22(0)*field_wp(1,i  ,j,k,isub)   + &
                        &        fac22(1)*field_wp(1,i+1,j,k,isub)   + &
                        &        fac22(2)*field_wp(1,i+2,j,k,isub))  + &
                        &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                        &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                        &2.0_mk*(     fac2+fac3)*field_wp(1,i,j,k,isub)          +&
                        & fac1x*(fac11(0)*field_wp(1,i  ,j,k,isub)   &
                        &               *field_up(1,i  ,j,k,isub) + &
                        &        fac11(1)*field_wp(1,i+1,j,k,isub)   &
                        &               *field_up(1,i+1,j,k,isub) + &
                        &        fac11(2)*field_wp(1,i+2,j,k,isub)   &
                        &               *field_up(1,i+2,j,k,isub))+ &
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k-1,isub))
                   
                   field_dwp(2,i,j,k,isub) = &
                        & fac2x*(fac22(0)*field_wp(2,i  ,j,k,isub)   + &
                        &        fac22(1)*field_wp(2,i+1,j,k,isub)   + &
                        &        fac22(2)*field_wp(2,i+2,j,k,isub))  + &
                        &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                        &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                        &2.0_mk*(     fac2+fac3)*field_wp(2,i,j,k,isub)          +&
                        & fac1x*(fac11(0)*field_wp(1,i  ,j,k,isub)   &
                        &               *field_up(2,i  ,j,k,isub) + &
                        &        fac11(1)*field_wp(1,i+1,j,k,isub)   &
                        &               *field_up(2,i+1,j,k,isub) + &
                        &        fac11(2)*field_wp(1,i+2,j,k,isub)   &
                        &               *field_up(2,i+2,j,k,isub))+ &
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k-1,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2x*(fac22(0)*field_wp(3,i  ,j,k,isub)   + &
                        &        fac22(1)*field_wp(3,i+1,j,k,isub)   + &
                        &        fac22(2)*field_wp(3,i+2,j,k,isub))  + &
                        &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                        &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                        &2.0_mk*(     fac2+fac3)*field_wp(3,i,j,k,isub)          +&
                        & fac1x*(fac11(0)*field_wp(1,i  ,j,k,isub)   &
                        &               *field_up(3,i  ,j,k,isub) + &
                        &        fac11(1)*field_wp(1,i+1,j,k,isub)   &
                        &               *field_up(3,i+1,j,k,isub) + &
                        &        fac11(2)*field_wp(1,i+2,j,k,isub)   &
                        &               *field_up(3,i+2,j,k,isub))+ &
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k-1,isub))

                END DO
             END DO
             
          CASE(2) ! upper YZ plane
             
             i=ndata(1,isubl)
             DO k=1,ndata(3,isubl)
                DO j=1,ndata(2,isubl)

                   field_dwp(1,i,j,k,isub) = &
                        & fac2x*(fac22(0)*field_wp(1,i  ,j,k,isub)   + &
                        &        fac22(1)*field_wp(1,i-1,j,k,isub)   + &
                        &        fac22(2)*field_wp(1,i-2,j,k,isub))  + &
                        &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))+&
                        &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                        &2.0_mk*(     fac2+fac3)*field_wp(1,i,j,k,isub)          -&
                        & fac1x*(fac11(0)*field_wp(1,i  ,j,k,isub)   &
                        &               *field_up(1,i  ,j,k,isub) + &
                        &        fac11(1)*field_wp(1,i-1,j,k,isub)   &
                        &               *field_up(1,i-1,j,k,isub) + &
                        &        fac11(2)*field_wp(1,i-2,j,k,isub)   &
                        &               *field_up(1,i-2,j,k,isub))+ &
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k-1,isub))
                   
                   field_dwp(2,i,j,k,isub) = &
                        & fac2x*(fac22(0)*field_wp(2,i  ,j,k,isub)   + &
                        &        fac22(1)*field_wp(2,i-1,j,k,isub)   + &
                        &        fac22(2)*field_wp(2,i-2,j,k,isub))  + &
                        &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))+&
                        &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                        &2.0_mk*(     fac2+fac3)*field_wp(2,i,j,k,isub)          -&
                        & fac1x*(fac11(0)*field_wp(1,i  ,j,k,isub)   &
                        &               *field_up(2,i  ,j,k,isub) + &
                        &        fac11(1)*field_wp(1,i-1,j,k,isub)   &
                        &               *field_up(2,i-1,j,k,isub) + &
                        &        fac11(2)*field_wp(1,i-2,j,k,isub)   &
                        &               *field_up(2,i-2,j,k,isub))+ &
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k-1,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2x*(fac22(0)*field_wp(3,i  ,j,k,isub)   + &
                        &        fac22(1)*field_wp(3,i-1,j,k,isub)   + &
                        &        fac22(2)*field_wp(3,i-2,j,k,isub))  + &
                        &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))+&
                        &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                        &2.0_mk*(     fac2+fac3)*field_wp(3,i,j,k,isub)          -&
                        & fac1x*(fac11(0)*field_wp(1,i  ,j,k,isub)   &
                        &               *field_up(3,i  ,j,k,isub) + &
                        &        fac11(1)*field_wp(1,i-1,j,k,isub)   &
                        &               *field_up(3,i-1,j,k,isub) + &
                        &        fac11(2)*field_wp(1,i-2,j,k,isub)   &
                        &               *field_up(3,i-2,j,k,isub))+ &
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k-1,isub))

                END DO
             END DO

          CASE(3) ! lower xz plane
             j = 1
             DO k=1,ndata(3,isubl)
                DO i=1,ndata(1,isubl)
                   
                   field_dwp(1,i,j,k,isub) = &
                        & fac2y*(fac22(0)*field_wp(1,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(1,i,j+1,k,isub)   + &
                        &        fac22(2)*field_wp(1,i,j+2,k,isub))  + &
                        &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                        &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                        &2.0_mk*(fac1     +fac3)*field_wp(1,i,j,k,isub)          +&
                        & fac1y*(fac11(0)*field_wp(2,i,j  ,k,isub)   &
                        &               *field_up(1,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(2,i,j+1,k,isub)   &
                        &               *field_up(1,i,j+1,k,isub) + &
                        &        fac11(2)*field_wp(2,i,j+2,k,isub)   &
                        &               *field_up(1,i,j+2,k,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k-1,isub))

                   field_dwp(2,i,j,k,isub) = &
                        & fac2y*(fac22(0)*field_wp(2,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(2,i,j+1,k,isub)   + &
                        &        fac22(2)*field_wp(2,i,j+2,k,isub))  + &
                        &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                        &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                        &2.0_mk*(fac1     +fac3)*field_wp(2,i,j,k,isub)          +&
                        & fac1y*(fac11(0)*field_wp(2,i,j  ,k,isub)   &
                        &               *field_up(2,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(2,i,j+1,k,isub)   &
                        &               *field_up(2,i,j+1,k,isub) + &
                        &        fac11(2)*field_wp(2,i,j+2,k,isub)   &
                        &               *field_up(2,i,j+2,k,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k-1,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2y*(fac22(0)*field_wp(3,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(3,i,j+1,k,isub)   + &
                        &        fac22(2)*field_wp(3,i,j+2,k,isub))  + &
                        &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                        &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                        &2.0_mk*(fac1     +fac3)*field_wp(3,i,j,k,isub)          +&
                        & fac1y*(fac11(0)*field_wp(2,i,j  ,k,isub)   &
                        &               *field_up(3,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(2,i,j+1,k,isub)   &
                        &               *field_up(3,i,j+1,k,isub) + &
                        &        fac11(2)*field_wp(2,i,j+2,k,isub)   &
                        &               *field_up(3,i,j+2,k,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k-1,isub))

                   
                   
                END DO

             END DO


          CASE(4) ! lower xz plane
             j = ndata(2,isubl)
             DO k=1,ndata(3,isubl)
                DO i=1,ndata(1,isubl)
                   
                   field_dwp(1,i,j,k,isub) = &
                        & fac2y*(fac22(0)*field_wp(1,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(1,i,j-1,k,isub)   + &
                        &        fac22(2)*field_wp(1,i,j-2,k,isub))  + &
                        &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                        &fac3*(field_wp(1,i,j,k+1,isub)+field_wp(1,i,j,k-1,isub))-&
                        &2.0_mk*(fac1     +fac3)*field_wp(1,i,j,k,isub)          -&
                        & fac1y*(fac11(0)*field_wp(2,i,j  ,k,isub)   &
                        &               *field_up(1,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(2,i,j-1,k,isub)   &
                        &               *field_up(1,i,j-1,k,isub) + &
                        &        fac11(2)*field_wp(2,i,j-2,k,isub)   &
                        &               *field_up(1,i,j-2,k,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(1,i,j,k-1,isub))

                   field_dwp(2,i,j,k,isub) = &
                        & fac2y*(fac22(0)*field_wp(2,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(2,i,j-1,k,isub)   + &
                        &        fac22(2)*field_wp(2,i,j-2,k,isub))  + &
                        &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                        &fac3*(field_wp(2,i,j,k+1,isub)+field_wp(2,i,j,k-1,isub))-&
                        &2.0_mk*(fac1     +fac3)*field_wp(2,i,j,k,isub)          -&
                        & fac1y*(fac11(0)*field_wp(2,i,j  ,k,isub)   &
                        &               *field_up(2,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(2,i,j-1,k,isub)   &
                        &               *field_up(2,i,j-1,k,isub) + &
                        &        fac11(2)*field_wp(2,i,j-2,k,isub)   &
                        &               *field_up(2,i,j-2,k,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(2,i,j,k-1,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2y*(fac22(0)*field_wp(3,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(3,i,j-1,k,isub)   + &
                        &        fac22(2)*field_wp(3,i,j-2,k,isub))  + &
                        &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                        &fac3*(field_wp(3,i,j,k+1,isub)+field_wp(3,i,j,k-1,isub))-&
                        &2.0_mk*(fac1     +fac3)*field_wp(3,i,j,k,isub)          -&
                        & fac1y*(fac11(0)*field_wp(2,i,j  ,k,isub)   &
                        &               *field_up(3,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(2,i,j-1,k,isub)   &
                        &               *field_up(3,i,j-1,k,isub) + &
                        &        fac11(2)*field_wp(2,i,j-2,k,isub)   &
                        &               *field_up(3,i,j-2,k,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                        & fac6*(field_wp(3,i,j,k+1,isub)*field_up(3,i,j,k+1,isub)-&
                        &       field_wp(3,i,j,k-1,isub)*field_up(3,i,j,k-1,isub))

                   
                   
                END DO

             END DO

          CASE(5)

             k=1
             DO j=1,ndata(2,isubl)
                DO i=1,ndata(1,isubl)
#ifdef __LOW_ORDER__                   
                   field_dwp(1,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(1,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(1,i,j,k+1,isub)   + &
                        &        fac22(2)*field_wp(1,i,j,k+2,isub))  + &
                        &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                        &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(1,i,j,k,isub)          +&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(1,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k+1,isub)   &
                        &               *field_up(1,i,j,k+1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k+2,isub)   &
                        &               *field_up(1,i,j,k+2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))

                   field_dwp(2,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(2,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(2,i,j,k+1,isub)   + &
                        &        fac22(2)*field_wp(2,i,j,k+2,isub))  + &
                        &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                        &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(2,i,j,k,isub)          +&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(2,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k+1,isub)   &
                        &               *field_up(2,i,j,k+1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k+2,isub)   &
                        &               *field_up(2,i,j,k+2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(3,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(3,i,j,k+1,isub)   + &
                        &        fac22(2)*field_wp(3,i,j,k+2,isub))  + &
                        &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                        &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(3,i,j,k,isub)          +&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(3,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k+1,isub)   &
                        &               *field_up(3,i,j,k+1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k+2,isub)   &
                        &               *field_up(3,i,j,k+2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))
#else
                   field_dwp(1,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(1,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(1,i,j,k+1,isub)   + &
                        &        fac22(2)*field_wp(1,i,j,k+2,isub)  + &
                        &        fac22(3)*field_wp(1,i,j,k+3,isub)   + &
                        &        fac22(4)*field_wp(1,i,j,k+4,isub))  + &
                        &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                        &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(1,i,j,k,isub)          +&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(1,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k+1,isub)   &
                        &               *field_up(1,i,j,k+1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k+2,isub)   &
                        &               *field_up(1,i,j,k+2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))

                   field_dwp(2,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(2,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(2,i,j,k+1,isub)   + &
                        &        fac22(2)*field_wp(2,i,j,k+2,isub)   + &
                        &        fac22(1)*field_wp(2,i,j,k+3,isub)   + &
                        &        fac22(2)*field_wp(2,i,j,k+4,isub))  + &
                        &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                        &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(2,i,j,k,isub)          +&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(2,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k+1,isub)   &
                        &               *field_up(2,i,j,k+1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k+2,isub)   &
                        &               *field_up(2,i,j,k+2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(3,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(3,i,j,k+1,isub)   + &
                        &        fac22(2)*field_wp(3,i,j,k+2,isub)  + &
                        &        fac22(3)*field_wp(3,i,j,k+3,isub)   + &
                        &        fac22(4)*field_wp(3,i,j,k+4,isub))  + &
                        &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                        &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(3,i,j,k,isub)          +&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(3,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k+1,isub)   &
                        &               *field_up(3,i,j,k+1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k+2,isub)   &
                        &               *field_up(3,i,j,k+2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))
                   
#endif
                END DO
             END DO

          CASE(6)

             k=ndata(3,isubl)
             DO j=1,ndata(2,isubl)
                DO i=1,ndata(1,isubl)

#ifdef __LOW_ORDER__
                   field_dwp(1,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(1,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(1,i,j,k-1,isub)   + &
                        &        fac22(2)*field_wp(1,i,j,k-2,isub))  + &
                        &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                        &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(1,i,j,k,isub)      -&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(1,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k-1,isub)   &
                        &               *field_up(1,i,j,k-1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k-2,isub)   &
                        &               *field_up(1,i,j,k-2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))

                   field_dwp(2,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(2,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(2,i,j,k-1,isub)   + &
                        &        fac22(2)*field_wp(2,i,j,k-2,isub))  + &
                        &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                        &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(2,i,j,k,isub)      -&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(2,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k-1,isub)   &
                        &               *field_up(2,i,j,k-1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k-2,isub)   &
                        &               *field_up(2,i,j,k-2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(3,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(3,i,j,k-1,isub)   + &
                        &        fac22(2)*field_wp(3,i,j,k-2,isub))  + &
                        &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                        &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(3,i,j,k,isub)      -&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(3,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k-1,isub)   &
                        &               *field_up(3,i,j,k-1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k-2,isub)   &
                        &               *field_up(3,i,j,k-2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))
#else
                   field_dwp(1,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(1,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(1,i,j,k-1,isub)   + &
                        &        fac22(2)*field_wp(1,i,j,k-2,isub)   + &
                        &        fac22(3)*field_wp(1,i,j,k-3,isub)   + &
                        &        fac22(4)*field_wp(1,i,j,k-4,isub))  + &
                        &fac1*(field_wp(1,i+1,j,k,isub)+field_wp(1,i-1,j,k,isub))+&
                        &fac2*(field_wp(1,i,j+1,k,isub)+field_wp(1,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(1,i,j,k,isub)      -&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(1,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k-1,isub)   &
                        &               *field_up(1,i,j,k-1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k-2,isub)   &
                        &               *field_up(1,i,j,k-2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(1,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(1,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(1,i,j-1,k,isub))

                   field_dwp(2,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(2,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(2,i,j,k-1,isub)   + &
                        &        fac22(2)*field_wp(2,i,j,k-2,isub)   + &
                        &        fac22(3)*field_wp(2,i,j,k-3,isub)   + &
                        &        fac22(4)*field_wp(2,i,j,k-4,isub))  + &
                        &fac1*(field_wp(2,i+1,j,k,isub)+field_wp(2,i-1,j,k,isub))+&
                        &fac2*(field_wp(2,i,j+1,k,isub)+field_wp(2,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(2,i,j,k,isub)      -&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(2,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k-1,isub)   &
                        &               *field_up(2,i,j,k-1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k-2,isub)   &
                        &               *field_up(2,i,j,k-2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(2,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(2,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(2,i,j-1,k,isub))

                   field_dwp(3,i,j,k,isub) = &
                        & fac2z*(fac22(0)*field_wp(3,i,j  ,k,isub)   + &
                        &        fac22(1)*field_wp(3,i,j,k-1,isub)   + &
                        &        fac22(2)*field_wp(3,i,j,k-2,isub)   + &
                        &        fac22(3)*field_wp(3,i,j,k-3,isub)   + &
                        &        fac22(4)*field_wp(3,i,j,k-4,isub))  + &
                        &fac1*(field_wp(3,i+1,j,k,isub)+field_wp(3,i-1,j,k,isub))+&
                        &fac2*(field_wp(3,i,j+1,k,isub)+field_wp(3,i,j-1,k,isub))-&
                        &2.0_mk*(fac1+fac2         )*field_wp(3,i,j,k,isub)      -&
                        & fac1z*(fac11(0)*field_wp(3,i,j  ,k,isub)   &
                        &               *field_up(3,i,j  ,k,isub) + &
                        &        fac11(1)*field_wp(3,i,j,k-1,isub)   &
                        &               *field_up(3,i,j,k-1,isub) + &
                        &        fac11(2)*field_wp(3,i,j,k-2,isub)   &
                        &               *field_up(3,i,j,k-2,isub))+ &
                        & fac4*(field_wp(1,i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                        &       field_wp(1,i-1,j,k,isub)*field_up(3,i-1,j,k,isub))+&
                        & fac5*(field_wp(2,i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                        &       field_wp(2,i,j-1,k,isub)*field_up(3,i,j-1,k,isub))

#endif                   
                END DO

             END DO
                
          END SELECT

       END DO

    END DO

  END SUBROUTINE wvic_dgammadt_wall
  

  !----------------------------------------------------------------------------!
  !  Extrapolate dgammadt to the ghost layer(s)
  !----------------------------------------------------------------------------!
  SUBROUTINE wvic_dgammadt_wall_ghosts

    USE module_wvic
    USE ppm_module_write
    
    INTEGER :: isub, i, j, k, isubl, iface, info
    REAL(mk) :: fac11(0:4), fac22(0:4), facx, facy, facz
    REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6
    REAL(mk) :: fac1x, fac1y, fac1z, fac2x, fac2y, fac2z
    REAL(mk) :: tim1s, tim1e, dxi, dyi, dzi
    CHARACTER(len=256) :: msg

    !--------------------------------------------------------------------------!
    !  Rationale
    !  loop over the faces, check whether they are wall, and if they are do
    !  either linear or quadratic interpolation to get the value at x=-h
    !--------------------------------------------------------------------------!

    DO isub = 1,nsublist

       DO iface = 1,6
          
          IF(.NOT. liswall( iface, isub ) ) CYCLE

          isubl = isublist(isub)
          
          SELECT CASE ( iface )
             
          CASE(1) ! lower YZ plane
             
             i=0
             DO k=0,ndata(3,isubl)+1
                DO j=0,ndata(2,isubl)+1
                   field_dwp(1,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(1,i+1,j,k,isub) - &
                        &        field_dwp(1,i+2,j,k,isub)
                   field_dwp(2,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(2,i+1,j,k,isub) - &
                        &        field_dwp(2,i+2,j,k,isub)
                   field_dwp(3,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(3,i+1,j,k,isub) - &
                        &        field_dwp(3,i+2,j,k,isub)
                END DO
             END DO
             
          CASE(2) ! upper YZ plane
             
             i=ndata(1,isubl)+1
             DO k=0,ndata(3,isubl)+1
                DO j=0,ndata(2,isubl)+1
                   field_dwp(1,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(1,i-1,j,k,isub) - &
                        &        field_dwp(1,i-2,j,k,isub)
                   field_dwp(2,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(2,i-1,j,k,isub) - &
                        &        field_dwp(2,i-2,j,k,isub)
                   field_dwp(3,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(3,i-1,j,k,isub) - &
                        &        field_dwp(3,i-2,j,k,isub)
                END DO
             END DO
             
          CASE(3) ! lower xz plane
             j=0
             DO k=0,ndata(3,isubl)+1
                DO i=0,ndata(1,isubl)+1
                   field_dwp(1,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(1,i,j+1,k,isub) - &
                        &        field_dwp(1,i,j+2,k,isub)
                   field_dwp(2,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(2,i,j+1,k,isub) - &
                        &        field_dwp(2,i,j+2,k,isub)
                   field_dwp(3,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(3,i,j+1,k,isub) - &
                        &        field_dwp(3,i,j+2,k,isub)
                END DO
             END DO
             
             
          CASE(4) ! lower xz plane
             j=ndata(2,isubl)+1
             DO k=0,ndata(3,isubl)+1
                DO i=0,ndata(1,isubl)+1
                   field_dwp(1,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(1,i,j-1,k,isub) - &
                        &        field_dwp(1,i,j-2,k,isub)
                   field_dwp(2,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(2,i,j-1,k,isub) - &
                        &        field_dwp(2,i,j-2,k,isub)
                   field_dwp(3,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(3,i,j-1,k,isub) - &
                        &        field_dwp(3,i,j-2,k,isub)
                END DO
             END DO
             
          CASE(5)
             
             k=0
             DO j=0,ndata(2,isubl)+1
                DO i=0,ndata(1,isubl)+1
                   field_dwp(1,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(1,i,j,k+1,isub) - &
                        &        field_dwp(1,i,j,k+2,isub)
                   field_dwp(2,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(2,i,j,k+1,isub) - &
                        &        field_dwp(2,i,j,k+2,isub)
                   field_dwp(3,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(3,i,j,k+1,isub) - &
                        &        field_dwp(3,i,j,k+2,isub)
                END DO
             END DO
             
          CASE(6)
             
             k=ndata(3,isubl)+1
             DO j=0,ndata(2,isubl)+1
                DO i=0,ndata(1,isubl)+1
                   field_dwp(1,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(1,i,j,k-1,isub) - &
                        &        field_dwp(1,i,j,k-2,isub)
                   field_dwp(2,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(2,i,j,k-1,isub) - &
                        &        field_dwp(2,i,j,k-2,isub)
                   field_dwp(3,i,j,k,isub) = &
                        & 2.0_mk*field_dwp(3,i,j,k-1,isub) - &
                        &        field_dwp(3,i,j,k-2,isub)
                END DO
             END DO
                
          END SELECT

       END DO

    END DO

  END SUBROUTINE wvic_dgammadt_wall_ghosts
  
