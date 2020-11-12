!-------------------------------------------------------------------------------
!* filename: wvic_compvel                                                    *!
!* project : wvic                                                            *!
!* purpose : compuate the gradient of the stream function                     *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 12:49:09 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!-------------------------------------------------------------------------------
! $Log: wvic_compvel.F,v $
! Revision 1.3  2006/09/27 09:30:21  pchatela
! Fixes, spectra calculation,
! most importantly: moved the u_infty out, so it does not kill the dgammadt
!
! Revision 1.2  2006/09/03 14:34:15  pchatela
! Added trivial routine to add U_infty in case of spectral velocities
! Fixed inverse transforms in the fft_velocities_bgw, the inverse transform should depend on the qty considered
!
! Revision 1.1.1.1  2006/07/25 15:13:46  menahel
! initial import
!
! Revision 1.3  2005/11/24 08:45:52  michaebe
! removed minus sign form curl -> we re going to conform to the std
!
! Revision 1.2  2005/11/21 17:36:36  michaebe
! major session. now let''s debug
!
! Revision 1.1  2005/09/28 11:40:28  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
!  Compute the velocites from the stream function
!-------------------------------------------------------------------------------
SUBROUTINE wvic_compvel (info)

  USE module_wvic
  USE ppm_module_write

  !-----------------------------------------------------------------------------
  !  Arguments
  !-----------------------------------------------------------------------------
  INTEGER , INTENT(out) :: info


  !-----------------------------------------------------------------------------
  !  Localities
  !-----------------------------------------------------------------------------
  INTEGER  :: i,j,k,isub,isubl
  INTEGER  :: maptype
  REAL(mk) :: tim1s, tim1e, fac1, fac2, fac3
  REAL(mk) :: phixy, phixz, phiyx, phiyz, phizx, phizy
  REAL(mk), DIMENSION(3) :: umax
  CHARACTER(len=256) :: msg
  REAL(mk), DIMENSION(-1:1) :: filter = (/-0.166666666666666_mk,&
                                         &1.333333333333333_mk,&
                                         &-0.166666666666666_mk/)
  REAL(mk), DIMENSION(3)    :: mom_tot, gmom_tot ! total momentum
  REAL(mk)                  :: wpabs, wpabs_tot
  REAL(mk)                  :: gwpabs_tot
  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()
  
  fac1 = 0.5_mk / dx
  fac2 = 0.5_mk / dy
  fac3 = 0.5_mk / dz
  
  IF(verbose) WRITE(msg,*) 'LBOUND of field_wps:',LBOUND(field_wps)
  IF(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)
  IF(verbose) WRITE(msg,*) 'UBOUND of field_wps:',UBOUND(field_wps)
  IF(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)
  
  umax = -HUGE(umax)
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)
        
        DO j=1,ndata(2,isubl)
           
           DO i=1,ndata(1,isubl)
              
              phixy = fac2*(field_wps(1,i,j+1,k,isub)-field_wps(1,i,j-1,k,isub))
              phixz = fac3*(field_wps(1,i,j,k+1,isub)-field_wps(1,i,j,k-1,isub))
              phiyx = fac1*(field_wps(2,i+1,j,k,isub)-field_wps(2,i-1,j,k,isub))
              phiyz = fac3*(field_wps(2,i,j,k+1,isub)-field_wps(2,i,j,k-1,isub))
              phizx = fac1*(field_wps(3,i+1,j,k,isub)-field_wps(3,i-1,j,k,isub))
              phizy = fac2*(field_wps(3,i,j+1,k,isub)-field_wps(3,i,j-1,k,isub))

              field_up(1,i,j,k,isub) = phizy - phiyz
              field_up(2,i,j,k,isub) = phixz - phizx
              field_up(3,i,j,k,isub) = phiyx - phixy

              IF(ABS(field_up(1,i,j,k,isub)).GT.umax(1)) umax(1) &
                   &=  ABS(field_up(1,i,j,k,isub))
              
              IF(ABS(field_up(2,i,j,k,isub)).GT.umax(2)) umax(2) &
                   &=  ABS(field_up(2,i,j,k,isub))
              
              IF(ABS(field_up(3,i,j,k,isub)).GT.umax(3)) umax(3) &
                   &=  ABS(field_up(3,i,j,k,isub))
              
              
           END DO

        END DO

     END DO

  END DO

  IF(itime.EQ.1.AND.g_istage.EQ.2) THEN
     
     IF(u_infty(1).GT.0.0_mk) THEN
        WRITE(msg,'(A,I3.3,A)') runtag(1:iruntag),itime,'.badie.u'
        OPEN(18,file=msg)
        DO i=1,ndata(1,1)
           WRITE(18,*) i,field_up(1:3,i,16,2,1)
        END DO
        CLOSE(18)
     END IF
     
     IF(u_infty(2).GT.0.0_mk) THEN
        WRITE(msg,'(A,I3.3,A)') runtag(1:iruntag),itime,'.badie.u'
        OPEN(18,file=msg)
        DO i=1,ndata(1,1)
           WRITE(18,*) i,field_up(1:3,16,i,2,1)
        END DO
        CLOSE(18)
     END IF
     
  END IF

  
  IF(MOD(itime,ndump).EQ.0) THEN
     PRINT*,'hello'
     WRITE(msg,'(A,A,I2.2,A,I4.4,A)') &
          &runtag(1:iruntag),'R',rank,'I',itime,'.velprof'
     OPEN(20,file=msg)
     WRITE(20,*) 1,0.0_mk,0.0_mk,field_wp(1,32,32,1,1),field_wp(2,32,32,1,1)
     DO k=2,ndata(3,1)-1
        WRITE(20,*) k,field_up(1,32,32,k,1),field_up(2,32,32,k,1),field_wp(1,32,32,k,1),field_wp(2,32&
             &,32,k,1)
     END DO
     WRITE(20,*) ndata(3,1),0.0_mk,0.0_mk,field_wp(1,32,32,ndata(3,1),1),field_wp(2,32&
          &,32,ndata(3,1),1)
     CLOSE(20)
  END IF

  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'umax ',umax
  if(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)  
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)

END SUBROUTINE wvic_compvel




!-------------------------------------------------------------------------------
!  Compute the velocites from the stream function
!  with a fourth order finite difference stencil
!-------------------------------------------------------------------------------
SUBROUTINE wvic_compvel_4th (info)

  USE module_wvic
  USE ppm_module_write

  !-----------------------------------------------------------------------------
  !  Arguments
  !-----------------------------------------------------------------------------
  INTEGER , INTENT(out) :: info


  !-----------------------------------------------------------------------------
  !  Localities
  !-----------------------------------------------------------------------------
  INTEGER  :: i,j,k,isub,isubl
  INTEGER  :: maptype
  REAL(mk) :: tim1s, tim1e, fac1, fac2, fac3, fac4, fac5, fac6
  REAL(mk) :: phixy, phixz, phiyx, phiyz, phizx, phizy
  REAL(mk), DIMENSION(3) :: umax
  CHARACTER(len=256) :: msg
  REAL(mk), DIMENSION(-1:1) :: filter = (/-0.166666666666666_mk,&
                                         &1.333333333333333_mk,&
                                         &-0.166666666666666_mk/)
  REAL(mk), DIMENSION(3)    :: mom_tot, gmom_tot ! total momentum
  REAL(mk)                  :: wpabs, wpabs_tot
  REAL(mk)                  :: gwpabs_tot
  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()
  
  fac1 = 8.0_mk / dx / 12.0_mk
  fac2 = 8.0_mk / dy / 12.0_mk
  fac3 = 8.0_mk / dz / 12.0_mk
  fac4 = 1.0_mk / dx / 12.0_mk
  fac5 = 1.0_mk / dy / 12.0_mk
  fac6 = 1.0_mk / dz / 12.0_mk
  
  
  IF(verbose) WRITE(msg,*) 'LBOUND of field_wps:',LBOUND(field_wps)
  IF(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)
  IF(verbose) WRITE(msg,*) 'UBOUND of field_wps:',UBOUND(field_wps)
  IF(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)
  
  umax = -HUGE(umax)
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)
        
        DO j=1,ndata(2,isubl)
           
           DO i=1,ndata(1,isubl)
              
              phixy=-fac2*(field_wps(1,i,j+1,k,isub)-field_wps(1,i,j-1,k,isub))&
                   &+fac5*(field_wps(1,i,j+2,k,isub)-field_wps(1,i,j-2,k,isub))
              phixz=-fac3*(field_wps(1,i,j,k+1,isub)-field_wps(1,i,j,k-1,isub))&
                   &+fac6*(field_wps(1,i,j,k+2,isub)-field_wps(1,i,j,k-2,isub))
              phiyx=-fac1*(field_wps(2,i+1,j,k,isub)-field_wps(2,i-1,j,k,isub))&
                   &+fac4*(field_wps(2,i+2,j,k,isub)-field_wps(2,i-2,j,k,isub))
              phiyz=-fac3*(field_wps(2,i,j,k+1,isub)-field_wps(2,i,j,k-1,isub))&
                   &+fac6*(field_wps(2,i,j,k+2,isub)-field_wps(2,i,j,k-2,isub))
              phizx=-fac1*(field_wps(3,i+1,j,k,isub)-field_wps(3,i-1,j,k,isub))&
                   &+fac4*(field_wps(3,i+2,j,k,isub)-field_wps(3,i-2,j,k,isub))
              phizy=-fac2*(field_wps(3,i,j+1,k,isub)-field_wps(3,i,j-1,k,isub))&
                   &+fac5*(field_wps(3,i,j+2,k,isub)-field_wps(3,i,j-2,k,isub))

              field_up(1,i,j,k,isub) = phizy - phiyz
              field_up(2,i,j,k,isub) = phixz - phizx
              field_up(3,i,j,k,isub) = phiyx - phixy

              IF(ABS(field_up(1,i,j,k,isub)).GT.umax(1)) umax(1) &
                   &=  ABS(field_up(1,i,j,k,isub))
              
              IF(ABS(field_up(2,i,j,k,isub)).GT.umax(2)) umax(2) &
                   &=  ABS(field_up(2,i,j,k,isub))
              
              IF(ABS(field_up(3,i,j,k,isub)).GT.umax(3)) umax(3) &
                   &=  ABS(field_up(3,i,j,k,isub))
              
              
           END DO

        END DO

     END DO

  END DO

  
  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'umax ',umax
  if(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)  
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_compvel',msg,info)

END SUBROUTINE wvic_compvel_4th


!-------------------------------------------------------------------------------
!  Compute the velocites spectrally
!  = already done, just need to compute diags
!-------------------------------------------------------------------------------
SUBROUTINE wvic_compvel_spec (info)

  USE module_wvic
  USE ppm_module_write

  !-----------------------------------------------------------------------------
  !  Arguments
  !-----------------------------------------------------------------------------
  INTEGER , INTENT(out) :: info


  !-----------------------------------------------------------------------------
  !  Localities
  !-----------------------------------------------------------------------------
  INTEGER  :: i,j,k,isub,isubl
  INTEGER  :: maptype
  REAL(mk) :: tim1s, tim1e
  REAL(mk), DIMENSION(3) :: umax
  CHARACTER(len=256) :: msg
  REAL(mk), DIMENSION(3)    :: mom_tot, gmom_tot ! total momentum
  REAL(mk)                  :: wpabs, wpabs_tot
  REAL(mk)                  :: gwpabs_tot
  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()
  
  umax = -HUGE(umax)
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)
        
        DO j=1,ndata(2,isubl)
           
           DO i=1,ndata(1,isubl)
              
              IF(ABS(field_up(1,i,j,k,isub)).GT.umax(1)) umax(1) &
                   &=  ABS(field_up(1,i,j,k,isub))
              
              IF(ABS(field_up(2,i,j,k,isub)).GT.umax(2)) umax(2) &
                   &=  ABS(field_up(2,i,j,k,isub))
              
              IF(ABS(field_up(3,i,j,k,isub)).GT.umax(3)) umax(3) &
                   &=  ABS(field_up(3,i,j,k,isub))
              
              
           END DO

        END DO

     END DO

  END DO

  
  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'umax ',umax
  if(verbose) CALL ppm_write(rank,'wvic_compvel_spec',msg,info)  
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_compvel_spec',msg,info)

END SUBROUTINE wvic_compvel_spec


!-------------------------------------------------------------------------------
!  Add the free-stream and compute diags
!-------------------------------------------------------------------------------
SUBROUTINE wvic_compvel_adduinfty (info)

  USE module_wvic
  USE ppm_module_write

  !-----------------------------------------------------------------------------
  !  Arguments
  !-----------------------------------------------------------------------------
  INTEGER , INTENT(out) :: info


  !-----------------------------------------------------------------------------
  !  Localities
  !-----------------------------------------------------------------------------
  INTEGER  :: i,j,k,isub,isubl
  INTEGER  :: maptype
  REAL(mk) :: tim1s, tim1e
  REAL(mk), DIMENSION(3) :: umax
  CHARACTER(len=256) :: msg
  REAL(mk), DIMENSION(3)    :: mom_tot, gmom_tot ! total momentum
  REAL(mk)                  :: wpabs, wpabs_tot
  REAL(mk)                  :: gwpabs_tot
  REAL(mk)                  :: ramp
  INCLUDE 'mpif.h'

  tim1s = MPI_WTIME()
  
  umax = -HUGE(umax)

  !-----------------------------------------------------------------------------
  !  Ramping up U_inf in case of onset flow
  !-----------------------------------------------------------------------------
  IF (((flow_case .EQ. 10) .or. (flow_case .EQ. 13) .or. (flow_case .EQ. 14)) &
     & .and. (time .lt. onset_ramptime) .and. (onset_ramptime .NE. 0.0_mk)) THEN
    ramp=time/onset_ramptime
    IF (ramp .GT. 1.0_mk) THEN
      ramp=1.0_mk
    ENDIF
  ELSE
    ramp=1.0_mk
  ENDIF

  DO isub=1,nsublist
    isubl = isublist(isub)
    
    DO k=1-ghostsize(3),ndata(3,isubl)+ghostsize(3)
        
      DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
        
        DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
              
          field_up(1,i,j,k,isub)=field_up(1,i,j,k,isub) + u_infty(1)*ramp
          field_up(2,i,j,k,isub)=field_up(2,i,j,k,isub) + u_infty(2)*ramp
          field_up(3,i,j,k,isub)=field_up(3,i,j,k,isub) + u_infty(3)*ramp

          IF(ABS(field_up(1,i,j,k,isub)).GT.umax(1)) umax(1) &
               &=  ABS(field_up(1,i,j,k,isub))
          
          IF(ABS(field_up(2,i,j,k,isub)).GT.umax(2)) umax(2) &
               &=  ABS(field_up(2,i,j,k,isub))
          
          IF(ABS(field_up(3,i,j,k,isub)).GT.umax(3)) umax(3) &
               &=  ABS(field_up(3,i,j,k,isub))
              
        END DO

      END DO

    END DO

  END DO

  
  tim1e = MPI_WTIME()
  if(verbose) WRITE(msg,*) 'umax ',umax
  if(verbose) CALL ppm_write(rank,'wvic_compvel_spec',msg,info)  
  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
  if(verbose) CALL ppm_write(rank,'wvic_compvel_spec',msg,info)

END SUBROUTINE wvic_compvel_adduinfty
