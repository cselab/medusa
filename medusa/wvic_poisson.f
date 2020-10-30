!------------------------------------------------------------------------------!
!* filename: wvic_poisson                                                     *!
!* project : ppm                                                              *!
!* purpose : solve poisson eq                                                 *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 15:10:29 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!------------------------------------------------------------------------------!
! $Log: wvic_poisson.F,v $
! Revision 1.3  2006/09/11 14:57:27  pchatela
! Fixed velocity computation with odd/even symmetries
! Added flag for adaptive time step
!
! Revision 1.2  2006/08/29 16:44:13  menahel
! made it call the symmetric ffts in case of trailvortex
!
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.6  2005/12/10 20:15:59  michaebe
! clean up
!
! Revision 1.5  2005/11/24 08:47:36  michaebe
! made the fishpack poisson equation conform to the std
!
! Revision 1.4  2005/11/21 17:36:39  michaebe
! major session. now let''s debug
!
! Revision 1.3  2005/11/11 14:04:25  michaebe
! clean up, additions, comments
!
! Revision 1.2  2005/11/10 10:25:01  michaebe
! added a fishpack based poisson solver
!
! Revision 1.1  2005/09/28 11:40:41  michaebe
! Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!  Solve the poisson equation for the streamfunction
!------------------------------------------------------------------------------!
SUBROUTINE wvic_poisson_fft (info)

  
  USE module_wvic
  USE ppm_module_write
  USE ppm_module_fdsolver_solve
  USE ppm_module_fft
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER , INTENT(out) :: info
  
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER, DIMENSION(3) :: ftopo_id
  INTEGER, DIMENSION(4) :: fmesh_id, t_topoid
  REAL(mk)              :: tim1s, tim1e 
  REAL(mk)              :: pim1s, pim1e
  INTEGER               :: i,j,k,isub,isubl
  REAL(mk)              :: scaler
  CHARACTER(len=256)    :: msg
  REAL(mk), DIMENSION(3)    :: mom_tot, gmom_tot ! total momentum
  REAL(mk)                  :: wpabs, wpabs_tot
  REAL(mk)                  :: gwpabs_tot
  INCLUDE 'mpif.h'
  
  tim1s = MPI_WTIME()
  CALL cpu_time(pim1s)
  
  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)
  !----------------------------------------------------------------------------!
  ! copy the stuff
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_up(1,i,j,k,isub) = field_wp(1,i,j,k,isub) 
              field_up(2,i,j,k,isub) = field_wp(2,i,j,k,isub) 
              field_up(3,i,j,k,isub) = field_wp(3,i,j,k,isub)
           END DO
        END DO
     END DO
  END DO
  IF(wvic_compvel_scheme.EQ.2) THEN
     IF(trailvortex) THEN
        CALL ppm_fft_velocities_bgw(field_up,mesh_id,topo_id,t_topoid,fmesh_id, & 
             & ghostsize, info)
     ELSE
        CALL ppm_fft_velocities(field_up,mesh_id,topo_id,t_topoid,fmesh_id, & 
             & ghostsize, info)
     END IF
  ELSE
     CALL ppm_fft_potential(field_up,mesh_id,topo_id,t_topoid,fmesh_id, & 
          & ghostsize, info)
     GOTO 2003
  END IF
  scaler = 1.0_mk/REAL((nx(1)-1)*(nx(2)-1)*(nx(3)-1),mk)
  IF(wvic_compvel_scheme.EQ.2) THEN
     
     IF(g_istage.EQ.1) THEN
        
        wpabs_tot = 0.0_mk
        mom_tot = 0.0_mk
        
        DO isub=1,nsublist
           isubl = isublist(isub)
           
           DO k=1,ndata(3,isubl)-1
              
              DO j=1,ndata(2,isubl)-1
                 
                 DO i=1,ndata(1,isubl)-1
                    wpabs      = SQRT(field_wp(1,i,j,k,isub)**2 +&
                         &            field_wp(2,i,j,k,isub)**2 +&
                         &            field_wp(3,i,j,k,isub)**2) * dx*dy*dz
                    wpabs_tot = wpabs_tot + wpabs
                    mom_tot(1) = mom_tot(1) + field_up(1,i,j,k,isub) * wpabs
                    mom_tot(2) = mom_tot(2) + field_up(2,i,j,k,isub) * wpabs
                    mom_tot(3) = mom_tot(3) + field_up(3,i,j,k,isub) * wpabs
                 END DO
                 
              END DO
              
           END DO
           
        END DO
        
        !----------------------------------------------------------------------!
        ! all reduce
        !----------------------------------------------------------------------!
        gmom_tot   = 0.0_mk
        gwpabs_tot = 0.0_mk
        CALL MPI_Allreduce(mom_tot,gmom_tot,3,mpi_prec,MPI_SUM,comm,info)
        CALL MPI_Allreduce(wpabs_tot,gwpabs_tot,1,mpi_prec,MPI_SUM,comm,info)
        
        mom_tot = gmom_tot / gwpabs_tot
        u_cmass = mom_tot
        !----------------------------------------------------------------------!
        ! killit
        !----------------------------------------------------------------------!
        
     ELSE
        
        mom_tot = u_cmass
        
     END IF
     
     
  ELSE

     DO isub=1,nsublist
        isubl = isublist(isub)
        DO k=1,ndata(3,isubl)
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)
                 field_wps(1,i,j,k,isub) = field_wps(1,i,j,k,isub)*scaler 
                 field_wps(2,i,j,k,isub) = field_wps(2,i,j,k,isub)*scaler 
                 field_wps(3,i,j,k,isub) = field_wps(3,i,j,k,isub)*scaler 
              END DO
           END DO
        END DO
     END DO
     
  END IF
2003 CONTINUE

  tim1e = MPI_WTIME()
  CALL cpu_time(pim1e)
  
  IF(verbose) WRITE(msg,*) 'took ',INT(1000.0*(tim1e-tim1s)),&
       &' [cpu:',INT(1000.0*(pim1e-pim1s)),']'
  IF(verbose) CALL ppm_write(rank,'wvic_poisson',msg,info)
  
END SUBROUTINE wvic_poisson_fft




!------------------------------------------------------------------------------!
!  Solve the poisson equation for the streamfunction
!------------------------------------------------------------------------------!
SUBROUTINE wvic_poisson_mg (info)

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_mg
  
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE wvic_alloc_field (vfield_up, info)
       USE module_wvic
       IMPLICIT NONE
       REAL (mk), DIMENSION (:, :, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field
  END INTERFACE
  
  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER , INTENT(out) :: info
  
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER               :: iteration, niterations, iterloc
  REAL(mk)              :: tim1s, tim1e 
  REAL(mk)              :: pim1s, pim1e
  LOGICAL,SAVE          :: firsttime = .TRUE.
  CHARACTER(len=256)    :: msg
  REAL(mk)              :: c1, c2, c3
  REAL(mk)              :: cresid1, cresid2, cresid3
  INTEGER               :: i,j,k,isub,isubl
  REAL(mk), DIMENSION(3):: l0rhs, gl0rhs

  
  INCLUDE 'mpif.h'
  
  c1 = 0.0_mk
  c2 = 0.0_mk
  c3 = 0.0_mk
  cresid1 = 0.0_mk;cresid2 = 0.0_mk;cresid3 = 0.0_mk
  IF(verbose) WRITE(msg,*) 'entered - a safe bet'
  IF(verbose) CALL ppm_write(rank,'wvic_poisson_mg',msg,info)
  tim1s = MPI_WTIME()

#ifndef __TIMING  
  IF(firsttime) THEN
     niterations = 5
     field_wps = 0.0_mk
     firsttime = .FALSE.
  ELSE
     niterations = 1
  END IF

!  CALL wvic_alloc_field(field_rhs,3,info)
  IF(info.NE.0) CALL ppm_write(rank,'wvic_poisson',&
       &'allocation faild, we''re dead',info)
  
  !----------------------------------------------------------------------------!
  ! split-copy the stuff
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!
  ! compute the circulation and subtract it
  !----------------------------------------------------------------------------!
   l0rhs = 0.0_mk
   DO isub=1,nsublist
      isubl = isublist(isub)
      DO k=1,ndata(3,isubl)-1
         DO j=1,ndata(2,isubl)-1
            DO i=1,ndata(1,isubl)-1
               l0rhs(1) = l0rhs(1) + field_wp(1,i,j,k,isub)
               l0rhs(2) = l0rhs(2) + field_wp(2,i,j,k,isub)
               l0rhs(3) = l0rhs(3) + field_wp(3,i,j,k,isub)
            END DO
         END DO
      END DO
   END DO

   CALL MPI_Allreduce(l0rhs,gl0rhs,3,mpi_prec,MPI_SUM,comm,info)
   l0rhs = gl0rhs / REAL((nx(1)-1)*(nx(2)-1)*(nx(3)-1),mk)

   DO isub=1,nsublist
      isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_rhs(1,i,j,k,isub) = -field_wp(1,i,j,k,isub) - l0rhs(1)
              field_rhs(2,i,j,k,isub) = -field_wp(2,i,j,k,isub) - l0rhs(2)
              field_rhs(3,i,j,k,isub) = -field_wp(3,i,j,k,isub) - l0rhs(3)
           END DO
        END DO
     END DO
  END DO

  DO iteration=1,niterations
     !-------------------------------------------------------------------------!
     ! solve
     !-------------------------------------------------------------------------!
     CALL ppm_mg_solv(field_wps,field_rhs,3,2,2,10,4,cresid1,info)
     iterloc = 1
     IF(niterations.LE.1) THEN
        cresid1 = cresid1 / max_vorticity 
        c1 = cresid1
        !----------------------------------------------------------------------!
        ! If target residual not reached, do another iteration
        !----------------------------------------------------------------------!
        IF(cresid1.GT.1.25e-3_mk) THEN
           CALL ppm_mg_solv(field_wps,field_rhs,3,2,2,10,4,cresid2,info)
           cresid2 = cresid2 / max_vorticity
           c2 = cresid2
           iterloc = 2
        END IF
        !----------------------------------------------------------------------!
        ! And once more
        !----------------------------------------------------------------------!
        IF(cresid2.GT.1.25e-3_mk) THEN
           CALL ppm_mg_solv(field_wps,field_rhs,3,2,2,10,4,cresid3,info)
           cresid3 = cresid3 / max_vorticity
           c3 = cresid3
           iterloc = 3
        END IF
        !----------------------------------------------------------------------!
        ! Complain if it didnt work out
        !----------------------------------------------------------------------!
        IF(cresid3.GT.1.25e-3_mk) THEN
           CALL ppm_write(rank,'wvic_poisson_mg',&
                &'Failed to reach target residual',info)
        END IF
        
        IF(.TRUE.) THEN
           WRITE(msg,*) 'scaled residual = ',c1,' ',c2,' ',c3,' using # ',iterloc
           CALL ppm_write(rank,'wvic_poisson',msg,info)
        END IF
     END IF
  END DO
  
  
  !----------------------------------------------------------------------------!
  !\/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/
#else
  IF(firsttime) THEN
     niterations = 2
     field_angel = 0.0_mk
     firsttime = .FALSE.
  ELSE
     niterations = 2
  END IF
  field_anger(1:3,:,:,:,:) = field_wp(1:3,:,:,:,:)
  DO iteration=1,niterations
     CALL ppm_mg_solv(field_angel,field_anger,3,1,1,10,4,cresid,info)
  END DO
#endif
  !/\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\
  !----------------------------------------------------------------------------!
  
  tim1e = MPI_WTIME()
  ! DEALLOCATE(field_rhs,stat=info)
  IF(info.NE.0) CALL ppm_write(rank,'wvic_poisson',&
       &' dealloc failed. We''re dead',info)
 !NULLIFY(field_rhs)
  
  IF(verbose) WRITE(msg,*) 'took ',INT(1000.0*(tim1e-tim1s))
  if(verbose) CALL ppm_write(rank,'wvic_poisson',msg,info)
END SUBROUTINE wvic_poisson_mg

