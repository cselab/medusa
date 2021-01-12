!-------------------------------------------------------------------------------
! WVIC_MOVEMENT
! 2007/2008
! Routines for moving objects, i.e. the step function
!
! wvic_solid_velocity updates the solid velocity field
!
! Johannes Tophoej Rasmussen betonarbejder@gmail.com
!-------------------------------------------------------------------------------

!JTR rename to wvic_solid_velocity_uniform
SUBROUTINE wvic_solid_velocity(info)

  USE module_wvic
  USE ppm_module_map_field_ghost
  USE ppm_module_data

  INTEGER, INTENT(out)  :: info
  INTEGER               :: i,j,k
  INTEGER               :: isub, isubl
  INTEGER               :: maptype
  REAL(mk)              :: pi

  pi = ACOS(-1.0_mk)

  IF ((SUM(ABS(harmonic_amplitude)) .NE. 0.0_mk) .AND. &
  &(SUM(ABS(harmonic_period)) .NE. 0.0_mk)) THEN
    IF (harmonic_amplitude(1) .NE. 0.0_mk) THEN
      u_solid(1) = harmonic_amplitude(1)*2*pi/harmonic_period(1) * &
                 & COS(2*pi/harmonic_period(1)*time+2*pi*harmonic_phase(1))
    END IF
    IF (harmonic_amplitude(2) .NE. 0.0_mk) THEN
      u_solid(2) = harmonic_amplitude(2)*2*pi/harmonic_period(2) * &
                 & COS(2*pi/harmonic_period(2)*time+2*pi*harmonic_phase(2))
    END IF
    IF (harmonic_amplitude(3) .NE. 0.0_mk) THEN
      u_solid(3) = harmonic_amplitude(3)*2*pi/harmonic_period(3) * &
                 & COS(2*pi/harmonic_period(3)*time+2*pi*harmonic_phase(3))
    END IF
  ENDIF


  !---------------------------------------------------------------------------
  ! Calculating solid velocity - so far uniform
  !---------------------------------------------------------------------------
  DO isub=1,nsublist
    isubl = isublist(isub)
      DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
          DO i=1,ndata(1,isubl)
            field_ubar(:,i,j,k,isub) = u_solid
          END DO !i
        END DO !j
      END DO !k
  END DO

  !---------------------------------------------------------------------------
  ! Ghosting - remove and update ghosts manually
  !---------------------------------------------------------------------------
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)

END SUBROUTINE wvic_solid_velocity


!JTR add description
! ?This an attempt to move the stepfunction by convection? usuccesful due to high gradients yes? That is abandoned?
SUBROUTINE wvic_advect_stepfunction(info)

  USE module_wvic
  USE ppm_module_map_field_ghost
  USE ppm_module_data



  INTEGER, INTENT(out)  :: info
  INTEGER               :: i,j,k
  INTEGER               :: isub, isubl
  INTEGER               :: maptype
  REAL(mk)              :: fac1,fac2,fac3

  fac1=1.0_mk/(dx*12.0_mk)
  fac2=1.0_mk/(dy*12.0_mk)
  fac3=1.0_mk/(dz*12.0_mk)

  !---------------------------------------------------------------------------
  !JTR TODO: make 3 point stencil
  !---------------------------------------------------------------------------
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO k=1,ndata(3,isubl)-1
      DO j=1,ndata(2,isubl)-1
        DO i=1,ndata(1,isubl)-1
          field_dwp(1,i,j,k,isub) = field_H(i,j,k,isub) - dt*( &
             & fac1 * &
             & (   -field_H(i+2,j,k,isub)*field_ubar(1,i+2,j,k,isub) + & 
             &  8.0*field_H(i+1,j,k,isub)*field_ubar(1,i+1,j,k,isub) &
             & -8.0*field_H(i-1,j,k,isub)*field_ubar(1,i-1,j,k,isub) + &
             &      field_H(i-2,j,k,isub)*field_ubar(1,i-2,j,k,isub)) + &
             & fac2 * &
             & (   -field_H(i,j+2,k,isub)*field_ubar(2,i,j+2,k,isub) + & 
             &  8.0*field_H(i,j+1,k,isub)*field_ubar(2,i,j+1,k,isub) &
             & -8.0*field_H(i,j-1,k,isub)*field_ubar(2,i,j-1,k,isub) + &
             &      field_H(i,j-2,k,isub)*field_ubar(2,i,j-2,k,isub)) + &
             & fac3 * &
             & (   -field_H(i,j,k+2,isub)*field_ubar(3,i,j,k+2,isub) + & 
             &  8.0*field_H(i,j,k+1,isub)*field_ubar(3,i,j,k+1,isub) &
             & -8.0*field_H(i,j,k-1,isub)*field_ubar(3,i,j,k-1,isub) + &
             &      field_H(i,j,k-2,isub)*field_ubar(3,i,j,k-2,isub)) )
        END DO !i
      END DO !j
    END DO !k
    DO k=1,ndata(3,isubl)-1
      DO j=1,ndata(2,isubl)-1
        DO i=1,ndata(1,isubl)-1
          field_H(i,j,k,isub) = field_dwp(1,i,j,k,isub)
        END DO !i
      END DO !j
    END DO !k
  END DO

!  ftopo_id = (/2,3,4/)
!  t_topoid = (/2,3,4,5/)
!  fmesh_id = (/2,3,4,5/)

  !---------------------------------------------------------------------------
  ! Ghosting
  !---------------------------------------------------------------------------
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)

END SUBROUTINE wvic_advect_stepfunction



SUBROUTINE wvic_advect_stepparticles(info)

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_interp_p2m
  USE ppm_module_map_part
  USE ppm_module_topo_check
  USE ppm_module_map_field_ghost



  INTEGER, INTENT(inout):: info
  INTEGER               :: i,j,k
  INTEGER               :: isub, isubl
  INTEGER               :: maptype
  INTEGER               :: mpart
  CHARACTER(len=256)       :: msg
  REAL(mk)              :: tmp,tmp2

  !-----------------------------------------------------
  !  Advect particles - implement rotation later
  !-----------------------------------------------------
  DO i=1,npc
    xpc(1,i) = xpc(1,i) + dt*u_solid(1)
    xpc(2,i) = xpc(2,i) + dt*u_solid(2)
    xpc(3,i) = xpc(3,i) + dt*u_solid(3)
  END DO

  !-----------------------------------------------------
  !  Catch particles leaving the domain
  !-----------------------------------------------------
  CALL wvic_pbc_stepfunc(info)

  !-----------------------------------------------------
  !  Map particles
  !-----------------------------------------------------
  maptype = ppm_param_map_partial
  CALL ppm_map_part(xpc,dime,npc,mpart,topo_id,maptype,info)
  maptype = ppm_param_map_push 
  CALL ppm_map_part(hpc,npc,mpart,topo_id,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_part(hpc,npc,mpart,topo_id,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_part(hpc,npc,mpart,topo_id,maptype,info)
  CALL ppm_map_part(xpc,dime,npc,mpart,topo_id,maptype,info)
  npc=mpart

  !-----------------------------------------------------
  !  Interpolate particle step function to mesh
  !-----------------------------------------------------
  !field_H = 0.0_mk
  CALL ppm_interp_p2m(xpc,npc,hpc,topo_id,mesh_id, &
       & ppm_param_rmsh_kernel_mp4,ghostsize,&
       & field_H,info)

!tmp=0.0_mk
!DO i=1,npc
!  tmp=tmp+xpc(3,i)
!END DO
!tmp2=0.0_mk
!  DO isub=1,nsublist
!    isubl = isublist(isub)
!    DO k=1,ndata(3,isubl)-1
!      DO j=1,ndata(2,isubl)-1
!        DO i=1,ndata(1,isubl)-1
!          tmp2=tmp2+(SUM(field_up(:,i,j,k,isub)))**2!!field_H(i,j,k,isub)
!        END DO !i
!      END DO !j
!    END DO !k
!  END DO
!
!WRITE(msg,*) 'JTR: rank,itime,npc', rank, itime, tmp2,'\n'
!WRITE(0,*) msg

  !-----------------------------------------------------
  !  Ghost stepfunction
  !-----------------------------------------------------
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
       & ghostsize,maptype,info)

END SUBROUTINE wvic_advect_stepparticles


#ifdef _JOHSSLUK


#endif


SUBROUTINE wvic_pbc_stepfunc (info)
  
  USE module_wvic
  USE ppm_module_data
  USE ppm_module_error
  USE ppm_module_write
  USE ppm_module_impose_part_bc

  !----------------------------------------------------------------------------!
  !  arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out) :: info
  
  !----------------------------------------------------------------------------!
  !  locals
  !----------------------------------------------------------------------------!
  INTEGER              :: kp,dispair
  CHARACTER(len=256)   :: msg,pnm
  REAL(mk), DIMENSION(3) :: len_physg
  
  WRITE(pnm,*) 'wvic_pbci_stepfunc'
  
  len_physg = max_physg - min_physg
  
!  DO dispair=1,2  !JTR why this? omitting...
     !----------------------------------------------------------------------!
     !  x direction
     !----------------------------------------------------------------------!
     DO kp=1,npc
        IF(xpc(1,kp).GE.max_physg(1)) THEN
           xpc(1,kp) = xpc(1,kp) - len_physg(1)
        ELSE
           IF(xpc(1,kp).LT.min_physg(1)) THEN
              xpc(1,kp) = xpc(1,kp) + len_physg(1)
           END IF
        END IF
     END DO

     !----------------------------------------------------------------------!
     !  y direction
     !----------------------------------------------------------------------!
     DO kp=1,npc
        IF(xpc(2,kp).GE.max_physg(2)) THEN
           xpc(2,kp) = xpc(2,kp) - len_physg(2)
        ELSE
           IF(xpc(2,kp).LT.min_physg(2)) THEN
              xpc(2,kp) = xpc(2,kp) + len_physg(2)
           END IF
        END IF
     END DO
     !----------------------------------------------------------------------!
     !  z direction
     !----------------------------------------------------------------------!
     DO kp=1,npc
        IF(xpc(3,kp).GE.max_physg(3)) THEN
           xpc(3,kp) = xpc(3,kp) - len_physg(3)
        ELSE
           IF(xpc(3,kp).LT.min_physg(3)) THEN
               xpc(3,kp) = xpc(3,kp) + len_physg(3)
           END IF
        END IF
     END DO
!  END DO

9999 CONTINUE  
  
END SUBROUTINE wvic_pbc_stepfunc


