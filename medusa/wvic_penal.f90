!-------------------------------------------------------------------------------
! WVIC_PENAL
! 2007/2008
! Impose penalization on dwp field
!
! Johannes Tophoej Rasmussen betonarbejder@gmail.com
!-------------------------------------------------------------------------------
! JRT - 9/5-2008: Adding implicit penalization
!
!-------------------------------------------------------------------------------

SUBROUTINE wvic_penalization_explicit (info)

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
  INTEGER                    :: isub, i , j , k, isubl
!  REAL(mk)                   :: lambda
  CHARACTER(len=256)         :: msg,filename
!  REAL(mk)                   :: tx,ty,tz
  REAL(mk)                   :: fac1,fac2,fac3,fac4,fac5,fac6
  REAL(mk), DIMENSION(3)     :: penalw, penalu, penal
  REAL(mk)                   :: maxpw,maxpu,maxrhs,magnitude
  REAL(mk), DIMENSION(3)     :: alphau,alphaw,alphaRHS,galphaRHS,galphausum,galphawsum
  REAL(mk)                   :: dV
  INTEGER                    :: ios
  
  INCLUDE 'mpif.h'




!  lambda=penalization_lambda/dt

#ifdef _showsmaximumdomegacontributionfromstretchinganddiffusion
  maxpw=0.0_mk
  maxpu=0.0_mk
  maxrhs=0.0_mk
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          magnitude=sqrt(field_dwp(1,i,j,k,isub)**2+field_dwp(2,i,j,k,isub)**2+field_dwp(3,i,j,k,isub)**2)
          IF (magnitude .GT. maxrhs) THEN
            maxrhs=magnitude
          END IF
        END DO
      END DO
    END DO
  END DO
#endif

  alphau=0.0_mk
  alphaw=0.0_mk
  alphaRHS=0.0_mk
  dV=dx*dy*dz

  IF(SUM(ghostsize).EQ.3) THEN
    fac1 = lambda * 0.5_mk / dx
    fac2 = lambda * 0.5_mk / dy
    fac3 = lambda * 0.5_mk / dz
  ELSE
     fac1 = lambda * 8.0_mk/dx/12.0_mk
     fac2 = lambda * 8.0_mk/dy/12.0_mk
     fac3 = lambda * 8.0_mk/dz/12.0_mk
     fac4 = -lambda * 1.0_mk/dx/12.0_mk
     fac5 = -lambda * 1.0_mk/dy/12.0_mk
     fac6 = -lambda * 1.0_mk/dz/12.0_mk
  END IF

  DO isub=1,nsublist
    isubl = isublist(isub)

    IF (penalization_clearrhs .EQV. .true.) THEN
      DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
          DO i=1,ndata(1,isubl)
             !------------------------------------------------------------------
             ! recording the contribution from RHS inside solid - as w moment
             ! F = - d(omega)/dt
             !------------------------------------------------------------------
#ifdef _forcefrompenalizationtermsucks
             alphaRHS(1)=alphaRHS(1) + 0.5_mk*dV*field_H(i,j,k,isub)* &
                  & (ty*(field_dwp(3,i,j,k,isub))-tz*(field_dwp(2,i,j,k,isub)))
             alphaRHS(2)=alphaRHS(2) + 0.5_mk*dV*field_H(i,j,k,isub)* &
                  & (tz*(field_dwp(1,i,j,k,isub))-tx*(field_dwp(3,i,j,k,isub)))
             alphaRHS(3)=alphaRHS(3) + 0.5_mk*dV*field_H(i,j,k,isub)* &
                  & (tx*(field_dwp(2,i,j,k,isub))-ty*(field_dwp(1,i,j,k,isub)))
#endif
             !------------------------------------------------------------------
             ! clear dwp contribution from RHS inside solid
             !------------------------------------------------------------------
             field_dwp(1,i,j,k,isub)=field_dwp(1,i,j,k,isub) &
                  & *(1.0_mk-field_H(i,j,k,isub))
             field_dwp(2,i,j,k,isub)=field_dwp(2,i,j,k,isub) &
                  & *(1.0_mk-field_H(i,j,k,isub))
             field_dwp(3,i,j,k,isub)=field_dwp(3,i,j,k,isub) &
                  & *(1.0_mk-field_H(i,j,k,isub))
          END DO
        END DO
      END DO
    END IF
 
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
!          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
!          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
!          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

          IF (penalization_oneterm .EQV. .TRUE.) THEN
            IF(SUM(ghostsize).EQ.3) THEN
              penal(1)=fac2*(-field_H(i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                      & (-field_H(i,j-1,k,isub)*field_up(3,i,j-1,k,isub)))-&
                      & fac3*(-field_H(i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                      & (-field_H(i,j,k-1,isub)*field_up(2,i,j,k-1,isub)))
              penal(2)=fac3*(-field_H(i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                      & (-field_H(i,j,k-1,isub)*field_up(1,i,j,k-1,isub)))-&
                      & fac1*(-field_H(i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                      & (-field_H(i-1,j,k,isub)*field_up(3,i-1,j,k,isub)))
              penal(3)=fac1*(-field_H(i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                      & (-field_H(i-1,j,k,isub)*field_up(2,i-1,j,k,isub)))-&
                      & fac2*(-field_H(i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                      & (-field_H(i,j-1,k,isub)*field_up(1,i,j-1,k,isub)))
          ELSE
              penal(1)=fac2*(-field_H(i,j+1,k,isub)*field_up(3,i,j+1,k,isub)-&
                      & (-field_H(i,j-1,k,isub)*field_up(3,i,j-1,k,isub)))-&
                      & fac3*(-field_H(i,j,k+1,isub)*field_up(2,i,j,k+1,isub)-&
                      & (-field_H(i,j,k-1,isub)*field_up(2,i,j,k-1,isub)))+&
                      & fac5*(-field_H(i,j+2,k,isub)*field_up(3,i,j+2,k,isub)-&
                      & (-field_H(i,j-2,k,isub)*field_up(3,i,j-2,k,isub)))-&
                      & fac6*(-field_H(i,j,k+2,isub)*field_up(2,i,j,k+2,isub)-&
                      & (-field_H(i,j,k-2,isub)*field_up(2,i,j,k-2,isub)))
              penal(2)=fac3*(-field_H(i,j,k+1,isub)*field_up(1,i,j,k+1,isub)-&
                      & (-field_H(i,j,k-1,isub)*field_up(1,i,j,k-1,isub)))-&
                      & fac1*(-field_H(i+1,j,k,isub)*field_up(3,i+1,j,k,isub)-&
                      & (-field_H(i-1,j,k,isub)*field_up(3,i-1,j,k,isub)))+&
                      & fac6*(-field_H(i,j,k+2,isub)*field_up(1,i,j,k+2,isub)-&
                      & (-field_H(i,j,k-2,isub)*field_up(1,i,j,k-2,isub)))-&
                      & fac4*(-field_H(i+2,j,k,isub)*field_up(3,i+2,j,k,isub)-&
                      & (-field_H(i-2,j,k,isub)*field_up(3,i-2,j,k,isub)))
              penal(3)=fac1*(-field_H(i+1,j,k,isub)*field_up(2,i+1,j,k,isub)-&
                      & (-field_H(i-1,j,k,isub)*field_up(2,i-1,j,k,isub)))-&
                      & fac2*(-field_H(i,j+1,k,isub)*field_up(1,i,j+1,k,isub)-&
                      & (-field_H(i,j-1,k,isub)*field_up(1,i,j-1,k,isub)))+&
                      & fac4*(-field_H(i+2,j,k,isub)*field_up(2,i+2,j,k,isub)-&
                      & (-field_H(i-2,j,k,isub)*field_up(2,i-2,j,k,isub)))-&
                      & fac5*(-field_H(i,j+2,k,isub)*field_up(1,i,j+2,k,isub)-&
                      & (-field_H(i,j-2,k,isub)*field_up(1,i,j-2,k,isub)))
            END IF
            field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub)+penal(1)
            field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub)+penal(2)
            field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub)+penal(3)
          ELSE
            penalw(1)=lambda*field_H(i,j,k,isub) * (- field_wp(1,i,j,k,isub))
            penalu(1)=lambda*(field_ubar(2,i,j,k,isub)*(-field_up(3,i,j,k,isub))&
                      & - field_ubar(3,i,j,k,isub)*(-field_up(2,i,j,k,isub)))
            field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub)+penalw(1)+penalu(1)
  
            penalw(2)=lambda*field_H(i,j,k,isub) * (- field_wp(2,i,j,k,isub)) 
            penalu(2)=lambda*(field_ubar(3,i,j,k,isub)*(-field_up(1,i,j,k,isub))&
                      & - field_ubar(1,i,j,k,isub)*(-field_up(3,i,j,k,isub)))
            field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub)+penalw(2)+penalu(2)
  
            penalw(3)=lambda*field_H(i,j,k,isub) * (- field_wp(3,i,j,k,isub))
            penalu(3)=lambda*(field_ubar(1,i,j,k,isub)*(-field_up(2,i,j,k,isub))&
                      & - field_ubar(2,i,j,k,isub)*(-field_up(1,i,j,k,isub)))
            field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub)+penalw(3)+penalu(3)
          END IF


#ifdef _showsmaximumdomegacontributionfromuandomegapenalization
           magnitude=sqrt(penalu(1)**2+penalu(2)**2+penalu(3)**2)
           IF (magnitude .GT. maxpu) THEN
             maxpu=magnitude
           END IF
           magnitude=sqrt(penalw(1)**2+penalw(2)**2+penalw(3)**2)
           IF (magnitude .GT. maxpw) THEN
             maxpw=magnitude
           END IF
#endif
#ifdef _forcefrompenalizationtermsucks
           !alpha integration: int X x w    (vorticity moments)
           alphau(1)=alphau(1) - 0.5_mk*dV*(ty*(penalu(3))-tz*(penalu(2)))
           alphau(2)=alphau(2) - 0.5_mk*dV*(tz*(penalu(1))-tx*(penalu(3)))
           alphau(3)=alphau(3) - 0.5_mk*dV*(tx*(penalu(2))-ty*(penalu(1)))
           alphaw(1)=alphaw(1) - 0.5_mk*dV*(ty*(penalw(3))-tz*(penalw(2)))
           alphaw(2)=alphaw(2) - 0.5_mk*dV*(tz*(penalw(1))-tx*(penalw(3)))
           alphaw(3)=alphaw(3) - 0.5_mk*dV*(tx*(penalw(2))-ty*(penalw(1)))
#endif

        END DO
      END DO
    END DO
  END DO

#ifdef _forcefrompenalizationtermsucks
  !Sum alpha from all nodes and write from node 0
  CALL MPI_Allreduce( alphau, galphausum, 3, mpi_prec,MPI_SUM,comm,info)
  CALL MPI_Allreduce( alphaw, galphawsum, 3, mpi_prec,MPI_SUM,comm,info)
  CALL MPI_Allreduce( alphaRHS, galphaRHS, 3, mpi_prec,MPI_SUM,comm,info)
  IF(rank.EQ.0) THEN
    alphausum=galphausum/(0.5_mk*(u_infty(1)**2+u_infty(2)**2+u_infty(3)**2))
    alphawsum=galphawsum/(0.5_mk*(u_infty(1)**2+u_infty(2)**2+u_infty(3)**2))
    alphaRHSsum=galphaRHS/(0.5_mk*(u_infty(1)**2+u_infty(2)**2+u_infty(3)**2))
  END IF
#endif

  
END SUBROUTINE wvic_penalization_explicit





SUBROUTINE wvic_penalization_implicit (info)

  USE module_wvic
  USE ppm_module_write
  USE ppm_module_interp_p2m
  USE ppm_module_data
  USE ppm_module_map_field_ghost
  IMPLICIT NONE

  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out)              :: info

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER                    :: isub, i , j , k, isubl
  INTEGER                    :: maptype
!  REAL(mk)                   :: lambda
  CHARACTER(len=256)         :: msg
  REAL(mk)                   :: fac1,fac2,fac3,fac4,fac5,fac6
  REAL(mk)                   :: Auy,Auz,Avx,Avz,Awx,Awy,Buy,Buz,Bvx,Bvz,Bwx,Bwy 
  REAL(mk)                   :: AAuy,AAuz,AAvx,AAvz,AAwx,AAwy
  REAL(mk)                   :: BBuy,BBuz,BBvx,BBvz,BBwx,BBwy 
  REAL(mk)                   :: tmp2
!JTMP
  REAL(mk)                 :: time1, time2
  REAL(mk)                 :: tx_center
  REAL(mk)                 :: ty_center
  REAL(mk)                 :: tz_center
  REAL(mk)                 :: tx, ty, tz
  CHARACTER(len=256)       :: filename
  INTEGER                  :: ios


  INCLUDE 'mpif.h'

  !----------------------------------------------------------------------------!
  ! Interpolate particle vorticity to mesh
  !----------------------------------------------------------------------------!
  CALL ppm_interp_p2m(xp,np,wp,lda,topo_id,mesh_id,&
       &              ppm_param_rmsh_kernel_mp4,&
       &              ghostsize, field_wp,info)

  !----------------------------------------------------------------------------!
  ! Ghost vorticity
  !----------------------------------------------------------------------------!
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)

  !----------------------------------------------------------------------------!
  ! Solve Poisson equation
  !----------------------------------------------------------------------------!
  CALL wvic_poisson_fft(info)

  !----------------------------------------------------------------------------!
  ! Get velocities and add U_infinity
  !----------------------------------------------------------------------------!
  SELECT CASE(wvic_compvel_scheme)
    CASE(0)
      CALL wvic_compvel(info)
    CASE(1)
      CALL wvic_compvel_4th(info)
    CASE(2)
      CALL wvic_compvel_spec(info)
  END SELECT

  CALL wvic_compvel_adduinfty(info)

  !JTR? is this ghost needed? It has been setup according to wvic_rhs_vort
  !----------------------------------------------------------------------------!
  ! Ghost vorticity
  !----------------------------------------------------------------------------!
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)

  !----------------------------------------------------------------------------!
  ! Ghost velocities
  !----------------------------------------------------------------------------!
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_up,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)

  !----------------------------------------------------------------------------!
  ! Calculate diagnostics
  !----------------------------------------------------------------------------!
  !JTR is moved outside again
  !CALL wvic_diagnostics (info)

  !-------------------------------------------------------------------------!
  ! Update solid velocity field - used for advection and penalization
  ! An array for storing the updated solid velocity field is needed.
  ! dwp is used for this which means that dwp cannot be correctly dumped now?
  !-------------------------------------------------------------------------!
  CALL cpu_time(time1) !JTR JTMP

  IF (flow_case .EQ. 11) THEN
    CALL wvic_stepfunc_medusa
!    u_solid = (/0.0_mk, -1.0_mk,0.0_mk/)
!    CALL wvic_solid_velocity(info) !REMOVE
  ELSE
    CALL wvic_solid_velocity(info)
  ENDIF

  !-------------------------------------------------------------------------!
  ! Advect the step function
  !-------------------------------------------------------------------------!
  IF (object_move .EQV. .TRUE.) THEN 
!JTR
write(*,*) 'flowers and sunshine'
    CALL wvic_advect_stepparticles(info)
  END IF

  CALL cpu_time(time2) !JTR JTMP
  CALL MPI_Reduce((time2-time1),time1,1,mpi_prec,MPI_MAX,0,comm,info)
  IF (rank .EQ. 0) THEN
    WRITE(msg,*) '\nitime=',itime,', Time used initialising and moving stepfunction: ',time1,'\n'
    WRITE(0,*) TRIM(msg)
  ENDIF

  !JTR? can this not be assumed constant? 
  CALL wvic_calculate_mass(info)
  IF (rank .EQ. 0) THEN
    WRITE(msg,*) '\nitime=',itime,', Medusa mass=',object_mass,'\n'
    WRITE(0,*) TRIM(msg)
    WRITE(msg,*) '\nitime=',itime,', Medusa center of mass',object_cmass,'\n'
    WRITE(0,*) TRIM(msg)
  ENDIF

!  lambda=penalization_lambda/dt

  IF(SUM(ghostsize).EQ.3) THEN
    fac1 = 0.5_mk / dx
    fac2 = 0.5_mk / dy
    fac3 = 0.5_mk / dz
  ELSE
     fac1 = 8.0_mk/dx/12.0_mk
     fac2 = 8.0_mk/dy/12.0_mk
     fac3 = 8.0_mk/dz/12.0_mk
     fac4 = -1.0_mk/dx/12.0_mk
     fac5 = -1.0_mk/dy/12.0_mk
     fac6 = -1.0_mk/dz/12.0_mk
  END IF

!tmp2=0.0_mk
!  DO isub=1,nsublist
!    isubl = isublist(isub)
!    DO k=1,ndata(3,isubl)-1
!      DO j=1,ndata(2,isubl)-1
!        DO i=1,ndata(1,isubl)-1
!          !tmp2=tmp2+(SUM(field_up(:,i,j,k,isub)-field_ubar(:,i,j,k,isub)))**2
!          tmp2=tmp2+field_H(i,j,k,isub)
!        END DO !i
!      END DO !j
!    END DO !k
!  END DO
!WRITE(msg,*) 'JTR: rank,itime,npc', rank, itime, npc,'\n'
!WRITE(0,*) msg


  !-------------------------------------------------------------------------!
  ! Rotation stuff/implicit penalization. Investigate and add comments
  !-------------------------------------------------------------------------!
  !JTR calculate lambda*dt outside loop
  DO isub=1,nsublist
    isubl = isublist(isub)

    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
!JTR use MIN .LT. 2  instead of SUM
          IF(SUM(ghostsize).EQ.3) THEN
        Auy   = (field_up(1,i,j-1,k,isub)+    lambda*dt*field_H(i,j-1,k,isub)* &
             & field_ubar(1,i,j-1,k,isub))/(1+lambda*dt*field_H(i,j-1,k,isub))
        Buy   = (field_up(1,i,j+1,k,isub)+    lambda*dt*field_H(i,j+1,k,isub)* &
             & field_ubar(1,i,j+1,k,isub))/(1+lambda*dt*field_H(i,j+1,k,isub))
        Auz   = (field_up(1,i,j,k-1,isub)+    lambda*dt*field_H(i,j,k-1,isub)* &
             & field_ubar(1,i,j,k-1,isub))/(1+lambda*dt*field_H(i,j,k-1,isub))
        Buz   = (field_up(1,i,j,k+1,isub)+    lambda*dt*field_H(i,j,k+1,isub)* &
             & field_ubar(1,i,j,k+1,isub))/(1+lambda*dt*field_H(i,j,k+1,isub))
        Avx   = (field_up(2,i-1,j,k,isub)+    lambda*dt*field_H(i-1,j,k,isub)* &
             & field_ubar(2,i-1,j,k,isub))/(1+lambda*dt*field_H(i-1,j,k,isub))
        Bvx   = (field_up(2,i+1,j,k,isub)+    lambda*dt*field_H(i+1,j,k,isub)* &
             & field_ubar(2,i+1,j,k,isub))/(1+lambda*dt*field_H(i+1,j,k,isub))
        Avz   = (field_up(2,i,j,k-1,isub)+    lambda*dt*field_H(i,j,k-1,isub)* &
             & field_ubar(2,i,j,k-1,isub))/(1+lambda*dt*field_H(i,j,k-1,isub))
        Bvz   = (field_up(2,i,j,k+1,isub)+    lambda*dt*field_H(i,j,k+1,isub)* &
             & field_ubar(2,i,j,k+1,isub))/(1+lambda*dt*field_H(i,j,k+1,isub))
        Awx   = (field_up(3,i-1,j,k,isub)+    lambda*dt*field_H(i-1,j,k,isub)* &
             & field_ubar(3,i-1,j,k,isub))/(1+lambda*dt*field_H(i-1,j,k,isub))
        Bwx   = (field_up(3,i+1,j,k,isub)+    lambda*dt*field_H(i+1,j,k,isub)* &
             & field_ubar(3,i+1,j,k,isub))/(1+lambda*dt*field_H(i+1,j,k,isub))
        Awy   = (field_up(3,i,j-1,k,isub)+    lambda*dt*field_H(i,j-1,k,isub)* &
             & field_ubar(3,i,j-1,k,isub))/(1+lambda*dt*field_H(i,j-1,k,isub))
        Bwy   = (field_up(3,i,j+1,k,isub)+    lambda*dt*field_H(i,j+1,k,isub)* &
             & field_ubar(3,i,j+1,k,isub))/(1+lambda*dt*field_H(i,j+1,k,isub))
            field_wp(1,i,j,k,isub) = fac2*(Bwy-Awy)-fac3*(Bvz-Avz)
            field_wp(2,i,j,k,isub) = fac3*(Buz-Auz)-fac1*(Bwx-Awx)
            field_wp(3,i,j,k,isub) = fac1*(Bvx-Avx)-fac2*(Buy-Auy)
          ELSE
        Auy   = (field_up(1,i,j-1,k,isub)+    lambda*dt*field_H(i,j-1,k,isub)* &
             & field_ubar(1,i,j-1,k,isub))/(1+lambda*dt*field_H(i,j-1,k,isub))
        Buy   = (field_up(1,i,j+1,k,isub)+    lambda*dt*field_H(i,j+1,k,isub)* &
             & field_ubar(1,i,j+1,k,isub))/(1+lambda*dt*field_H(i,j+1,k,isub))
        Auz   = (field_up(1,i,j,k-1,isub)+    lambda*dt*field_H(i,j,k-1,isub)* &
             & field_ubar(1,i,j,k-1,isub))/(1+lambda*dt*field_H(i,j,k-1,isub))
        Buz   = (field_up(1,i,j,k+1,isub)+    lambda*dt*field_H(i,j,k+1,isub)* &
             & field_ubar(1,i,j,k+1,isub))/(1+lambda*dt*field_H(i,j,k+1,isub))
        Avx   = (field_up(2,i-1,j,k,isub)+    lambda*dt*field_H(i-1,j,k,isub)* &
             & field_ubar(2,i-1,j,k,isub))/(1+lambda*dt*field_H(i-1,j,k,isub))
        Bvx   = (field_up(2,i+1,j,k,isub)+    lambda*dt*field_H(i+1,j,k,isub)* &
             & field_ubar(2,i+1,j,k,isub))/(1+lambda*dt*field_H(i+1,j,k,isub))
        Avz   = (field_up(2,i,j,k-1,isub)+    lambda*dt*field_H(i,j,k-1,isub)* &
             & field_ubar(2,i,j,k-1,isub))/(1+lambda*dt*field_H(i,j,k-1,isub))
        Bvz   = (field_up(2,i,j,k+1,isub)+    lambda*dt*field_H(i,j,k+1,isub)* &
             & field_ubar(2,i,j,k+1,isub))/(1+lambda*dt*field_H(i,j,k+1,isub))
        Awx   = (field_up(3,i-1,j,k,isub)+    lambda*dt*field_H(i-1,j,k,isub)* &
             & field_ubar(3,i-1,j,k,isub))/(1+lambda*dt*field_H(i-1,j,k,isub))
        Bwx   = (field_up(3,i+1,j,k,isub)+    lambda*dt*field_H(i+1,j,k,isub)* &
             & field_ubar(3,i+1,j,k,isub))/(1+lambda*dt*field_H(i+1,j,k,isub))
        Awy   = (field_up(3,i,j-1,k,isub)+    lambda*dt*field_H(i,j-1,k,isub)* &
             & field_ubar(3,i,j-1,k,isub))/(1+lambda*dt*field_H(i,j-1,k,isub))
        Bwy   = (field_up(3,i,j+1,k,isub)+    lambda*dt*field_H(i,j+1,k,isub)* &
             & field_ubar(3,i,j+1,k,isub))/(1+lambda*dt*field_H(i,j+1,k,isub))

        AAuy  = (field_up(1,i,j-2,k,isub)+    lambda*dt*field_H(i,j-2,k,isub)* &
             & field_ubar(1,i,j-2,k,isub))/(1+lambda*dt*field_H(i,j-2,k,isub))
        BBuy  = (field_up(1,i,j+2,k,isub)+    lambda*dt*field_H(i,j+2,k,isub)* &
             & field_ubar(1,i,j+2,k,isub))/(1+lambda*dt*field_H(i,j+2,k,isub))
        AAuz  = (field_up(1,i,j,k-2,isub)+    lambda*dt*field_H(i,j,k-2,isub)* &
             & field_ubar(1,i,j,k-2,isub))/(1+lambda*dt*field_H(i,j,k-2,isub))
        BBuz  = (field_up(1,i,j,k+2,isub)+    lambda*dt*field_H(i,j,k+2,isub)* &
             & field_ubar(1,i,j,k+2,isub))/(1+lambda*dt*field_H(i,j,k+2,isub))
        AAvx  = (field_up(2,i-2,j,k,isub)+    lambda*dt*field_H(i-2,j,k,isub)* &
             & field_ubar(2,i-2,j,k,isub))/(1+lambda*dt*field_H(i-2,j,k,isub))
        BBvx  = (field_up(2,i+2,j,k,isub)+    lambda*dt*field_H(i+2,j,k,isub)* &
             & field_ubar(2,i+2,j,k,isub))/(1+lambda*dt*field_H(i+2,j,k,isub))
        AAvz  = (field_up(2,i,j,k-2,isub)+    lambda*dt*field_H(i,j,k-2,isub)* &
             & field_ubar(2,i,j,k-2,isub))/(1+lambda*dt*field_H(i,j,k-2,isub))
        BBvz  = (field_up(2,i,j,k+2,isub)+    lambda*dt*field_H(i,j,k+2,isub)* &
             & field_ubar(2,i,j,k+2,isub))/(1+lambda*dt*field_H(i,j,k+2,isub))
        AAwx  = (field_up(3,i-2,j,k,isub)+    lambda*dt*field_H(i-2,j,k,isub)* &
             & field_ubar(3,i-2,j,k,isub))/(1+lambda*dt*field_H(i-2,j,k,isub))
        BBwx  = (field_up(3,i+2,j,k,isub)+    lambda*dt*field_H(i+2,j,k,isub)* &
             & field_ubar(3,i+2,j,k,isub))/(1+lambda*dt*field_H(i+2,j,k,isub))
        AAwy  = (field_up(3,i,j-2,k,isub)+    lambda*dt*field_H(i,j-2,k,isub)* &
             & field_ubar(3,i,j-2,k,isub))/(1+lambda*dt*field_H(i,j-2,k,isub))
        BBwy  = (field_up(3,i,j+2,k,isub)+    lambda*dt*field_H(i,j+2,k,isub)* &
             & field_ubar(3,i,j+2,k,isub))/(1+lambda*dt*field_H(i,j+2,k,isub))

            field_wp(1,i,j,k,isub) = fac2*(Bwy-Awy)-fac3*(Bvz-Avz)+ &
                                   & fac5*(BBwy-AAwy)-fac6*(BBvz-AAvz)
            field_wp(2,i,j,k,isub) = fac3*(Buz-Auz)-fac1*(Bwx-Awx)+ &
                                   & fac6*(BBuz-AAuz)-fac4*(BBwx-AAwx)
            field_wp(3,i,j,k,isub) = fac1*(Bvx-Avx)-fac2*(Buy-Auy)+ &
                                   & fac4*(BBvx-AAvx)-fac5*(BBuy-AAuy)
            END IF

        END DO
      END DO
    END DO
  END DO

  !----------------------------------------------------------------------------!
  ! Vorticity field has been updated on field_wp.
  ! For the remeshing particle vorticity is interpolated to the field_wp - 
  ! since field_wp already contains the new vorticity field the interpolation
  ! from particles to mesh in wvic_tvdrk3 should be deactivated in case
  ! penalization is done implicitly
  !----------------------------------------------------------------------------!
  !----------------------------------------------------------------------------!
  ! Ghost vorticity - this may be unnecessary - no I think it is nescessary
  !----------------------------------------------------------------------------!
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)

  !JTMP JTR
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))
  tz_center = 0.5_mk*(min_physg(3) + max_physg(3))
!JTR what is this?
 IF (rank .EQ. -10) THEN  
    WRITE(filename,'(A,I3.3,A)') 'vort',itime
    OPEN(14,file=filename,iostat=ios,position='append',status='new')
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO k=ndata(3,isubl),ndata(3,isubl)!1-ghostsize(3),ndata(3,isubl)+ghostsize(3)!
      DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)!ndata(2,isubl),ndata(2,isubl)!
        DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center
          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center
          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center
          WRITE(14,*) tx,ty,tz,field_wp(:,i,j,k,isub)!field_ubar(:,i,j,k,isub)! 
        END DO
      END DO
    END DO
  END DO
    CLOSE(14)
 ENDIF

END SUBROUTINE wvic_penalization_implicit

