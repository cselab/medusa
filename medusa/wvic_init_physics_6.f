!-------------------------------------------------------------------------------
!          Homogeneous decaying turbulence test case
!-------------------------------------------------------------------------------
!* filename: wvic_init_physics_6                                              *!
!* project : ppm                                                              *!
!* purpose : initial condition for homogeneous periodic turbulence            *!
!*         :                                                                  *!
!* author  : Philippe Chatelain, Michael Bergdorf                             *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 09:45:17 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log $
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! = WVIC Vorticity =
!-------------------------------------------------------------------------------
SUBROUTINE wvic_init_physics_6

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_mg
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part
  USE ppm_module_fdsolver_solve
  USE ppm_module_fft
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! interfaces
  INTERFACE
     SUBROUTINE wvic_alloc_field_s (vfield_up, info)
       USE module_wvic
       IMPLICIT NONE
       REAL (mk), DIMENSION (:, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_s
  END INTERFACE
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! localities: noise stuff
  INTEGER, PARAMETER       :: mkd = KIND(1.0d0)
  INTEGER(mkd)             :: an, a32
  REAL(mk), DIMENSION(500) :: aky,bky,akz,bkz
  INTEGER                  :: nk,kk
  REAL(mk), DIMENSION(3)   :: RING_CENTER
  REAL(mk)                 :: rnk, amp
  INTEGER, DIMENSION(3)    :: ilo, ihi
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! localities: geometry stuff
  REAL(mk)                 :: maxvortl, maxvortg,scalf
  REAL(mk), DIMENSION(3)   :: sumvortl, sumvortg
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk)                 :: ink_a0, ink_b0
  REAL(mk)                 :: tz_ink1, ty_ink1
  REAL(mk)                 :: tz_ink2, ty_ink2
  REAL(mk)                 :: phixy, phixz, phiyx, phiyz, phizx, phizy
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! localities: stuff stuff
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: ig,jg,kg,i,j,k,isub,isubl,kmax
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! localities: mg stuff
  INTEGER,  DIMENSION(6)     :: ibcdef
  REAL(mk), DIMENSION(1,1,1) :: ibcvalue  
  REAL(mk)                     :: cresid
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! localities: helmholtz stuff
  REAL(mk)                     :: fac1,fac2,fac3,fac4,fac5,fac6,ldiv
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid
  
  !-----------------------------------------------------------------------------
  ! localities: timing counters
  INTEGER  :: c1,c2,c3,c4,c5,c6
  
  INCLUDE 'mpif.h'
  !-----------------------------------------------------------------------------
  ! = WHAT HAS TO BE DONE? =
  ! - generate noisy vorticity
  ! - reproject vorticity to get a divergence free field
  !-----------------------------------------------------------------------------

  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)

  nk =  wvic_noise_nmodes
  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_init_physics_6',msg,info)
  END IF
  
  !-----------------------------------------------------------------------------
  ! jump start the random number generator
  !-----------------------------------------------------------------------------
  a32 = 0_mkd
  a32 = IBSET(a32,32)
  an  = 1_mkd
  an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
  
  !-----------------------------------------------------------------------------
  ! Domain size
  !-----------------------------------------------------------------------------
  length = max_physg - min_physg
  
  !=================================================================
  ! VORTICITY PARAMETERS
  !-----------------------------------------------------------------
  IF (target_re.GT.0.0_mk) THEN
     max_vorticity = nu*target_re*(MINVAL(length))**2
     IF(rank.EQ.0) THEN
        WRITE(msg,*) 'From box Reynolds number, set vorticity amplitude to ', max_vorticity
        CALL ppm_write(rank,'wvic_init_physics_6',msg,info)
     END IF
  ELSEIF (wvic_noise_amp.GT.0.0_mk) THEN
     max_vorticity = wvic_noise_amp
     IF(rank.EQ.0) THEN
        WRITE(msg,*) 'From noise amplitude, set vorticity amplitude to ', max_vorticity
        CALL ppm_write(rank,'wvic_init_physics_6',msg,info)
     END IF
  END IF
  
  !-----------------------------------------------------
  ! check if we are restarting
  !-----------------------------------------------------
  field_wp = 0.0_mk
  IF(netcdf_restart) THEN
     CALL wvic_netcdf2field(info)
     time  = netcdf_time
     itime = netcdf_itime
     dt    = netcdf_dt
     IF(rank.EQ.0) THEN
        WRITE(msg,*) 'restarted from netcdf, itime = ',netcdf_itime
        CALL ppm_write(rank,'wvic_init_physics_6',msg,info)
     END IF
     GOTO 1122
  END IF
  
  !-----------------------------------------------------------------------------
  ! Careful, we want an initial condition independent of the nr of procs
  !-----------------------------------------------------------------------------
  DO isub=1,nsublist
     isubl = isublist(isub)
     !--------------------------------------------------------------------------
     ! Re-jump start the random number generator
     !--------------------------------------------------------------------------
     a32 = 0_mkd
     a32 = IBSET(a32,32)
     an  = 1_mkd
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     !--------------------------------------------------------------------------
     ! Loop over all points
     !--------------------------------------------------------------------------
     kg = 1
     jg = 1
     ig = 1
     DO 
        IF (kg.GE.(istart(3,isubl)+ndata(3,isubl))) THEN
           EXIT
        ELSEIF ( (ig.GE.istart(1,isubl)).AND. &
             &   (jg.GE.istart(2,isubl)).AND. &
             &   (kg.GE.istart(3,isubl)).AND. &
             &   (ig.LT.(istart(1,isubl)+ndata(1,isubl))).AND. &
             &   (jg.LT.(istart(2,isubl)+ndata(2,isubl))) ) THEN
           i = ig - istart(1,isubl) + 1
           j = jg - istart(2,isubl) + 1
           k = kg - istart(3,isubl) + 1
           field_wp(1,i,j,k,isub) = max_vorticity * REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
           an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
           field_wp(2,i,j,k,isub) = max_vorticity * REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
           an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
           field_wp(3,i,j,k,isub) = max_vorticity * REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
           an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
        ELSE
           ! Advance RNG
           an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
           an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
           an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
        END IF
        ! Advance global indices
        ig = ig + 1
        IF (ig.GT.nx(1)) THEN
           ig = 1
           jg = jg + 1
        END IF
        IF (jg.GT.nx(2)) THEN
           jg = 1
           kg = kg + 1
        END IF
     END DO
  END DO
  
  !-----------------------------------------------------
  !  compute velocity from this
  !-----------------------------------------------------
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
  
  CALL ppm_fft_velocities(field_up,&
       & mesh_id,topo_id,t_topoid,fmesh_id,ghostsize,info)
  !-----------------------------------------------------
  !  Set total momentum to zero
  !-----------------------------------------------------
  sumvortl = 0.0_mk
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)-1
        DO j=1,ndata(2,isubl)-1
           DO i=1,ndata(1,isubl)-1
              sumvortl(1) = sumvortl(1) + field_up(1,i,j,k,isub)
              sumvortl(2) = sumvortl(2) + field_up(2,i,j,k,isub)
              sumvortl(3) = sumvortl(3) + field_up(3,i,j,k,isub)
           END DO
        END DO
     END DO
  END DO
  CALL MPI_Allreduce(sumvortl,sumvortg,3,mpi_prec,MPI_SUM,comm,info)
  sumvortg = sumvortg/((nx(1))*(nx(2)-1)*(nx(3)-1))
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_up(1,i,j,k,isub) = field_up(1,i,j,k,isub) - sumvortg(1)
              field_up(2,i,j,k,isub) = field_up(2,i,j,k,isub) - sumvortg(2)
              field_up(3,i,j,k,isub) = field_up(3,i,j,k,isub) - sumvortg(3)
           END DO
        END DO
     END DO
  END DO
  !-----------------------------------------------------
  !  take curl of velocity
  !-----------------------------------------------------
  CALL wvic_ghost(wvic_prm_velocity,info)
  fac1 = 8.0_mk / dx / 12.0_mk
  fac2 = 8.0_mk / dy / 12.0_mk
  fac3 = 8.0_mk / dz / 12.0_mk
  fac4 = 1.0_mk / dx / 12.0_mk
  fac5 = 1.0_mk / dy / 12.0_mk
  fac6 = 1.0_mk / dz / 12.0_mk
  maxvortl = 0.0_mk
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              ! d u / d y
              phixy= fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))&
                   &-fac5*(field_up(1,i,j+2,k,isub)-field_up(1,i,j-2,k,isub))
              ! d u / d z
              phixz= fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))&
                   &-fac6*(field_up(1,i,j,k+2,isub)-field_up(1,i,j,k-2,isub))
              ! d v / d x
              phiyx= fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))&
                   &-fac4*(field_up(2,i+2,j,k,isub)-field_up(2,i-2,j,k,isub))
              ! d v / d z
              phiyz= fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))&
                   &-fac6*(field_up(2,i,j,k+2,isub)-field_up(2,i,j,k-2,isub))
              ! d w / d x
              phizx= fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))&
                   &-fac4*(field_up(3,i+2,j,k,isub)-field_up(3,i-2,j,k,isub))
              ! d w / d y
              phizy= fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))&
                   &-fac5*(field_up(3,i,j+2,k,isub)-field_up(3,i,j-2,k,isub))
              
              field_wp(1,i,j,k,isub) = phizy - phiyz
              field_wp(2,i,j,k,isub) = phixz - phizx
              field_wp(3,i,j,k,isub) = phiyx - phixy
              maxvortl = MAX(maxvortl,SUM(field_wp(:,i,j,k,isub)**2))
           END DO
        END DO
     END DO
  END DO
  
#if 0
  sumvortl = 0.0_mk
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)-1
        DO j=1,ndata(2,isubl)-1
           DO i=1,ndata(1,isubl)-1
              sumvortl(1) = sumvortl(1) + field_wp(1,i,j,k,isub)
              sumvortl(2) = sumvortl(2) + field_wp(2,i,j,k,isub)
              sumvortl(3) = sumvortl(3) + field_wp(3,i,j,k,isub)
           END DO
        END DO
     END DO
  END DO
  CALL MPI_Allreduce(sumvortl,sumvortg,3,mpi_prec,MPI_SUM,comm,info)
  sumvortg = sumvortg/((nx(1))*(nx(2)-1)*(nx(3)-1))
  maxvortl = 0.0_mk
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_wp(1,i,j,k,isub) = field_wp(1,i,j,k,isub) - sumvortg(1)
              field_wp(2,i,j,k,isub) = field_wp(2,i,j,k,isub) - sumvortg(2)
              field_wp(3,i,j,k,isub) = field_wp(3,i,j,k,isub) - sumvortg(3)
              maxvortl = MAX(maxvortl,SUM(field_wp(:,i,j,k,isub)**2))
           END DO
        END DO
     END DO
  END DO
#endif
  !-----------------------------------------------------
  !  Re-scale to achieve target re or maxvort
  !-----------------------------------------------------
  maxvortl = SQRT(maxvortl)
  CALL MPI_Allreduce(maxvortl,maxvortg,1,mpi_prec,MPI_MAX,comm,info)
  scalf = max_vorticity/maxvortg
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_wp(1,i,j,k,isub) = scalf*field_wp(1,i,j,k,isub) 
              field_wp(2,i,j,k,isub) = scalf*field_wp(2,i,j,k,isub) 
              field_wp(3,i,j,k,isub) = scalf*field_wp(3,i,j,k,isub) 
           END DO
        END DO
     END DO
  END DO
  
  
  
9212 CONTINUE
1122 CONTINUE
  !-----------------------------------------------------------------------------
  ! allocate auxilliary fields
  CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,topo_id,mesh_id,&
       & (/cutoff,HUGE(cutoff)/),info,resetpos=.TRUE.)

  WRITE(msg,*) ' created ',np,' vortobots'
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_physics_6',msg,info)
  !-----------------------------------------------------------------------------
  ! all set
  !-----------------------------------------------------------------------------
  ! we dont like plot3d no more
  CALL system_clock(c1,c2,c3)
  CALL wvic_field2netcdf(info)
  CALL system_clock(c4,c5,c6)
  IF(rank.eq.0) THEN
      write(msg,*) 'Writing netcdf took: ',REAL(c4-c1)/REAL(c5)
      CALL ppm_write(rank,'wvic_init_physics_6',msg,info)
  END IF
  
END SUBROUTINE wvic_init_physics_6