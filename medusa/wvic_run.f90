!-------------------------------------------------------------------------------
!* filename: wvic_run                                                          *!
!* project : ppm                                                              *!
!* purpose : run that jack                                                    *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Aug 18 10:36:19 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_run.F,v $
!  Revision 1.2  2006/10/25 11:35:50  menahel
!  re-enabled .rst dumping and added NRESTART parameter
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.1  2005/09/28 11:40:44  michaebe
!  Fork from ppm_pvc
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  run the simulation
!-------------------------------------------------------------------------------
SUBROUTINE wvic_run (niter, info)

  USE module_wvic
  !--- ppm general
  USE ppm_module_data
  USE ppm_module_write
  !--- ode solver
  USE ppm_module_ode_step
  USE ppm_module_ode_init
  USE ppm_module_ode_alldone
  USE ppm_module_ode_create_ode
  USE ppm_module_ode_map_pop
  USE ppm_module_ode_map_push
  USE ppm_module_ode_start
  !--- remeshing
  USE ppm_module_rmsh_remesh
  USE ppm_module_rmsh_comp_weights
  USE ppm_module_rmsh_create_part
  USE ppm_module_interp_p2m
  !--- mapping
  USE ppm_module_map_part
  USE ppm_module_topo_check

  !-----------------------------------------------------------------------------
  !  Interfaces
  !-----------------------------------------------------------------------------
  INTERFACE
     SUBROUTINE wvic_m2p (vxp, vup, vfield_up, info)
       USE module_wvic
       USE ppm_module_data_rmsh
       USE ppm_module_write
       REAL (mk), DIMENSION (:, :), POINTER :: vxp
       REAL (mk), DIMENSION (:, :), POINTER :: vup
       REAL (mk), DIMENSION (:, :, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_m2p
  END INTERFACE

  INTERFACE
     SUBROUTINE wvic_p2m (vxp, vup, info)
       USE module_wvic
       USE ppm_module_rmsh_comp_weights
       USE ppm_module_rmsh_remesh
       USE ppm_module_write
       REAL (mk), DIMENSION (:, :), POINTER :: vxp
       REAL (mk), DIMENSION (:, :), POINTER :: vup
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_p2m
  END INTERFACE
  

  
  
  !-----------------------------------------------------------------------------
  !  Arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(inout) :: info
  INTEGER                :: niter ! howmany
#ifdef __asdfkjas
  !-----------------------------------------------------------------------------
  !  Local/ODE
  !-----------------------------------------------------------------------------
  REAL(mk),  DIMENSION(4)                  :: ode_w_time, ode_x_time
  INTEGER                                  :: ode_w_ibfr, ode_x_ibfr
  INTEGER                                  :: ode_w_id  , ode_x_id
  LOGICAL                                  :: ode_adaptive
  REAL(mk),  DIMENSION(:,:)  , POINTER     :: ode_w_bfr,  ode_x_bfr
  INTEGER ,  EXTERNAL                      :: wvic_rhs_vort, wvic_rhs_loca
  INTEGER                                  :: ode_w_nstage, ode_x_nstage
  INTEGER                                  :: ode_w_scheme, ode_x_scheme
  INTEGER                                  :: istage
  LOGICAL ,  DIMENSION(:,:),   POINTER     :: ode_lpack
  INTEGER ,  DIMENSION(:,:),   POINTER     :: ode_ipack
  REAL(mk),  DIMENSION(:,:),   POINTER     :: ode_rpack

  !-----------------------------------------------------------------------------
  !  Local/MAPPING
  !-----------------------------------------------------------------------------
  INTEGER                                  :: maptype

  !-----------------------------------------------------------------------------
  !  Local/MISC
  !-----------------------------------------------------------------------------
  CHARACTER(len=256) :: msg
  INTEGER            :: istat, p, isub, i,j,k, mpart
  LOGICAL            :: ok
  REAL(mk)           :: tim1s, tim1e, timr1, timr2
  LOGICAL            :: abort = .FALSE.
  INCLUDE 'mpif.h'
  REAL(mk) :: cpu1, cpu2
  INTEGER :: clockspersec, t1_dum, t2_dum
  !-----------------------------------------------------
  ! time stuff
  !-----------------------------------------------------
  ! system clock
  INTEGER :: t1_run_clock, t1_p2m_clock, t1_mgs_clock
  INTEGER :: t2_run_clock, t2_p2m_clock, t2_mgs_clock
  ! mpi wtime
  REAL(mk)    :: t1_run_wtime, t1_p2m_wtime, t1_mgs_wtime
  REAL(mk)    :: t2_run_wtime, t2_p2m_wtime, t2_mgs_wtime
  ! cpu time
  REAL(mk)    :: t1_run_cpu, t1_p2m_cpu, t1_mgs_cpu
  REAL(mk)    :: t2_run_cpu, t2_p2m_cpu, t2_mgs_cpu
  REAL(mk)    :: runtime_now, runtime_tot
  
  
  !-----------------------------------------------------------------------------
  !  Nullify all empty pointers for the sake of certainty
  !-----------------------------------------------------------------------------
  NULLIFY(wx,wy,wz)
  NULLIFY(ode_x_bfr, ode_w_bfr, ode_lpack, ode_ipack, ode_rpack)
  ALLOCATE(ode_lpack(1,1))
  ode_lpack(1,1) = .TRUE.
  runtime_now = 0.0
  !-----------------------------------------------------------------------------
  !  Initialize the Ode solver
  !-----------------------------------------------------------------------------
  CALL ppm_ode_init (info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_run','Initialization of ode solver failed',&
          & info)
     GOTO 9999
  END IF


  !-----------------------------------------------------------------------------
  !  Create the modes
  !-----------------------------------------------------------------------------
  ode_adaptive = .FALSE.
  ode_w_id = -1
  ode_x_id = -1
  !-----------------------------------------------------
  ! mid point rk
  ! ode_w_scheme = ppm_param_ode_scheme_midrk2
  ! ode_x_scheme = ppm_param_ode_scheme_midrk2
  !-----------------------------------------------------

  !-----------------------------------------------------
  ! trapezoidal
  ode_w_scheme = ppm_param_ode_scheme_tvdrk3
  ode_x_scheme = ppm_param_ode_scheme_tvdrk3
  !-----------------------------------------------------
  
  CALL ppm_ode_create_ode(ode_w_id,ode_w_ibfr,ode_w_nstage,&
       & ode_w_scheme, ode_w_scheme, ode_adaptive, info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_run','failed to create a mode for wp',info)
     GOTO 9999
  END IF

  CALL ppm_ode_create_ode(ode_x_id,ode_x_ibfr,ode_x_nstage,&
       & ode_x_scheme, ode_x_scheme, ode_adaptive, info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_run','failed to create a mode for xp',info)
     GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  !  Allocate memory for the buffers
  !-----------------------------------------------------------------------------
  ALLOCATE(ode_x_bfr(ode_x_ibfr*dime,np),&
       &   ode_w_bfr(ode_w_ibfr*lda,np), stat=istat)
  IF(istat.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_run','failed to allocate buffers',info)
     GOTO 9999
  END IF

  ode_x_bfr = 0.0_mk ; ode_w_bfr = 0.0_mk

  !-----------------------------------------------------------------------------
  !  reset the timer
  !-----------------------------------------------------------------------------
  ode_w_time(1) = time             ! start time
  ode_w_time(2) = tend
  ode_w_time(3) = time             ! current time
  ode_w_time(4) = dt
  ode_x_time    = ode_w_time



  !-----------------------------------------------------------------------------
  !  Start the ode solver (not start "solving")
  !-----------------------------------------------------------------------------
  CALL ppm_ode_start(info)
  IF(info.NE.0) THEN
     if(verbose) CALL ppm_write(rank,'wvic_run','failed to start ode solver',info)
     GOTO 9999
  END IF
  
  
  !-----------------------------------------------------------------------------
  !  Entering the loops of timne
  !-----------------------------------------------------------------------------
  DO WHILE(.NOT.ppm_ode_alldone(info))
     !-----------------------------------------------------
     !  start clock
     !-----------------------------------------------------
     CALL SYSTEM_CLOCK(t1_run_clock,clockspersec,t1_dum)
     t1_run_wtime = MPI_WTIME()
     CALL cpu_time(t1_run_cpu)

     !--------------------------------------------------------------------
     !  Loop over stages
     !-------------------------------------------------------------------
     DO istage=1,MAX(ode_x_nstage,ode_w_nstage)
        g_istage = istage
        !------------------------------------------------------------------
        !  call a step for vorticity
        !-------------------------------------------------------------------
        CALL ppm_ode_step(ode_w_id,xp,wp,dwp,lda,np,ode_w_bfr,istage,&
             & ode_w_time, wvic_rhs_vort, info=info)
        IF(info.NE.0) THEN
           IF(verbose) THEN
              CALL ppm_write(rank,'wvic_run','ode step for wp failed',info)
           END IF
           GOTO 9999
        END IF
        !--------------------------------------------------------------------
        !  call a step for the locations
        !------------------------------------------------------------------
        CALL ppm_ode_step(ode_x_id,xp,xp,up,dime,np,ode_x_bfr,istage,&
             & ode_x_time, wvic_rhs_loca, info=info, lpackdata=ode_lpack)
        IF(info.NE.0) THEN
           IF(verbose) THEN
              CALL ppm_write(rank,'wvic_run','ode step for xp failed',info)
           END IF
           GOTO 9999
        END IF
        
        IF(verbose) CALL ppm_write(rank,'wvic_run','going to map now',info)
        
        !------------------------------------------------------------------
        !  have to remap the particles
        !-------------------------------------------------------------------
        CALL wvic_pbc(info)
        
        !-----------------------------------------------------
        !  measure mapping time
        !-----------------------------------------------------
        CALL SYSTEM_CLOCK(t1_p2m_clock,clockspersec,t1_dum)
        t1_p2m_wtime = MPI_WTIME()
        CALL cpu_time(t1_p2m_cpu)
        
        maptype = ppm_param_map_partial
        CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
        ! ===   push   ====
        maptype = ppm_param_map_push
        CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
        CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
        CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
        ! ===   push bfr   ====
        CALL ppm_ode_map_push(ode_w_id,ode_w_bfr,lda,np,mpart,info)
        CALL ppm_ode_map_push(ode_x_id,ode_x_bfr,dime,np,mpart,info)
        ! ===   send   ====
        maptype = ppm_param_map_send
        CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
        ! ===   pop bfr    ====
        CALL ppm_ode_map_pop(ode_x_id,ode_x_bfr,dime,np,mpart,info)
        CALL ppm_ode_map_pop(ode_w_id,ode_w_bfr,lda,np,mpart,info)
        ! ===   pop   ====
        maptype = ppm_param_map_pop
        CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
        CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
        CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
        CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
        np = mpart
        CALL SYSTEM_CLOCK(t2_p2m_clock,clockspersec,t2_dum)
        t2_p2m_wtime = MPI_WTIME()
        CALL cpu_time(t2_p2m_cpu)
        WRITE(msg,*) itime,' ',REAL(t2_p2m_clock-t1_p2m_clock)&
             &/REAL(clockspersec),&
             & ' ',(t2_p2m_wtime-t1_p2m_wtime),' ',(t2_p2m_cpu-t1_p2m_cpu)
        CALL ppm_write(rank,'wvic_part_maps',msg,info)
        
           
        IF(info.NE.0) THEN
           IF(verbose) THEN
              CALL ppm_write(rank,'wvic_run','map partial failed',info)
              GOTO 9999
           END IF
        END IF
        !----------------------------------------------------------------
        !  check that particles have been mapped correctly
        !---------------------------------------------------------------
        ok = .FALSE.
        CALL ppm_topo_check(xp,np,ok,info)
        IF(.NOT.ok) THEN
           IF(verbose) THEN
              CALL ppm_write(rank,'wvic_run','mapping is wrong.'&
                   & ,info)
           END IF
           GOTO 9999
        END IF
     END DO

        
     !--------------------------------------------------------------------------
     !  do some diagnostics
     !--------------------------------------------------------------------------
     
     itime = itime + 1
     time  = ode_w_time(3)
     CALL wvic_diagnostics (info)
     ode_w_time(4) = dt
     ode_x_time(4) = dt
     
     !--------------------------------------------------------------------------
     !  Time to remesh them guys (dont tell Jens)
     !--------------------------------------------------------------------------
     IF(ASSOCIATED(wx)) THEN
        DEALLOCATE(wx,wy,wz,stat=istat)
        IF(istat.NE.0) THEN
           if(verbose) CALL ppm_write(rank,'wvic_run','failed to dealloc weights',info)
        END IF
        NULLIFY(wx,wy,wz)
     ELSE
        NULLIFY(wx,wy,wz)
     END IF

     IF(.NOT.fast_rmsh) THEN     
        WRITE(0,*)'JTRrun' !JTR
        CALL ppm_rmsh_comp_weights(xp,np,topo_id,mesh_id,krnl,info, &
             & wx1_user=wx,wx2_user=wy,wx3_user=wz)
        IF(info.NE.0) THEN
           IF(verbose) CALL ppm_write(rank,'wvic_run','(rmsh) comp.weights. failed',info)
           GOTO 9999
        END IF
        
        CALL ppm_rmsh_remesh(xp,np,wp,lda,topo_id,mesh_id,krnl,ghostsize, &
             & info, field_up=field_wp, wx1_user=wx, wx2_user=wy, wx3_user=wz)!,&
        !          & lrnrm=.false.)
        IF(info.NE.0) THEN
           IF(verbose) CALL ppm_write(rank,'wvic_run','(rmsh) remesh. failed',info)
           GOTO 9999
        END IF
        
        CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,&
             & topo_id,mesh_id,(/cutoff,huge(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
        IF(info.NE.0) THEN
           IF(verbose) CALL ppm_write(rank,'wvic_run','(rmsh) create parts. failed',info)
           GOTO 9999
        END IF
        !--------------------------------------------------------------------------
        !  kill  the weights
        !--------------------------------------------------------------------------
        IF(ASSOCIATED(wx))  DEALLOCATE(wx,wy,wz,stat=istat)
        IF(istat.NE.0) THEN
           IF(verbose) CALL ppm_write(rank,'wvic_run','failed to dealloc weights',info)
        END IF
        NULLIFY(wx,wy,wz)
     ELSE

        !-----------------------------------------------------
        ! time the remeshing
        !-----------------------------------------------------
        CALL SYSTEM_CLOCK(t1_p2m_clock,clockspersec,t1_dum)
        t1_p2m_wtime = MPI_WTIME()
        CALL cpu_time(t1_p2m_cpu)
        
        CALL ppm_interp_p2m(xp,np,wp,lda,topo_id,mesh_id,krnl,ghostsize, &
             & info, field_up=field_wp)

        CALL SYSTEM_CLOCK(t2_p2m_clock,clockspersec,t2_dum)
        t2_p2m_wtime = MPI_WTIME()
        CALL cpu_time(t2_p2m_cpu)
        WRITE(msg,*) itime,' ',REAL(t2_p2m_clock-t1_p2m_clock)&
             &/REAL(clockspersec),&
             & ' ',(t2_p2m_wtime-t1_p2m_wtime),' ',(t2_p2m_cpu-t1_p2m_cpu)
        CALL ppm_write(rank,'wvic_p2m',msg,info)

        CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,&
             & topo_id,mesh_id,(/cutoff,huge(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
        
        
     END IF

     
     !--------------------------------------------------------------------------
     !  now hack the ode solver so we dont have to restart it
     !--------------------------------------------------------------------------
     DEALLOCATE(ode_w_bfr,ode_x_bfr,up,dwp)
     IF(istat.NE.0) THEN
        IF(verbose) CALL ppm_write(rank,'wvic_run','failed to dealloc buffers',info)
     END IF
     ALLOCATE(ode_x_bfr(dime*ode_x_ibfr,np), dwp(lda,np),  &
          &   ode_w_bfr(lda*ode_w_ibfr,np), up(dime,np), stat=istat)
     IF(istat.NE.0) THEN
        IF(verbose) CALL ppm_write(rank,'wvic_run','failed to allocate new buffers',info)
     END IF
     up = 0.0_mk



     CALL cpu_time(cpu2)
     tim1e = MPI_WTIME()

     IF(verbose) WRITE(msg,*) itime,' time ',ode_w_time(3),' took ',INT(1000.0*(tim1e-tim1s))
     IF(verbose) CALL ppm_write(rank,'wvic_run',msg,info)

     !-----------------------------------------------------
     !  stop timers
     !-----------------------------------------------------
     CALL SYSTEM_CLOCK(t2_run_clock,clockspersec,t2_dum)
     t2_run_wtime = MPI_WTIME()
     CALL cpu_time(t2_run_cpu)

     !-----------------------------------------------------
     !  write timings
     !-----------------------------------------------------
     WRITE(msg,*) itime,' ',REAL(t2_run_clock-t1_run_clock)/REAL(clockspersec),&
          & ' ',(t2_run_wtime-t1_run_wtime),' ',(t2_run_cpu-t1_run_cpu)
     
     CALL ppm_write(rank,'wvic_run',msg,info)
     IF(MOD(itime,ndump).EQ.0) THEN
        CALL wvic_field2netcdf(info)
     END IF
     IF(MOD(itime,nrestart).EQ.0) THEN
        CALL wvic_dump_restart(info)
     END IF

     runtime_now = runtime_now + (t2_run_wtime-t1_run_wtime)
     
     
     !--------------------------------------------------------------------------
     !  check for abort file
     !--------------------------------------------------------------------------
     CALL wvic_check_abort(abort,info)
     ! run for 90 minutes maximum
!     IF(runtime_now.GT.60.0) THEN
!        abort = .TRUE.
!     END IF
     CALL MPI_BCast(abort,1,MPI_LOGICAL,0,comm,info)
     !-----------------------------------------------------
     !  root broad casts
     !-----------------------------------------------------
     IF(abort) THEN
        CALL ppm_write(rank,'wvic_run',&
             &'Abort file found. Terminizing.',info)
        EXIT
     END IF

     IF(itime.EQ.niter) EXIT
  END DO



  CALL wvic_dump_restart(info)
  CALL MPI_Barrier(comm,info)
  DEALLOCATE(ode_lpack)
  
9999 CONTINUE
#endif
  IF(info.NE.0) CALL wvic_died
  
END SUBROUTINE wvic_run
