!------------------------------------------------------------------------------!
!* filename: wvic_tvdrk3                                                      *!
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
!  $Log: wvic_tvdrk3.F,v $
!  Revision 1.13  2006/10/23 08:19:11  pchatela
!  LES models
!  Bugfixes in KE spectra and factor for Parseval identity
!  Removed the reset of noise in init_tvphysics_0 and_1
!
!  Revision 1.12  2006/10/16 17:24:21  pchatela
!  Projection to solenoidal vorticity field for the periodic case of tbphysics_0
!  spatial diags have to be computed for Flow_case 5 and 6 in the periodic case
!
!  Revision 1.11  2006/10/16 15:38:46  pchatela
!  Added divergence free (analytical) wvic_init_tvphysics_1.F
!
!  Revision 1.10  2006/10/03 16:11:58  pchatela
!  Added spectra in the z direction (misnamed kx spectrum...)
!  Added spatial diagnostics, like kinetic energy, enstrophy, circulation
!  as functions of z, dumped at the frequency ndump
!
!  Revision 1.9  2006/09/27 09:30:21  pchatela
!  Fixes, spectra calculation,
!  most importantly: moved the u_infty out, so it does not kill the dgammadt
!
!  Revision 1.8  2006/09/16 00:22:03  pchatela
!  Implemented the kinetic energy spectrum, dumped into an ascii file.
!
!  Revision 1.7  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.6  2006/09/01 15:47:03  pchatela
!  Added option to dump velocity and streamfunction fields
!
!  Revision 1.5  2006/08/30 17:42:20  pchatela
!  Bugfixes in initial conditions
!  Added dump at t=0
!
!  Revision 1.4  2006/08/29 13:42:55  menahel
!  uncommented call to enforcetv and added a conditional call to
!  ppm_interp_p2m_reno
!
!  Revision 1.3  2006/08/24 11:26:05  menahel
!  commented out the enforcetv0 call
!
!  Revision 1.2  2006/08/23 09:54:50  pchatela
!  Added calls to wvic_enforcetv0 in the case of the trailing vortex
!  formation
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.5  2006/02/06 16:35:33  pchatela
!  Added call to wvic_diags_mom
!
!  Revision 1.4  2005/11/21 17:36:39  michaebe
!  major session. now let''s debug
!
!  Revision 1.3  2005/11/11 17:16:35  michaebe
!  removed p2m m2p interfaces
!
!  Revision 1.2  2005/11/11 14:04:25  michaebe
!  clean up, additions, comments
!
!  Revision 1.1  2005/09/28 11:40:45  michaebe
!  Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!  run the simulation
!------------------------------------------------------------------------------!
SUBROUTINE wvic_tvdrk3 (niter, info)

  USE module_wvic
  USE wvic_module_io
  USE ppm_module_fft
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
  USE ppm_module_map_field_ghost !JTR
  USE MPI
  
  !----------------------------------------------------------------------------!
  !  Interfaces
  !----------------------------------------------------------------------------!
  INTERFACE
     SUBROUTINE wvic_alloc_field (vfield_up, info)
       USE module_wvic
       REAL (mk), DIMENSION (:, :, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field
     SUBROUTINE wvic_alloc_field_4 (vfield_up, info)
       USE module_wvic
       REAL (mk), DIMENSION (:, :, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_4
  END INTERFACE

  INTERFACE
     FUNCTION wvic_rhs_loca (vdummy, vxp, vup, vdime, vnp, &
          & ipack, lpack, rpack, info)
       USE module_wvic
       USE ppm_module_write
       INTEGER, INTENT(in)                        :: vdime, vnp
       REAL(mk),DIMENSION(:,:), POINTER           :: vdummy, vxp, vup
       REAL(MK),DIMENSION(:,:), POINTER,OPTIONAL  :: rpack
       INTEGER, DIMENSION(:,:), POINTER,OPTIONAL  :: ipack
       LOGICAL, DIMENSION(:,:), POINTER,OPTIONAL  :: lpack  
       INTEGER, INTENT(inout)            :: info
       INTEGER                                    :: wvic_rhs_loca
     END FUNCTION wvic_rhs_loca
  END INTERFACE

  INTERFACE
     FUNCTION wvic_rhs_vort ( vxp, vwp, vdwp, vdime, vnp, ipack, &
          &lpack, rpack, rkstep, info)
       USE module_wvic
       USE ppm_module_write
       INTEGER, INTENT(in)               :: vdime, vnp
       REAL(mk),DIMENSION(:,:),POINTER     :: vxp, vwp, vdwp
       REAL(MK),DIMENSION(:,:), POINTER,OPTIONAL  :: rpack
       INTEGER, DIMENSION(:,:), POINTER,OPTIONAL  :: ipack
       LOGICAL, DIMENSION(:,:), POINTER,OPTIONAL  :: lpack  
       INTEGER, INTENT(inout)            :: info
       INTEGER, INTENT(in)               :: rkstep
       INTEGER                                    :: wvic_rhs_vort
     END FUNCTION wvic_rhs_vort
  END INTERFACE
  
  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(inout) :: info
  INTEGER                :: niter ! howmany

  !----------------------------------------------------------------------------!
  !  Local/ODE
  !----------------------------------------------------------------------------!
  INTEGER                                  :: istage
  REAL(mk), PARAMETER :: onethird         = 0.3333333333333332_mk
  REAL(mk), PARAMETER :: twothird         = 0.6666666666666667_mk
  REAL(mk), PARAMETER :: fiveninth        = 0.5555555555555556_mk
  REAL(mk), PARAMETER :: fifteensixteenth = 0.9375_mk
  REAL(mk), PARAMETER :: onefivethree     = 1.1953125_mk
  REAL(mk), PARAMETER :: eightfifteenth   = 0.533333333333332_mk
  
  !----------------------------------------------------------------------------!
  !  Local/MAPPING
  !----------------------------------------------------------------------------!
  INTEGER                                  :: maptype

  !----------------------------------------------------------------------------!
  !  Local/MISC
  !----------------------------------------------------------------------------!
  CHARACTER(len=256) :: msg
  INTEGER            :: istat, p, isub, i,j,k, mpart,l,kp, isubl
  LOGICAL            :: ok
  REAL(mk)           :: tim1s, tim1e, timr1, timr2
  LOGICAL            :: abort = .FALSE.
  REAL(mk) :: cpu1, cpu2
  INTEGER :: clockspersec, t1_dum, t2_dum

  !----------------------------------------------------------------------------!
  ! time stuff
  !----------------------------------------------------------------------------!
  ! system clock
  INTEGER     :: t1_run_clock, t1_p2m_clock, t1_mgs_clock
  INTEGER     :: t2_run_clock, t2_p2m_clock, t2_mgs_clock
  ! mpi wtime
  REAL(mk)    :: t1_run_wtime, t1_p2m_wtime, t1_mgs_wtime
  REAL(mk)    :: t2_run_wtime, t2_p2m_wtime, t2_mgs_wtime
  ! cpu time
  REAL(mk)    :: t1_run_cpu, t1_p2m_cpu, t1_mgs_cpu
  REAL(mk)    :: t2_run_cpu, t2_p2m_cpu, t2_mgs_cpu
  REAL(mk)    :: runtime_now, runtime_tot
  LOGICAL     :: ink 
  !----------------------------------------------------------------------------!
  ! some dummy stuff
  !----------------------------------------------------------------------------!
  REAL(mk), DIMENSION(:,:), POINTER :: rpack
  INTEGER,  DIMENSION(:,:), POINTER :: ipack
  LOGICAL,  DIMENSION(:,:), POINTER :: lpack
  REAL(mk), DIMENSION(3)            :: len_phys

  !-----------------------------------------------------
  !  memory measurements
  !-----------------------------------------------------
  INTEGER    :: eye
  INTEGER*4  :: fragments
  INTEGER*8  :: total_free, largest_free, total_used
  INTEGER    :: heap_info,printmemsize

  !-----------------------------------------------------
  !  FFT stuff
  !-----------------------------------------------------
  INTEGER, DIMENSION(3) :: ftopo_id
  INTEGER, DIMENSION(4) :: fmesh_id, t_topoid
  
  !----------------------------------------------------------------------------!
  !  reset the timer
  !----------------------------------------------------------------------------!
  runtime_now = 0.0
  len_phys = max_physg - min_physg
  !----------------------------------------------------------------------------!
  !  allocate the buffers
  !----------------------------------------------------------------------------!
  NULLIFY(ipack,rpack,lpack)
  !----------------------------------------------------------------------------!
  !  Allocate temporary particle arrays (for multistep schemes)
  !----------------------------------------------------------------------------!
  ALLOCATE(xp0(dime,np),stat=istat)
  IF(istat.NE.0) THEN
     CALL ppm_write(rank,'wvic_tvdrk3','failed to allocate xp0',info)
  END IF
  ALLOCATE(wp0(lda,np),stat=istat)
  IF(istat.NE.0) THEN
     CALL ppm_write(rank,'wvic_tvdrk3','failed to allocate wp0',info)
  END IF
  
  IF(MOD(itime,ndump).EQ.0) THEN
        CALL MPI_Barrier(comm,info)
        CALL wvic_field2netcdf(info)
  END IF

  !----------------------------------------------------------------------------!
  ! = INTEGRATE =
  !----------------------------------------------------------------------------!
  DO WHILE(time.LT.tend)
     !-------------------------------------------------------------------------!
     !  start clock
     !-------------------------------------------------------------------------!
     CALL SYSTEM_CLOCK(t1_run_clock,clockspersec,t1_dum)
     t1_run_wtime = MPI_WTIME()
     CALL cpu_time(t1_run_cpu)
     !PPP----------------------------------------------------------------------!
     ! Here we do implicit penalization
     ! - diagnostics are calculated here as well
     ! - the step function is advected
     ! - solid velocity is computed
     !-------------------------------------------------------------------------!
     g_istage = 1 !JTR what is this?
     IF (penalization_implicit .EQV. .TRUE.) THEN
       CALL wvic_penalization_implicit(info)  !JTR
     END IF

! JTR I Have to interpolate the vorticity to the particles - and then back to the mesh in rhs_vort. Why? Jens might understand
!This is also done at the end of the iteration (still it is before the implicit penalization)
     CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,&
       & topo_id,mesh_id,(/cutoff,HUGE(cutoff)/),info,&
       & resetpos=.TRUE.,cutoff_weights=cow)
     !-------------------------------------------------------------------------!
     ! Vorticity field is ready for RK3
     ! !Vorticity should not be interpolated to mesh in rhs_vort first run
     !PPP----------------------------------------------------------------------!


     !-------------------------------------------------------------------------!
     !  Williamsons RK3
     !-------------------------------------------------------------------------!
!g_istage = 1 moved from here
     !-------------------------------------------------------------------------!
     ! Compute up and dwp
     !-------------------------------------------------------------------------!
     CALL measure(1,1)
     istat = wvic_rhs_vort(xp,wp,dwp,dime,np,ipack,lpack,rpack,1,info)
     CALL measure(2,1)
     istat = wvic_rhs_loca(xp,xp,up, dime,np,ipack,lpack,rpack,info)
     CALL measure(3,1)
!JTR clean up after the other guys. Or maybe I did this? Either way...
!explicit... what then?
!     IF ( ( (flow_case.EQ.5).OR. (flow_case.EQ.6)  .OR. (flow_case.EQ.7)).OR. (trailvortex) ) THEN
!        CALL wvic_diagnostics_trail (info)
!     ELSE
!        CALL wvic_diagnostics (info)
!     END IF
     CALL wvic_diagnostics (info)

     CALL measure(4,1)
     !-------------------------------------------------------------------------!
     ! Optional dumps
     !-------------------------------------------------------------------------

     IF(MOD(itime,ndump).EQ.0) THEN
	WRITE(msg,'(A)') 'H'
        CALL wvic_io_dump_field_vtk(field_H,ndata,msg, itime, info)
	WRITE(msg,'(A)') 'vel'
        CALL wvic_io_dump_field_vtk(field_up,ndata,msg, itime, info)
	WRITE(msg,'(A)') 'vrt'
        CALL wvic_io_dump_field_vtk(field_wp,ndata,msg, itime, info)
     END IF
     IF(MOD(itime,ndump).EQ.0) THEN
        CALL MPI_Barrier(comm,info)
        IF (dumpvel) THEN

!  DO isub=1,nsublist !JTR - old and abandoned. Used to output mask
!    DO k=1,ndata(3,isubl)
!      DO j=1,ndata(2,isubl)
!        DO i=1,ndata(1,isubl)
!          field_up(1,i,j,k,isub) = field_H(1,j,k,isub)
!          field_up(2,i,j,k,isub) = 0.0_mk
!          field_up(3,i,j,k,isub) = 0.0_mk
!        END DO !i
!      END DO !j
!    END DO !k
!  END DO !isub


!            WRITE(msg,'(A)') 'vel'
!            CALL wvic_io_dump_field_netcdf(field_up,ndata,msg, itime, info)
!            WRITE(msg,'(A)') 'dwp'
!            CALL wvic_io_dump_field_netcdf(field_dwp,ndata,msg, itime, info)
!            WRITE(msg,'(A)') 'H'
!            CALL wvic_io_dump_field_netcdf(field_H,ndata,msg, itime, info)
            ! use lvl instead of H
            WRITE(msg,'(A)') 'H'
            CALL wvic_io_dump_field_vtk(field_H,ndata,msg, itime, info)
            WRITE(msg,'(A)') 'vel'
            CALL wvic_io_dump_field_vtk(field_up,ndata,msg, itime, info)
            WRITE(msg,'(A)') 'vrt'
            CALL wvic_io_dump_field_vtk(field_wp,ndata,msg, itime, info)
        END IF
        IF (dumppsi) THEN
            WRITE(msg,'(A)') 'psi'
            CALL wvic_io_dump_field_netcdf(field_wps,ndata,msg, itime, info)
        END IF
        IF (dumpkespec) THEN
           ftopo_id = (/2,3,4/)
           t_topoid = (/2,3,4,5/)
           fmesh_id = (/2,3,4,5/)
           field_wp = field_up
           WRITE(msg,'(A,A,I5.5,A)') runtag(1:iruntag),'I',itime,'.kespec'
           CALL ppm_fft_kespectrum(field_wp,mesh_id,topo_id,t_topoid,fmesh_id, & 
             & ghostsize, info,msg)
        END IF
        IF (dumpkespeckx) THEN
           ftopo_id = (/2,3,4/)
           t_topoid = (/2,3,4,5/)
           fmesh_id = (/2,3,4,5/)
           field_wp = field_up
           WRITE(msg,'(A,A,I5.5,A)') runtag(1:iruntag),'I',itime,'.kespeckx'
           CALL ppm_fft_kespectrumkx(field_up,mesh_id,topo_id,t_topoid,fmesh_id, & 
             & ghostsize, info,msg)
        END IF
     END IF
     !-------------------------------------------------------------------------!
     ! store right hand side
     !-------------------------------------------------------------------------!
     DO kp=1,np
        xp0(1,kp) = up(1,kp)
        xp0(2,kp) = up(2,kp)
        xp0(3,kp) = up(3,kp)
        wp0(1,kp) = dwp(1,kp)
        wp0(2,kp) = dwp(2,kp)
        wp0(3,kp) = dwp(3,kp)
     END DO
     !-------------------------------------------------------------------------!
     ! euler step
     !-------------------------------------------------------------------------!
!print *, 'max_xp', maxval(xp)
     DO kp=1,np
        xp(1,kp) = xp(1,kp) + onethird*dt * up(1,kp)
        xp(2,kp) = xp(2,kp) + onethird*dt * up(2,kp)
        xp(3,kp) = xp(3,kp) + onethird*dt * up(3,kp)
        wp(1,kp) = wp(1,kp) + onethird*dt * dwp(1,kp)
        wp(2,kp) = wp(2,kp) + onethird*dt * dwp(2,kp)
        wp(3,kp) = wp(3,kp) + onethird*dt * dwp(3,kp)
     END DO
!print *, 'max_xp2', maxval(xp)
!print *, 'onethird', onethird
!print *, 'dt', dt
!print *, 'max_up', maxval(up)
!    write(6,*) 'time',time,'proc',rank,'memory',printmemsize()/1024/1024,'after euler'
!WRITE(0,*) 'JTRstop' !JTR
!call mpi_finalize(info)
!stop
!call wvic_died

     !-------------------------------------------------------------------------!
     ! time to map them
     !-------------------------------------------------------------------------!
     CALL wvic_pbc(info)
     CALL measure(5,1)
#ifndef __SMALL_FOOTPRINT__
     maptype = ppm_param_map_partial
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     ! ===   push   ====
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#else
     maptype = ppm_param_map_partial

     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(xp); NULLIFY(xp)
      maptype = ppm_param_map_send
     !-----------------------------------------------------
     !  measure memory here
     !-----------------------------------------------------
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(up); NULLIFY(up)
      maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(wp); NULLIFY(wp)
      maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(xp0); NULLIFY(xp0)
      maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(wp0,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(wp0); NULLIFY(wp0)
      maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp0,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(dwp); NULLIFY(dwp)
      maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#endif
!    write(6,*) 'time',time,'proc',rank,'memory',printmemsize()/1024/1024,'after small print'
     CALL measure(6,1)
     g_istage = 2
!     WRITE(msg,*) 'Williamson-',g_istage
!     IF(rank.EQ.0)      CALL ppm_write(rank,'wvic_tvdrk3',msg,info)
     !-------------------------------------------------------------------------!
     ! Compute up and dwp
     !-------------------------------------------------------------------------!
     CALL measure(1,2)
     istat = wvic_rhs_vort(xp,wp,dwp,dime,np,ipack,lpack,rpack,0,info)
     CALL measure(2,2)
     istat = wvic_rhs_loca(xp,xp,up, dime,np,ipack,lpack,rpack,info)
     CALL measure(3,2)
     CALL measure(4,2)
     !-------------------------------------------------------------------------!
     ! store q2
     !-------------------------------------------------------------------------!
     DO kp=1,np
        xp0(1,kp) = up(1,kp)  - fiveninth * xp0(1,kp)
        xp0(2,kp) = up(2,kp)  - fiveninth * xp0(2,kp)
        xp0(3,kp) = up(3,kp)  - fiveninth * xp0(3,kp)
        wp0(1,kp) = dwp(1,kp) - fiveninth * wp0(1,kp)
        wp0(2,kp) = dwp(2,kp) - fiveninth * wp0(2,kp)
        wp0(3,kp) = dwp(3,kp) - fiveninth * wp0(3,kp)
     END DO
     !-------------------------------------------------------------------------!
     ! do an euler step
     !-------------------------------------------------------------------------!
     DO kp=1,np
        xp(1,kp) = xp(1,kp) + fifteensixteenth*dt * xp0(1,kp)
        xp(2,kp) = xp(2,kp) + fifteensixteenth*dt * xp0(2,kp)
        xp(3,kp) = xp(3,kp) + fifteensixteenth*dt * xp0(3,kp)
        wp(1,kp) = wp(1,kp) + fifteensixteenth*dt * wp0(1,kp)
        wp(2,kp) = wp(2,kp) + fifteensixteenth*dt * wp0(2,kp)
        wp(3,kp) = wp(3,kp) + fifteensixteenth*dt * wp0(3,kp)
     END DO
     !-------------------------------------------------------------------------!
     ! time to map them
     !-------------------------------------------------------------------------!
     CALL wvic_pbc(info)
     CALL measure(5,2)
!     write(6,*) 'time',time,'proc',rank,'memory',printmemsize()/1024/1024,'after will1'
#ifndef __SMALL_FOOTPRINT__
     maptype = ppm_param_map_partial
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     ! ===   push   ====
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#else
     maptype = ppm_param_map_partial
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(xp); NULLIFY(xp)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(up); NULLIFY(up)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(wp); NULLIFY(wp)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(xp0); NULLIFY(xp0)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(xp0,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(wp0,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(wp0); NULLIFY(wp0)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp0,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(dwp); NULLIFY(dwp)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#endif
!    write(6,*) 'time',time,'proc',rank,'memory',printmemsize()/1024/1024,'after small2'
     CALL measure(6,2)

     g_istage = 3
!     WRITE(msg,*) 'Williamson-',g_istage
!     IF(rank.EQ.0)      CALL ppm_write(rank,'wvic_tvdrk3',msg,info)
     !-------------------------------------------------------------------------!
     ! Compute up and dwp
     !-------------------------------------------------------------------------!
     CALL measure(1,3)
     istat = wvic_rhs_vort(xp,wp,dwp,dime,np,ipack,lpack,rpack,0,info)
     CALL measure(2,3)
     istat = wvic_rhs_loca(xp,xp,up, dime,np,ipack,lpack,rpack,info)
     CALL measure(3,3)
     CALL measure(4,3)
     !-------------------------------------------------------------------------!
     ! compute q3
     !-------------------------------------------------------------------------!
     DO kp=1,np
        xp0(1,kp) = up(1,kp)  - onefivethree * xp0(1,kp)
        xp0(2,kp) = up(2,kp)  - onefivethree * xp0(2,kp)
        xp0(3,kp) = up(3,kp)  - onefivethree * xp0(3,kp)
        wp0(1,kp) = dwp(1,kp) - onefivethree * wp0(1,kp)
        wp0(2,kp) = dwp(2,kp) - onefivethree * wp0(2,kp)
        wp0(3,kp) = dwp(3,kp) - onefivethree * wp0(3,kp)
     END DO
     !-------------------------------------------------------------------------!
     ! final euler
     !-------------------------------------------------------------------------!
     DO kp=1,np
        xp(1,kp) = xp(1,kp) + eightfifteenth*dt * xp0(1,kp)
        xp(2,kp) = xp(2,kp) + eightfifteenth*dt * xp0(2,kp)
        xp(3,kp) = xp(3,kp) + eightfifteenth*dt * xp0(3,kp)
        wp(1,kp) = wp(1,kp) + eightfifteenth*dt * wp0(1,kp)
        wp(2,kp) = wp(2,kp) + eightfifteenth*dt * wp0(2,kp)
        wp(3,kp) = wp(3,kp) + eightfifteenth*dt * wp0(3,kp)
     END DO

     !-------------------------------------------------------------------------!
     ! Integrate medusa motion, move particles - explicit Euler
     !-------------------------------------------------------------------------!
     IF (flow_case .EQ. 11) THEN
       CALL wvic_medusa_move
     ENDIF

     !-------------------------------------------------------------------------!
     ! final map
     !-------------------------------------------------------------------------!
     !-------------------------------------------------------------------------!
     ! time to map them
     !-------------------------------------------------------------------------!
     CALL wvic_pbc(info)
     CALL measure(5,3)
!    write(6,*) 'time',time,'proc',rank,'memory',printmemsize()/1024/1024,'after will2'
#ifndef __SMALL_FOOTPRINT__
     maptype = ppm_param_map_partial
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     ! ===   push   ====
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#else
     maptype = ppm_param_map_partial
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(xp); NULLIFY(xp)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(up); NULLIFY(up)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(wp); NULLIFY(wp)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     DEALLOCATE(dwp); NULLIFY(dwp)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#endif
     CALL measure(6,3)
     
     time  = time + dt
     itime = itime + 1
     !-------------------------------------------------------------------------!
     ! time the remeshing
     !-------------------------------------------------------------------------!
     CALL SYSTEM_CLOCK(t1_p2m_clock,clockspersec,t1_dum)
     t1_p2m_wtime = MPI_WTIME()
     CALL cpu_time(t1_p2m_cpu)

 
     IF(trailvortex) THEN
        CALL ppm_interp_p2m_reno(xp,np,wp,lda,topo_id,mesh_id,krnl,ghostsize, &
             & field_wp,info)
     ELSE
        CALL ppm_interp_p2m(xp,np,wp,lda,topo_id,mesh_id,krnl,ghostsize, &
             & field_wp,info)
     END IF

     IF (trailvortex) THEN
        SELECT CASE(flow_case)
        CASE(0)
           !-----------------------------------------------------
           !  4 straight tubes
           !-----------------------------------------------------
           CALL wvic_enforcetv0
        CASE(1)
           !-----------------------------------------------------
           !  4 tubes with divergence-free perturbation 
           !-----------------------------------------------------
           CALL wvic_enforcetv1
        CASE(2)
           !-----------------------------------------------------
           !  4 tubes with divergence free perturbation at x=0
           !-----------------------------------------------------
           CALL wvic_enforcetv2
        END SELECT
     END IF
     
     IF (nkilldiv.GT.0) THEN
        IF (MOD(itime,nkilldiv).EQ.0) THEN 
        IF (.NOT.trailvortex) THEN
           ftopo_id = (/2,3,4/)
           t_topoid = (/2,3,4,5/)
           fmesh_id = (/2,3,4,5/)
           CALL ppm_fft_solenoidal(field_wp,mesh_id,topo_id,t_topoid,fmesh_id, &
             & ghostsize, info)
        END IF
        END IF
     END IF
     
     CALL measure(7,3)
     CALL SYSTEM_CLOCK(t2_p2m_clock,clockspersec,t2_dum)
     t2_p2m_wtime = MPI_WTIME()
     CALL cpu_time(t2_p2m_cpu)
     
     CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,&
          & topo_id,mesh_id,(/cutoff,HUGE(cutoff)/),info,&
          & resetpos=.TRUE.,cutoff_weights=cow)
     CALL measure(8,3)
! JTR? Now that vorticity is on both then new particles and the mesh why
! interpolate them back to the mesh at the beginning of the loop?
     !-------------------------------------------------------------------------!
     ! Reboot the other values, too
     !-------------------------------------------------------------------------!
     DEALLOCATE(up,dwp,xp0,wp0)
     ALLOCATE(up(dime,np),dwp(lda,np),xp0(dime,np),wp0(lda,np),stat=istat)
     IF(istat.NE.0) THEN
        CALL ppm_write(rank,'wvic_tvdrk3','failed to allocate dwp...',info)
     END IF
     up = 0.0_mk; dwp = 0.0_mk;
     
     CALL cpu_time(cpu2)
     tim1e = MPI_WTIME()
     CALL measure(9,3)
     !-------------------------------------------------------------------------!
     !  stop timers
     !-------------------------------------------------------------------------!
     CALL SYSTEM_CLOCK(t2_run_clock,clockspersec,t2_dum)
     t2_run_wtime = MPI_WTIME()
     CALL cpu_time(t2_run_cpu)
     !-------------------------------------------------------------------------!
     !  write timings
     !-------------------------------------------------------------------------!
     WRITE(msg,*) itime,' ',REAL(t2_run_clock-t1_run_clock)/REAL(clockspersec),&
          & ' ',(t2_run_wtime-t1_run_wtime),' ',(t2_run_cpu-t1_run_cpu)
     
     IF(rank.EQ.0.AND.MOD(itime,ndump).EQ.0) CALL ppm_write(rank,'wvic_tvdrk3'&
          &,msg,info)

     IF(MOD(itime,ndump).EQ.0) THEN
        CALL MPI_Barrier(comm,info)
        CALL wvic_field2netcdf(info)
     END IF
     
     runtime_now = runtime_now + (t2_run_wtime-t1_run_wtime)
     
     !-------------------------------------------------------------------------!
     !  check for abort file
     !-------------------------------------------------------------------------!
     CALL wvic_check_abort(abort,info)
     IF(rank.EQ.0) THEN
        IF(INT(runtime_now).GT.tot) THEN
           abort = .TRUE.
        END IF
     END IF
     IF(itime.EQ.niter) THEN
        abort = .TRUE.
     END IF
     CALL MPI_BCast(abort,1,MPI_LOGICAL,0,comm,info)
     !-------------------------------------------------------------------------!
     !  root broad casts
     !-------------------------------------------------------------------------!
     IF(abort) THEN
        CALL ppm_write(rank,'wvic_tvdrk3',&
             &'Abort file found. Terminizing.',info)
        CALL wvic_dump_restart(info)
        CALL MPI_Barrier(comm,info)
        EXIT
     END IF
        
     IF(itime.EQ.niter) EXIT
     
  END DO !end of iteration
  
  

  CALL MPI_Barrier(comm,info)
  
9999 CONTINUE
  IF(info.NE.0) CALL wvic_died
END SUBROUTINE wvic_tvdrk3


SUBROUTINE measure(tag,which)
  USE module_wvic
  INTEGER    :: tag,which
  
#ifdef __WITH_MEMORY_MEASUREMENT__
  INTEGER    :: eye
  INTEGER*4  :: fragments
  INTEGER*8  :: total_free, largest_free, total_used
  INTEGER    :: heap_info
  CHARACTER(len=256) :: filename


  WRITE(filename,'(A,I1.1,A)') 'memory',which,'.dat'
  OPEN(43,file=filename(1:LEN_TRIM(filename)),position='append',status='unknown')
  eye = heap_info(fragments, total_free, largest_free, total_used)
  WRITE(43,'(I,A,I,A,I)') itime,' ', tag,' ', total_used
  CLOSE(43)
#endif
END SUBROUTINE measure

  
  








