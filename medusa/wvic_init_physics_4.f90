!  this is the rings for the CG lab
!-------------------------------------------------------------------------------
!* filename: pvc_vorticity                                                    *!
!* project : ppm                                                              *!
!* purpose : initial condition for elliptical vortex ring                     *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 09:45:17 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_init_physics_4.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.4  2005/03/19 00:11:54  michaebe
!  adopted to the new create_part
!
!  Revision 1.3  2005/01/06 12:28:40  michaebe
!  all the changes needed for the tubes
!
!  Revision 1.2  2004/12/03 12:46:25  michaebe
!  forgot to deallocate potential and divergence
!
!  Revision 1.1  2004/12/03 12:08:59  michaebe
!  initial implementation of crow instability stuff
!
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! = PVC Vorticity =
! Generate the crow instability initial condition which is two antiparallel
! vortices perturbed by noise
!-------------------------------------------------------------------------------
SUBROUTINE wvic_init_physics_4

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
  REAL(mk)                 :: a0, b0, gamma
  REAL(mk)                 :: tx_center, tx_vtx1, tx_vtx2
  REAL(mk)                 :: ty_center, ty_vtx1, ty_vtx2, r_vtx1, r_vtx2
  REAL(mk)                 :: tz_center, tz_vtx1, tz_vtx2
  REAL(mk)                 :: fac1, ra0, tx, ty, tz, noise, typ, tzp, omega_x
  REAL(mk)                 :: ra0t, fac1t
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
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
  REAL(mk)                     :: fac2,fac3,fac4,fac5,fac6,ldiv
  REAL(mk), DIMENSION(:,:,:,:), POINTER :: divergence, potential
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid
  !-----------------------------------------------------------------------------
  ! localities: Vortex  Ring stuff
  REAL(mk) :: rad1t, rad1r, theta1, radstr, rad1sq, radstrength, vrrad, rletrad
  REAL(mk) :: rpos, pertu,rad1sqTILDA,vsrad,ringlet_radius
  INTEGER  :: imirr,jmirr,kmirr,c1,c2,c3,c4,c5,c6
  !-----------------------------------------------------------------------------
  ! = WHAT HAS TO BE DONE? =
  ! - generate noise (angle of dislocation and amplitute of dislocation)
  ! - calculate the position of the vortices using b0 (intervortex distance)
  ! - initialize the vortices using a0 (size of vortices)
  ! - reproject vorticity to get a divergence free field
  !-----------------------------------------------------------------------------

  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)

  nk =  wvic_noise_nmodes
  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_init_physics_4',msg,info)
  END IF

  !-----------------------------------------------------------------------------
  ! jump start the random number generator
  !-----------------------------------------------------------------------------
  a32 = 0_mkd
  a32 = IBSET(a32,32)
  an  = 1_mkd
  an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
  !-----------------------------------------------------------------------------
  ! print 1st seed as debug
  !-----------------------------------------------------------------------------
  DO i=1,nk

     aky(i) = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     bky(i) = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     akz(i)   = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     bkz(i)   = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)

  END DO

  !=================================================================
  ! VORTEX RING PARAMETERS
  !-----------------------------------------------------------------
  !  z position of ring, middle of the box
  !  outer radius of ring
  !  roughly 2.6cm
  vrrad  = 0.013
  !  inner radius of ring
  !  used to be vrrad / 3.5
  IF(tube_radius.GT.0.0_mk) THEN
     vsrad = vrrad / tube_radius
  ELSE
     vsrad  = vrrad / 2.379
  END IF
  !  circulation
  IF(target_re.GT.0.0_mk) THEN
     gamma = nu * target_re
  ELSE
     gamma  = 0.8415_mk
  END IF
  !*****************************************************************
  !-----------------------------------------------------------------------------
  ! compute the random amplitudes
  !-----------------------------------------------------------------------------
  nk = wvic_noise_nmodes
  IF(nk.GT.0) THEN
     rnk = 1.0_mk/(2.0_mk*REAL(nk,mk))
  ELSE
     rnk = 0.0_mk
  END IF

  amp = wvic_noise_amp * rnk * vrrad
  

  !  derived parameters
  radstr = 1.0_MK/vsrad**2;
  length = max_physg - min_physg
  
  !-----------------------------------------------------
  ! check if we are restarting
  !-----------------------------------------------------
  IF(netcdf_restart) THEN
     CALL wvic_netcdf2field(info)
     time  = netcdf_time
     itime = netcdf_itime
     dt    = netcdf_dt
     GOTO 1122
  END IF
  
  !-----------------------------------------------------------------------------
  ! create vortices
  !-----------------------------------------------------------------------------
  max_vorticity = gamma*radstr/M_PI
  
  RING_CENTER = 0.5_mk*(min_physg+max_physg)
  RING_CENTER(3) = min_physg(3)
#include "wvic_put_a_ring.inc"
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
  !  take curl of velocity
  !-----------------------------------------------------
  CALL wvic_ghost(wvic_prm_velocity,info)
  fac1 = 8.0_mk / dx / 12.0_mk
  fac2 = 8.0_mk / dy / 12.0_mk
  fac3 = 8.0_mk / dz / 12.0_mk
  fac4 = 1.0_mk / dx / 12.0_mk
  fac5 = 1.0_mk / dy / 12.0_mk
  fac6 = 1.0_mk / dz / 12.0_mk
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
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_physics_4',msg,info)
  CALL system_clock(c1,c2,c3)
  CALL wvic_field2netcdf(info)
  CALL system_clock(c4,c5,c6)
  IF(rank.eq.0) THEN
      write(msg,*) 'Writing netcdf took: ',REAL(c4-c1)/REAL(c5)
      CALL ppm_write(rank,'wvic_init_physics_4',msg,info)
  END IF
  
END SUBROUTINE wvic_init_physics_4
