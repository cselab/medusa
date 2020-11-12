!------------------------------------------------------------------------------!
!* filename: wvic_init_sphere                                                 *!
!* project : ppm                                                              *!
!* purpose : initial condition for helical vortex instability                 *!
!*             &                                                              *!
!* purpose : impose zero vorticity and (velocity) inside object               *!
!*         :                                                                  *!
!* author  : Johannes Tophoej Rasmussen                                       *!
!*         : Technical University of Denmark (DTU)                            *!
!*         :                                                                  *!
!* date    : Sep 2007                                                         *!
!* please return to <johs@johs.nu>                                            *!
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!  initialize vorticity field for stokes flow around sphere
!------------------------------------------------------------------------------!
SUBROUTINE wvic_init_sphere_stokes

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft

  !----------------------------------------------------------------------------!
  ! interfaces
  INTERFACE
     SUBROUTINE wvic_alloc_field_s (vfield_up, info)
       USE module_wvic
       REAL (mk), DIMENSION (:, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_s
     SUBROUTINE wvic_alloc_field(vfield_up, ilda, info)
       USE module_wvic
       REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: vfield_up
       INTEGER                   , INTENT(out) :: info
       INTEGER                   , INTENT(in ) :: ilda
     END SUBROUTINE wvic_alloc_field
  END INTERFACE
  !----------------------------------------------------------------------------!
  REAL(MK), EXTERNAL :: rtbis, rtflsp, rtsec, zriddr, zbrent, sisec, sisec2


  !----------------------------------------------------------------------------!
  ! localities: noise stuff
!  INTEGER, PARAMETER       :: mkd = KIND(1.0d0)
!  INTEGER(mkd)             :: an, a32 !&
!  REAL(mk), DIMENSION(500) :: aky,bky,akz,bkz !&
!  INTEGER                  :: nk,kk !&
!  REAL(mk)                 :: rnk, amp !&
!  INTEGER, DIMENSION(3)    :: ilo, ihi !&
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! phi is the phase between the helices
  ! Rs - radius of sphere
  ! R, theta, phi - spherical coordinates
  ! Uinf - freestream velocity (x-direction only)
  ! vort - local vorticity along phi unit vector
  !----------------------------------------------------------------------------!
  REAL(mk)                 :: tx_center
  REAL(mk)                 :: ty_center
  REAL(mk)                 :: tz_center
  REAL(mk)                 :: tx, ty, tz
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
  REAL(mk), DIMENSION(7)   :: offsetx, offsety, offsetz
  REAL(mk)                 :: Rs,Uinf
  REAL(mk)                 :: R,theta,phi,invR
  REAL(mk)                 :: vort
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: i,j,k,o,pp,isub,isubl,iv
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: mg stuff
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! TODO!??!!: - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_init_sphere_stokes',msg,info)
     STOP ! brutal
  END IF

  isubl = isublist(1)
!  ilo = istart(:,isubl) !&
!  ihi = ilo + ndata(:,isubl)-1 !&

  !-----------------------------------------------------
  !                          Restart from netcdf file
  !-----------------------------------------------------
  IF(netcdf_restart) THEN
     time  = netcdf_time
     itime = netcdf_itime
     dt    = netcdf_dt
     CALL wvic_netcdf2field(info)
     IF(rank.EQ.0) THEN
        WRITE(msg,*) 'restarted from netcdf, itime = ',netcdf_itime
        CALL ppm_write(rank,'wvic_init_sphere',msg,info)
     END IF
     GOTO 1122
  END IF

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! ty/tx/tz_center is the center of the helix
  !----------------------------------------------------------------------------!
  length    = max_physg - min_physg 
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))
  tz_center = 0.5_mk*(min_physg(3) + max_physg(3))
  Rs = sphere_radius
  Uinf = u_infty(3) !!ER DET DEN RIGTIGE RETNING?
  !----------------------------------------------------------------------------!
  ! determining offset for neighbouring boxes (periodic bc's)
  !----------------------------------------------------------------------------!
  offsetx(1)=0_mk
  offsety(1)=0_mk
  offsetz(1)=0_mk
  offsetx(2)=length(1)
  offsety(2)=0_mk
  offsetz(2)=0_mk
  offsetx(3)=-length(1)
  offsety(3)=0_mk
  offsetz(3)=0_mk
  offsetx(4)=0_mk
  offsety(4)=length(2)
  offsetz(4)=0_mk
  offsetx(5)=0_mk
  offsety(5)=-length(2)
  offsetz(5)=0_mk
  offsetx(6)=0_mk
  offsety(6)=0_mk
  offsetz(6)=length(3)
  offsetx(7)=0_mk
  offsety(7)=0_mk
  offsetz(7)=-length(3)

  !----------------------------------------------------------------------------!
  ! calculating predefined, multiply used, variables for each helix
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! create vortices
  !----------------------------------------------------------------------------!
  max_vorticity = Uinf / Rs
  !----------------------------------------------------------------------------!
  ! For each subdomain assigned to the processor...
  ! isubl is the index of the currently treated subdomain
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
    isubl = isublist(isub)
    !-------------------------------------------------------------------------!
    ! loop through the domain. outer loop is the z-layers
    ! i,j,k are counters for x, y, z
    !-------------------------------------------------------------------------!
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          !-------------------------------------------------------------------!
          ! tx, ty, tz are the coordinates (closest to O) of the cell
          ! tx/ty/tz_center is subtracted since the helix revolves around (0,0)
          !-------------------------------------------------------------------!
          twp = 0.0_mk
          DO o=1,1
            tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center + offsetx(o)
            ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center + offsety(o)
            tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center + offsetz(o)

            r = sqrt(tx**2 + ty**2 + tz**2)

            IF (r .GT. Rs*0.5_mk) THEN 
              theta = ATAN2(ty,tx) 
              phi = ATAN2(sqrt(tx**2 + ty**2),tz) !ACOS(tz/r) =
              vort = -1.5_mk * Uinf*SIN(phi)*Rs/r**2
            ELSE
              vort = 0.0_mk
              theta = 0.0_mk
              phi = 0.0_mk
            END IF

            twp(1) =  twp(1) - sin(theta)*vort
            twp(2) =  twp(2) + cos(theta)*vort
            twp(3) =  twp(3)

          END DO !o

          !-----------------------------------------------------
          !  SET FIELD VORTICITY VALUES
          !-----------------------------------------------------
          field_wp(1,i,j,k,isub) = twp(1)
          field_wp(2,i,j,k,isub) = twp(2)
          field_wp(3,i,j,k,isub) = twp(3)
        END DO !i
      END DO !j
    END DO !k
  END DO !isub

  !---- Vorticity has been initialized on the field.
  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)
  CALL ppm_fft_solenoidal(field_wp,mesh_id,topo_id,t_topoid,fmesh_id, & 
     & ghostsize, info)



  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)

  CALL wvic_imposezero_sphere

1122 CONTINUE

  !----------------------------------------------------------------------------!
  ! create particles with style
  !----------------------------------------------------------------------------!
  CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,topo_id,mesh_id,&
       & (/cutoff,HUGE(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
  WRITE(msg,*) ' created ',np,' particles'
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_sphere_stokes',msg,info)
  !----------------------------------------------------------------------------!
  ! all set
  !----------------------------------------------------------------------------!

END SUBROUTINE wvic_init_sphere_stokes


!==============================================================================!
!------------------------------------------------------------------------------!
!  initialize vorticity field for onset flow
!------------------------------------------------------------------------------!
SUBROUTINE wvic_init_onset

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft

  !----------------------------------------------------------------------------!
  ! interfaces
  INTERFACE
     SUBROUTINE wvic_alloc_field_s (vfield_up, info)
       USE module_wvic
       REAL (mk), DIMENSION (:, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_s
     SUBROUTINE wvic_alloc_field(vfield_up, ilda, info)
       USE module_wvic
       REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: vfield_up
       INTEGER                   , INTENT(out) :: info
       INTEGER                   , INTENT(in ) :: ilda
     END SUBROUTINE wvic_alloc_field
  END INTERFACE
  !----------------------------------------------------------------------------!
  REAL(MK), EXTERNAL :: rtbis, rtflsp, rtsec, zriddr, zbrent, sisec, sisec2

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! R, theta, phi - spherical coordinates
  ! Uinf - freestream velocity (x-direction only)
  ! vort - local vorticity along phi unit vector
  !----------------------------------------------------------------------------!
!  REAL(mk)                 :: tx_center
!  REAL(mk)                 :: ty_center
!  REAL(mk)                 :: tz_center
!  REAL(mk)                 :: tx, ty, tz
!  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
!  REAL(mk), DIMENSION(8)   :: offsetx, offsety, offsetz
!  REAL(mk)                 :: Uinf
!  REAL(mk)                 :: R,theta,phi
!  REAL(mk)                 :: vort
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: i,j,k,o,pp,isub,isubl,iv
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: mg stuff
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! TODO!??!!: - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_init_onset',msg,info)
     STOP ! brutal
  END IF

  isubl = isublist(1)

  !-----------------------------------------------------
  !                          Restart from netcdf file
  !-----------------------------------------------------
  IF(netcdf_restart) THEN
     time  = netcdf_time
     itime = netcdf_itime
     dt    = netcdf_dt
     CALL wvic_netcdf2field(info)
     IF(rank.EQ.0) THEN
        WRITE(msg,*) 'restarted from netcdf, itime = ',netcdf_itime
        CALL ppm_write(rank,'wvic_init_onset',msg,info)
     END IF
     GOTO 1123
  END IF

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! ty/tx/tz_center is the center of the helix
  !----------------------------------------------------------------------------!
!  length    = max_physg - min_physg !!skal denne bruges overhovedet? 
!  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))
!  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))
!  tz_center = 0.5_mk*(min_physg(3) + max_physg(3))
!  Uinf = u_infty(3)
  !----------------------------------------------------------------------------!
  ! determining offset for neighbouring boxes (periodic bc's)
  !----------------------------------------------------------------------------!
!  offsetx(1)=0_mk
!  offsety(1)=0_mk
!  offsetz(1)=0_mk
!  offsetx(2)=length(1)
!  offsety(2)=0_mk
!  offsetz(2)=0_mk
!  offsetx(3)=-length(1)
!  offsety(3)=0_mk
!  offsetz(3)=0_mk
!  offsetx(4)=0_mk
!  offsety(4)=length(2)
!  offsetz(4)=0_mk
!  offsetx(5)=0_mk
!  offsety(5)=-length(2)
!  offsetz(5)=0_mk
!  offsetx(6)=0_mk
!  offsety(6)=0_mk
!  offsetz(6)=length(3)
!  offsetx(7)=0_mk
!  offsety(7)=0_mk
!  offsetz(7)=-length(3)

  !----------------------------------------------------------------------------!
  ! calculating predefined, multiply used, variables for each helix
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! create vortices
  !----------------------------------------------------------------------------!
  max_vorticity = 0.0E0_mk
  !----------------------------------------------------------------------------!
  ! For each subdomain assigned to the processor...
  ! isubl is the index of the currently treated subdomain
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
    isubl = isublist(isub)
    !-------------------------------------------------------------------------!
    ! loop through the domain. outer loop is the z-layers
    ! i,j,k are counters for x, y, z
    !-------------------------------------------------------------------------!
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          !-------------------------------------------------------------------!
          ! tx, ty, tz are the coordinates (closest to O) of the cell
          ! tx/ty/tz_center is subtracted since the helix revolves around (0,0)
          !-------------------------------------------------------------------!
          twp = 0.0_mk
!            tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center + offsetx(o)
!            ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center + offsety(o)
!            tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center + offsetz(o)

!            R = sqrt(tx**2 + ty**2 + tz**2)
!            theta = ATAN2(ty,tx)
!            phi = ACOS(tz/R)
 

            twp(1) =  0.0_mk
            twp(2) =  0.0_mk
            twp(3) =  0.0_mk

          !-----------------------------------------------------
          !  SET FIELD VORTICITY VALUES
          !-----------------------------------------------------
          field_wp(1,i,j,k,isub) = twp(1)
          field_wp(2,i,j,k,isub) = twp(2)
          field_wp(3,i,j,k,isub) = twp(3)
        END DO !i
      END DO !j
    END DO !k
  END DO !isub

  !---- Vorticity has been initialized on the field.
  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)
  CALL ppm_fft_solenoidal(field_wp,mesh_id,topo_id,t_topoid,fmesh_id, & 
     & ghostsize, info)



  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)

1123 CONTINUE

  !----------------------------------------------------------------------------!
  ! create particles with style
  !----------------------------------------------------------------------------!
  CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,topo_id,mesh_id,&
       & (/cutoff,HUGE(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
  WRITE(msg,*) ' created ',np,' particles'
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_onset',msg,info)
  !----------------------------------------------------------------------------!
  ! all set
  !----------------------------------------------------------------------------!

END SUBROUTINE wvic_init_onset


!==============================================================================!
!------------------------------------------------------------------------------!
!  initialize vorticity field for Poiseulle flow
!------------------------------------------------------------------------------!
SUBROUTINE wvic_init_poiseulle

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft

  !----------------------------------------------------------------------------!
  ! interfaces
  INTERFACE
     SUBROUTINE wvic_alloc_field_s (vfield_up, info)
       USE module_wvic
       REAL (mk), DIMENSION (:, :, :, :), POINTER :: vfield_up
       INTEGER, INTENT (Out) :: info
     END SUBROUTINE wvic_alloc_field_s
     SUBROUTINE wvic_alloc_field(vfield_up, ilda, info)
       USE module_wvic
       REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: vfield_up
       INTEGER                   , INTENT(out) :: info
       INTEGER                   , INTENT(in ) :: ilda
     END SUBROUTINE wvic_alloc_field
  END INTERFACE
  !----------------------------------------------------------------------------!
  REAL(MK), EXTERNAL :: rtbis, rtflsp, rtsec, zriddr, zbrent, sisec, sisec2

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! R, theta, phi - spherical coordinates
  ! Uinf - freestream velocity (x-direction only)
  ! vort - local vorticity along phi unit vector
  !----------------------------------------------------------------------------!
!  REAL(mk)                 :: tx_center
  REAL(mk)                 :: ty_center
!  REAL(mk)                 :: tz_center
!  REAL(mk)                 :: tx, ty, tz
  REAL(mk)                 :: ty, ty2
!  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: i,j,k,o,pp,isub,isubl,iv
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: mg stuff
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! TODO!??!!: - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_init_poiseulle',msg,info)
     STOP ! brutal
  END IF

  isubl = isublist(1)

  !-----------------------------------------------------
  !                          Restart from netcdf file
  !-----------------------------------------------------
  IF(netcdf_restart) THEN
     time  = netcdf_time
     itime = netcdf_itime
     dt    = netcdf_dt
     CALL wvic_netcdf2field(info)
     IF(rank.EQ.0) THEN
        WRITE(msg,*) 'restarted from netcdf, itime = ',netcdf_itime
        CALL ppm_write(rank,'wvic_init_poiseulle',msg,info)
     END IF
     GOTO 1123
  END IF

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! ty/tx/tz_center is the center of the helix
  !----------------------------------------------------------------------------!
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))

  !----------------------------------------------------------------------------!
  ! calculating predefined, multiply used, variables for each helix
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! create vortices
  !----------------------------------------------------------------------------!
  max_vorticity = 0.0E0_mk
  !----------------------------------------------------------------------------!
  ! For each subdomain assigned to the processor...
  ! isubl is the index of the currently treated subdomain
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
    isubl = isublist(isub)
    !-------------------------------------------------------------------------!
    ! loop through the domain. outer loop is the z-layers
    ! i,j,k are counters for x, y, z
    !-------------------------------------------------------------------------!
    DO j=1,ndata(2,isubl)
      !-------------------------------------------------------------------!
      ! tx, ty, tz are the coordinates (closest to O) of the cell
      ! tx/ty/tz_center is subtracted since the helix revolves around (0,0)
      !-------------------------------------------------------------------!
      twp = 0.0_mk
      ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center

      IF (ty .GT. 0.0_mk) THEN
        ty2 = ty - 0.5_mk*(max_physg(2) - min_physg(2)) 
      ELSE
        ty2 = ty + 0.5_mk*(max_physg(2) - min_physg(2)) 
      END IF

      twp(1) = -1.5_mk * (max_physg(2) - min_physg(2)) * u_infty(3) * ty2 &
             & / ((max_physg(2)-min_physg(2))*0.5_mk - thickness05init) ** 3
      twp(2) = 0.0_mk
      twp(3) = 0.0_mk
      DO k=1,ndata(3,isubl)
        DO i=1,ndata(1,isubl)
            !-----------------------------------------------------
            !  SET FIELD VORTICITY VALUES
            !-----------------------------------------------------
            field_wp(1,i,j,k,isub) = twp(1)*(1.0_mk - field_H(i,j,k,isub))
            field_wp(2,i,j,k,isub) = twp(2)
            field_wp(3,i,j,k,isub) = twp(3)
        END DO !i
      END DO !j
    END DO !k
  END DO !isub

  !---- Vorticity has been initialized on the field.
  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)
  CALL ppm_fft_solenoidal(field_wp,mesh_id,topo_id,t_topoid,fmesh_id, & 
     & ghostsize, info)



  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
  maptype = ppm_param_map_init
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id,&
       & ghostsize,maptype,info)
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(field_wp,lda,topo_id,mesh_id, &
       & ghostsize,maptype,info)

1123 CONTINUE

  !----------------------------------------------------------------------------!
  ! create particles with style
  !----------------------------------------------------------------------------!
  CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,topo_id,mesh_id,&
       & (/cutoff,HUGE(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
  WRITE(msg,*) ' created ',np,' particles'
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_poiseulle',msg,info)
  !----------------------------------------------------------------------------!
  ! all set
  !----------------------------------------------------------------------------!



END SUBROUTINE wvic_init_poiseulle


SUBROUTINE wvic_imposezero_sphere (info)

  USE module_wvic
  USE ppm_module_write

  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER, INTENT(out)              :: info

  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER :: isub, i , j , k, isubl
  REAL(mk) :: tim1s, tim1e
  CHARACTER(len=256) :: msg

  INCLUDE 'mpif.h'

!??  tim1s = MPI_WTIME()

  DO isub=1,nsublist
     isubl = isublist(isub)

     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)

              field_wp(1,i,j,k,isub) = field_wp(1,i,j,k,isub) * &
              & (1-field_H(i,j,k,isub))
              field_wp(2,i,j,k,isub) = field_wp(2,i,j,k,isub) * &
              & (1-field_H(i,j,k,isub))
              field_wp(3,i,j,k,isub) = field_wp(3,i,j,k,isub) * &
              & (1-field_H(i,j,k,isub))
!              field_up(1,i,j,k,isub) = field_up(1,i,j,k,isub) * &
!              & (1-field_H(i,j,k,isub))
!              field_up(2,i,j,k,isub) = field_up(2,i,j,k,isub) * &
!              & (1-field_H(i,j,k,isub))
!              field_up(3,i,j,k,isub) = field_up(3,i,j,k,isub) * &
!              & (1-field_H(i,j,k,isub))

           END DO

        END DO

     END DO

  END DO

!??  tim1e = MPI_WTIME()
!??  if(verbose) WRITE(msg,*) 'took ',1000.0*(tim1e-tim1s)
!??  if(verbose) CALL ppm_write(rank,'wvic_dgammadt',msg,info)

END SUBROUTINE wvic_imposezero_sphere

