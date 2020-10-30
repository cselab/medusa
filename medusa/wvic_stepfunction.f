!-------------------------------------------------------------------------------
! WVIC_INIT_STDOBJ
! 2007/2008
! Initialises mollified step function \tilde{\chi}_S analytically for basic
! geometries
!
! Johannes Tophoej Rasmussen betonarbejder@gmail.com
!-------------------------------------------------------------------------------

SUBROUTINE wvic_init_stepfunc
  USE module_wvic
  USE ppm_module_rmsh_create_part

  REAL(mk) :: pi
  pi = ACOS(-1.0_mk)

  CALL wvic_solid_velocity(info)

  If ((SUM(ABS(harmonic_amplitude)) .NE. 0.0_mk) &
  & .AND. (SUM(ABS(harmonic_period)) .NE. 0.0_mk)) THEN
    object_move = .true.
    IF (rank .EQ. 0) THEN
      WRITE(UNIT=0,*) 'Harmonic oscillation activated'
    ENDIF
  ELSEIf ((SUM(u_solid**2)) .EQ. 0.0_mk) THEN
    object_move = .false.
    IF (rank .EQ. 0) THEN
      WRITE(UNIT=0,*) 'Static object, particle step function disabled'
    ENDIF
  ELSE
    object_move = .true.
  ENDIF


  SELECT CASE(flow_case)
     CASE(9)
       !-----------------------------------------------------
       ! Poiseulle flow between two parallel plates
       !-----------------------------------------------------
     CASE(10)
       !-----------------------------------------------------
       ! Onset flow past sphere
       !-----------------------------------------------------
     CASE(11)
       !-----------------------------------------------------
       ! Medusa
       !-----------------------------------------------------
!       CALL wvic_stepfunc_medusa !REMOVE
     CASE(12)
       !-----------------------------------------------------
       ! Cylinder array - Dong
       !-----------------------------------------------------
     CASE(13)
       !-----------------------------------------------------
       ! Harmonic oscilations of plate
       !-----------------------------------------------------
       IF (harmonic_amplitude(1) .NE. 0.0_mk) THEN
         object_offset(1) = object_offset(1) + harmonic_amplitude(1)* &
                      & SIN(2*pi/harmonic_period(1)*time+2*pi*harmonic_phase(1))
         u_solid(1) = harmonic_amplitude(1)*2*pi/harmonic_period(1) * &
               & COS(2*pi/harmonic_period(1)*time+2*pi*harmonic_phase(1))
       ENDIF 
       IF (harmonic_amplitude(2) .NE. 0.0_mk) THEN
         object_offset(2) = object_offset(2) + harmonic_amplitude(2)* & 
                      & SIN(2*pi/harmonic_period(2)*time+2*pi*harmonic_phase(2))
         u_solid(2) = harmonic_amplitude(2)*2*pi/harmonic_period(2) * &
               & COS(2*pi/harmonic_period(2)*time+2*pi*harmonic_phase(2))
       ENDIF 
       IF (harmonic_amplitude(3) .NE. 0.0_mk) THEN
         object_offset(3) = object_offset(3) + harmonic_amplitude(3)* &
                      & SIN(2*pi/harmonic_period(3)*time+2*pi*harmonic_phase(3))
         u_solid(3) = harmonic_amplitude(3)*2*pi/harmonic_period(3) * &
               & COS(2*pi/harmonic_period(3)*time+2*pi*harmonic_phase(3))
       ENDIF 
     CASE(14)
       !-----------------------------------------------------
       ! STL
       !-----------------------------------------------------
  END SELECT


!JTR mass and center of mass could also be calculated here

END SUBROUTINE wvic_init_stepfunc


SUBROUTINE wvic_stepfunc_sphere

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft
  IMPLICIT NONE

  REAL(MK), EXTERNAL :: stepfunction1
  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! tx/ty/tz_center = coordinates of the center of the sphere
  ! tx/ty/tz = coordinates of treated cell
  ! length = size of the domain
  ! !for later use:! offsetx/y/z = location of periodic/adjacent domains
  ! Rs = radius of sphere
  ! level = levelset/epsilon
  ! epsilon = 1/epsilon (epsilon ~ dx)
  !----------------------------------------------------------------------------!
  REAL(mk)                 :: tx_center
  REAL(mk)                 :: ty_center
  REAL(mk)                 :: tz_center
  REAL(mk)                 :: tx, ty, tz
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk)                 :: Rs, theta, phi, level, epsilon
  REAL(mk)                 :: tmp, move,tmp2
  CHARACTER(len=32)        :: filename
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  ! i,j,k = counter for spatial loop
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  INTEGER                  :: maptype, info
  INTEGER                  :: i,j,k,isub,isubl
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff ??
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! tx/ty/tz_center is the center of the sphere
  ! ra = 1/(2 sigma^2)
  !----------------------------------------------------------------------------!
  length    = max_physg - min_physg !!skal denne bruges overhovedet?
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))+object_offset(1)
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))+object_offset(2)
  tz_center = 0.5_mk*(min_physg(3) + max_physg(3))+object_offset(3)
  epsilon = 1.0_mk / (sqrt(dx**2 + dy**2 + dz**2)*stepfunction_band)
  Rs = sphere_radius
  !----------------------------------------------------------------------------!
  ! For each subdomain assigned to the processor...
  ! isubl is the index of the currently treated subdomain
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
    isubl = isublist(isub)
    !-------------------------------------------------------------------------!
    ! do for each Z-layer of the subdomain...
    ! i,j,k are counters for x, y, z
    !-------------------------------------------------------------------------!
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          !-------------------------------------------------------------------!
          ! The j and i loops covers all cells in the z-layer
          ! tx, ty, tz are the coordinates in a coordinate system with origo in
          ! the center of the sphere (tx/ty/tz_center)
          !-------------------------------------------------------------------!
          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center
          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center
          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center

          level = (sqrt(tx**2 + ty**2 + tz**2) - Rs)
          IF (step_function .EQ. 0) THEN
!            theta = ATAN2(ty,tx)
!            phi = ATAN2(sqrt(tx**2 + ty**2),tz)

            field_H(i,j,k,isub) = -0.5_mk * TANH(level*epsilon) + 0.5_mk
          ELSEIF (step_function .EQ. 1) THEN
            field_H(i,j,k,isub) = stepfunction1(level*epsilon)
          END IF 
        END DO !i
      END DO !j
    END DO !k
  END DO !isub

!JTR mass and center of mass should be calculated here - or moved to just after this routine - to insert the calculation more generally

  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
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

END SUBROUTINE wvic_stepfunc_sphere



SUBROUTINE wvic_stepfunc_plate

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  CHARACTER(len=256) :: msg
  INTEGER, PARAMETER :: md = kind(2.0d0)

  INTEGER                    :: i,tr
  REAL(mk), DIMENSION(3)     :: vertex1,vertex2,vertex3,vecu,vecv,vecw
  REAL(mk), DIMENSION(12,3,3):: blockdata
  REAL(mk), DIMENSION(3)     :: vert1,vert2,vert3,vert4,vert5,vert6,vert7,vert8
  REAL(mk), DIMENSION(3)     :: t_center
  INCLUDE 'mpif.h'


  t_center = 0.5_mk*(max_physg + min_physg) + object_offset 

  tri_count = 12
  !-----------------------------------------------------------------------------
  ! Allocate vector dotproduct arrays
  !-----------------------------------------------------------------------------
  ALLOCATE(tri_norm(tri_count,3))
  ALLOCATE(tri_base(tri_count,3))
  ALLOCATE(tri_vecu(tri_count,3))
  ALLOCATE(tri_vecv(tri_count,3))
  ALLOCATE(tri_vecw(tri_count,3))
  ALLOCATE(tri_denom(tri_count))
  ALLOCATE(tri_udotv(tri_count))
  ALLOCATE(tri_udotu(tri_count))
  ALLOCATE(tri_vdotv(tri_count))
  ALLOCATE(tri_wdotw(tri_count))

  vert1= (/-0.5_mk*block_wx,-0.5_mk*block_wy,-0.5_mk*block_wz/)
  vert2= (/-0.5_mk*block_wx,-0.5_mk*block_wy, 0.5_mk*block_wz/)
  vert3= (/ 0.5_mk*block_wx,-0.5_mk*block_wy,-0.5_mk*block_wz/)
  vert4= (/ 0.5_mk*block_wx,-0.5_mk*block_wy, 0.5_mk*block_wz/)
  vert5= (/-0.5_mk*block_wx, 0.5_mk*block_wy,-0.5_mk*block_wz/)
  vert6= (/-0.5_mk*block_wx, 0.5_mk*block_wy, 0.5_mk*block_wz/)
  vert7= (/ 0.5_mk*block_wx, 0.5_mk*block_wy,-0.5_mk*block_wz/)
  vert8= (/ 0.5_mk*block_wx, 0.5_mk*block_wy, 0.5_mk*block_wz/)

  vert1 = vert1 + t_center 
  vert2 = vert2 + t_center 
  vert3 = vert3 + t_center 
  vert4 = vert4 + t_center 
  vert5 = vert5 + t_center 
  vert6 = vert6 + t_center 
  vert7 = vert7 + t_center 
  vert8 = vert8 + t_center 

  tri_norm(1,:)= (/-1.0_mk,0.0_mk,0.0_mk/)
  tri_norm(2,:)= (/-1.0_mk,0.0_mk,0.0_mk/)
  tri_norm(3,:)= (/0.0_mk,0.0_mk,1.0_mk/)
  tri_norm(4,:)= (/0.0_mk,0.0_mk,1.0_mk/)
  tri_norm(5,:)= (/1.0_mk,0.0_mk,0.0_mk/)
  tri_norm(6,:)= (/1.0_mk,0.0_mk,0.0_mk/)
  tri_norm(7,:)= (/0.0_mk,0.0_mk,-1.0_mk/)
  tri_norm(8,:)= (/0.0_mk,0.0_mk,-1.0_mk/)
  tri_norm(9,:)= (/0.0_mk,1.0_mk,0.0_mk/)
  tri_norm(10,:)= (/0.0_mk,1.0_mk,0.0_mk/)
  tri_norm(11,:)= (/0.0_mk,-1.0_mk,0.0_mk/)
  tri_norm(12,:)= (/0.0_mk,-1.0_mk,0.0_mk/)

  blockdata(1,1,:) = vert1
  blockdata(1,2,:) = vert2
  blockdata(1,3,:) = vert6

  blockdata(2,1,:) = vert1
  blockdata(2,2,:) = vert6
  blockdata(2,3,:) = vert5

  blockdata(3,1,:) = vert2
  blockdata(3,2,:) = vert4
  blockdata(3,3,:) = vert8

  blockdata(4,1,:) = vert2
  blockdata(4,2,:) = vert8
  blockdata(4,3,:) = vert6

  blockdata(5,1,:) = vert4
  blockdata(5,2,:) = vert3
  blockdata(5,3,:) = vert7

  blockdata(6,1,:) = vert4
  blockdata(6,2,:) = vert7
  blockdata(6,3,:) = vert8

  blockdata(7,1,:) = vert3
  blockdata(7,2,:) = vert1
  blockdata(7,3,:) = vert5

  blockdata(8,1,:) = vert3
  blockdata(8,2,:) = vert5
  blockdata(8,3,:) = vert7

  blockdata(9,1,:) = vert5
  blockdata(9,2,:) = vert6
  blockdata(9,3,:) = vert8

  blockdata(10,1,:) = vert5
  blockdata(10,2,:) = vert8
  blockdata(10,3,:) = vert7

  blockdata(11,1,:) = vert1
  blockdata(11,2,:) = vert2
  blockdata(11,3,:) = vert4

  blockdata(12,1,:) = vert1
  blockdata(12,2,:) = vert4
  blockdata(12,3,:) = vert3

  DO tr=1,tri_count
    vertex1 = blockdata(tr,1,:)
    vertex2 = blockdata(tr,2,:)
    vertex3 = blockdata(tr,3,:)
    tri_base(tr,:) = vertex1
    vecu = vertex2 - vertex1
    vecv = vertex3 - vertex1
    vecw = vertex3 - vertex2
    tri_vecu(tr,:) = vecu
    tri_vecv(tr,:) = vecv
    tri_vecw(tr,:) = vecw
    tri_udotu(tr) = (vecu(1)*vecu(1)+vecu(2)*vecu(2)+vecu(3)*vecu(3))
    tri_vdotv(tr) = (vecv(1)*vecv(1)+vecv(2)*vecv(2)+vecv(3)*vecv(3))
    tri_wdotw(tr) = (vecw(1)*vecw(1)+vecw(2)*vecw(2)+vecw(3)*vecw(3))
    tri_udotv(tr) = (vecu(1)*vecv(1)+vecu(2)*vecv(2)+vecu(3)*vecv(3))
    tri_denom(tr)  = 1.0_mk/(tri_udotu(tr)*tri_vdotv(tr)-tri_udotv(tr)**2)
  END DO    

  !BOUNDING BOX	!
  bndminx = -0.5_mk*block_wx + t_center(1)
  bndmaxx =  0.5_mk*block_wx + t_center(1)
  bndminy = -0.5_mk*block_wy + t_center(2)
  bndmaxy =  0.5_mk*block_wy + t_center(2)
  bndminz = -0.5_mk*block_wz + t_center(3)
  bndmaxz =  0.5_mk*block_wz + t_center(3)

END SUBROUTINE wvic_stepfunc_plate



SUBROUTINE wvic_stepfunc_cylinderarray

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft
  IMPLICIT NONE

  REAL(MK), EXTERNAL :: stepfunction1
  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! tx/ty/tz_center = coordinates of the center of the sphere
  ! tx/ty/tz = coordinates of treated cell
  ! length = size of the domain
  ! !for later use:! offsetx/y/z = location of periodic/adjacent domains
  ! Rc = radius of cylinder
  ! level = levelset/epsilon
  ! epsilon = 1/epsilon (epsilon ~ dx)
  !----------------------------------------------------------------------------!
  REAL(mk)                           :: tx_center
  REAL(mk)                           :: ty_center
  REAL(mk)                           :: tz_center
  REAL(mk)                           :: tx, ty, tz
  REAL(mk), DIMENSION(3)             :: length
  REAL(mk), DIMENSION(8)             :: offsetx, offsety, offsetz
  REAL(mk)                           :: Rc, theta, phi, level, epsilon
  REAL(mk)                           :: tmp, move
  REAL(mk)                           :: tmpH
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  ! i,j,k = counter for spacial loop
  ! o loops adjacent domains due to periodic BC's
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: i,j,k,o,no,isub,isubl
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff ??
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! tx/ty/tz_center is the center of the sphere
  ! ra = 1/(2 sigma^2)
  !----------------------------------------------------------------------------!
  length    = max_physg - min_physg
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))+object_offset(1)
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))+object_offset(2)
  tz_center = 0.5_mk*(min_physg(3) + max_physg(3))+object_offset(3)
  epsilon = 1.0_mk / (sqrt(dx**2 + dy**2 + dz**2)*stepfunction_band)
  Rc = cylinder_radius
  IF (cylinder_no .EQ. 1) THEN
    no=1
    offsetx(1)= 0.0_mk
    offsety(1)= 0.0_mk
    offsetz(1)= 0.0_mk
  ELSEIF (cylinder_no .EQ. 2) THEN
    no=8
    offsetx(1)= 0.0_mk
    offsety(1)= 0.25*length(2)
    offsetz(1)= -0.25*length(3)
    offsetx(2)= 0.0_mk
    offsety(2)= -0.25*length(2)
    offsetz(2)= 0.25*length(3)
    offsetx(3)= 0.0_mk
    offsety(3)= 0.75*length(2)
    offsetz(3)= -0.75*length(3)
    offsetx(4)= 0.0_mk
    offsety(4)= 0.75*length(2)
    offsetz(4)= 0.25*length(3)
    offsetx(5)= 0.0_mk
    offsety(5)= 0.25*length(2)
    offsetz(5)= 0.75*length(3)
    offsetx(6)= 0.0_mk
    offsety(6)= -0.75*length(2)
    offsetz(6)= 0.75*length(3)
    offsetx(7)= 0.0_mk
    offsety(7)= -0.75*length(2)
    offsetz(7)= -0.25*length(3)
    offsetx(8)= 0.0_mk
    offsety(8)= -0.25*length(2)
    offsetz(8)= -0.75*length(3)
  END IF
  !----------------------------------------------------------------------------!
  ! For each subdomain assigned to the processor...
  ! isubl is the index of the currently treated subdomain
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
    isubl = isublist(isub)
    !-------------------------------------------------------------------------!
    ! do for each Z-layer of the subdomain...
    ! i,j,k are counters for x, y, z
    !-------------------------------------------------------------------------!
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          !-------------------------------------------------------------------!
          ! The j and i loops covers all cells in the z-layer
          ! tx, ty, tz are the coordinates in a coordinate system with origo in
          ! the center of the sphere (tx/ty/tz_center)
          !-------------------------------------------------------------------!
          tmpH = 0.0_mk
          DO o=1,no
            tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center + offsetx(o)
            ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center + offsety(o)
            tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center + offsetz(o)

            level = (sqrt(ty**2 + tz**2) - Rc)
!            theta = ATAN2(ty,tz)
            IF (step_function .EQ. 0) THEN
              tmpH = tmpH + ( -0.5_mk * TANH(level*epsilon) + 0.5_mk )
            ELSEIF (step_function .EQ. 1) THEN
              tmpH = tmpH + stepfunction1(level*epsilon)
            END IF 
          END DO !o
          field_H(i,j,k,isub) = tmpH 
        END DO !i
      END DO !j
    END DO !k
  END DO !isub

!JTR mass and center of mass should be calculated here - or moved to just after this routine - to insert the calculation more generally

!  tmp=0.0_mk
!  DO isub=1,nsublist
!    isubl = isublist(isub)
!    DO k=1,ndata(3,isubl)
!      DO j=1,ndata(2,isubl)
!        DO i=1,ndata(1,isubl)
!          tmp=tmp+field_H(i,j,k,isub)
!        END DO !i
!      END DO !j
!    END DO !k
!  END DO !isub
!WRITE(msg,*) 'JTR: rank', rank, tmp,'\n'
!WRITE(UNIT=0,*) msg
!stop

  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
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

END SUBROUTINE wvic_stepfunc_cylinderarray



SUBROUTINE wvic_stepfunc_complex

!------------------------------------------------------------------------------
! OBSOLETE & ABANDONED. TESTS DEVIATION FROM ANALYTIC GRADIENT FIELD TO FD
!------------------------------------------------------------------------------

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft
  IMPLICIT NONE

  INCLUDE 'mpif.h' !JTR - fjern naar mpi ikke er noedvendig

  INTEGER, PARAMETER :: md = kind(2.0d0)

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! tx/ty/tz_center = coordinates of the center of the sphere
  ! tx/ty/tz = coordinates of treated cell
  ! length = size of the domain
  ! !for later use:! offsetx/y/z = location of periodic/adjacent domains
  ! Rs = radius of sphere
  ! level = levelset/epsilon
  ! epsilon = 1/epsilon (epsilon ~ dx)
  !----------------------------------------------------------------------------!
  REAL(mk)                 :: tx_center
  REAL(mk)                 :: ty_center
  REAL(mk)                 :: tz_center
  REAL(mk)                 :: tx, ty, tz
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk)                 :: Pg, Rg, dR, theta, phi, Rs, level, epsilon
  REAL(mk)                 :: tmp, move
  REAL(mk)                 :: dvnc, gdvnc
  REAL(mk)                 :: lobjectvolume,gobjectvolume,dv !JTR slet
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  ! i,j,k = counter for spacial loop
  ! o loops adjacent domains due to periodic BC's
  !----------------------------------------------------------------------------!
  CHARACTER(len=8)         :: file
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: i,j,k,o,isub,isubl
  REAL(mk)                 :: fac1,fac2,fac3,fac4
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! tx/ty/tz_center is the center of the sphere
  ! ra = 1/(2 sigma^2)
  !----------------------------------------------------------------------------!
  length    = max_physg - min_physg !!skal denne bruges overhovedet?
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))
  tz_center = 0.5_mk*(min_physg(3) + max_physg(3))
  epsilon = 1.0_mk / (sqrt(dx**2 + dy**2 + dz**2)*stepfunction_band)
  Rs = sphere_radius
  Rg = 1.0_mk
  Pg = 6.0_mk
  dR = 0.3_mk
  dvnc = 0.0_mk
  fac1=1.0_mk/(dx*840.0_mk)
  fac2=1.0_mk/(dy*840.0_mk)
  fac3=1.0_mk/(dz*840.0_mk)
  fac4=840.0_mk/820.0_mk
  fac1=1.0_mk/(dx*2.0_mk)
  fac2=1.0_mk/(dy*2.0_mk)
  fac3=1.0_mk/(dz*2.0_mk)
  fac1=1.0_mk/(dx*12.0_mk)
  fac2=1.0_mk/(dy*12.0_mk)
  fac3=1.0_mk/(dz*12.0_mk)
  fac1=1.0_mk/(dx*60.0_mk)
  fac2=1.0_mk/(dy*60.0_mk)
  fac3=1.0_mk/(dz*60.0_mk)

  DO isub=1,nsublist
  !----------------------------------------------------------------------------!
  ! initialize heaviside - non-smooth
  !----------------------------------------------------------------------------!
    isubl = isublist(isub)
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          !-------------------------------------------------------------------!
          ! The j and i loops covers all cells in the z-layer
          ! tx, ty, tz are the coordinates in a coordinate system with origo in
          ! the center of the sphere (tx/ty/tz_center)
          !-------------------------------------------------------------------!
          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center
          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center
          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center

          !note sign on level. this may alter ubar
          level = (sqrt(tx**2 + ty**2 + tz**2) - Rs)
          tmp = -0.5_mk * (1.0_mk - TANH(level*epsilon)**2)*epsilon

          theta = ATAN2(ty,tx)
          phi = ATAN2(sqrt(tx**2 + ty**2),tz)

          field_H(i,j,k,isub) = -0.5_mk * TANH(level*epsilon) + 0.5_mk
          field_ubar(1,i,j,k,isub) = COS(theta)*SIN(phi)*tmp
          field_ubar(2,i,j,k,isub) = SIN(theta)*SIN(phi)*tmp
          field_ubar(3,i,j,k,isub) =            COS(phi)*tmp
        END DO !i
      END DO !j
    END DO !k

#ifdef _johssluk
  !----------------------------------------------------------------------------!
  ! smooth heaviside
  !----------------------------------------------------------------------------!
  DO o=1,11 !the number of time to diffuse the step function
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
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          field_ubar(1,i,j,k,isub) = 0.0_mk
          field_ubar(1,i,j,k,isub) = field_ubar(1,i,j,k,isub) + 0.28_mk * &
                              & field_H(i,j,k,isub)
          field_ubar(1,i,j,k,isub) = field_ubar(1,i,j,k,isub) + 0.12_mk * &
                              & field_H(i+1,j,k,isub)
          field_ubar(1,i,j,k,isub) = field_ubar(1,i,j,k,isub) + 0.12_mk * &
                              & field_H(i-1,j,k,isub)
          field_ubar(1,i,j,k,isub) = field_ubar(1,i,j,k,isub) + 0.12_mk * &
                              & field_H(i,j+1,k,isub)
          field_ubar(1,i,j,k,isub) = field_ubar(1,i,j,k,isub) + 0.12_mk * &
                              & field_H(i,j-1,k,isub)
          field_ubar(1,i,j,k,isub) = field_ubar(1,i,j,k,isub) + 0.12_mk * &
                              & field_H(i,j,k+1,isub)
          field_ubar(1,i,j,k,isub) = field_ubar(1,i,j,k,isub) + 0.12_mk * &
                              & field_H(i,j,k-1,isub)

          !field_H(i,j,k,isub) = 0.0_mk
          !field_ubar(1,i,j,k,isub) = COS(theta)*SIN(phi)*tmp
          !field_ubar(2,i,j,k,isub) = SIN(theta)*SIN(phi)*tmp
          !field_ubar(3,i,j,k,isub) =            COS(phi)*tmp
        END DO !i
      END DO !j
    END DO !k
    !--------------------------------------------------------------------------!
    ! save smoothened step function to vorticity field(1) to allow ghosting
    !--------------------------------------------------------------------------!
    DO k=1,ndata(3,isubl)-1
      DO j=1,ndata(2,isubl)-1
        DO i=1,ndata(1,isubl)-1
          field_H(i,j,k,isub) = field_ubar(1,i,j,k,isub)
        END DO !i
      END DO !j
    END DO !k
  END DO !o
  !----------------------------------------------------------------------------!
  ! calculate gradient of step function field
  !----------------------------------------------------------------------------!
#endif
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
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
#ifdef _johsstopST3
          dvnc = dvnc + sqrt( &
               & ( field_ubar(1,i,j,k,isub) - (field_H(i+1,j,k,isub)- &
               & field_H(i-1,j,k,isub))*0.5_mk/dx )**2 + &
               & ( field_ubar(2,i,j,k,isub) - (field_H(i,j+1,k,isub)- &
               & field_H(i,j-1,k,isub))*0.5_mk/dy )**2 + &
               & ( field_ubar(3,i,j,k,isub) - (field_H(i,j,k+1,isub)- &
               & field_H(i,j,k-1,isub))*0.5_mk/dz )**2 )
#endif
#ifdef _johsstopST5
          dvnc = dvnc + sqrt( &
               & ( field_ubar(1,i,j,k,isub) - &
               & (-field_H(i+2,j,k,isub) + 8.0*field_H(i+1,j,k,isub) &
               &  -8.0*field_H(i-1,j,k,isub) + field_H(i-2,j,k,isub) &
               &  )/(dx*12.0) )**2 + &
               & ( field_ubar(2,i,j,k,isub) - &
               & (-field_H(i,j+2,k,isub) + 8.0*field_H(i,j+1,k,isub) &
               &  -8.0*field_H(i,j-1,k,isub) + field_H(i,j-2,k,isub) &
               &  )/(dy*12.0) )**2 + &
               & ( field_ubar(3,i,j,k,isub) - &
               & (-field_H(i,j,k+2,isub) + 8.0*field_H(i,j,k+1,isub) &
               &  -8.0*field_H(i,j,k-1,isub) + field_H(i,j,k-2,isub) &
               &  )/(dz*12.0) )**2)
#endif
#ifdef _johsstopST7
          dvnc = dvnc + sqrt( &
             & (field_ubar(1,i,j,k,isub) - fac1 * &
             & (         field_H(i+3,j,k,isub) -  9.0_mk*field_H(i+2,j,k,isub) &
             & + 45.0_mk*field_H(i+1,j,k,isub) - 45.0_mk*field_H(i-1,j,k,isub) &
         & +  9.0_mk*field_H(i-2,j,k,isub) -         field_H(i-3,j,k,isub)))**2 &
           & + (field_ubar(2,i,j,k,isub) - fac2 * &
             & (         field_H(i,j+3,k,isub) -  9.0_mk*field_H(i,j+2,k,isub) &
             & + 45.0_mk*field_H(i,j+1,k,isub) - 45.0_mk*field_H(i,j-1,k,isub) &
             & +  9.0_mk*field_H(i,j-2,k,isub) -     field_H(i,j-3,k,isub)))**2 &
           & + (field_ubar(3,i,j,k,isub) - fac3 * &
             & (         field_H(i,j,k+3,isub) -  9.0_mk*field_H(i,j,k+2,isub) &
             & + 45.0_mk*field_H(i,j,k+1,isub) - 45.0_mk*field_H(i,j,k-1,isub) &
         & +  9.0_mk*field_H(i,j,k-2,isub) -         field_H(i,j,k-3,isub)))**2 )
#endif
#ifdef _johsstopST9
          dvnc = dvnc + sqrt( &
     & ( field_ubar(1,i,j,k,isub) - fac1 * &
             & (  - fac4*field_H(i+4,j,k,isub) + 32.0_mk*field_H(i+3,j,k,isub) &
             & -168.0_mk*field_H(i+2,j,k,isub) +672.0_mk*field_H(i+1,j,k,isub) &
             & -672.0_mk*field_H(i-1,j,k,isub) +168.0_mk*field_H(i-2,j,k,isub) &
             & - 32.0_mk*field_H(i-3,j,k,isub) +    fac4*field_H(i-4,j,k,isub)))**2 + &
     & ( field_ubar(2,i,j,k,isub) - fac2 * &
             & (  - fac4*field_H(i,j+4,k,isub) + 32.0_mk*field_H(i,j+3,k,isub) &
             & -168.0_mk*field_H(i,j+2,k,isub) +672.0_mk*field_H(i,j+1,k,isub) &
             & -672.0_mk*field_H(i,j-1,k,isub) +168.0_mk*field_H(i,j-2,k,isub) &
             & - 32.0_mk*field_H(i,j-3,k,isub) +    fac4*field_H(i,j-4,k,isub)))**2 + &
      & ( field_ubar(3,i,j,k,isub) - fac3 * &
             & (  - fac4*field_H(i,j,k+4,isub) + 32.0_mk*field_H(i,j,k+3,isub) &
             & -168.0_mk*field_H(i,j,k+2,isub) +672.0_mk*field_H(i,j,k+1,isub) &
             & -672.0_mk*field_H(i,j,k-1,isub) +168.0_mk*field_H(i,j,k-2,isub) &
             & - 32.0_mk*field_H(i,j,k-3,isub) +    fac4*field_H(i,j,k-4,isub)))**2 )
#endif
#ifdef _johsstopTOTAL
!          dvnc = dvnc + sqrt( &
!               & ( field_ubar(1,i,j,k,isub) )**2 + &
!               & ( field_ubar(2,i,j,k,isub) )**2 + &
!               & ( field_ubar(3,i,j,k,isub) )**2 )
#endif

#ifndef _johsstopST7
        field_ubar(1,i,j,k,isub) = fac1 * &
             & (         field_H(i+3,j,k,isub) -  9.0_mk*field_H(i+2,j,k,isub) &
             & + 45.0_mk*field_H(i+1,j,k,isub) - 45.0_mk*field_H(i-1,j,k,isub) &
             & +  9.0_mk*field_H(i-2,j,k,isub) -         field_H(i-3,j,k,isub))
        field_ubar(2,i,j,k,isub) = fac2 * &
             & (         field_H(i,j+3,k,isub) -  9.0_mk*field_H(i,j+2,k,isub) &
             & + 45.0_mk*field_H(i,j+1,k,isub) - 45.0_mk*field_H(i,j-1,k,isub) &
             & +  9.0_mk*field_H(i,j-2,k,isub) -         field_H(i,j-3,k,isub))
        field_ubar(3,i,j,k,isub) = fac3 * &
             & (         field_H(i,j,k+3,isub) -  9.0_mk*field_H(i,j,k+2,isub) &
             & + 45.0_mk*field_H(i,j,k+1,isub) - 45.0_mk*field_H(i,j,k-1,isub) &
             & +  9.0_mk*field_H(i,j,k-2,isub) -         field_H(i,j,k-3,isub))
#endif
#ifdef _johsstopST5
        field_ubar(1,i,j,k,isub) = fac1 * &
             & (-field_H(i+2,j,k,isub) + 8.0*field_H(i+1,j,k,isub) &
             &  -8.0*field_H(i-1,j,k,isub) + field_H(i-2,j,k,isub))
        field_ubar(2,i,j,k,isub) = fac2 * &
             & (-field_H(i,j+2,k,isub) + 8.0*field_H(i,j+1,k,isub) &
             &  -8.0*field_H(i,j-1,k,isub) + field_H(i,j-2,k,isub))
        field_ubar(3,i,j,k,isub) = fac3 * &
             & (-field_H(i,j,k+2,isub) + 8.0*field_H(i,j,k+1,isub) &
             &  -8.0*field_H(i,j,k-1,isub) + field_H(i,j,k-2,isub))
#endif
#ifdef _johsstopST3
        field_ubar(1,i,j,k,isub) = fac1 * &
             & (field_H(i+1,j,k,isub) - field_H(i-1,j,k,isub))
        field_ubar(2,i,j,k,isub) = fac2 * &
             & (field_H(i,j+1,k,isub) - field_H(i,j-1,k,isub))
        field_ubar(3,i,j,k,isub) = fac3 * &
             & (field_H(i,j,k+1,isub) - field_H(i,j,k-1,isub))
#endif

        END DO !i
      END DO !j
    END DO !k
    dv=dx*dy*dz
    lobjectvolume=0.0_mk
    DO k=1,ndata(3,isubl)-1
      DO j=1,ndata(2,isubl)-1
        DO i=1,ndata(1,isubl)-1 !HER SKAL NYE GRAENSER PAA
          lobjectvolume = lobjectvolume + field_H(i,j,k,isub)*dv
        END DO !i
      END DO !j
    END DO !k
    CALL MPI_Reduce(lobjectvolume,gobjectvolume,1,mpi_prec,MPI_SUM,0,comm,info)
    IF (rank .EQ. 0) THEN
      WRITE(msg,*) '\n Objectvolume:', gobjectvolume,'\n'
      WRITE(*,*) msg
      WRITE(unit=0,*) msg
    END IF
  END DO !isub

  CALL MPI_Allreduce(dvnc,gdvnc,1,mpi_prec,MPI_SUM,comm,info)
  IF(rank.EQ.0) THEN

    WRITE(msg,*) 'aDeviance from analytical heaviside gradient to FD: ', gdvnc
    WRITE(UNIT=0,*) msg

  END IF
  !---- Vorticity has been initialized on the field.
  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)

  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
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
END SUBROUTINE wvic_stepfunc_complex

!---------------------------------------------------------------------------
! Smoothens heaviside step function by convolution
!---------------------------------------------------------------------------
SUBROUTINE wvic_stepfunc_smoothen

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_fdsolver_solve
  USE ppm_module_map
  USE ppm_module_fft
  IMPLICIT NONE

  INTEGER, PARAMETER :: md = kind(2.0d0)

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! tx/ty/tz_center = coordinates of the center of the sphere
  ! tx/ty/tz = coordinates of treated cell
  ! length = size of the domain
  ! !for later use:! offsetx/y/z = location of periodic/adjacent domains
  ! Rs = radius of sphere
  ! level = levelset/epsilon
  ! epsilon = 1/epsilon (epsilon ~ dx)
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: i,j,k,o,isub,isubl
  REAL(mk)                 :: fac1,fac2,fac3,fac4
  REAL(mk)                 :: faca,facb,facc,facd
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

!  faca = 0.3658829878_mk
!  facb = 0.07283746515_mk
!  facc = 0.01449998088_mk
!  facd = 0.002886556324_mk
!  faca = 0.529692587072715_mk !9ny - ikke god
!  facb = 0.025970041213382_mk
!  facc = 0.019613309578477_mk
!  facd = 0.009890931338159_mk
!  faca = 0.540521097167045_mk !9
!  facb = 0.029217651871221_mk
!  facc = 0.009328655041086_mk
!  facd = 0.021528641389074_mk
!  faca = 0.735872025121785_mk !19
!  facb = 0.028398324003127_mk
!  facc = 0.0_mk
!  facd = 0.011717253857432_mk
  faca = 0.707835037010006_mk !17
  facb = 0.030792896870216_mk
  facc = 0.000359821140528_mk
  facd = 0.012886216010295_mk
!  faca = 0.833942167852334_mk !31
!  facb = 0.017094811880438_mk
!  facc = 0.000000000000000_mk
!  facd = 0.007936120108130_mk

  IF (stl_stencil .EQ. 0) THEN
    fac1=1.0_mk/(dx*2.0_mk)
    fac2=1.0_mk/(dy*2.0_mk)
    fac3=1.0_mk/(dz*2.0_mk)
  END IF
  IF (stl_stencil .EQ. 1) THEN
    fac1=1.0_mk/(dx*12.0_mk)
    fac2=1.0_mk/(dy*12.0_mk)
    fac3=1.0_mk/(dz*12.0_mk)
  END IF
  IF (stl_stencil .EQ. 2) THEN
    fac1=1.0_mk/(dx*60.0_mk)
    fac2=1.0_mk/(dy*60.0_mk)
    fac3=1.0_mk/(dz*60.0_mk)
  END IF
  IF (stl_stencil .EQ. 3) THEN
    fac1=1.0_mk/(dx*840.0_mk)
    fac2=1.0_mk/(dy*840.0_mk)
    fac3=1.0_mk/(dz*840.0_mk)
    fac4=840.0_mk/820.0_mk
  END IF

  DO isub=1,nsublist
    isubl = isublist(isub)

    !--------------------------------------------------------------------------!
    ! smooth heaviside by convolution
    !--------------------------------------------------------------------------!
    DO o=1,stl_mollify !the number of time to smoothen the step function

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
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          field_ubar(1,i,j,k,isub) = &
          ! a: inner point
            & faca *   field_H(i,j,k,isub) + &
          ! b: faces (6)
            & facb * ( field_H(i+1,j,k,isub) + &
            &          field_H(i-1,j,k,isub) + &
            &          field_H(i,j+1,k,isub) + &
            &          field_H(i,j-1,k,isub) + &
            &          field_H(i,j,k+1,isub) + &
            &          field_H(i,j,k-1,isub) ) + &
          ! c: edges (12)
            & facc * ( field_H(i+1,j+1,k,isub) + &
            &          field_H(i+1,j-1,k,isub) + &
            &          field_H(i+1,j,k+1,isub) + &
            &          field_H(i+1,j,k-1,isub) + &
            &          field_H(i-1,j+1,k,isub) + &
            &          field_H(i-1,j-1,k,isub) + &
            &          field_H(i-1,j,k+1,isub) + &
            &          field_H(i-1,j,k-1,isub) + &
            &          field_H(i,j+1,k+1,isub) + &
            &          field_H(i,j+1,k-1,isub) + &
            &          field_H(i,j-1,k+1,isub) + &
            &          field_H(i,j-1,k-1,isub) ) + &
          ! d: corners (8)
            & facd * ( field_H(i+1,j+1,k+1,isub) + &
            &          field_H(i+1,j+1,k-1,isub) + &
            &          field_H(i+1,j-1,k+1,isub) + &
            &          field_H(i+1,j-1,k-1,isub) + &
            &          field_H(i-1,j+1,k+1,isub) + &
            &          field_H(i-1,j+1,k-1,isub) + &
            &          field_H(i-1,j-1,k+1,isub) + &
            &          field_H(i-1,j-1,k-1,isub) )
        END DO !i
      END DO !j
    END DO !k
    !--------------------------------------------------------------------------!
    ! save smoothened step function to vorticity field(1) to allow ghosting
    !--------------------------------------------------------------------------!
    DO k=1,ndata(3,isubl)
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          field_H(i,j,k,isub) = field_ubar(1,i,j,k,isub)
        END DO !i
      END DO !j
    END DO !k
  END DO !o

  !----------------------------------------------------------------------------!
  ! calculate gradient of step function field
  !----------------------------------------------------------------------------!
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
  DO k=1,ndata(3,isubl)
    DO j=1,ndata(2,isubl)
      DO i=1,ndata(1,isubl)
        IF (stl_stencil .EQ. 3) THEN
          field_ubar(1,i,j,k,isub) = fac1 * &
             & (  - fac4*field_H(i+4,j,k,isub) + 32.0_mk*field_H(i+3,j,k,isub) &
             & -168.0_mk*field_H(i+2,j,k,isub) +672.0_mk*field_H(i+1,j,k,isub) &
             & -672.0_mk*field_H(i-1,j,k,isub) +168.0_mk*field_H(i-2,j,k,isub) &
             & - 32.0_mk*field_H(i-3,j,k,isub) +    fac4*field_H(i-4,j,k,isub))
          field_ubar(2,i,j,k,isub) = fac2 * &
             & (  - fac4*field_H(i,j+4,k,isub) + 32.0_mk*field_H(i,j+3,k,isub) &
             & -168.0_mk*field_H(i,j+2,k,isub) +672.0_mk*field_H(i,j+1,k,isub) &
             & -672.0_mk*field_H(i,j-1,k,isub) +168.0_mk*field_H(i,j-2,k,isub) &
             & - 32.0_mk*field_H(i,j-3,k,isub) +    fac4*field_H(i,j-4,k,isub))
          field_ubar(3,i,j,k,isub) = fac3 * &
             & (  - fac4*field_H(i,j,k+4,isub) + 32.0_mk*field_H(i,j,k+3,isub) &
             & -168.0_mk*field_H(i,j,k+2,isub) +672.0_mk*field_H(i,j,k+1,isub) &
             & -672.0_mk*field_H(i,j,k-1,isub) +168.0_mk*field_H(i,j,k-2,isub) &
             & - 32.0_mk*field_H(i,j,k-3,isub) +    fac4*field_H(i,j,k-4,isub))
        END IF
        IF (stl_stencil .EQ. 2) THEN
          field_ubar(1,i,j,k,isub) = fac1 * &
             & (         field_H(i+3,j,k,isub) -  9.0_mk*field_H(i+2,j,k,isub) &
             & + 45.0_mk*field_H(i+1,j,k,isub) - 45.0_mk*field_H(i-1,j,k,isub) &
             & +  9.0_mk*field_H(i-2,j,k,isub) -         field_H(i-3,j,k,isub))
          field_ubar(2,i,j,k,isub) = fac2 * &
             & (         field_H(i,j+3,k,isub) -  9.0_mk*field_H(i,j+2,k,isub) &
             & + 45.0_mk*field_H(i,j+1,k,isub) - 45.0_mk*field_H(i,j-1,k,isub) &
             & +  9.0_mk*field_H(i,j-2,k,isub) -         field_H(i,j-3,k,isub))
          field_ubar(3,i,j,k,isub) = fac3 * &
             & (         field_H(i,j,k+3,isub) -  9.0_mk*field_H(i,j,k+2,isub) &
             & + 45.0_mk*field_H(i,j,k+1,isub) - 45.0_mk*field_H(i,j,k-1,isub) &
             & +  9.0_mk*field_H(i,j,k-2,isub) -         field_H(i,j,k-3,isub))
        END IF
        IF (stl_stencil .EQ. 1) THEN
          field_ubar(1,i,j,k,isub) = fac1 * &
             & (-field_H(i+2,j,k,isub) + 8.0*field_H(i+1,j,k,isub) &
             &  -8.0*field_H(i-1,j,k,isub) + field_H(i-2,j,k,isub))
          field_ubar(2,i,j,k,isub) = fac2 * &
             & (-field_H(i,j+2,k,isub) + 8.0*field_H(i,j+1,k,isub) &
             &  -8.0*field_H(i,j-1,k,isub) + field_H(i,j-2,k,isub))
          field_ubar(3,i,j,k,isub) = fac3 * &
             & (-field_H(i,j,k+2,isub) + 8.0*field_H(i,j,k+1,isub) &
             &  -8.0*field_H(i,j,k-1,isub) + field_H(i,j,k-2,isub))
        END IF
        IF (stl_stencil .EQ. 0) THEN
          field_ubar(1,i,j,k,isub) = fac1 * &
             & (field_H(i+1,j,k,isub) - field_H(i-1,j,k,isub))
          field_ubar(2,i,j,k,isub) = fac2 * &
             & (field_H(i,j+1,k,isub) - field_H(i,j-1,k,isub))
          field_ubar(3,i,j,k,isub) = fac3 * &
             & (field_H(i,j,k+1,isub) - field_H(i,j,k-1,isub))
        END IF
      END DO !i
    END DO !j
  END DO !k

  END DO !isub

  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)
  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
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

END SUBROUTINE wvic_stepfunc_smoothen


!------------------------------------------------------------------------------
! Mass and center of mass calculation
!------------------------------------------------------------------------------
SUBROUTINE wvic_calculate_mass(info)
  USE module_wvic
  USE ppm_module_data
  USE ppm_module_map_field_ghost
  IMPLICIT NONE


  INTEGER, INTENT(out)  :: info
  INTEGER               :: i,j,k
  INTEGER               :: isub, isubl
  INTEGER               :: maptype
  REAL(mk)              :: sum_mass,gsum_mass,dv
  REAL(mk)              :: sum_cmassx,gsum_cmassx,sum_cmassy,gsum_cmassy, &
                         & sum_cmassz,gsum_cmassz
  REAL(mk)              :: tx,ty,tz
  CHARACTER(len=256)       :: msg
  INCLUDE 'mpif.h'

  dv=dx*dy*dz
  sum_mass = 0.0_mk
  sum_cmassx = 0.0_mk
  sum_cmassy = 0.0_mk
  sum_cmassz = 0.0_mk
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO k=1,ndata(3,isubl)-1
      DO j=1,ndata(2,isubl)-1
        DO i=1,ndata(1,isubl)-1
          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

          sum_mass = sum_mass + field_H(i,j,k,isub);
          sum_cmassx = sum_cmassx + field_H(i,j,k,isub)*tx;
          sum_cmassy = sum_cmassy + field_H(i,j,k,isub)*ty;
          sum_cmassz = sum_cmassz + field_H(i,j,k,isub)*tz;

        END DO !i
      END DO !j
    END DO !k
  END DO
  sum_mass=sum_mass*dv
  sum_cmassx=sum_cmassx*dv
  sum_cmassy=sum_cmassy*dv
  sum_cmassz=sum_cmassz*dv

!  WRITE (msg,*) 'rank',rank,'mass',sum_mass,'\n'
!  WRITE (unit=0,*) msg
  CALL MPI_AllReduce(sum_mass,gsum_mass,1,mpi_prec,MPI_SUM,comm,info)
  CALL MPI_AllReduce(sum_cmassx,gsum_cmassx,1,mpi_prec,MPI_SUM,comm,info)
  CALL MPI_AllReduce(sum_cmassy,gsum_cmassy,1,mpi_prec,MPI_SUM,comm,info)
  CALL MPI_AllReduce(sum_cmassz,gsum_cmassz,1,mpi_prec,MPI_SUM,comm,info)
!
  object_mass=gsum_mass
  object_cmass(1)=gsum_cmassx/gsum_mass
  object_cmass(2)=gsum_cmassy/gsum_mass
  object_cmass(3)=gsum_cmassz/gsum_mass

END SUBROUTINE wvic_calculate_mass


!------------------------------------------------------------------------------
! Step functions
!------------------------------------------------------------------------------
SUBROUTINE init_stepfunction1 (info)
  USE module_wvic
  INTEGER, INTENT(out)              :: info

  step1A =  step1_interval*0.5_mk + step1_offset
  step1B =  step1_interval*0.5_mk * step1_linearfraction + step1_offset
  step1C = -step1_interval*0.5_mk * step1_linearfraction + step1_offset
  step1E = -step1_interval*0.5_mk + step1_offset

  step1vara = 1.0_mk/ (step1C * step1B - step1C * step1A - step1B ** 2 + &
            & step1A ** 2 + step1E * step1B - step1E * step1A)
  step1varb = -2.0_mk / (step1C * step1B - step1C * step1A - step1B ** 2 + &
            & step1A ** 2 + step1E * step1B - step1E * step1A) * step1A
  step1varc = 1.0_mk/ (step1C * step1B - step1C * step1A - step1B ** 2 + &
              & step1A ** 2 + step1E * step1B - step1E * step1A) * step1A ** 2
  step1vard = 1.0_mk/ (step1A - step1C - step1E + step1B) / (-step1C + step1E)
  step1vare = -2.0_mk / (step1A - step1C - step1E + step1B) / (-step1C + &
            & step1E) * step1E
  step1varj = (-step1C * step1A + step1E * step1A + step1C ** 2 - step1C * &
            & step1B + step1E * step1B) / (step1A - step1C - step1E + step1B) &
            & / (-step1C + step1E)
  step1vark = -2.0_mk / (step1A - step1C - step1E + step1B)
  step1varl = (step1A + step1B) / (step1A - step1C - step1E + step1B)


END SUBROUTINE init_stepfunction1


FUNCTION stepfunction1(xin)
  USE module_wvic
  IMPLICIT NONE
  REAL(mk), INTENT(IN) :: xin
  REAL(mk)             :: x
  REAL(mk) :: stepfunction1

  x =xin

  IF (x .GE. step1A) THEN
    stepfunction1 = 0.0_mk
  ELSEIF (x .GE. step1B) THEN
    stepfunction1 = step1vara*x**2 + step1varb*x + step1varc
  ELSEIF (x .GE. step1C) THEN
    stepfunction1 = step1vark*x + step1varl
  ELSEIF (x .GE. step1E) THEN
    stepfunction1 = step1vard*x**2 + step1vare*x + step1varj
  ELSE
    stepfunction1 = 1.0_mk
  END IF

END FUNCTION stepfunction1

