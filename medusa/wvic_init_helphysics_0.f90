!------------------------------------------------------------------------------!
!* filename: wvic_init_helphysics_0                                           *!
!* project : ppm                                                              *!
!* purpose : initial condition for helical vortex instability                 *!
!*         :                                                                  *!
!* author  : Michael Bergdorf  Philippe Chatelain                             *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 09:45:17 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log $
!
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic INIT HELICAL VORTICES PHYSICS 0 =
!  initialize vorticity here
!------------------------------------------------------------------------------!
SUBROUTINE wvic_init_helphysics_0

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
  INTEGER, PARAMETER       :: mkd = KIND(1.0d0)
  INTEGER(mkd)             :: an, a32 !&
  REAL(mk), DIMENSION(500) :: aky,bky,akz,bkz !&
  INTEGER                  :: nk,kk !&
  REAL(mk)                 :: rnk, amp !&
  INTEGER, DIMENSION(3)    :: ilo, ihi !&
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  ! phi is the phase between the helices
  ! u0 - coordinates of helix
  ! e1, e2, e3 - frenet frame unit vectors = tangent, curve, crossprod
  ! fxy - xy-coordinates of point in frenet frame
  ! vort - strength of vorticity
  ! T - period/pitch
  ! r2 - distance from helix to point squared
  ! Rh - radius of helix, temporary variable
  !----------------------------------------------------------------------------!
  REAL(mk)                 :: ty_center
  REAL(mk)                 :: tx_center
  REAL(mk)                 :: tx, ty, tz
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
  REAL(mk), DIMENSION(helvortex_n*2) :: fac, ra, nrm, radii, phases
  REAL(mk), DIMENSION(5)   :: offsetx, offsety
  REAL(mk)                 :: phi, T
  REAL(mk)                 :: u0x, u0y, u0z
  REAL(mk)                 :: e1x, e1y, e1z
  REAL(mk)                 :: e2x, e2y, e2z
  REAL(mk)                 :: e3x, e3y, e3z
  REAL(mk)                 :: fx, fy, r2
  REAL(mk)                 :: helvortex_periodz ! evt fra Ctrl
  REAL(mk)                 :: vort,Rh
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  ! ss01 & ss02, sss0, s0 - 1st & 2nd root of g' , root of g'' , root of g
  ! sr, sl - interval for bisection
  ! imark, iroot - number of marks and roots to and from the bisection
  ! smarks - the marks to the bisection
  ! gleft,gright - interval to search for roots
  ! nrm, coss, sins - supporting variables
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  CHARACTER(LEN=32)        :: filename
  INTEGER                  :: info, maptype
  INTEGER                  :: i,j,k,o,pp,isub,isubl,iv
  REAL(mk)                 :: sr, sl !&
  INTEGER                  :: smark
  REAL(mk), DIMENSION(10)   :: smarks !?hvormange
  REAL(mk)                 :: sss0, ss01, ss02, s0
  REAL(mk)                 :: gleft, gright
  REAL(mk)                 :: coss, sins
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: mg stuff
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  INTEGER, DIMENSION(4)        :: t_topoid

  COMMON / varsisec / tx, ty, tz, phi, T, Rh

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! - generate noise (angle of dislocation and amplitute of dislocation)
  ! - calculate the position of the vortices using b0 (intervortex distance)
  ! - initialize the vortices using a0 (size of vortices)
  ! - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)
  !pi = ACOS(-1.0_mk)

  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_init_helphysics_0',msg,info)
     STOP ! brutal
  END IF

  isubl = isublist(1)
  ilo = istart(:,isubl) !&
  ihi = ilo + ndata(:,isubl)-1 !&


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
        CALL ppm_write(rank,'wvic_init_physics_5',msg,info)
     END IF
     GOTO 1122
  END IF

  CALL stopwatch("STAR")
  !----------------------------------------------------------------------------!
  ! compute the random amplitudes
  !----------------------------------------------------------------------------!
!!Hvad er dette!?

  IF (helvortex_period.EQ.0.0_mk) THEN
    helvortex_period = max_physg(3) - min_physg(3)
  END IF

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! ty/tx_center is the center of the helix
  ! fac = gamma/(2 pi sigma^2)
  ! ra = 1/(2 sigma^2)
  !----------------------------------------------------------------------------!
  length    = max_physg - min_physg !!skal denne bruges overhovedet? 
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))
  helvortex_periodz=1
  T=helvortex_period
  gleft=-T*helvortex_periodz
  gright=T*(2+helvortex_periodz)
  !----------------------------------------------------------------------------!
  ! determining offset for neighbouring boxes (periodic bc's)
  !----------------------------------------------------------------------------!
  offsetx(1)=0_mk
  offsety(1)=0_mk
  offsetx(2)=length(1)
  offsety(2)=0_mk
  offsetx(3)=-length(1)
  offsety(3)=0_mk
  offsetx(4)=0_mk
  offsety(4)=length(2)
  offsetx(5)=0_mk
  offsety(5)=-length(2)

  !----------------------------------------------------------------------------!
  ! calculating predefined, multiply used, variables for each helix
  !----------------------------------------------------------------------------!
  DO iv=1,helvortex_n
    fac(iv*2-1) = helvortex_gamma/(helvortex_a1**2_mk * M_PI )
    fac(iv*2) = helvortex_ratio * helvortex_gamma/(helvortex_a2**2_mk * M_PI )
    ra(iv*2-1) = 1.0_mk/(helvortex_a1**2_mk )
    ra(iv*2) = 1.0_mk/(helvortex_a2**2_mk )
    radii(iv*2-1) = helvortex_r1
    radii(iv*2) = helvortex_r2
    nrm(iv*2-1) = 1.0_mk/sqrt(4*helvortex_r1**2*M_PI**2+T**2);
    nrm(iv*2) = 1.0_mk/sqrt(4*helvortex_r2**2*M_PI**2+T**2);
    phases(iv*2-1) = 2.0_mk * M_PI / helvortex_n * REAL(iv-1,mk)
    phases(iv*2) = 2.0_mk * M_PI / helvortex_n * REAL(iv-1,mk)
  END DO
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! create vortices
  !----------------------------------------------------------------------------!
  max_vorticity = MAX(fac(1),fac(2))
  !----------------------------------------------------------------------------!
  ! For each subdomain assigned to the processor...
  ! isubl is the index of the currently treated subdomain
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
    isubl = isublist(isub)
    write(*,*) 'starting subdomain ', isubl
    DO k=1,ndata(3,isubl)
    !-------------------------------------------------------------------------!
    ! do for each Z-layer of the subdomain...
    ! i,j,k are counters for x, y, z
    !-------------------------------------------------------------------------!
      DO j=1,ndata(2,isubl)
        DO i=1,ndata(1,isubl)
          !-------------------------------------------------------------------!
          ! The two loops covers all cells in the z-layer
          ! tx, ty, tz are the coordinates (closest to O) of the cell
          ! tx/ty_center is subtracted since the helix revolves around (0,0)
          !-------------------------------------------------------------------!
          twp = 0.0_mk
          DO o=1,5
            tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center + offsetx(o)
            ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center + offsety(o)
            tz = min_sub(3,isubl) + REAL(k-1,mk)*dz
  
            DO iv=1,helvortex_n*2
              Rh=radii(iv)
              phi=phases(iv)
              ! solutions are in interval -3/2T : 1/2T
              sss0= T * (-phi + ATAN2(ty , tx)) /(M_PI * 2_mk);
              smark=1
              smarks(smark)=gleft
              
              IF ((SISEC2(sss0)*SISEC2(sss0+T/2_mk) .LT. 0) .and. (ty .ne. 0.0) .and. (tx .ne. 0.0)) THEN
                ! solutions are in interval -2T : 1T
                ss01=zriddr(sisec2,sss0-T/2_mk,sss0,helvortex_precision_isec2,2)
                ss02=zriddr(sisec2,sss0,sss0+T/2_mk,helvortex_precision_isec2,3)
                ! to cover the domain only the above interval has to be added
                ! 2 periods forward only
                ! to cover the domain and an extra period forward and backward
                ! 1 period has to be added backward and 3 forward
                ! -to the above the endpoints have to be added
                ! -and they need to be filtered/sorted if not in desired domain
                DO pp=-helvortex_periodz,(2+helvortex_periodz) !has been tested in MATLAB
                  IF ((ss01+REAL(pp,mk)*T-gleft)*(ss01+REAL(pp,mk)*T-gright) .LT. 0) THEN
                  ! the point is gleft < point < gright              
                    smark=smark+1
                    smarks(smark)=ss01+pp*T
                  END IF
                  IF ((ss02+REAL(pp*T,mk)-gleft)*(ss02+REAL(pp*T,mk)-gright) .LT. 0) THEN
                  ! the point is gleft < point < gright              
                    smark=smark+1
                    smarks(smark)=ss02+pp*T
                  END IF
                END DO
              END IF
              smark=smark+1
              smarks(smark)=gright
  
              DO pp=1,(smark-1)
                IF (SISEC(smarks(pp))*SISEC(smarks(pp+1)) .LT. 0) THEN
                  s0=zriddr(SISEC,smarks(pp),smarks(pp+1),helvortex_precision_isec1,1)!,outstring)
                  coss=COS( 2_mk * M_PI * s0 / T + phi )
                  sins=SIN( 2_mk * M_PI * s0 / T + phi )
                  u0x = Rh * coss
                  u0y = Rh * sins
                  u0z = s0;
                  e1x = -2_mk * M_PI * Rh * sins * nrm(iv)
                  e1y =  2_mk * M_PI * Rh * coss * nrm(iv)
                  e1z = T * nrm(iv)
                  e2x = - coss
                  e2y = - sins
                  e2z = 0;
                  e3x = -T*sins*nrm(iv)
                  e3y =  T*coss*nrm(iv)
                  e3z = -2_mk*Rh*M_PI*nrm(iv)
                  r2 = (tx-u0x)**2+(ty-u0y)**2+(tz-u0z)**2
                  fx = ((tx-u0x)*e2x+(ty-u0y)*e2y)
                  fy = ((tx-u0x)*e3x+(ty-u0y)*e3y+(tz-u0z)*e3z)
                  !magnitude of vorticity
                  vort = fac(iv) * (sqrt(16_mk*Rh**2*M_PI**4*(fx**2-2_mk*fx*Rh+Rh**2)+&
                  4_mk*T**2*M_PI**2*(fx**2-2_mk*fx*Rh+2_mk*Rh**2+fy**2)+T**4) &
                  * nrm(iv)**2) * EXP( -r2 * ra(iv))

                  !vort = fac(iv) * EXP( -r2 * ra(iv))

                  twp(1) =  twp(1) + e1x*vort
                  twp(2) =  twp(2) + e1y*vort
                  twp(3) =  twp(3) + e1z*vort
                END IF
              END DO

            END DO !iv
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
  write(*,*) "finished subdomain, ", isubl
  CALL stopwatch("STOP")

 !---- Vorticity has been initialized on the field.
  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)
  CALL ppm_fft_solenoidal(field_wp,mesh_id,topo_id,t_topoid,fmesh_id, & 
     & ghostsize, info)

  !----------------------------------------------------------------------------!
  ! get ghosts for the new vorticity
  !----------------------------------------------------------------------------!
  CALL ppm_write(rank,'wvic_init_helphysics_0','ghosting',info)
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
  CALL ppm_write(rank,'wvic_init_helphysics_0','ghst complete',info)

1122 CONTINUE

  !----------------------------------------------------------------------------!
  ! create particles with style
  !----------------------------------------------------------------------------!
  CALL ppm_rmsh_create_part(xp,np,wp,lda,field_wp,topo_id,mesh_id,&
       & (/cutoff,HUGE(cutoff)/),info,resetpos=.TRUE.,cutoff_weights=cow)
  WRITE(msg,*) ' created ',np,' particles'
  IF(rank.EQ.0) CALL ppm_write(rank,'wvic_init_helphysics_0',msg,info)
END SUBROUTINE wvic_init_helphysics_0

!-------purely debugging---------------
! IF ((SISEC2(sss0)*SISEC2(sss0-T/2_mk) .GT. 0)) THEN
!   write(filename,'(a,i4.4)') 'debug',rank
!   open(10,file=filename,position='append')
!   write(10,*) 'error, T=' , T , ' phi=' , phi , ' R=' , Rh , ' ty=' , ty , ' tx=' , tx , ' tz=' , tz , ' sss0=' , sss0 , ' sisec2(-)=', sisec2(sss0-T/2_mk)  , ' sisec2()=', sisec2(sss0) , ' sisec2(-)=', sisec2(sss0+T/2_mk)
!   close(10)
! END IF

