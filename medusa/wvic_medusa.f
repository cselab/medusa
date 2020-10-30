SUBROUTINE wvic_readmedusa

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  INTEGER                  :: info

  CHARACTER(len=256) :: msg
  INTEGER, PARAMETER :: md = kind(2.0d0)

  INTEGER                  :: ilenmedusa,iline
  INTEGER                  :: iUnit
  INTEGER                  :: ILEN,ios,filenumber
  CHARACTER(LEN=256)       :: cbuf, tmpbuf, filein
  LOGICAL                  :: lExist
  REAL(mk), DIMENSION(4)   :: pointvel
  REAL(mk), DIMENSION(2)   :: point,pointold,panel,vel
  REAL(mk)                 :: panelsq

  !-----------------------------------------------------------------------------
  ! Define frame time step -TODO
  !-----------------------------------------------------------------------------
  medusa_dt      = medusa_period/REAL(medusa_frames,mk)
  medu_mass_orig = medu_mass_orig * medusa_scale**3
  medu_move_vel  = 0.0_mk
  medu_move_pos  = 0.0_mk
  !-----------------------------------------------------------------------------
  ! Definition of file unit
  !-----------------------------------------------------------------------------
  iUnit = 30

  !-----------------------------------------------------------------------------
  ! open medusa outline file
  !-----------------------------------------------------------------------------
  ilenmedusa = LEN_TRIM(medusa_file)
  WRITE(filein,'(A,I5.5)') medusa_file(1:ilenmedusa), 0
  INQUIRE(FILE=filein(1:(ilenmedusa+5)), EXIST=lExist)
  IF(.NOT. lExist) THEN
    WRITE(unit=0,'(2A)')'\nFile does not exist: ',filein(1:(ilenmedusa+5)), '\n'
    Info = 1
    call mpi_finalize(info)
    stop
    call wvic_died
    GOTO 9999
  END IF

  OPEN(iUnit, FILE=filein(1:(ilenmedusa+5)), IOSTAT=ios, ACTION='READ')
  IF(ios .NE. 0) THEN
    WRITE(unit=0,'(3A)')'\nFailed to open medusa file: ', &
     & filein(1:(ilenmedusa+5)), '\n'
    Info = 1
    GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  ! Count points
  !-----------------------------------------------------------------------------
!JTR medusa_count is not intuitive. consider using medusa_npoints
  medusa_count = 0
  DO
    READ(iUnit,'(A)',END=50,ERR=200) cbuf
    medusa_count = medusa_count + 1
  END DO

50  CLOSE(iUnit)

  !-----------------------------------------------------------------------------
  ! Allocate vectors, dotproduct arrays
  !-----------------------------------------------------------------------------
  ALLOCATE(medusa_points_array(medusa_frames,medusa_count,2))
  ALLOCATE(medusa_vel_array(medusa_frames,medusa_count,2))
  ALLOCATE(medusa_points(medusa_count,2))
  ALLOCATE(medusa_vel(medusa_count,2))
  ALLOCATE(medusa_panel(medusa_count-1,2))
  ALLOCATE(medusa_panel_squared(medusa_count-1))
  ALLOCATE(medusa_normal(medusa_count-1,2))
  ALLOCATE(medusa_maxx(medusa_frames))
  ALLOCATE(medusa_maxy(medusa_frames))
  ALLOCATE(medusa_miny(medusa_frames))

  !-----------------------------------------------------------------------------
  ! read file
  !-----------------------------------------------------------------------------
  DO filenumber=1,medusa_frames
    WRITE(filein,'(A,I5.5)') medusa_file(1:ilenmedusa), (filenumber-1)
    OPEN(iUnit, FILE=filein(1:(ilenmedusa+5)), IOSTAT=ios, ACTION='READ')
    IF(ios .NE. 0) THEN
      WRITE(unit=0,'(3A)')'\nFailed to open file:',filein(1:(ilenmedusa+5)),'\n'
      Info = 1
      GOTO 9999
    END IF

    DO iline=1,medusa_count 
      READ(iUnit,'(A)',END=100,ERR=200) cbuf
      cbuf = ADJUSTL(cbuf)
      ILEN = LEN_TRIM(cbuf)

      CALL UpperCase(cbuf,ILEN)
      tmpbuf = cbuf(1:ILEN)
      tmpbuf = ADJUSTL(tmpbuf)
      READ(tmpbuf,*,iostat=ios,err=200) pointvel
      vel    = pointvel(3:4)*medusa_scale/(2.0_mk*medusa_dt)
      point  = pointvel(1:2)*medusa_scale

      !------------------------------------------------------------------------
      ! Determine x and y bounds of each frame
      !------------------------------------------------------------------------
      IF (iline .EQ. 1) THEN
        medusa_maxx(filenumber) = point(1)
        medusa_maxy(filenumber) = point(2)
        medusa_miny(filenumber) = point(2)
      END IF
      medusa_maxx(filenumber) = MAX(medusa_maxx(filenumber),point(1))
      medusa_maxy(filenumber) = MAX(medusa_maxy(filenumber),point(2))
      medusa_miny(filenumber) = MIN(medusa_miny(filenumber),point(2))

      !JTR WTF? are we moving point(1) slightly to the left to close the curve?
      IF (iline .EQ. medusa_count) THEN
        point(1) = point(1) - 0.0001*medusa_scale
      END IF
      medusa_points_array(filenumber,iline,:) = point
      medusa_vel_array(filenumber,iline,:)    = vel
    END DO !iline
  CLOSE(iUnit)
  END DO !filenumber

  GOTO 9999

200 CONTINUE
  WRITE(*,'(A,I5,2A)') '\nError reading line: ',iline,                      &
       &                 ' of file: ',filein(1:(ilenmedusa+5))
  ILEN = LEN_TRIM(cbuf)
  WRITE(*,'(A)') cbuf(1:ILEN)
  Info = 1
  GOTO 9999

  !-----------------------------------------------------------------------------
  !  End of file
  !-----------------------------------------------------------------------------
100 Info = 0

  !-----------------------------------------------------------------------------
  !  Close file
  !-----------------------------------------------------------------------------

  IF (rank .EQ. 0) THEN
    WRITE(msg,*) '\nrank:',rank,'Succesfully read', medusa_count, &
      & ' points from ', medusa_file(1:ilenmedusa),'\n'
    WRITE(unit=0,*) msg
    WRITE(msg,*)'\nmaxx:',medusa_maxx,'maxy:',medusa_maxy,'miny:', &
      & medusa_miny,'\n'
    WRITE(unit=0,*) msg
  END IF

  !-----------------------------------------------------------------------------
  !  Return
  !-----------------------------------------------------------------------------
9999 CONTINUE


END SUBROUTINE wvic_readmedusa


SUBROUTINE wvic_stepfunc_medusa

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
  ! level = levelset/epsilon
  ! epsilon = 1/epsilon (epsilon ~ dx)
  !----------------------------------------------------------------------------!
  !JTR clean up variables
  REAL(mk)                 :: tx_center
  REAL(mk)                 :: ty_center
  REAL(mk)                 :: tz_center
  REAL(mk)                 :: tx, ty, tz
  REAL(mk)                 :: px, py
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(2)   :: vecP,panel
  REAL(mk)                 :: level, epsilon
  REAL(mk)                 :: dist,pdotpanel,panel_part
  REAL(mk)                 :: cosa,sina,cosangxz,sinangxz
  REAL(mk)                 :: tst1,tst2
  REAL(mk)                 :: pola,polb,polc
  REAL(mk)                 :: x1,y1,x2,y2,vx1,vy1,vx2,vy2
  REAL(mk)                 :: panel_slope, prev_point
  INTEGER                  :: panel_xing,lastinout,p,ang,np1,np1cand,np2
  INTEGER                  :: npointshalf
  INTEGER                  :: ibmin,ibmax,jbmin,jbmax,kbmin,kbmax
  INTEGER                  :: imin,imax,jmin,jmax,kmin,kmax
  LOGICAL                  :: inside,offofpanel,offofpanel_cand
  REAL(mk)                 :: ubarmax
  !----------------------------------------------------------------------------!
  INCLUDE 'mpif.h'

  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  ! i,j,k = counter for spatial loop
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: msg
  CHARACTER(len=512)       :: msgl
  INTEGER                  :: maptype, info
  INTEGER                  :: i,j,k,isub,isubl
  !----------------------------------------------------------------------------!
  CHARACTER(len=256)       :: filename
  INTEGER                  :: ios

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff ??
  INTEGER, DIMENSION(3)    :: ftopo_id
  INTEGER, DIMENSION(4)    :: fmesh_id
  INTEGER, DIMENSION(4)    :: t_topoid

  ftopo_id = (/2,3,4/)
  fmesh_id = (/2,3,4,5/)
  t_topoid = (/2,3,4,5/)

  !----------------------------------------------------------------------------!
  ! Interpolate between medusa frames
  !----------------------------------------------------------------------------!
  CALL wvic_interpolate_medusa
  
  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  ! tx/ty/tz_center is the center of the sphere
  ! ra = 1/(2 sigma^2)
  !----------------------------------------------------------------------------!
  field_H     = 0.0_mk
  field_ubar  = 0.0_mk
  length      = max_physg - min_physg
  tx_center   = 0.5_mk*(min_physg(1) + max_physg(1))+object_offset(1)
  ty_center   = 0.5_mk*(min_physg(2) + max_physg(2))+object_offset(2)
  tz_center   = 0.5_mk*(min_physg(3) + max_physg(3))+object_offset(3)
  !JTR get rid of stepfunction_band
  epsilon     = 1.0_mk / (sqrt(dx**2 + dy**2 + dz**2)*stepfunction_band)
  !JTR what if medusa_count is not divisable by 2 - it is ALWAYS! but do a check
  npointshalf = (medusa_count)/2
  !----------------------------------------------------------------------------!
  ! For each subdomain assigned to the processor...
  ! isubl is the index of the currently treated subdomain
  !JTR? what is ibmin/max... ijkminmax?
  ! ibmin - is local lower bound (index) in x-direction
  ! jbmin - is local lower bound (index) in y-direction
  ! kbmin - is local lower bound (index) in z-direction
  ! ibmax - is local upper bound (index) in x-direction
  ! jbmax - is local upper bound (index) in y-direction
  ! kbmax - is local upper bound (index) in z-direction
  ! the + in jbmin is ok as mminy just gives the lower bound in the relative 
  ! coordinates who are added to the centre of the domain. From this the local 
  ! minimum is subtracted, the result divided by dx/dy/dz etc. and rounded down
  ! to get the index of the field array. ibmin and kbmin has a - (minus) due to 
  ! mmax being the UPPER bound.
  ! Finally imin...kmax is created to ensure that the field array indicies do 
  ! not extend the local array bounds
  ! 
  ! If the ghost layer is included no ghosting should be nescessary
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
    isubl  = isublist(isub)
    i      = ceiling(sqrt(3.0_mk)*(step1_interval*0.5_mk + step1_offset))
    !TODO : change ibmin to imin etc.
    ibmin  = FLOOR( (tx_center-mmaxx-min_sub(1,isubl))/dx+0.5_mk)-1-i
    jbmin  = FLOOR( (ty_center+mminy-min_sub(2,isubl))/dy+0.5_mk)-1-i!+ is ok
    kbmin  = FLOOR( (tz_center-mmaxx-min_sub(3,isubl))/dz+0.5_mk)-1-i
    ibmax  = FLOOR( (tx_center+mmaxx-min_sub(1,isubl))/dx+0.5_mk)+1+i
    jbmax  = FLOOR( (ty_center+mmaxy-min_sub(2,isubl))/dy+0.5_mk)+1+i
    kbmax  = FLOOR( (tz_center+mmaxx-min_sub(3,isubl))/dz+0.5_mk)+1+i
    imin   = max(ibmin,1-ghostsize(1))
    jmin   = max(jbmin,1-ghostsize(2))
    kmin   = max(kbmin,1-ghostsize(3))
    imax   = min(ibmax,ndata(1,isubl)+ghostsize(1))
    jmax   = min(jbmax,ndata(2,isubl)+ghostsize(2))
    kmax   = min(kbmax,ndata(3,isubl)+ghostsize(3))
    !-------------------------------------------------------------------------!
    ! do for each Z-layer of the subdomain...
    ! i,j,k are counters for x, y, z
    !-------------------------------------------------------------------------!
    !JTR money to be made here: That is the local bounds are not being heeded
    !JTR did anyone say pencil?
!    DO k=1-ghostsize(3),ndata(3,isubl)+ghostsize(3)!kmin,kmax !
!      DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)!jmin,jmax !
!        DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)!imin,imax !
    WRITE(msg,*) '\nindex:',rank,imin,imax,jmin,jmax,kmin,kmax
    WRITE(UNIT=0,*) msg
    WRITE(msg,*) '\nbounds:',rank,mmaxx,mminy,mmaxy
    WRITE(UNIT=0,*) msg

    DO k=kmin,kmax,1 
      DO j=jmin,jmax,1 
        DO i=imin,imax,1 
          !-------------------------------------------------------------------!
          ! The j and i loops covers all cells in the z-layer
          ! tx, ty, tz are the coordinates in a coordinate system with origo in
          ! the center of the sphere (tx/ty/tz_center)
          !-------------------------------------------------------------------!
          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center
          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center
          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center

          !JTR? We are getting cylindrical coordinates here?
          !JTR if this is not done in pencils, do it!
          !PEEENCIIIL!!!!!!!!!
          px = sqrt(tx**2 + tz**2)
          py = ty
          !---------------------------------------------------------------------
          ! Compute projection components to project radial velocities along
          ! the x and z directions
          !---------------------------------------------------------------------
          sinangxz = ATAN2(tz,tx)
          tst1     = sinangxz
          cosangxz = COS(sinangxz)
          sinangxz = SIN(sinangxz)
          IF (rank .eq. -10) THEN
            WRITE(unit=0,*) '\nCOORDS ',tx,ty,tz,px,py,tst1,cosangxz,sinangxz
          END IF
          IF (px .EQ. 0.0_mk) THEN
            !JTR consider printing cosangxz,sinangxz to see what ATAN2 gives
            !points coinciding with rotational axis should have no radial vel.
            cosangxz = 0.0_mk
            sinangxz = 0.0_mk
          ENDIF

          !---------------------------------------------------------------------
          ! DETERMINE IF THE POINT IS INSIDE OR OUTSIDE THE MEDUSA
          !---------------------------------------------------------------------
          panel_xing = 0
          prev_point = -1.0_mk
          DO p=2,medusa_count
            !-------------------------------------------------------------------
            ! Check if px is over a medusa point (then continue to next point)
            ! But if p is the last point the comparison is always done
            !-------------------------------------------------------------------
            IF ((px .EQ. medusa_points(p,1)) .AND. &
              & (p .NE. medusa_count)) THEN
              CYCLE
            !-------------------------------------------------------------------
            ! Else check if the point is between the two points -
            ! and if so check if the point is above or below the panel(s)
            ! Also if p is the last point the comparison is always done
            ! JTR - this might be done just by changing LT to LE...
            !-------------------------------------------------------------------
            ELSEIF ((px-prev_point)*(px-medusa_points(p,1)) .LE. 0.0_mk) THEN
              !-----------------------------------------------------------------
              ! More explanation here. The slope should also be calculated in
              ! the interpolation routine !JTR!!
              ! the slope is calculated backwards (from p to p-1)
              !-----------------------------------------------------------------
              panel_slope = (medusa_points(p-1,2)-medusa_points(p,2)) / &
                            (medusa_points(p-1,1)-medusa_points(p,1))
              IF (panel_slope*(px-medusa_points(p,1)) + medusa_points(p,2) &
                   .GT. py) THEN
                panel_xing  = panel_xing  + 1
              END IF
            END IF
            prev_point = medusa_points(p,1)
          END DO

          IF (MOD(panel_xing,2) .NE. 0) THEN
            inside = .true.
          ELSE
            inside = .false.
          END IF

          !---------------------------------------------------------------------
          ! DETERMINE MINIMUM DISTANCE
          !---------------------------------------------------------------------
          ! vectors: vecP (basepoint to P), panel (panel), t (panel unit vect)
          ! P just off panel => vecP * t              < len(panel)     =>
          !                     vecP * t * len(panel) < len(panel)^2
          ! and  t * len(panel) = panel
          !                     vecP * panel < len(panel)^2
          ! Also                vecP * panel > 0
          !---------------------------------------------------------------------
          level      = SUM(length**2) 
          offofpanel = .FALSE.
          DO p=1,medusa_count-1
            offofpanel_cand = .FALSE.
            vecP(1)    = px - medusa_points(p,1)
            vecP(2)    = py - medusa_points(p,2)
            Pdotpanel  = vecP(1)*medusa_panel(p,1) + vecP(2)*medusa_panel(p,2)
            !-------------------------------------------------------------------
            ! Point is closest to endpoint
            IF (Pdotpanel .GE. medusa_panel_squared(p)) THEN
              np1cand = p+1
              vecP    = vecP - medusa_panel(np1cand,:)
              dist    = SUM(vecP**2)
            !-------------------------------------------------------------------
            ! Point is closest to basepoint
            ELSEIF (Pdotpanel .LE. 0.0_mk) THEN
              np1cand = p
              dist    = SUM(vecP**2)
            !-------------------------------------------------------------------
            ! Point is closest to panel
            ELSE
              dist = (vecP(1)*medusa_normal(p,1)+vecP(2)*medusa_normal(p,2))**2

              offofpanel_cand = .TRUE.
              np1cand         = p
              !this may not be correct. maybe p should also be used:
              !IF (Pdotpanel .LE. 0.5_mk*medusa_panel_squared(p)) THEN
                np1cand = p
              !ELSE
                !np1cand = p+1
              !END IF
            ENDIF

            IF (dist .LT. level) THEN
              level = dist
              np1   = np1cand
              IF (offofpanel_cand) THEN
                panel_part = Pdotpanel/medusa_panel_squared(p)
                offofpanel = .TRUE.
              ELSE
                offofpanel = .FALSE.
              END IF
            END IF
          END DO

          !---------------------------------------------------------------------
          ! ASSIGN VELOCITY TO FIELD
          !---------------------------------------------------------------------
          ! For points interior to the medusa AND points not too close to the
          ! boundary. Points very close to the boundary cause trouble due to
          ! truncation. If a very small or great medusa scaling is used care
          ! should be taken - the level should be compared to medusa_scaling
          ! No 'ERROR's have been seen for abs(level) greater than 0.28E-11 and
          ! the medusa is of the order O(medusa_scaling) so 1E-8 should be ok
          !---------------------------------------------------------------------
          IF (inside .AND. (level .GT. 1.0E-8_mk)) THEN
            !-------------------------------------------------------------------
            ! Determine points for interpolation (See sketch further below)
            !-------------------------------------------------------------------
            IF (np1 .GT. npointshalf) THEN
              np1=npointshalf*2+1-np1
            ENDIF
            IF (np1 .EQ. npointshalf) THEN
              np1=np1-1
            END IF
            np2 = 2*npointshalf+1-np1-1

            !-------------------------------------------------------------------
            ! IDENTIFY THE CELL CONTAINING THE POINT
            !
            !    ---+------------+---
            !       | \         /|\
            !       |  \ vp      |
            !    v1 |   \        | v2
            !       |    _|      |
            !      \|/           |
            !    ---+------------+---
            !
            ! if -^v1 * vp .GE. 0 the cell we are looking for is to the left
            ! if -^v1 * vp .GT. 0 the cell we are looking for is to the right
            ! the first dot product is stored in tst1 the second in tst2
            !-------------------------------------------------------------------
            tst1 = (px - medusa_points(np1,1))* &
            &(medusa_points(np2,2)+medusa_panel(np2,2)-medusa_points(np1,2))-&
            &(py - medusa_points(np1,2))* &
            &(medusa_points(np2,1)+medusa_panel(np2,1)-medusa_points(np1,1))

            DO WHILE (tst1 .GE. 0.0_mk)
              !WE ENCOUNTERED THE BUTCHER MOVE LEFT
              np1 = np1-1
              np2 = 2*npointshalf+1-np1-1
              IF (np1 .LT. 1) THEN
                !Dump Medusa outline
                WRITE(filename,'(A,I5.5,A,I5.5)') 'outlineR',rank,'I',itime
                OPEN(14,file=filename,iostat=ios,position='append',status='new')
                DO p=1,medusa_count
                  WRITE(14,*) medusa_points(p,:), medusa_vel(p,:)
                END DO
                CLOSE(14)
                WRITE(msgl,*) '\nERROR in wvic_stepfunction_medusa, np1 < 1:', &
                  & 'rank',rank,'ijk',i,j,k,'xyz,pxy,',tx,ty,tz,px,py, &
                  & 'tst1',tst1,'np1,np2',np1,np2,'level',level, &
                  & "panels (x1x1'x2x2')", &
                  & medusa_points(np1+1,1),&
                  & medusa_points(np1+1,2),&
                  & medusa_points(np1+1,1)+medusa_panel(np1+1,1),&
                  & medusa_points(np1+1,2)+medusa_panel(np1+1,2),&
                  & medusa_points(np2-1,1),&
                  & medusa_points(np2-1,2),&
                  & medusa_points(np2-1,1)+medusa_panel(np2-1,1),&
                  & medusa_points(np2-1,2)+medusa_panel(np2-1,2),'\n'
                WRITE(unit=0,*) TRIM(msgl)
                STOP
              ENDIF
              tst1 = (px - medusa_points(np1,1))* & 
              &(medusa_points(np2,2)+medusa_panel(np2,2)-medusa_points(np1,2))-&
              &(py - medusa_points(np1,2))* & 
              &(medusa_points(np2,1)+medusa_panel(np2,1)-medusa_points(np1,1))
            END DO

            tst2 = (px - medusa_points(np2,1))* & 
            &(medusa_points(np1,2)+medusa_panel(np1,2)-medusa_points(np2,2))-&
            &(py - medusa_points(np2,2))* & 
            &(medusa_points(np1,1)+medusa_panel(np1,1)-medusa_points(np2,1))
            DO WHILE (tst2 .GT. 0.0_mk)
              !NO TOWARDS THE STAIRS MOVE RIGHT
              np1 = np1+1
              np2 = 2*npointshalf-np1
              IF (np1 .GE. npointshalf) THEN
                !Dump Medusa outline
                WRITE(filename,'(A,I5.5,A,I5.5)') 'outlineR',rank,'I',itime
                OPEN(14,file=filename,iostat=ios,position='append',status='new')
                DO p=1,medusa_count
                  WRITE(14,*) medusa_points(p,:), medusa_vel(p,:)
                END DO
                CLOSE(14)
                WRITE(msgl,*) '\nERROR in wvic_stepfunction_medusa, np1 >= ', &
                  & 'np/2: rank',rank,'ijk',i,j,k,'xyz,pxy,',tx,ty,tz,px,py, &
                  & 'tst1',tst1,'tst2',tst2,'np1,np2',np1,np2,'level',level, &
                  & "panels (x1x1'x2x2')", &
                  & medusa_points(np1-1,1),&
                  & medusa_points(np1-1,2),&
                  & medusa_points(np1-1,1)+medusa_panel(np1-1,1),&
                  & medusa_points(np1-1,2)+medusa_panel(np1-1,2),&
                  & medusa_points(np2+1,1),&
                  & medusa_points(np2+1,2),&
                  & medusa_points(np2+1,1)+medusa_panel(np2+1,1),&
                  & medusa_points(np2+1,2)+medusa_panel(np2+1,2),'\n'
                WRITE(unit=0,*) TRIM(msgl)
                STOP
                STOP
              ENDIF
              tst2 = (px - medusa_points(np2,1))* & 
              &(medusa_points(np1,2)+medusa_panel(np1,2)-medusa_points(np2,2))-&
              &(py - medusa_points(np2,2))* & 
              &(medusa_points(np1,1)+medusa_panel(np1,1)-medusa_points(np2,1))
            END DO

            !This could be a spot for testing tst1 and tst2 again
            tst1 = (px - medusa_points(np1,1))* & 
            &(medusa_points(np2,2)+medusa_panel(np2,2)-medusa_points(np1,2))-&
            &(py - medusa_points(np1,2))* & 
            &(medusa_points(np2,1)+medusa_panel(np2,1)-medusa_points(np1,1))
            tst2 = (px - medusa_points(np2,1))* & 
            &(medusa_points(np1,2)+medusa_panel(np1,2)-medusa_points(np2,2))-&
            &(py - medusa_points(np2,2))* & 
            &(medusa_points(np1,1)+medusa_panel(np1,1)-medusa_points(np2,1))

            IF ((tst1 .GT. 0.0_mk) .OR. (tst2 .GT. 0.0_mk)) THEN
                WRITE(msg,*) '\nERROR tst1 or tst2 positive\n'
                WRITE(unit=0,*) msg
            END IF

            !INNER SANCTUM REACHED, TAKE HEED AND BEAR WITNESS TO THE TRUTH
            !-------------------------------------------------------------------
            ! DETERMINE INTERPOLATION VARIABLES a AND b
            ! Minimise array access, store in x1,y1,...
            ! determine coefficients to polynomial and solve to get 
            ! interpolation variable 'b' between medusa surfaces(cf Maple sheet)
            ! pola,polb,polc are constants in the 2nd degree polynomial.
            ! Order of the points here and in interp_coord_aligned_fortran.mw:
            !            (vx1)
            !           ----->
            !        --x1----x1'---
            !        (np1)
            !               (np2)
            !        --x2'---x2----
            !           <-----
            !            (vx2)
            !-------------------------------------------------------------------
            x1  = medusa_points(np1,1)
            y1  = medusa_points(np1,2)
            x2  = medusa_points(np2,1)
            y2  = medusa_points(np2,2)
            vx1 = medusa_panel(np1,1)
            vy1 = medusa_panel(np1,2)
            vx2 = medusa_panel(np2,1)
            vy2 = medusa_panel(np2,2)

            !pola = x1*vy2-vx1*y1+vx1*y2+vx1*vy2-vx2*y1+vx2*y2+x1*vy1 - &
            !     & vx2*vy1-x2*vy1-x2*vy2
            !polb = -vx1*py+px*vy2+2.0_mk*vx1*y1+vx2*vy1+px*vy1-x1*vy2 + &
            !     & vx2*y1-vx1*y2-vx1*vy2+x2*vy1-vx2*py-2*x1*vy1
            !polc = -px*vy1-vx1*y1+vx1*py+x1*vy1
            pola = -x2*vy2-vx2*vy1+x1*vy2+x1*vy1+vx2*y2-vx2*y1-vx1*y1-x2*vy1+ &
                 & vx1*y2+vx1*vy2
            polb = px*vy1-x1*vy2-vx1*y2+vx2*vy1+px*vy2+vx2*y1-vx1*py-vx1*vy2- &
                 & 2*x1*vy1-vx2*py+x2*vy1+2*vx1*y1
            polc = x1*vy1-px*vy1-vx1*y1+vx1*py

            !JTR consider testing determinant and/or pola
            !-------------------------------------------------------------------
            ! Solve polynamial for b and check for a solution between 0 and 1
            ! Temporarily store a and b in x1 and x2 respectively
            !-------------------------------------------------------------------
            tst1 = (-polb + sqrt(polb**2 - 4.0_mk*pola*polc))/(2.0_mk*pola)
            tst2 = (-polb - sqrt(polb**2 - 4.0_mk*pola*polc))/(2.0_mk*pola)

            !JTR maybe add test on the discriminant
            !JTR at some point test if perhaps only tst1 or tst2 is used!!
            IF ((tst1 .GE. 0.0_mk) .AND. (tst1 .LE. 1.0_mk)) THEN
              !should be ok when the vectors point in opposite directions
              x1 = (x1-tst1*x1+tst1*x2+tst1*vx2-px)/(-vx1+tst1*vx1+tst1*vx2)
              x2 = tst1
            ELSEIF ((tst2 .GE. 0.0_mk) .AND. (tst2 .LE. 1.0_mk)) THEN
              x1 = (x1-tst2*x1+tst2*x2+tst2*vx2-px)/(-vx1+tst2*vx1+tst2*vx2)
              x2 = tst2
            ELSE
              !-----------------------------------------------------------------
              ! I am not sure why I do this / if it has been a problem at all.
              ! But I solve the polynomial for a = 0
              !-----------------------------------------------------------------
              WRITE(msgl,*) '\nINFO: switching to linear equation: itime=', &
                & itime,', rank ',rank,'ijk',i,j,k,'pol abc',pola,polb,polc, &
                & 'xyz',tx,ty,tz,'level',level,'roots',tst1,tst2,'coords', &
                & x1,y1,x2,y2,vx1,vy1,vx2,vy2,px,py,'\n'
              WRITE(unit=0,*) TRIM(msgl)
              tst1 = -polc/polb
              tst2 = tst1
              IF ((tst1 .GE. 0.0_mk) .AND. (tst1 .LE. 1.0_mk)) THEN
                x1 = (x1-tst1*x1+tst1*x2+tst1*vx2-px)/(-vx1+tst1*vx1+tst1*vx2)
                x2 = tst1
              ELSE
                !PROBLEMS
                !JTR clean up in all this - remove all this output
                WRITE(msgl,*) '\nERROR in interpolation variables: rank',rank, &
                  & 'ijk',i,j,k,'xyz',tx,ty,tz,level,'pol abc',pola,polb,polc, &
                  & 'tst1',tst1,'ubar',field_ubar(:,i,j,k,isub),'\n'
                WRITE(unit=0,*) TRIM(msgl)
              ENDIF
            ENDIF
            !LAZARUS HAS BEEN KILLED GO FINISH OFF DIABLO: set ubar from x1 & x2
            !-------------------------------------------------------------------
            ! Interpolate velocity from boundary. ...symmetry around y-axis?
            !-------------------------------------------------------------------
            field_ubar(1,i,j,k,isub) = cosangxz * ( &
               & (1-x2)*((1-x1)*medusa_vel(np1,1) + x1*medusa_vel(np1+1,1)) + &
               & x2*((1-x1)*medusa_vel(np2+1,1)+x1*medusa_vel(np2,1)))
            field_ubar(2,i,j,k,isub) = &
               & (1-x2)*((1-x1)*medusa_vel(np1,2) + x1*medusa_vel(np1+1,2)) + &
               & x2*((1-x1)*medusa_vel(np2+1,2)+x1*medusa_vel(np2,2))
            field_ubar(3,i,j,k,isub) = sinangxz * ( &
               & (1-x2)*((1-x1)*medusa_vel(np1,1) + x1*medusa_vel(np1+1,1)) + &
               & x2*((1-x1)*medusa_vel(np2+1,1)+x1*medusa_vel(np2,1)))
          ELSE
            ! This may be insufficient for wide penalization intervals
            IF (offofpanel) THEN
              field_ubar(1,i,j,k,isub) = cosangxz * (medusa_vel(np1,1)*&
                & (1.0_mk-panel_part) + medusa_vel(np1+1,1)*panel_part)
              field_ubar(2,i,j,k,isub) = (medusa_vel(np1,2)*&
                & (1.0_mk-panel_part) + medusa_vel(np1+1,2)*panel_part)
              field_ubar(3,i,j,k,isub) = sinangxz * (medusa_vel(np1,1)*&
                & (1.0_mk-panel_part) + medusa_vel(np1+1,1)*panel_part)
            ELSE
              field_ubar(1,i,j,k,isub) = cosangxz * medusa_vel(np1,1)
              field_ubar(2,i,j,k,isub) = medusa_vel(np1,2)
              field_ubar(3,i,j,k,isub) = sinangxz * medusa_vel(np1,1)
            END IF 
          ENDIF

          !---------------------------------------
          ! Sign distance according to in/out
          !---------------------------------------
          level = sqrt(level)
          IF (inside) THEN
            level = sign(level,-1.0_mk)
          ENDIF

          IF (step_function .EQ. 0) THEN
            field_H(i,j,k,isub) = -0.5_mk * TANH(level*epsilon) + 0.5_mk
          ELSEIF (step_function .EQ. 1) THEN
            field_H(i,j,k,isub) = stepfunction1(level*epsilon)
          END IF 

        END DO !i
      END DO !j
    END DO !k
  END DO !isub



  ftopo_id = (/2,3,4/)
  t_topoid = (/2,3,4,5/)
  fmesh_id = (/2,3,4,5/)
  !JTR shouldn't I wait with this?
!  CALL ppm_fft_solenoidal(field_ubar,mesh_id,topo_id,t_topoid,fmesh_id, &
!       & ghostsize, info)

  !JTR if I do not use the solenoidal correction it shouldn't be necessary to 
  !ghost ubar
!  maptype = ppm_param_map_init
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_ghost_get
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_push
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_send
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_pop
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)

!JTMP
  IF (rank .EQ. -10) THEN  
    WRITE(filename,'(A,I5.5,A)') 'outline',itime
    OPEN(14,file=filename,iostat=ios,position='append',status='new')
    DO i=1,medusa_count
      WRITE(14,*) medusa_points(i,:), medusa_vel(i,:)
    END DO
    CLOSE(14)

    WRITE(filename,'(A,I3.3,I3.3,A)') 'ubar',rank,itime
    OPEN(14,file=filename,iostat=ios,position='append',status='new')
    DO isub=1,nsublist
      isubl = isublist(isub)
      DO k=ndata(3,isubl),ndata(3,isubl)
        DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
          DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
            tx = min_sub(1,isubl) + REAL(i-1,mk)*dx - tx_center
            ty = min_sub(2,isubl) + REAL(j-1,mk)*dy - ty_center
            tz = min_sub(3,isubl) + REAL(k-1,mk)*dz - tz_center
            WRITE(14,'7(E)') tx,ty,tz,field_ubar(:,i,j,k,isub),field_H(i,j,k,isub)!
          END DO
        END DO
      END DO
    END DO
    CLOSE(14)

  WRITE(msg,*) '\nrank:',rank , min_sub(1,isubl), min_sub(2,isubl),  min_sub(3,isubl)
  WRITE(unit=0,*) msg
  ENDIF


  !JTR this should eventually be deleted
  ubarmax = 0.0_mk
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO k=ndata(3,isubl),ndata(3,isubl)
      DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
        DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
          ubarmax = MAX(ubarmax, SUM(field_ubar(1:3,i,j,k,isub)**2))
        END DO
      END DO
    END DO
  END DO
  CALL MPI_Reduce(ubarmax,tx,1,mpi_prec,MPI_MAX,0,comm,info)
  IF (rank.eq.0) THEN
    WRITE (msg,*) '\nitime=',itime,', ubarmax=',ubarmax,'\n'
    WRITE (unit=0,*) TRIM(msg)
  END IF



END SUBROUTINE wvic_stepfunc_medusa


SUBROUTINE wvic_interpolate_medusa
  USE module_wvic
  IMPLICIT NONE

  CHARACTER(len=256)       :: msg
  REAL(mk), DIMENSION(2)   :: point1,point2,panel
  REAL(mk)                 :: fractiont,panelsq,epsilon
  INTEGER                  :: p,frame1,frame2

  !-----------------------------------------------------------------------------
  ! Interpolate between frames
  !-----------------------------------------------------------------------------
  frame1=MOD(FLOOR(time/medusa_dt),medusa_frames)+1
  IF (frame1 .EQ. medusa_frames) THEN
    frame2 = 1
  ELSE
    frame2 = frame1+1
  ENDIF
  fractiont = MOD(time,medusa_dt)/medusa_dt

  medusa_points = (1.0_mk - fractiont)*medusa_points_array(frame1,:,:) + &
                & fractiont*medusa_points_array(frame2,:,:)
  medusa_vel    = (1.0_mk - fractiont)*medusa_vel_array   (frame1,:,:) + &
                & fractiont*medusa_vel_array   (frame2,:,:) 

  !-----------------------------------------------------------------------------
  ! Add medusa translation velocity and integrated position
  !-----------------------------------------------------------------------------
  DO p=1,medusa_count
    medusa_points(p,2) = medusa_points(p,2) + medu_move_pos
    medusa_vel(p,2)    =    medusa_vel(p,2) + medu_move_vel
  END DO

  !-----------------------------------------------------------------------------
  ! Calculate max/min bounds for medusa, taking into account the width of the
  ! step function and its offset AND the movement of the medusa. Strictly
  ! speaking only half the step function interval needs to be considered
  !-----------------------------------------------------------------------------
  epsilon     = (sqrt(dx**2 + dy**2 + dz**2)*stepfunction_band)
  mmaxx = MAX(medusa_maxx(frame1),medusa_maxx(frame2))+epsilon*(step1_interval &
        & + step1_offset)
  mmaxy = MAX(medusa_maxy(frame1),medusa_maxy(frame2))+epsilon*(step1_interval &
        & + step1_offset) + medu_move_pos 
  mminy = MAX(medusa_miny(frame1),medusa_miny(frame2))-epsilon*(step1_interval &
        & + step1_offset) + medu_move_pos


  !-----------------------------------------------------------------------------
  ! Initialize interpolated frames, i.e. calculate panels, normals and 
  ! their squared lengths
  !-----------------------------------------------------------------------------
  !JTR cleanup arrays. It seems medusa_panel is not used
  DO p=1,medusa_count-1
    point1                  = medusa_points(p,:)
    point2                  = medusa_points(p+1,:)
    panel                   = point2 - point1
    panelsq                 = SUM(panel**2)
    medusa_panel(p,:)       = panel
    medusa_panel_squared(p) = panelsq
    medusa_normal(p,1)      = -panel(2)/sqrt(panelsq)
    medusa_normal(p,2)      = panel(1)/sqrt(panelsq)
  END DO

END SUBROUTINE wvic_interpolate_medusa

!JTR rename to wvic_medusa_integrate_motion ?
SUBROUTINE wvic_medusa_move

  USE module_wvic
  USE ppm_module_rmsh_create_part
  USE ppm_module_interp_p2m
  USE ppm_module_data
  USE ppm_module_map_field_ghost
  USE ppm_module_map_part

  IMPLICIT NONE

  REAL(mk), DIMENSION(3) :: len_physg
  INTEGER                :: i,kp
  INTEGER                :: info
  INTEGER                :: mpart
  INTEGER                :: maptype
  LOGICAL                :: alldone
  CHARACTER(len=256)     :: msg

  len_physg = max_physg - min_physg


!  CALL ppm_rmsh_create_part(xpc,npc,hpc,field_H,topo_id,mesh_id,&
!     & (/-10.0_mk,1.01_mk/),info,resetpos=.TRUE.,field_wp = field_ubar,wp&
!     & = ubarp, lda2 = dime)

!JTR medusa_velocity is not intuitive. Use medusa_cm_velocity or medusa_total_velocity
  medu_move_vel = medu_move_vel + dt*medu_thrusty/medu_mass_orig
  medu_move_pos = medu_move_pos + medu_move_vel*dt
IF (rank.EQ. 0) THEN
WRITE(msg,*) '\nJTR3',rank,medu_move_pos,medu_move_vel,medu_thrusty,medu_mass_orig,'\n' !JTR
WRITE(UNIT=0,*) msg
ENDIF 

  !-----------------------------------------------------------------------------
  ! Move all particles backward 
  ! (This is when not moving medusa but instead moving the flow/solution)
  !-----------------------------------------------------------------------------
  !DO kp=1,np
  !  xp(2,kp) = xp(2,kp) - medu_move_vel*dt
  !END DO

  !-----------------------------------------------------------------------------
  ! Catch particles outside domain
  !-----------------------------------------------------------------------------
!  alldone = .FALSE.
!  DO WHILE (alldone .EQ. .FALSE.)
!  alldone = .TRUE.
!    DO kp=1,np
!      IF(xp(2,kp).GE.max_physg(2)) THEN
!        xp(2,kp) = xp(2,kp) - len_physg(2)
!        alldone   = .FALSE.
!      ELSE
!        IF(xp(2,kp).LT.min_physg(2)) THEN
!          xp(2,kp) = xp(2,kp) + len_physg(2)
!          alldone   = .FALSE.
!        END IF
!      END IF
!    END DO
!  END DO

  !-----------------------------------------------------------------------------
  ! Map the particles
  !-----------------------------------------------------------------------------
!  maptype = ppm_param_map_partial
!  CALL ppm_map_part(xpc,dime,npc,mpart,topo_id,maptype,info)
!  maptype = ppm_param_map_push
!  CALL ppm_map_part(hpc,npc,mpart,topo_id,maptype,info)
!  maptype = ppm_param_map_send
!  CALL ppm_map_part(hpc,npc,mpart,topo_id,maptype,info)
!  maptype = ppm_param_map_pop
!  CALL ppm_map_part(hpc,npc,mpart,topo_id,maptype,info)
!  CALL ppm_map_part(xpc,dime,npc,mpart,topo_id,maptype,info)
!  npc=mpart

  ! Is this nescessary or does interp_p2m set fields to zero?
!  field_H = 0.0_mk
!  field_ubar = 0.0_mk

!  CALL ppm_interp_p2m(xpc,npc,hpc,topo_id,mesh_id, &
!       & ppm_param_rmsh_kernel_mp4,ghostsize,&
!       & field_H,info)
!  CALL ppm_interp_p2m(xpc,npc,ubarp,lda,topo_id,mesh_id,&
!       &              ppm_param_rmsh_kernel_mp4,&
!       &              ghostsize, field_ubar,info)

  !-----------------------------------------------------------------------------
  ! Ghost the fields - Can this be done manually? (cf. wvic_solid_velocity)
  !-----------------------------------------------------------------------------
!  maptype = ppm_param_map_init
!  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_ghost_get
!  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_push
!  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_send
!  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_pop
!  CALL ppm_map_field_ghost(field_H,topo_id,mesh_id, &
!       & ghostsize,maptype,info)

!  maptype = ppm_param_map_init
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_ghost_get
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_push
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_send
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)
!  maptype = ppm_param_map_pop
!  CALL ppm_map_field_ghost(field_ubar,lda,topo_id,mesh_id,&
!       & ghostsize,maptype,info)


END SUBROUTINE wvic_medusa_move

