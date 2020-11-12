!-------------------------------------------------------------------------------
! WVIC_READSTL
! 2007/2008
! Read stl file, import and treat facets
!
! 060608 added reading medusa contour
!
! Johannes Tophoej Rasmussen betonarbejder@gmail.com
!-------------------------------------------------------------------------------

SUBROUTINE wvic_readstl(info)

  USE module_wvic
  USE ppm_module_write
  USE MPI


  INTEGER, INTENT(inout)                :: info

  CHARACTER(len=156) :: msg
  INTEGER, PARAMETER :: md = kind(2.0d0)

  INTEGER                  :: i,tr
  INTEGER                  :: vertex,state
  INTEGER                  :: ilenstl,iline
  INTEGER                  :: iUnit
  INTEGER                  :: ILEN,ios
  CHARACTER(LEN=256)       :: cbuf, tmpbuf
  CHARACTER(LEN=256)       :: cvalue,carg
  LOGICAL                  :: lExist
  REAL(mk), DIMENSION(3)   :: vertex1,vertex2,vertex3,vecu,vecv,vecw

  !-----------------------------------------------------------------------------
  ! Definition of file unit
  !-----------------------------------------------------------------------------
  iUnit = 30

  !-----------------------------------------------------------------------------
  ! open STL file
  !-----------------------------------------------------------------------------
  ilenstl = LEN_TRIM(stlfile)         ! length of filename string
  INQUIRE(FILE=stlfile(1:ilenstl), EXIST=lExist)
  IF(.NOT. lExist) THEN
    WRITE(*,'(2A)')'No such STL file: ',stlfile(1:ilenstl), '\n'
     Info = 1
     call mpi_finalize(info)
     stop
!     GOTO 9999
  END IF

  OPEN(iUnit, FILE=stlfile(1:ilenstl), IOSTAT=ios, ACTION='READ')
  IF(ios .NE. 0) THEN
    WRITE(*,'(2A)')'Failed to open STL file: ',stlfile(1:ilenstl), '\n'
     Info = 1
     GOTO 9999
  END IF

  !-----------------------------------------------------------------------------
  ! Count facets
  !-----------------------------------------------------------------------------
  tri_count = 0
  DO
    READ(iUnit,'(A)',END=50,ERR=200) cbuf
    cbuf=ADJUSTL(cbuf)
    ILEN = LEN_TRIM(cbuf)
    CALL UpperCase(cbuf,5)
    IF('FACET' .EQ. cbuf(1:5)) THEN
      tri_count = tri_count + 1
    END IF
  END DO

50 REWIND(iUnit)
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

  !-----------------------------------------------------------------------------
  ! scan file
  !-----------------------------------------------------------------------------
  state = 0
  iline = 0
  tr=0
  vertex = 1
  DO
    iline = iline + 1        ! increment line
    READ(iUnit,'(A)',END=100,ERR=200) cbuf
    !  IF (iline .EQ. 1) THEN
    !    WRITE(msg,*) '\n rank:',rank,'+',cbuf ,'+',ILEN,'\n'
    !    WRITE(0,*) msg
    !  END IF
    cbuf=ADJUSTL(cbuf)
    ILEN = LEN_TRIM(cbuf)

    !--------------------------------------------------------------------------
    !  Skip comment or empty lines
    !--------------------------------------------------------------------------
    IF(ILEN .GT. 0 .AND. cbuf(1:1) .NE. '#') THEN
      CALL UpperCase(cbuf,ILEN)
      IF('SOLID' .EQ. cbuf(1:5)) THEN

        IF(state .NE. 0) THEN
          WRITE(*,*)'STL input: misplaced SOLID, line ',iline, '\n'
          Info = 1
          GOTO 9999
        ELSE
          state=1
          CYCLE
        END IF
      END IF
      IF('FACET NORMAL' .EQ. cbuf(1:12)) THEN
        IF(state .NE. 1) THEN
          WRITE(*,*)'STL input: misplaced FACET, line ',iline, '\n'
          Info = 1
          GOTO 9999
        ELSE
          tr=tr+1
          tmpbuf=cbuf(13:ILEN)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) tri_norm(tr,:)
          state=2
          CYCLE
        END IF
      END IF
      IF('VERTEX' .EQ. cbuf(1:6)) THEN
        IF(state .NE. 2) THEN
          WRITE(*,*)'STL input: misplaced VERTEX, line ',iline, '\n'
          Info = 1
          GOTO 9999
        ELSE IF (vertex .EQ. 1) THEN
          tmpbuf=cbuf(7:ILEN)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) vertex1
          vertex1 = vertex1 * stl_scale
          vertex1 = vertex1 + stl_translate
          IF (tr .EQ. 1) THEN
            bndminx = vertex1(1)
            bndmaxx = vertex1(1)
            bndminy = vertex1(2)
            bndmaxy = vertex1(2)
            bndminz = vertex1(3)
            bndmaxz = vertex1(3)
          END IF
          bndminx = MIN(bndminx,vertex1(1))
          bndmaxx = MAX(bndmaxx,vertex1(1))
          bndminy = MIN(bndminy,vertex1(2))
          bndmaxy = MAX(bndmaxy,vertex1(2))
          bndminz = MIN(bndminz,vertex1(3))
          bndmaxz = MAX(bndmaxz,vertex1(3))
          vertex=2
          CYCLE
        ELSE IF (vertex .EQ. 2) THEN
          tmpbuf=cbuf(7:ILEN)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) vertex2
          vertex2 = vertex2 * stl_scale
          vertex2 = vertex2 + stl_translate
          bndminx = MIN(bndminx,vertex2(1))
          bndmaxx = MAX(bndmaxx,vertex2(1))
          bndminy = MIN(bndminy,vertex2(2))
          bndmaxy = MAX(bndmaxy,vertex2(2))
          bndminz = MIN(bndminz,vertex2(3))
          bndmaxz = MAX(bndmaxz,vertex2(3))
          vertex=3
          CYCLE
        ELSE IF (vertex .EQ. 3) THEN
          tmpbuf=cbuf(7:ILEN)
          tmpbuf=ADJUSTL(tmpbuf)
          READ(tmpbuf,*,iostat=ios,err=200) vertex3
          vertex3 = vertex3 * stl_scale
          vertex3 = vertex3 + stl_translate
          bndminx = MIN(bndminx,vertex3(1))
          bndmaxx = MAX(bndmaxx,vertex3(1))
          bndminy = MIN(bndminy,vertex3(2))
          bndmaxy = MAX(bndmaxy,vertex3(2))
          bndminz = MIN(bndminz,vertex3(3))
          bndmaxz = MAX(bndmaxz,vertex3(3))
          vertex=1
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

          CYCLE
        END IF
      END IF
      IF('ENDFACET' .EQ. cbuf(1:8)) THEN
        IF(state .NE. 2) THEN
          WRITE(*,*)'STL input: misplaced ENDFACET, line ',iline, '\n'
          Info = 1
          GOTO 9999
        ELSE
          state=1
          CYCLE
        END IF
      END IF
      IF('ENDSOLID' .EQ. cbuf(1:8)) THEN
        IF((state .NE. 1) .AND. (vertex .NE. 1)) THEN
          WRITE(*,*)'STL input: misplaced ENDSOLID, line ',iline, '\n'
          Info = 1
          GOTO 9999
        ELSE
          state=0
          CYCLE
        END IF
      END IF
    END IF

  END DO   




  GOTO 9999



200 CONTINUE
  WRITE(*,'(A,I5,2A)') 'Error reading line: ',iline,                      &
       &                 ' of file: ',stlfile(1:ilenstl)
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
  CLOSE(iUnit)

  !-----------------------------------------------------------------------------
  !  Do last checks
  !-----------------------------------------------------------------------------
  IF (tri_count .NE. tr) THEN
    WRITE(*,*)' Error: Bad triangle count ', tri_count, ' /= ',tr ,'\n'
    Info = 1
    GOTO 9999
  END IF
  IF (tri_count .LT. 4) THEN
    WRITE(*,*)' Error: At least 4 facets required\n'
    Info = 1
    GOTO 9999
  END IF

  IF ( (bndminx .LT. min_physg(1)) .OR. (bndmaxx .GT. max_physg(1)) .OR. &
     & (bndminy .LT. min_physg(2)) .OR. (bndmaxy .GT. max_physg(2)) .OR. &
     & (bndminz .LT. min_physg(3)) .OR. (bndmaxz .GT. max_physg(3)) ) THEN
    IF (stl_check_bounding .EQV. .true.) THEN 
      WRITE(*,*)' Error: STL-object exceeds  computational domain:\n'
      WRITE(*,*)' minx: ', bndminx, ' miny: ', bndminy, ' minz: ', bndminz, '\n'
      WRITE(*,*)' maxx: ', bndmaxx, ' maxy: ', bndmaxy, ' maxz: ', bndmaxz, '\n'
      Info = 1
      GOTO 9999
      ! THIS SHOULD TERMINATE THE PROGRAM
    ELSE
      WRITE(*,*)' Warning: Object exceeds domain\n'
    END IF 
  END IF

  IF (rank .EQ. 0) THEN
    WRITE(msg,*) '\nrank:',rank,'Succesfully read', tri_count, ' facets from ', stlfile(1:ilenstl),'\n'
    WRITE(*,*) msg
    WRITE(msg,*)' minx:', bndminx, 'bndminy:', bndminy, 'minz:', bndminz, '\n',&
              & ' maxx:', bndmaxx, 'bndmaxy:', bndmaxy, 'maxz:', bndmaxz, '\n'
    WRITE(6,*) msg
  END IF

  !-----------------------------------------------------------------------------
  !  Return
  !-----------------------------------------------------------------------------
9999 CONTINUE

END SUBROUTINE wvic_readstl


