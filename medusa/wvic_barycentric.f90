!-------------------------------------------------------------------------------
! WVIC_BARYCENTRIC
! 2008
! Subroutine for initialising the step function \chi_S based on triangulated 
! surface. Calculates intersections between (a line from the point to a 
! reference point) and (the triangles). Calculates barycentric coordinates and
! determines whether the line intersects the triangle. If so the direction of
! interesection is stored and after checking all triangles it can be determined
! if the point is interior or exterior.
!
! Johannes Tophoej Rasmussen betonarbejder@gmail.com
!-------------------------------------------------------------------------------
SUBROUTINE wvic_barycentric

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map_field_ghost
  USE ppm_module_rmsh_create_part

  REAL(MK), EXTERNAL :: stepfunction1
  CHARACTER(len=156) :: msg
  INTEGER, PARAMETER :: md = kind(2.0d0)

  INTEGER                  :: i,j,k,op,op2,tr,isub,isubl,maptype
  INTEGER                  :: s1,s2,s3
  INTEGER                  :: multicheck,lastlevel,inout,intersections,info,ios
  REAL(mk), DIMENSION(3)   :: rdir,vecT0,vecT,vecP,vecS,vece
  REAL(mk), DIMENSION(8,3) :: fixpoints
  REAL(mk)                 :: rdotn,Pdotu,Pdotv,a,b,c
  REAL(mk)                 :: lobjectvolume,gobjectvolume,dv,epsilon
  INTEGER                  :: ibmin,ibmax,jbmin,jbmax,kbmin,kbmax
  INTEGER                  :: imin,imax,jmin,jmax,kmin,kmax
  REAL(mk)                 :: minimumdist,diag
  REAL(mk)                 :: int1,int2,int3,int4,int5,int6
  INCLUDE 'mpif.h'



  fixpoints(1,1) = min_physg(1)
  fixpoints(1,2) = min_physg(2)
  fixpoints(1,3) = min_physg(3)

  fixpoints(2,1) = max_physg(1)
  fixpoints(2,2) = min_physg(2)
  fixpoints(2,3) = max_physg(3)

  fixpoints(3,1) = min_physg(1)
  fixpoints(3,2) = max_physg(2)
  fixpoints(3,3) = max_physg(3)

  fixpoints(4,1) = max_physg(1)
  fixpoints(4,2) = max_physg(2)
  fixpoints(4,3) = min_physg(3)

  fixpoints(5,1) = min_physg(1)
  fixpoints(5,2) = min_physg(2)
  fixpoints(5,3) = max_physg(3)

  fixpoints(6,1) = min_physg(1)
  fixpoints(6,2) = max_physg(2)
  fixpoints(6,3) = min_physg(3)

  fixpoints(7,1) = max_physg(1)
  fixpoints(7,2) = max_physg(2)
  fixpoints(7,3) = max_physg(3)

  fixpoints(8,1) = max_physg(1)
  fixpoints(8,2) = min_physg(2)
  fixpoints(8,3) = min_physg(3)

  field_H = 0.0_mk
  epsilon = 1.0_mk / (sqrt(dx**2 + dy**2 + dz**2))
  diag = (max_physg(1) - min_physg(1))**2 + &
       & (max_physg(2) - min_physg(2))**2 + &
       & (max_physg(3) - min_physg(3))**2
  DO isub=1,nsublist
    isubl = isublist(isub)
    !---------------------------------------------------------------------------
    ! determine boundaries for the cells to be evaluated
    ! the boundaries is expanded by the number of cells equal to 50% of 
    ! step1_interval cell diagonals
    ! ijkbmin will only exceed ijkbmax when a part of the CV is in the subdomain
    ! so that DO loops only initiate when CVs are present in the subdomain:
    !---------------------------------------------------------------------------
    i = ceiling(sqrt(3.0_mk)*(step1_interval*0.5_mk + &
      & step1_offset))
    ibmin=FLOOR( (bndminx-min_sub(1,isubl))/dx+0.5_mk)-1-i
    jbmin=FLOOR( (bndminy-min_sub(2,isubl))/dy+0.5_mk)-1-i
    kbmin=FLOOR( (bndminz-min_sub(3,isubl))/dz+0.5_mk)-1-i
    ibmax=FLOOR( (bndmaxx-min_sub(1,isubl))/dx+0.5_mk)+1+i
    jbmax=FLOOR( (bndmaxy-min_sub(2,isubl))/dy+0.5_mk)+1+i
    kbmax=FLOOR( (bndmaxz-min_sub(3,isubl))/dz+0.5_mk)+1+i
    imin=max(ibmin,1)
    jmin=max(jbmin,1)
    kmin=max(kbmin,1)
    imax=min(ibmax,ndata(1,isubl))
    jmax=min(jbmax,ndata(2,isubl))
    kmax=min(kbmax,ndata(3,isubl))
    DO k=kmin,kmax 
      DO j=jmin,jmax 
	DO i=imin,imax 
	  vecT(1) = min_sub(1,isubl) + REAL(i-1,mk)*dx
	  vecT(2) = min_sub(2,isubl) + REAL(j-1,mk)*dy
	  vecT(3) = min_sub(3,isubl) + REAL(k-1,mk)*dz

	  multicheck=0

	  DOfixpoints: DO op=1,8
	    inout=0
	    intersections=0
	    
	    rdir = fixpoints(op,:) - vecT
            minimumdist = diag
	    DO tr=1,tri_count
	      rdotn = (tri_norm(tr,1)*rdir(1) + tri_norm(tr,2)*rdir(2) + &
		    & tri_norm(tr,3)*rdir(3))
	      vecS=tri_base(tr,:)-vecT

              !----------------------------------------------------------------
              ! meassure distance
              !----------------------------------------------------------------
              vecP = tri_norm(tr,:) * &
                   & ( tri_norm(tr,1)*vecS(1) + tri_norm(tr,2)*vecS(2) + &
                   & tri_norm(tr,3)*vecS(3)) - vecS

	      Pdotu = vecP(1)*tri_vecu(tr,1) + vecP(2)*tri_vecu(tr,2) + &
		    & vecP(3)*tri_vecu(tr,3)
	      Pdotv = vecP(1)*tri_vecv(tr,1) + vecP(2)*tri_vecv(tr,2) + &
		    & vecP(3)*tri_vecv(tr,3)

	      a = (Pdotu*tri_vdotv(tr) - Pdotv*tri_udotv(tr)) * tri_denom(tr)
	      b = (Pdotv*tri_udotu(tr) - Pdotu*tri_udotv(tr)) * tri_denom(tr)

	      IF ((a .GE. 0.0_mk) .AND. (b .GE. 0.0_mk) .AND. & 
                                                  & ((a+b) .LE. 1.0_mk)) THEN
                !get distance normal to triangle
                 a = (tri_norm(tr,1)*vecS(1) + tri_norm(tr,2)*vecS(2) + &
                   & tri_norm(tr,3)*vecS(3))**2
              ELSE
                !get distance to edge or corner
                a = vecP(1)*tri_vecu(tr,1) + vecP(2)*tri_vecu(tr,2) + &
                  & vecP(3)*tri_vecu(tr,3)
                b = vecP(1)*tri_vecv(tr,1) + vecP(2)*tri_vecv(tr,2) + &
                  & vecP(3)*tri_vecv(tr,3)
                c = (vecP(1)-tri_vecu(tr,1))*tri_vecw(tr,1) + &
                  & (vecP(2)-tri_vecu(tr,2))*tri_vecw(tr,2) + &
                  & (vecP(3)-tri_vecu(tr,3))*tri_vecw(tr,3)
  
                IF (a .GT. tri_udotu(tr)) THEN
                   !a = tri_udotu(tr)
                   vece=-vecS-tri_vecu(tr,:)
                   a = SUM(vece**2)
                ELSEIF (a .LT. 0.0_mk) THEN
                   !a = 0
                   !vece=-vecS
                   a = SUM(vecS**2)
                ELSE
                   vece=-vecS-tri_vecu(tr,:)*a/tri_udotu(tr)
                   a = SUM(vece**2)
                   !this may not be the cheapest way to do this
                ENDIF 

                IF (b .GT. tri_vdotv(tr)) THEN
                   !b = tri_vdotv(tr)
                   vece=-vecS-tri_vecv(tr,:)
                   b = SUM(vece**2)
                ELSEIF (b .LT. 0.0_mk) THEN
                   !b = 0
                   !vece=-vecS
                   b = SUM(vecS**2)
                ELSE
                   vece=-vecS-tri_vecv(tr,:)*b/tri_vdotv(tr)
                   b = SUM(vece**2)
                ENDIF 
  
                IF (c .GT. tri_wdotw(tr)) THEN
                   !c = tri_wdotw(tr)
                   vece=-vecS-tri_vecv(tr,:)
                   c = SUM(vece**2)
                ELSEIF (c .LT. 0.0_mk) THEN
                   !c = 0
                   vece=-vecS-tri_vecu(tr,:)
                   c = SUM(vecS**2)
                ELSE
                   vece=-vecS-tri_vecu(tr,:)-tri_vecw(tr,:)*c/tri_wdotw(tr)
                   c = SUM(vece**2)
                ENDIF

                a = min(a,b)
                a = min(a,c)

              ENDIF

              minimumdist = min(minimumdist,a)

              IF (minimumdist .EQ. 0.0_mk) THEN
                EXIT dofixpoints 
              END IF

	      !-----------------------------------------------------------------
	      ! Checking if the triangle is aligned with the line of sight
	      ! If so, loop through the remaining fixpoints until a LOS that is
	      ! not parallel to the triangle is found, calculate P and
	      ! restore rdir. Otherwise calculate P
	      !-----------------------------------------------------------------
	      IF (rdotn .EQ. 0.0_mk) THEN
		DO op2=(op+1),8
		  rdir = fixpoints(op2,:) - vecT
		  rdotn = (tri_norm(tr,1)*rdir(1) + tri_norm(tr,2)*rdir(2) + &
			& tri_norm(tr,3)*rdir(3))
		  IF (rdotn .NE. 0.0_mk) THEN
		    EXIT
		  END IF
		END DO
		c = (tri_norm(tr,1)*vecS(1) + tri_norm(tr,2)*vecS(2) + &
		  & tri_norm(tr,3)*vecS(3))/rdotn
                IF (c .LT. 0.0_mk) THEN
                  CYCLE 
                END IF
		vecP = rdir * c - vecS
		rdir = fixpoints(op,:) - vecT
	      ELSE
		c = (tri_norm(tr,1)*vecS(1) + tri_norm(tr,2)*vecS(2) + &
		  & tri_norm(tr,3)*vecS(3))/rdotn
                IF (c .LT. 0.0_mk) THEN
                  CYCLE
                END IF
		vecP = rdir * c - vecS
	      END IF

	      Pdotu = vecP(1)*tri_vecu(tr,1) + vecP(2)*tri_vecu(tr,2) + &
		    & vecP(3)*tri_vecu(tr,3)
	      Pdotv = vecP(1)*tri_vecv(tr,1) + vecP(2)*tri_vecv(tr,2) + &
		    & vecP(3)*tri_vecv(tr,3)

	      a = (Pdotu*tri_vdotv(tr) - Pdotv*tri_udotv(tr)) * tri_denom(tr)
	      b = (Pdotv*tri_udotu(tr) - Pdotu*tri_udotv(tr)) * tri_denom(tr)

  
	      IF (a .LT. 0.0_mk) THEN
		CYCLE
	      END IF
	      IF (b .LT. 0.0_mk) THEN
		CYCLE
	      END IF
  
	      IF ((a+b) .GT. 1.0_mk) THEN
		CYCLE
	      END IF

	      !-----------------------------------------------------------------
	      ! The line is intersecting the triangle. Update intersections
	      !-----------------------------------------------------------------
	      intersections=intersections+1

	      IF (rdotn .GT. 0.0_mk) THEN
		inout=inout+1;
	      ELSE
		inout=inout-1;
	      END IF

	    END DO !tr triangles

	    !-------------------------------------------------------------------
	    ! Test number of intersections against intersection directions
	    !-------------------------------------------------------------------
	    IF (stl_check_intersections .EQV. .true.) THEN
            IF (((MOD(intersections,2) .EQ. 0) .AND. (inout .NE. 0))) THEN
	      !diverging normals at intersections
              IF (stl_nonverbose .NEQV. .TRUE.) THEN
		WRITE(msg,*) '\nrank:',rank, 'ijk:', &
                & i,j,k, 'normals diverge (a), intersections:',intersections,&
                & 'inout:',inout
                WRITE(6,*) msg
              ENDIF
	      CYCLE
	    END IF
	    IF (MOD(intersections,2) .NE. 0) THEN
	      IF ((stl_min_phys_inside .EQV. .true.) .AND. (inout .NE. -1)) THEN
                IF (stl_nonverbose .NEQV. .TRUE.) THEN
		  WRITE(msg,*) '\nrank:',rank, 'ijk:', &
                  & i,j,k, 'normals diverge (b), intersections:',intersections,&
                  & 'inout:',inout
		  WRITE(6,*) msg
                ENDIF
	        CYCLE
	      END IF
	      IF ((stl_min_phys_inside .EQV. .false.) .AND. (inout .NE. 1)) THEN
                IF (stl_nonverbose .NEQV. .TRUE.) THEN
		  WRITE(msg,*) '\nrank:',rank, 'ijk:', &
                  & i,j,k, 'normals diverge (c), intersections:',intersections,&
                  & 'inout:',inout
                  WRITE(6,*) msg
                ENDIF
	        CYCLE
	      END IF
	    END IF
            END IF

	    !-------------------------------------------------------------------
	    ! Double check in/out status with an additional reference point
	    !-------------------------------------------------------------------
            IF (stl_double_check .EQV. .true.) THEN
              IF ((multicheck .EQ. 0) .AND. (op .NE. 8)) THEN
                lastlevel = inout
                multicheck = 1
                CYCLE
              ELSE IF (lastlevel .EQ. inout) THEN
                multicheck = 2
                EXIT
              ELSE IF (op .EQ. 8) THEN
                WRITE(msg,*) '\nrank:',rank,'unable to achive multicheck, ijk',&
                           & i,j,k
                WRITE(0,*) msg 
                call mpi_finalize(info)
                stop
                call wvic_died
                EXIT
              ELSE
                multicheck = 0
                IF (stl_nonverbose .NEQV. .TRUE.) THEN
                  WRITE(msg,*) '\nrank:',rank, &
                             & ' multicheck failed, recalculate, ijk',i,j,k
                  WRITE(6,*) msg 
                ENDIF
                CYCLE
              END IF
            END IF
            EXIT
          END DO dofixpoints!op

	  !------------------------------------------------------------------
	  ! Update field_H. The intersections have been tested or passed
	  !------------------------------------------------------------------
          minimumdist = sqrt(minimumdist)
	  IF ((inout .GT. 0) .AND. (stl_min_phys_inside .EQV. .false.)) THEN
            minimumdist = sign(minimumdist,-1.0_mk)
	  ELSEIF ((inout .GE. 0) .AND. (stl_min_phys_inside .EQV. .true.)) THEN
            minimumdist = sign(minimumdist,-1.0_mk)
	  END IF

          field_H(i,j,k,isub) = stepfunction1(minimumdist*epsilon)

	END DO !i
      END DO !j
    END DO !k

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

    dv=dx*dy*dz
    lobjectvolume=0.0_mk
    DO k=1,ndata(3,isubl)-1
      DO j=1,ndata(2,isubl)-1
	DO i=1,ndata(1,isubl)-1 
          lobjectvolume = lobjectvolume + field_H(i,j,k,isub)
	END DO !i
      END DO !j
    END DO !k
    lobjectvolume = lobjectvolume * dv

  END DO !isub

  CALL MPI_Reduce(lobjectvolume,gobjectvolume,1,mpi_prec,MPI_SUM,0,comm,info)
  IF (rank .EQ. 0) THEN
    WRITE(msg,*) '\nSTL object, volume ' ,gobjectvolume
    WRITE(0,*) msg
  END IF


END SUBROUTINE wvic_barycentric

