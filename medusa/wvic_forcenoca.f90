!-------------------------------------------------------------------------------
! WVIC_FORCENOCA
! 2007/2008
! Subroutine for measuring forces in arbitrary control volumes based on the
! "impulse equation" expression from F. Noca et al.,
! "A comparison of methods for evaluating time-dependent fluid dynamic forces
! on bodies using only velocity fields and their derivatives."
!
! Johannes Tophoej Rasmussen betonarbejder@gmail.com
!-------------------------------------------------------------------------------

SUBROUTINE wvic_forcenoca (tforcenoca,gnocadimensions,info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE



  INTEGER , INTENT(inout) :: info
  CHARACTER(len=512) :: msg
  INTEGER, PARAMETER :: md = kind(2.0d0)

  INTEGER  :: i,j,k,ni,isub,isubl
  REAL(mk) :: tx, ty, tz
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9
  REAL(mk) :: dv
  INTEGER, DIMENSION(3) :: ncimin, ncimax, cncimin, cncimax
  REAL(mk) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, dsu, dsv, dsw
  REAL(mk) :: duxy, dvxy, duxz, dwxz, dvyz, dwyz
  REAL(mk) :: duxx, duyy, duzz, dvxx, dvyy, dvzz, dwxx, dwyy, dwzz
  !-----------------------------------------------------------------------------
  !  sforcenoca/gsforcenoca - surface contribution
  !  vforcenoca/gvforcenoca - volume contribution (vorticity moments)
  !  tforcenoca - total noca force
  !  (g)forcenoca(W,E,S,N,B,T) - surface contribution from individual surfaces
  !  (l,s)nocadimensions - diagnostics: no. cells, volume and surface area
  !-----------------------------------------------------------------------------
  REAL(mk), DIMENSION(3*nforces) :: sforcenoca, tforcenoca, gsforcenoca
  REAL(md), DIMENSION(3*nforces) :: vforcenoca, gvforcenoca
  REAL(md), DIMENSION(3*nforces) :: sforcenocaW, gsforcenocaW, &
       sforcenocaE, gsforcenocaE, sforcenocaS, &
       gsforcenocaS, sforcenocaN, gsforcenocaN, &
       sforcenocaB, gsforcenocaB, sforcenocaT, gsforcenocaT
  REAL(mk), DIMENSION(3*nforces) :: lnocadimensions,gnocadimensions
  REAL(mk) :: o1, o2, o3, u1, u2, u3

  INCLUDE 'mpif.h'


  !-----------------------------------------------------------------------------
  !  force a la Noca (Impulse Equation)
  !-----------------------------------------------------------------------------
  fac1 = 0.5_mk / dx
  fac2 = 0.5_mk / dy
  fac3 = 0.5_mk / dz
  fac4 = 1.0_mk / dx**2
  fac5 = 1.0_mk / dy**2
  fac6 = 1.0_mk / dz**2
  fac7 = dy*dz
  fac8 = dx*dz
  fac9 = dx*dy
  dv = dx*dy*dz
  vforcenoca = 0.0_mk
  sforcenoca = 0.0_mk
  sforcenocaW = 0.0_mk
  sforcenocaE = 0.0_mk
  sforcenocaS = 0.0_mk
  sforcenocaN = 0.0_mk
  sforcenocaB = 0.0_mk
  sforcenocaT = 0.0_mk
  lnocadimensions = 0.0_mk
  DO ni=1,nforces
  DO isub=1,nsublist
    isubl = isublist(isub)
    !---------------------------------------------------------------------------
    ! Determining local loop indexes. Check for indexes outside subdomain
    ! nci  - noca (counting) index
    ! cnci - checked noca index
    !---------------------------------------------------------------------------
    ncimin(1)=FLOOR( (forcecv(1+6*(ni-1))-min_sub(1,isubl)) / dx + 0.5_mk)+1
    ncimin(2)=FLOOR( (forcecv(2+6*(ni-1))-min_sub(2,isubl)) / dy + 0.5_mk)+1
    ncimin(3)=FLOOR( (forcecv(3+6*(ni-1))-min_sub(3,isubl)) / dz + 0.5_mk)+1
    ncimax(1)=FLOOR( (forcecv(4+6*(ni-1))-min_sub(1,isubl)) / dx + 0.5_mk)
    ncimax(2)=FLOOR( (forcecv(5+6*(ni-1))-min_sub(2,isubl)) / dy + 0.5_mk)
    ncimax(3)=FLOOR( (forcecv(6+6*(ni-1))-min_sub(3,isubl)) / dz + 0.5_mk)
    ! ncimax will only exceed ncimin when a part of the CV is in the subdomain
    ! so that DO loops only initiate when CVs are present in the subdomain:
    cncimin(1)=max(ncimin(1),1)
    cncimin(2)=max(ncimin(2),1)
    cncimin(3)=max(ncimin(3),1)
    cncimax(1)=min(ncimax(1),ndata(1,isubl)-1) 
    cncimax(2)=min(ncimax(2),ndata(2,isubl)-1)
    cncimax(3)=min(ncimax(3),ndata(3,isubl)-1)
    !---------------------------------------------------------------------------
    !  Noca volume integral 
    !---------------------------------------------------------------------------
    DO k=cncimin(3),cncimax(3)
      DO j=cncimin(2),cncimax(2)
        DO i=cncimin(1),cncimax(1)
          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

          lnocadimensions(1+3*(ni-1)) = lnocadimensions(1+3*(ni-1)) + dv
          lnocadimensions(2+3*(ni-1)) = lnocadimensions(2+3*(ni-1)) + 1.0_mk

          vforcenoca(1+3*(ni-1)) = vforcenoca(1+3*(ni-1)) &
                & - 0.5_mk*dv*(ty*field_wp(3,i,j,k,isub) &
                & - tz*field_wp(2,i,j,k,isub))
          vforcenoca(2+3*(ni-1)) = vforcenoca(2+3*(ni-1)) &
                & - 0.5_mk*dv*(tz*field_wp(1,i,j,k,isub) &
                & - tx*field_wp(3,i,j,k,isub))
          vforcenoca(3+3*(ni-1)) = vforcenoca(3+3*(ni-1)) &
                & - 0.5_mk*dv*(tx*field_wp(2,i,j,k,isub) &
                & - ty*field_wp(1,i,j,k,isub))

        END DO
      END DO
    END DO
    !---------------------------------------------------------------------------
    !  Noca gamma_imp surface integral: W 
    !---------------------------------------------------------------------------
  IF ((ncimin(1) .GE. 1) .and. (ncimin(1) .LE. (ndata(1,isubl)-1))) THEN
    DO k=cncimin(3),cncimax(3)
      DO j=cncimin(2),cncimax(2)
        i=ncimin(1)
        lnocadimensions(3+3*(ni-1)) = lnocadimensions(3+3*(ni-1)) + dy*dz

        o1 = field_wp(1,i,j,k,isub)
        o2 = field_wp(2,i,j,k,isub)
        o3 = field_wp(3,i,j,k,isub)

        u1 = field_up(1,i,j,k,isub)
        u2 = field_up(2,i,j,k,isub)
        u3 = field_up(3,i,j,k,isub)

        tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
        ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
        tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

        dux = fac1*(field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub))
        dvx = fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))
        dwx = fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
        
        duy = fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
!        dvy = fac2*(field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub))
!        dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))

        duz = fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))
!        dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
!        dwz = fac3*(field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub))

        duxy = fac1*fac2*(field_up(1,i+1,j+1,k,isub)-field_up(1,i+1,j-1,k,isub)&
             & -field_up(1,i-1,j+1,k,isub)+field_up(1,i-1,j-1,k,isub))
        dvxy = fac1*fac2*(field_up(2,i+1,j+1,k,isub)-field_up(2,i+1,j-1,k,isub)&
             & -field_up(2,i-1,j+1,k,isub)+field_up(2,i-1,j-1,k,isub))
!n/u    dwxy = fac1*fac2*(field_up(3,i+1,j+1,k,isub)-field_up(3,i+1,j-1,k,isub)&
!            & -field_up(3,i-1,j+1,k,isub)+field_up(3,i-1,j-1,k,isub))

        duxz = fac1*fac3*(field_up(1,i+1,j,k+1,isub)-field_up(1,i+1,j,k-1,isub)&
             & -field_up(1,i-1,j,k+1,isub)+field_up(1,i-1,j,k-1,isub))
!n/u    dvxz = fac1*fac3*(field_up(2,i+1,j,k+1,isub)-field_up(2,i+1,j,k-1,isub)&
!            & -field_up(2,i-1,j,k+1,isub)+field_up(2,i-1,j,k-1,isub))
        dwxz = fac1*fac3*(field_up(3,i+1,j,k+1,isub)-field_up(3,i+1,j,k-1,isub)&
             & -field_up(3,i-1,j,k+1,isub)+field_up(3,i-1,j,k-1,isub))

!n/u    duyz = fac2*fac3*(field_up(1,i,j+1,k+1,isub)-field_up(1,i,j+1,k-1,isub)&
!            & -field_up(1,i,j-1,k+1,isub)+field_up(1,i,j-1,k-1,isub))
        dvyz = fac2*fac3*(field_up(2,i,j+1,k+1,isub)-field_up(2,i,j+1,k-1,isub)&
             & -field_up(2,i,j-1,k+1,isub)+field_up(2,i,j-1,k-1,isub))
        dwyz = fac2*fac3*(field_up(3,i,j+1,k+1,isub)-field_up(3,i,j+1,k-1,isub)&
             & -field_up(3,i,j-1,k+1,isub)+field_up(3,i,j-1,k-1,isub))

        duxx = fac4*(field_up(1,i+1,j,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i-1,j,k,isub))
        duyy = fac5*(field_up(1,i,j+1,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j-1,k,isub))
        duzz = fac6*(field_up(1,i,j,k+1,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j,k-1,isub))

        dvxx = fac4*(field_up(2,i+1,j,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i-1,j,k,isub))
        dvyy = fac5*(field_up(2,i,j+1,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j-1,k,isub))
        dvzz = fac6*(field_up(2,i,j,k+1,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j,k-1,isub))

        dwxx = fac4*(field_up(3,i+1,j,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i-1,j,k,isub))
        dwyy = fac5*(field_up(3,i,j+1,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j-1,k,isub))
        dwzz = fac6*(field_up(3,i,j,k+1,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j,k-1,isub))

        sforcenocaW(1+3*(ni-1)) = sforcenocaW(1+3*(ni-1)) - fac7 * (&
          & -0.5D0 * u1 ** 2 + 0.5D0 * u2 ** 2 + 0.5D0 * u3 ** 2 - 0.5D0 * u1 &
          & * (ty * o3 - tz * o2) + 0.5D0 * o1 * (ty * u3 - tz * u2) + 0.5D0 &
          & * nu * (dvxx + duxy + 2.0_mk * dvyy + dwyz + dvzz) * ty + 0.5D0 * &
          & nu * (dwxx + duxz + dwyy + dvyz + 2.0_mk * dwzz) * tz + 0.2D1 * &
          & nu * dux)
        sforcenocaW(2+3*(ni-1)) = sforcenocaW(2+3*(ni-1)) - fac7 * (&
          & -u1 * u2 - 0.5D0 * u1 * (tz * o1 - tx * o3) + 0.5D0 * o1 * (tz * &
          & u1 - tx * u3) - 0.5D0 * tx * nu * (dvxx + duxy + 2.0_mk * dvyy + &
          & dwyz + dvzz) + nu * (dvx + duy))
        sforcenocaW(3+3*(ni-1)) = sforcenocaW(3+3*(ni-1)) - fac7 * (&
          & -u1 * u3 - 0.5D0 * u1 * (tx * o2 - ty * o1) + 0.5D0 * o1 * (tx * &
          & u2 - ty * u1) - 0.5D0 * tx * nu * (dwxx + duxz + dwyy + dvyz + &
          & 2.0_mk * dwzz) + nu * (dwx + duz))


      END DO
    END DO
  END IF
  !---------------------------------------------------------------------------
  !  Noca gamma_imp surface integral: E
  !---------------------------------------------------------------------------
  IF ((ncimax(1) .GE. 1) .and. (ncimax(1) .LE. (ndata(1,isubl)-1))) THEN
    DO k=cncimin(3),cncimax(3)
      DO j=cncimin(2),cncimax(2)
        i=ncimax(1)
        lnocadimensions(3+3*(ni-1)) = lnocadimensions(3+3*(ni-1)) + dy*dz

        o1 = field_wp(1,i,j,k,isub)
        o2 = field_wp(2,i,j,k,isub)
        o3 = field_wp(3,i,j,k,isub)

        u1 = field_up(1,i,j,k,isub)
        u2 = field_up(2,i,j,k,isub)
        u3 = field_up(3,i,j,k,isub)

        tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
        ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
        tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

        dux = fac1*(field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub))
        dvx = fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))
        dwx = fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
        
        duy = fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))

        duz = fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))

        duxy = fac1*fac2*(field_up(1,i+1,j+1,k,isub)-field_up(1,i+1,j-1,k,isub)&
             & -field_up(1,i-1,j+1,k,isub)+field_up(1,i-1,j-1,k,isub))
        dvxy = fac1*fac2*(field_up(2,i+1,j+1,k,isub)-field_up(2,i+1,j-1,k,isub)&
             & -field_up(2,i-1,j+1,k,isub)+field_up(2,i-1,j-1,k,isub))

        duxz = fac1*fac3*(field_up(1,i+1,j,k+1,isub)-field_up(1,i+1,j,k-1,isub)&
             & -field_up(1,i-1,j,k+1,isub)+field_up(1,i-1,j,k-1,isub))
        dwxz = fac1*fac3*(field_up(3,i+1,j,k+1,isub)-field_up(3,i+1,j,k-1,isub)&
             & -field_up(3,i-1,j,k+1,isub)+field_up(3,i-1,j,k-1,isub))

        dvyz = fac2*fac3*(field_up(2,i,j+1,k+1,isub)-field_up(2,i,j+1,k-1,isub)&
             & -field_up(2,i,j-1,k+1,isub)+field_up(2,i,j-1,k-1,isub))
        dwyz = fac2*fac3*(field_up(3,i,j+1,k+1,isub)-field_up(3,i,j+1,k-1,isub)&
             & -field_up(3,i,j-1,k+1,isub)+field_up(3,i,j-1,k-1,isub))

        duxx = fac4*(field_up(1,i+1,j,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i-1,j,k,isub))
        duyy = fac5*(field_up(1,i,j+1,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j-1,k,isub))
        duzz = fac6*(field_up(1,i,j,k+1,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j,k-1,isub))

        dvxx = fac4*(field_up(2,i+1,j,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i-1,j,k,isub))
        dvyy = fac5*(field_up(2,i,j+1,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j-1,k,isub))
        dvzz = fac6*(field_up(2,i,j,k+1,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j,k-1,isub))

        dwxx = fac4*(field_up(3,i+1,j,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i-1,j,k,isub))
        dwyy = fac5*(field_up(3,i,j+1,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j-1,k,isub))
        dwzz = fac6*(field_up(3,i,j,k+1,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j,k-1,isub))

        sforcenocaE(1+3*(ni-1)) = sforcenocaE(1+3*(ni-1)) + fac7 * (&
          & -0.5D0 * u1 ** 2 + 0.5D0 * u2 ** 2 + 0.5D0 * u3 ** 2 - 0.5D0 * u1 &
          & * (ty * o3 - tz * o2) + 0.5D0 * o1 * (ty * u3 - tz * u2) + 0.5D0 &
          & * nu * (dvxx + duxy + 2.0_mk * dvyy + dwyz + dvzz) * ty + 0.5D0 * &
          & nu * (dwxx + duxz + dwyy + dvyz + 2.0_mk * dwzz) * tz + 0.2D1 * &
          & nu * dux)
        sforcenocaE(2+3*(ni-1)) = sforcenocaE(2+3*(ni-1)) + fac7 * (&
          & -u1 * u2 - 0.5D0 * u1 * (tz * o1 - tx * o3) + 0.5D0 * o1 * (tz * &
          & u1 - tx * u3) - 0.5D0 * tx * nu * (dvxx + duxy + 2.0_mk * dvyy + &
          & dwyz + dvzz) + nu * (dvx + duy))
        sforcenocaE(3+3*(ni-1)) = sforcenocaE(3+3*(ni-1)) + fac7 * (&
          & -u1 * u3 - 0.5D0 * u1 * (tx * o2 - ty * o1) + 0.5D0 * o1 * (tx * &
          & u2 - ty * u1) - 0.5D0 * tx * nu * (dwxx + duxz + dwyy + dvyz + &
          & 2.0_mk * dwzz) + nu * (dwx + duz))

      END DO
    END DO
  END IF
  !---------------------------------------------------------------------------
  !  Noca gamma_imp surface integral: S
  !---------------------------------------------------------------------------
  IF ((ncimin(2) .GE. 1) .and. (ncimin(2) .LE. (ndata(2,isubl)-1))) THEN
    DO k=cncimin(3),cncimax(3)
      DO i=cncimin(1),cncimax(1)
        j=ncimin(2)
        lnocadimensions(3+3*(ni-1)) = lnocadimensions(3+3*(ni-1)) + dx*dz

        o1 = field_wp(1,i,j,k,isub)
        o2 = field_wp(2,i,j,k,isub)
        o3 = field_wp(3,i,j,k,isub)

        u1 = field_up(1,i,j,k,isub)
        u2 = field_up(2,i,j,k,isub)
        u3 = field_up(3,i,j,k,isub)

        tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
        ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
        tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

        dvx = fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))
        
        duy = fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
        dvy = fac2*(field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub))
        dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))

        dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))

        duxy = fac1*fac2*(field_up(1,i+1,j+1,k,isub)-field_up(1,i+1,j-1,k,isub)&
             & -field_up(1,i-1,j+1,k,isub)+field_up(1,i-1,j-1,k,isub))
        dvxy = fac1*fac2*(field_up(2,i+1,j+1,k,isub)-field_up(2,i+1,j-1,k,isub)&
             & -field_up(2,i-1,j+1,k,isub)+field_up(2,i-1,j-1,k,isub))

        duxz = fac1*fac3*(field_up(1,i+1,j,k+1,isub)-field_up(1,i+1,j,k-1,isub)&
             & -field_up(1,i-1,j,k+1,isub)+field_up(1,i-1,j,k-1,isub))
        dwxz = fac1*fac3*(field_up(3,i+1,j,k+1,isub)-field_up(3,i+1,j,k-1,isub)&
             & -field_up(3,i-1,j,k+1,isub)+field_up(3,i-1,j,k-1,isub))

        dvyz = fac2*fac3*(field_up(2,i,j+1,k+1,isub)-field_up(2,i,j+1,k-1,isub)&
             & -field_up(2,i,j-1,k+1,isub)+field_up(2,i,j-1,k-1,isub))
        dwyz = fac2*fac3*(field_up(3,i,j+1,k+1,isub)-field_up(3,i,j+1,k-1,isub)&
             & -field_up(3,i,j-1,k+1,isub)+field_up(3,i,j-1,k-1,isub))

        duxx = fac4*(field_up(1,i+1,j,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i-1,j,k,isub))
        duyy = fac5*(field_up(1,i,j+1,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j-1,k,isub))
        duzz = fac6*(field_up(1,i,j,k+1,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j,k-1,isub))

        dvxx = fac4*(field_up(2,i+1,j,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i-1,j,k,isub))
        dvyy = fac5*(field_up(2,i,j+1,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j-1,k,isub))
        dvzz = fac6*(field_up(2,i,j,k+1,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j,k-1,isub))

        dwxx = fac4*(field_up(3,i+1,j,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i-1,j,k,isub))
        dwyy = fac5*(field_up(3,i,j+1,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j-1,k,isub))
        dwzz = fac6*(field_up(3,i,j,k+1,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j,k-1,isub))

        sforcenocaS(1+3*(ni-1)) = sforcenocaS(1+3*(ni-1)) - fac8 * (&
          & -u1 * u2 - 0.5D0 * u2 * (ty * o3 - tz * o2) + 0.5D0 * o2 * (ty * &
          & u3 - tz * u2) - 0.5D0 * ty * nu * (2.0_mk * duxx + dvxy + &
          & duyy + dwxz + duzz) + nu * (dvx + duy))
        sforcenocaS(2+3*(ni-1)) = sforcenocaS(2+3*(ni-1)) - fac8 * (&
          & 0.5D0 * u1 ** 2 - 0.5D0 * u2 ** 2 + 0.5D0 * u3 ** 2 - 0.5D0 * u2 &
          & * (tz * o1 - tx * o3) + 0.5D0 * o2 * (tz * u1 - tx * u3) + 0.5D0 &
          & * nu * (2.0_mk * duxx + dvxy + duyy + dwxz + duzz) * tx + &
          & 0.5D0 * nu * (dwxx + duxz + dwyy + dvyz + 2.0_mk * dwzz) * tz &
          & + 0.2D1 * nu * dvy)
        sforcenocaS(3+3*(ni-1)) = sforcenocaS(3+3*(ni-1)) - fac8 * (&
          & -u2 * u3 - 0.5D0 * u2 * (tx * o2 - ty * o1) + 0.5D0 * o2 * (tx * &
          & u2 - ty * u1) - 0.5D0 * ty * nu * (dwxx + duxz + dwyy + dvyz &
          & + 2.0_mk * dwzz) + nu * (dwy + dvz))

      END DO
    END DO
  END IF
  !---------------------------------------------------------------------------
  !  Noca gamma_imp surface integral: N
  !---------------------------------------------------------------------------
  IF ((ncimax(2) .GE. 1) .and. (ncimax(2) .LE. (ndata(2,isubl)-1))) THEN
    DO k=cncimin(3),cncimax(3)
      DO i=cncimin(1),cncimax(1)
        j=ncimax(2)
        lnocadimensions(3+3*(ni-1)) = lnocadimensions(3+3*(ni-1)) + dx*dz

        o1 = field_wp(1,i,j,k,isub)
        o2 = field_wp(2,i,j,k,isub)
        o3 = field_wp(3,i,j,k,isub)

        u1 = field_up(1,i,j,k,isub)
        u2 = field_up(2,i,j,k,isub)
        u3 = field_up(3,i,j,k,isub)

        tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
        ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
        tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

        dvx = fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))
        
        duy = fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
        dvy = fac2*(field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub))
        dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))

        dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))

        duxy = fac1*fac2*(field_up(1,i+1,j+1,k,isub)-field_up(1,i+1,j-1,k,isub)&
             & -field_up(1,i-1,j+1,k,isub)+field_up(1,i-1,j-1,k,isub))
        dvxy = fac1*fac2*(field_up(2,i+1,j+1,k,isub)-field_up(2,i+1,j-1,k,isub)&
             & -field_up(2,i-1,j+1,k,isub)+field_up(2,i-1,j-1,k,isub))

        duxz = fac1*fac3*(field_up(1,i+1,j,k+1,isub)-field_up(1,i+1,j,k-1,isub)&
             & -field_up(1,i-1,j,k+1,isub)+field_up(1,i-1,j,k-1,isub))
        dwxz = fac1*fac3*(field_up(3,i+1,j,k+1,isub)-field_up(3,i+1,j,k-1,isub)&
             & -field_up(3,i-1,j,k+1,isub)+field_up(3,i-1,j,k-1,isub))

        dvyz = fac2*fac3*(field_up(2,i,j+1,k+1,isub)-field_up(2,i,j+1,k-1,isub)&
             & -field_up(2,i,j-1,k+1,isub)+field_up(2,i,j-1,k-1,isub))
        dwyz = fac2*fac3*(field_up(3,i,j+1,k+1,isub)-field_up(3,i,j+1,k-1,isub)&
             & -field_up(3,i,j-1,k+1,isub)+field_up(3,i,j-1,k-1,isub))

        duxx = fac4*(field_up(1,i+1,j,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i-1,j,k,isub))
        duyy = fac5*(field_up(1,i,j+1,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j-1,k,isub))
        duzz = fac6*(field_up(1,i,j,k+1,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j,k-1,isub))

        dvxx = fac4*(field_up(2,i+1,j,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i-1,j,k,isub))
        dvyy = fac5*(field_up(2,i,j+1,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j-1,k,isub))
        dvzz = fac6*(field_up(2,i,j,k+1,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j,k-1,isub))

        dwxx = fac4*(field_up(3,i+1,j,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i-1,j,k,isub))
        dwyy = fac5*(field_up(3,i,j+1,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j-1,k,isub))
        dwzz = fac6*(field_up(3,i,j,k+1,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j,k-1,isub))

        sforcenocaN(1+3*(ni-1)) = sforcenocaN(1+3*(ni-1)) + fac8 * (&
          & -u1 * u2 - 0.5D0 * u2 * (ty * o3 - tz * o2) + 0.5D0 * o2 * (ty * &
          & u3 - tz * u2) - 0.5D0 * ty * nu * (2.0_mk * duxx + dvxy + &
          & duyy + dwxz + duzz) + nu * (dvx + duy))
        sforcenocaN(2+3*(ni-1)) = sforcenocaN(2+3*(ni-1)) + fac8 * (&
          & 0.5D0 * u1 ** 2 - 0.5D0 * u2 ** 2 + 0.5D0 * u3 ** 2 - 0.5D0 * u2 &
          & * (tz * o1 - tx * o3) + 0.5D0 * o2 * (tz * u1 - tx * u3) + 0.5D0 &
          & * nu * (2.0_mk * duxx + dvxy + duyy + dwxz + duzz) * tx + &
          & 0.5D0 * nu * (dwxx + duxz + dwyy + dvyz + 2.0_mk * dwzz) * tz &
          & + 0.2D1 * nu * dvy)
        sforcenocaN(3+3*(ni-1)) = sforcenocaN(3+3*(ni-1)) + fac8 * (&
          & -u2 * u3 - 0.5D0 * u2 * (tx * o2 - ty * o1) + 0.5D0 * o2 * (tx * &
          & u2 - ty * u1) - 0.5D0 * ty * nu * (dwxx + duxz + dwyy + dvyz &
          & + 2.0_mk * dwzz) + nu * (dwy + dvz))

      END DO
    END DO
  END IF
  !---------------------------------------------------------------------------
  !  Noca gamma_imp surface integral: B
  !---------------------------------------------------------------------------
  IF ((ncimin(3) .GE. 1) .and. (ncimin(3) .LE. (ndata(3,isubl)-1))) THEN
    DO j=cncimin(2),cncimax(2)
      DO i=cncimin(1),cncimax(1)
        k=ncimin(3)
        lnocadimensions(3+3*(ni-1)) = lnocadimensions(3+3*(ni-1)) + dx*dy

        o1 = field_wp(1,i,j,k,isub)
        o2 = field_wp(2,i,j,k,isub)
        o3 = field_wp(3,i,j,k,isub)

        u1 = field_up(1,i,j,k,isub)
        u2 = field_up(2,i,j,k,isub)
        u3 = field_up(3,i,j,k,isub)

        tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
        ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
        tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

        dwx = fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
        
        dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))

        duz = fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))
        dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
        dwz = fac3*(field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub))

        duxy = fac1*fac2*(field_up(1,i+1,j+1,k,isub)-field_up(1,i+1,j-1,k,isub)&
             & -field_up(1,i-1,j+1,k,isub)+field_up(1,i-1,j-1,k,isub))
        dvxy = fac1*fac2*(field_up(2,i+1,j+1,k,isub)-field_up(2,i+1,j-1,k,isub)&
             & -field_up(2,i-1,j+1,k,isub)+field_up(2,i-1,j-1,k,isub))

        duxz = fac1*fac3*(field_up(1,i+1,j,k+1,isub)-field_up(1,i+1,j,k-1,isub)&
             & -field_up(1,i-1,j,k+1,isub)+field_up(1,i-1,j,k-1,isub))
        dwxz = fac1*fac3*(field_up(3,i+1,j,k+1,isub)-field_up(3,i+1,j,k-1,isub)&
             & -field_up(3,i-1,j,k+1,isub)+field_up(3,i-1,j,k-1,isub))

        dvyz = fac2*fac3*(field_up(2,i,j+1,k+1,isub)-field_up(2,i,j+1,k-1,isub)&
             & -field_up(2,i,j-1,k+1,isub)+field_up(2,i,j-1,k-1,isub))
        dwyz = fac2*fac3*(field_up(3,i,j+1,k+1,isub)-field_up(3,i,j+1,k-1,isub)&
             & -field_up(3,i,j-1,k+1,isub)+field_up(3,i,j-1,k-1,isub))

        duxx = fac4*(field_up(1,i+1,j,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i-1,j,k,isub))
        duyy = fac5*(field_up(1,i,j+1,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j-1,k,isub))
        duzz = fac6*(field_up(1,i,j,k+1,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j,k-1,isub))

        dvxx = fac4*(field_up(2,i+1,j,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i-1,j,k,isub))
        dvyy = fac5*(field_up(2,i,j+1,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j-1,k,isub))
        dvzz = fac6*(field_up(2,i,j,k+1,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j,k-1,isub))

        dwxx = fac4*(field_up(3,i+1,j,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i-1,j,k,isub))
        dwyy = fac5*(field_up(3,i,j+1,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j-1,k,isub))
        dwzz = fac6*(field_up(3,i,j,k+1,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j,k-1,isub))

        sforcenocaB(1+3*(ni-1)) = sforcenocaB(1+3*(ni-1)) - fac9 * (&
          & -u1 * u3 - 0.5D0 * u3 * (ty * o3 - tz * o2) + 0.5D0 * o3 * (ty * &
          & u3 - tz * u2) - 0.5D0 * tz * nu * (2.0_mk * duxx + dvxy + duyy + &
          & dwxz + duzz) + nu * (dwx + duz))
        sforcenocaB(2+3*(ni-1)) = sforcenocaB(2+3*(ni-1)) - fac9 * (&
          & -u2 * u3 - 0.5D0 * u3 * (tz * o1 - tx * o3) + 0.5D0 * o3 * (tz * &
          & u1 - tx * u3) - 0.5D0 * tz * nu * (dvxx + duxy + 2.0_mk * dvyy + &
          & dwyz + dvzz) + nu * (dwy + dvz))
        sforcenocaB(3+3*(ni-1)) = sforcenocaB(3+3*(ni-1)) - fac9 * (&
          & 0.5D0 * u1 ** 2 + 0.5D0 * u2 ** 2 - 0.5D0 * u3 ** 2 - 0.5D0 * u3 &
          & * (tx * o2 - ty * o1) + 0.5D0 * o3 * (tx * u2 - ty * u1) + 0.5D0 &
          & * nu * (2.0_mk * duxx + dvxy + duyy + dwxz + duzz) * tx + 0.5D0 * &
          & nu * (dvxx + duxy + 2.0_mk * dvyy + dwyz + dvzz) * ty + 0.2D1 * &
          & nu * dwz)

      END DO
    END DO
  END IF
  !---------------------------------------------------------------------------
  !  Noca gamma_imp surface integral: T
  !---------------------------------------------------------------------------
  IF ((ncimax(3) .GE. 1) .and. (ncimax(3) .LE. (ndata(3,isubl)-1))) THEN
    DO j=cncimin(2),cncimax(2)
      DO i=cncimin(1),cncimax(1)
        k=ncimax(3)
        lnocadimensions(3+3*(ni-1)) = lnocadimensions(3+3*(ni-1)) + dx*dy

        o1 = field_wp(1,i,j,k,isub)
        o2 = field_wp(2,i,j,k,isub)
        o3 = field_wp(3,i,j,k,isub)

        u1 = field_up(1,i,j,k,isub)
        u2 = field_up(2,i,j,k,isub)
        u3 = field_up(3,i,j,k,isub)

        tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
        ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
        tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

        dwx = fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
        
        dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))

        duz = fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))
        dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
        dwz = fac3*(field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub))

        duxy = fac1*fac2*(field_up(1,i+1,j+1,k,isub)-field_up(1,i+1,j-1,k,isub)&
             & -field_up(1,i-1,j+1,k,isub)+field_up(1,i-1,j-1,k,isub))
        dvxy = fac1*fac2*(field_up(2,i+1,j+1,k,isub)-field_up(2,i+1,j-1,k,isub)&
             & -field_up(2,i-1,j+1,k,isub)+field_up(2,i-1,j-1,k,isub))

        duxz = fac1*fac3*(field_up(1,i+1,j,k+1,isub)-field_up(1,i+1,j,k-1,isub)&
             & -field_up(1,i-1,j,k+1,isub)+field_up(1,i-1,j,k-1,isub))
        dwxz = fac1*fac3*(field_up(3,i+1,j,k+1,isub)-field_up(3,i+1,j,k-1,isub)&
             & -field_up(3,i-1,j,k+1,isub)+field_up(3,i-1,j,k-1,isub))

        dvyz = fac2*fac3*(field_up(2,i,j+1,k+1,isub)-field_up(2,i,j+1,k-1,isub)&
             & -field_up(2,i,j-1,k+1,isub)+field_up(2,i,j-1,k-1,isub))
        dwyz = fac2*fac3*(field_up(3,i,j+1,k+1,isub)-field_up(3,i,j+1,k-1,isub)&
             & -field_up(3,i,j-1,k+1,isub)+field_up(3,i,j-1,k-1,isub))

        duxx = fac4*(field_up(1,i+1,j,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i-1,j,k,isub))
        duyy = fac5*(field_up(1,i,j+1,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j-1,k,isub))
        duzz = fac6*(field_up(1,i,j,k+1,isub)-2.0_mk*field_up(1,i,j,k,isub) &
             & +field_up(1,i,j,k-1,isub))

        dvxx = fac4*(field_up(2,i+1,j,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i-1,j,k,isub))
        dvyy = fac5*(field_up(2,i,j+1,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j-1,k,isub))
        dvzz = fac6*(field_up(2,i,j,k+1,isub)-2.0_mk*field_up(2,i,j,k,isub) &
             & +field_up(2,i,j,k-1,isub))

        dwxx = fac4*(field_up(3,i+1,j,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i-1,j,k,isub))
        dwyy = fac5*(field_up(3,i,j+1,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j-1,k,isub))
        dwzz = fac6*(field_up(3,i,j,k+1,isub)-2.0_mk*field_up(3,i,j,k,isub) &
             & +field_up(3,i,j,k-1,isub))

        sforcenocaT(1+3*(ni-1)) = sforcenocaT(1+3*(ni-1)) + fac9 * (&
          & -u1 * u3 - 0.5D0 * u3 * (ty * o3 - tz * o2) + 0.5D0 * o3 * (ty * &
          & u3 - tz * u2) - 0.5D0 * tz * nu * (2.0_mk * duxx + dvxy + duyy + &
          & dwxz + duzz) + nu * (dwx + duz))
        sforcenocaT(2+3*(ni-1)) = sforcenocaT(2+3*(ni-1)) + fac9 * (&
          & -u2 * u3 - 0.5D0 * u3 * (tz * o1 - tx * o3) + 0.5D0 * o3 * (tz * &
          & u1 - tx * u3) - 0.5D0 * tz * nu * (dvxx + duxy + 2.0_mk * dvyy + &
          & dwyz + dvzz) + nu * (dwy + dvz))
        sforcenocaT(3+3*(ni-1)) = sforcenocaT(3+3*(ni-1)) + fac9 * (&
          & 0.5D0 * u1 ** 2 + 0.5D0 * u2 ** 2 - 0.5D0 * u3 ** 2 - 0.5D0 * u3 &
          & * (tx * o2 - ty * o1) + 0.5D0 * o3 * (tx * u2 - ty * u1) + 0.5D0 &
          & * nu * (2.0_mk * duxx + dvxy + duyy + dwxz + duzz) * tx + 0.5D0 * &
          & nu * (dvxx + duxy + 2.0_mk * dvyy + dwyz + dvzz) * ty + 0.2D1 * &
          & nu * dwz)

      END DO
    END DO
  END IF

  END DO
  END DO !ni

sforcenoca = REAL(sforcenocaW+sforcenocaE+sforcenocaS+sforcenocaN+ &
           & sforcenocaB+sforcenocaT,mk)
CALL MPI_Allreduce(vforcenoca,gvforcenoca,3*nforces, MPI_DOUBLE_PRECISION, &
     & MPI_SUM,comm,info)
CALL MPI_Allreduce(sforcenoca,gsforcenoca,3*nforces,mpi_prec,MPI_SUM,comm,info)
CALL MPI_Allreduce(lnocadimensions,gnocadimensions, 3*nforces, mpi_prec, &
     & MPI_SUM,comm,info)

  IF(rank.EQ.0) THEN
     tforcenoca = ((REAL(gvforcenoca,mk)-vforcenocaold)/dt + gsforcenoca)
     vforcenocaold=REAL(gvforcenoca,mk)
  END IF

END SUBROUTINE wvic_forcenoca

