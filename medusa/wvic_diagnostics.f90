!-------------------------------------------------------------------------------
!* filename: wvic_diagnostics                                                *!
!* project : ppm                                                              *!
!* purpose : compute some diagnostics                                         *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Aug 11 13:47:11 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_diagnostics.F,v $
!  Revision  2008/06/03 12:00:00  johannes
!  Removed double precision energy and enstrophy, extensive cleanup.
!  Structured output and added force calculations
!
!  Revision 1.5  2006/10/18 12:19:59  pchatela
!  Removed the free-stream from kinetic energy measurements!
!
!  Revision 1.4  2006/10/04 09:24:45  pchatela
!  Spatial evolution of the center of axial vorticity
!
!  Revision 1.3  2006/10/03 16:11:58  pchatela
!  Added spectra in the z direction (misnamed kx spectrum...)
!  Added spatial diagnostics, like kinetic energy, enstrophy, circulation
!  as functions of z, dumped at the frequency ndump
!
!  Revision 1.2  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.1  2005/09/28 11:40:31  michaebe
!  Fork from ppm_pvc
!
!-------------------------------------------------------------------------------

SUBROUTINE wvic_diagnostics (info)

  USE module_wvic
  USE ppm_module_write

  INTEGER , INTENT(inout) :: info
  INTEGER, PARAMETER :: md = kind(1.0d0)
  CHARACTER(len=512) :: filename,msg
  INTEGER  :: i,isub,isubl,j,k,ni
  REAL(mk) :: lens ! local enstrophy
  REAL(mk) :: gens ! global enstrophy
!  REAL(mk), DIMENSION(3) :: force ! global moment of viscosity
  REAL(mk) :: tx, ty, tz
  REAL(mk) :: lensobj ! local object enstrophy
  REAL(mk) :: gensobj ! global object enstrophy
  REAL(mk) :: leng ! local kinetic energy      
  REAL(mk) :: geng ! global kinetic energy      
  REAL(mk) :: lfou ! local fourier number
  REAL(mk) :: gfou ! global fourier number
  REAL(mk) :: ldiv,ldiv1 ! local divergence of vorticity
  REAL(mk) :: gdiv ! global divergence of vorticity
  REAL(mk) :: ludiv,ludiv1 ! local divergence of velocity
  REAL(mk) :: gudiv ! global divergence of velocity
  REAL(mk) :: gshr, lshr
  REAL(mk) :: gvor, lvor 
  REAL(mk) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, dsu, dsv, dsw
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9
  REAL(mk) :: dv, umax, midplane
  REAL(mk) :: lcfl, gcfl
  REAL(mk) :: newdt,newlambda
  REAL(mk), DIMENSION(3) :: u_wcmass, rot, erot
  REAL(mk) :: lrot, grot  ! local, global error
  INTEGER  :: ios

  REAL(mk) :: duxy, dvxy, duxz, dwxz, dvyz, dwyz 
  REAL(mk) :: duxx, duyy, duzz, dvxx, dvyy, dvzz, dwxx, dwyy, dwzz 
  REAL(mk) :: u1, u2, u3
  REAL(md), DIMENSION(3) :: ldp,gdp
  REAL(md), DIMENSION(3) :: ldpt,gdpt
  REAL(mk), DIMENSION(3) :: forcepenal, gforcepenal 
  REAL(mk), DIMENSION(3) :: forcewmoment, lforcewmoment, gforcewmoment
  REAL(mk), DIMENSION(3*nforces) :: forcenoca
  REAL(mk), DIMENSION(3*nforces) :: nocadimensions

  INCLUDE 'mpif.h'

  dv = dx*dy*dz

  fac1 = 0.5_mk / dx
  fac2 = 0.5_mk / dy
  fac3 = 0.5_mk / dz
  IF(wvic_compvel_scheme.GT.0) THEN
     fac1 = 8.0_mk/dx/12.0_mk
     fac2 = 8.0_mk/dy/12.0_mk
     fac3 = 8.0_mk/dz/12.0_mk
     fac4 = -1.0_mk/dx/12.0_mk
     fac5 = -1.0_mk/dy/12.0_mk
     fac6 = -1.0_mk/dz/12.0_mk
  END IF
  
  !-----------------------------------------------------------------------------
  !  Set initial values. Ghost
  !-----------------------------------------------------------------------------
  lens = 0.0_mk
  lensobj = 0.0_mk
  leng = 0.0_mk
  ldiv = 0.0_mk
  ludiv = 0.0_mk
  umax = -HUGE(umax)
  lshr = -HUGE(lshr)
  CALL wvic_ghost(wvic_prm_vorticity,info)
  CALL wvic_ghost(wvic_prm_velocity,info)
  lrot = 0.0_mk
  lvor = -HUGE(lvor)
  forcepenal = 0.0_mk
  lforcewmoment = 0.0_mk
!  leng_d = 0.0_md
!  lens_d = 0.0_mk
  midplane = (min_physg(1)+max_physg(1))*0.5_mk
  DO isub=1,nsublist
    isubl = isublist(isub)
    
    DO k=1,ndata(3,isubl)-1
        
      DO j=1,ndata(2,isubl)-1
           
        DO i=1,ndata(1,isubl)-1

          tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
          ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
          tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

          !-------------------------------------------------------------------
          !  Compute  "omega - curl(u)"
          !-------------------------------------------------------------------
          IF(SUM(ghostsize).EQ.3) THEN
            rot(1)=fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))-&
                 & fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
            rot(2)=fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))-&
                 & fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
            rot(3)=fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))-&
                 & fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
          ELSE
            rot(1)=fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))-&
                 & fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))+&
                 & fac5*(field_up(3,i,j+2,k,isub)-field_up(3,i,j-2,k,isub))-&
                 & fac6*(field_up(2,i,j,k+2,isub)-field_up(2,i,j,k-2,isub))

            rot(2)=fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))-&
                 & fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))+&
                 & fac6*(field_up(1,i,j,k+2,isub)-field_up(1,i,j,k-2,isub))-&
                 & fac4*(field_up(3,i+2,j,k,isub)-field_up(3,i-2,j,k,isub))
            
            rot(3)=fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))-&
                 & fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))+&
                 & fac4*(field_up(2,i+2,j,k,isub)-field_up(2,i-2,j,k,isub))-&
                 & fac5*(field_up(1,i,j+2,k,isub)-field_up(1,i,j-2,k,isub))
          END IF
          erot = rot - field_wp(1:3,i,j,k,isub)
          lrot = lrot + SUM(erot**2)*dv

          !-----------------------------------------------------
          ! Compute divergence of vorticity and velocity
          !-----------------------------------------------------
          ldiv1 =fac1*(field_wp(1,i+1,j,k,isub)-field_wp(1,i-1,j,k,isub))+&
               & fac2*(field_wp(2,i,j+1,k,isub)-field_wp(2,i,j-1,k,isub))+&
               & fac3*(field_wp(3,i,j,k+1,isub)-field_wp(3,i,j,k-1,isub))
          ldiv = ldiv + field_H(i,j,k,isub)*ldiv1**2 * dv !JTR OBS JTMP remove field_H
          ludiv1 =fac1*(field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub))+&
               & fac2*(field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub))+&
               & fac3*(field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub))
          ludiv = ludiv + field_H(i,j,k,isub)*ludiv1**2 * dv !JTR OBS JTMP remove field_H

          !-----------------------------------------------------
          ! do it on the particles !JTR?
          !-----------------------------------------------------
          u_wcmass(1) = field_up(1,i,j,k,isub)
          u_wcmass(2) = field_up(2,i,j,k,isub)
          u_wcmass(3) = field_up(3,i,j,k,isub)
          !-----------------------------------------------------
          ! Compute kinetic energy wo/cmass
          !-----------------------------------------------------
          leng = leng + dv*SUM(u_wcmass**2)
          !-----------------------------------------------------
          ! Compute domain enstrophy and enstrophy inside object
          !-----------------------------------------------------
          lensobj = lensobj + dv*SUM(field_wp(1:3,i,j,k,isub)**2) &
                    & *field_H(i,j,k,isub)
          lens = lens + dv*SUM(field_wp(1:3,i,j,k,isub)**2)

!          IF(mk.NE.md) THEN
!             leng_d = leng_d + DBLE(SUM(u_wcmass**2))
!             lens_d = lens_d + DBLE(SUM(field_wp(1:3,i,j,k,isub)**2))
!          END IF

          !-----------------------------------------------------
          ! Maxmimum velcity       
          !-----------------------------------------------------
          umax = MAX(umax,MAXVAL(ABS(field_up(:,i,j,k,isub))))

          !-----------------------------------------------------
          ! Compute shear of flow  
          !-----------------------------------------------------
          dux = fac1*(field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub))
          dvx = fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))
          dwx = fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
              
          duy = fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
          dvy = fac2*(field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub))
          dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))

          duz = fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))
          dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
          dwz = fac3*(field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub))

          dsu = ABS(dux)+ABS(duy)+ABS(duz)
          dsv = ABS(dvx)+ABS(dvy)+ABS(dvz)
          dsw = ABS(dwx)+ABS(dwy)+ABS(dwz)

          lshr = MAX(lshr,dsu**2 + dsv**2 + dsw**2)

          !-----------------------------------------------------
          ! Compute maximum vorticity (magnitude)
          !-----------------------------------------------------
          lvor = MAX(lvor, SUM(field_wp(1:3,i,j,k,isub)**2))

          !-----------------------------------------------------
          ! Compute force a la Vasilyev / Angot
          ! (this is actually the force on the fluid - heed signs)
          ! At this point the new velocity field has already been calculated
          ! from the vorticity field. It is therefore not the correct 
          ! velocity field that is being used in this expression to calculate
          ! the penalization forces. ubar is ok.
          ! HENCE move this to a previous point!!!JTR
          !-----------------------------------------------------
          forcepenal(1) = forcepenal(1) + lambda*field_H(i,j,k,isub)* &
                     & (field_ubar(1,i,j,k,isub)-field_up(1,i,j,k,isub))
          forcepenal(2) = forcepenal(2) + lambda*field_H(i,j,k,isub)* &
                     & (field_ubar(2,i,j,k,isub)-field_up(2,i,j,k,isub))
          forcepenal(3) = forcepenal(3) + lambda*field_H(i,j,k,isub)* &
                     & (field_ubar(3,i,j,k,isub)-field_up(3,i,j,k,isub))

          !-----------------------------------------------------
          ! Compute vorticity moments
          !-----------------------------------------------------
          lforcewmoment(1) = lforcewmoment(1) - &
               & 0.5_mk*(ty*field_wp(3,i,j,k,isub) - tz*field_wp(2,i,j,k,isub))
          lforcewmoment(2) = lforcewmoment(2) - & 
               & 0.5_mk*(tz*field_wp(1,i,j,k,isub) - tx*field_wp(3,i,j,k,isub))
          lforcewmoment(3) = lforcewmoment(3) - &
               & 0.5_mk*(tx*field_wp(2,i,j,k,isub) - ty*field_wp(1,i,j,k,isub))

        END DO

      END DO

    END DO

  END DO

  !-----------------------------------------------------
  !  Post integration calculations
  !-----------------------------------------------------
  lforcewmoment = dv*lforcewmoment
  forcepenal = dv*forcepenal



  !-----------------------------------------------------------------------------
  !  cfl number
  !-----------------------------------------------------------------------------
  lcfl = umax / (dx/dt)
  !-----------------------------------------------------
  !  fourier number
  !-----------------------------------------------------
  lfou = nu * dt / (dx**2) 
  
  !-----------------------------------------------------------------------
  !  Integrate pressure along xyz-axis   
  !-----------------------------------------------------------------------
  ldp = 0.0_md
  ldpt = 0.0_md
  fac1 = 0.5_mk / dx
  fac2 = 0.5_mk / dy
  fac3 = 0.5_mk / dz
  fac4 = 1.0_mk / dx**2
  fac5 = 1.0_mk / dy**2
  fac6 = 1.0_mk / dz**2
  fac7 = dy*dz
  fac8 = dx*dz
  fac9 = dx*dy
  !-----------------------------------------------------------------------
  ! x-axis
  !-----------------------------------------------------------------------
  IF ((min_physg(2) .eq. min_sub(2,isubl)) .AND. &
     & (min_physg(3) .eq. min_sub(3,isubl))) THEN
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO i=1,ndata(1,isubl)-1
      j=1
      k=1
      u1 = field_up(1,i,j,k,isub)
      u2 = field_up(2,i,j,k,isub)
      u3 = field_up(3,i,j,k,isub)
      duxx = fac4*(field_up(1,i+1,j,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
           & +field_up(1,i-1,j,k,isub))
      duyy = fac5*(field_up(1,i,j+1,k,isub)-2.0_mk*field_up(1,i,j,k,isub) &
           & +field_up(1,i,j-1,k,isub))
      duzz = fac6*(field_up(1,i,j,k+1,isub)-2.0_mk*field_up(1,i,j,k,isub) &
           & +field_up(1,i,j,k-1,isub))
      dux = fac1*(field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub))
      duy = fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
      duz = fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))

      ldp(1) = ldp(1) + dx*(nu * (duxx+duyy+duzz) - (u1*dux + u2*duy + u3*duz))
      ldpt(1) = ldpt(1) - dx*(u1-u1old(i,isub))/dt 
      u1old(i,isub)=u1
    END DO
  END DO
  END IF
  !-----------------------------------------------------------------------
  ! y-axis
  !-----------------------------------------------------------------------
  IF ((min_physg(1) .eq. min_sub(1,isubl)) .AND. &
     & (min_physg(3) .eq. min_sub(3,isubl))) THEN
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO j=1,ndata(2,isubl)-1
      k=1
      i=1
      u1 = field_up(1,i,j,k,isub)
      u2 = field_up(2,i,j,k,isub)
      u3 = field_up(3,i,j,k,isub)
      dvxx = fac4*(field_up(2,i+1,j,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
           & +field_up(2,i-1,j,k,isub))
      dvyy = fac5*(field_up(2,i,j+1,k,isub)-2.0_mk*field_up(2,i,j,k,isub) &
           & +field_up(2,i,j-1,k,isub))
      dvzz = fac6*(field_up(2,i,j,k+1,isub)-2.0_mk*field_up(2,i,j,k,isub) &
           & +field_up(2,i,j,k-1,isub))
      dvx = fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))
      dvy = fac2*(field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub))
      dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
      ldp(2) = ldp(2) + dy*(nu * (dvxx+dvyy+dvzz) - (u1*dvx + u2*dvy + u3*dvz))
      ldpt(2) = ldpt(2) - dy*(u2-u2old(j,isub))/dt 
      u2old(j,isub)=u2
    END DO
  END DO
  END IF
  !-----------------------------------------------------------------------------
  ! z-axis
  !-----------------------------------------------------------------------------
  IF ((min_physg(1) .eq. min_sub(1,isubl)) .AND. &
     & (min_physg(2) .eq. min_sub(2,isubl))) THEN
  DO isub=1,nsublist
    isubl = isublist(isub)
    DO k=1,ndata(3,isubl)-1
      j=1 
      i=1 
      u1 = field_up(1,i,j,k,isub)
      u2 = field_up(2,i,j,k,isub)
      u3 = field_up(3,i,j,k,isub)
      dwxx = fac4*(field_up(3,i+1,j,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
           & +field_up(3,i-1,j,k,isub))
      dwyy = fac5*(field_up(3,i,j+1,k,isub)-2.0_mk*field_up(3,i,j,k,isub) &
           & +field_up(3,i,j-1,k,isub))
      dwzz = fac6*(field_up(3,i,j,k+1,isub)-2.0_mk*field_up(3,i,j,k,isub) &
           & +field_up(3,i,j,k-1,isub))
      dwx = fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
      dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))
      dwz = fac3*(field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub))
      ldp(3) = ldp(3) + dz*(nu * (dwxx+dwyy+dwzz) - (u1*dwx + u2*dwy + u3*dwz))
      ldpt(3) = ldpt(3) - dz*(u3-u3old(k,isub))/dt 
      u3old(k,isub)=u3
    END DO
  END DO
  END IF

  !-----------------------------------------------------------------------------
  !  Output velocity profile   !JTR obsolete
  !-----------------------------------------------------------------------------
  IF ((flow_case .EQ. 9) .AND. (0 .EQ. 1)) THEN
    IF ((min_physg(1) .eq. min_sub(1,isubl)) .AND. &
       & (min_physg(3) .eq. min_sub(3,isubl))) THEN
      DO isub=1,nsublist
        isubl = isublist(isub)
        i=1
        k=1
        WRITE(filename,'(A,I3.3,A)') 'poiseuille',rank,'.dat'
        OPEN(14,file=filename,iostat=ios,position='append',status='unknown')
        ! write coordinates on first line of file
        IF (itime .EQ. 0) THEN
          WRITE(14,'((E12.5,A),$)') time , ' '
          DO j=1,ndata(2,isubl)-1
            ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
            WRITE(14,'((E12.5,A),$)') ty, ' '
          END DO
          WRITE(14,*)
        END IF
        ! write time and velocities - time is first figure on each line
        WRITE(14,'((E12.5,A),$)') time , ' '
        DO j=1,ndata(2,isubl)-1
          WRITE(14,'((E12.5,A),$)') field_up(3,i,j,k,isub) , ' '
        END DO
        WRITE(14,*)
        CLOSE(14)
      END DO
    END IF
  END IF
 
  !-----------------------------------------------------------------------------
  ! Calculate forces via Noca
  !-----------------------------------------------------------------------------
  CALL wvic_forcenoca (forcenoca,nocadimensions,info)

!  leng_d = leng_d * DBLE(dv)
!  lens_d = lens_d * DBLE(dv)
  !-----------------------------------------------------------------------------
  !  send it to the root
  !-----------------------------------------------------------------------------
  CALL MPI_Reduce(lens,gens,1, mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(ldp,gdp,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(ldpt,gdpt,3,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(forcepenal,gforcepenal,3,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lforcewmoment,gforcewmoment,3,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lensobj,gensobj,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(leng,geng,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(ldiv,gdiv,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(ludiv,gudiv,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lcfl,gcfl,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lshr,gshr,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lvor,gvor,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lfou,gfou,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lrot,grot,1,mpi_prec,MPI_SUM,0,comm,info)
!  IF(md.NE.mk) THEN
!     CALL MPI_Reduce(lens_d,gens_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
!     CALL MPI_Reduce(leng_d,geng_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
!  ELSE
!     gens_d = gens
!     geng_d = geng
!  END IF

  !-------------------------------------------------------------------------
  ! Get solid mass/volume
  !-------------------------------------------------------------------------
  CALL wvic_calculate_mass(info)

  !-------------------------------------------------------------------------
  ! Force from medusa
  !-------------------------------------------------------------------------
  IF (flow_case .EQ. 11) THEN
    medu_thrusty = 0.0_mk
    IF ((rank .EQ. 0) .AND. (itime .NE. 0)) THEN
      medu_thrusty = -(max_physg(1)-min_physg(1)) * &
                   & (max_physg(3)-min_physg(3))*REAL(gdp(2)+gdpt(2))
    END IF
  CALL MPI_BCast(medu_thrusty,1,mpi_prec,0,comm,info)
  END IF
     
i = LEN_TRIM(runtag)
  IF(rank.EQ.0) THEN
     forcewmoment = (gforcewmoment - gforcewmomentold)/dt
     gforcewmomentold = gforcewmoment
!     !-------------------------------------------------------------------------
!     ! convert pressure drop to average pressure gradient
!     !-------------------------------------------------------------------------
!     gdp(1)  = DBLE(gdp(1)/(max_physg(1)-min_physg(1)))
!     gdpt(1) = DBLE(gdpt(1)/(max_physg(1)-min_physg(1)))
!     gdp(2)  = DBLE(gdp(2)/(max_physg(2)-min_physg(2)))
!     gdpt(2) = DBLE(gdpt(2)/(max_physg(2)-min_physg(2)))
!     gdp(3)  = DBLE(gdp(3)/(max_physg(3)-min_physg(3)))
!     gdpt(3) = DBLE(gdpt(3)/(max_physg(3)-min_physg(3)))
!
     !--------------------------------------------------------------------------
     ! Determine adaptive time step - brodcast it after appending output .dat
     !--------------------------------------------------------------------------
!JTR!!
!     IF (dt_adapt) THEN
       newdt = dt_max !MIN(dt_max,dt)
       !-----------------------------------------------------------------------
       ! Time constraint from 
       !-----------------------------------------------------------------------
       IF(gvor .NE. 0.0_mk) THEN
          newdt = MIN(newdt,0.25_mk/SQRT(gvor)*0.5_mk)
       END IF
       ! As long as this time step is an ok initial timestep leave this be
       !-----------------------------------------------------------------------
       ! Time constraint from 
       !-----------------------------------------------------------------------
       newdt = MIN(newdt,0.48_mk*dx**2/(nu*3.0_mk))
       !-----------------------------------------------------------------------
       ! Time constraint from flow shear
       !-----------------------------------------------------------------------
       IF(gshr .NE. 0.0_mk) THEN
         newdt = MIN(newdt,0.5_mk/SQRT(gshr))
       END IF

       newlambda = penalization_lambda/newdt
       newdt = newdt*dt_adapt_fraction

       IF ((dt_rampsteps .NE. 0) .and. (dt_rampsteps .GT. itime)) THEN
         newdt = newdt*(REAL(itime,mk)/REAL(dt_rampsteps,mk)*(1-dt_rampstart) &
               & +dt_rampstart)
       END IF
!     END IF
 
     WRITE(filename,'(A,A)') runtag(1:iruntag),'-diag.dat'
     OPEN(14,file=filename,iostat=ios,position='append',status='unknown')
     IF (itime .EQ. 0) THEN
        WRITE(14, '(A)') '# itime, time, kinetic energy w/o mass' // &
          'global enstrophy, enstrophy in object, dt,' // &
          'max(omega), omega-rot(u), div(omega),' // &
          'div(u), norm(u_solid)**2,' // &
          'force_penal_x, force_penal_y, force_penal_z,' // &
          'force_wmoment_x, force_wmoment_y, force_wmoment_z,' // &
          'Dp_x, Dp_y, Dp_z,' // &
          '(forcenoca&nocadimensions) x nforces'
     END IF
     WRITE(14,'(I16,A, 19(E12.5,A),$)') &
          & itime,' ',time,' ',geng,' ', & !3
          & gens,' ',gensobj,' ',dt,' ', & !6
          & SQRT(gvor),' ',grot,' ',gdiv,' ', & !9
          & gudiv,' ',medu_thrusty,' ',& !11
          & gforcepenal(1),' ',gforcepenal(2),' ',gforcepenal(3),' ', & !14
          & medu_move_pos,' ',medu_move_vel,' ',0.0_mk,' ',&!17
          & REAL(gdp(1)+gdpt(1)),' ',REAL(gdp(2)+gdpt(2)),' ', & !19
          & REAL(gdp(3)+gdpt(3)),' ' !20
     DO ni=1,nforces 
       WRITE(14,'(6(E12.5,A),$)') &
          & ,forcenoca(1+3*(ni-1)),' ' ,forcenoca(2+3*(ni-1)),' ' &
          & ,forcenoca(3+3*(ni-1)),' ' ,nocadimensions(1+3*(ni-1)),' ' &
          & ,nocadimensions(2+3*(ni-1)),' ' ,nocadimensions(3+3*(ni-1)),' '
     END DO
     WRITE(14,*)
     CLOSE(14)

  END IF
  !-----------------------------------------------------------------------------
  ! Broadcast new dt and lambda
  !-----------------------------------------------------------------------------
  IF (dt_adapt) THEN
    CALL MPI_BCast(newdt,1,mpi_prec,0,comm,info)
    CALL MPI_BCast(newlambda,1,mpi_prec,0,comm,info)
    dt=newdt
    lambda = newlambda
!    CALL MPI_BCast(max_vorticity,1,mpi_prec,0,comm,info)
  ELSE
    CALL MPI_BCast(newlambda,1,mpi_prec,0,comm,info)
    lambda = newlambda
  END IF

2002 CONTINUE
!  CALL MPI_BCast(max_vorticity,1,mpi_prec,0,comm,info)
END SUBROUTINE wvic_diagnostics


SUBROUTINE wvic_diagnostics_trail(info)

  USE module_wvic
  USE ppm_module_write

  INTEGER , INTENT(inout) :: info
  INTEGER, PARAMETER :: md = kind(1.0d0)
  CHARACTER(len=256) :: filename,msg
  INTEGER  :: i,isub,isubl,j,k,ldu
  REAL(mk) :: lcir ! local circlulation
  REAL(mk) :: gcir ! global circulateion
  REAL(mk) :: lens ! local circlulation
  REAL(mk) :: gens ! global circulateion
  REAL(mk) :: leng ! local circlulation
  REAL(mk) :: geng ! global circulateion
  REAL(mk) :: lfou ! local fourier number
  REAL(mk) :: gfou ! global fourier number
  REAL(mk) :: ldiv,ldiv1 ! local circlulation
  REAL(mk) :: gdiv ! global circulateion
  
  ! Arrays for the spatial evolutions
  REAL(mk), DIMENSION(:), POINTER :: lcirofx, gcirofx
  REAL(mk), DIMENSION(:), POINTER :: lensofx, gensofx
  REAL(mk), DIMENSION(:), POINTER :: lengofx, gengofx
  
  REAL(mk) :: gshr, lshr
  REAL(mk) :: gvor, lvor, grcir,rcir
  REAL(mk) :: dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, dsu, dsv, dsw
  REAL(mk) :: fac1, fac2, fac3, fac4, fac5, fac6
  REAL(mk) :: tx, ty, tz
  REAL(mk) :: dv, umax, xplane(3), midplane
  REAL(mk) :: lcfl, gcfl
  REAL(mk) :: newdt
  REAL(mk), DIMENSION(3) :: u_wcmass, rot, erot
  REAL(mk) :: lrot, grot
  INTEGER  :: ios

  !-----------------------------------------------------
  ! some double precision stuff
  !-----------------------------------------------------
  REAL(md) :: lens_d, gens_d, leng_d, geng_d, tmp1_d, tmp2_d
  REAL(md) :: lrcirr_d, grcirr_d, lrcirl_d, grcirl_d
  REAL(md) :: lrcir25r_d, lrcir25l_d, grcir25r_d, grcir25l_d
  REAL(md) :: lrcir26r_d, lrcir26l_d, grcir26r_d, grcir26l_d
  REAL(md), DIMENSION(:), POINTER :: lcirofx_d, gcirofx_d
  REAL(md), DIMENSION(:), POINTER :: lensofx_d, gensofx_d
  REAL(md), DIMENSION(:), POINTER :: lengofx_d, gengofx_d
  REAL(md), DIMENSION(:), POINTER :: lxomegaofx_d, gxomegaofx_d
  REAL(md), DIMENSION(:), POINTER :: lyomegaofx_d, gyomegaofx_d
  
  INCLUDE 'mpif.h'

  dv = dx*dy*dz

  fac1 = 0.5_mk / dx
  fac2 = 0.5_mk / dy
  fac3 = 0.5_mk / dz
  IF(wvic_compvel_scheme.GT.0) THEN
     fac1 = 8.0_mk/dx/12.0_mk
     fac2 = 8.0_mk/dy/12.0_mk
     fac3 = 8.0_mk/dz/12.0_mk
     fac4 = -1.0_mk/dx/12.0_mk
     fac5 = -1.0_mk/dy/12.0_mk
     fac6 = -1.0_mk/dz/12.0_mk
  END IF
  
  !-----------------------------------------------------------------------------
  !  compute the local circulation
  !-----------------------------------------------------------------------------
  lcir = SUM(wp(1:3,1:np))*dv
  !-----------------------------------------------------------------------------
  !  compute the enstropy <-- this is BAD, dont compare bananas and kiwis
  !-----------------------------------------------------------------------------
  !  lens = SUM(wp(1:3,1:np)**2)*dv
  lens = 0.0_mk
  lrcir25r_d = 0.0_md
  lrcir26r_d = 0.0_md
  lrcir25l_d = 0.0_md
  lrcir26l_d = 0.0_md
  lrcirl_d = 0.0_md
  lrcirr_d = 0.0_md
  
  !-----------------------------------------------------------------------------
  !  compute the energy, and the divergence of the vorticity
  !  Wed Dec 29 23:30:01 CET 2004: compute energy and enstrophy in double
  !  precision
  !-----------------------------------------------------------------------------
  leng = 0.0_mk
  ldiv = 0.0_mk
  umax = -HUGE(umax)
  lshr = -HUGE(lshr)
  CALL wvic_ghost(wvic_prm_vorticity,info)
  CALL wvic_ghost(wvic_prm_velocity,info)
  lrot = 0.0_mk
  lvor = -HUGE(lvor)
  leng_d = 0.0_md
  lens_d = 0.0_mk
  rcir   = 0.0_mk
  !-----------------------------------------------------------------------------
  !  Compute spatial evolutions of energy, enstrophy and circulation
  !-----------------------------------------------------------------------------
  ldu = nx(3)
  ALLOCATE(lcirofx_d(ldu),gcirofx_d(ldu))
  ALLOCATE(lensofx_d(ldu),gensofx_d(ldu))
  ALLOCATE(lengofx_d(ldu),gengofx_d(ldu))
  ALLOCATE(lxomegaofx_d(ldu),gxomegaofx_d(ldu))
  ALLOCATE(lyomegaofx_d(ldu),gyomegaofx_d(ldu))
  lcirofx_d = 0.0_md
  gcirofx_d = 0.0_md
  lensofx_d = 0.0_md
  gensofx_d = 0.0_md
  lengofx_d = 0.0_md
  gengofx_d = 0.0_md
  lxomegaofx_d = 0.0_md
  gxomegaofx_d = 0.0_md
  lyomegaofx_d = 0.0_md
  gyomegaofx_d = 0.0_md
  
  midplane = (min_physg(1)+max_physg(1))*0.5_mk
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)-1
        
        DO j=1,ndata(2,isubl)-1
           
           DO i=1,ndata(1,isubl)-1
           
              tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
              ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
              tz = min_sub(3,isubl) + REAL(k-1,mk)*dz

              IF((i+istart(1,isubl)-1).GT.(nx(1)-1)/2 &
                   & .AND. (j+istart(2,isubl)-1).EQ.(nx(2)-1)/2) THEN
                 lrcirl_d = lrcirl_d + DBLE(dx*dy*field_wp(2,i,j,k,isub))
                 IF(field_wp(2,i,j,k,isub).GT.2.5_mk) THEN
                    lrcir25l_d = lrcir25l_d + &
                         & DBLE(dx*dy*field_wp(2,i,j,k,isub))
                 ELSEIF(field_wp(2,i,j,k,isub).LT.-2.5_mk) THEN
                    lrcir26l_d = lrcir26l_d + &
                         & DBLE(dx*dy*field_wp(2,i,j,k,isub))
                 END IF
                    
              ELSEIF( (j+istart(2,isubl)-1).EQ.(nx(2)-1)/2) THEN
                 lrcirr_d = lrcirr_d + DBLE(dx*dy*field_wp(2,i,j,k,isub))
                 IF(field_wp(2,i,j,k,isub).LT.-2.5_mk) THEN
                    lrcir25r_d = lrcir25r_d + &
                         & DBLE(dx*dy*field_wp(2,i,j,k,isub))
                 ELSEIF(field_wp(2,i,j,k,isub).GT.2.5_mk) THEN
                    lrcir26r_d = lrcir26r_d + &
                         & DBLE(dx*dy*field_wp(2,i,j,k,isub))
                 END IF
              END IF
              
              IF(SUM(ghostsize).EQ.3) THEN
                 rot(1)=fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))-&
                      & fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
                 rot(2)=fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))-&
                      & fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
                 rot(3)=fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))-&
                      & fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
              ELSE
                 rot(1)=fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))-&
                      & fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))+&
                      & fac5*(field_up(3,i,j+2,k,isub)-field_up(3,i,j-2,k,isub))-&
                      & fac6*(field_up(2,i,j,k+2,isub)-field_up(2,i,j,k-2,isub))

                 rot(2)=fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))-&
                      & fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))+&
                      & fac6*(field_up(1,i,j,k+2,isub)-field_up(1,i,j,k-2,isub))-&
                      & fac4*(field_up(3,i+2,j,k,isub)-field_up(3,i-2,j,k,isub))
                 
                 rot(3)=fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))-&
                      & fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))+&
                      & fac4*(field_up(2,i+2,j,k,isub)-field_up(2,i-2,j,k,isub))-&
                      & fac5*(field_up(1,i,j+2,k,isub)-field_up(1,i,j-2,k,isub))
              END IF
              erot = rot - field_wp(1:3,i,j,k,isub)
              lrot = lrot + SUM(erot**2)*dv
              ldiv1 =fac1*(field_wp(1,i+1,j,k,isub)-field_wp(1,i-1,j,k,isub))+&
                   & fac2*(field_wp(2,i,j+1,k,isub)-field_wp(2,i,j-1,k,isub))+&
                   & fac3*(field_wp(3,i,j,k+1,isub)-field_wp(3,i,j,k-1,isub))
              ldiv = ldiv + ldiv1**2 * dv
              !-----------------------------------------------------
              ! do it on the particles
              !-----------------------------------------------------
              u_wcmass(1) = field_up(1,i,j,k,isub)-u_infty(1)
              u_wcmass(2) = field_up(2,i,j,k,isub)-u_infty(2)
              u_wcmass(3) = field_up(3,i,j,k,isub)-u_infty(3)
              !-----------------------------------------------------
              ! kinetic energy wo/cmass
              !-----------------------------------------------------
              leng = leng + dv*SUM(u_wcmass**2)
              lens = lens + dv*SUM(field_wp(1:3,i,j,k,isub)**2)
              IF(mk.NE.md) THEN
                 leng_d = leng_d + DBLE(SUM(u_wcmass**2))
                 lens_d = lens_d + DBLE(SUM(field_wp(1:3,i,j,k,isub)**2))
              END IF
              
              !-----------------------------------------------------
              ! Spatial evolutions of kinetic energy, 
              ! circulation through half plane,
              ! position of center of vorticity...
              lensofx_d(k+istart(3,isubl)-1) = lensofx_d(k+istart(3,isubl)-1) + &
                                            & DBLE(SUM(field_wp(1:3,i,j,k,isub)**2))
              lengofx_d(k+istart(3,isubl)-1) = lengofx_d(k+istart(3,isubl)-1) + &
                                            & DBLE(SUM(u_wcmass**2))
              IF ( (j+istart(2,isubl)-1).LT.(nx(2)-1)/2) THEN
                  lcirofx_d(k+istart(3,isubl)-1) = lcirofx_d(k+istart(3,isubl)-1) + &
                                            & DBLE(dx*dy*field_wp(3,i,j,k,isub))
                  lxomegaofx_d(k+istart(3,isubl)-1) = lxomegaofx_d(k+istart(3,isubl)-1) + &
                                            & DBLE(dx*dy*tx*field_wp(3,i,j,k,isub))
                  lyomegaofx_d(k+istart(3,isubl)-1) = lyomegaofx_d(k+istart(3,isubl)-1) + &
                                            & DBLE(dx*dy*ty*field_wp(3,i,j,k,isub))
              END IF
              umax = MAX(umax,MAXVAL(ABS(field_up(:,i,j,k,isub))))

              dux = fac1*(field_up(1,i+1,j,k,isub)-field_up(1,i-1,j,k,isub))
              dvx = fac1*(field_up(2,i+1,j,k,isub)-field_up(2,i-1,j,k,isub))
              dwx = fac1*(field_up(3,i+1,j,k,isub)-field_up(3,i-1,j,k,isub))
              
              duy = fac2*(field_up(1,i,j+1,k,isub)-field_up(1,i,j-1,k,isub))
              dvy = fac2*(field_up(2,i,j+1,k,isub)-field_up(2,i,j-1,k,isub))
              dwy = fac2*(field_up(3,i,j+1,k,isub)-field_up(3,i,j-1,k,isub))

              duz = fac3*(field_up(1,i,j,k+1,isub)-field_up(1,i,j,k-1,isub))
              dvz = fac3*(field_up(2,i,j,k+1,isub)-field_up(2,i,j,k-1,isub))
              dwz = fac3*(field_up(3,i,j,k+1,isub)-field_up(3,i,j,k-1,isub))

              dsu = ABS(dux)+ABS(duy)+ABS(duz)
              dsv = ABS(dvx)+ABS(dvy)+ABS(dvz)
              dsw = ABS(dwx)+ABS(dwy)+ABS(dwz)

              lshr = MAX(lshr,dsu**2 + dsv**2 + dsw**2)

              lvor = MAX(lvor, SUM(field_wp(1:3,i,j,k,isub)**2))
           END DO


        END DO

     END DO

  END DO

  !-----------------------------------------------------------------------------
  !  cfl number
  !-----------------------------------------------------------------------------
  lcfl = umax / (dx/dt)
  !-----------------------------------------------------
  !  fourier number
  !-----------------------------------------------------
  lfou = nu * dt / (dx**2) 
  
  leng_d = leng_d * DBLE(dv)
  lens_d = lens_d * DBLE(dv)
  lengofx_d = lengofx_d * DBLE(dx*dy)
  lensofx_d = lensofx_d * DBLE(dx*dy)
  !-----------------------------------------------------------------------------
  !  send it to the root
  !-----------------------------------------------------------------------------
  CALL MPI_Reduce(lcir,gcir,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lens,gens,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(leng,geng,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(ldiv,gdiv,1,mpi_prec,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lcfl,gcfl,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lshr,gshr,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lvor,gvor,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lfou,gfou,1,mpi_prec,MPI_MAX,0,comm,info)
  CALL MPI_Reduce(lrot,grot,1,mpi_prec,MPI_SUM,0,comm,info)
  !JTR IS THIS redundant?
  CALL MPI_Reduce(lrcirl_d,grcirl_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lrcirr_d,grcirr_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lrcir25l_d,grcir25l_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lrcir25r_d,grcir25r_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lrcir26l_d,grcir26l_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lrcir26r_d,grcir26r_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  IF(md.NE.mk) THEN
     CALL MPI_Reduce(lens_d,gens_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
     CALL MPI_Reduce(leng_d,geng_d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  ELSE
     gens_d = gens
     geng_d = geng
  END IF
  
  ldu = nx(3)
  CALL MPI_Reduce(lengofx_d,gengofx_d,ldu,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lensofx_d,gensofx_d,ldu,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lcirofx_d,gcirofx_d,ldu,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lxomegaofx_d,gxomegaofx_d,ldu,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)
  CALL MPI_Reduce(lyomegaofx_d,gyomegaofx_d,ldu,MPI_DOUBLE_PRECISION,MPI_SUM,0,comm,info)

  ! adapt time step -- this is stupid
  ! -
  !
  ! WRITE(msg,*) 'gshr = ',gshr
  !   CALL ppm_write(rank,'wvic_diagnostics',msg,info)
  !   IF(SQRT(gshr).GT.0.0_mk) THEN
  !      newdt = MIN(dt_max,0.25_mk/SQRT(gvor)*0.5_mk)
  !      newdt = MIN(newdt,0.48_mk*dx**2/(nu*3.0_mk))
  !      newdt = MIN(newdt,0.5_mk/SQRT(gshr))
  !   END IF
  !   ! find minimal dt for the whole domain
  !   CALL MPI_AllReduce(newdt,dt,1,mpi_prec,MPI_MIN,comm,info)
  ! max_vorticity = SQRT(gvor)
  !
  IF (dt_adapt) THEN
     IF(rank.EQ.0) THEN
        IF(SQRT(gshr).GT.0.0_mk) THEN
           IF(gvor .NE. 0.0_mk) THEN
              newdt = MIN(dt_max,0.25_mk/SQRT(gvor)*0.5_mk)
           END IF
           newdt = MIN(newdt,0.48_mk*dx**2/(nu*3.0_mk))
           IF(gshr .NE. 0.0_mk) THEN
              newdt = MIN(newdt,0.5_mk/SQRT(gshr))
           END IF
        END IF
      END IF
      CALL MPI_BCast(newdt,1,mpi_prec,0,comm,info)
      dt = newdt
      CALL MPI_BCast(max_vorticity,1,mpi_prec,0,comm,info)
  END IF
  
  IF(rank.EQ.0) THEN
     WRITE(filename,'(A,A)') runtag(1:iruntag),'-diag.dat'
     OPEN(14,file=filename,iostat=ios,position='append',status='unknown')
     ! timestep, time, circ, rcirc, energy, enstrophy, dt, maxvorticity, rot
     WRITE(14,'(I16,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,&
      A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,&
      E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5)') &
          & itime,' ',time,' ',gcir,' ',grcirr_d,' ',grcirl_d,' ',&
          & grcir25r_d,' ',grcir25l_d,' ',&
          & grcir26r_d,' ',grcir26l_d,' ',&
          &geng_d,' ',gens_d,' ',&
          & dt,' ',SQRT(gvor),' ',grot,' ',u_cmass(1),' ',u_cmass(2),' ',u_cmass(3)
     CLOSE(14)
     IF (MOD(itime,ndump).EQ.0) THEN 
         DO k=1,ldu-1
            gxomegaofx_d(k) = gxomegaofx_d(k) / gcirofx_d(k)
            gyomegaofx_d(k) = gyomegaofx_d(k) / gcirofx_d(k)
         END DO
         WRITE(filename,'(A,A,I5.5,A)') runtag(1:iruntag),'-spatialdiagI',itime,'.dat'
         OPEN(15,file=filename,iostat=ios,status='unknown',action='WRITE')
         ldu = nx(3)
         WRITE(15,*)  (/ (i*dz,i=0,ldu-1) /)
         WRITE(15,*)  gengofx_d(1:ldu)
         WRITE(15,*)  gensofx_d(1:ldu)
         WRITE(15,*)  gcirofx_d(1:ldu)
         WRITE(15,*)  gxomegaofx_d(1:ldu)
         WRITE(15,*)  gyomegaofx_d(1:ldu)
         CLOSE(15,iostat=info)
      END IF
  END IF

  !-----------------------------------------------------------------------------
  !  Deallocation
  !-----------------------------------------------------------------------------
  DEALLOCATE(lengofx_d,gengofx_d,stat=info)
  DEALLOCATE(lensofx_d,gensofx_d,stat=info)
  DEALLOCATE(lcirofx_d,gcirofx_d,stat=info)
  DEALLOCATE(lxomegaofx_d,gxomegaofx_d,stat=info)
  DEALLOCATE(lyomegaofx_d,gyomegaofx_d,stat=info)

2002 CONTINUE
  CALL MPI_BCast(max_vorticity,1,mpi_prec,0,comm,info)


END SUBROUTINE wvic_diagnostics_trail
  
  
