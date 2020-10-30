!------------------------------------------------------------------------------!
!* filename: wvic_iswall                                                      *!
!* project : ppm                                                              *!
!* purpose : Enforce the no-slip condition for stationary walls               *!
!*           This routines assumes that the iswall array has been             *!
!*           initialized                                                      *!
!*         :                                                                  *!
!* author  : Philippe Chatelain / Michael Bergdorf                            *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Aug 18 10:36:19 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_enforcewall.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.8  2005/11/24 08:49:41  pchatela
!  Bugfix: wrong sign for the factx-y-z in Thom's formula
!
!  Revision 1.7  2005/11/21 17:36:37  michaebe
!  major session. now let''s debug
!
!  Revision 1.6  2005/11/15 11:03:19  michaebe
!  csmtcs
!
!  Revision 1.5  2005/11/14 22:45:34  pchatela
!  Added the contributions of the potentianl field
!  Added extrapolation of velocity field to velocity field bcs
!  * linear for components parallel to wall, could be smarter and use
!  wall vorticity value...
!  * quadratic for normal component
!
!  Revision 1.4  2005/10/14 09:03:00  michaebe
!  implemented enforcewall_velocity
!
!  Revision 1.3  2005/10/09 13:38:00  michaebe
!  renamed routine, added one for velocity (to be implemented)
!
!  Revision 1.2  2005/10/09 12:49:00  michaebe
!  cosmetics
!
!  Revision 1.1  2005/10/07 14:33:53  pchatela
!  Generation of "Wall boundary conditions vs subdomain" flags
!  Enforcement of wall boundary conditions a la Thom
!  Added entries in Makefile
!
!  Revision 1.1  2005/09/28 11:40:44  michaebe
!  Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!  Enforce the no slip condition using a Thom formula...
!------------------------------------------------------------------------------!
SUBROUTINE wvic_enforcewall_vorticity (info)

  USE module_wvic
  !--- ppm general
  USE ppm_module_data
  USE ppm_module_data_mesh
  USE ppm_module_write
  USE ppm_module_error
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
  IMPLICIT NONE
  
  !----------------------------------------------------------------------------!
  !  Arguments
  !----------------------------------------------------------------------------!
  INTEGER,  INTENT(INOUT)                :: info
  
  INTEGER                                :: i,j,k,topoidint, meshidint
  !INTEGER,  DIMENSION(:,:)     , POINTER :: istart, ndata
  INTEGER                                :: nsubs, isub, isubl, top, topm1
  INTEGER,  DIMENSION(3)                 :: ldu,ldl
  INTEGER,  DIMENSION(3)                 :: Nm
  REAL(MK)                               :: factx,facty,factz
  REAL(MK)                               :: potx, poty, potz

#ifdef __TOPOUNKNOWN__
  !----------------------------------------------------------------------------!
  !  Get the internal topoid
  !----------------------------------------------------------------------------!
  topoidint = ppm_internal_topoid(topo_id)
  !----------------------------------------------------------------------------!
  !  Get the internal meshid
  !----------------------------------------------------------------------------!
  meshidint = ppm_meshid(topoidint)%internal(mesh_id)
  !----------------------------------------------------------------------------!
  !  Get istart
  !----------------------------------------------------------------------------!
  istart => ppm_cart_mesh(meshidint,topoidint)%istart
  
  !----------------------------------------------------------------------------!
  !  Assignment of the useful arrays/scalar
  !----------------------------------------------------------------------------!
  Nm(1:dime) = ppm_cart_mesh(meshidint,topoidint)%Nm
  nsubs = ppm_nsublist(topoidint)
#else
  !----------------------------------------------------------------------------!
  !  All the above information is already known...
  !----------------------------------------------------------------------------!
#endif

  

  
  !----------------------------------------------------------------------------!
  ! Initialize the fators
  !----------------------------------------------------------------------------!
  factx = 2.0/dx**2
  facty = 2.0/dy**2
  factz = 2.0/dz**2
  
  !----------------------------------------------------------------------------!
  ! Loop over the subdomains of this proc
  !----------------------------------------------------------------------------!
  DO isub = 1,nsublist
     isubl = isublist(isub)
     !-------------------------------------------------------------------------!
     ! Dimension 1
     !-------------------------------------------------------------------------!
     IF (iswall(1,1,isub).EQ.1) THEN
	poty = 2.0*u_infty(3)/dx
	potz = 2.0*u_infty(2)/dx
        DO k=1,ndata(3,isubl)
           DO j=1,ndata(2,isubl)
              field_wp(1,1,j,k,isub) = 0.0_MK
              field_wp(2,1,j,k,isub) = -factx*field_wps(2,2,j,k,isub) - poty	
              field_wp(3,1,j,k,isub) = -factx*field_wps(3,2,j,k,isub) + potz
           END DO
        END DO
     END IF
     IF (iswall(2,1,isub).EQ.1) THEN
	poty = 2.0*u_infty(3)/dx
	potz = 2.0*u_infty(2)/dx
        top = ndata(1,isubl)
        topm1 = top-1
        DO k=1,ndata(3,isubl)
           DO j=1,ndata(2,isubl)
              field_wp(1,top,j,k,isub) = 0.0_MK
              field_wp(2,top,j,k,isub) = -factx*field_wps(2,topm1,j,k,isub) &
                   & + poty
              field_wp(3,top,j,k,isub) = -factx*field_wps(3,topm1,j,k,isub) &
                   & - potz
           END DO
        END DO
     END IF
     !-------------------------------------------------------------------------!
     ! Dimension 2
     !-------------------------------------------------------------------------!
     IF (iswall(1,2,isub).EQ.1) THEN
		potx = 2.0*u_infty(3)/dy
		potz = 2.0*u_infty(1)/dy
        DO k=1,ndata(3,isubl)
           DO i=1,ndata(1,isubl)
              field_wp(1,i,1,k,isub) = -facty*field_wps(1,i,2,k,isub) &
                   & + potx
              field_wp(2,i,1,k,isub) = 0.0_MK
              field_wp(3,i,1,k,isub) = -facty*field_wps(3,i,2,k,isub) &
                   & - potz
           END DO
        END DO
     END IF
     IF (iswall(2,2,isub).EQ.1) THEN
        potx = 2.0*u_infty(3)/dy
        potz = 2.0*u_infty(1)/dy
        top = ndata(2,isubl)
        topm1 = top-1
        DO k=1,ndata(3,isubl)
           DO i=1,ndata(1,isubl)
              field_wp(1,i,top,k,isub) = -facty*field_wps(1,i,topm1,k,isub) &
                   & - potx
              field_wp(2,i,top,k,isub) = 0.0_MK
              field_wp(3,i,top,k,isub) = -facty*field_wps(3,i,topm1,k,isub) &
                   & + potz
           END DO
        END DO
     END IF
     !-------------------------------------------------------------------------!
     ! Dimension 3
     !-------------------------------------------------------------------------!
     IF (iswall(1,3,isub).EQ.1) THEN
		potx = 2.0*u_infty(2)/dz
		poty = 2.0*u_infty(1)/dz
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_wp(1,i,j,1,isub) = -factz*field_wps(1,i,j,2,isub) - potx
              field_wp(2,i,j,1,isub) = -factz*field_wps(2,i,j,2,isub) + poty
              field_wp(3,i,j,1,isub) = 0.0_MK
           END DO
        END DO
     END IF
     IF (iswall(2,3,isub).EQ.1) THEN
        potx = 2.0*u_infty(2)/dz
        poty = 2.0*u_infty(1)/dz
        top = ndata(3,isubl)
        topm1 = top-1
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_wp(1,i,j,top,isub) = -factz*field_wps(1,i,j,topm1,isub) &
                   & + potx
              field_wp(2,i,j,top,isub) = -factz*field_wps(2,i,j,topm1,isub) &
                   & - poty
              field_wp(3,i,j,top,isub) = 0.0_MK
           END DO
        END DO
     END IF
  END DO
  
  !----------------------------------------------------------------------------!
  !  Return 
  !----------------------------------------------------------------------------!
9999 CONTINUE
  
  RETURN
  
END SUBROUTINE wvic_enforcewall_vorticity






!------------------------------------------------------------------------------!
!  Enforce boundary condition on velocity
!------------------------------------------------------------------------------!
SUBROUTINE wvic_enforcewall_velocity (info)

  USE module_wvic
  IMPLICIT NONE
  
  INTEGER,  INTENT(INOUT)                :: info
  
  INTEGER                                :: i,j,k
  !INTEGER,  DIMENSION(:,:)     , POINTER :: istart, ndata
  INTEGER                                :: nsubs, isub, isubl, top, topm1,topp1
  INTEGER,  DIMENSION(3)                 :: ldu,ldl
  INTEGER,  DIMENSION(3)                 :: Nm
  REAL(MK)                               :: factx,facty,factz

  !----------------------------------------------------------------------------!
  ! Loop over the subdomains of this proc
  !----------------------------------------------------------------------------!
  DO isub = 1,nsublist
     isubl = isublist(isub)
     !-------------------------------------------------------------------------!
     ! Dimension 1
     !-------------------------------------------------------------------------!
     IF (iswall(1,1,isub).EQ.1) THEN
        DO k=1,ndata(3,isubl)
           DO j=1,ndata(2,isubl)
              field_up(1,1,j,k,isub) = 0.0_MK
              field_up(2,1,j,k,isub) = 0.0_MK
              field_up(3,1,j,k,isub) = 0.0_MK
              field_up(1,0,j,k,isub) =  field_up(1,2,j,k,isub)
              field_up(2,0,j,k,isub) = -field_up(2,2,j,k,isub)
              field_up(3,0,j,k,isub) = -field_up(3,2,j,k,isub)
           END DO
        END DO
     END IF
     IF (iswall(2,1,isub).EQ.1) THEN
        top = ndata(1,isubl)
        topm1 = top-1
	topp1 = top+1
        DO k=1,ndata(3,isubl)
           DO j=1,ndata(2,isubl)
              field_up(1,top,j,k,isub) = 0.0_MK
              field_up(2,top,j,k,isub) = 0.0_MK
              field_up(3,top,j,k,isub) = 0.0_MK
              field_up(1,topp1,j,k,isub) =  field_up(1,topm1,j,k,isub)
              field_up(2,topp1,j,k,isub) = -field_up(2,topm1,j,k,isub)
              field_up(3,topp1,j,k,isub) = -field_up(3,topm1,j,k,isub)
           END DO
        END DO
     END IF
     !-------------------------------------------------------------------------!
     ! Dimension 2
     !-------------------------------------------------------------------------!
     IF (iswall(1,2,isub).EQ.1) THEN
        DO k=1,ndata(3,isubl)
           DO i=1,ndata(1,isubl)
              field_up(1,i,1,k,isub) = 0.0_MK
              field_up(2,i,1,k,isub) = 0.0_MK
              field_up(3,i,1,k,isub) = 0.0_MK
              field_up(1,i,0,k,isub) = -field_up(1,i,2,k,isub)
              field_up(2,i,0,k,isub) =  field_up(2,i,2,k,isub)
              field_up(3,i,0,k,isub) = -field_up(3,i,2,k,isub)
           END DO
        END DO
     END IF
     IF (iswall(2,2,isub).EQ.1) THEN
        top = ndata(2,isubl)
        topm1 = top-1
	topp1 = top+1
        DO k=1,ndata(3,isubl)
           DO i=1,ndata(1,isubl)
              field_up(1,i,top,k,isub) = 0.0_MK
              field_up(2,i,top,k,isub) = 0.0_MK
              field_up(3,i,top,k,isub) = 0.0_MK
              field_up(1,i,topp1,k,isub) = -field_up(1,i,topm1,k,isub)
              field_up(2,i,topp1,k,isub) =  field_up(2,i,topm1,k,isub)
              field_up(3,i,topp1,k,isub) = -field_up(3,i,topm1,k,isub)
           END DO
        END DO
     END IF
     !-------------------------------------------------------------------------!
     ! Dimension 3
     !-------------------------------------------------------------------------!
     IF (iswall(1,3,isub).EQ.1) THEN
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_up(1,i,j,1,isub) = 0.0_MK
              field_up(2,i,j,1,isub) = 0.0_MK
              field_up(3,i,j,1,isub) = 0.0_MK
              field_up(1,i,j,0,isub) = -field_up(1,i,j,2,isub)
              field_up(2,i,j,0,isub) = -field_up(2,i,j,2,isub)
              field_up(3,i,j,0,isub) =  field_up(3,i,j,2,isub)
           END DO
        END DO
     END IF
     IF (iswall(2,3,isub).EQ.1) THEN
        top = ndata(3,isubl)
        topm1 = top-1
        topp1 = top+1
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              field_up(1,i,j,top,isub) = 0.0_MK
              field_up(2,i,j,top,isub) = 0.0_MK
              field_up(3,i,j,top,isub) = 0.0_MK
              field_up(1,i,j,topp1,isub) = -field_up(1,i,j,topm1,isub)
              field_up(2,i,j,topp1,isub) = -field_up(2,i,j,topm1,isub)
              field_up(3,i,j,topp1,isub) =  field_up(3,i,j,topm1,isub)
           END DO
        END DO
     END IF
  END DO
  
  !----------------------------------------------------------------------------!
  !  Return 
  !----------------------------------------------------------------------------!
9999 CONTINUE
  
  RETURN

END SUBROUTINE wvic_enforcewall_velocity
