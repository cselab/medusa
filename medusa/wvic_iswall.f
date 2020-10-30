!------------------------------------------------------------------------------!
!* filename: wvic_iswall                                                      *!
!* project : ppm                                                              *!
!* purpose : generate an array with flags saying if the boundary              *!
!*           of a subdomain is a wall                                         *!
!*         :                                                                  *!
!* author  : Philippe Chatelain / Michael Bergdorf                            *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Aug 18 10:36:19 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_iswall.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.3  2005/11/21 17:36:38  michaebe
!  major session. now let''s debug
!
!  Revision 1.2  2005/10/09 13:07:32  michaebe
!  corrections (nsub versus nsublist) and cosmetics
!
!  Revision 1.1  2005/10/07 14:33:54  pchatela
!  Generation of "Wall boundary conditions vs subdomain" flags
!  Enforcement of wall boundary conditions a la Thom
!  Added entries in Makefile
!
!  Revision 1.1  2005/09/28 11:40:44  michaebe
!  Fork from ppm_pvc
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
!  Generate iswall arrays
!------------------------------------------------------------------------------!
  SUBROUTINE wvic_iswall (info)

  USE module_wvic
  !--- ppm general
  USE ppm_module_data
  USE ppm_module_data_mesh
  USE ppm_module_write
  USE ppm_module_alloc
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
  
  INTEGER                                :: topoidint, meshidint
  !INTEGER,  DIMENSION(:,:)     , POINTER :: istart, ndata
  INTEGER                                :: nsubs, isub, d
  INTEGER,  DIMENSION(3)                 :: ldu,ldl
  INTEGER                                :: iopt
  
  !----------------------------------------------------------------------------!
  !  Alloc memory for the iswall array
  !----------------------------------------------------------------------------!
  iopt   = ppm_param_alloc_fit
  ldu(1) = 2
  ldu(2) = dime
  ldu(3) = nsublist
  CALL ppm_alloc(iswall,ldu,iopt,info)
  IF (info .NE. 0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'wvic_iswall',     &
          &        'iswall flags ISWALL',__LINE__,info)
     GOTO 9999
  ENDIF
  !----------------------------------------------------------------------------!
  !  Kill in time
  !----------------------------------------------------------------------------!
  iopt   = ppm_param_alloc_fit
  ldu(1) = 6
  ldu(2) = nsublist
  CALL ppm_alloc(liswall,ldu,iopt,info)
  IF (info .NE. 0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'wvic_iswall',     &
          &        'iswall flags ISWALL',__LINE__,info)
     GOTO 9999
  ENDIF
  !----------------------------------------------------------------------------!
  ! Loop over the subdomains of this proc
  !----------------------------------------------------------------------------!
  DO isub = 1,nsublist
     !-------------------------------------------------------------------------!
     ! Loop over the dimensions
     !-------------------------------------------------------------------------!
     DO d=1,dime
        iswall(:,d,isub) = 0
        IF ((istart(d,isub).EQ.1).AND.(bcdef((d-1)*2 + 1).EQ.ppm_param_bcdef_dirichlet)) THEN
           iswall(1,d,isub) = 1
        ENDIF
        IF (((istart(d,isub)+ndata(d,isublist(isub))-1).EQ.Nx(dime)).AND.(bcdef((d-1)*2 + 2).EQ.ppm_param_bcdef_dirichlet)) THEN
           iswall(2,d,isub) = 1
        ENDIF
        liswall(2*d-1,isub) = (iswall(1,d,isub).EQ.1)
        liswall(2*d  ,isub) = (iswall(2,d,isub).EQ.1)
     END DO
  END DO
  
  !----------------------------------------------------------------------------!
  !  Return 
  !----------------------------------------------------------------------------!
9999 CONTINUE
  
     RETURN

  END SUBROUTINE wvic_iswall
        
