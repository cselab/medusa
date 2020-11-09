#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_mg_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine initializes the solver for 
      !                 2D and 3D problems
      !
      !  Input        :  equation   (I)  :  KIND OF EQUATION TO BE SOLVED 
      !                                     FOR THE MOMENT ONLY POISSON
      !                  order      (I)  :  ORDER OF FINITE DIFFERENCES
      !                                     NOW SECOND. THE GHOSTSIZE IS 
      !                                     AUTOMATICALLY ADJUSTED 
      !                  smoother   (I)  :  NOW GAUSS-SEIDEL
      !
      !                  [lda]     (I)   : LEADING DIMENSION, ONLY TO BE
      !                                    GIVEN FOR VECTOR CASES
      !                
      !                  ibcdef     (I)  : ARRAY OF BOUNDARY CONDITION 
      !
      !
      !                  bcvalue   (F)   : ARRAY WHERE THE VALUES OF THE BC
      !                                    ARE STORED.IN CASE OF PERIODIC 
      !                                    JUST GIVE ANY KIND OF VALUE
      !
      !                  EPSU      (F)   : STOPPING CRITERIUM. DETAIL:SHOULD
      !                                    BE SCALED WITH THE MAXIMUM VALUE           !                                    OF THE RHS.
      !                  
      !                  limlev    (I)    :Number of levels that the user 
      !                                    wants to coarse.
      !                
      !                  wcycle    (L)    : TRUE if the user wants W-cycle.
      !                                    OTHERWISE FALSE
      !                  lprint    (L)    : TRUE IF YOU WANT TO DUMP OUT
      !                                     INFORMATION
      !                  
      !                   omega     (F)    : relaxation parameter for SOR
      !
      !  
      !  Input/output :     
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Remarks      :  PLEASE PAY ATTENTION THAT IN ORDER TO DIVIDE 
      !                  FURTHER A MESH IT SHOULD BE DIVISIBLE WITH 2.
      !                  IF YOU WANT TO SOLVE DIFFERENT EQUATIONS 
      !                  THE WHOLE MACHINERY SHOULD BE CALLED TWICE.
      !                  ALSO THE SOLVER IS NOW PROGRAMMED FOR THE POISSON
      !                  PROBLEM. A FUTURE IMPROVEMENT WOULD BE
      !                  TO USE A GENERAL STENCIL.      
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mg_init.f,v $
      !  Revision 1.2  2006/08/22 15:54:37  pchatela
      !  Added a hopefully appropriate scaling factor in the comparisons against
      !  lmyeps
      !
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.7  2006/02/03 09:34:03  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.6  2005/12/08 12:43:16  kotsalie
      !  commiting dirichlet
      !
      !  Revision 1.5  2005/01/04 09:47:45  kotsalie
      !  ghostsize=2 for scalar case
      !
      !  Revision 1.4  2004/10/29 15:59:14  kotsalie
      !  RED BLACK SOR FOR 3d vec case. 2d will soon follow.
      !
      !  Revision 1.3  2004/09/28 14:04:49  kotsalie
      !  Changes concerning 4th order finite differences
      !
      !  Revision 1.2  2004/09/23 09:38:30  kotsalie
      !  added details in the header
      !
      !  Revision 1.1  2004/09/22 18:27:09  kotsalie
      !  MG new version
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if __DIM == __SFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_sca_s(equation,iorder,smoother,ibcdef,&
     &          bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_sca_d(equation,iorder,smoother,ibcdef,&
     &                          bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#elif  __MESH_DIM  == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_sca_s(equation,iorder,smoother,ibcdef,&
     &                         bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_sca_d(equation,iorder,smoother,ibcdef,&
     &                         bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_vec_s(equation,iorder,smoother,lda,ibcdef,&
     &    bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_2d_vec_d(equation,iorder,smoother,lda,ibcdef,&
     &   bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#elif  __MESH_DIM  == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_vec_s(equation,iorder,smoother,lda,ibcdef,&
     &              bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_init_3d_vec_d(equation,iorder,smoother,lda,ibcdef,&
     &                    bcvalue,EPSU,mesh_id,limlev,wcycle,lprint,omega,info)
#endif
#endif
#endif

        !----------------------------------------------------------------------
        !  Includes
        !----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  Modules 
        !-----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_data_mesh
        USE ppm_module_data_mg
        USE ppm_module_alloc
        USE ppm_module_mg_alloc
        USE ppm_module_error
        USE ppm_module_mesh_derive
        USE ppm_module_substart
        USE ppm_module_substop 

#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !-----------------------------------------------------------------------  
        !  Arguments     
        !-----------------------------------------------------------------------
        INTEGER, INTENT(IN)                                :: equation
        INTEGER, INTENT(IN)                                :: iorder  
        INTEGER, INTENT(IN)                                :: smoother
#if __DIM == __VFIELD
        INTEGER,              INTENT(IN)                   ::  lda  
#endif
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
        INTEGER,DIMENSION(:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:)                            ::  bcvalue
#elif __MESH_DIM == __3D
        INTEGER,DIMENSION(:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:,:)                          ::  bcvalue
#endif  
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
        INTEGER,DIMENSION(:,:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:,:)                            ::  bcvalue
#elif __MESH_DIM == __3D
        INTEGER,DIMENSION(:,:)                               ::  ibcdef
        REAL(MK),DIMENSION(:,:,:,:)                          ::  bcvalue
#endif
#endif

        INTEGER,  INTENT(IN)                               :: mesh_id
        REAL(MK),INTENT(IN)                                :: EPSU
        INTEGER,INTENT(IN)                                 :: limlev
        LOGICAL,INTENT(IN)                                 :: wcycle
        LOGICAL,INTENT(IN)                                 :: lprint
        REAL(MK),INTENT(IN)                                :: omega
        INTEGER, INTENT(OUT)                               :: info
        !--------------------------------------------------------------------
        !  Local variables 
        !-----------------------------------------------------------------------
        REAL(MK)                             :: t0
        REAL(MK)                             :: lmyeps
        INTEGER                              :: meshid,mlev,isub 
        INTEGER                              :: idom
        INTEGER                              ::  count,ilda,iface
        INTEGER                              :: i,j,k
        INTEGER                              :: kk
#if __MESH_DIM == __2D
        INTEGER                              :: dir
#endif
        INTEGER                              :: iter1,iter2,ix,iy
        INTEGER                              :: ipoint,jpoint 
        INTEGER                              :: newmeshid,lmesh_id
        INTEGER , DIMENSION(1)               :: ldu1
        INTEGER , DIMENSION(2)               :: ldu2,ldl2 ,direc
        INTEGER , DIMENSION(3)               :: ldu3,ldl3 
#if __MESH_DIM == __3D
        INTEGER                              :: dir1,dir2,jj,iz
        INTEGER , DIMENSION(4)               :: ldu4,ldl4
#endif
        INTEGER , DIMENSION(ppm_dim)         :: Nml 
        REAL(MK), DIMENSION(ppm_dim)         :: min_phys,max_phys
        REAL(MK), DIMENSION(ppm_dim,ppm_nsubs(ppm_field_topoid)) &
             &                               :: min_sub,max_sub
        INTEGER                              :: iopt,topoid


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_2d_sca_s),DIMENSION(:,:),POINTER :: mgfield => null()
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_2d_sca_d),DIMENSION(:,:),POINTER :: mgfield => null()
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_3d_sca_s),DIMENSION(:,:),POINTER :: mgfield => null()
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_3d_sca_d),DIMENSION(:,:),POINTER :: mgfield => null()
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_2d_vec_s),DIMENSION(:,:),POINTER :: mgfield => null()
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_2d_vec_d),DIMENSION(:,:),POINTER :: mgfield => null()
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_3d_vec_s),DIMENSION(:,:),POINTER :: mgfield => null()
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_3d_vec_d),DIMENSION(:,:),POINTER :: mgfield => null()
#endif
#endif
#endif



#if __DIM == __SFIELD
#if __MESH_DIM == __2D
REAL(MK),DIMENSION(:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:),POINTER :: terr
#elif __MESH_DIM == __3D
REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:,:),POINTER :: terr
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:,:),POINTER :: terr
#elif __MESH_DIM == __3D
REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
REAL(MK),DIMENSION(:,:,:,:),POINTER :: terr
#endif
#endif


        !-----------------------------------------------------------------------       
        !  Externals 
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  Initialize 
        !-----------------------------------------------------------------------

        CALL substart('ppm_mg_init',t0,info)


        !-----------------------------------------------------------------------  
        !  Check arguments
        !-----------------------------------------------------------------------
          IF (ppm_debug.GT.0) THEN
#if __DIM == __VFIELD
          IF (lda.LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_poiss_mg_init',  &
     &            'lda must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
          IF (EPSU.LE.0.0_MK) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_poiss_mg_init',  &
    &            'EPSU must be >0',__LINE__,info)
             GOTO 9999
          ENDIF
        ENDIF

        !---------------------------------------------------------------------
        ! Definition of necessary variables and allocation of arrays
        !---------------------------------------------------------------------
#if __DIM == __SFIELD
        vecdim = 1
#elif __DIM == __VFIELD
        vecdim = lda
#endif
        w_cycle=wcycle
        l_print=lprint
                                                   
        topoid = ppm_field_topoid 
        nsubs  = ppm_nsublist(topoid) 
        !PRINT *,'nsub:',nsubs 
        meshid = ppm_meshid(topoid)%internal(mesh_id)
        lmesh_id = mesh_id

                                           
#if    __KIND == __SINGLE_PRECISION
        min_phys(:)=ppm_min_physs(:,topoid)
        max_phys(:)=ppm_max_physs(:,topoid)
        min_sub(:,:)=ppm_min_subs(:,:,topoid)
        max_sub(:,:)=ppm_max_subs(:,:,topoid)
        EPSU_s = EPSU
        omega_s=omega
        lmyeps=ppm_myepss 
#elif  __KIND == __DOUBLE_PRECISION
        min_phys(:)=ppm_min_physd(:,topoid)
        max_phys(:)=ppm_max_physd(:,topoid)
        min_sub(:,:)=ppm_min_subd(:,:,topoid)
        max_sub(:,:)=ppm_max_subd(:,:,topoid)
        EPSU_d = EPSU
        omega_d=omega
        lmyeps=ppm_myepsd 
#endif
#if __MESH_DIM == __2D
        Nml(1) = ppm_cart_mesh(meshid,topoid)%Nm(1) 
        Nml(2) = ppm_cart_mesh(meshid,topoid)%Nm(2) 
        maxlev = INT(log10(Nml(1)*Nml(2)*REAL(ppm_nproc,MK))/log10(2.0_MK))
        IF (maxlev.GT.limlev) THEN
         maxlev=limlev
        ENDIF 
#if __KIND == __SINGLE_PRECISION
        dx_s = (max_phys(1)-min_phys(1))/(Nml(1)-1) 
        dy_s = (max_phys(2)-min_phys(2))/(Nml(2)-1) 
        rdx2_s  = 1/(dx_s*dx_s)
        rdy2_s  = 1/(dy_s*dy_s) 
#elif __KIND == __DOUBLE_PRECISION
        dx_d = (max_phys(1)-min_phys(1))/(Nml(1)-1) 
        dy_d = (max_phys(2)-min_phys(2))/(Nml(2)-1) 

        rdx2_d  = 1/(dx_d*dx_d)
        rdy2_d  = 1/(dy_d*dy_d) 

#endif
#elif __MESH_DIM == __3D
        Nml(1) = ppm_cart_mesh(meshid,topoid)%Nm(1) 
        Nml(2) = ppm_cart_mesh(meshid,topoid)%Nm(2) 
        Nml(3) = ppm_cart_mesh(meshid,topoid)%Nm(3) 
        maxlev = INT(log10(Nml(1)*Nml(2)*Nml(3)* &
     &           REAL(ppm_nproc,MK))/log10(2.0_MK))

        IF (maxlev.GT.limlev) THEN
         maxlev=limlev
        ENDIF 
#if __KIND == __SINGLE_PRECISION
        dx_s = (max_phys(1)-min_phys(1))/(Nml(1)-1) 
        dy_s = (max_phys(2)-min_phys(2))/(Nml(2)-1) 
        dz_s = (max_phys(3)-min_phys(3))/(Nml(3)-1) 
        rdx2_s = 1/(dx_s*dx_s)
        rdy2_s = 1/(dy_s*dy_s) 
        rdz2_s = 1/(dz_s*dz_s)
#elif __KIND == __DOUBLE_PRECISION
        dx_d = (max_phys(1)-min_phys(1))/(Nml(1)-1) 
        dy_d = (max_phys(2)-min_phys(2))/(Nml(2)-1) 
        dz_d = (max_phys(3)-min_phys(3))/(Nml(3)-1) 
        rdx2_d = 1/(dx_d*dx_d)
        rdy2_d = 1/(dy_d*dy_d) 
        rdz2_d = 1/(dz_d*dz_d)
#endif
#endif


#if __DIM == __SFIELD
!Print *,'SField'
        iopt = ppm_param_alloc_fit    
        ldu2(1) = nsubs
        ldu2(2) = 2*ppm_dim
        CALL ppm_alloc(bcdef_sca,ldu2,iopt,info)
        IF (info .NE. 0) THEN 
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'Boundary condiotions',__LINE__,info)
           GOTO 9999
        ENDIF
        bcdef_sca(:,:)=0

        !---------------------------------------------------------
        !MICHAEL
        !--------------------------------------------------------
        DO isub=1,nsubs 
         idom=ppm_isublist(isub,topoid)
         !----------------------------------------------------------------------
         !  compare the west boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(1,idom)-min_phys(1)) .LT. lmyeps*(max_sub(1,i)-min_sub(1,i))) THEN
            bcdef_sca(isub,1)=ibcdef(1)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the east boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(1,idom)-max_phys(1)) .LT. lmyeps*(max_sub(1,i)-min_sub(1,i))) THEN
            bcdef_sca(isub,2)=ibcdef(2)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the south boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(2,idom)-min_phys(2)) .LT. lmyeps*(max_sub(2,i)-min_sub(2,i))) THEN
            bcdef_sca(isub,3)=ibcdef(3)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the north boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(2,idom)-max_phys(2)) .LT. lmyeps*(max_sub(2,i)-min_sub(2,i))) THEN
            bcdef_sca(isub,4)=ibcdef(4)
         ENDIF
	 !----------------------------------------------------------------------
         !  compare the south boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(3,idom)-min_phys(3)) .LT. lmyeps*(max_sub(3,i)-min_sub(3,i))) THEN
            bcdef_sca(isub,5)=ibcdef(5)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the north boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(3,idom)-max_phys(3)) .LT. lmyeps*(max_sub(3,i)-min_sub(3,i))) THEN
            bcdef_sca(isub,6)=ibcdef(6)
         ENDIF
        ENDDO  
        lperiodic=.TRUE.  
        
        DO isub=1,nsubs  
         DO i=1,2*ppm_dim
          IF (bcdef_sca(isub,i).NE.ppm_param_bcdef_periodic) THEN
           lperiodic=.FALSE.  
           EXIT  
          ENDIF 
         ENDDO  
        ENDDO 

#elif __DIM == __VFIELD
!print *,'Vfiedl'
        iopt = ppm_param_alloc_fit
        ldu3(1) = vecdim
        ldu3(2) = nsubs
        ldu3(3) = 2*ppm_dim
        CALL ppm_alloc(bcdef_vec,ldu3,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'Boundary condiotions',__LINE__,info)
           GOTO 9999
        ENDIF

	bcdef_vec(:,:,:)=0
	Do isub=1,nsubs
	  idom=ppm_isublist(isub,topoid)
	  Do ilda=1,vecdim
	    
	    !----------------------------------------------------------------------
         !  compare the west boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(1,idom)-min_phys(1)) .LT. lmyeps*(max_sub(1,i)-min_sub(1,i))) THEN
            bcdef_vec(ilda,isub,1)=ibcdef(ilda,1)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the east boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(1,idom)-max_phys(1)) .LT. lmyeps*(max_sub(1,i)-min_sub(1,i))) THEN
            bcdef_vec(ilda,isub,2)=ibcdef(ilda,2)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the south boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(2,idom)-min_phys(2)) .LT. lmyeps*(max_sub(2,i)-min_sub(2,i))) THEN
            bcdef_vec(ilda,isub,3)=ibcdef(ilda,3)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the north boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(2,idom)-max_phys(2)) .LT. lmyeps*(max_sub(2,i)-min_sub(2,i))) THEN
            bcdef_vec(ilda,isub,4)=ibcdef(ilda,4)
         ENDIF
	 !----------------------------------------------------------------------
         !  compare the south boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(3,idom)-min_phys(3)) .LT. lmyeps*(max_sub(3,i)-min_sub(3,i))) THEN
            bcdef_vec(ilda,isub,5)=ibcdef(ilda,5)
         ENDIF

         !----------------------------------------------------------------------
         !  compare the north boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(3,idom)-max_phys(3)) .LT. lmyeps*(max_sub(3,i)-min_sub(3,i))) THEN
            bcdef_vec(ilda,isub,6)=ibcdef(ilda,6)
         ENDIF
	enddo
	enddo
	lperiodic=.TRUE.
	Do isub=1,nsubs
	  DO i=1,2*ppm_dim
	   DO ilda=1,vecdim
	    IF (bcdef_vec(ilda,isub,i).NE.ppm_param_bcdef_periodic) Then
		lperiodic=.FALSE.
		EXIT
	    ENDIF
	   ENDDO
	  ENDDO
	ENDDO
#endif


        iopt = ppm_param_alloc_fit    
        ldu1(1) = ppm_dim
        CALL ppm_alloc(ghostsize,ldu1,iopt,info)
        IF (info .NE. 0) THEN 
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'ghostsize',__LINE__,info)
           GOTO 9999
        ENDIF

        IF (iorder.EQ.ppm_param_order_2) THEN
         ghostsize(:)=1
         order=iorder
        ELSEIF (iorder.EQ.ppm_param_order_4) THEN
         ghostsize(:)=2
         order=ppm_param_order_4
        ENDIF 

        iopt = ppm_param_alloc_fit    
        ldu1(1) = ppm_dim
        CALL ppm_alloc(factor,ldu1,iopt,info)
        IF (info .NE. 0) THEN 
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &                   'factor',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit    
        ldu1(1) = maxlev
        CALL ppm_alloc(meshid_g,ldu1,iopt,info)
        IF (info .NE. 0) THEN 
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                  'meshid_g',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit    
        ldu1(1) = maxlev
        CALL ppm_alloc(mesh_id_g,ldu1,iopt,info)
        IF (info .NE. 0) THEN 
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
      &                  'mesh_id_g',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit
        ldu3(1) = ppm_dim
        ldu3(2) = nsubs
        ldu3(3) = maxlev
        CALL ppm_alloc(start,ldu3,iopt,info)
        IF (info .NE. 0) THEN   
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',  &
     &             'starting indices when updating the field',__LINE__,info)
           GOTO 9999
        ENDIF


        iopt = ppm_param_alloc_fit
        ldu3(1) = ppm_dim
        ldu3(2) = nsubs
        ldu3(3) = maxlev
        CALL ppm_alloc(stop,ldu3,iopt,info)
        IF (info .NE. 0) THEN   
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'stopping indices when updating the field',__LINE__,info)
           GOTO 9999
        ENDIF


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_sca_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_2d_sca_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_sca_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF
        mgfield => mgfield_2d_sca_d
#endif

#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_sca_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
                &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_sca_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_sca_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_sca_d 

#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_vec_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_2d_vec_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_2d_vec_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF
        mgfield => mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_vec_s,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
                &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_vec_s

#elif __KIND == __DOUBLE_PRECISION
        iopt = ppm_param_alloc_fit
        ldu2(1) = nsubs
        ldu2(2) = maxlev
        CALL ppm_mg_alloc(mgfield_3d_vec_d,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',   &
     &        'Multigrid fields used on the different levels',__LINE__,info)
           GOTO 9999
        ENDIF

        mgfield => mgfield_3d_vec_d 

#endif
#endif
#endif



        iopt = ppm_param_alloc_fit
        ldu2(1) = 2*ppm_dim
        ldu2(2) = nsubs
        CALL ppm_alloc(lboundary,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the boundary alloc.',__LINE__,info)
           GOTO 9999
        ENDIF

        iopt = ppm_param_alloc_fit
        ldu2(1) = ppm_dim
        ldu2(2) = maxlev
        CALL ppm_alloc(max_node,ldu2,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with a maximum number alloc.',__LINE__,info)
           GOTO 9999
        ENDIF
        max_node(:,:)=0

        lboundary(:,:)=.FALSE.
        start(:,:,:)=1


        !-----------------------------------------------------------------------
        ! Derive coarser meshes 
        !-----------------------------------------------------------------------

        DO mlev=1,maxlev
            
          
#if __MESH_DIM == __2D

           !--------------------------------------------------------------------
           ! Go through the subs, define the stopping indices on each mesh,
           ! check and store if it is on the boundary, allocate the 
           ! multigrid fields, pass the boundary values.
           !--------------------------------------------------------------------
           DO i=1,nsubs
              idom=ppm_isublist(i,topoid)

               stop(:,i,mlev)= ppm_cart_mesh(meshid,topoid)%nnodes(:,idom)

              DO j=1,ppm_dim
                 IF (max_node(j,mlev).LT.stop(j,i,mlev)) THEN
                    max_node(j,mlev)=stop(j,i,mlev)  
                 ENDIF
              ENDDO


              !-----------------------------------------------------------------
              ! Allocate the function correction, the restricted errors,
              ! the residuals and the values on the boundary on each level.
              !----------------------------------------------------------------
#if __DIM == __SFIELD
              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF            
                 
              tuc=>mgfield(i,mlev)%uc
              tuc(:,:)=0.0_MK 
              
              iopt = ppm_param_alloc_fit
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
         &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

               mgfield(i,mlev)%fc(:,:)=0.0_MK

              iopt = ppm_param_alloc_fit
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)

              CALL ppm_alloc(mgfield(i,mlev)%err,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              terr=>mgfield(i,mlev)%err  
              terr(:,:)=0.0_MK
           

              !------------------------------------------------------------------
              !MICHAEL
              !------------------------------------------------------------------ 
              !ALLOCATE THE BCVALUE(IT IS A TYPE!!)
             !PRINT *,'LPERIODIC:',lperiodic
             IF (.NOT.lperiodic) THEN
              iopt = ppm_param_alloc_fit
              ldu1(1) = 2*ppm_dim
              CALL ppm_mg_alloc(mgfield(i,mlev)%bcvalue,ldu1,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
         
              !ALLOCATE THE PBCVALUE 
            
             DO iface=1,2*ppm_dim 
              iopt = ppm_param_alloc_fit
              IF (iface.EQ.1.OR.iface.EQ.2) THEN
               ldu1(1) = max_node(2,mlev)
              ELSE
               ldu1(1) = max_node(1,mlev)
              ENDIF
	      
	      
              CALL ppm_alloc(mgfield(i,mlev)%bcvalue(iface)%pbcvalue,ldu1,iopt,info)
	      IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF 
             ENDDO


             DO iface=1,2*ppm_dim
               IF (iface.EQ.1.OR.iface.EQ.2) THEN 
	       	direc(1)=2
		elseif (iface.EQ.3.OR.iface.EQ.4) then
	  	direc(1)=1
		endif
                 DO ipoint=1,max_node(direc(1),mlev)
              
                  IF (mlev.EQ.1) THEN                         
		     mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint)=bcvalue(iface,ipoint)
		          ELSE
                   IF(bcdef_sca(i,iface).EQ.ppm_param_bcdef_neumann) THEN 
                    mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint)=&
               &            mgfield(i,mlev-1)%bcvalue(iface)%pbcvalue(2*ipoint-1) 
                   ELSE
                   !NO CORRECTIONS FOR THE DIRICHLET  
                   ! mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint)=0.0_MK
		  
			         mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint)=0.0_MK
                   ENDIF  
                  ENDIF
                 ENDDO
               \
             ENDDO!faces 
         ENDIF!lperiodic
#elif __DIM == __VFIELD

              iopt = ppm_param_alloc_fit
              ldl3(1) = 1
              ldl3(2) = 1-ghostsize(1)
              ldl3(3) = 1-ghostsize(2)
              ldu3(1) = vecdim
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
             
                
              tuc=>mgfield(i,mlev)%uc
              tuc(:,:,:)=0.0_MK

              iopt = ppm_param_alloc_fit
              ldu3(1) = vecdim  
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

               mgfield(i,mlev)%fc(:,:,:)=0.0_MK
	 	
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1
              ldl3(2) = 1-ghostsize(1)
              ldl3(3) = 1-ghostsize(2)
              ldu3(1) = vecdim
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              terr=>mgfield(i,mlev)%err  
              terr(:,:,:)=0.0_MK
	      
#endif
              iopt = ppm_param_alloc_fit  
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%mask_red,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the mask  alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              iopt = ppm_param_alloc_fit  
              ldl2(1) = 1-ghostsize(1)
              ldl2(2) = 1-ghostsize(2)
              ldu2(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu2(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              CALL ppm_alloc(mgfield(i,mlev)%mask_black,ldl2,ldu2,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with mask alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              !---------------------------------------------------------------- 
              !Filling the mask for communication (red black) 
              !----------------------------------------------------------------
              DO iy=1-ghostsize(2),ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
                 DO ix=1-ghostsize(1),ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)

                    IF (MOD(ix+iy,2).EQ.0) THEN  

                       mgfield(i,mlev)%mask_red(ix,iy)=.TRUE.      
                       mgfield(i,mlev)%mask_black(ix,iy)=.FALSE.      

                    ELSE 

                       mgfield(i,mlev)%mask_red(ix,iy)   = .FALSE.      
                       mgfield(i,mlev)%mask_black(ix,iy) = .TRUE.      

                    ENDIF
                 ENDDO
              ENDDO

           ENDDO!DO 1,nsubs



#elif __MESH_DIM == __3D 


           DO i=1,nsubs

              idom=ppm_isublist(i,topoid)
              stop(:,i,mlev) = ppm_cart_mesh(meshid,topoid)%nnodes(:,idom)

              DO j=1,ppm_dim
                 IF (max_node(j,mlev).LT.stop(j,i,mlev)) THEN
                    max_node(j,mlev)=stop(j,i,mlev)  
                 ENDIF
              ENDDO

              IF (ppm_subs_bc(1,idom,topoid).EQ.1) THEN

                 lboundary(1,i)=.TRUE.          


              ELSEIF (ppm_subs_bc(3,idom,topoid).EQ.1) THEN

                 lboundary(3,i)=.TRUE.

              ELSEIF (ppm_subs_bc(2,idom,topoid).EQ.1) THEN

                 lboundary(2,i)=.TRUE.  


              ELSEIF (ppm_subs_bc(4,idom,topoid).EQ.1) THEN

                 lboundary(4,i)=.TRUE.

              ELSEIF (ppm_subs_bc(5,idom,topoid).EQ.1) THEN

                 lboundary(5,i)=.TRUE. 


              ELSEIF (ppm_subs_bc(6,idom,topoid).EQ.1) THEN

                 lboundary(6,i)=.TRUE.


              ENDIF


              !----------------------------------------------------------------
              ! Allocate the function correction, the restricted errors and the 
              !residuals on each level.
              !----------------------------------------------------------------

#if __DIM == __SFIELD
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
                
              tuc=>mgfield(i,mlev)%uc
              tuc(:,:,:)=0.0_MK              

              iopt = ppm_param_alloc_fit
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              mgfield(i,mlev)%fc(:,:,:)=0.0_MK

              
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
	      terr=>mgfield(i,mlev)%err  
              terr(:,:,:)=0.0_MK 
	      
	      !ALLOCATE THE BCVALUE(IT IS A TYPE!!)
             !PRINT *,'LPERIODIC:',lperiodic
             IF (.NOT.lperiodic) THEN
              iopt = ppm_param_alloc_fit
              ldu1(1) = 2*ppm_dim
              CALL ppm_mg_alloc(mgfield(i,mlev)%bcvalue,ldu1,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
	       !ALLOCATE THE PBCVALUE 
             
	      
             DO iface=1,2*ppm_dim 
              iopt = ppm_param_alloc_fit
              IF (iface.EQ.1.OR.iface.EQ.2) THEN
               ldu2(1) = max_node(2,mlev)
	       ldu2(2)= max_node(3,mlev)
              ELSEif (iface.EQ.3.OR. iface.EQ.4) then	      
               ldu2(1) = max_node(1,mlev)
	       ldu2(2)=max_node(3,mlev)
	      else
	       ldu2(1)=max_node(1,mlev)
	       ldu2(2)=max_node(2,mlev)
              ENDIF
	      
              CALL ppm_alloc(mgfield(i,mlev)%bcvalue(iface)%pbcvalue,ldu2,iopt,info)
	      !Print *,size(mgfield(i,mlev)%bcvalue(iface)%pbcvalue,1)
	      IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF 
             ENDDO
		
              
	      DO iface=1,2*ppm_dim
               IF (iface.EQ.1.OR.iface.EQ.2) THEN  
	       		direc(1)=2
			direc(2)=3
		elseif (iface.EQ.3.OR.iface.EQ.4) THEN
			direc(1)=1
			direc(2)=3
		else
			direc(1)=1
			direc(2)=2
		endif
                 DO ipoint=1,max_node(direc(1),mlev)
                  DO jpoint=1,max_node(direc(2),mlev)
                   IF (mlev.EQ.1) THEN                         
                      mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint)=bcvalue(iface,ipoint,jpoint)
                      

		   ELSE
                    IF(bcdef_sca(i,iface).EQ.ppm_param_bcdef_neumann) THEN 
                           mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint)=&
               &           mgfield(i,mlev-1)%bcvalue(iface)%pbcvalue(2*ipoint-1,2*jpoint-1) 
                    ELSE
                   !NO CORRECTIONS FOR THE DIRICHLET  
                   
		
			mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint)=0.0_MK
!			If (mlev.EQ.5) then
!				Print *,ipoint,jpoint,mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint)
!			endif
                   ENDIF  
                  ENDIF
                 ENDDO
		enddo 	       		          
             ENDDO!faces 
         endif !lperiodic
#elif __DIM == __VFIELD

              iopt = ppm_param_alloc_fit
              ldl4(1) = 1
              ldl4(2) = 1-ghostsize(1)
              ldl4(3) = 1-ghostsize(2)
              ldl4(4) = 1-ghostsize(3)
              ldu4(1) = vecdim
              ldu4(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu4(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu4(4) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%uc,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the function corr. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
                
              tuc=>mgfield(i,mlev)%uc
              tuc(:,:,:,:)=0.0_MK              
 

              iopt = ppm_param_alloc_fit
              ldu4(1) = vecdim
              ldu4(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)
              ldu4(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)
              ldu4(4) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)
              CALL ppm_alloc(mgfield(i,mlev)%fc,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the restricted err. alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              mgfield(i,mlev)%fc(:,:,:,:)=0.0_MK


              iopt = ppm_param_alloc_fit
              ldl4(1) = 1
              ldl4(2) = 1-ghostsize(1)
              ldl4(3) = 1-ghostsize(2)
              ldl4(4) = 1-ghostsize(3)
              ldu4(1) = vecdim
              ldu4(2) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu4(3) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu4(4) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%err,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the residual alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              terr=>mgfield(i,mlev)%err  
              terr(:,:,:,:)=0.0_MK
	      
	       !ALLOCATE THE BCVALUE(IT IS A TYPE!!)
             !PRINT *,'LPERIODIC:',lperiodic
             IF (.NOT.lperiodic) THEN
              iopt = ppm_param_alloc_fit
	      ldu1=2*ppm_dim
	      !ldu2(1)=vecdim
              !ldu2(2) = 2*ppm_dim

	      !allocate(mgfield(i,mlev)%bcvalue(3,6))
              CALL ppm_mg_alloc(mgfield(i,mlev)%bcvalue,ldu1,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF
	       !ALLOCATE THE PBCVALUE 

             DO iface=1,2*ppm_dim 
              iopt = ppm_param_alloc_fit
	      ldu3(1)=vecdim
              IF (iface.EQ.1.OR.iface.EQ.2) THEN
	       
               ldu3(2) = max_node(2,mlev)
	       ldu3(3)= max_node(3,mlev)
              ELSEif (iface.EQ.3.OR. iface.EQ.4) then	      
               ldu3(2) = max_node(1,mlev)
	       ldu3(3)=max_node(3,mlev)
	      else
	       ldu3(2)=max_node(1,mlev)
	       ldu3(3)=max_node(2,mlev)
              ENDIF
	
	      CALL ppm_alloc(mgfield(i,mlev)%bcvalue(iface)%pbcvalue,ldu3,iopt,info)


	      IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the BOUNDARY alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF 
             ENDDO
	      	     
	      DO iface=1,2*ppm_dim
                 IF (iface.EQ.1.OR.iface.EQ.2) THEN  
	       	     	direc(1)=2
		         	direc(2)=3
	             elseif (iface.EQ.3.OR.iface.EQ.4) THEN
			        direc(1)=1
			        direc(2)=3
		         else
			        direc(1)=1
			        direc(2)=2
		         endif

               DO ipoint=1,max_node(direc(1),mlev)
                  DO jpoint=1,max_node(direc(2),mlev)
		           DO ilda=1,vecdim

                   IF (mlev.EQ.1) THEN                         
                       !PRINT *,'bcef',ilda,iface,bcvalue(ilda,iface,ipoint,jpoint)	
			        mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ilda,ipoint,jpoint)=bcvalue(ilda,iface,ipoint,jpoint)
            !PRINT *,ipoint,jpoint,iface,ilda,bcvalue(ilda,iface,ipoint,jpoint),mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ilda,ipoint,jpoint)
            
	!	Print*,'Boundariu',bcvalue(ilda,iface,ipoint,jpoint)	
		           ELSE	    
                    IF(bcdef_vec(ilda,i,iface).EQ.ppm_param_bcdef_neumann) THEN 
                           mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ilda,ipoint,jpoint)=&
               &           mgfield(i,mlev-1)%bcvalue(iface)%pbcvalue(ilda,2*ipoint-1,2*jpoint-1) 
                    ELSE
                   !NO CORRECTIONS FOR THE DIRICHLET  

		
		        	mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ilda,ipoint,jpoint)=0.0_MK
!			If (mlev.EQ.5) then
!				Print *,ipoint,jpoint,mgfield(i,mlev)%bcvalue(iface)%pbcvalue(ipoint,jpoint)
!			endif
                   ENDIF  
                  ENDIF
		 ENDDO
                 ENDDO
		enddo 	       		          
             ENDDO!faces 
	   ENDIF !lperiodic
	      
#endif


              iopt = ppm_param_alloc_fit  
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%mask_red,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
     &        'Problem with the mask  alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF

              iopt = ppm_param_alloc_fit  
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1-ghostsize(3)
              ldu3(1) = ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)
              ldu3(2) = ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
              ldu3(3) = ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
              CALL ppm_alloc(mgfield(i,mlev)%mask_black,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_poiss_mg_init',    &
                      &        'Problem with mask alloc.',__LINE__,info)
                 GOTO 9999
              ENDIF


              !----------------------------------------------------------------
              !Filling the mask for communication (red black) 
              !-----------------------------------------------------------------
              DO iz=1-ghostsize(3),&
               ppm_cart_mesh(meshid,topoid)%nnodes(3,idom)+ghostsize(3)
                                      
                 DO iy=1-ghostsize(2),&
                   & ppm_cart_mesh(meshid,topoid)%nnodes(2,idom)+ghostsize(2)
                    DO ix=1-ghostsize(1),&
                   &  ppm_cart_mesh(meshid,topoid)%nnodes(1,idom)+ghostsize(1)

                       IF (MOD(ix+iy+iz,2).EQ.0) THEN  

                          mgfield(i,mlev)%mask_red(ix,iy,iz)=.TRUE.      
                          mgfield(i,mlev)%mask_black(ix,iy,iz)=.FALSE.      
                          !mgfield(i,mlev)%mask_black(ix,iy,iz)=.TRUE.      

                       ELSE 

                          mgfield(i,mlev)%mask_red(ix,iy,iz)   = .FALSE.      
                          mgfield(i,mlev)%mask_black(ix,iy,iz) = .TRUE.      

                      ENDIF
                    ENDDO
                 ENDDO
              ENDDO



           ENDDO!DO i=1,nsubs

#endif


           factor(:)=2
           mesh_id_g(mlev)=lmesh_id
           meshid_g(mlev)=meshid
           newmeshid=-1

           IF (mlev.LT.maxlev) THEN 
	    !Print *,'dfj',meshid,ppm_param_mesh_coarsen,factor,newmeshid
            CALL ppm_mesh_derive(topoid,meshid,ppm_param_mesh_coarsen,factor,&
     &                          newmeshid,info)

           
            lmesh_id = newmeshid
            meshid = ppm_meshid(topoid)%internal(lmesh_id)

           ENDIF 
	   
        ENDDO!DO mlev=1,maxlev 
	   	      

        !----------------------------------------------------------------------
        !  Return 
        !----------------------------------------------------------------------
9999    CONTINUE
        CALL substop('ppm_mg_init',t0,info)
        RETURN
#if    __DIM       == __SFIELD
#if    __MESH_DIM  == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_sca_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if    __MESH_DIM  == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_2d_vec_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_init_3d_vec_d
#endif
#endif
#endif
