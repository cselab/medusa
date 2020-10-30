      !-------------------------------------------------------------------------
      !     Subroutine   :                ppm_hamjac_ext_3d
      !-------------------------------------------------------------------------
      !     
      !     Purpose      : Solve Hamilton-Jacobi for Gowas extension
      !      
      !     Input        : 
      !                    
      !     Input/Output : 
      !                    
      !     Output       : 
      !      
      !     Remarks      : 
      !                    
      !     
      !     References   :
      !     
      !     Revisions    :
      !-------------------------------------------------------------------------
      !     $Log: ppm_hamjac_ext_3d.f,v $
      !     Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !     initial import
      !
      !     Revision 1.2  2005/08/25 16:48:50  ivos
      !     Fixed format string. pgf90 barked.
      !
      !     Revision 1.1  2005/07/25 00:34:00  ivos
      !     Initial check-in.
      !
      !-------------------------------------------------------------------------
      !     Parallel Particle Mesh Library (PPM)
      !     Institute of Computational Science
      !     ETH Zentrum, Hirschengraben 84
      !     CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------


#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_hamjac_ext_3ds (phi, psi, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_ext_3dd (phi, psi, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_hamjac_ext_3dsv (phi, psi, lda, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_ext_3ddv (phi, psi, lda, tol, maxstep, &
           &                     topo_id, mesh_id, ghostsize, info)
#endif
#endif

        USE ppm_module_data
        USE ppm_module_data_mesh
        USE ppm_module_error
        USE ppm_module_write
        USE ppm_module_substart
        USE ppm_module_alloc
        USE ppm_module_substop
        USE ppm_module_map
        IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION       
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

        !-----------------------------------------------------
        !  Arguments
        !-----------------------------------------------------
        REAL(MK), DIMENSION(:,:,:,:), POINTER :: phi
#if   __MODE == __SCA
        REAL(mk), DIMENSION(:,:,:,:), POINTER :: psi
#elif __MODE == __VEC
        REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: psi
        INTEGER, INTENT(in)                   :: lda
#endif        
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(3), INTENT(in)     :: ghostsize
        INTEGER, INTENT(out)                  :: info
        INTEGER, INTENT(in)                   :: maxstep
        REAL(mk), INTENT(in)                  :: tol

        !-----------------------------------------------------
        !  Aliases
        !-----------------------------------------------------
        INTEGER, DIMENSION(:), POINTER        :: isublist
#if   __MODE == __SCA
        REAL(mk), DIMENSION(:,:,:,:), POINTER :: tpsi
#elif __MODE == __VEC
        REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: tpsi
#endif        
        INTEGER                               :: nsublist
        INTEGER, DIMENSION(:,:), POINTER      :: ndata
        INTEGER                               :: topoid,meshid
        REAL(mk), DIMENSION(:,:), POINTER     :: min_phys, max_phys
        
        !-----------------------------------------------------
        !  Local variables
        !-----------------------------------------------------
        INTEGER                               :: isub,isubl,i,j,k,maptype
        INTEGER                               :: istep,iopt
#if   __MODE == __SCA
        INTEGER                               :: ldl(4), ldu(4)
#elif __MODE == __VEC        
        INTEGER                               :: ldl(5), ldu(5)
#endif
        INTEGER                               :: ndata_max(3)
        REAL(mk)                              :: t0, res
        CHARACTER(LEN=ppm_char)               :: cbuf
        REAL(mk)                              :: dx(3), len_phys(3)
        !-----------------------------------------------------
        !  Initialisation
        !-----------------------------------------------------
        CALL substart('ppm_hamjac_ext_3d',t0,info)
        
        !-----------------------------------------------------
        !  Get the mesh data
        !-----------------------------------------------------
        topoid = ppm_internal_topoid(topo_id)
        meshid = ppm_meshid(topoid)%internal(mesh_id)
        nsublist = ppm_nsublist(topoid)
        ndata    => ppm_cart_mesh(meshid,topoid)%nnodes
        !  COMMENT Thu May 26 19:39:51 PDT 2005:  experimental
        isublist => ppm_isublist(:,topoid)
#if    __KIND == __SINGLE_PRECISION
        min_phys => ppm_min_physs
        max_phys => ppm_max_physs
#elif  __KIND == __DOUBLE_PRECISION       
        min_phys => ppm_min_physd
        max_phys => ppm_max_physd
#endif


        !-----------------------------------------------------
        !  RATIONALE Thu May 26 20:51:19 PDT 2005:
        !  loop ghostmap doit. easy.
        !-----------------------------------------------------


        !-----------------------------------------------------
        !  allocate temporary storage
        !-----------------------------------------------------
#if __MODE == __SCA
        ldl(1:3) = 1 - ghostsize(1:3); ldl(4) = 1
#elif __MODE == __VEC
        ldl(1) = 1
        ldl(2:4) = 1- ghostsize(1:3); ldl(5) = 1
#endif        
        
        ndata_max(1) = MAXVAL(ndata(1,1:nsublist))
        ndata_max(2) = MAXVAL(ndata(2,1:nsublist))
        ndata_max(3) = MAXVAL(ndata(3,1:nsublist))
#if __MODE == __SCA
        ldu(1)   = ndata_max(1) + ghostsize(1)
        ldu(2)   = ndata_max(2) + ghostsize(2)
        ldu(3)   = ndata_max(3) + ghostsize(3)
        ldu(4)   = nsublist
#elif __MODE == __VEC
        ldu(1)   = lda
        ldu(2)   = ndata_max(1) + ghostsize(1)
        ldu(3)   = ndata_max(2) + ghostsize(2)
        ldu(4)   = ndata_max(3) + ghostsize(3)
        ldu(5)   = nsublist
#endif
        
        iopt     = ppm_param_alloc_fit
        CALL ppm_alloc(tpsi,ldl,ldu,iopt,info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_hamjac_ext_3d', &
                &        'temp storage for donowatocalit',__LINE__,info)
           GOTO 9999
        END IF

        len_phys(1) = max_phys(1,topoid) - min_phys(1,topoid)
        len_phys(2) = max_phys(2,topoid) - min_phys(2,topoid)
        len_phys(3) = max_phys(3,topoid) - min_phys(3,topoid)
        dx(1)       = len_phys(1)/REAL(ppm_cart_mesh(meshid,topoid)%nm(1)-1,mk)
        dx(2)       = len_phys(2)/REAL(ppm_cart_mesh(meshid,topoid)%nm(2)-1,mk)
        dx(3)       = len_phys(3)/REAL(ppm_cart_mesh(meshid,topoid)%nm(3)-1,mk)
        
        
        !--- ready to blast
        maptype = ppm_param_map_init
        CALL ppm_map_field_ghost(phi,topo_id,mesh_id,ghostsize,maptype,info)

        !--- COMMENT Thu May 26 21:05:23 PDT 2005:  simple euler here, do TVD
        !--- map the gowas
        maptype = ppm_param_map_ghost_get
        CALL ppm_map_field_ghost(phi,topo_id,mesh_id,ghostsize,maptype,info)
        maptype = ppm_param_map_push
        CALL ppm_map_field_ghost(phi,topo_id,mesh_id,ghostsize,maptype,info)
        maptype = ppm_param_map_send
        CALL ppm_map_field_ghost(phi,topo_id,mesh_id,ghostsize,maptype,info)
        maptype = ppm_param_map_pop
        CALL ppm_map_field_ghost(phi,topo_id,mesh_id,ghostsize,maptype,info)
           !--- map the function
        DO istep=1,maxstep
#if   __MODE == __SCA
           maptype = ppm_param_map_ghost_get
           CALL ppm_map_field_ghost(psi,topo_id,mesh_id,ghostsize,maptype,info)
           maptype = ppm_param_map_push
           CALL ppm_map_field_ghost(psi,topo_id,mesh_id,ghostsize,maptype,info)
           maptype = ppm_param_map_send
           CALL ppm_map_field_ghost(psi,topo_id,mesh_id,ghostsize,maptype,info)
           maptype = ppm_param_map_pop
           CALL ppm_map_field_ghost(psi,topo_id,mesh_id,ghostsize,maptype,info)
#elif __MODE == __VEC
           maptype = ppm_param_map_ghost_get
           CALL ppm_map_field_ghost(psi,lda,topo_id,mesh_id,ghostsize,maptype,info)
           maptype = ppm_param_map_push
           CALL ppm_map_field_ghost(psi,lda,topo_id,mesh_id,ghostsize,maptype,info)
           maptype = ppm_param_map_send
           CALL ppm_map_field_ghost(psi,lda,topo_id,mesh_id,ghostsize,maptype,info)
           maptype = ppm_param_map_pop
           CALL ppm_map_field_ghost(psi,lda,topo_id,mesh_id,ghostsize,maptype,info)
#endif           
           !     IF (ppm_debug .GT. 0) THEN
           !     ENDIF
#if __MODE == __SCA
           CALL ppm_hamjac_ext_step (phi,psi,tpsi,res,topo_id,mesh_id&
                &,                  ghostsize,info)
#elif __MODE == __VEC
           CALL ppm_hamjac_ext_step (phi,psi,lda,tpsi,res,topo_id,mesh_id&
                &,                  ghostsize,info)
#endif
           IF(ppm_debug.GT.0) THEN
              WRITE(cbuf,'(A,I4,A,E12.5)') 'Iteration ',istep,' Residual: ',res
              CALL ppm_write(ppm_rank,'ppm_hamjac_ext_3d',cbuf,info)
           END IF
           DO isub=1,nsublist
              isubl = isublist(isub)
              DO k=1,ndata(3,isubl);DO j=1,ndata(2,isubl);DO i=1,ndata(1,isubl)
                 IF(ABS(phi(i,j,k,isub)).GT.6.0_mk*dx(1)) CYCLE
#if __MODE == __SCA
                 psi(i,j,k,isub) = tpsi(i,j,k,isub)
#elif __MODE == __VEC
                 psi(1:lda,i,j,k,isub) = tpsi(1:lda,i,j,k,isub)
#endif                 
              END DO; END DO; END DO
           END DO
           IF(res.LT.tol) GOTO 666
        END DO

        info = ppm_error_warning
        CALL ppm_error(ppm_err_converge,'ppm_hamjac_ext_3d', &
             &         'failed to reach target residual',__LINE__,info)

666     CONTINUE

        iopt = ppm_param_dealloc
        CALL ppm_alloc(tpsi,ldl,ldu,iopt,info)
        IF(info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_hamjac_ext_3d', &
                &        'temp storage for donowatocalit not freed',__LINE__,info)
           GOTO 9999
        END IF


9999    CONTINUE

#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_3ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_3dd 
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_3dsv
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_3ddv
#endif
#endif


        
           

        
        

