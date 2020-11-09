      !-------------------------------------------------------------------------
      !     Subroutine   :              ppm_hamjac_ext_step_3d
      !-------------------------------------------------------------------------
      !     
      !     Purpose      : Extension
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
      !     $Log: ppm_hamjac_ext_step_3d.f,v $
      !     Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !     initial import
      !
      !     Revision 1.3  2006/03/27 08:17:33  michaebe
      !     bug fix
      !
      !     Revision 1.2  2005/08/12 14:38:00  ivos
      !     bugfix: index bounds in loop corrected.
      !
      !     Revision 1.1  2005/07/25 00:34:01  ivos
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
      SUBROUTINE ppm_hamjac_ext_step_3ds (phi, psi, tpsi, res, &
           &                          topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_ext_step_3dd (phi, psi, tpsi, res, &
           &                          topo_id, mesh_id, ghostsize, info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_hamjac_ext_step_3dsv (phi, psi, lda, tpsi, res, &
           &                          topo_id, mesh_id, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_hamjac_ext_step_3ddv (phi, psi, lda, tpsi, res, &
           &                          topo_id, mesh_id, ghostsize, info)
#endif
#endif

        USE ppm_module_data
        USE ppm_module_data_mesh
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_substop
        
        
#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION       
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

        !-----------------------------------------------------
        !  Arguments
        !-----------------------------------------------------
        REAL(MK), DIMENSION(:,:,:,:),   POINTER :: phi
#if   __MODE == __SCA
        REAL(mk), DIMENSION(:,:,:,:  ), POINTER :: psi, tpsi
#elif __MODE == __VEC
        REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: psi, tpsi
        INTEGER,  INTENT(in)                    :: lda
#endif
        INTEGER, INTENT(in)                   :: topo_id, mesh_id
        INTEGER, DIMENSION(3), INTENT(in)     :: ghostsize
        INTEGER, INTENT(inout)                :: info
        REAL(mk),INTENT(out)                  :: res

        !-----------------------------------------------------
        !  Aliases
        !-----------------------------------------------------
        INTEGER, DIMENSION(:),    POINTER     :: isublist
        INTEGER                               :: nsublist
        INTEGER, DIMENSION(:,:),  POINTER     :: ndata
        INTEGER                               :: topoid, meshid
        REAL(mk), DIMENSION(:,:), POINTER     :: min_phys, max_phys
        
        !-----------------------------------------------------
        !  standard stuff
        !-----------------------------------------------------
        INTEGER                               :: isub,isubl,i,j,k
        REAL(mk)                              :: len_phys(3)
        !-----------------------------------------------------
        !  WENO stuff
        !-----------------------------------------------------
        REAL(mk) :: oneg(3), opos(3), wenoeps, wenotau, pbs, n(3)
        REAL(mk) :: laps(-1:1,3), rpos(3), rneg(3), dx(3), dxi(3)
        REAL(mk) :: phip(3), phin(3), phimid(3), rms, sij
#if    __MODE == __SCA
        REAL(MK) :: dphi_dt
#elif  __MODE == __VEC
        REAL(mk) :: dphi_dt(10)
#endif        
        INTEGER  :: ilap
        INTEGER, PARAMETER, DIMENSION(3,3) :: offs &
             & = RESHAPE((/2,1,0,1,0,-1,0,-1,-2/),(/3,3/))
        REAL(mk) :: t0, dxavg

        dphi_dt = 0.0_mk
        CALL substart('ppm_hamjac_ext_step_3d',t0,info)
        
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

        len_phys(1) = max_phys(1,topoid) - min_phys(1,topoid)
        len_phys(2) = max_phys(2,topoid) - min_phys(2,topoid)
        len_phys(3) = max_phys(3,topoid) - min_phys(3,topoid)
        dx(1)       = len_phys(1)/REAL(ppm_cart_mesh(meshid,topoid)%nm(1)-1,mk)
        dx(2)       = len_phys(2)/REAL(ppm_cart_mesh(meshid,topoid)%nm(2)-1,mk)
        dx(3)       = len_phys(3)/REAL(ppm_cart_mesh(meshid,topoid)%nm(3)-1,mk)
        dxavg       = SUM(dx(1:3))/3.0_mk
        dxi(1)      = 1.0_mk/dx(1)
        dxi(2)      = 1.0_mk/dx(2)
        dxi(3)      = 1.0_mk/dx(3)
        wenoeps = 1.0e-6_mk
        wenotau = 0.5_mk*MINVAL(dx)
        rms = -HUGE(rms)

        DO isub=1,nsublist
           isubl = isublist(isub)
           DO k=1,ndata(3,isubl)
              DO j=1,ndata(2,isubl)
                 DO i=1,ndata(1,isubl)

                    IF(ABS(phi(i,j,k,isub)).GT.7.0_mk*dx(1)) CYCLE
                    IF(ABS(phi(i,j,k,isub)).LT.dx(1)) THEN
#if   __MODE == __SCA
                       tpsi(i,j,k,isub) = psi(i,j,k,isub)
#elif __MODE == __VEC
                       tpsi(1:lda,i,j,k,isub) = psi(1:lda,i,j,k,isub)
#endif
                    ELSE
                    phimid(1) = phi(i+1,j,k,isub)-phi(i-1,j,k,isub)
                    phimid(2) = phi(i,j+1,k,isub)-phi(i,j-1,k,isub)
                    phimid(3) = phi(i,j,k+1,isub)-phi(i,j,k-1,isub)
                    sij       = phi(i,j,k,isub) &
                         &       /SQRT(phi(i,j,k,isub)**2+dxavg**2)
                    n         = phimid / SQRT(SUM(phimid**2))
#if   __MODE == __SCA
                    dphi_dt    = &
                         & MAX(n(1)*sij,0.0_mk)*dxi(1)*          &
                         & (psi(i,j,k,isub)-psi(i-1,j,k,isub)) + &
                         & MIN(n(1)*sij,0.0_mk)*dxi(1)*          &
                         & (psi(i+1,j,k,isub)-psi(i,j,k,isub)) + &
                         & MAX(n(2)*sij,0.0_mk)*dxi(2)*          &
                         & (psi(i,j,k,isub)-psi(i,j-1,k,isub)) + &
                         & MIN(n(2)*sij,0.0_mk)*dxi(2)*          &
                         & (psi(i,j+1,k,isub)-psi(i,j,k,isub)) + &
                         & MAX(n(3)*sij,0.0_mk)*dxi(3)*          &
                         & (psi(i,j,k,isub)-psi(i,j,k-1,isub)) + &
                         & MIN(n(3)*sij,0.0_mk)*dxi(3)*          &
                         & (psi(i,j,k+1,isub)-psi(i,j,k,isub)) 
                    
                    tpsi(i,j,k,isub) = psi(i,j,k,isub) - wenotau * dphi_dt
                    rms = MAX(rms,ABS(dphi_dt))
#elif __MODE == __VEC
                    dphi_dt(1:lda)    = &
                         & MAX(n(1)*sij,0.0_mk)*dxi(1)*          &
                         & (psi(1:lda,i,j,k,isub)-psi(1:lda,i-1,j,k,isub)) + &
                         & MIN(n(1)*sij,0.0_mk)*dxi(1)*          &
                         & (psi(1:lda,i+1,j,k,isub)-psi(1:lda,i,j,k,isub)) + &
                         & MAX(n(2)*sij,0.0_mk)*dxi(2)*          &
                         & (psi(1:lda,i,j,k,isub)-psi(1:lda,i,j-1,k,isub)) + &
                         & MIN(n(2)*sij,0.0_mk)*dxi(2)*          &
                         & (psi(1:lda,i,j+1,k,isub)-psi(1:lda,i,j,k,isub)) + &
                         & MAX(n(3)*sij,0.0_mk)*dxi(3)*          &
                         & (psi(1:lda,i,j,k,isub)-psi(1:lda,i,j,k-1,isub)) + &
                         & MIN(n(3)*sij,0.0_mk)*dxi(3)*          &
                         & (psi(1:lda,i,j,k+1,isub)-psi(1:lda,i,j,k,isub)) 
                    
                    tpsi(1:lda,i,j,k,isub) = psi(1:lda,i,j,k,isub) &
                         & - wenotau * dphi_dt(1:lda)

                    rms = MAX(rms,SUM(ABS(dphi_dt)))
#endif                    
                    END IF
                    
                    

                 END DO

              END DO

           END DO

        END DO

        res = rms

        CALL substop('ppm_hamjac_ext_step_3d',t0,info)
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_step_3ds 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_step_3dd 
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_step_3dsv 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_hamjac_ext_step_3ddv 
#endif
#endif
      
                    


                    
           
           



        
        
        
        
