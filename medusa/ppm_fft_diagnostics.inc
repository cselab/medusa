      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fft_diagnostics
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !  Input        : 
      !  Input/output : 
      !  Output       : info       (I) return status. =0 if no error.
      !                   
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fft_diagnostics.F,v $
      !  Revision 1.4  2006/10/23 08:19:11  pchatela
      !  LES models
      !  Bugfixes in KE spectra and factor for Parseval identity
      !  Removed the reset of noise in init_tvphysics_0 and_1
      !
      !  Revision 1.3  2006/10/03 16:11:58  pchatela
      !  Added spectra in the z direction (misnamed kx spectrum...)
      !  Added spatial diagnostics, like kinetic energy, enstrophy, circulation
      !  as functions of z, dumped at the frequency ndump
      !
      !  Revision 1.2  2006/09/27 09:30:21  pchatela
      !  Fixes, spectra calculation,
      !  most importantly: moved the u_infty out, so it does not kill the dgammadt
      !
      !  Revision 1.1  2006/09/16 00:22:02  pchatela
      !  Implemented the kinetic energy spectrum, dumped into an ascii file.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#define __OLDSCHOOL 1
#if __KIND == __DOUBLE_PRECISION
#error DOUBLE PRECISION NOT IMPLEMENTED [MB]
#endif
#if   __WHAT == __KESPECTRUM
#undef __PPM_FFT_GENERIC
#define __PPM_FFT_GENERIC ppm_fft_kespectrum
#elif __WHAT == __KESPECTRUM_KX
#undef __PPM_FFT_GENERIC
#define __PPM_FFT_GENERIC ppm_fft_kespectrumkx
#endif
      SUBROUTINE __PPM_FFT_GENERIC (field_rhs, mesh_id, topo_id, &
           &                        t_topo_id, t_mesh_id,        &
           &                        ghostsize, info,filename) 

        
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_mktopo
        USE ppm_module_map_field
        USE ppm_module_data
        USE ppm_module_data_mesh
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_mesh_define
        USE ppm_module_write
        USE MPI
        USE, INTRINSIC :: ISO_C_BINDING
        INCLUDE 'fftw3.f03'
        INTEGER, PARAMETER :: mk = ppm_kind_single
        !-----------------------------------------------------------------------
        !  arguments
        !-----------------------------------------------------------------------
        REAL(mk), DIMENSION(:,:,:,:,:), POINTER       :: field_rhs  ! rrhs
        INTEGER                       , INTENT(in   ) :: mesh_id    ! mesh id
        INTEGER                       , INTENT(in   ) :: topo_id    ! mesh id
        INTEGER,  DIMENSION(3)        , INTENT(in   ) :: ghostsize
        INTEGER                       , INTENT(  out) :: info
        CHARACTER(len=256)            , INTENT(in   ) :: filename
        !-----------------------------------------------------------------------
        
        !-----------------------------------------------------------------------
        !  localities
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(4), INTENT(in)             :: t_topo_id  ! tmp topo
        INTEGER, DIMENSION(4), INTENT(in)             :: t_mesh_id  ! tmp mesh
        INTEGER :: isub, isubl, i, j, k, iopt, istat
        INTEGER, DIMENSION(5)                         :: ldl, ldu
        INTEGER                                       :: topoid, meshid
        CHARACTER(len=256)                            :: msg, pnm
        REAL(mk)                                      :: t0
        COMPLEX(mk), DIMENSION(:,:,:,:,:) , POINTER   :: cfld1, cfld2
        INTEGER, DIMENSION(3)                         :: t_ghostsize
        INTEGER                                       :: ct1,ct2,ct3,ct4,ct5,ct6
        !-----------------------------------------------------------------------
        !  topology stuff
        !-----------------------------------------------------------------------
        INTEGER, DIMENSION(6)                         :: bcdef
        INTEGER                                       :: decomp, assignm
        REAL(mk), DIMENSION(:,:)      , POINTER       :: xp
        INTEGER                                       :: np
        REAL(mk), DIMENSION(3)                        :: min_phys
        REAL(mk), DIMENSION(3)                        :: max_phys
        INTEGER , DIMENSION(3)                        :: nx
        INTEGER , DIMENSION(:,:)      , POINTER       :: ndata, istart
        INTEGER , DIMENSION(:,:)      , POINTER       :: ndata_xy
        INTEGER , DIMENSION(:  )      , POINTER       :: isublist_xy
        INTEGER                                       :: nsubs_xy
        INTEGER , DIMENSION(:  )      , POINTER       :: isublist,sub2proc
        INTEGER                                       :: nsublist, nsubs
        REAL(mk), DIMENSION(:,:)      , POINTER       :: min_sub, max_sub
        REAL(mk), DIMENSION(:  )      , POINTER       :: sub_cost
        INTEGER                                       :: u_topo_id
        INTEGER                                       :: v_mesh_id, u_mesh_id
        INTEGER                                       :: maptype
        INTEGER                                       :: dbg_meshid, dbg_topoid
        INTEGER                                       :: ndata_old
        LOGICAL                                       :: verbose = .false.
        INTEGER                                       :: mpi_prec
        !-----------------------------------------------------
        !  fft stuff
        !-----------------------------------------------------
        INTEGER                                       :: idist, odist, rank
        INTEGER                                       :: istride, ostride
        INTEGER,  DIMENSION(2)                        :: fnx, iembed, oembed
        INTEGER                                       :: howmany
        TYPE(C_PTR)                                   :: plan
        COMPLEX(mk), DIMENSION(:,:), POINTER          :: c_i, c_o

        !-----------------------------------------------------
        !  poisson stuff
        !-----------------------------------------------------
        INTEGER , DIMENSION(2)                        :: lda_end
        REAL(mk)                                      :: pi2_lx,pi2_ly,pi2_lz
        REAL(mk)                                      :: kinv
#if   __WHAT == __KESPECTRUM
        COMPLEX(mk)                                   :: kx,ky,kz
        REAL(mk)                                      :: rkx,rky,rkz,ivol,isqrtvol
        REAL(mk)                                      :: sigx,sigy,sigz
        REAL(mk), DIMENSION(3)                        :: dx, dxi
        ! local and global arrays for kinetic energy spectrum arrays
        REAL(mk), DIMENSION(:), POINTER               :: kespecl, kespecg
        INTEGER                                       :: ike, ikemin, ikemax, nkespec
        REAL(mk)                                      :: ksqrt, &
             ksq, dirac, dirac_eps, dirac_epsi, kmax_kespec, &
             dkespec, dkespeci, ke, symfac
        INTEGER                                       :: kefile
#elif __WHAT == __KESPECTRUM_KX
        COMPLEX(mk)                                   :: kx,ky,kz
        COMPLEX(mk), DIMENSION(3)                     :: uhat
        REAL(mk)                                      :: rkx,rky,rkz,ivol,isqrtvol
        REAL(mk)                                      :: sigx,sigy,sigz
        REAL(mk), DIMENSION(3)                        :: dx, dxi
        ! local and global arrays for kinetic energy spectrum arrays
        REAL(mk), DIMENSION(:), POINTER               :: kespecl, kespecg
        INTEGER                                       :: ike, nkespec
        REAL(mk)                                      :: ksqrt, ksq, kmax_kespec, dkespec, dkespeci, ke, symfac,dkxdky
        INTEGER                                       :: kefile
#endif
        INTEGER                                       :: ig,jg,kg
        REAL(mk), DIMENSION(3)                        :: length
        
#if   __WHAT == __KESPECTRUM
        WRITE(pnm,*) 'ppm_fft_kespectrum'
#elif  __WHAT == __KESPECTRUM_KX
        WRITE(pnm,*) 'ppm_fft_kespectrumkx'
#endif

#if __KIND == __DOUBLE_PRECISION
        mpi_prec = MPI_DOUBLE_PRECISION
#else
        mpi_prec = MPI_REAL
#endif
        !-----------------------------------------------------------------------
        !  START
        !-----------------------------------------------------------------------
        CALL substart(pnm,t0,info)


        !-----------------------------------------------------------------------
        !  create topology     xy slab
        !-----------------------------------------------------------------------
        topoid  = ppm_field_topoid
        meshid  = ppm_meshid(topoid)%internal(mesh_id)
        min_phys(1:3) = ppm_min_physs(1:3,ppm_field_topoid)
        max_phys(1:3) = ppm_max_physs(1:3,ppm_field_topoid)
        length(1) = max_phys(3)-min_phys(3)
        length(2) = max_phys(2)-min_phys(2)
        length(3) = max_phys(1)-min_phys(1)
        nx(1:3) = ppm_cart_mesh(meshid,topoid)%nm(1:3)
        ivol = 1.0_mk/REAL((nx(1)-1)*(nx(2)-1)*(nx(3)-1))
        isqrtvol = SQRT(ivol)
#if   (__WHAT == __KESPECTRUM) ||  (__WHAT == __KESPECTRUM_KX)
        dx(1) = length(3)/REAL(nx(1)-1,mk)
        dx(2) = length(2)/REAL(nx(2)-1,mk)
        dx(3) = length(1)/REAL(nx(3)-1,mk)
        dxi(1) = 1.0_mk / dx(1)
        dxi(2) = 1.0_mk / dx(2)
        dxi(3) = 1.0_mk / dx(3)
#endif
        WRITE(msg,*) 'grid dimensions = ',nx
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        np      = -1
        t_ghostsize = 0
        assignm = ppm_param_assign_internal
        decomp  = ppm_param_decomp_xy_slab
        bcdef   = ppm_param_bcdef_periodic
        u_topo_id = t_topo_id(1)
        u_mesh_id = t_mesh_id(1)
        IF(Verbose) WRITE(msg,*) 'DBG--creating xy_slab topology'
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
                
        CALL ppm_mktopo (xp, np, nx, decomp, assignm, min_phys, max_phys, &
             &           bcdef, t_ghostsize, u_topo_id, u_mesh_id,    &
             &           min_sub, max_sub, sub_cost, sub2proc, nsubs,   &
             &           isublist, nsublist, istart, ndata, info)
        IF(Verbose) WRITE(msg,*) 'DBG--ndata = ',ndata
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        
        ALLOCATE(ndata_xy(3,4096))
        ALLOCATE(isublist_xy(4096))
        ndata_xy(1:3,isublist(1:nsublist)) = ndata(1:3,isublist(1:nsublist))
        isublist_xy(1:nsubs) = isublist(1:nsubs)
        nsubs_xy = nsubs
        WRITE(msg,*) ndata_xy(1:3,isublist(1:nsublist))
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        !-----------------------------------------------------------------------
        WRITE(msg,*) 'xy_slab created'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)

        !-----------------------------------------------------------------------
        !  map data from user topo to xy_slab topology
        !  hopefully this baby is gonna get its ghosts busted by the mapping
        !-----------------------------------------------------------------------
        maptype = ppm_param_map_global
        WRITE(msg,*) 'before map ',UBOUND(field_rhs)
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        WRITE(msg,*) 'before map topo/mesh/mesh',t_topo_id,mesh_id,t_mesh_id(1)
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        CALL ppm_map_field(field_rhs, 3, t_topo_id(1), mesh_id, t_mesh_id(1),&
             &             t_ghostsize, maptype, info)
        maptype = ppm_param_map_push
        CALL ppm_map_field(field_rhs, 3, t_topo_id(1), mesh_id, t_mesh_id(1),&
             &             t_ghostsize, maptype, info)
        maptype = ppm_param_map_send
        CALL ppm_map_field(field_rhs, 3, t_topo_id(1), mesh_id, t_mesh_id(1),&
             &             t_ghostsize, maptype, info)
        maptype = ppm_param_map_pop
        CALL ppm_map_field(field_rhs, 3, t_topo_id(1), mesh_id, t_mesh_id(1),&
             &             t_ghostsize, maptype, info)
        WRITE(msg,*) 'after map ',UBOUND(field_rhs)
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  get new internal ids
        !-----------------------------------------------------------------------
        WRITE(msg,*) 'retrieving new topo specs'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        topoid = ppm_field_topoid
        meshid = ppm_meshid(topoid)%internal(t_mesh_id(1))
        
        WRITE(msg,*) 'allocating complex representation'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        
        !-----------------------------------------------------------------------
        !  allocate complex array (TODO:check if size is ok)
        !-----------------------------------------------------------------------
        ndata_old = MAXVAL(ndata(1,isublist(1:nsublist)))
        ldu(1) = 3
        ldu(2) = MAXVAL(ndata(1,isublist(1:nsublist))-1)/2+2
        ldu(3) = MAXVAL(ndata(2,isublist(1:nsublist)))
        ldu(4) = MAXVAL(ndata(3,isublist(1:nsublist)))
        ldu(5) = nsublist
        iopt   = ppm_param_alloc_fit
        CALL ppm_alloc(cfld1, ldu, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,pnm,'crhs',__LINE__,info)
           GOTO 9999
        END IF
        WRITE(msg,*) 'sizeof(cfld1) = ',UBOUND(cfld1)
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        WRITE(msg,*) 'sizeof(flrhs) = ',UBOUND(field_rhs)
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        !-----------------------------------------------------------------------

        
        !-----------------------------------------------------------------------
        !  2D FFT - real->cfld1
        !  USING FFTW:
        !   a 3-fold 2d fft 
        !-----------------------------------------------------------------------
        DO isub=1,nsublist
           isubl = isublist(isub)
           WRITE(msg,*) 'XYFFT sub[',isub,'/',nsublist,']'
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
           DO k=1,ndata(3,isublist(isub))
              rank      = 2
              fnx(1)    = ndata(1,isublist(isub))-1
              fnx(2)    = ndata(2,isublist(isub))-1
              iembed(1) = UBOUND(field_rhs,2)
              iembed(2) = UBOUND(field_rhs,3)
              oembed(1) = UBOUND(cfld1,2)
              oembed(2) = UBOUND(cfld1,3)
              istride   = 3
              ostride   = 3
              idist     = 1
              odist     = 1
              howmany   = 3
              CALL sfftw_plan_many_dft_r2c(plan, rank, fnx(1), howmany, &
                   & field_rhs(1,1,1,k,isub), iembed(1), istride, idist, &
                   & cfld1    (1,1,1,k,isub), oembed(1), ostride, odist, &
                   & FFTW_ESTIMATE)
              CALL sfftw_execute(plan)
              CALL sfftw_destroy_plan(plan)
              DO j=1,ndata(2,isubl)
                 cfld1(1,(ndata(1,isubl)-1)/2+2,j,k,isub) = cfld1(1,1,j,k,isub)
                 cfld1(2,(ndata(1,isubl)-1)/2+2,j,k,isub) = cfld1(2,1,j,k,isub)
                 cfld1(3,(ndata(1,isubl)-1)/2+2,j,k,isub) = cfld1(3,1,j,k,isub)
              END DO
              DO i=1,(ndata(1,isubl)-1)/2+2
                 cfld1(1,i,ndata(2,isubl),k,isub) = cfld1(1,i,1,k,isub)
                 cfld1(2,i,ndata(2,isubl),k,isub) = cfld1(2,i,1,k,isub)
                 cfld1(3,i,ndata(2,isubl),k,isub) = cfld1(3,i,1,k,isub)
              END DO
              cfld1(1,(ndata(1,isubl)-1)/2+2,ndata(2,isubl),k,isub) = &
                   & cfld1(1,1,1,k,isub)
              cfld1(2,(ndata(1,isubl)-1)/2+2,ndata(2,isubl),k,isub) = &
                   & cfld1(2,1,1,k,isub)
              cfld1(3,(ndata(1,isubl)-1)/2+2,ndata(2,isubl),k,isub) = &
                   & cfld1(3,1,1,k,isub)
              

           END DO
           !-----------------------------------------------------
           ! copy margin for ghosts
           !-----------------------------------------------------

        END DO

        !-----------------------------------------------------------------------
        
        iopt = ppm_param_dealloc
        CALL ppm_alloc(field_rhs,ldu,iopt,info)

        
        !-----------------------------------------------------------------------
        !  make a mesh conforming with the complex field
        !-----------------------------------------------------------------------
        nx(1) = (nx(1)-1)/2 + 1
        u_mesh_id = t_mesh_id(2)
        u_topo_id = t_topo_id(1)
        CALL ppm_mesh_define(u_topo_id,nx,min_phys,max_phys,              &
             &               u_mesh_id, istart, ndata, info)
        WRITE(msg,*) 'created complex mesh: nx = ',nx
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        !-----------------------------------------------------------------------

        
        !-----------------------------------------------------------------------
        !  decompose domain into z pencils
        !-----------------------------------------------------------------------
        decomp = ppm_param_decomp_zpencil
        u_topo_id = t_topo_id(2)
        u_mesh_id = t_mesh_id(3)
        CALL ppm_mktopo(xp,np,nx,decomp,assignm,min_phys,max_phys,           &
             &          bcdef, t_ghostsize, u_topo_id, u_mesh_id,            &
             &          min_sub, max_sub, sub_cost, sub2proc, nsubs,             &
             &          isublist, nsublist, istart, ndata, info)
        !-----------------------------------------------------------------------

        IF(ppm_internal_topoid(t_topo_id(1)).NE.ppm_field_topoid) THEN
           WRITE(msg,*) '"from" topology is not current topology (?)'
           IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
           STOP
        END IF
        !-----------------------------------------------------------------------
        !  map data from xy slab to z pencil
        !-----------------------------------------------------------------------
        dbg_topoid = ppm_internal_topoid(t_topo_id(1))
        dbg_meshid = ppm_meshid(dbg_topoid)%internal(t_mesh_id(2))
        WRITE(msg,*) '"from" mesh specs: ',             &
             & ppm_cart_mesh(dbg_meshid,dbg_topoid)%nm(1:3)
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        dbg_topoid = ppm_internal_topoid(t_topo_id(2))
        dbg_meshid = ppm_meshid(dbg_topoid)%internal(t_mesh_id(3))
        WRITE(msg,*) '  "to" mesh specs: ',             &
             & ppm_cart_mesh(dbg_meshid,dbg_topoid)%nm(1:3)
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        
        
        maptype = ppm_param_map_global
        CALL ppm_map_field(cfld1, 3, t_topo_id(2),t_mesh_id(2),t_mesh_id(3),&
             &             t_ghostsize, maptype, info)
        maptype = ppm_param_map_push
        CALL ppm_map_field(cfld1, 3, t_topo_id(2),t_mesh_id(2),t_mesh_id(3),&
             &             t_ghostsize, maptype, info)
        maptype = ppm_param_map_send
        CALL ppm_map_field(cfld1, 3, t_topo_id(2),t_mesh_id(2),t_mesh_id(3),&
             &             t_ghostsize, maptype, info)
        maptype = ppm_param_map_pop
        CALL ppm_map_field(cfld1, 3, t_topo_id(2),t_mesh_id(2),t_mesh_id(3),&
             &             t_ghostsize, maptype, info)

        !-----------------------------------------------------------------------
        WRITE(msg,*) 'mapped to z-pencil topology'
        IF(verbose) CALL  ppm_write(ppm_rank,pnm,msg,info)
        
        !-----------------------------------------------------------------------
        !  now transpose data
        !-----------------------------------------------------------------------
        ldu(1) = 3
        ldu(2) = MAXVAL(ndata(3,isublist(1:nsublist)))
        ldu(3) = MAXVAL(ndata(2,isublist(1:nsublist)))
        ldu(4) = MAXVAL(ndata(1,isublist(1:nsublist)))
        ldu(5) = nsublist
        iopt   = ppm_param_alloc_fit
        CALL ppm_alloc(cfld2,ldu,iopt,info)

        !-----------------------------------------------------------------------
        ! transpose!
        !-----------------------------------------------------------------------
        DO isub=1,nsublist
           isubl = isublist(isub)
           DO k=1,ndata(1,isubl)
              DO j=1,ndata(2,isubl)
                 DO i=1,ndata(3,isubl)
                    cfld2(1,i,j,k,isub) = cfld1(1,k,j,i,isub)
                    cfld2(2,i,j,k,isub) = cfld1(2,k,j,i,isub)
                    cfld2(3,i,j,k,isub) = cfld1(3,k,j,i,isub)
                 END DO
              END DO
           END DO
        END DO
        
        !-----------------------------------------------------------------------
        ! reshape cfld1
        !-----------------------------------------------------------------------
        iopt = ppm_param_dealloc
        CALL ppm_alloc(cfld1,ldu,iopt,info)
        iopt = ppm_param_alloc_fit
        CALL ppm_alloc(cfld1,ldu,iopt,info)
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  now transform into the z direction, cfld2->cfld1
        !  USING FFTW:
        !   a nx*ny-fold 1d fft 
        !-----------------------------------------------------------------------
#ifdef __OLDSCHOOL
        isub = 1
        iopt = ppm_param_alloc_fit
        ldl(1) = 3
        ldl(2) = ndata(3,(isublist(isub)))-1
        CALL ppm_alloc(c_i,ldl,iopt,info)
        CALL ppm_alloc(c_o,ldl,iopt,info)
        DO isub=1,nsublist
           isubl = isublist(isub)
           WRITE(msg,*) ' ZFFT sub[',isub,'/',nsublist,']'
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
           rank      = 1
           fnx(1)    = ndata(3,isublist(isub))-1
           iembed(1) = -1!UBOUND(cfld2,2)*3
           oembed(1) = -1!UBOUND(cfld1,2)*3
           idist     = 1
           odist     = 1
           howmany   = 3
           CALL sfftw_plan_many_dft(plan, rank, fnx(1), howmany, &
                   & c_i, iembed(1), istride, idist, &
                   & c_o, iembed(1), ostride, odist, &
                   & FFTW_FORWARD, FFTW_ESTIMATE)
           DO k=1,ndata(1,isublist(isub))
              DO j=1,ndata(2,isublist(isub))
                 DO i=1,ndata(3,isublist(isub))-1
                    c_i(1,i) = cfld2(1,i,j,k,isub)
                    c_i(2,i) = cfld2(2,i,j,k,isub)
                    c_i(3,i) = cfld2(3,i,j,k,isub)
                 END DO
                 CALL sfftw_execute(plan)
                 DO i=1,ndata(3,isublist(isub))-1
                    cfld1(1,i,j,k,isub) = c_o(1,i)
                    cfld1(2,i,j,k,isub) = c_o(2,i)
                    cfld1(3,i,j,k,isub) = c_o(3,i)
                 END DO
                 cfld1(1,ndata(3,isublist(isub)),j,k,isub) = c_o(1,1)
                 cfld1(2,ndata(3,isublist(isub)),j,k,isub) = c_o(2,1)
                 cfld1(3,ndata(3,isublist(isub)),j,k,isub) = c_o(3,1)
                 
              END DO
           END DO
           CALL sfftw_destroy_plan(plan)

        END DO
#else
        iopt = ppm_param_alloc_fit
        ldl(1) = ndata(3,(isublist(isub)))-1
        ldl(2) = 3
        CALL ppm_alloc(c_i,ldl,iopt,info)
        CALL ppm_alloc(c_o,ldl,iopt,info)
        DO isub=1,nsublist
           WRITE(msg,*) ' ZFFT sub[',isub,'/',nsublist,']'
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
           rank      = 1
           fnx(1)    = ndata(3,isublist(isub))-1
           iembed(1) = -1!UBOUND(cfld2,2)*3
           oembed(1) = -1!UBOUND(cfld1,2)*3
           idist     = fnx(1)
           odist     = fnx(1)
           istride   = 1
           ostride   = 1
           howmany   = 3
           CALL sfftw_plan_many_dft(plan, rank, fnx(1), howmany, &
                   & c_i, iembed(1), istride, idist, &
                   & c_o, iembed(1), ostride, odist, &
                   & FFTW_FORWARD, FFTW_ESTIMATE)
           DO k=1,ndata(1,isublist(isub))
              DO j=1,ndata(2,isublist(isub))
                 DO i=1,ndata(3,isublist(isub))-1
                    c_i(i,1) = cfld2(1,i,j,k,isub)
                    c_i(i,2) = cfld2(2,i,j,k,isub)
                    c_i(i,3) = cfld2(3,i,j,k,isub)
                 END DO
                 CALL sfftw_execute(plan)
                 DO i=1,ndata(3,isublist(isub))-1
                    cfld1(1,i,j,k,isub) = c_o(i,1)
                    cfld1(2,i,j,k,isub) = c_o(i,2)
                    cfld1(3,i,j,k,isub) = c_o(i,3)
                 END DO
              END DO
           END DO
           CALL sfftw_destroy_plan(plan)

        END DO
#endif
        !-----------------------------------------------------------------------
        ! DO OPERATIONS HERE: (cfld1->cfld1)
        !                      KESPECTRUM
        !                      VELOCITIES
        !                      WHATNOT
        ! check istart!
        !-----------------------------------------------------------------------
#if __WHAT == __KESPECTRUM
        pi2_Lz = 2.0_MK*ppm_pi_s/length(1)
        pi2_Ly = 2.0_MK*ppm_pi_s/length(2)
        pi2_Lx = 2.0_MK*ppm_pi_s/length(3)
        lda_end(1)=(Nx(3)-1)/2 +1
        lda_end(2)=(Nx(2)-1)/2 +1
                
        ! Allocate local and global ke spectrum arrays
        nkespec = MIN(lda_end(1),lda_end(2),nx(1))
        kmax_kespec = ppm_pi_s * MINVAL(dxi)
        dkespec = kmax_kespec/REAL(nkespec-1,mk)
        dkespeci = 1.0_mk/dkespec
        ldl(1)  = 0
        ldu(1)  = nkespec-1
        iopt = ppm_param_alloc_fit
        CALL ppm_alloc(kespecl,ldl,ldu,iopt,info)
        CALL ppm_alloc(kespecg,ldl,ldu,iopt,info)
        kespecl = 0.0_mk
        kespecg = 0.0_mk

        DO isub=1,nsublist
           isubl = isublist(isub)
           DO k=1,ndata(1,isublist(isub))-1
              kg = k + istart(1,isubl)-2
              sigz = 1.0_mk
              symfac = 1.0_mk
              IF (kg.EQ.0) symfac = 1.0_mk
              DO j=1,ndata(2,isublist(isub))-1
                 jg = j + istart(2,isubl)-2
                 sigy = 1.0_mk
                 IF(jg.GE.lda_end(2)) THEN
                    jg = nx(2)-jg-1
                    sigy = -1.0_mk
                 END IF
                 DO i=1,ndata(3,isublist(isub))-1
                    ig = i + istart(3,isubl)-2
                    sigx = 1.0_mk
                    IF(ig.GE.(lda_end(1))) THEN
                        ig = nx(3)-ig-1
                        sigx = -1.0_mk
                    END IF
                    IF(ig.EQ.0.AND.jg.EQ.0.AND.kg.EQ.0) THEN
                       ke = 0.0_mk
                    ELSE
                       rkz=pi2_Lx*REAL(ig,MK)*sigx
                       rky=pi2_Ly*REAL(jg,MK)*sigy
                       rkx=pi2_Lz*REAL(kg,MK)
                       ! Contribution of this point to energy
                       ke =  REAL(cfld1(1,i,j,k,isub)*CONJG(cfld1(1,i,j,k,isub)) + &
                       &          cfld1(2,i,j,k,isub)*CONJG(cfld1(2,i,j,k,isub)) + &
                       &          cfld1(3,i,j,k,isub)*CONJG(cfld1(3,i,j,k,isub)),mk)
                       
                       ksq = (rkx*rkx + rky*rky + rkz*rkz)
                       ksqrt = SQRT(ksq)
                       !kinv = (1.0_mk/ksq)*ivol
                       !ke =  REAL(cfld1(1,i,j,k,isub)*CONJG(cfld1(1,i,j,k,isub)) + &
                       !&          cfld1(2,i,j,k,isub)*CONJG(cfld1(2,i,j,k,isub)) + &
                       !&          cfld1(3,i,j,k,isub)*CONJG(cfld1(3,i,j,k,isub)),mk) * kinv
                       ! Computation of index in kespec array
                       ike = INT(ksqrt*dkespeci)
                       ! Smoothed dirac epsilon (Engquist, 2005) for consistent Dirac
                       ! integration
                       dirac_eps = (ABS(rkx)+ABS(rky)+ABS(rkz))/ksqrt * dkespec
                       dirac_epsi = 1.0_mk/dirac_eps
                       ! Starting and ending points for the contributions
                       ikemin = MAX(INT((ksqrt-dirac_eps)*dkespeci)+1,1)
                       ikemax = MIN(INT((ksqrt+dirac_eps)*dkespeci),nkespec-1)
                       !WRITE(*,*) 'k ', rkx, rky, rkz, ' ksqrt ', ksqrt, ke 
                       DO ike=ikemin,ikemax
                          ! MAX it with 0, just to remove the (tiny) negative contributions
                          dirac = MAX(0.0_mk,(1.0 - ABS(ksqrt-ike*dkespec)*dirac_epsi)*dirac_epsi)
                          !WRITE(*,*) ' dirac ', dirac, ' modified eps', dirac_eps 
                          ! CAVEATs
                          ! 1) Shouldnt there be a dk volume factor somewhere?
                          ! 2) Using the fact that this is the transform of a real function
                          !    we have been only considering the positive freqs of x
                          !    cfld1(-k) = CONJG(clfld1(k)) -> contributions are the same
                          !    There is a factor of 2.0_mk for all the points corresponding to
                          !    kg.GT.0 and 1.0_mk for the frequencies in the symmetry plane
                          ! 3) Why transform back cfld1??? Dont know... Should be removed probably
                          !   Unless you want to get a correlation function in physical space...
                          kespecl(ike) = kespecl(ike) + symfac*dirac*ke
                       END DO
                    ENDIF
                 END DO
              END DO
           END DO
        END DO
        
        ! We re done we reduce all these guys into the global array
        CALL MPI_Reduce(kespecl,kespecg,nkespec,mpi_prec,MPI_SUM,0,ppm_comm,info)
        IF (ppm_rank.EQ.0) THEN
           ! Factor for Parseval identity
           kespecg = 0.5*ppm_pi_s * (0.5_mk/ppm_pi_s)**3 * (dx(1)*dx(2)*dx(3))**2 * pi2_Lx*pi2_Ly*pi2_Lz * kespecg
           
           kefile = 20
           OPEN(kefile,file=filename,iostat=info,action='WRITE')
           WRITE(kefile,*) (/ (i*dkespec,i=0,nkespec-1) /)
           WRITE(kefile,*) kespecg(0:nkespec-1)
           CLOSE(kefile,iostat=info)
        END IF
        IF (ppm_rank.EQ.0) THEN
           WRITE(msg,'(A,F10.6)') 'Total energy (summation over shells) ', SUM(kespecg)*dkespec
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        END IF
        iopt = ppm_param_dealloc
        CALL ppm_alloc(kespecl,ldl,ldu,iopt,info)
        CALL ppm_alloc(kespecg,ldl,ldu,iopt,info)
        
#elif __WHAT == __KESPECTRUM_KX

        pi2_Lz = 2.0_MK*ppm_pi_s/length(1)
        pi2_Ly = 2.0_MK*ppm_pi_s/length(2)
        pi2_Lx = 2.0_MK*ppm_pi_s/length(3)
        lda_end(1)=(Nx(3)-1)/2 +1
        lda_end(2)=(Nx(2)-1)/2 +1
                
        ! Allocate local and global ke spectrum arrays
        nkespec = lda_end(1)
        kmax_kespec = ppm_pi_s * MINVAL(dxi)
        dkespec = kmax_kespec/REAL(nkespec-1,mk)
        dkespeci = 1.0_mk/dkespec
        dkxdky = pi2_Lx * pi2_Ly
        
        ldl(1)  = 0
        ldu(1)  = nkespec-1
        iopt = ppm_param_alloc_fit
        CALL ppm_alloc(kespecl,ldl,ldu,iopt,info)
        CALL ppm_alloc(kespecg,ldl,ldu,iopt,info)
        kespecl = 0.0_mk
        kespecg = 0.0_mk

        DO isub=1,nsublist
           isubl = isublist(isub)
           DO k=1,ndata(1,isublist(isub))-1
              kg = k + istart(1,isubl)-2
              sigz = 1.0_mk
              symfac = 1.0_mk
              IF (kg.EQ.0) symfac = 1.0_mk
              DO j=1,ndata(2,isublist(isub))-1
                 jg = j + istart(2,isubl)-2
                 sigy = 1.0_mk
                 IF(jg.GE.lda_end(2)) THEN
                    jg = nx(2)-jg-1
                    sigy = -1.0_mk
                 END IF
                 DO i=1,ndata(3,isublist(isub))-1
                    ig = i + istart(3,isubl)-2
                    sigx = 1.0_mk
                    IF(ig.GE.(lda_end(1))) THEN
                        ig = nx(3)-ig-1
                        sigx = -1.0_mk
                    END IF
                    IF(ig.EQ.0.AND.jg.EQ.0.AND.kg.EQ.0) THEN
                       ke = 0.0_mk
                    ELSE
                       rkz=pi2_Lx*REAL(ig,MK)*sigx
                       rky=pi2_Ly*REAL(jg,MK)*sigy
                       rkx=pi2_Lz*REAL(kg,MK)
                       
                       ke =  REAL(cfld1(1,i,j,k,isub)*CONJG(cfld1(1,i,j,k,isub)) + &
                       &          cfld1(2,i,j,k,isub)*CONJG(cfld1(2,i,j,k,isub)) + &
                       &          cfld1(3,i,j,k,isub)*CONJG(cfld1(3,i,j,k,isub)),mk)
                    
                       
                       ! Contribution of this point to energy
                       !ksq = (rkx*rkx + rky*rky + rkz*rkz)
                       !ksqrt = SQRT(ksq)
                       !kinv = (1.0_mk/ksq)*ivol
                       
                       !kx=CMPLX(0.0_mk,rkx)
                       !ky=CMPLX(0.0_mk,rky)
                       !kz=CMPLX(0.0_mk,rkz)
                       !-----------------------------------------------------
                       ! take the curl of the inverse nabla
                       !kinv = 1.0_mk/(rkx*rkx + rky*rky + rkz*rkz)*isqrtvol
                       !uhat(1)=kinv*(ky*cfld1(3,i,j,k,isub)-kz*cfld1(2,i,j,k,isub))
                       !uhat(2)=kinv*(kz*cfld1(1,i,j,k,isub)-kx*cfld1(3,i,j,k,isub))
                       !uhat(3)=kinv*(kx*cfld1(2,i,j,k,isub)-ky*cfld1(1,i,j,k,isub))
                       !ke =  REAL(uhat(1)*CONJG(uhat(1)) + &
                       !           uhat(2)*CONJG(uhat(2)) + &
                       !           uhat(3)*CONJG(uhat(3)),mk)
                       
                       !ke =  REAL(cfld1(1,i,j,k,isub)*CONJG(cfld1(1,i,j,k,isub)) + &
                       !&          cfld1(2,i,j,k,isub)*CONJG(cfld1(2,i,j,k,isub)) + &
                       !&          cfld1(3,i,j,k,isub)*CONJG(cfld1(3,i,j,k,isub)),mk) * kinv
                       ! Computation of index in kespec array
                       ike = ig
                       ! CAVEATs
                       ! 1) Shouldnt there be a dk volume factor somewhere?
                       ! 2) Using the fact that this is the transform of a real function
                       !    we have been only considering the positive freqs of x
                       !    cfld1(-k) = CONJG(clfld1(k)) -> contributions are the same
                       !    There is a factor of 2.0_mk for all the points corresponding to
                       !    kg.GT.0 and 1.0_mk for the frequencies in the symmetry plane
                       ! 3) Why transform back cfld1??? Dont know... Should be removed probably
                       !   Unless you want to get a correlation function in physical space...
                       !kespecl(ike) = kespecl(ike) + dkxdky*symfac*ke
                       kespecl(ike) = kespecl(ike) + symfac*ke
                    ENDIF
                 END DO
              END DO
           END DO
        END DO
        
        ! We re done we reduce all these guys into the global array
        CALL MPI_Reduce(kespecl,kespecg,nkespec,mpi_prec,MPI_SUM,0,ppm_comm,info)
        IF (ppm_rank.EQ.0) THEN
           kefile = 20
           ! Factor for Parseval identity
           kespecg = 0.5*ppm_pi_s * (0.5_mk/ppm_pi_s)**3 * (dx(1)*dx(2)*dx(3))**2 * dkxdky * kespecg
           OPEN(kefile,file=filename,iostat=info,action='WRITE')
           WRITE(kefile,*) (/ (i*dkespec,i=0,nkespec-1) /)
           WRITE(kefile,*) kespecg(0:nkespec-1)
           CLOSE(kefile,iostat=info)
        END IF
        IF (ppm_rank.EQ.0) THEN
           WRITE(msg,'(A,F10.8)') 'Total energy (summation over kx planes) ', SUM(kespecg)*dkespec
           IF (verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        END IF
        iopt = ppm_param_dealloc
        CALL ppm_alloc(kespecl,ldl,ldu,iopt,info)
        CALL ppm_alloc(kespecg,ldl,ldu,iopt,info)
#endif
        
      




        !-----------------------------------------------------
        !
        !
        !
        !
        !-----------------------------------------------------




        
        !-----------------------------------------------------------------------
        ! transform back in the z direction
        ! cfld1 -> cfld2
        !-----------------------------------------------------------------------
#ifdef __OLDSCHOOL
        DO isub=1,nsublist
           WRITE(msg,*) ' ZFFT sub[',isub,'/',nsublist,']'
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
           rank      = 1
           fnx(1)    = ndata(3,isublist(isub))-1
           iembed(1) = -1!UBOUND(cfld2,2)*3
           oembed(1) = -1!UBOUND(cfld1,2)*3
           idist     = 1
           odist     = 1
           howmany   = 3
           CALL sfftw_plan_many_dft(plan, rank, fnx(1), howmany, &
                   & c_i, iembed(1), istride, idist, &
                   & c_o, iembed(1), ostride, odist, &
                   & FFTW_BACKWARD, FFTW_ESTIMATE)
           DO k=1,ndata(1,isublist(isub))
              DO j=1,ndata(2,isublist(isub))
                 DO i=1,ndata(3,isublist(isub))-1
                    c_i(1,i) = cfld1(1,i,j,k,isub)
                    c_i(2,i) = cfld1(2,i,j,k,isub)
                    c_i(3,i) = cfld1(3,i,j,k,isub)
                 END DO
                 CALL sfftw_execute(plan)
                 DO i=1,ndata(3,isublist(isub))-1
                    cfld2(1,i,j,k,isub) = c_o(1,i)
                    cfld2(2,i,j,k,isub) = c_o(2,i)
                    cfld2(3,i,j,k,isub) = c_o(3,i)
                 END DO
                 cfld2(1,ndata(3,isublist(isub)),j,k,isub) = c_o(1,1)
                 cfld2(2,ndata(3,isublist(isub)),j,k,isub) = c_o(2,1)
                 cfld2(3,ndata(3,isublist(isub)),j,k,isub) = c_o(3,1)
              END DO
           END DO
           CALL sfftw_destroy_plan(plan)

        END DO
#else
        DO isub=1,nsublist
           WRITE(msg,*) ' ZFFT sub[',isub,'/',nsublist,']'
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
           rank      = 1
           fnx(1)    = ndata(3,isublist(isub))-1
           iembed(1) = -1!UBOUND(cfld2,2)*3
           oembed(1) = -1!UBOUND(cfld1,2)*3
           idist     = ndata(3,isublist(isub))-1
           odist     = ndata(3,isublist(isub))-1
           howmany   = 3
           istride   = 1
           ostride   = 1
           CALL sfftw_plan_many_dft(plan, rank, fnx(1), howmany, &
                   & c_i, iembed(1), istride, idist, &
                   & c_o, iembed(1), ostride, odist, &
                   & FFTW_BACKWARD, FFTW_ESTIMATE)
           DO k=1,ndata(1,isublist(isub))
              DO j=1,ndata(2,isublist(isub))
                 DO i=1,ndata(3,isublist(isub))-1
                    c_i(i,1) = cfld1(1,i,j,k,isub)
                    c_i(i,2) = cfld1(2,i,j,k,isub)
                    c_i(i,3) = cfld1(3,i,j,k,isub)
                 END DO
                 CALL sfftw_execute(plan)
                 DO i=1,ndata(3,isublist(isub))-1
                    cfld2(1,i,j,k,isub) = c_o(i,1)
                    cfld2(2,i,j,k,isub) = c_o(i,2)
                    cfld2(3,i,j,k,isub) = c_o(i,3)
                 END DO
              END DO
           END DO
           CALL sfftw_destroy_plan(plan)

        END DO

#endif
        
        !-----------------------------------------------------------------------
        ! reshape cfld1
        !-----------------------------------------------------------------------
        iopt =ppm_param_dealloc
        CALL ppm_alloc(c_i,ldu,iopt,info)
        CALL ppm_alloc(c_o,ldu,iopt,info)
        CALL ppm_alloc(cfld1,ldu,iopt,info)
        ldu(1) = 3
        ldu(2) = MAXVAL(ndata(1,isublist(1:nsublist)))
        ldu(3) = MAXVAL(ndata(2,isublist(1:nsublist)))
        ldu(4) = MAXVAL(ndata(3,isublist(1:nsublist)))
        ldu(5) = nsublist
        iopt   = ppm_param_alloc_fit
        CALL ppm_alloc(cfld1,ldu,iopt,info)


        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        ! transpose back
        !-----------------------------------------------------------------------
        DO isub=1,nsublist
           isubl = isublist(isub)
           DO k=1,ndata(3,isubl)
              DO j=1,ndata(2,isubl)
                 DO i=1,ndata(1,isubl)
                    cfld1(1,i,j,k,isub) = cfld2(1,k,j,i,isub)
                    cfld1(2,i,j,k,isub) = cfld2(2,k,j,i,isub)
                    cfld1(3,i,j,k,isub) = cfld2(3,k,j,i,isub)
                 END DO
              END DO
           END DO
        END DO
        !-----------------------------------------------------------------------

        
        !-----------------------------------------------------------------------
        ! map the stuff from the z to the xy slab topology
        !-----------------------------------------------------------------------
        IF(ppm_internal_topoid(t_topo_id(2)).NE.ppm_field_topoid) THEN
           WRITE(msg,*) '"from" topology is not current topology (?)'
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        END IF

        dbg_topoid = ppm_field_topoid
        dbg_meshid = ppm_meshid(dbg_topoid)%internal(t_mesh_id(3))
        WRITE(msg,*) '"from" mesh specs: ',             &
             & ppm_cart_mesh(dbg_meshid,dbg_topoid)%nm(1:3)
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        dbg_topoid = ppm_internal_topoid(t_topo_id(1))
        dbg_meshid = ppm_meshid(dbg_topoid)%internal(t_mesh_id(2))
        WRITE(msg,*) '  "to" mesh specs: ',             &
             & ppm_cart_mesh(dbg_meshid,dbg_topoid)%nm(1:3)
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)

        maptype = ppm_param_map_global
        CALL ppm_map_field(cfld1, 3, t_topo_id(1),t_mesh_id(3),t_mesh_id(2),&
             &             t_ghostsize, maptype, info)
        if(verbose) WRITE(msg,*) 'ppm_param_map_global'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        maptype = ppm_param_map_push
        CALL ppm_map_field(cfld1, 3, t_topo_id(1),t_mesh_id(3),t_mesh_id(2),&
             &             t_ghostsize, maptype, info)
        if(verbose) WRITE(msg,*) 'ppm_param_map_push'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        maptype = ppm_param_map_send
        CALL ppm_map_field(cfld1, 3, t_topo_id(1),t_mesh_id(3),t_mesh_id(2),&
             &             t_ghostsize, maptype, info)
        if(verbose) WRITE(msg,*) 'ppm_param_map_send'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        maptype = ppm_param_map_pop
        CALL ppm_map_field(cfld1, 3, t_topo_id(1),t_mesh_id(3),t_mesh_id(2),&
             &             t_ghostsize, maptype, info)
        if(verbose) WRITE(msg,*) 'ppm_param_map_pop'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        !-----------------------------------------------------------------------
        WRITE(msg,*) 'mapped back to xy-pencil'
        if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        
        !-----------------------------------------------------------------------
        ! dealloc cfld2
        !-----------------------------------------------------------------------
        iopt = ppm_param_dealloc
        CALL ppm_alloc(cfld2,ldu,iopt,info)
        !-----------------------------------------------------------------------
        ldu(2) = MAXVAL(ndata_xy(1,isublist_xy(1:nsublist)))
        ldu(3) = MAXVAL(ndata_xy(2,isublist_xy(1:nsublist)))
        ldu(4) = MAXVAL(ndata_xy(3,isublist_xy(1:nsublist)))
        ldu(5) = nsublist
        iopt = ppm_param_alloc_fit
        CALL ppm_alloc(field_rhs,ldu,iopt,info)
        WRITE(msg,*) 'recreated field_rhs w/size ',UBOUND(field_rhs)
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        isublist(1:nsubs_xy) = isublist_xy(1:nsubs_xy)
        ndata(1:3,isublist(1:nsublist)) = ndata_xy(1:3,isublist(1:nsublist))
        !-----------------------------------------------------------------------
        ! 2d ifft in xy
        !  USING FFTW:
        !   a 3-fold 2d fft 
        !-----------------------------------------------------------------------
        DO isub=1,nsublist
           isubl = isublist(isub)
           WRITE(msg,*) 'iXYFFT sub[',isub,'/',nsublist,']'
           if(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
           DO k=1,ndata(3,isubl)
              rank      = 2
              ! have to save ndata_old for all isubs!!!!!! TODO
              fnx(1)    = ndata(1,isubl)-1
              fnx(2)    = ndata(2,isubl)-1
              oembed(1) = UBOUND(field_rhs,2)
              oembed(2) = UBOUND(field_rhs,3)
              iembed(1) = UBOUND(cfld1,2)
              iembed(2) = UBOUND(cfld1,3)
              idist     = 1
              odist     = 1
              istride   = 3
              ostride   = 3
              howmany   = 3
              CALL sfftw_plan_many_dft_c2r(plan, rank, fnx(1), howmany, &
                   & cfld1    (1,1,1,k,isub), iembed(1), istride, idist, &
                   & field_rhs(1,1,1,k,isub), oembed(1), ostride, odist, &
                   & FFTW_ESTIMATE)
              CALL sfftw_execute(plan)
              CALL sfftw_destroy_plan(plan)
              DO j=1,ndata(2,isubl)
                 field_rhs(1,ndata(1,isubl),j,k,isub) = field_rhs(1,1,j,k,isub)
                 field_rhs(2,ndata(1,isubl),j,k,isub) = field_rhs(2,1,j,k,isub)
                 field_rhs(3,ndata(1,isubl),j,k,isub) = field_rhs(3,1,j,k,isub)
              END DO
              DO i=1,ndata(1,isubl)
                 field_rhs(1,i,ndata(2,isubl),k,isub) = field_rhs(1,i,1,k,isub)
                 field_rhs(2,i,ndata(2,isubl),k,isub) = field_rhs(2,i,1,k,isub)
                 field_rhs(3,i,ndata(2,isubl),k,isub) = field_rhs(3,i,1,k,isub)
              END DO
              field_rhs(1,ndata(1,isubl),ndata(2,isubl),k,isub) = &
                   & field_rhs(1,1,1,k,isub)
              field_rhs(2,ndata(1,isubl),ndata(2,isubl),k,isub) = &
                   & field_rhs(2,1,1,k,isub)
              field_rhs(3,ndata(1,isubl),ndata(2,isubl),k,isub) = &
                   & field_rhs(3,1,1,k,isub)
              
              
           END DO
        END DO
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        ! map back to user topology
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        ! put ghostlayer back on
        !-----------------------------------------------------------------------
        WRITE(msg,*) 'before map field_rhs is ',UBOUND(field_rhs)
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)
        iopt = ppm_param_dealloc
        CALL ppm_alloc(cfld1,ldu,iopt,info)
        iopt = ppm_param_alloc_fit_preserve
        
        ldl(1) = 1; ldl(5) = 1;
        
        ldl(2) = 1-ghostsize(1)
        ldl(3) = 1-ghostsize(2)
        ldl(4) = 1-ghostsize(3)

        ldu(2) = ldu(2) + ghostsize(1)
        ldu(3) = ldu(3) + ghostsize(2)
        ldu(4) = ldu(4) + ghostsize(3)
        CALL ppm_alloc(field_rhs, ldl, ldu, iopt, info)
        
        maptype = ppm_param_map_global
        CALL ppm_map_field(field_rhs, 3, topo_id, t_mesh_id(1),mesh_id,&
             &             ghostsize, maptype, info)
        maptype = ppm_param_map_push
        CALL ppm_map_field(field_rhs, 3, topo_id, t_mesh_id(1),mesh_id,&
             &             ghostsize, maptype, info)
        maptype = ppm_param_map_send
        CALL ppm_map_field(field_rhs, 3, topo_id, t_mesh_id(1),mesh_id,&
             &             ghostsize, maptype, info)
        maptype = ppm_param_map_pop
        CALL ppm_map_field(field_rhs, 3, topo_id, t_mesh_id(1),mesh_id,&
             &             ghostsize, maptype, info)
        WRITE(msg,*) 'after  map field_rhs is ',UBOUND(field_rhs)
        IF(verbose) CALL ppm_write(ppm_rank,pnm,msg,info)

        !-----------------------------------------------------------------------
9999    CONTINUE

        !-----------------------------------------------------------------------
        ! clean up
        !-----------------------------------------------------------------------
        ! done
        !-----------------------------------------------------------------------
        iopt = ppm_param_dealloc
        CALL ppm_alloc(ndata,ldu,iopt,info)
        CALL ppm_alloc(isublist,ldu,iopt,info)
        CALL ppm_alloc(min_sub,ldu,iopt,info)
        CALL ppm_alloc(max_sub,ldu,iopt,info)
        CALL ppm_alloc(istart,ldu,iopt,info)
        CALL ppm_alloc(sub2proc,ldu,iopt,info)
        CALL ppm_alloc(sub_cost,ldu,iopt,info)
        
        
        DEALLOCATE(isublist_xy,ndata_xy)
        CALL substop(pnm,t0,info)
      END SUBROUTINE __PPM_FFT_GENERIC

!#undef __PPM_FFT_GENERIC
!#undef __WHAT
        
        
