!------------------------------------------------------------------------------!
!* filename: wvic_enforcetv2                                                  *!
!* project : ppm                                                              *!
!* purpose : initial condition for crow instability                           *!
!*         :                                                                  *!
!* author  : Michael Bergdorf  Philippe Chatelain                             *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Fri Dec  3 09:45:17 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: $
!
!
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
! = wvic ENFORCE TRAILING VORTICES CONDITION 2
!  initialize vorticity here
!------------------------------------------------------------------------------!
SUBROUTINE wvic_enforcetv2

  USE module_wvic
  USE ppm_module_data
  USE ppm_module_write
  USE ppm_module_map
  USE MPI
  USE, INTRINSIC :: ISO_C_BINDING
  INCLUDE 'fftw3.f03'

  !----------------------------------------------------------------------------!
  ! localities: noise stuff
  INTEGER, PARAMETER       :: mkd = KIND(1.0d0)
  INTEGER(mkd)             :: an, a32
  REAL(mk), DIMENSION(500) :: aky,bky,akz,bkz
  INTEGER                  :: nk,kk
  REAL(mk)                 :: rnk, amp
  INTEGER, DIMENSION(3)    :: ilo, ihi
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: geometry stuff
  REAL(mk)                 :: a0, b0, gamma
  REAL(mk)                 :: ty_center, ty_vtx1, ty_vtx2, r_vtx1, r_vtx2
  REAL(mk)                 :: tx_center, tz_vtx1, tz_vtx2
  REAL(mk)                 :: ty_vtx3, ty_vtx4, tx_vtx3, tx_vtx4
  REAL(mk)                 :: tx_vtx1, tx_vtx2, r_vtx3, r_vtx4
  REAL(mk)                 :: fac1, ra1, tx, ty, tz, noise, txp, typ, tzp, omega_x
  REAL(mk)                 :: dnoise, dtxp, dtyp
  REAL(mk)                 :: fac2, ra2
  REAL(mk), DIMENSION(3)   :: length
  REAL(mk), DIMENSION(4)   :: twp
  REAL(mk)                 :: ink_a0, ink_b0
  REAL(mk)                 :: tz_ink1, ty_ink1
  REAL(mk)                 :: tz_ink2, ty_ink2, pi
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  ! localities: stuff stuff
  CHARACTER(len=256)       :: msg
  INTEGER                  :: info, maptype
  INTEGER                  :: ig,jg,kg,i,j,k,isub,isubl,kmax,kl
  REAL(mk), PARAMETER :: onethird         = 0.3333333333333332_mk
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: mg stuff
  INTEGER,  DIMENSION(6)     :: ibcdef
  INTEGER,  DIMENSION(1,6)   :: ibcdef_v
  REAL(mk), DIMENSION(1,1,1) :: ibcvalue
  REAL(mk), DIMENSION(1,1,1,1) :: ibcvalue_v
  REAL(mk)                     :: cresid
  !----------------------------------------------------------------------------!

  !----------------------------------------------------------------------------!
  ! localities: helmholtz stuff
  REAL(mk)                     :: fac3,fac4,fac5,fac6,ldiv
  REAL(mk), DIMENSION(:,:,:,:), POINTER :: divergence, potential
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: divergence_v, potential_v
  
  INTEGER, DIMENSION(3)        :: ftopo_id
  INTEGER, DIMENSION(4)        :: fmesh_id
  
  !----------------------------------------------------------------------------!
  ! localities: Div free condition stuff
  REAL(mk), DIMENSION(:,:), POINTER :: upstream_planel, upstream_plane
  REAL(mk), DIMENSION(:,:,:), POINTER :: upstream_wpl, upstream_wp
  INTEGER, DIMENSION(:,:), POINTER     :: istartl
  INTEGER                           :: send_size, recv_size

!define SOLVE_WFISHPACK
#define SOLVE_WFFTW
#ifdef SOLVE_WFISHPACK
  REAL(mk)                          :: a,b,c,d, pertrbfp, lambda
  REAL(mk), DIMENSION(:), POINTER   :: workfp
  REAL(mk), DIMENSION(1)            :: dum
  INTEGER                           :: idimf, nfp, mfp, mbdcnd, nbdcnd
#endif

#ifdef SOLVE_WFFTW
  COMPLEX(mk), DIMENSION(:,:), POINTER :: upstream_planec
  INTEGER                                       :: idist, odist, frank
  INTEGER                                       :: istride, ostride
  INTEGER,  DIMENSION(2)                        :: fnx, iembed, oembed
  INTEGER                                       :: howmany
  TYPE(C_PTR)                                   :: plan
  REAL(mk)                                      :: pi2_lx,pi2_ly
  REAL(mk)                                      :: kinv
  INTEGER , DIMENSION(2)                        :: lda_end
  REAL(mk)                                      :: kx,ky,ivol
#endif

  INTEGER                           :: ip1, im1, jp1, jm1
  INTEGER :: c1,c2,c3,c4,c5,c6
  ! MPI comm status
  INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: status
  INTEGER :: dbgfile

  !----------------------------------------------------------------------------!
  ! = WHAT HAS TO BE DONE? =
  ! - generate noise (angle of dislocation and amplitute of dislocation)
  ! - calculate the position of the vortices using b0 (intervortex distance)
  ! - initialize the vortices using a0 (size of vortices)
  ! - reproject vorticity to get a divergence free field
  !----------------------------------------------------------------------------!

   WRITE(msg,'(A,I5.5,A,I1.1)') ' Enforcement at ',itime,' and ',g_istage
   IF ((verbose).AND.(rank.EQ.0)) CALL ppm_write(rank,'wvic_enforcetv2',msg,info)

  info = 0
  IF(nsublist.GT.1) THEN
     WRITE(msg,*) 'random init only for 1 sub/cpu'
     CALL ppm_write(rank,'wvic_enforcetv2',msg,info)
     STOP ! brutal
  END IF

  isubl = isublist(1)
  ilo = istart(:,isubl)
  ihi = ilo + ndata(:,isubl)-1

  !----------------------------------------------------------------------------!
  ! jump start the random number generator
  !----------------------------------------------------------------------------!
  a32 = 0_mkd
  a32 = IBSET(a32,32)
  an  = 1_mkd
  an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
   !----------------------------------------------------------------------------!
  ! compute the random amplitudes
  !----------------------------------------------------------------------------!
  IF (wvic_noise_basemode.EQ.0.0_mk) THEN
     wvic_noise_basemode = max_physg(3) - min_physg(3)
  END IF
  
  nk = wvic_noise_nmodes
  rnk = 1.0_mk/(2.0_mk*REAL(nk,mk))
  amp = wvic_noise_amp * rnk * trailvortex_b1
  DO i=1,nk

     aky(i) = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     bky(i) = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     akz(i)   = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)
     bkz(i)   = REAL(MOD(an/65536_mkd,32768_mkd),mk)/REAL(32769,mk)
     an  = MOD((1103515245_mkd*an + 12345_mkd),a32)

  END DO

  !----------------------------------------------------------------------------!
  ! derived parameters - vorticity
  !----------------------------------------------------------------------------!
  ty_center = 0.5_mk*(min_physg(2) + max_physg(2))
  ty_vtx1   = ty_center - 0.5_mk * trailvortex_b1
  ty_vtx2   = ty_center + 0.5_mk * trailvortex_b1
  ty_vtx3   = ty_center - 0.5_mk * trailvortex_b2
  ty_vtx4   = ty_center + 0.5_mk * trailvortex_b2
  tx_center = 0.5_mk*(min_physg(1) + max_physg(1))
  tx_vtx1   = tx_center
  tx_vtx2   = tx_center
  tx_vtx3   = tx_center + trailvortex_z12
  tx_vtx4   = tx_center + trailvortex_z12
  
  fac1      = trailvortex_gamma/(trailvortex_a1**2 * M_PI)
  ra1       = 1.0_mk/trailvortex_a1**2
  fac2      = trailvortex_r * trailvortex_gamma/(trailvortex_a2**2 * M_PI)
  ra2       = 1.0_mk/trailvortex_a2**2
  length    = max_physg - min_physg
  !----------------------------------------------------------------------------!

  
  !----------------------------------------------------------------------------!
  ! First impose omegaz on upstream plane
  !----------------------------------------------------------------------------!
  DO isub=1,nsublist
     isubl = isublist(isub)
	 IF (min_sub(3,isubl).LT.(min_physg(3)+0.5*dz)) THEN
        DO k=1-ghostsize(3),1
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)
			  tx = min_sub(1,isubl) + REAL(i-1,mk)*dx
              ty = min_sub(2,isubl) + REAL(j-1,mk)*dy
			  tz = min_sub(3,isubl) + REAL(k-1,mk)*dz
              
			  !-----------------------------------------------------
              !  noise
              !-----------------------------------------------------
              omega_x = 2.0_mk*M_PI*(tz-u_infty(3)*time)/wvic_noise_basemode
                            noise   = 0.0_mk
              DO kk=1,nk
                 noise = noise &
                      & + SIN(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*akz(kk))) &
                      & + COS(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*bkz(kk)))
              END DO
              txp = noise * amp
              noise   = 0.0_mk
              DO kk=1,nk
                 noise = noise &
                      & + SIN(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*aky(kk))) &
                      & + COS(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*bky(kk)))
              END DO
              typ = noise * amp
              tzp = 0.0_mk
              
              twp = 0.0_mk
              
              !----------------------------------------------------------------!
              ! compute vorticity - vtx1 in ...
              !----------------------------------------------------------------!
              ! center box
			  IF(nk.EQ.1) typ = 0.0_mk
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tx+txp-tx_vtx1)**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! left neighbor
              r_vtx1 = (ty+typ-ty_vtx1-length(2))**2 + (tx+txp-tx_vtx1)**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! right neighbor
              r_vtx1 = (ty+typ-ty_vtx1+length(2))**2 + (tx+txp-tx_vtx1)**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! bottom neighbor
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tx+txp-tx_vtx1-length(1))**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              ! top neighbor
              r_vtx1 = (ty+typ-ty_vtx1)**2 + (tx+txp-tx_vtx1+length(1))**2
              twp(3) =  twp(3) + fac1 * EXP(-r_vtx1*ra1)
              !----------------------------------------------------------------!
              ! compute vorticity - vtx2 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx2 = (ty+txp-ty_vtx2)**2 + (tx+typ-tx_vtx2)**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! left neighbor
              r_vtx2 = (ty+txp-ty_vtx2-length(2))**2 + (tx+typ-tx_vtx2)**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! right neighbor
              r_vtx2 = (ty+txp-ty_vtx2+length(2))**2 + (tx+typ-tx_vtx2)**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! bottom neighbor
              r_vtx2 = (ty+txp-ty_vtx2)**2 + (tx+typ-tx_vtx2-length(1))**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              ! top neighbor
              r_vtx2 = (ty+txp-ty_vtx2)**2 + (tx+typ-tx_vtx2+length(1))**2
              twp(3) =  twp(3) - fac1 * EXP(-r_vtx2*ra1)
              !----------------------------------------------------------------!
              ! compute vorticity - vtx3 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx3 = (ty-txp-ty_vtx3)**2 + (tx-typ-tx_vtx3)**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! left neighbor
              r_vtx3 = (ty-txp-ty_vtx3-length(2))**2 + (tx-typ-tx_vtx3)**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! right neighbor
              r_vtx3 = (ty-txp-ty_vtx3+length(2))**2 + (tx-typ-tx_vtx3)**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! bottom neighbor
              r_vtx3 = (ty-txp-ty_vtx3)**2 + (tx-typ-tx_vtx3-length(1))**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              ! top neighbor
              r_vtx3 = (ty-txp-ty_vtx3)**2 + (tx-typ-tx_vtx3+length(1))**2
              twp(3) =  twp(3) + fac2 * EXP(-r_vtx3*ra2)
              !----------------------------------------------------------------!
              ! compute vorticity - vtx4 in ...
              !----------------------------------------------------------------!
              ! center box
              r_vtx4 = (ty-typ-ty_vtx4)**2 + (tx-txp-tx_vtx4)**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! left neighbor
              r_vtx4 = (ty-typ-ty_vtx4-length(2))**2 + (tx-txp-tx_vtx4)**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! right neighbor
              r_vtx4 = (ty-typ-ty_vtx4+length(2))**2 + (tx-txp-tx_vtx4)**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! bottom neighbor
              r_vtx4 = (ty-typ-ty_vtx4)**2 + (tx-txp-tx_vtx4-length(1))**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)
              ! top neighbor
              r_vtx4 = (ty-typ-ty_vtx4)**2 + (tx-txp-tx_vtx4+length(1))**2
              twp(3) =  twp(3) - fac2 * EXP(-r_vtx4*ra2)

              !-----------------------------------------------------
              !  SET FIELD VORTICITY VALUES
              !-----------------------------------------------------
              field_wp(1,i,j,k,isub) = twp(1)
              field_wp(2,i,j,k,isub) = twp(2)
              field_wp(3,i,j,k,isub) = twp(3)
              !-----------------------------------------------------
              END DO
           END DO
        END DO
     END IF
  END DO
  

  IF (upstream_rank.NE.MPI_UNDEFINED) THEN
     isub=1
     isubl = isublist(isub)
     
     ALLOCATE(upstream_planel(ndata(1,isubl),ndata(2,isubl)),stat=info)
     ALLOCATE(upstream_wpl(2,ndata(1,isubl),ndata(2,isubl)),stat=info)
     fac3 = -1.5_mk / dz
     fac4 =  2.0_mk / dz
     fac5 = -0.5_mk / dz
     DO k=1,1
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)
                 upstream_planel(i,j) = fac3 * field_wp(3,i,j,k,isub) + &
                                    &   fac4 * field_wp(3,i,j,k+1,isub) + &
                                    &   fac5 * field_wp(3,i,j,k+2,isub)
           END DO
        END DO
     END DO
     
     IF (upstream_rank.EQ.0) THEN
         ALLOCATE(upstream_plane(nx(1),nx(2)),stat=info)
         ALLOCATE(upstream_wp(2,nx(1),nx(2)),stat=info)
         ALLOCATE(istartl(3,n_upstream))
         
         upstream_plane(istart(1,isubl):istart(1,isubl)+ndata(1,isubl)-1, &
         & istart(2,isubl):istart(2,isubl)+ndata(2,isubl)-1) = &
         & upstream_planel(1:ndata(1,isubl),1:ndata(2,isubl))
         istartl(:,1) = istart(:,isubl)
         DO i=2,n_upstream
            recv_size = ndata(1,isubl)*ndata(2,isubl) ! Assuming they're all the same size!!
            CALL MPI_RECV(istartl(:,i),3,MPI_INTEGER, i-1, 0, upstream_comm, status,info)
            CALL MPI_RECV(upstream_planel,recv_size,mpi_prec,i-1,0,upstream_comm,status,info)
            upstream_plane(istartl(1,i):istartl(1,i)+ndata(1,isubl)-1,istartl(2,i):istartl(2,i)+ndata(2,isubl)-1) = &
            & upstream_planel
         END DO
         
         !dbgfile = 76
         !WRITE(msg,'(A,I2.2,A,I1.1)') 'upstrsrc',itime,'_',g_istage
         !OPEN(dbgfile,file=msg,iostat=info,action='WRITE')
         !WRITE(*,*) 'info = ', info
         !WRITE(dbgfile,*) upstream_plane
           
           
#ifdef SOLVE_WFISHPACK         
         a = min_physg(1); b = max_physg(1)
         c = min_physg(2); d = max_physg(2)
         nfp = nx(1)-1
         mfp = nx(2)-1
         mbdcnd = 0
         nbdcnd = 0
         lambda = 0.0
         idimf = mfp+1
         
         ALLOCATE(workfp(4*(nfp+1)+(13+INT(LOG(REAL(nfp+1))/LOG(2.0)))*(mfp+1)))
         CALL SYSTEM_CLOCK(c1,c2,c3)
         CALL hwscrt(a,b,mfp,mbdcnd,dum,dum,c,d,nfp,nbdcnd,dum,dum,&
                   & lambda,upstream_plane,idimf,pertrbfp,info)!,workfp)
         DEALLOCATE(workfp)
         CALL SYSTEM_CLOCK(c4,c5,c6)
         !PRINT*,'info:',info,', time',REAL(c4-c1)/REAL(c5)
#endif

#ifdef SOLVE_WFFTW
         ALLOCATE(upstream_planec((nx(1)-1)/2 +2,nx(2)),stat=info)

         length(1) = max_physg(1)-min_physg(1)
         length(2) = max_physg(2)-min_physg(2)

         frank      = 2
         fnx(1)    = nx(1)-1
         fnx(2)    = nx(2)-1
         iembed(1) = UBOUND(upstream_plane,1)
         iembed(2) = UBOUND(upstream_plane,2)
         oembed(1) = UBOUND(upstream_planec,1)
         oembed(2) = UBOUND(upstream_planec,2)
         istride   = 1
         ostride   = 1
         idist     = 1
         odist     = 1
         CALL sfftw_plan_dft_r2c(plan, frank, fnx(1), &
               & upstream_plane(1,1), &
               & upstream_planec(1,1), &
               & FFTW_ESTIMATE)
         CALL sfftw_execute(plan)
         CALL sfftw_destroy_plan(plan)
         DO j=1,nx(2)
            upstream_planec((nx(1)-1)/2+2,j) = upstream_planec(1,j)
         END DO
         DO i=1,nx(1)
            upstream_planec(i,nx(2)) = upstream_planec(i,1)
         END DO
         upstream_planec((nx(1)-1)/2+2,nx(2)) = &
                   & upstream_planec(1,1)
         
         pi2_Lx = 2.0_MK*ppm_pi_s/length(1)
         pi2_Ly = 2.0_MK*ppm_pi_s/length(2)
         ivol = 1.0_mk/REAL((nx(1)-1)*(nx(2)-1))
         lda_end(1)=(nx(1)-1)/2 +1
         lda_end(2)=(nx(2)-1)/2 +1
         DO j=1,nx(2)
            jg = j-1
            DO i=1,nx(1)
               ig = i-1
               IF(ig.EQ.0.AND.jg.EQ.0) THEN
                  upstream_planec(i,j) = 0.0_mk
               ELSE
                  IF(ig.GE.(lda_end(1))) THEN
                     ig = nx(1)-ig-1
                  END IF
                  IF(jg.GE.(lda_end(2))) THEN
                     jg = nx(2)-jg-1
                  END IF
                  kx=pi2_Lx*REAL(ig,MK)  
                  ky=pi2_Ly*REAL(jg,MK)
                  kinv = -1.0_mk/(kx*kx + ky*ky)*ivol
                  upstream_planec(i,j)=kinv&
                            &*upstream_planec(i,j)
                ENDIF
            END DO                 
        END DO
         
#endif         
         !WRITE(dbgfile,*) upstream_plane
         !CLOSE(dbgfile,iostat=info)
         
         fac3 = -0.5 / dx
         fac4 = -0.5 / dy
         DO j=1,nx(2)
            DO i=1,nx(1)
                 ip1 = MOD(i,nx(1)-1)+1
                 im1 = MOD(i-3+nx(1),nx(1)-1)+1
                 jp1 = MOD(j,nx(2)-1)+1
                 jm1 = MOD(j-3+nx(2),nx(2)-1)+1
                 
                 upstream_wp(1,i,j) = fac3 * (upstream_plane(ip1,j) - &
                                           & upstream_plane(im1,j))
                 upstream_wp(2,i,j) = fac4 * (upstream_plane(i,jp1) - &
                                           & upstream_plane(i,jm1))
           END DO
        END DO
        WRITE(msg,*) ' Found the transversal vorticities ', &
             MINVAL(upstream_wp(1,:,:)), &
             MAXVAL(upstream_wp(1,:,:)), &
             MINVAL(upstream_wp(2,:,:)), &
             MAXVAL(upstream_wp(2,:,:))
        IF ((verbose).AND.(rank.EQ.0)) &
             CALL ppm_write(rank,'wvic_enforcetv2',msg,info)
        DO i=2,n_upstream
            send_size = ndata(1,isubl)*ndata(2,isubl)*2 ! Assuming they're all the same size!!
            upstream_wpl = upstream_wp(:,istartl(1,i):istartl(1,i)+ndata(1,isubl)-1,istartl(2,i):istartl(2,i)+ndata(2,isubl)-1)
            CALL MPI_SEND(upstream_wpl,send_size,mpi_prec,i-1,0,upstream_comm,info)
         END DO
         upstream_wpl = upstream_wp(:,istartl(1,1):istartl(1,1)+ndata(1,isubl)-1,istartl(2,1):istartl(2,1)+ndata(2,isubl)-1)
         DEALLOCATE(upstream_plane,upstream_wp,istartl)
     ELSE
         send_size = ndata(1,isubl)*ndata(2,isubl)
         CALL MPI_SEND(istart(:,isubl),3, MPI_INTEGER, 0, 0, upstream_comm,info)
         CALL MPI_SEND(upstream_planel,send_size, mpi_prec, 0, 0, upstream_comm,info)
         
         recv_size = ndata(1,isubl)*ndata(2,isubl)*2
         CALL MPI_RECV(upstream_wpl,recv_size,mpi_prec,0,0,upstream_comm,status,info)
     END IF
     ! TO DO ASSIGN upstream x and y to field_wp
     DO k=1,1
           DO j=1,ndata(2,isubl)
              DO i=1,ndata(1,isubl)
                 field_wp(1:2,i,j,k,isub) = upstream_wpl(1:2,i,j)
           END DO
        END DO
     END DO
     DEALLOCATE(upstream_wpl,upstream_planel)
  END IF
END SUBROUTINE wvic_enforcetv2
