      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_fmm_traverse
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Recursive routine to traverse the tree.
      !                 Called by the ppm_fmm_expansion subroutine.
      !
      !  Input        : root         (I) index of the root box.
      !                 order        (I) expansion order 
      !                 prec         (I) not used dummy argument
      !                                  to determine precision
      !
      !  Input/output :     
      !
      !  Output       : 
      !                 info         (I) return status. 0 upon success
      !
      !  Remarks      :  
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_traverse.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.12  2005/09/19 13:03:30  polasekb
      !  code cosmetics
      !
      !  Revision 1.11  2005/09/05 06:27:59  polasekb
      !  checking if box is on proc
      !
      !  Revision 1.10  2005/08/23 14:11:55  polasekb
      !  fixed the sign for the Dnm
      !
      !  Revision 1.9  2005/08/08 13:37:08  polasekb
      !  nullify some data pointers
      !
      !  Revision 1.8  2005/08/04 16:03:17  polasekb
      !  moved data allocation to init
      !
      !  Revision 1.7  2005/07/29 12:35:54  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.6  2005/07/27 14:59:57  polasekb
      !  now computing centerofbox from children
      !
      !  Revision 1.5  2005/07/25 14:38:19  polasekb
      !  changed call to subroutine,
      !  now saving constants in module_data_fmm file
      !
      !  Revision 1.4  2005/07/21 13:17:03  polasekb
      !  bugfix in do-loop
      !
      !  Revision 1.3  2005/06/02 14:24:01  polasekb
      !  removed variable totalmass
      !  bugfix calling cart2sph
      !
      !  Revision 1.2  2005/05/27 12:44:28  polasekb
      !  removed some debug output
      !
      !  Revision 1.1  2005/05/27 07:59:30  polasekb
      !  initial implementation
      !  TODO: remove debug outputs
      !
      !  Revision 0  2004/12/02 15:59:14 polasekb
      !  start
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE ppm_fmm_traverse_s(root,order,prec,info)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE ppm_fmm_traverse_d(root,order,prec,info)
#endif

!      RECURSIVE SUBROUTINE ppm_fmm_traverse(root,order,info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop 
      USE ppm_module_util_cart2sph
      USE ppm_module_write

      
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      USE MPI
#endif
      

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: root
      INTEGER                 , INTENT(IN   ) :: order
      REAL(MK)                , INTENT(IN   ) :: prec !dummy arg for prec.
      INTEGER                 , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      LOGICAL                              :: onlocalproc
      INTEGER                              :: m,n,l,j,i,iopt
      INTEGER                              :: fir,las,box,istat
      INTEGER                              :: first,last
      REAL(MK)                             :: sine,cosine,val,prod
      REAL(MK)                             :: angle,reci,t0
      REAL(MK)                             :: dx,dy,dz,tmp  
      REAL(MK),DIMENSION(:  ),POINTER      :: box_rho => null()
      REAL(MK),DIMENSION(:  ),POINTER      :: box_theta => null()
      REAL(MK),DIMENSION(:  ),POINTER      :: box_phi => null()
      REAL(MK),DIMENSION(:  ),POINTER      :: fracfac => null()
      REAL(MK),DIMENSION(:  ),POINTER      :: totalmass => null()
      REAL(MK),DIMENSION(:  ),POINTER      :: radius => null()
      REAL(MK),DIMENSION(:,:),POINTER      :: Pnm => null()
      REAL(MK),DIMENSION(:,:),POINTER      :: Anm => null()
      REAL(MK),DIMENSION(:,:),POINTER      :: sqrtfac => null()
      REAL(MK),DIMENSION(:,:),POINTER      :: centerofbox => null()
      COMPLEX(MK)                          :: csum
      COMPLEX(MK),PARAMETER                :: CI=(0.0_MK,1.0_MK)
      COMPLEX(MK),DIMENSION(:,:  ),POINTER :: Inner => null()
      COMPLEX(MK),DIMENSION(:,:  ),POINTER :: Ynm => null()
      COMPLEX(MK),DIMENSION(:,:,:),POINTER :: expansion => null()
      CHARACTER(LEN=ppm_char)              :: cbuf

      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------

      CALL substart('ppm_fmm_traverse',t0,info)
      !-------------------------------------------------------------------------
      ! Checking precision and pointing tree data to correct variables
      !-------------------------------------------------------------------------
        
#if     __KIND == __SINGLE_PRECISION
      centerofbox  => centerofbox_s
      totalmass  => totalmass_s
      expansion    => expansion_s
      radius    => radius_s
      Anm          => Anm_s
      sqrtfac      => sqrtfac_s
      fracfac      => fracfac_s
      Ynm          => Ynm_s
      Pnm          => Pnm_s
      box_rho      => rho_s
      box_theta    => theta_s
      box_phi      => phi_s
      Inner        => Inner_s
#else
      centerofbox  => centerofbox_d
      totalmass  => totalmass_d
      expansion    => expansion_d
      radius    => radius_d
      Anm          => Anm_d
      sqrtfac      => sqrtfac_d
      fracfac      => fracfac_d
      Ynm          => Ynm_d
      Pnm          => Pnm_d
      box_rho      => rho_d
      box_theta    => theta_d
      box_phi      => phi_d
      Inner        => Inner_d
#endif
      !-------------------------------------------------------------------------
      ! Checking if current root is a leaf,
      ! if yes, no shifting has to be done
      !-------------------------------------------------------------------------

      IF (nchld(root) .EQ. 0) THEN 
          IF (ppm_debug .GT. 0) THEN
            CALL ppm_write(ppm_rank,'ppm_fmm_traverse','leaf box',info)
          ENDIF

      !-------------------------------------------------------------------------
      ! Current root has more children, call routine recursively
      !-------------------------------------------------------------------------

      ELSE
         !----------------------------------------------------------------------
         ! Check if box on local proc., if not, no computation for this box
         !----------------------------------------------------------------------
         onlocalproc = .FALSE.
         IF (nbpl(blevel(root)) .LT. ppm_nproc) THEN
            !level topology not defined
            onlocalproc = .TRUE.
         ELSE
           DO i=1,ppm_nsublist(ppm_internal_topoid(blevel(root)))
              IF (ppm_boxid(ppm_isublist(i,ppm_internal_topoid(blevel(root))), &
                  & blevel(root)) .EQ. root) THEN
                 onlocalproc = .TRUE.
                 EXIT
              ENDIF
           ENDDO
         ENDIF
         IF (onlocalproc) THEN
         DO i=1,nchld(root)
            CALL ppm_fmm_traverse(child(i,root),order,prec,info)
            IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_fmm_traverse', &
              &              'Called traverse',info)
            ENDIF
         ENDDO

        !-----------------------------------------------------------------------
        ! Now all children have been called recursively, computation can start
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  Compute centerofbox and totalmass
        !-----------------------------------------------------------------------

        fir = child(1,root)
        las = child((nchld(root)),root)
        totalmass(root) = SUM(totalmass(fir:las))
        centerofbox(1,root) = SUM(centerofbox(1,fir:las)*totalmass(fir:las))/ &
        &                     totalmass(root)
        centerofbox(2,root) = SUM(centerofbox(2,fir:las)*totalmass(fir:las))/ &
        &                     totalmass(root)
        centerofbox(3,root) = SUM(centerofbox(3,fir:las)*totalmass(fir:las))/ &
        &                     totalmass(root)
        !-----------------------------------------------------------------------
        !  Init Variables
        !-----------------------------------------------------------------------
        box_rho   = (0.0_MK)
        box_phi   = (0.0_MK)
        box_theta = (0.0_MK)

        !-----------------------------------------------------------------------
        !  Compute spherical coordinates of child boxes
        !-----------------------------------------------------------------------
        CALL ppm_util_cart2sph(centerofbox(1,fir:las),centerofbox(2,fir:las), &
              & centerofbox(3,fir:las),nchld(root), &
              & centerofbox(1,root),centerofbox(2,root),centerofbox(3,root), &
              & box_rho(1:nchld(root)),box_theta(1:nchld(root)), &
              & box_phi(1:nchld(root)),info)
        IF (info .NE. 0) THEN
            CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_traverse', &
                 & 'Failed calling util_cart2sph',__LINE__,info)
        ENDIF
        DO i=1,nchld(root)
           first = lhbx(1,child(i,root))
           last  = lhbx(2,child(i,root))
           IF (last-first+1 .EQ. 0) CYCLE
           !--------------------------------------------------------------------
           !  Init Variables
           !--------------------------------------------------------------------
           Ynm       = (0.0_MK,0.0_MK)
           Pnm       = (0.0_MK,0.0_MK)
           Inner     = (0.0_MK,0.0_MK)
           !--------------------------------------------------------------------
           !  Compute radius
           !--------------------------------------------------------------------
           
           box = child(i,root)
           dx = centerofbox(1,box) - centerofbox(1,root)
           dy = centerofbox(2,box) - centerofbox(2,root)
           dz = centerofbox(3,box) - centerofbox(3,root)
           tmp = SQRT(dx**2 + dy**2 + dz**2) + radius(box)
           IF (tmp .GT. radius(root)) THEN
              radius(root) = tmp
           ENDIF

           !--------------------------------------------------------------------
           ! Compute Legendre polynomial
           !--------------------------------------------------------------------
           reci      = 1.0_MK/box_rho(i)
           sine      = SIN(box_theta(i))
           cosine    = COS(box_theta(i))
           val       = -sine
           prod      = 1.0_MK
           DO m=0,order
             Pnm(m,m) = fracfac(m)*prod
             prod     = prod * val
           ENDDO
           DO m=0,order-1
              Pnm(m+1,m) = cosine*REAL(2*m + 1,MK)*Pnm(m,m)
           ENDDO
           DO n=2,order
              val = cosine*REAL(2*n-1,MK)
              DO m=0,n-1
                 Pnm(n,m)=(val*Pnm(n-1,m)-REAL(n+m-1,MK)* &
                          Pnm(n-2,m))/REAL(n-m,MK)
              ENDDO
           ENDDO
           IF (ppm_debug .GT. 0) THEN
             CALL ppm_write(ppm_rank,'ppm_fmm_traverse','Computed Pnm',info)
           ENDIF
           !--------------------------------------------------------------------
           ! Compute Ynm(n,m) and Ynm(n,-m)
           !--------------------------------------------------------------------
           DO n=0,order
              m = 0
              angle = REAL(m,MK)*box_phi(i)
              Ynm(n,m) = sqrtfac(n,m)*Pnm(n,m)* &
                         & CMPLX(COS(angle),SIN(angle))
              DO m=1,n
                 angle     = REAL(m,MK)*box_phi(i)
                 Ynm(n,m)  = sqrtfac(n,m)*Pnm(n,m)* &
                            & CMPLX(COS(angle),SIN(angle))
                 Ynm(n,-m) = CONJG(Ynm(n,m))
              ENDDO
           ENDDO
           IF (ppm_debug .GT. 0) THEN
             CALL ppm_write(ppm_rank,'ppm_fmm_traverse','Computed Ynm',info)
           ENDIF
           !--------------------------------------------------------------------
           ! Compute Inner expansion
           !--------------------------------------------------------------------
           DO n=0,order
              DO m=-n,n
                 Inner(n,m)=CI**(-ABS(m))*Anm(n,m)*box_rho(i)**n*Ynm(n,m)
              ENDDO
           ENDDO
           IF (ppm_debug .GT. 0) THEN
             CALL ppm_write(ppm_rank,'ppm_fmm_traverse','Computed Inner',info)
           ENDIF

           !--------------------------------------------------------------------
           ! Compute Dnm(n,m) = expansion coefficient
           !--------------------------------------------------------------------
           DO l=0,order
             DO j=-l,l
                csum = (0.0_MK,0.0_MK)
                DO n=0,l !order
                   DO m=MAX(j+n-l,-n),MIN(j+l-n,n) !-n,n
                      csum = csum + (-1)**(l-n)*Inner((l-n),(j-m)) &
                      &      *expansion(child(i,root),n,m)
                   ENDDO
                ENDDO
                ! add to expansion in fmm_module_data file
                expansion(root,l,j) = expansion(root,l,j) + csum
             ENDDO
           ENDDO
           IF (ppm_debug .GT. 0) THEN
             CALL ppm_write(ppm_rank,'ppm_fmm_traverse','Computed Dnm',info)
           ENDIF
         ENDDO
      ENDIF !on local proc
      ENDIF

      !-------------------------------------------------------------------------
      !  Nullify data pointers
      !-------------------------------------------------------------------------

      NULLIFY(centerofbox)
      NULLIFY(radius)
      NULLIFY(expansion)
      NULLIFY(Anm)
      NULLIFY(sqrtfac)
      NULLIFY(fracfac)
      NULLIFY(Ynm)
      NULLIFY(Pnm)
      NULLIFY(Inner)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999    CONTINUE
        CALL substop('ppm_fmm_traverse',t0,info)
        RETURN

#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fmm_traverse_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fmm_traverse_d
#endif
      

