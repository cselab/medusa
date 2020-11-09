      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_fmm_expansion
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Compute the expansions of the leaf boxes of the created
      !                 tree structure.
      !                 Calls ppm_fmm_traverse, which traverses the tree and 
      !                 shifts the expansions up the tree.
      !                 
      !
      !  Input        : xp(:,:)      (F) field particle positions
      !                 wp(:)        (F) field particle strengths
      !                 min_dom(:)   (F) the minimum coordinate of the
      !                                  domain
      !                 max_dom(:)   (F) the maximum coordinate of the
      !                                  domain      
      !                 Nm(:)        (I) number of grid points in the
      !                                  global mesh. (0,0,0) if there is
      !                                  no mesh.
      !                 Np           (I) the number of field particles.
      !                 order        (I) expansion order
      !
      !  Input/output :     
      !
      !  Output       :
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      :  
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_expansion.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.16  2006/02/03 09:39:12  ivos
      !  Changed the expansion loop to allow vectorization.
      !
      !  Revision 1.15  2005/09/19 13:03:27  polasekb
      !  code cosmetics
      !
      !  Revision 1.14  2005/09/11 11:43:59  polasekb
      !  moved mapping and second tree call to init
      !
      !  Revision 1.13  2005/08/25 13:52:28  polasekb
      !  mapping of boxpart added
      !
      !  Revision 1.12  2005/08/23 14:22:43  polasekb
      !  corrected to xpunord and wpunord
      !
      !  Revision 1.11  2005/08/23 14:11:38  polasekb
      !  fixed the formula for the Cnm
      !
      !  Revision 1.10  2005/08/11 15:13:15  polasekb
      !  now using maxboxcost from the data file
      !
      !  Revision 1.9  2005/08/08 13:33:41  polasekb
      !  init some more variables
      !
      !  Revision 1.8  2005/08/04 16:01:01  polasekb
      !  moved some data allocation to init
      !
      !  Revision 1.7  2005/07/29 12:35:27  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.6  2005/07/27 15:01:27  polasekb
      !  changed tree variables
      !
      !  Revision 1.5  2005/07/25 14:39:54  polasekb
      !  adapted to new constants saved in the
      !  ppm_module_data_fmm file
      !  adapted function call to ppm_fmm_traverse
      !
      !  Revision 1.4  2005/07/21 13:20:25  polasekb
      !  nullify lpdx, lhbx pointers
      !
      !  Revision 1.3  2005/06/02 14:18:39  polasekb
      !  removed variable totalmass
      !  changed some comments
      !
      !  Revision 1.2  2005/05/27 12:45:29  polasekb
      !  removed debug output
      !
      !  Revision 1.1  2005/05/27 07:56:40  polasekb
      !  initial implementation
      !
      !  Revision 0  2004/11/16 15:59:14 polasekb
      !  start
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fmm_expansion_s(xpunord,wpunord,min_dom,max_dom,Nm,Np, &
                 & order,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fmm_expansion_d(xpunord,wpunord,min_dom,max_dom,Nm,Np, &
                 & order,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_fmm_traverse
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

      
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK),DIMENSION(:,:),POINTER       :: xpunord
      REAL(MK),DIMENSION(:)  ,POINTER       :: wpunord
      REAL(MK),DIMENSION(:)  ,INTENT(IN   ) :: min_dom,max_dom
      INTEGER ,DIMENSION(:  ),INTENT(IN   ) :: Nm
      INTEGER                ,INTENT(INOUT) :: Np
      INTEGER                ,INTENT(IN   ) :: order
      INTEGER                ,INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------

      LOGICAL ,DIMENSION(3)                 :: fixed
      INTEGER                               :: iopt,i,j,k,m,n
      INTEGER                               :: root,isub,topoid
      INTEGER                               :: nsublist,in_topoid
      INTEGER                               :: first,last,box,nrpart
      INTEGER ,DIMENSION(2)                 :: ldu2 
      REAL(MK)                              :: sine,cosine,val,prod 
      REAL(MK)                              :: angle,reci 
      REAL(MK)                              :: t0,x0,y0,z0
      REAL(MK)                              :: dx,dy,dz,dist
      REAL(MK),DIMENSION(:)  , POINTER      :: rho,theta,phi,radius 
      REAL(MK),DIMENSION(:)  , POINTER      :: fac,fracfac,boxcost 
      REAL(MK),DIMENSION(:,:), POINTER      :: Anm,Pnm,sqrtfac,xp 
      REAL(MK),DIMENSION(:,:), POINTER      :: centerofbox
      REAL(MK),DIMENSION(:,:), POINTER      :: min_box,max_box
      REAL(MK),DIMENSION(:,:), POINTER      :: min_sub,max_sub
      COMPLEX(MK),PARAMETER                 :: CI=(0.0_MK,1.0_MK)
      COMPLEX(MK)                           :: temp
      COMPLEX(MK),DIMENSION(:,:),  POINTER  :: Cnm,Ynm
      COMPLEX(MK),DIMENSION(:,:,:),POINTER  :: expansion
      CHARACTER(LEN=ppm_char)               :: cbuf

      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------

      CALL substart('ppm_fmm_expansion',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN  
         DO i=1,ppm_dim
            IF (order .LT. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_expansion',   &
      &               'expansion order must be >= 0 !',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      ! Checking precision and pointing tree data to correct variables
      !-------------------------------------------------------------------------
        
#if   __KIND == __SINGLE_PRECISION
      min_box      => min_box_s
      max_box      => max_box_s
      boxcost      => boxcost_s
      centerofbox  => centerofbox_s
      radius       => radius_s
      expansion    => expansion_s
      Anm          => Anm_s
      sqrtfac      => sqrtfac_s
      fracfac      => fracfac_s
      fac          => fac_s
      Ynm          => Ynm_s
      Pnm          => Pnm_s
      Cnm          => Cnm_s
      rho          => rho_s
      theta        => theta_s
      phi          => phi_s
#else
      min_box      => min_box_d
      max_box      => max_box_d
      boxcost      => boxcost_d
      centerofbox  => centerofbox_d
      radius       => radius_d
      expansion    => expansion_d
      Anm          => Anm_d
      sqrtfac      => sqrtfac_d
      fracfac      => fracfac_d
      fac          => fac_d
      Ynm          => Ynm_d
      Pnm          => Pnm_d
      Cnm          => Cnm_d
      rho          => rho_d
      theta        => theta_d
      phi          => phi_d
#endif

      !-------------------------------------------------------------------------
      ! Allocating xp particles sorted according to tree
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit

      ldu2(1) = ppm_dim
      ldu2(2) = Np
      CALL ppm_alloc(xp,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expansion', &
      &       'error allocating xp',__LINE__,info)
      GOTO 9999
      ENDIF  
      
      !-------------------------------------------------------------------------
      ! Init variables
      !-------------------------------------------------------------------------

      fac    = (0.0_MK)
      rho    = (0.0_MK)
      phi    = (0.0_MK)
      theta  = (0.0_MK)
      topoid = nlevel

      !-------------------------------------------------------------------------
      ! Computing some constants for expansion coefficients
      ! and
      ! compute factorial function up to order*2
      !-------------------------------------------------------------------------

      fac(0) = 1
      DO i=1,order*2
         fac(i) = fac(i-1)*REAL(i,MK)
      ENDDO
      DO n=0,order
         DO m=-n,n
            sqrtfac(n,m) = SQRT(fac(n-ABS(m))/fac(n+ABS(m)))
            Anm(n,m)     = (-1)**n/SQRT(fac(n-m)*fac(n+m))
         ENDDO
      ENDDO

      DO m=0,order
         fracfac(m) = fac(2*m)/(2.**m*fac(m))
      ENDDO
      !-------------------------------------------------------------------------
      ! Computing the expansions of leafs -> particle-particle
      !-------------------------------------------------------------------------

      in_topoid = ppm_internal_topoid(topoid)
      DO i=1,ppm_nsublist(in_topoid)
         !----------------------------------------------------------------------
         !initializing arrays
         !----------------------------------------------------------------------
         Cnm = (0.0_MK,0.0_MK)
         Ynm = (0.0_MK,0.0_MK)
         Pnm = (0.0_MK,0.0_MK)
         
         isub = ppm_isublist(i,in_topoid)
         box = ppm_boxid(isub,topoid)
         first = lhbx(1,box)
         last  = lhbx(2,box)
         !----------------------------------------------------------------------
         ! Computing new array with particle order from tree
         !----------------------------------------------------------------------

         IF (ppm_dim .GT. 2) THEN
             DO j=first,last
                 xp(1,j) = xpunord(1,lpdx(j))
                 xp(2,j) = xpunord(2,lpdx(j))
                 xp(3,j) = xpunord(3,lpdx(j))
             ENDDO
         ELSE
             DO j=first,last
                 xp(1,j) = xpunord(1,lpdx(j))
                 xp(2,j) = xpunord(2,lpdx(j))
             ENDDO
         ENDIF

         nrpart = last-first + 1
         CALL ppm_util_cart2sph(xp(1,first:last), &
              & xp(2,first:last),xp(3,first:last),nrpart, &
              & centerofbox(1,box),centerofbox(2,box),centerofbox(3,box), &
              & rho(first:last),theta(first:last),phi(first:last),info)
         IF (info .NE. 0) THEN
            CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_expansion', &
                 & 'Failed calling util_cart2sph',__LINE__,info)
         ENDIF
         DO j=first,last !loop over particles in leaf
            cosine = COS(theta(j))
            sine   = SIN(theta(j))
            !-------------------------------------------------------------------
            !  Recurrence for Pnm
            !-------------------------------------------------------------------
            val  = -sine
            prod = 1.0_MK
            DO m=0,order
               Pnm(m,m) = fracfac(m)*prod
               prod     = prod * val
            ENDDO
            DO m=0,order-1
               Pnm(m+1,m) = cosine*REAL(2*m + 1,MK)*Pnm(m,m)
            ENDDO
            DO n=2,order
               val = cosine*REAL(2*n-1,MK)
               DO m=0,n-2
                  Pnm(n,m)=(val*Pnm(n-1,m) - & 
                  &        DBLE(n+m-1)*Pnm(n-2,m))/REAL(n-m,MK)
               ENDDO
            ENDDO
            !-------------------------------------------------------------------
            !  Compute Ynm(n,m) and Ynm(n,-m)
            !-------------------------------------------------------------------
            DO n=0,order
               m = 0
               angle    = REAL(m,MK)*phi(j)
               Ynm(n,m) = sqrtfac(n,m)*Pnm(n,m) &
               &           *CMPLX(COS(angle),SIN(angle))
               DO m=1,n
                  angle     = REAL(m,MK)*phi(j)
                  Ynm(n,m)  = sqrtfac(n,m) &
                  &           *Pnm(n,m)*CMPLX(COS(angle),SIN(angle))
                  Ynm(n,-m) = CONJG(Ynm(n,m))
               ENDDO
            ENDDO
            !-------------------------------------------------------------------
            !  Compute Cnm(n,m) - the expansion coefficients
            !-------------------------------------------------------------------
            prod = 1.0_MK
            val  = rho(j)
            DO n=0,order
               DO m=0,n
                  Cnm(n,m) = Cnm(n,m) + (wpunord(lpdx(j))*prod*Ynm(n,m)* &
                  &          Anm(n,m))/((-1)**n)
               ENDDO
               prod = prod * val
            ENDDO
         ENDDO !particles in one sub
         DO n=0,order
            DO m=1,n
               Cnm(n,-m) = CONJG(Cnm(n,m))
            ENDDO
         ENDDO
         DO m=1,order
            temp = CI**(-m)
            DO n=m,order
               Cnm(n,m) = Cnm(n,m)*temp
               Cnm(n,-m)= Cnm(n,-m)*temp
            ENDDO
         ENDDO
         ! saving expansion in fmm_module_data file
         expansion(box,0:order,-order:order) = Cnm(0:order,-order:order)
      ENDDO !loop over all subs

      !-------------------------------------------------------------------------
      !  Nullify data pointers
      !-------------------------------------------------------------------------

      NULLIFY(min_box)
      NULLIFY(max_box)
      NULLIFY(boxcost)
      NULLIFY(centerofbox)
      NULLIFY(radius)
      NULLIFY(expansion)
      NULLIFY(Anm)
      NULLIFY(sqrtfac)
      NULLIFY(fracfac)
      NULLIFY(Ynm)
      NULLIFY(Pnm)
      NULLIFY(Cnm)

      !-------------------------------------------------------------------------
      ! Traversing the tree and shifting the expansions upwards
      !-------------------------------------------------------------------------

      IF (parent(1) .EQ. ppm_param_undefined) THEN
         root = 1
      ELSE
         DO i=1,nbox
            IF (parent(i) .EQ. ppm_param_undefined) THEN
               root = i
               EXIT
            ENDIF
         ENDDO
      ENDIF

      CALL ppm_fmm_traverse(root,order,t0,info)
      IF (info.NE.0) THEN
         CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_expansion', &
              &         'traversing tree failed',__LINE__,info)
      ENDIF
      IF (ppm_debug .GT. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_fmm_expansion','traversed tree',info)
      ENDIF     

      !-------------------------------------------------------------------------
      ! deallocate local data
      !-------------------------------------------------------------------------
      ldu2 = (0)
      CALL ppm_alloc(xp,ldu2,ppm_param_dealloc,info)
      IF (info .NE. 0) THEN
          WRITE(cbuf,'(A,I3,A)') 'for ',info,'error while dealloc'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expansion',cbuf,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_fmm_expansion',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fmm_expansion_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fmm_expansion_d
#endif
