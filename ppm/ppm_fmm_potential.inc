      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_fmm_potential
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Compute the potential of the target points using
      !                 the expansion coefficients
      !                 In the parallel case:
      !                 Calls ppm_fmm_pretraverse and ppm_fmm_expchange
      !                 Maps the target points onto the leaf topolgy
      !                 
      !  Input        : xp(:,:)      (F) the field points
      !                 wp(:)        (F) the strength of the field points
      !                 Np           (I) the number of field points.
      !                 tp(:,:)      (F) the target points
      !                 Ntp          (I) the number of target points
      !                 order        (I) (maximum) expansion order
      !                 theta        (F) acceptance factor
      !
      !  Input/output :     
      !
      !  Output       : potential(:) (F) the multipole expansion potential
      !                                  for each point
      !                                  size 1:Np
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      :  
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_potential.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.25  2005/09/19 13:03:29  polasekb
      !  code cosmetics
      !
      !  Revision 1.24  2005/09/12 09:14:10  hiebers
      !  added mapping of target points
      !
      !  Revision 1.23  2005/09/11 18:05:31  polasekb
      !  (final?) corrected version
      !  (also works parallel :-)
      !
      !  Revision 1.22  2005/09/10 07:50:04  polasekb
      !  changed init of stack for parallel version
      !
      !  Revision 1.21  2005/09/05 06:23:57  polasekb
      !  corrected variable initialisation
      !
      !  Revision 1.20  2005/08/29 15:18:00  polasekb
      !  bugfix when computing direct way
      !
      !  Revision 1.19  2005/08/25 13:52:10  polasekb
      !  mapping corrected from wp to wpunord
      !
      !  Revision 1.18  2005/08/23 14:35:20  polasekb
      !  changed call to ppm_fmm_pretraverse
      !
      !  Revision 1.17  2005/08/23 14:30:05  polasekb
      !  changed call to ppm_fmm_expchange
      !
      !  Revision 1.16  2005/08/23 14:16:50  polasekb
      !  corrected wpunord and wp
      !
      !  Revision 1.15  2005/08/11 15:12:04  polasekb
      !  fixed indices of Outer
      !
      !  Revision 1.14  2005/08/11 13:32:34  polasekb
      !  added variable theta
      !
      !  Revision 1.13  2005/08/08 13:35:52  polasekb
      !  deallocate some local variables
      !
      !  Revision 1.12  2005/08/04 16:04:51  polasekb
      !  removed some data allocation
      !
      !  Revision 1.11  2005/07/29 14:06:52  polasekb
      !  changed check of eqalness to ppm_myeps
      !
      !  Revision 1.10  2005/07/29 12:36:13  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.9  2005/07/27 14:59:24  polasekb
      !  now using constants from data file
      !
      !  Revision 1.8  2005/07/25 15:01:20  polasekb
      !  bugfix with indices
      !
      !  Revision 1.7  2005/07/25 14:40:31  polasekb
      !  adapted computation of the potential
      !
      !  Revision 1.6  2005/07/21 12:42:07  polasekb
      !  bugfix in allocating an array
      !
      !  Revision 1.5  2005/07/21 08:26:09  polasekb
      !  changed function call, now different target 
      !  points and field points can be
      !  specified by the user
      !
      !  Revision 1.4  2005/06/02 19:18:51  polasekb
      !  corrected syntax error
      !
      !  Revision 1.3  2005/05/30 09:37:01  polasekb
      !  bugfix: corrected call to ppm_util_cart2sph
      !
      !  Revision 1.2  2005/05/27 12:42:48  polasekb
      !  initialized further arrays
      !
      !  Revision 1.1  2005/05/27 08:01:23  polasekb
      !  initial implementation
      !  TODO: remove debug output
      !
      !  Revision 0  2005/01/16 15:59:14 polasekb
      !  start
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fmm_potential_s(xpunord,wpunord,Np,tp,Ntp,order,theta, &
                 &                   potential,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fmm_potential_d(xpunord,wpunord,Np,tp,Ntp,order,theta, &
                 &                   potential,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_fmm_expchange
      USE ppm_module_fmm_pretraverse
      USE ppm_module_map
      USE ppm_module_map_part_get_sub
      USE ppm_module_mktopo
      USE ppm_module_substart
      USE ppm_module_substop 
      USE ppm_module_topo
      USE ppm_module_topo_box2subs
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
      REAL(MK), DIMENSION(:,:), POINTER       :: xpunord
      REAL(MK), DIMENSION(:  ), POINTER       :: wpunord
      INTEGER                 , INTENT(INOUT) :: Np
      REAL(MK), DIMENSION(:,:), POINTER       :: tp
      INTEGER                 , INTENT(INOUT) :: Ntp
      INTEGER                 , INTENT(IN   ) :: order
      REAL(MK)                , INTENT(IN   ) :: theta
#if   __KIND == __SINGLE_PRECISION
      REAL(ppm_kind_single),DIMENSION(:),POINTER :: potential
#elif __KIND == __DOUBLE_PRECISION
      REAL(ppm_kind_double),DIMENSION(:),POINTER :: potential
#endif
      INTEGER                 , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------

      LOGICAL                              :: check,drct,OK
      INTEGER                              :: i,j,k,cnt,iopt,m,n
      INTEGER                              :: pcount,ccount
      INTEGER                              :: mapt,Mpart,root,istat
      INTEGER                              :: first,last,level
      INTEGER                              :: stackpointer,curbox
      INTEGER                              :: pexp,topoid,isymm
      INTEGER ,DIMENSION(1)                :: ldu1
      INTEGER ,DIMENSION(2)                :: ldu2 
      INTEGER ,DIMENSION(:  ), POINTER     :: part_subtop => null()
      INTEGER ,DIMENSION(:  ), POINTER     :: newlpdx => null()
      INTEGER ,DIMENSION(:  ), POINTER     :: stack => null()
      INTEGER ,DIMENSION(:  ), POINTER     :: wp => null()
      INTEGER ,DIMENSION(:,:), POINTER     :: coeff_subtop => null()
      REAL(MK)                             :: t0,ghostsize,cutoff
      REAL(MK)                             :: dx,dy,dz,dist,rad
      REAL(MK)                             :: thetap,eps,angle,reci 
      REAL(MK)                             :: sine,cosine,val,prod 
      REAL(MK),DIMENSION(1)                :: curboxrho,curboxphi,curboxtheta
      REAL(MK),DIMENSION(:  ),     POINTER :: fracfac => null()
      REAL(MK),DIMENSION(:  ),     POINTER :: radius => null()
      REAL(MK),DIMENSION(:  ),     POINTER :: boxcost => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: sqrtfac => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: xp => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: Pnm => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: Anm => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: centerofbox => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: min_box => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: max_box => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: min_sub => null()
      REAL(MK),DIMENSION(:,:),     POINTER :: max_sub => null()
      COMPLEX(MK),DIMENSION(:,:),  POINTER :: Ynm => null()
      COMPLEX(MK),DIMENSION(:,:),  POINTER :: Outer => null()
      COMPLEX(MK),DIMENSION(:,:,:),POINTER :: expansion => null()
      COMPLEX(MK),PARAMETER                :: CI=(0.0_MK,1.0_MK)
      CHARACTER(LEN=ppm_char)              :: cbuf

      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------

      CALL substart('ppm_fmm_potential',t0,info)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN  
         DO i=1,ppm_dim
            IF (order .EQ. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_potential',   &
      &               'expansion order must be > 0 !',__LINE__,info)
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
      sqrtfac      => sqrtfac_s
      fracfac      => fracfac_s
      Anm          => Anm_s
      eps          = ppm_myepss
      Ynm          => Ynm_s
      Pnm          => Pnm_s
      Outer        => Cnm_s
#else
      min_box      => min_box_d
      max_box      => max_box_d
      boxcost      => boxcost_d
      centerofbox  => centerofbox_d
      radius       => radius_d
      expansion    => expansion_d
      sqrtfac      => sqrtfac_d
      fracfac      => fracfac_d
      Anm          => Anm_d
      eps          = ppm_myepsd
      Ynm          => Ynm_d
      Pnm          => Pnm_d
      Outer        => Cnm_d
#endif

      IF (ppm_nproc .GT. 1) THEN
      !choose lowest level of the tree as topology
      topoid = nlevel
      
      !-------------------------------------------------------------------------
      ! Calling fmm_expchange to communicate all expansions
      !-------------------------------------------------------------------------
      
      CALL ppm_fmm_expchange(order,t0,info)
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to call expchange.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! Mapping target points
      !-------------------------------------------------------------------------
      mapt = ppm_param_map_global
      CALL ppm_map_part(tp,ppm_dim,Ntp,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF
      mapt = ppm_param_map_send
      CALL ppm_map_part(tp,ppm_dim,Ntp,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF
      mapt = ppm_param_map_pop
      CALL ppm_map_part(tp,ppm_dim,Ntp,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF

      Ntp = Mpart
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &          'Done mapping target points.',info)
          WRITE(cbuf,'(A,I6)') 'Local number of target points now:',Ntp
          CALL ppm_write(ppm_rank,'ppm_fmm_potential',cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that particles have been mapped correctly
      !-------------------------------------------------------------------------
      ! If processor has no more particles, end
      IF (Np .EQ. 0) THEN
         GOTO 9999
      ENDIF

      CALL ppm_topo_check(tp,Ntp,OK,info)
      IF (info .NE. 0) THEN
         CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
         &    'Failed to check topology.',info)
      ENDIF
      IF (.NOT.OK) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'target points not mapped correctly!',info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      ! Calling pretraversal routine to build communication lists
      !-------------------------------------------------------------------------

      CALL ppm_fmm_pretraverse(tp,Ntp,topoid,theta,coeff_subtop,ccount,&
      &                       part_subtop,pcount,info)
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to call pretraverse.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! Calling map_ghost routine to get particles on local processor
      !-------------------------------------------------------------------------
      
      CALL ppm_map_part_get_sub(part_subtop,pcount,topoid,xpunord,Np,info)
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to call map_part_get_sub.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! Mapping particles, wheights and boxpart as ghosts
      !-------------------------------------------------------------------------
      isymm  = 0
      cutoff = 1.0_MK

      mapt = ppm_param_map_push
      CALL ppm_map_part_ghost(wpunord,Np,Mpart,isymm,cutoff,mapt,info) !strengths
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      CALL ppm_map_part_ghost(boxpart,Np,Mpart,isymm,cutoff,mapt,info) !boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF      
      mapt = ppm_param_map_send
      CALL ppm_map_part_ghost(boxpart,Np,Mpart,isymm,cutoff,mapt,info)   ! send
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to send particles.',info)
          GOTO 9999
      ENDIF
      mapt = ppm_param_map_pop
      CALL ppm_map_part_ghost(boxpart,Np,Mpart,isymm,cutoff,mapt,info)  !boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF      
      CALL ppm_map_part_ghost(wpunord,Np,Mpart,isymm,cutoff,mapt,info) !strengths
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to pop strengths.',info)
          GOTO 9999
      ENDIF !positions
      CALL ppm_map_part_ghost(xpunord,ppm_dim,Np,Mpart,isymm,cutoff,mapt,info)
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to pop positions.',info)
          GOTO 9999
      ENDIF

      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &          'Done mapping ghost particles.',info)
          WRITE(cbuf,'(A,I6)') 'Received ghost particles:',Mpart-Np
          CALL ppm_write(ppm_rank,'ppm_fmm_potential',cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Sorting new (ghost) particles into array
      !-------------------------------------------------------------------------

      ldu1(1) = Mpart
      CALL ppm_alloc(newlpdx,ldu1,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_potential', &
      &       'error allocating newlpdx',__LINE__,info)
      GOTO 9999
      ENDIF

      newlpdx(1:Np) = lpdx(1:Np)

      cnt = Np
      DO i=1,nbox
         check = .TRUE.
         DO j=Np+1,Mpart
            IF (boxpart(j) .EQ. i) THEN
               cnt = cnt + 1
               IF (((lhbx(1,i) .EQ. 1) .OR. (lhbx(1,i) .EQ. Np+1)) .AND. &
                    (check)) THEN
                  check = .FALSE.
                  lhbx(1,i) = cnt
                  lhbx(2,i) = cnt
                  newlpdx(cnt) = j
               ELSE
                  lhbx(2,i) = lhbx(2,i) + 1
                  newlpdx(cnt) = j
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      !-------------------------------------------------------------------------
      ! ATTENTION: now tree lists lhbx and lpdx are not valid anymore
      ! can be extended later
      !-------------------------------------------------------------------------
      ENDIF !ppm_nproc .GT. 1

      !-------------------------------------------------------------------------
      ! Allocating array for potentials of target points
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit

      ldu1(1) = Ntp
      CALL ppm_alloc(potential,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_potential', &
      &       'error allocating potential',__LINE__,info)
      GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! Initialize arrays
      !-------------------------------------------------------------------------
      potential = (0.0_MK)

      !-------------------------------------------------------------------------
      ! allocating and initializing further arrays
      !-------------------------------------------------------------------------

      istat = 0
      ldu1(1) = nbox
      CALL ppm_alloc(stack,ldu1,iopt,info)
      istat = istat + info
 
      ldu2(1) = 3
      IF (ppm_nproc .GT. 1) THEN
         ldu2(2) = Mpart
      ELSE
         ldu2(2) = Np
      ENDIF
      CALL ppm_alloc(xp,ldu2,iopt,info)
      istat = istat + info
      IF (ppm_nproc .GT. 1) THEN
         ldu1(1) = Mpart
      ELSE
         ldu1(1) = Np
      ENDIF
      CALL ppm_alloc(wp,ldu1,iopt,info)
      istat = istat + info
      IF (istat .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_potential', &
      &       'error allocating variables',__LINE__,info)
      GOTO 9999
      ENDIF

      stack = (0)
      wp    = (0.0_MK)
      xp    = (0.0_MK,0.0_MK)

      !-------------------------------------------------------------------------
      ! Finding the root of the tree (serial)
      ! Finding top level of tree (parallel)
      !-------------------------------------------------------------------------

      IF (ppm_nproc .GT. 1) THEN        
         ! finding top level
         DO i=1,nlevel
           IF (nbpl(i) .GE. ppm_nproc) THEN
              level = i
              EXIT
           ENDIF
         ENDDO
      ELSE
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
      ENDIF

      !-------------------------------------------------------------------------
      ! ordering particles according to tree
      !-------------------------------------------------------------------------
      IF (ppm_nproc .GT. 1) THEN
          DO i=1,nbox
             first = lhbx(1,i)
             last  = lhbx(2,i)
             DO j=first,last
                xp(:,j) = xpunord(:,newlpdx(j))
                wp(j)   = wpunord(newlpdx(j))
             ENDDO
          ENDDO
      ELSE
          DO i=1,nbox
             first = lhbx(1,i)
             last  = lhbx(2,i)
             DO j=first,last
                xp(:,j) = xpunord(:,lpdx(j))
                wp(j)   = wpunord(lpdx(j))
             ENDDO
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      ! Computing the potential for the target points
      !-------------------------------------------------------------------------
      DO i=1,Ntp
         IF (ppm_nproc .GT. 1) THEN        
!          init stack parallel
           stackpointer = 1
           cnt = 0
           DO j=1,nbox
             IF (blevel(j) .EQ. level) THEN
               stack(stackpointer) = j 
               stackpointer = stackpointer + 1
               cnt = cnt +1
               IF (cnt .EQ. nbpl(level)) THEN
                  EXIT
               ENDIF
             ENDIF
           ENDDO
         ELSE
!         init stack serial
          stackpointer = 1
          stack(stackpointer) = root
          stackpointer = stackpointer + 1
         ENDIF
         DO WHILE (stackpointer .GT. 1)
           curbox = stack(stackpointer-1)
           !pop top box
           stackpointer = stackpointer - 1
           dx = tp(1,i) - centerofbox(1,curbox)
           dy = tp(2,i) - centerofbox(2,curbox)
           dz = tp(3,i) - centerofbox(3,curbox)
           dist = SQRT(dx*dx + dy*dy + dz*dz)
           !--------------------------------------------------------------------
           ! Checking Barnes-Hut Criterium
           !--------------------------------------------------------------------
           drct = .FALSE.
           IF (radius(curbox) .EQ. 0) THEN
              !only one particle in box, do direct computation
              drct = .TRUE.
           ENDIF
           IF ((dist/(2*radius(curbox)) .GT. theta) .AND. (.NOT. drct)) THEN
             !-----------------------------------------------------------------
             !  far enough, compute part-box interaction
             !-----------------------------------------------------------------
             IF (ppm_debug .GT. 0) THEN
                  CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                  &    'far enough',info)
             ENDIF
             !------------------------------------------------------------------
             ! TODO:
             ! Computing the expansion order pexp according to Wang
             !------------------------------------------------------------------
             !thetap = 0.75_MK+0.2_MK*(order-pexp)+0.05_MK*(order-pexp)**2
             pexp = order
             !DO WHILE((thetap .LE. dist/radius(curbox)).AND.(pexp .GT. 3))
             !   pexp = pexp - 1
             !   thetap = 0.75_MK+0.2_MK*(order-pexp)+0.05_MK*(order-pexp)**2
             !ENDDO

             CALL ppm_util_cart2sph(tp(1,i:i),tp(2,i:i),tp(3,i:i),1, &
                 & centerofbox(1,curbox),centerofbox(2,curbox), &
                 & centerofbox(3,curbox), &
                 & curboxrho,curboxtheta,curboxphi,info)
             IF (info .NE. 0) THEN
                CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_potential', &
                     & 'Failed calling util_cart2sph',__LINE__,info)
             ENDIF
             !compute expansion
             !------------------------------------------------------------------
             !      Compute Legendre polynomial box-particle interaction
             !------------------------------------------------------------------
             ! init variables
             Ynm       = (0.0_MK,0.0_MK)
             Outer     = (0.0_MK,0.0_MK)
             Pnm       = (0.0_MK,0.0_MK)

             reci      = 1.0_MK/curboxrho(1)
             sine      = SIN(curboxtheta(1))
             cosine    = COS(curboxtheta(1))
             val       = -sine
             prod      = 1.0_MK
             DO m=0,pexp
               Pnm(m,m) = fracfac(m)*prod
               prod     = prod * val
             ENDDO
             DO m=0,pexp-1
                Pnm(m+1,m) = cosine*REAL(2*m + 1,MK)*Pnm(m,m)
             ENDDO
             DO n=2,pexp
                val = cosine*REAL(2*n-1,MK)
                DO m=0,n-1
                   Pnm(n,m)=(val*Pnm(n-1,m)-REAL(n+m-1,MK)* &
                                      Pnm(n-2,m))/REAL(n-m,MK)
                ENDDO
             ENDDO
             !------------------------------------------------------------------
             !      Compute Ynm(n,m) and Ynm(n,-m)
             !------------------------------------------------------------------
             DO n=0,pexp
                m = 0
                angle = REAL(m,MK)*curboxphi(1)
                Ynm(n,m) = sqrtfac(n,m)*Pnm(n,m)*CMPLX(COS(angle),SIN(angle))
                DO m=1,n
                   angle     = REAL(m,MK)*curboxphi(1)
                   Ynm(n,m)  = sqrtfac(n,m)*Pnm(n,m)*   &
                               CMPLX(COS(angle),SIN(angle))
                   Ynm(n,-m) = CONJG(Ynm(n,m))
                ENDDO
             ENDDO

             !------------------------------------------------------------------
             !    Compute the Outer expansion
             !------------------------------------------------------------------
             prod = 1.0_MK
             DO n=0,pexp
                prod = prod * curboxrho(1)
                DO m=-n,n
                   Outer(n,m) = (-1)**n*CI**ABS(m)*Ynm(n,m)/(Anm(n,m)* & 
                   &            prod)
                ENDDO
             ENDDO

             !------------------------------------------------------------------
             !    Evaluate potential, using multipole expansion coefficients
             !------------------------------------------------------------------
             DO n=0,pexp
                DO m=-n,n
                   potential(i)=potential(i) + expansion(curbox,n,m) &
                    &           *Outer(n,-m)
                ENDDO
             ENDDO

           ELSE
             !-----------------------------------------------------------------
             !  not far enough, push childern
             !-----------------------------------------------------------------
             IF (ppm_debug .GT. 0) THEN
                CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                &    'not far enough',info)
             ENDIF
             IF (nchld(curbox) .GT. 0) THEN
                IF (ppm_debug .GT. 0) THEN
                   CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                   &    'push children',info)
                ENDIF
               DO j=1,nchld(curbox)
                  stack(stackpointer) = child(j,curbox)
                  stackpointer = stackpointer + 1
               ENDDO
             ELSE
             !------------------------------------------------------------------
             !  no children, direct computation
             !------------------------------------------------------------------
              IF (ppm_debug .GT. 0) THEN
                 CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                 &    'no children',info)
              ENDIF
              first = lhbx(1,curbox)
              last  = lhbx(2,curbox)
              DO j=first,last !loop over particles in leaf
                  !-------------------------------------------------------------
                  !      Evaluate potential, direct method
                  !-------------------------------------------------------------
                  !IF (i .EQ. j) THEN
                  IF ((ABS(tp(1,i) - xp(1,j)) .LE. eps) .AND. &
                  &   (ABS(tp(2,i) - xp(2,j)) .LE. eps) .AND. &
                  &   (ABS(tp(3,i) - xp(3,j)) .LE. eps)) THEN
                  ELSE
   	            dx = xp(1,j) - tp(1,i)
	            dy = xp(2,j) - tp(2,i)
	            dz = xp(3,j) - tp(3,i)
	            rad = SQRT(dx*dx + dy*dy + dz*dz)                  
                    potential(i) = potential(i) + wp(j)/rad
                  ENDIF
              ENDDO
             ENDIF
           ENDIF
         ENDDO 
      ENDDO

      !-------------------------------------------------------------------------
      !  Nullify data pointers
      !-------------------------------------------------------------------------

      NULLIFY(min_box)
      NULLIFY(max_box)
      NULLIFY(boxcost)
      NULLIFY(centerofbox)
      NULLIFY(radius)

      !-------------------------------------------------------------------------
      !  Deallocate local data
      !-------------------------------------------------------------------------

      istat     = 0
      ldu1(1)   = 0
      ldu2(1:2) = 0
      CALL ppm_alloc(newlpdx,ldu1,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(stack,ldu1,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(xp,ldu2,ppm_param_dealloc,info)
      istat = istat + info
      IF (istat .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expansion', &
      &       'error deallocating newlpdx',__LINE__,info)
      GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_fmm_potential',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fmm_potential_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fmm_potential_d
#endif
