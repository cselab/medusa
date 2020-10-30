      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_fmm_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Initialisation of FMM. Calls the tree to get the tree 
      !                 information and saves this in the ppm_module_data_fmm 
      !                 file. 
      !                 Maps the particles to the leafs of the tree.
      !                 Computes then the center of the boxes and the radius
      !                 of the leaf boxes and saves this in the same file.
      !
      !  Input        : xp(:,:)      (F) the field points
      !                 wp(:)        (F) field particle strenghts
      !                 Np           (I) the number of field points.
      !                 Nm(:)        (I) number of grid points in the
      !                                  global mesh. (0,0,0) if there is
      !                                  no mesh. If a mesh is present, the
      !                                  box boundaries will be aligned
      !                                  with mesh planes.
      !                 order        (I) expansion order
      !                 min_dom(:)   (F) the minimum coordinate of the
      !                                  domain
      !                 max_dom(:)   (F) the maximum coordinate of the
      !                                  domain
      !                 maxboxcost   (F) the maximum number of particles
      !                                  allowed in a box
      !  Input/output :     
      !
      !  Output       : nrofbox     (I) the total number of boxes
      !                 info        (I) return status. 0 upon success.
      !
      !  Remarks      :  
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_init.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.17  2005/09/19 13:03:28  polasekb
      !  code cosmetics
      !
      !  Revision 1.16  2005/09/12 13:30:33  polasekb
      !  added ppm_subid
      !
      !  Revision 1.15  2005/09/11 18:05:30  polasekb
      !  (final?) corrected version
      !  (also works parallel :-)
      !
      !  Revision 1.14  2005/09/11 11:43:39  polasekb
      !  moved mapping and second tree call to init
      !
      !  Revision 1.13  2005/08/30 08:48:30  polasekb
      !  added timing for tree
      !
      !  Revision 1.12  2005/08/25 13:51:49  polasekb
      !  corrected data allocation of theta,phi,rho
      !
      !  Revision 1.11  2005/08/11 15:12:53  polasekb
      !  added argument maxboxcost
      !
      !  Revision 1.10  2005/08/08 13:34:25  polasekb
      !  removec fmm_prec
      !
      !  Revision 1.9  2005/08/04 16:00:41  polasekb
      !  moved some allocation to init
      !
      !  Revision 1.8  2005/07/29 12:35:05  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.7  2005/07/27 14:58:26  polasekb
      !  added new argument wp to subroutine call
      !
      !  Revision 1.6  2005/07/25 15:01:57  polasekb
      !  adapted some tree coefficients
      !
      !  Revision 1.5  2005/07/25 13:39:20  polasekb
      !  bugfix in array indices
      !
      !  Revision 1.4  2005/07/21 13:21:32  polasekb
      !  removed nullify
      !
      !  Revision 1.3  2005/06/02 13:54:55  polasekb
      !  removed totalmass
      !
      !  Revision 1.2  2005/05/30 09:37:24  polasekb
      !  correctet computing of centerofbox
      !
      !  Revision 1.1  2005/05/27 07:53:40  polasekb
      !  Initial Implementation
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

#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fmm_init_s(xp,wp,Np,Nm,order,min_dom,max_dom,maxboxcost, &
      &          nrofbox,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fmm_init_d(xp,wp,Np,Nm,order,min_dom,max_dom,maxboxcost, &
      &          nrofbox,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_tree
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_map
      USE ppm_module_topo_box2subs
      USE ppm_module_topo      
      USE ppm_module_substart
      USE ppm_module_substop 
      USE ppm_module_write

      IMPLICIT NONE
      
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      
      
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), POINTER       :: xp
      REAL(MK), DIMENSION(:  ), POINTER       :: wp
      INTEGER                 , INTENT(INOUT) :: Np
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      INTEGER                 , INTENT(IN   ) :: order
      REAL(MK)                , INTENT(IN   ) :: maxboxcost
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_dom,max_dom
      INTEGER                 , INTENT(  OUT) :: nrofbox,info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------

      LOGICAL                             :: pruneboxes,OK
      LOGICAL,DIMENSION(3)                :: fixed
      INTEGER,DIMENSION(1)                :: ldu1
      INTEGER,DIMENSION(2)                :: ldu2
      INTEGER,DIMENSION(3)                :: ldu,ldl
      INTEGER,DIMENSION(6)                :: bcdef      
      INTEGER,DIMENSION(:),    POINTER    :: box2proc,boxid 
      INTEGER,DIMENSION(:),    POINTER    :: subs2proc,isublist
      INTEGER,DIMENSION(:),    POINTER    :: new_subs2proc
      INTEGER                             :: treetype,minboxes
      INTEGER                             :: maxlevels,level
      INTEGER                             :: iopt,i,j
      INTEGER                             :: mapt,Mpart
      INTEGER                             :: box,first,last,nrpbox
      INTEGER                             :: nsubs,topoid,istat
      INTEGER                             :: decomp,assig,nsublist   
      REAL(MK)                            :: maxvariance
      REAL(MK)                            :: t0, ghostsize,tmp
      REAL(MK),DIMENSION(3  )             :: minboxsize
      REAL(MK),DIMENSION(3  )             :: diagvec
      REAL(MK),DIMENSION(3,2)             :: weights
      REAL(MK),DIMENSION(:),   POINTER    :: treewp,cost,boxcost 
      REAL(MK),DIMENSION(:),   POINTER    :: radius,totalmass
      REAL(MK),DIMENSION(:,:), POINTER    :: centerofbox,treepart
      REAL(MK),DIMENSION(:,:), POINTER    :: min_box,max_box
      REAL(MK),DIMENSION(:,:), POINTER    :: min_sub,max_sub      
      CHARACTER(LEN=ppm_char)             :: cbuf

      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------

      CALL substart('ppm_fmm_init',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN  
         DO i=1,ppm_dim
            IF (min_dom(i) .GT. max_dom(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_init',   &
      &               'min_dom must be <= max_dom !',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      ! Setting tree variables 
      !-------------------------------------------------------------------------
      treetype     = ppm_param_tree_oct
      minboxes     = ppm_nproc
      pruneboxes   = .FALSE.
      minboxsize   = (/0.001_MK,0.001_MK,0.001_MK/)
      maxvariance  = -1.0_MK
      !maxboxcost   = input argument
      fixed        = (/.FALSE.,.FALSE.,.FALSE./)
      weights(:,1) = (/1.0_MK,0.0_MK,0.0_MK/)
      weights(:,2) = (/0.0_MK,0.0_MK,1.0_MK/)
      maxlevels    = 20

#if   __KIND == __SINGLE_PRECISION
      !-------------------------------------------------------------------------
      ! Calling the tree single (1)
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,             &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_s,max_box_s,lhbx,lpdx,boxcost_s,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_s
      max_box => max_box_s
      boxcost => boxcost_s
#else
      !-------------------------------------------------------------------------
      ! Calling the tree double (2)
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,             &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_d,max_box_d,lhbx,lpdx,boxcost_d,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_d
      max_box => max_box_d
      boxcost => boxcost_d
#endif

      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init',  &
       &       'Calling tree failed.',__LINE__,info)
          GOTO 9999
      ENDIF
      nrofbox = nbox
!     output information!
      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init', &
         &    'calling tree successful',info)
         WRITE (cbuf,'(A,I)') 'nbox = ',nbox
         CALL ppm_write(ppm_rank,'ppm_fmm_init',cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Allocating boxpart
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit
      ldu1(1) = Np
      CALL ppm_alloc(boxpart,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating boxpart',__LINE__,info)
      GOTO 9999
      ENDIF

      boxpart = (0)

      !-------------------------------------------------------------------------
      ! Storing which particle is in which leaf-box
      !-------------------------------------------------------------------------
      IF (Np .GT. 0) THEN
        DO i=1,nbox
          IF (nchld(i) .EQ. 0) THEN
             DO j=lhbx(1,i),lhbx(2,i)
                boxpart(lpdx(j)) = i
             ENDDO
           ENDIF
        ENDDO
      ELSE
        boxpart = (0)
      ENDIF

      !-------------------------------------------------------------------------
      ! make top level topology seperate
      !-------------------------------------------------------------------------

      DO i=1,nlevel
        IF (nbpl(i) .GE. ppm_nproc) THEN
           level = i
           EXIT
        ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      ! Allocating ppm_boxid, box2proc and cost
      !-------------------------------------------------------------------------

      ldu2(1) = nbox
      ldu2(2) = nlevel
      CALL ppm_alloc(ppm_boxid,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating ppm_boxid',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_subid,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating ppm_boxid',__LINE__,info)
      GOTO 9999
      ENDIF

      ldu1(1) = nbox
      CALL ppm_alloc(box2proc,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating box2proc',__LINE__,info)
      GOTO 9999
      ENDIF

      ldu1(1) = MAXVAL(nbpl(:))
      CALL ppm_alloc(cost,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating box2proc',__LINE__,info)
      GOTO 9999
      ENDIF

      ppm_boxid = (0,0)
      ppm_subid = (0,0)
      box2proc  = (0)
      cost      = (1.0_MK)
      decomp     = ppm_param_decomp_user_defined
      assig      = ppm_param_assign_internal
      ghostsize  = 0.0_MK
      bcdef(1:6) = ppm_param_bcdef_freespace
      CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,max_sub, &
           &                 nsubs,info,boxid,-level,blevel,child)
      IF (info.NE.0) THEN
         CALL ppm_write(ppm_rank,'ppm_fmm_init', &
         &    'topo_box2subs failed',info)
      ENDIF

      CALL ppm_mktopo(decomp,assig,min_dom,max_dom,bcdef,ghostsize, &
           &          level,min_sub,max_sub,cost,subs2proc,nsubs,isublist, &
           &          nsublist,info)
      IF (info.NE.0) THEN
         CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init', &
              &         'mktopo failed',__LINE__,info)
      ENDIF
      ppm_boxid(1:nsubs,level) = boxid(1:nsubs)
      DO j=1,SIZE(boxid)
         ppm_subid(boxid(j),level) = j
      ENDDO

      DO i=1,nsubs
         box2proc(boxid(i)) = subs2proc(i)
      ENDDO

      !-------------------------------------------------------------------------
      ! Looping over levels of tree and registering as topolgy
      !-------------------------------------------------------------------------

      ! topoid of level i = i

      assig      = ppm_param_assign_user_defined
      DO i=level+1,nlevel

        topoid = i

        !-----------------------------------------------------------------------
        ! Calling subroutine to get subs
        !-----------------------------------------------------------------------

        CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,max_sub, &
          &                      nsubs,info,boxid,-topoid,blevel,child)
        IF (info.NE.0) THEN
           CALL ppm_write(ppm_rank,'ppm_fmm_init', &
           &    'topo_box2subs failed',info)
        ENDIF

        !-----------------------------------------------------------------------
        ! Allocating new subs2proc
        !-----------------------------------------------------------------------

        iopt = ppm_param_alloc_grow
        ldu1(1) = nsubs
        CALL ppm_alloc(new_subs2proc,ldu1,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
        &       'error allocating new_subs2proc',__LINE__,info)
        GOTO 9999
        ENDIF

        DO j=1,nsubs
           new_subs2proc(j) = box2proc(parent(boxid(j)))
           box2proc(boxid(j)) = new_subs2proc(j)
        ENDDO

        !-----------------------------------------------------------------------
        ! Calling subroutine to get topology
        !-----------------------------------------------------------------------
        CALL ppm_mktopo(decomp,assig,min_dom,max_dom,bcdef,ghostsize, & 
        &          topoid,min_sub,max_sub,cost,new_subs2proc,nsubs, &  
        &          isublist,nsublist,info)
        IF (info.NE.0) THEN
           CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init', &
           &         'mktopo failed',__LINE__,info)
        ENDIF
        ppm_boxid(1:nsubs,topoid) = boxid(1:nsubs)
        DO j=1,SIZE(boxid)
           ppm_subid(boxid(j),topoid)   = j
        ENDDO
      ENDDO 

      !-------------------------------------------------------------------------
      ! Mapping for the lowest level (leafs of tree)
      !-------------------------------------------------------------------------

      topoid = nlevel

      IF (ppm_nproc .GT. 1) THEN

      !-------------------------------------------------------------------------
      !  Map the particles onto the finest tree topology = topoid=nlevel
      !-------------------------------------------------------------------------

      mapt = ppm_param_map_global
      CALL ppm_map_part(xp,ppm_dim,Np,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF
      mapt = ppm_param_map_push
      CALL ppm_map_part(wp,Np,Mpart,topoid,mapt,info)   ! strengths
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      CALL ppm_map_part(boxpart,Np,Mpart,topoid,mapt,info)   ! boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      mapt = ppm_param_map_send
      CALL ppm_map_part(wp,Np,Mpart,topoid,mapt,info)   ! send
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to send particles.',info)
          GOTO 9999
      ENDIF
      mapt = ppm_param_map_pop
      CALL ppm_map_part(boxpart,Np,Mpart,topoid,mapt,info)   ! boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      CALL ppm_map_part(wp,Np,Mpart,topoid,mapt,info)   ! strengths
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to pop strengths.',info)
          GOTO 9999
      ENDIF
     CALL ppm_map_part(xp,ppm_dim,Np,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to pop positions.',info)
          GOTO 9999
      ENDIF
      Np = Mpart
      IF (ppm_debug .GT. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_fmm_init', &
        &             'Done mapping particles.',info)
        WRITE(cbuf,'(A,I)') 'Local number of particles is now: ',Np
        CALL ppm_write(ppm_rank,'ppm_fmm_init',cbuf,info)
      ENDIF                                                                
      !-------------------------------------------------------------------------
      !  Check that particles have been mapped correctly
      !-------------------------------------------------------------------------
      ! If processor has no more particles, end
      IF (Np .EQ. 0) THEN
         GOTO 9999
      ENDIF

      CALL ppm_topo_check(xp,Np,OK,info)
      IF (info .NE. 0) THEN
         CALL ppm_write(ppm_rank,'ppm_fmm_init', &
         &    'Failed to check topology.',info)
      ENDIF
      IF (.NOT.OK) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Particles not mapped correctly!',info) 
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      ! Recall tree to get the correct lpdx and lhbx arrays
      ! (after mapping particles changed)
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      !-------------------------------------------------------------------------
      ! Calling the tree single (2)
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,        &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_s,max_box_s,lhbx,lpdx,boxcost_s,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_s
      max_box => max_box_s
      boxcost => boxcost_s
#else
      !-------------------------------------------------------------------------
      ! Calling the tree double (2)
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,        &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_d,max_box_d,lhbx,lpdx,boxcost_d,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_d
      max_box => max_box_d
      boxcost => boxcost_d
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init',  &
       &       'Calling tree (2) failed.',__LINE__,info)
          GOTO 9999
      ENDIF
      nrofbox = nbox
!     ouput information!
      IF (ppm_debug.GT.0) THEN
         CALL ppm_write(ppm_rank,'ppm_fmm_init', &
         &    'calling tree (2) successful',info)
         WRITE (cbuf,'(A,I)') 'nbox = ',nbox
         CALL ppm_write(ppm_rank,'ppm_fmm_init',cbuf,info)
      ENDIF

      ENDIF ! ppm_nrpco .GT. 1

       
      !-------------------------------------------------------------------------
      ! Allocating all data (single/double)
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Allocating expansion,sqrtfac,Anm,fracfac (single/double)
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit

#if   __KIND == __SINGLE_PRECISION
      ldl(1)  = 1
      ldl(2)  = 0
      ldl(3)  = -order
      ldu(1)  = nbox
      ldu(2)  = order
      ldu(3)  = order
      CALL ppm_alloc(expansion_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating expansion',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(sqrtfac_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating sqrtfac',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(Anm_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Anm',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldu(1) = order
      CALL ppm_alloc(fracfac_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating frac_fac',__LINE__,info)
      GOTO 9999
      ENDIF
      ! Init array
      expansion_s = (0.0_MK,0.0_MK)
      Anm_s       = (0.0_MK,0.0_MK)
      sqrtfac_s   = (0.0_MK,0.0_MK)
      fracfac_s   = (0.0_MK)
#else
      ldl(1)  = 1
      ldl(2)  = 0
      ldl(3)  = -order
      ldu(1)  = nbox
      ldu(2)  = order
      ldu(3)  = order
      CALL ppm_alloc(expansion_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating expansion',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(sqrtfac_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating sqrtfac',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(Anm_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Anm',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldu(1) = order
      CALL ppm_alloc(fracfac_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating fracfac',__LINE__,info)
      GOTO 9999
      ENDIF
      ! Init array
      expansion_d = (0.0_MK,0.0_MK)
      Anm_d       = (0.0_MK,0.0_MK)
      sqrtfac_d   = (0.0_MK,0.0_MK)
      fracfac_d   = (0.0_MK)
#endif

      !-------------------------------------------------------------------------
      ! Allocating multipole coefficient variables      
      ! Pnm, Cnm, Ynm, fac, rho, theta, phi      
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      ldl(1) = 0
      ldu(1) = 2*order
      CALL ppm_alloc(fac_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating fac',__LINE__,info)
      GOTO 9999
      ENDIF
      ldu1(1) = ppm_nproc*Np
      CALL ppm_alloc(rho_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating rho',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(theta_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating theta',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(phi_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating phi',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(Cnm_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Cnm',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(Ynm_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Ynm',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(Pnm_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Pnm',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldl(2) = -2*order
      ldu(1) = 2*order
      ldu(2) = 2*order
      CALL ppm_alloc(Inner_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Inner',__LINE__,info)
      GOTO 9999
      ENDIF
      fac_s   = (0.0_MK)
      rho_s   = (0.0_MK)
      phi_s   = (0.0_MK)
      theta_s = (0.0_MK)
      Inner_s = (0.0_MK,0.0_MK)
#else
      ldl(1) = 0
      ldu(1) = 2*order
      CALL ppm_alloc(fac_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating fac',__LINE__,info)
      GOTO 9999
      ENDIF
      ldu1(1) = ppm_nproc*Np
      CALL ppm_alloc(rho_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating rho',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(theta_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating theta',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(phi_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating phi',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(Cnm_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Cnm',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(Ynm_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Ynm',__LINE__,info)
      GOTO 9999
      ENDIF
      CALL ppm_alloc(Pnm_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Pnm',__LINE__,info)
      GOTO 9999
      ENDIF
      ldl(1) = 0
      ldl(2) = -2*order
      ldu(1) = 2*order
      ldu(2) = 2*order
      CALL ppm_alloc(Inner_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Inner',__LINE__,info)
      GOTO 9999
      ENDIF
      fac_d   = (0.0_MK)
      rho_d   = (0.0_MK)
      phi_d   = (0.0_MK)
      theta_d = (0.0_MK)
      Inner_d = (0.0_MK,0.0_MK)
#endif

      !-------------------------------------------------------------------------
      ! Allocating diagnoal (single/double)
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit

#if   __KIND == __SINGLE_PRECISION
      ldu1(1) = nbox
      CALL ppm_alloc(radius_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating radius',__LINE__,info)
      GOTO 9999
      ENDIF
#else
      ldu1(1) = nbox
      CALL ppm_alloc(radius_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating radius',__LINE__,info)
      GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      ! Allocating centerofbox (single/double)
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      ldu2(1) = ppm_dim
      ldu2(2) = nbox
      CALL ppm_alloc(centerofbox_s,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating centerofbox',__LINE__,info)
      GOTO 9999
      ENDIF
#else
      ldu2(1) = ppm_dim
      ldu2(2) = nbox
      CALL ppm_alloc(centerofbox_d,ldu2,iopt,info)
      IF (info .NE. 0) THEN 
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &        'error allocating centerofbox',__LINE__,info)
      GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      ! Allocating totalmass (single/double)
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      ldu1(1) = nbox
      CALL ppm_alloc(totalmass_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating totalmass',__LINE__,info)
      GOTO 9999
      ENDIF
#else
      ldu1(1) = nbox
      CALL ppm_alloc(totalmass_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN 
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &        'error allocating totalmass',__LINE__,info)
      GOTO 9999
      ENDIF
#endif

      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init','alloc data successful' &
      &      ,info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Checking precision and pointing to the correct variables
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      centerofbox  => centerofbox_s
      totalmass    => totalmass_s
      radius       => radius_s
      maxboxcost_s = maxboxcost
#else
      centerofbox  => centerofbox_d
      totalmass    => totalmass_d
      radius       => radius_d
      maxboxcost_d = maxboxcost
#endif

      !-------------------------------------------------------------------------
      ! Allocating treepart and treewp
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit

      ldu2(1) = ppm_dim
      ldu2(2) = Np
      CALL ppm_alloc(treepart,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating treepart',__LINE__,info)
      GOTO 9999
      ENDIF
      
      ldu1(1) = Np
      CALL ppm_alloc(treewp,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating treewp',__LINE__,info)
      GOTO 9999
      ENDIF

      centerofbox(1,:) = (0.0_MK)
      centerofbox(2,:) = (0.0_MK)
      centerofbox(3,:) = (0.0_MK)

      DO box=1,nbox
         IF (nchld(box) .NE. 0) CYCLE
         first = lhbx(1,box)
         last  = lhbx(2,box)
         !----------------------------------------------------------------------
         ! Computing new array with particle order from tree
         !----------------------------------------------------------------------
         nrpbox = last-first+1
         DO j=first,last
            treepart(:,j) = xp(:,lpdx(j))
            treewp(j)   = wp(lpdx(j))
         ENDDO
         totalmass(box) = SUM(ABS(treewp(first:last)))
         !----------------------------------------------------------------------
         ! Computing the centers of the leaf boxes
         ! and necessary information (nr. of part. in box)
         !----------------------------------------------------------------------
         IF (nrpbox .GT. 0) THEN
            centerofbox(1,box)=SUM(treepart(1,first:last)* &
            &                  ABS(treewp(first:last)))/ &
            &                  totalmass(box)              
            centerofbox(2,box)=SUM(treepart(2,first:last)* &
            &                  ABS(treewp(first:last)))/ &
            &                  totalmass(box)               
            centerofbox(3,box)=SUM(treepart(3,first:last)* &
            &                  ABS(treewp(first:last)))/ &
            &                  totalmass(box)               
         ELSE
            centerofbox(1,box) = 0.5_MK*(max_box(1,box) + min_box(1,box))
            centerofbox(2,box) = 0.5_MK*(max_box(2,box) + min_box(2,box))
            centerofbox(3,box) = 0.5_MK*(max_box(3,box) + min_box(3,box))
         ENDIF
      ENDDO
      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init','computed centers',info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Computing the radius of the leaf boxes
      !-------------------------------------------------------------------------

      radius = (0.0_MK)

      DO box=1,nbox
         IF (nchld(box) .NE. 0) CYCLE
         first = lhbx(1,box)
         last  = lhbx(2,box)
         DO j=first,last
           diagvec(1) = treepart(1,j) - centerofbox(1,box)
           diagvec(2) = treepart(2,j) - centerofbox(2,box)
           diagvec(3) = treepart(3,j) - centerofbox(3,box)
           tmp = SQRT(diagvec(1)**2 + diagvec(2)**2 + diagvec(3)**2)
           IF (tmp .GT. radius(box)) THEN
              radius(box) = tmp
           ENDIF
         ENDDO
      ENDDO

      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init','computed radius',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  deallocate local variables
      !-------------------------------------------------------------------------

      ldu1(1)=0
      ldu2(1:2) = 0
      istat = 0
      iopt = ppm_param_dealloc
      CALL ppm_alloc(treepart,ldu2,iopt,info)
      istat=istat+info
      CALL ppm_alloc(treewp,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(box2proc,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(cost,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(new_subs2proc,ldu1,iopt,info)
      istat=istat+info

      IF (istat .NE. 0) THEN
          WRITE(cbuf,'(A,I3,A)') 'for ',istat,'error while dealloc'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fmm_init',cbuf,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Nullify data pointers
      !-------------------------------------------------------------------------
      NULLIFY(min_box)
      NULLIFY(max_box)
      NULLIFY(boxcost)
      NULLIFY(centerofbox)
      NULLIFY(totalmass)
      NULLIFY(radius)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_fmm_init',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fmm_init_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fmm_init_d
#endif
