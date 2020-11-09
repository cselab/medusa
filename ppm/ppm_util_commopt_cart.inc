      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_util_commopt
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Determine an approximately optimal communication
      !                 sequence for each processor to SendRecv data with
      !                 its neighbors only. Such that: no conflicts occur
      !                 (i.e. A wants to send to B, but B is currently busy
      !                 receiving from C) and that a minimum number of
      !                 communication rounds are needed. This is done by
      !                 using the Vizing approximation to the minimal edge
      !                 coloring problem of the processor topology graph.
      !
      !  Input        : topoid     (I) topology ID to be optimized
      !                                (internal numbering)
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Routines     : vizing_coloring (libvizing)
      !
      !  Remarks      : 
      !
      !  References   : V.G. Vizing, On an estimate of the chromatic class
      !                 of a p-graph. Discret. Analiz. 3, 1964, pp.25-30.
      !                 (In Russian).
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_commopt_cart.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.15  2004/10/01 16:33:39  ivos
      !  cosmetics.
      !
      !  Revision 1.14  2004/10/01 16:09:13  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.13  2004/08/31 13:29:59  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.12  2004/07/26 07:42:32  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.11  2004/03/19 08:24:22  hiebers
      !  minor change at including files
      !
      !  Revision 1.10  2004/01/23 17:22:13  ivos
      !  Cleanup: (1) updated header, (2) inserted ppm_write and ppm_error, (3)
      !  inserted checks after every allocate, (4) added argument checking.
      !
      !  Revision 1.9  2003/12/17 14:00:57  ivos
      !  bug fix: size check for ilinks changed from nlinks to 2*nlinks.
      !
      !  Revision 1.8  2003/12/16 13:38:48  ivos
      !  bug fix: length of ilinks set to 2*nlinks (because each link has 2 
      !  nodes). Eliminated SIGSEGV in vizing_coloring.
      !
      !  Revision 1.7  2003/12/16 08:48:16  ivos
      !  Replaced grow_preserve with fit_preserve when shrinking lists before
      !  calling C++.
      !
      !  Revision 1.6  2003/12/11 14:01:05  ivos
      !  Bugfixes: (1) MPI bug: ranks .GT. 0 tried to access data that was not
      !  allocated on them. Fixed by placing an IF(rank.EQ.0) around it. (2)
      !  ppm_icommseq now always includes the processor itself as its first entry
      !  and ppm_ncommseq is increased by 1 as map_part_send needs it that way.
      !
      !  Revision 1.5  2003/12/09 11:27:52  hiebers
      !  merged
      !
      !  Revision 1.4  2003/12/09 11:18:41  ivos
      !  Bugfix: status is now in an #ifdef block for __MPI.
      !
      !  Revision 1.3  2003/12/09 09:35:39  ivos
      !  Changed INTENT of info to INOUT.
      !
      !  Revision 1.2  2003/12/09 08:56:50  ivos
      !  First complete version of the communication optimizer utility.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_util_commopt_cart(topoid,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_write
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_check_topoid
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      USE MPI
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! The ID of the topology which is to be optimized
      INTEGER                 , INTENT(IN   ) :: topoid
      ! return status
      INTEGER                 , INTENT(INOUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3,8)               :: lcoords
      LOGICAL, DIMENSION(:), POINTER        :: ichecked
      INTEGER, DIMENSION(3)                 :: coords, lcoords2
      REAL,    DIMENSION(3,26)              :: displ
      INTEGER, DIMENSION(3,26,8)            :: commsequence
      REAL,    DIMENSION(3,3)               :: rotxy,rotyz
      REAL,    DIMENSION(3)                 :: rdispl,dsplc
      REAL,    DIMENSION(4)                 :: angles
      REAL                                  :: M_PI
      INTEGER                               :: idir,icpu,iangle,jangle,lrank
      INTEGER                               :: idirect,lidir
      !-----------------------------------------------------
      REAL(ppm_kind_double)                 :: t0
      INTEGER, DIMENSION(2)                 :: ldu
      INTEGER                               :: iopt
      ! number of neighborhood relations in total
      INTEGER                               :: nlinks
      ! DISTINCT neighbor pairs, i.e. 1<->2 and 2<->1 is the same and only
      ! listed once. even indices are first points, odd ones second ones of
      ! the same edges.
      INTEGER, DIMENSION(:  ) , POINTER     :: ilinks
      ! optimal edge coloring determined. sequence of triples (p1,p2,c),...
      ! with p1 and p2 being the 2 vertices of each edge and c its color.
      INTEGER, DIMENSION(:  ) , POINTER     :: optres
      ! number of neighbors of all every CPU. index: MPI rank
      INTEGER, DIMENSION(:  ) , POINTER     :: nneighprocs
      ! all neighbors of all processors. 1st index: neighbor nr., 2nd:
      ! processor rank
      INTEGER, DIMENSION(:,:) , POINTER     :: ineighprocs
      INTEGER                               :: i,j,maxneigh,isize,ii,isin,k
      ! processor ranks
      INTEGER                               :: p1,p2,isb,jsb,ksb
      ! min and max of assigned colors
      INTEGER                               :: mincolor,maxcolor
      LOGICAL                               :: valid
#ifdef __MPI
      ! MPI comm status
      INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: status
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      !-----------------------------------------------------
      !  memory measurements
      !-----------------------------------------------------
      INTEGER    :: eye
      INTEGER*4  :: fragments
      INTEGER*8  :: total_free, largest_free, total_used
      INTEGER    :: heap_info
      
      CHARACTER(len=256) :: mesg
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_commopt',t0,info)
      M_PI = ACOS(-1.0)
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_commopt',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if there are more than 1 processor. If not, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc .LT. 2) THEN
          !---------------------------------------------------------------------
          !  Allocate memory for communication protocols
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_grow_preserve
          ldu(1) = 1
          ldu(2) = ppm_max_topoid
          CALL ppm_alloc(ppm_icommseq,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'communication sequence PPM_ICOMMSEQ',__LINE__,info)
              GOTO 9999
          ENDIF
          iopt   = ppm_param_alloc_grow_preserve
          ldu(1) = ppm_max_topoid
          CALL ppm_alloc(ppm_ncommseq,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'number of comm.rounds PPM_NCOMMSEQ',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Set the trivial protocol: only myself
          !---------------------------------------------------------------------
          ppm_ncommseq(topoid) = 1
          ppm_icommseq(1,topoid) = ppm_rank
          GOTO 9999
      END IF

#ifdef __MPI
      !-----------------------------------------------------
      !  Allocate the checked array
      !-----------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nproc
      CALL ppm_alloc(ichecked,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
              &            'ichecked array',__LINE__,info)
         GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      icpu = 1
      DO k=1,2
         DO j=1,2
            DO i=1,2
               lcoords(:,icpu) = (/i-1,j-1,k-1/)
               icpu = icpu + 1
            END DO
         END DO
      END DO
      rdispl    = (/1.0,0.0,0.0/)
      angles(1) = 0.0
      angles(2) = 0.25 * M_PI
      angles(3) = 0.5  * M_PI
      angles(4) = 0.75 * M_PI
      idir = 0
      DO jangle=1,4
         DO iangle=1,4
            rotxy(:,1) = (/ COS(angles(jangle)), SIN(angles(jangle)),0.0/)
            rotxy(:,2) = (/-SIN(angles(jangle)), COS(angles(jangle)),0.0/)
            rotxy(:,3) = (/ 0.0, 0.0, 1.0 /)
            rotyz(:,1) = (/ 1.0, 0.0, 0.0 /)
            rotyz(:,2) = (/0.0, COS(angles(iangle)),SIN(angles(iangle))/)
            rotyz(:,3) = (/0.0,-SIN(angles(iangle)),COS(angles(iangle))/)
            dsplc = MATMUL(MATMUL(rotyz,rotxy),rdispl)
            idir = idir + 1
            IF(idir.LT.4) CYCLE
            displ(:,2*(idir-3)-1) = dsplc
            displ(:,2*(idir-3)  ) = -dsplc
         END DO
      END DO
      !-----------------------------------------------------
      !  truncate
      DO idir=1,26
         DO i=1,3
            IF(ABS(displ(i,idir)).LT.0.05) THEN
               displ(i,idir) = 0.0
            END IF
         END DO
      END DO
      !-----------------------------------------------------
      !  modulate
      DO icpu=1,8
         DO idir=1,26
            IF(ABS(displ(1,idir)).GT.0.0) THEN
               idirect = 1
            END IF
            IF(ABS(displ(2,idir)).GT.0.0) THEN
               idirect = 2
            END IF
            IF(ABS(displ(3,idir)).GT.0.0) THEN
               idirect = 3
            END IF
            IF(MOD(lcoords(idirect,icpu),2).EQ.1) THEN
               IF(displ(1,idir).LT.0.0) THEN
                  commsequence(1,idir,icpu) = -1
               ELSEIF(displ(1,idir).GT.0.0) THEN
                  commsequence(1,idir,icpu) =  1
               ELSE
                  commsequence(1,idir,icpu) =  0
               END IF
               IF(displ(2,idir).LT.0.0) THEN
                  commsequence(2,idir,icpu) = -1
               ELSEIF(displ(2,idir).GT.0.0) THEN
                  commsequence(2,idir,icpu) =  1
               ELSE
                  commsequence(2,idir,icpu) =  0
               END IF
               IF(displ(3,idir).LT.0.0) THEN
                  commsequence(3,idir,icpu) = -1
               ELSEIF(displ(3,idir).GT.0.0) THEN
                  commsequence(3,idir,icpu) =  1
               ELSE
                  commsequence(3,idir,icpu) =  0
               END IF

            ELSE
               IF(displ(1,idir).LT.0.0) THEN
                  commsequence(1,idir,icpu) =  1
               ELSEIF(displ(1,idir).GT.0.0) THEN
                  commsequence(1,idir,icpu) = -1
               ELSE
                  commsequence(1,idir,icpu) =  0
               END IF
               IF(displ(2,idir).LT.0.0) THEN
                  commsequence(2,idir,icpu) =  1
               ELSEIF(displ(2,idir).GT.0.0) THEN
                  commsequence(2,idir,icpu) = -1
               ELSE
                  commsequence(2,idir,icpu) =  0
               END IF
               IF(displ(3,idir).LT.0.0) THEN
                  commsequence(3,idir,icpu) =  1
               ELSEIF(displ(3,idir).GT.0.0) THEN
                  commsequence(3,idir,icpu) = -1
               ELSE
                  commsequence(3,idir,icpu) =  0
               END IF
            END IF
         END DO
      END DO
      !-----------------------------------------------------
      !   we have got the basic pattern now
      !   now evaluate it for this rank
      CALL MPI_CART_COORDS(ppm_comm,ppm_rank,3,coords,info)
      !-----------------------------------------------------
      ichecked = .FALSE.
      !   allocate
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = ppm_max_topoid
      CALL ppm_alloc(ppm_ncommseq,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
              &        'number of comm. sequences PPM_NCOMMSEQ',__LINE__,info)
         GOTO 9999
      ENDIF
      ppm_ncommseq(topoid) = MIN(27,ppm_nproc)
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = ppm_ncommseq(topoid)   
      ldu(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_icommseq,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
              &        'final communication sequence PPM_ICOMMSEQ',__LINE__,info)
         GOTO 9999
      ENDIF
      !-----------------------------------------------------
      !  what is my icpu
      icpu = 1*MOD(coords(1),2)+2*MOD(coords(2),2)+4*MOD(coords(3),2)+1
      ppm_icommseq(1,topoid) = ppm_rank
      lidir = 1
      DO idir=1,26
         lcoords2 = coords + commsequence(:,idir,icpu)
         !WRITE(mesg,*) coords,'/',lcoords2
         !CALL ppm_write(ppm_rank,'commopt',mesg,info)
         lrank = -1
         CALL MPI_CART_RANK(ppm_comm,lcoords2,lrank,info)
         !WRITE(mesg,*) 'rank =>',lrank,'(',info,')'
         !CALL ppm_write(ppm_rank,'commopt',mesg,info)
         IF(ichecked(lrank+1)) CYCLE
         ppm_icommseq(lidir+1,topoid) = lrank
         ichecked(lrank+1) = .TRUE.
         lidir = lidir + 1
      END DO
      ppm_ncommseq(topoid) = lidir
      ! IF(ppm_rank.EQ.0) THEN
!          WRITE(mesg,*) 'ppm_ncommseq is ',ppm_ncommseq(topoid),lidir
!          CALL ppm_write(ppm_rank,'ppm_util_commopt_cart',mesg,info)
!          DO idir=1,ppm_ncommseq(topoid)
!             WRITE(mesg,*) ' entry ',ppm_icommseq(idir,topoid)
!             CALL ppm_write(ppm_rank,'ppm_util_commopt_cart',mesg,info)
!          END DO
!       END IF

      !  done
      !-----------------------------------------------------
      !-------------------------------------------------------------------------
      !  Everybody gets the memory needed
      !-------------------------------------------------------------------------
!       WRITE(mesg,*) 'ppm_icommseq ',ppm_icommseq(1:9,topoid)
!       CALL ppm_write(ppm_rank,'commopt',mesg,info)
!       WRITE(mesg,*) 'ppm_icommseq ',ppm_icommseq(10:18,topoid)
!       CALL ppm_write(ppm_rank,'commopt',mesg,info)
!       WRITE(mesg,*) 'ppm_icommseq ',ppm_icommseq(19:27,topoid)
!       CALL ppm_write(ppm_rank,'commopt',mesg,info)
      
          
      !-------------------------------------------------------------------------
      !  Mark this topology as done
      !-------------------------------------------------------------------------
      ppm_isoptimized(topoid) = .TRUE.

      !-------------------------------------------------------------------------
      !  Deallocate temporary storage
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ichecked,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
              &            'number of comm.rounds PPM_NCOMMSEQ',__LINE__,info)
         GOTO 9999
      ENDIF

#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_commopt',t0,info)
      RETURN
    END SUBROUTINE ppm_util_commopt_cart
