#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_fmm_expchange
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine exchanges the expansion coefficents,
      !                 the radius and the centerofboxes 
      !                 for the fmm module.
      !
      !  Input        : order        (I) : expansion order
      !                 prec         (F) : dummy to determine precision
      !                                    
      !  Input/output : info         (I) : return status, 0 on success
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_expchange.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:19  menahel
      !  initial import
      !
      !  Revision 1.10  2005/09/19 13:03:28  polasekb
      !  code cosmetics
      !
      !  Revision 1.9  2005/09/11 18:05:30  polasekb
      !  (final?) corrected version
      !  (also works parallel :-)
      !
      !  Revision 1.8  2005/09/11 11:44:20  polasekb
      !  also communicating radius and centerofbox
      !
      !  Revision 1.7  2005/08/30 08:48:16  polasekb
      !  removed debug output
      !
      !  Revision 1.6  2005/08/25 14:16:02  polasekb
      !  corrected size of send/recv buffers
      !  exchanged sendrank/recvrank
      !
      !  Revision 1.4  2005/08/23 14:30:28  polasekb
      !  now making difference between single/double precision
      !
      !  Revision 1.3  2005/08/23 07:56:24  polasekb
      !  added #ifdef __MPI where needed
      !
      !  Revision 1.2  2005/08/23 07:49:26  polasekb
      !  corrected error output
      !
      !  Revision 1.1  2005/05/27 07:57:54  polasekb
      !  initial implementation
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fmm_expchange_s(order,prec,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fmm_expchange_d(order,prec,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
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
      INTEGER                 , INTENT(IN   ) :: order
      REAL(MK)                , INTENT(IN   ) :: prec
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                              :: i,j,k,isub,m,n
      INTEGER                              :: nsend,nrecv,curtopoid,box
      INTEGER                              :: sendrank,recvrank
      INTEGER                              :: nsendexp,nrecvexp
      INTEGER                              :: nsendrad,nrecvrad
      INTEGER                              :: nsendcen,nrecvcen
      INTEGER                              :: iopt,level,cnt
      INTEGER                              :: tag1,tag2,istat
      INTEGER, DIMENSION(1)                :: ldu1
      INTEGER, DIMENSION(2)                :: ldu2
      INTEGER, DIMENSION(3)                :: ldu3
      REAL(MK)                             :: t0
      REAL(MK),DIMENSION(:  ),     POINTER :: radius,recvrad,sendrad
      REAL(MK),DIMENSION(:,:),     POINTER :: recvcen,sendcen,centerofbox
      COMPLEX(MK),DIMENSION(:,:,:),POINTER :: expansion,recvexp,sendexp
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE)  :: status
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_fmm_expchange',t0,info)

      !-------------------------------------------------------------------------
      !  pointing to correct variables (single/double)
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      expansion   => expansion_s
      centerofbox => centerofbox_s
      radius      => radius_s
#else
      expansion   => expansion_d
      centerofbox => centerofbox_d
      radius      => radius_d
#endif

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ppm_nsendlist = ppm_nproc
      ppm_nrecvlist = ppm_nproc
      ldu1(1)        = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu1,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange',     &
     &        'send list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_irecvlist,ldu1,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange',     &
     &        'receive list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      tag1   = 100
      tag2   = 200

      !-------------------------------------------------------------------------
      !  compute top level topology
      !-------------------------------------------------------------------------

      DO i=1,nlevel
        IF (nbpl(i) .GE. ppm_nproc) THEN
           level = i
           EXIT
        ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !   Compute nsendexp,nsendcen,nsendrad (nr of exp.,cen.,rad. to be sent)
      !-------------------------------------------------------------------------
      nsendexp = 0
      DO i=level,nlevel
         !----------------------------------------------------------------------
         !  Get the ppm internal topoid
         !----------------------------------------------------------------------
         curtopoid = ppm_internal_topoid(i)
         nsendexp = nsendexp + ppm_nsublist(curtopoid)
         nsendcen = nsendexp
         nsendrad = nsendexp
      ENDDO
      !-------------------------------------------------------------------------
      !   Set up own lists for sending to other processors
      !-------------------------------------------------------------------------
      ldu3(1) = nsendexp
      ldu3(2) = order+1
      ldu3(3) = 2*order+1
      CALL ppm_alloc(sendexp,ldu3,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
      &       'error allocating sendexp',__LINE__,info)
      GOTO 9999
      ENDIF 

      ldu2(1) = 3
      ldu2(2) = nsendcen
      CALL ppm_alloc(sendcen,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
      &       'error allocating sendcen',__LINE__,info)
      GOTO 9999
      ENDIF 

      ldu1(1) = nsendrad
      CALL ppm_alloc(sendrad,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
      &       'error allocating sendrad',__LINE__,info)
      GOTO 9999
      ENDIF 
      
      !-------------------------------------------------------------------------
      !  loop over all topologies, levels of tree
      !-------------------------------------------------------------------------
      cnt = 0
      DO i=level,nlevel
         !----------------------------------------------------------------------
         !  Get the ppm internal topoid
         !----------------------------------------------------------------------
         curtopoid = ppm_internal_topoid(i)
         !----------------------------------------------------------------------
         !   loop over local subs and store in sendexp,sendcen,sendrad
         !----------------------------------------------------------------------
         DO j=1,ppm_nsublist(curtopoid)
            box = ppm_boxid(ppm_isublist(j,curtopoid),i)
            cnt = cnt + 1
            sendexp(cnt,:,:) = expansion(box,:,:)
            sendcen(:,cnt)   = centerofbox(:,box)
            sendrad(cnt)     = radius(box)
         ENDDO
      ENDDO
      
      !-------------------------------------------------------------------------
      !   Initialize the buffers
      !-------------------------------------------------------------------------
      sendrank           = ppm_rank - 1
      recvrank           = ppm_rank + 1
      ppm_nsendlist      = 0
      ppm_nrecvlist      = 0

      !-------------------------------------------------------------------------
      !  Since we skip the local processor entirely, increment the pointers once
      !-------------------------------------------------------------------------
      sendrank                     = sendrank + 1
      recvrank                     = recvrank - 1
      ppm_nsendlist                = ppm_nsendlist + 1
      ppm_isendlist(ppm_nsendlist) = sendrank
      ppm_nrecvlist                = ppm_nrecvlist + 1
      ppm_irecvlist(ppm_nrecvlist) = recvrank

      !-------------------------------------------------------------------------
      !  Loop over all processors but skip the processor itself
      !-------------------------------------------------------------------------
      DO i=2,ppm_nproc
         !----------------------------------------------------------------------
         !  compute the next processor
         !----------------------------------------------------------------------
         sendrank = sendrank + 1
         IF (sendrank.GT.ppm_nproc-1) sendrank = sendrank - ppm_nproc 
         recvrank = recvrank - 1
         IF (recvrank.LT.          0) recvrank = recvrank + ppm_nproc 

         !----------------------------------------------------------------------
         !  Store the processor to which we will send to
         !----------------------------------------------------------------------
         ppm_nsendlist                = ppm_nsendlist + 1 
         ppm_isendlist(ppm_nsendlist) = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will recv from
         !----------------------------------------------------------------------
         ppm_nrecvlist                = ppm_nrecvlist + 1 
         ppm_irecvlist(ppm_nrecvlist) = recvrank
  
         !----------------------------------------------------------------------
         !  reset counter for nr of exp.coeff.,rad,centers to be received
         !----------------------------------------------------------------------
         nrecvexp                     = 0
         nrecvcen                     = 0
         nrecvrad                     = 0
         
         !----------------------------------------------------------------------
         !  loop over all topologies and check all subs 
         !----------------------------------------------------------------------
         DO j=level,nlevel
            !-------------------------------------------------------------------
            !  Get the ppm internal topoid
            !-------------------------------------------------------------------
            curtopoid = ppm_internal_topoid(j)
            DO isub=1,ppm_nsubs(curtopoid)
              !-----------------------------------------------------------------
              !  Check if they belong to the processor from where we will rcv
              !  data
              !-----------------------------------------------------------------
              IF (ppm_subs2proc(isub,curtopoid) .EQ. recvrank) THEN
                 !--------------------------------------------------------------
                 !  If yes, increase counter for correct allocation
                 !--------------------------------------------------------------
                 nrecvexp = nrecvexp + 1
                 nrecvcen = nrecvcen + 1
                 nrecvrad = nrecvrad + 1
              ELSE
                 !--------------------------------------------------------------
                 !  will be exchanged in another round
                 !--------------------------------------------------------------
              ENDIF 
           ENDDO
         ENDDO
         !----------------------------------------------------------------------
         !  Allocate our recv-arrays, exp.,rad., centers to be received
         !----------------------------------------------------------------------
         !iopt   = ppm_param_alloc_grow
         iopt   = ppm_param_alloc_fit

         ldu3(1) = nrecvexp
         ldu3(2) = order+1
         ldu3(3) = 2*order+1
         CALL ppm_alloc(recvexp,ldu3,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
         &       'error allocating recvexp',__LINE__,info)
         GOTO 9999
         ENDIF

         ldu2(1) = 3
         ldu2(2) = nrecvcen
         CALL ppm_alloc(recvcen,ldu2,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
         &       'error allocating recvcen',__LINE__,info)
         GOTO 9999
         ENDIF

         ldu1(1) = nrecvrad
         CALL ppm_alloc(recvrad,ldu1,iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
         &       'error allocating recvrad',__LINE__,info)
         GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  receive the expansions
         !----------------------------------------------------------------------
         nsend = nsendexp*(order+1)*(2*order+1)
         nrecv = nrecvexp*(order+1)*(2*order+1)
#ifdef __MPI
#if   __KIND == __SINGLE_PRECISION
         CALL MPI_SendRecv(sendexp,nsend,MPI_COMPLEX,sendrank,tag2, &
     &                     recvexp,nrecv,MPI_COMPLEX,recvrank,tag2, &
     &                     ppm_comm,status,info)
#else
         CALL MPI_SendRecv(sendexp,nsend,MPI_DOUBLE_COMPLEX,sendrank,tag2,&
     &                     recvexp,nrecv,MPI_DOUBLE_COMPLEX,recvrank,tag2,&
     &                     ppm_comm,status,info)
#endif
#endif
         !----------------------------------------------------------------------
         !  receive the centers
         !----------------------------------------------------------------------
         nsend = nsendcen*3
         nrecv = nrecvcen*3
#ifdef __MPI
         CALL MPI_SendRecv(sendcen,nsend,ppm_mpi_kind,sendrank,tag2, &
     &                     recvcen,nrecv,ppm_mpi_kind,recvrank,tag2, &
     &                     ppm_comm,status,info)
#endif
         !----------------------------------------------------------------------
         !  receive the radius
         !----------------------------------------------------------------------
         nsend = nsendrad
         nrecv = nrecvrad
#ifdef __MPI
         CALL MPI_SendRecv(sendrad,nsend,ppm_mpi_kind,sendrank,tag2,&
     &                     recvrad,nrecv,ppm_mpi_kind,recvrank,tag2,&
     &                     ppm_comm,status,info)
#endif
         !----------------------------------------------------------------------
         !  and store the received data
         !  loop over all topologies and check all subs
         !----------------------------------------------------------------------
         cnt = 0
         DO j=level,nlevel
            !-------------------------------------------------------------------
            !  Get the ppm internal topoid
            !-------------------------------------------------------------------
            curtopoid = ppm_internal_topoid(j)
            DO isub=1,ppm_nsubs(curtopoid)
              !-----------------------------------------------------------------
              !  Check if the sub belongs to the processor from where we 
              !  received data
              !-----------------------------------------------------------------
              IF (ppm_subs2proc(isub,curtopoid) .EQ. recvrank) THEN
                 box = ppm_boxid(isub,j)
                 !--------------------------------------------------------------
                 !  If yes, store
                 !--------------------------------------------------------------
                 cnt = cnt + 1
                 expansion(box,:,:) = recvexp(cnt,:,:)
                 centerofbox(:,box) = recvcen(:,cnt)
                 radius(box)        = recvrad(cnt)
              ENDIF
           ENDDO
         ENDDO
      ENDDO ! end loop over nproc

      !-------------------------------------------------------------------------
      !  Deallocate the memory for the lists
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      ldu3(1) = 0
      ldu3(2) = 0
      ldu3(3) = 0

      istat   = 0
      CALL ppm_alloc(recvexp,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(recvcen,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(recvrad,ldu3,iopt,info)
      istat = istat + info
      IF (istat.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expchange',     &
     &       'recvexp',__LINE__,info)
      ENDIF

      istat   = 0
      CALL ppm_alloc(sendexp,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(sendexp,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(sendexp,ldu3,iopt,info)
      istat = istat + info
      IF (istat.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expchange',     &
     &       'sendexp',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fmm_expchange',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fmm_expchange_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fmm_expchange_d
#endif
