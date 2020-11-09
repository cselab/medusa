#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_field_send
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs the actual send/recv of the 
      !                 mesh blocks and all pushed data.
      !
      !  Input        : 
      !
      !  Input/output : 
      ! 
      !  Output       : info    (I) return status. 0 upon success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data. 
      !
      !                 Using non-blocking communication (ISend,Irecv)
      !                 could possibly speed-up this routine since only a
      !                 Send/Recv/SendRecv is done where needed (we know a
      !                 priori how much data we are going to receive).
      !                 Skipped actions could then already loop to the next
      !                 item and would not have to wait for others.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_field_send_noblock.f,v $
      !  Revision 1.4  2006/09/04 13:29:37  pchatela
      !  Bugfix: send buffer offsets fixed
      !  Removed the wrong non blocking snipplets
      !
      !  Revision 1.3  2006/08/31 16:16:47  pchatela
      !  Bugfix in the "Waitany" version
      !  -> two lists of requests, one for sends, the other for recvs
      !
      !  Revision 1.2  2006/08/24 11:27:50  pchatela
      !  Added the WAIT_ANY non-blocking loop for the single-prec case
      !
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.15  2004/11/11 15:26:16  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.14  2004/10/14 10:32:42  ivos
      !  Local buffers are now never smaller than 1 (MPI_SendRecv could have
      !  problems otherwise).
      !
      !  Revision 1.13  2004/10/01 16:33:36  ivos
      !  cosmetics.
      !
      !  Revision 1.12  2004/10/01 16:09:06  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.11  2004/09/22 10:43:24  ivos
      !  send/recv-buffers are now deallocated after their use to save
      !  memory. If this is a performance issue and memory is not a
      !  problem, we can remove this again.
      !
      !  Revision 1.10  2004/08/30 09:27:04  ivos
      !  fix: replaced hard-coded tag1=600 by tag1=k in order to avoid
      !  communication-round cross talk.
      !
      !  Revision 1.9  2004/08/25 16:13:33  michaebe
      !  Bugfix for #0000024
      !
      !  Revision 1.8  2004/07/26 07:42:43  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/04/05 08:52:11  ivos
      !  bugfix: added missing 3d case when counting size of on-processor
      !  data (thanks to simone for reporting this).
      !
      !  Revision 1.6  2004/02/25 11:52:40  hiebers
      !  bug fix in write format
      !
      !  Revision 1.5  2004/02/24 15:42:17  ivos
      !  bugfix: mpi tag was wrong. now is the same for all cases, but different
      !  in each communication round. Communication does not deadlock any more
      !  on an odd number of processors.
      !
      !  Revision 1.4  2004/02/23 12:19:01  ivos
      !  Several bugs fixed. Tested on 2 processors with a scalar field.
      !  Added debug output in several places.
      !
      !  Revision 1.3  2004/02/20 16:24:55  ivos
      !  bugfix: Ndata was not calculated and ppm_nrecvbuffer thus initialized
      !  wrong.
      !
      !  Revision 1.2  2004/02/20 09:10:14  ivos
      !  bugfix: confused single and double buffer declaration fixed.
      !
      !  Revision 1.1  2004/02/17 16:12:32  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_field_send_noblock(info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      USE MPI
#endif
      integer, parameter :: mk = ppm_kind_double 
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER              , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ibuffer,jbuffer,bdim,offs
      INTEGER               :: iopt,tag1,Ndata,msend,mrecv,allsend,allrecv
      CHARACTER(ppm_char)   :: mesg
      REAL(ppm_kind_double) :: t0
#ifdef __MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE) :: commstat
#endif
      !-----------------------------------------------------
      !  Non-blocking
      !-----------------------------------------------------
      INTEGER               :: kbuffer, kreq, rindex
      INTEGER, DIMENSION(:), POINTER :: klist => null()
      INTEGER, DIMENSION(:), POINTER :: isendoff => null()
      INTEGER, DIMENSION(:), POINTER :: irecvoff => null()
      INTEGER, DIMENSION(:), POINTER :: request => null()
      ! 3 lines for new separate request arrays
      INTEGER               :: kreqrecv, kreqsend
      INTEGER, DIMENSION(:), POINTER :: sendreqs => null()
      INTEGER, DIMENSION(:), POINTER :: recvreqs => null()
      INTEGER, DIMENSION(:), POINTER :: klistsend => null()
      INTEGER, DIMENSION(:), POINTER :: klistrecv => null()
      REAL(ppm_kind_double), DIMENSION(:), POINTER :: isendd => null()
      REAL(ppm_kind_double), DIMENSION(:), POINTER :: irecvd => null()
      REAL(ppm_kind_single), DIMENSION(:), POINTER :: isends => null()
      REAL(ppm_kind_single), DIMENSION(:), POINTER :: irecvs => null()
	  
	  LOGICAL               :: rflag
	  REAL(ppm_kind_double) :: nonblock_time0, nonblock_timer, nonblock_timeout
      
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_field_send_noblock',t0,info)

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit 
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(nsend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'send counter NSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'particle send counter PSEND',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'receive counter NRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(precv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'particle receive counter PRECV',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nrecvlist 
      ldu(2) = ppm_buffer_set 
      CALL ppm_alloc(pp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'work buffer PP',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = ppm_nsendlist 
      ldu(2) = ppm_buffer_set 
      CALL ppm_alloc(qq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'work buffer QQ',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Count the size of the buffer that will not be send
      !-------------------------------------------------------------------------
      Ndata = 0
      IF (ppm_dim .LT. 3) THEN
          !---------------------------------------------------------------------
          !  access mesh blocks belonging to the 1st processor
          !---------------------------------------------------------------------
          DO j=ppm_psendbuffer(1),ppm_psendbuffer(2)-1
             !------------------------------------------------------------------
             !  Get the number of mesh points in this block
             !------------------------------------------------------------------
             Ndata = Ndata + (ppm_mesh_isendblksize(1,j)*    &
     &                ppm_mesh_isendblksize(2,j))
          ENDDO
      ELSE
          !---------------------------------------------------------------------
          !  access mesh blocks belonging to the 1st processor
          !---------------------------------------------------------------------
          DO j=ppm_psendbuffer(1),ppm_psendbuffer(2)-1
             !------------------------------------------------------------------
             !  Get the number of mesh points in this block
             !------------------------------------------------------------------
             Ndata = Ndata + (ppm_mesh_isendblksize(1,j)*    &
     &                ppm_mesh_isendblksize(2,j)*ppm_mesh_isendblksize(3,j))
          ENDDO
      ENDIF
      ibuffer = 0
      DO j=1,ppm_buffer_set
         bdim     = ppm_buffer_dim(j)
         ibuffer  = ibuffer  + bdim*Ndata
      ENDDO 

      !-------------------------------------------------------------------------
      !  Initialize the buffer counters
      !-------------------------------------------------------------------------
      ppm_nrecvbuffer = ibuffer
      nsend(1)        = ibuffer
      nrecv(1)        = ibuffer
      psend(1)        = Ndata
      precv(1)        = Ndata
      mrecv           = -1
      msend           = -1

      !-------------------------------------------------------------------------
      !  Count the size of the buffers that will be sent
      !-------------------------------------------------------------------------
      DO k=2,ppm_nsendlist
          nsend(k) = 0
          psend(k) = 0
          Ndata    = 0

          !---------------------------------------------------------------------
          !  Number of mesh points to be sent off to the k-th processor in
          !  the sendlist
          !---------------------------------------------------------------------
          IF (ppm_dim .LT. 3) THEN
              DO i=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_isendblksize(1,i)*    &
     &                ppm_mesh_isendblksize(2,i))
              ENDDO
          ELSE
              DO i=ppm_psendbuffer(k),ppm_psendbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_isendblksize(1,i)*    &
     &                ppm_mesh_isendblksize(2,i)*ppm_mesh_isendblksize(3,i))
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Store the number of mesh points in psend
          !---------------------------------------------------------------------
          psend(k) = Ndata

          !---------------------------------------------------------------------
          !  Store the size of the data to be sent
          !---------------------------------------------------------------------
          DO j=1,ppm_buffer_set
              nsend(k) = nsend(k) + (ppm_buffer_dim(j)*Ndata)
          ENDDO

          !---------------------------------------------------------------------
          !  Find the maximum buffer length (for the allocate)
          !---------------------------------------------------------------------
          msend = MAX(msend,nsend(k))
          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I9)') 'msend = ',msend
              CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Count the size of the buffers that will be received
      !  (For meshes we know this a priori. This is different from
      !  particles)
      !-------------------------------------------------------------------------
      DO k=2,ppm_nrecvlist
          nrecv(k) = 0
          precv(k) = 0
          Ndata    = 0

          !---------------------------------------------------------------------
          !  Number of mesh points to be received from the k-th processor in
          !  the recvlist
          !---------------------------------------------------------------------
          IF (ppm_dim .LT. 3) THEN
              DO i=ppm_precvbuffer(k),ppm_precvbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_irecvblksize(1,i)*    &
     &                ppm_mesh_irecvblksize(2,i))
              ENDDO
          ELSE
              DO i=ppm_precvbuffer(k),ppm_precvbuffer(k+1)-1
                  Ndata = Ndata + (ppm_mesh_irecvblksize(1,i)*    &
     &                ppm_mesh_irecvblksize(2,i)*ppm_mesh_irecvblksize(3,i))
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Store the number of mesh points in precv
          !---------------------------------------------------------------------
          precv(k) = Ndata

          !---------------------------------------------------------------------
          !  Store the size of the data to be received
          !---------------------------------------------------------------------
          DO j=1,ppm_buffer_set
              nrecv(k) = nrecv(k) + (ppm_buffer_dim(j)*Ndata)
          ENDDO

          !---------------------------------------------------------------------
          !  Find the maximum buffer length (for the allocate)
          !---------------------------------------------------------------------
          mrecv = MAX(mrecv,nrecv(k))
          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I9)') 'mrecv = ',mrecv
              CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Increment the total receive buffer count
          !---------------------------------------------------------------------
          ppm_nrecvbuffer = ppm_nrecvbuffer + nrecv(k)

          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I9)') 'ppm_nrecvbuffer = ',ppm_nrecvbuffer
              CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate the memory for the copy of the particle buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldu(1) = ppm_nrecvbuffer 
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         CALL ppm_alloc(ppm_recvbufferd,ldu,iopt,info)
      ELSE
         CALL ppm_alloc(ppm_recvbuffers,ldu,iopt,info)
      ENDIF 
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'global receive buffer PPM_RECVBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the smaller send and receive buffer
      !-------------------------------------------------------------------------
      ! only allocate if there actually is anything to be sent/recvd with
      ! other processors. Otherwise mrecv and msend would both still be -1
      ! (as initialized above) and the alloc would throw a FATAL. This was
      ! Bug ID 000012.
      IF (ppm_nrecvlist .GT. 1) THEN
          iopt   = ppm_param_alloc_grow
          ldu(1) = MAX(mrecv,1)
          IF (ppm_kind.EQ.ppm_kind_double) THEN
             CALL ppm_alloc(recvd,ldu,iopt,info)
          ELSE
             CALL ppm_alloc(recvs,ldu,iopt,info)
          ENDIF 
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &            'local receive buffer RECV',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      IF (ppm_nsendlist .GT. 1) THEN
          ldu(1) = MAX(msend,1)
          IF (ppm_kind.EQ.ppm_kind_double) THEN
             CALL ppm_alloc(sendd,ldu,iopt,info)
          ELSE
             CALL ppm_alloc(sends,ldu,iopt,info)
          ENDIF 
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &            'local send buffer SEND',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Sum of all mesh points that will be sent and received
      !-------------------------------------------------------------------------
      allsend = SUM(psend(1:ppm_nsendlist))
      allrecv = SUM(precv(1:ppm_nrecvlist))

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main send 
      !  buffer 
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I9)') 'ppm_buffer_set=',ppm_buffer_set
          CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
      ENDIF
      bdim = 0
      offs = 0
      DO k=1,ppm_buffer_set
         offs    = offs + allsend*bdim
         qq(1,k) = offs + 1
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nsendlist
            qq(j,k) = qq(j-1,k) + psend(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9)') 'qq(j,k)=',qq(j,k)
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Compute the pointer to the position of the data in the main receive 
      !  buffer 
      !-------------------------------------------------------------------------
      bdim = 0
      offs = 0
      DO k=1,ppm_buffer_set
         offs    = offs + allrecv*bdim
         pp(1,k) = offs + 1
         bdim    = ppm_buffer_dim(k)
         DO j=2,ppm_nrecvlist
            pp(j,k) = pp(j-1,k) + precv(j-1)*bdim
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9)') 'pp(j,k)=',pp(j,k)
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
                WRITE(mesg,'(A,I9,A,I4)') 'precv(j-1)=',precv(j-1),', bdim=',bdim
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  First copy the on processor data - which is in the first buffer
      !-------------------------------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1
            jbuffer = qq(1,k) - 1
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1
               jbuffer                  = jbuffer + 1
               ppm_recvbufferd(ibuffer) = ppm_sendbufferd(jbuffer)
            ENDDO 
         ENDDO 
      ELSE
         DO k=1,ppm_buffer_set
            ibuffer = pp(1,k) - 1
            jbuffer = qq(1,k) - 1
            DO j=1,psend(1)*ppm_buffer_dim(k)
               ibuffer                  = ibuffer + 1
               jbuffer                  = jbuffer + 1
               ppm_recvbuffers(ibuffer) = ppm_sendbuffers(jbuffer)
            ENDDO 
         ENDDO 
      ENDIF 

      !-----------------------------------------------------
      !  Allocate buffers for non-blocking communication
      !-----------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = 2*ppm_nsendlist
      CALL ppm_alloc(klist,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'request array',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = 2*ppm_nsendlist
      CALL ppm_alloc(request,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'request array',__LINE__,info)
         GOTO 9999
      ENDIF
      
      ! 4 new arrays for separate request lists
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(klistsend,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'k list for send requests',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(klistrecv,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'k list for recv requests',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(sendreqs,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'send requests array',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(recvreqs,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'recv requests array',__LINE__,info)
         GOTO 9999
      ENDIF
      
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(isendoff,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'isendoff array',__LINE__,info)
         GOTO 9999
      ENDIF
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(irecvoff,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
              &            'irecvoff array',__LINE__,info)
         GOTO 9999
      ENDIF

      IF(ppm_kind.EQ.ppm_kind_double) THEN
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nsend(2:ppm_nsendlist))
         CALL ppm_alloc(isendd,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'isendd array',__LINE__,info)
            GOTO 9999
         ENDIF
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nrecv(2:ppm_nsendlist))
         CALL ppm_alloc(irecvd,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'irecvd array',__LINE__,info)
            GOTO 9999
         ENDIF
      ELSE
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nsend(2:ppm_nsendlist))
         CALL ppm_alloc(isends,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'isendd array',__LINE__,info)
            GOTO 9999
         ENDIF
         iopt   = ppm_param_alloc_fit
         ldu(1) = SUM(nrecv(2:ppm_nsendlist))
         CALL ppm_alloc(irecvs,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
                 &            'irecvd array',__LINE__,info)
            GOTO 9999
         ENDIF
      END IF
      
   
      !-----------------------------------------------------
      !  loop over cpus in isendlist, skip the first guy which is the local
      !  processor, put the stuff in one (bigish) send buffer and store the
      !  index
      !-----------------------------------------------------
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         kbuffer = 0
         DO k=2,ppm_nsendlist
            IF(psend(k).GT.0) THEN
                isendoff(k) = kbuffer+1
            ELSE
                isendoff(k) = -1
            END IF
            !-----------------------------------------------------
            !  collect each buffer set
            !-----------------------------------------------------
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  kbuffer = kbuffer + 1
                  jbuffer = jbuffer + 1
                  isendd(kbuffer) = ppm_sendbufferd(jbuffer)
               END DO
            END DO
         END DO
         kbuffer = 1
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  construct irecvoff
            !-----------------------------------------------------
            IF(precv(k).GT.0) THEN
               irecvoff(k) = kbuffer
               kbuffer = kbuffer + nrecv(k)
            ELSE
               irecvoff(k) = -1
            END IF
         END DO
         
#ifdef __MPI
         kreqsend = 0
         kreqrecv = 0
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  perform the actual send/recv as needed,
            !  skip zero-communication rounds (see ppm_map_field_send.f)
            !-----------------------------------------------------
            IF(ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
               tag1 = k
               IF(psend(k).GT.0) THEN
                  kreqsend = kreqsend + 1
                  CALL MPI_ISend(isendd(isendoff(k)),nsend(k),ppm_mpi_kind, &
                       & ppm_isendlist(k),tag1,ppm_comm,sendreqs(kreqsend),info)
                  klistsend(kreqsend) = k
               END IF
               IF(precv(k).GT.0) THEN
                  kreqrecv = kreqrecv + 1
                  CALL MPI_IRecv(irecvd(irecvoff(k)),nrecv(k),ppm_mpi_kind, &
                       & ppm_irecvlist(k),tag1,ppm_comm,recvreqs(kreqrecv),info)
                  klistrecv(kreqrecv) = k
               END IF
            ELSE
               !-----------------------------------------------------
               ! do nothing
               !-----------------------------------------------------
            END IF
         END DO
         IF (ppm_debug .GT. 1) THEN
            WRITE(mesg,'(A,I9,I9)') 'send/recv requests= ', kreqsend, kreqrecv
            CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
#else
#error not implemented for usage without MPI. take ppm_map_field_send.f instead
#endif
         
      ELSE
         !-----------------------------------------------------
         !  single precision case
         !-----------------------------------------------------
         kbuffer = 0
         DO k=2,ppm_nsendlist
            IF(psend(k).GT.0) THEN
                isendoff(k) = kbuffer+1
            ELSE
                isendoff(k) = -1
            END IF
            !-----------------------------------------------------
            !  collect each buffer set
            !-----------------------------------------------------
            DO j=1,ppm_buffer_set
               jbuffer = qq(k,j) - 1
               DO i=1,psend(k)*ppm_buffer_dim(j)
                  kbuffer = kbuffer + 1
                  jbuffer = jbuffer + 1
                  isends(kbuffer) = ppm_sendbuffers(jbuffer)
               END DO
               IF ((psend(k).GT.0).AND.(ppm_debug .GT. 1)) THEN
                  jbuffer = qq(k,j)
                  WRITE(mesg,'(A,I9,I9,F10.4,A,I9)') &
                       'Initial send hash ', &
                       k, &
                       j, &
                       SUM(ppm_sendbuffers( &
                       jbuffer:jbuffer+psend(k)*ppm_buffer_dim(j)-1)), &
                       ' to ',&
                       ppm_isendlist(k)
                  CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
               ENDIF
            END DO
         END DO
         kbuffer = 1
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  construct irecvoff
            !-----------------------------------------------------
            IF(precv(k).GT.0) THEN
               irecvoff(k) = kbuffer
               kbuffer = kbuffer + nrecv(k)
            ELSE
               irecvoff(k) = -1
            END IF
         END DO
         
#ifdef __MPI
         kreqsend = 0
         kreqrecv = 0
         DO k=2,ppm_nsendlist
            !-----------------------------------------------------
            !  perform the actual send/recv as needed,
            !  skip zero-communication rounds (see ppm_map_field_send.f)
            !-----------------------------------------------------
            IF(ppm_isendlist(k).GE.0.AND.ppm_irecvlist(k).GE.0) THEN
               tag1 = k
               IF(psend(k).GT.0) THEN
                  kreqsend = kreqsend + 1
                  CALL MPI_ISend(isends(isendoff(k)), &
                       nsend(k), &
                       ppm_mpi_kind, &
                       ppm_isendlist(k), &
                       tag1, &
                       ppm_comm, &
                       sendreqs(kreqsend), &
                       info)
                  klistsend(kreqsend) = k
               END IF
               IF(precv(k).GT.0) THEN
                  kreqrecv = kreqrecv + 1
                  CALL MPI_IRecv(irecvs(irecvoff(k)), &
                       nrecv(k), &
                       ppm_mpi_kind, &
                       ppm_irecvlist(k), &
                       tag1,ppm_comm,recvreqs(kreqrecv), &
                       info)
                  klistrecv(kreqrecv) = k
               END IF
            ELSE
               !-----------------------------------------------------
               ! do nothing
               !-----------------------------------------------------
            END IF
            END DO
         IF (ppm_debug .GT. 1) THEN
            WRITE(mesg,'(A,I9,I9)') 'send/recv requests= ', kreqsend, kreqrecv
            CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
#else
#error not implemented for usage without MPI. take ppm_map_field_send.f instead
#endif
      ENDIF

      IF(ppm_kind.EQ.ppm_kind_double) THEN
#undef __NONBLOCKDBG
#ifdef __NONBLOCKDBG
         CALL MPI_Barrier(ppm_comm,info)
         nonblock_time0 = MPI_Wtime(info)
         WRITE(*,*) ppm_rank, ': nonblocking ddebug started'
         nonblock_timeout = 20.0_mk
         DO
            CALL MPI_TestAny(kreq,request,rindex,rflag,MPI_STATUSES_IGNORE,info)
            IF (rflag.EQV..TRUE.) THEN
               IF (rindex.EQ.MPI_UNDEFINED) THEN
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
                  WRITE(*,*) ppm_rank, ': all non-blocking stuff completed in ', nonblock_timer
                  EXIT
               ELSE
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
               ENDIF
            ELSE
               nonblock_timer = MPI_Wtime(info) - nonblock_time0
               IF (nonblock_timer.GT.nonblock_timeout) THEN
                  WRITE(*,*) ppm_rank, ': non-blocking stuff timed out'
                  DO k=1,kreq
                     CALL MPI_Test(request(k),rflag,MPI_STATUS_IGNORE,info)
                     IF (.NOT.rflag) THEN
                        WRITE(*,*) ppm_rank, ': was still waiting for recv from ', ppm_irecvlist(klist(k))
                     ENDIF
                  ENDDO
                  CALL MPI_Abort(ppm_comm,info)
               ENDIF
            ENDIF
         ENDDO
         
#endif
         DO 
            CALL MPI_WaitAny(kreqrecv,recvreqs,rindex,MPI_STATUSES_IGNORE,info)
            IF(rindex.EQ.MPI_UNDEFINED) EXIT
            k = klistrecv(rindex)
            kbuffer = irecvoff(k)-1
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9,A,I9)') 'IRecv ', rindex, 'completed, k = ', k
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
            ENDIF
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  jbuffer = jbuffer + 1
                  kbuffer = kbuffer + 1
                  ppm_recvbufferd(jbuffer) = irecvd(kbuffer)
               END DO
            END DO
         END DO
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All IRecv completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
         CALL MPI_Waitall(kreqsend,sendreqs,MPI_STATUSES_IGNORE,info)
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All ISend completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
#ifdef __OLDSTYLE         
         !-----------------------------------------------------
         !  as we have just done non-blocking stuff we need to sync the cpus
         !-----------------------------------------------------
         IF(kreq.GT.0) CALL MPI_Waitall(kreq,request,MPI_STATUSES_IGNORE,info)
         !-----------------------------------------------------
         !  now copy the data back into the real receive buffer
         !-----------------------------------------------------
         kbuffer = 0
         DO k=2,ppm_nsendlist
            kbuffer = irecvoff(k)-1
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  jbuffer = jbuffer + 1
                  kbuffer = kbuffer + 1
                  ppm_recvbufferd(jbuffer) = irecvd(kbuffer)
               END DO
            END DO
         END DO
#endif         

         
         !-----------------------------------------------------
         !  deallocate non-blocking things
         !-----------------------------------------------------
         iopt = ppm_param_dealloc
         CALL ppm_alloc(klist,ldu,iopt,info)
         CALL ppm_alloc(request,ldu,iopt,info)
         ! 4 new arrays
         CALL ppm_alloc(klistsend,ldu,iopt,info)
         CALL ppm_alloc(sendreqs,ldu,iopt,info)
         CALL ppm_alloc(klistrecv,ldu,iopt,info)
         CALL ppm_alloc(recvreqs,ldu,iopt,info)
         
         CALL ppm_alloc(isendoff,ldu,iopt,info)
         CALL ppm_alloc(irecvoff,ldu,iopt,info)
         CALL ppm_alloc(isendd,ldu,iopt,info)
         CALL ppm_alloc(irecvd,ldu,iopt,info)
      ELSE
	  
#ifdef __NONBLOCKDBG
         CALL MPI_Barrier(ppm_comm,info)
         nonblock_time0 = MPI_Wtime(info)
         WRITE(*,*) ppm_rank, ': nonblocking ddebug started'
         nonblock_timeout = 20.0_mk
         DO
            CALL MPI_TestAny(kreq,request,rindex,rflag,MPI_STATUSES_IGNORE,info)
            IF (rflag.EQV..TRUE.) THEN
               IF (rindex.EQ.MPI_UNDEFINED) THEN
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
                  WRITE(*,*) ppm_rank, ': all non-blocking stuff completed in ', nonblock_timer
                  EXIT
               ELSE
                  nonblock_timer = MPI_Wtime(info) - nonblock_time0
               ENDIF
            ELSE
               nonblock_timer = MPI_Wtime(info) - nonblock_time0
               IF (nonblock_timer.GT.nonblock_timeout) THEN
                  WRITE(*,*) ppm_rank, ': non-blocking stuff timed out'
                  DO k=1,kreq
                     CALL MPI_Test(request(k),rflag,MPI_STATUS_IGNORE,info)
                     IF (.NOT.rflag) THEN
                        WRITE(*,*) ppm_rank, ': was still waiting for recv from ', ppm_irecvlist(klist(k))
                     ENDIF
                  ENDDO
                  CALL MPI_Abort(ppm_comm,info)
               ENDIF
            ENDIF
         ENDDO
#endif
         
         DO 
            CALL MPI_WaitAny(kreqrecv,recvreqs,rindex,MPI_STATUSES_IGNORE,info)
            IF(rindex.EQ.MPI_UNDEFINED) EXIT
            k = klistrecv(rindex)
            kbuffer = irecvoff(k)-1
            IF (ppm_debug .GT. 1) THEN
                WRITE(mesg,'(A,I9,A,I9)') 'IRecv ', rindex, 'completed, k = ', k
                CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
            ENDIF
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  jbuffer = jbuffer + 1
                  kbuffer = kbuffer + 1
                  ppm_recvbuffers(jbuffer) = irecvs(kbuffer)
               END DO
               IF ((precv(k).GT.0).AND.(ppm_debug .GT. 1)) THEN
                  jbuffer = pp(k,j)
                  WRITE(mesg,'(A,I9,I9,F10.4,A,I9)') &
                       'Final recv hash ',&
                       k, &
                       j, &
                       SUM(ppm_recvbuffers( &
                       jbuffer:jbuffer+precv(k)*ppm_buffer_dim(j)-1)), &
                       ' from ', &
                       ppm_irecvlist(k)
                  CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
               ENDIF
            END DO
         END DO
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All IRecv completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
         CALL MPI_Waitall(kreqsend,sendreqs,MPI_STATUSES_IGNORE,info)
         IF (ppm_debug .GT. 1) THEN
             WRITE(mesg,'(A)') 'All ISend completed'
             CALL ppm_write(ppm_rank,'ppm_map_field_send_noblock',mesg,info)
         ENDIF
#ifdef __OLDSTYLE 
         !-----------------------------------------------------
         !  as we have just done non-blocking stuff we need to sync the cpus
         !-----------------------------------------------------
         IF(kreq.GT.0) CALL MPI_Waitall(kreq,request,MPI_STATUSES_IGNORE,info)
         !-----------------------------------------------------
         !  now copy the data back into the real receive buffer
         !-----------------------------------------------------
         kbuffer = 0
         DO k=2,ppm_nsendlist
            kbuffer = irecvoff(k)-1
            DO j=1,ppm_buffer_set
               jbuffer = pp(k,j) - 1
               DO i=1,precv(k)*ppm_buffer_dim(j)
                  jbuffer = jbuffer + 1
                  kbuffer = kbuffer + 1
                  ppm_recvbuffers(jbuffer) = irecvs(kbuffer)
               END DO
            END DO
         END DO
#endif
         !-----------------------------------------------------
         !  deallocate non-blocking things
         !-----------------------------------------------------
         iopt = ppm_param_dealloc
		 CALL ppm_alloc(klist,ldu,iopt,info)
         CALL ppm_alloc(request,ldu,iopt,info)
         ! 4 new arrays
         CALL ppm_alloc(klistsend,ldu,iopt,info)
         CALL ppm_alloc(sendreqs,ldu,iopt,info)
         CALL ppm_alloc(klistrecv,ldu,iopt,info)
         CALL ppm_alloc(recvreqs,ldu,iopt,info)
         
         CALL ppm_alloc(isendoff,ldu,iopt,info)
         CALL ppm_alloc(irecvoff,ldu,iopt,info)
         CALL ppm_alloc(isends,ldu,iopt,info)
         CALL ppm_alloc(irecvs,ldu,iopt,info)
      END IF
         
      !-------------------------------------------------------------------------
      !  Deallocate the send buffer to save memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
      ELSE
          CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
      ENDIF
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_map_field_send_noblock',     &
     &        'send buffer PPM_SENDBUFFER',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nsend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'send counter NSEND',__LINE__,info)
      ENDIF
      CALL ppm_alloc(nrecv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'receive counter NRECV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(psend,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'particle send counter PSEND',__LINE__,info)
      ENDIF
      CALL ppm_alloc(precv,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'particle receive counter PRECV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(   pp,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'work array PP',__LINE__,info)
      ENDIF
      CALL ppm_alloc(   qq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'work array QQ',__LINE__,info)
      ENDIF
      CALL ppm_alloc(recvd,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local receive buffer RECVD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(recvs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local receive buffer RECVS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sendd,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local send buffer SENDD',__LINE__,info)
      ENDIF
      CALL ppm_alloc(sends,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_field_send_noblock',     &
     &        'local send buffer SENDS',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_field_send_noblock',t0,info)
      RETURN
      END SUBROUTINE ppm_map_field_send_noblock
