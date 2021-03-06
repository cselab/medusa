#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_map_part_ring
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine maps the particles onto the given
      !                 processors without respect to the topology (i.e. every
      !                 processor gets as many particles as it can handle
      !                 according to ppm_proc_speed).
      !
      !  Input        : xp(:,:)   (F) the position of the particles
      !                 Npart     (I) the number of particles (on the
      !                               local processor)
      !
      !  Input/output : 
      !
      !  Output       : info      (I)  return status. 0 on success.
      !
      !  Remarks      : The first part of the buffer contains the on processor
      !                 data.
      !
      !                 The storing of ppm_map_type is not used (yet)
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_map_part_ring.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.10  2004/11/11 15:26:18  ivos
      !  Moved allocatable work arrays to module.
      !
      !  Revision 1.9  2004/10/01 16:09:08  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.8  2004/07/26 07:42:46  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.7  2004/07/13 10:40:09  oingo
      !  Cosmetics
      !
      !  Revision 1.6  2004/06/10 16:20:02  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.5  2004/05/17 15:48:01  oingo
      !  Changed the check for Npart from .LE. 0 to .LT. 0
      !
      !  Revision 1.4  2004/04/28 11:51:16  oingo
      !  Added a better distribution scheme for the particles that have not
      !  been associated to a processor due to truncation errors by using the
      !  ppm_proc_speed array
      !
      !  Revision 1.3  2004/04/20 12:00:59  oingo
      !  Bugfix: corrected some out-of-bound issues
      !  Added : consistency check for the excess list after distributing the
      !          particles to the processors
      !
      !  Revision 1.2  2004/04/13 15:16:32  oingo
      !  ppm_map_type is now ppm_param_map_global
      !
      !  Revision 1.1  2004/04/07 11:08:14  oingo
      !  Initial release of ppm_map_part_ring module.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ring_s(xp,Npart,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_map_part_ring_d(xp,Npart,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
#ifdef __MPI
      USE MPI
#endif
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      INTEGER                 , INTENT(IN   ) :: Npart
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,ipart
      INTEGER               :: iopt,iset,ibuffer,dim,totalpart
      INTEGER               :: ipos,max_val,max_exc,min_val,min_exc,temp_val
      CHARACTER(ppm_char)   :: mesg
      REAL(MK)              :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_ring',t0,info)
      dim = ppm_dim

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_map_part_ring',  &
     &            'Npart must be >=0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls (not used yet)
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_global

      !-------------------------------------------------------------------------
      !  Allocate memory for particle per processor lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nproc
      CALL ppm_alloc(plist_des,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',  &
     &        'particle per processor desired list PLIST_DES',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(plist_act,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',  &
     &        'particle per processor actual list PLIST_ACT',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(plist_exc,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',  &
     &        'particle per processor excess list PLIST_EXC',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Allocate memory for particle to/from processor lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nproc - 1
      CALL ppm_alloc(srlist1,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',  &
     &        'particle to/from processor list 1 SRLIST1',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(srlist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',  &
     &        'particle to/from processor list 2 SRLIST2',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find out how many particles are on each processor
      !-------------------------------------------------------------------------
      plist_act(:) = 0
      plist_des(:) = 0

      plist_act(ppm_rank + 1) = Npart

#ifdef __MPI
      CALL MPI_AllReduce(plist_act,plist_des,ppm_nproc,MPI_INTEGER, &
     &                   MPI_SUM,ppm_comm,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_mpi_fail,'ppm_map_part_ring', &
     &                  'MPI_AllReduce failed',__LINE__,info)
         GOTO 9999
      ENDIF
#else
      plist_des(ppm_rank + 1) = Npart
#endif

      totalpart = 0
      DO j = 1,ppm_nproc
         plist_act(j) = plist_des(j)
         totalpart = totalpart + plist_act(j)
      ENDDO

      !-------------------------------------------------------------------------
      !  Calculate the desired number of particles per processor based
      !  on ppm_proc_speed and correct a possible truncation error
      !-------------------------------------------------------------------------
      ipart = 0
      DO j = 1,ppm_nproc
         plist_des(j) = INT(REAL(totalpart,ppm_kind_double) * &
     &                      ppm_proc_speed(j - 1))
         ipart = ipart + plist_des(j)
      ENDDO

      IF (totalpart .NE. ipart) THEN
          j = 1
          DO
             IF (j .GT. ppm_nproc) j = 1
             plist_des(j) = plist_des(j) + 1
             ipart = ipart + 1
             IF (ipart .EQ. totalpart) EXIT
             j = j + 1
          ENDDO
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Calculate the particle excess (too much or too little particles) of
      !  each processor
      !-------------------------------------------------------------------------
      DO j = 1,ppm_nproc
         plist_exc(j) = plist_act(j) - plist_des(j)
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Now lets do some magic. Here we have to figure out how many particles
      !  have to be send to whom or received from. srlist1 contains the
      !  processors to send to or receive from, srlist2 contains the amount of
      !  particles
      !-------------------------------------------------------------------------
      ipos = 0
      
      DO
         max_exc = 1
         max_val = plist_exc(max_exc)
         min_exc = 1
         min_val = plist_exc(min_exc)

         !----------------------------------------------------------------------
         !  First find the maximum and minimum excess in the list. The index
         !  of plist_exc corresponds to the processor number + 1.
         !----------------------------------------------------------------------
         DO j = 1,ppm_nproc
            IF (plist_exc(j) .GT. max_val) THEN
                max_val = plist_exc(j)
                max_exc = j
            ELSEIF (plist_exc(j) .LT. min_val) THEN
                min_val = plist_exc(j)
                min_exc = j
            ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  If max_exc and min_exc are equal (e.g. maximum and minimum excess
         !  are the same) there is nothing more to shuffle around and we are
         !  finished.
         !----------------------------------------------------------------------
         IF (max_exc .EQ. min_exc) EXIT

         !----------------------------------------------------------------------
         !  Calculate how many particles can be send/received
         !----------------------------------------------------------------------
         IF (ABS(min_val) .LT. max_val) THEN
             temp_val = ABS(min_val)
         ELSE
             temp_val = max_val
         ENDIF

         !----------------------------------------------------------------------
         !  Fill the lists only if this processor has to send or receive
         !  something. srlist1(ipos) is the processor to send to or receive
         !  from. srlist2(ipos) is the number of particles to send to or
         !  receive from the processor in srlist1(ipos).
         !----------------------------------------------------------------------
         IF (((max_exc - 1) .EQ. ppm_rank) .OR. &
        &    ((min_exc - 1) .EQ. ppm_rank)) THEN
             IF ((max_exc - 1) .EQ. ppm_rank) THEN
                 srlist1(ipos + 1) = min_exc - 1
             ELSE
                 srlist1(ipos + 1) = max_exc - 1
             ENDIF

             srlist2(ipos + 1) = temp_val

             !------------------------------------------------------------------
             !  Count the entries in the srlists
             !------------------------------------------------------------------
             ipos = ipos + 1
         ENDIF

         !----------------------------------------------------------------------
         !  Update the excess list to the new situation. In the end all entries
         !  in plist_exc should be zero.
         !----------------------------------------------------------------------
         plist_exc(max_exc) = plist_exc(max_exc) - temp_val
         plist_exc(min_exc) = plist_exc(min_exc) + temp_val
      ENDDO
      
      !-------------------------------------------------------------------------
      !  Consistency check: all entries in plist_exc should be zero
      !-------------------------------------------------------------------------
      DO j = 1,ppm_nproc
         IF (plist_exc(j) .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_map_incomp,'ppm_map_part_ring', &
     &           'Could not distribute all particles to processors', &
     &           __LINE__,info)
            GOTO 9999
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Calculate the particle excess again since it has been destroyed
      !  during the above caluculation and to figure out if this processor
      !  has to send, receive or do nothing
      !-------------------------------------------------------------------------
      DO j = 1,ppm_nproc
         plist_exc(j) = plist_act(j) - plist_des(j)
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for the pointer to the buffer; we only need to talk
      !  to (ipos + 1) processors (including this one)
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ipos + 2
      CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',     &
     &        'particle send buffer PPM_PSENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = Npart
      CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',     &
     &        'buffer-to-particles map PPM_BUFFER2PART',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the number of buffer entries (this is the first)
      !-------------------------------------------------------------------------
      ppm_buffer_set = 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the field registers that holds the dimension and 
      !  type of the data
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
          GOTO 9999
      ENDIF

      ppm_buffer_dim(ppm_buffer_set)  = dim
      IF (ppm_kind .EQ. ppm_kind_single) THEN
          ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
      ELSE
          ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
      ENDIF

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the buffer
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = dim * Npart
      IF (ppm_kind .EQ. ppm_kind_double) THEN
         CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
      ELSE
         CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
      ENDIF
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',     &
     &        'global send buffer PPM_SENDBUFFER',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      ppm_nsendlist = ipos + 1
      ppm_nrecvlist = ipos + 1

      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',     &
     &        'send list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = ppm_nrecvlist
      CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_map_part_ring',     &
     &        'receive list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Fill the sendbuffer with the particles that will stay on this
      !  processor
      !-------------------------------------------------------------------------
      ppm_psendbuffer(1) = 1
      iset               = 0
      ibuffer            = 0
      
      ppm_isendlist(1) = ppm_rank
      ppm_irecvlist(1) = ppm_rank
      
      IF (plist_exc(ppm_rank + 1) .LE. 0) THEN
          temp_val = Npart
      ELSE
          temp_val = plist_des(ppm_rank + 1)
      ENDIF
      
      DO ipart = 1,temp_val
         !----------------------------------------------------------------------
         !  Store the id of the particle
         !----------------------------------------------------------------------
         ppm_buffer2part(ipart) = ipart

         !----------------------------------------------------------------------
         !  Store the particle
         !----------------------------------------------------------------------
         DO i = 1,dim
            ibuffer                  = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
            ppm_sendbuffers(ibuffer) = xp(i,ipart)
#else
            ppm_sendbufferd(ibuffer) = xp(i,ipart)
#endif
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Lets see if we have anything more to send
      !-------------------------------------------------------------------------
      DO j = 2,ipos + 1
         ppm_psendbuffer(j) = temp_val + 1

         !----------------------------------------------------------------------
         !  Only if we have to send something add it to the sendbuffer
         !----------------------------------------------------------------------
         IF (plist_exc(ppm_rank + 1) .GT. 0) THEN
             DO k = 1,srlist2(j - 1)
                temp_val = temp_val + 1

                !---------------------------------------------------------------
                !  Store the id of the particle
                !---------------------------------------------------------------
                ppm_buffer2part(temp_val) = temp_val

                !---------------------------------------------------------------
                !  Store the particle
                !---------------------------------------------------------------
                DO i = 1, dim
                   ibuffer                  = ibuffer + 1
#if    __KIND == __SINGLE_PRECISION
                   ppm_sendbuffers(ibuffer) = xp(i,temp_val)
#else
                   ppm_sendbufferd(ibuffer) = xp(i,temp_val)
#endif
                ENDDO
             ENDDO
         ENDIF

         ppm_isendlist(j) = srlist1(j - 1)
         ppm_irecvlist(j) = srlist1(j - 1)
      ENDDO

      !-------------------------------------------------------------------------
      !  Everything has to go!
      !-------------------------------------------------------------------------
      ppm_psendbuffer(ipos + 2) = Npart + 1

      !-------------------------------------------------------------------------
      !  Store the current size of the buffer
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ibuffer

      !-------------------------------------------------------------------------
      !  Deallocate the memory for the lists
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(srlist1,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ring',     &
     &        'particle to/from processor list 1 SRLIST1',__LINE__,info)
      ENDIF
      CALL ppm_alloc(srlist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ring',     &
     &        'particle to/from processor list 2 SRLIST2',__LINE__,info)
      ENDIF
      CALL ppm_alloc(plist_des,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ring',     &
     &        'particle per processor desired list PLIST_DES',__LINE__,info)
      ENDIF
      CALL ppm_alloc(plist_act,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ring',     &
     &        'particle per processor actual list PLIST_ACT',__LINE__,info)
      ENDIF
      CALL ppm_alloc(plist_exc,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ring',     &
     &        'particle per processor excess list PLIST_EXC',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_ring',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ring_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_map_part_ring_d
#endif
