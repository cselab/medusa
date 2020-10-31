      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_topo_subs2proc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine assigns the subdomains to the processors
      !                 using the cost assigned to each subdomain. It assigns
      !                 the domains recursively, by assigning domains and their
      !                 physical neighbours. At the end all processors will 
      !                 know who has what.
      !
      !  Input        : cost(:)     (F) the estimated cost associated 
      !                                 with the subdomains
      !                 nneigh(:)   (I) the number of neighbours of a
      !                                 subdomain
      !                 ineigh(:,:) (I) pointers to these neighbours
      !                 nsubs       (I) total number of subdomains
      !
      !  Input/output :                                            
      !
      !  Output       : sub2proc(:) (I) full list of processor 
      !                                 affiliation of subdomains
      !                 isublist(:) (I) list of subdomains assigned to
      !                                 the local processors
      !                 nsublist    (I) the number of subdomains assigned
      !                 info        (I) return status 
      !
      !  Remarks      : i am not sure we need the truncation handling
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_subs2proc.f,v $
      !  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
      !  initial import
      !
      !  Revision 1.3  2005/08/25 07:41:28  hiebers
      !  Bugfix: Completed former bugfix to ensure that each proc really gets one
      !  sub
      !
      !  Revision 1.2  2005/08/24 16:22:31  walther
      !  Bugfix: make sure each proc. gets at least one sub. Introduced the
      !  variable nassigned to count the number of currently assigned subs.
      !
      !  Revision 1.1  2004/07/26 14:12:22  ivos
      !  Renamed from the versions without topo to resolve name conflict
      !  with global ppm_subs2proc data array.
      !
      !  Revision 1.1  2004/07/26 08:56:44  ivos
      !  Renamed.
      !
      !  Revision 1.14  2004/07/26 07:42:38  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.13  2004/03/04 13:45:19  walther
      !  Bug fix: totalcost was not used.
      !
      !  Revision 1.12  2004/02/25 14:48:09  ivos
      !  bugfix: ppm_proc_speed is always double and needs to be converted to
      !  single for the single version of this routine.
      !
      !  Revision 1.11  2004/02/25 14:36:08  ivos
      !  Changed indices: ppm_proc_speed is now ppm_proc_speed(0:nproc-1).
      !
      !  Revision 1.10  2004/02/25 14:21:14  walther
      !  Removed a commented line and added the ppm_proc_speed(1:ppm_nproc)
      !
      !  Revision 1.9  2004/02/12 14:58:09  ivos
      !  bugfix: Now also allocates and returns ineighlist if nproc=1.
      !  Added argument checks and replaced a print with ppm_write.
      !
      !  Revision 1.8  2004/02/02 12:42:30  walther
      !  Had problems with sticky tags ... moved the #include to column 1 ?!
      !
      !  Revision 1.7  2004/01/23 17:24:14  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.6  2004/01/08 17:57:25  walther
      !  Previous release was not enough; rewrote the routine to fix it.
      !
      !  Revision 1.5  2004/01/08 10:13:21  walther
      !  Two bug fixes: 1) before we could end up on an island; this has 
      !  been fixed and 2) the truncation handling towards the end of the 
      !  routine now has a do loop, hence assuming more than one subdomains 
      !  can be left over from the general assignment.
      !
      !  Revision 1.4  2004/01/06 14:08:10  ivos
      !  Changed alloc_grow to alloc_fit for local subs2proc array.
      !
      !  Revision 1.3  2003/12/09 14:44:23  walther
      !  Added comments in header.
      !
      !  Revision 1.2  2003/12/09 13:42:49  hiebers
      !  added deallocation of temporary variables
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_subs2proc_s(cost,nneigh,ineigh,nsubs,sub2proc, &
     &                                  isublist,nsublist,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_subs2proc_d(cost,nneigh,ineigh,nsubs,sub2proc, &
     &                                  isublist,nsublist,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: cost
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: nneigh
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: ineigh
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc,isublist
      INTEGER                 , INTENT(IN   ) :: nsubs
      INTEGER                 , INTENT(  OUT) :: nsublist
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK):: costsum,totalcost,t0
      INTEGER , DIMENSION(:), POINTER :: list => null()
      LOGICAL , DIMENSION(:), POINTER :: not_assigned => null()
      LOGICAL , DIMENSION(:), POINTER :: not_listed => null()
      INTEGER , DIMENSION(ppm_dim) :: ldc
      INTEGER :: i,j,ii,jj,k,iopt
      INTEGER :: istat,isize,rank,nlist,ilist,isub,jsub,nassigned
      CHARACTER(ppm_char) :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_subs2proc',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (nsubs .LE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_subs2proc', &
     &           'nsubs must be > 0',__LINE__,info)
            GOTO 9999
         ENDIF
         DO isub=1,nsubs
             IF (nneigh(isub) .LT. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_topo_subs2proc', &
     &               'nneigh must be >= 0 for all subs',__LINE__,info)
                GOTO 9999
             ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Make sure we have enough memory for the sub2proc 
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubs
      CALL ppm_alloc(sub2proc,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)') 'allocating ',nsubs,' sub2procs failed'
          CALL ppm_error(ppm_err_alloc,'ppm_topo_subs2proc',            &
     &                   mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we only have one processor, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc .EQ. 1) THEN
         !----------------------------------------------------------------------
         !  Allocate the isublist array
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             WRITE (mesg,'(A,I10,A)') 'allocating ',nsubs,' isublist failed'
             CALL ppm_error(ppm_err_alloc,'ppm_topo_subs2proc',            &
     &           mesg,__LINE__,info)
             GOTO 9999
         ENDIF
         !----------------------------------------------------------------------
         !  Assign all subdomains to the processor and return
         !----------------------------------------------------------------------
         nsublist = nsubs
         DO isub=1,nsubs
            isublist(isub) = isub
            sub2proc(isub) = ppm_rank
         ENDDO
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  and allocate some memory for the lists of logicals and initialize 
      !  them to true, since none of the subdomains have been assigned yet.
      !  the list not_assigned indicate if the subdomain has been assigned
      !  to a processor and the list: not_listed holds the subdomains that
      !  have been pushed onto the stack (list(:)) of subs to be assigned - 
      !  we need this stack since one new subdomain will have several 
      !  neighbours and the way subs are assigned to procs requires this stack 
      !-------------------------------------------------------------------------
      CALL ppm_alloc(not_assigned,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)')                                 &
     &    'allocating local array (1) of size ',nsubs,' failed'
          CALL ppm_error(ppm_err_alloc,'ppm_topo_subs2proc',     &
     &                   mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      CALL ppm_alloc(not_listed  ,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)')                                 &
     &    'allocating local array (2) of size ',nsubs,' failed'
          CALL ppm_error(ppm_err_alloc,'ppm_topo_subs2proc',     &
     &                   mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the lists to true
      !-------------------------------------------------------------------------
      not_assigned = .TRUE. 
      not_listed   = .TRUE.

      !-------------------------------------------------------------------------
      !  Allocate some memory for the stack (list(:))
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit 
      isize  = MAXVAL(nneigh)
      ldc(1) = isize
      CALL ppm_alloc(list,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)')                                 &
     &    'allocating local array (3) of size ',isize,' failed'
          CALL ppm_error(ppm_err_alloc,'ppm_topo_subs2proc',     &
     &                   mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the total cost and the average cost per processor
      !-------------------------------------------------------------------------
      totalcost   = 0.0_MK
      DO i=1,nsubs
         totalcost = totalcost + cost(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Push the first subdomain onto the stack 
      !-------------------------------------------------------------------------
      isub             = 1
      nlist            = 1
      ilist            = 0
      list(nlist)      = isub
      nassigned        = 0

      !-------------------------------------------------------------------------
      !  loop over the processors
      !-------------------------------------------------------------------------
      DO rank=0,ppm_nproc-1
         !----------------------------------------------------------------------
         !  initialize the cost of the current processor and reset the stack
         !  to point to the next sub of the previous stack
         !----------------------------------------------------------------------
         costsum          = 0.0_MK
         isub             = list(ilist+1)
         list(1)          = isub
         not_listed       = .TRUE.
         not_listed(isub) = .FALSE.
         ilist            = 0
         nlist            = 1

         !----------------------------------------------------------------------
         !  keep including subdomains until the cost exceeds the average cost
         !  but make sure you get at least one sub and leave at least one sub 
         !  to the remaining procs
         !  the factor 0.5 helps splitting the load more equal amongst the procs
         !----------------------------------------------------------------------
         DO WHILE (nassigned+1.EQ.nsubs-(ppm_nproc-rank-1).OR.                    &
     &            ( nassigned.LT.nsubs-(ppm_nproc-rank-1).AND.                    &
     &      costsum+0.5_MK*cost(isub).LE.totalcost*REAL(ppm_proc_speed(rank),MK)))
            !-------------------------------------------------------------------
            !  Pop the list and assign the next subdomain to the processor: rank
            !-------------------------------------------------------------------
            ilist              = ilist + 1
            isub               = list(ilist)
            costsum            = costsum + cost(isub)
            sub2proc(isub)     = rank
            not_assigned(isub) = .FALSE.

            !-------------------------------------------------------------------
            !  Increment the counter for the number of assigned subs
            !-------------------------------------------------------------------
            nassigned          = nassigned + 1

            !-------------------------------------------------------------------
            !  Reallocate the list to hold the number of neighbours of the sub
            !-------------------------------------------------------------------
            IF (nlist + nneigh(isub).GT.isize) THEN
               iopt   = ppm_param_alloc_grow_preserve
               isize  = nlist + nneigh(isub)
               ldc(1) = isize
               CALL ppm_alloc(list,ldc,iopt,info)
               IF (info.NE.0) THEN
                   info = ppm_error_fatal
                   WRITE (mesg,'(A,I10,A)') &
     &             'allocating local array (3) of size ',isize,' failed'
                   CALL ppm_error(ppm_err_alloc,'ppm_topo_subs2proc',     &
     &                            mesg,__LINE__,info)
                   GOTO 9999
               ENDIF
            ENDIF

            !-------------------------------------------------------------------
            !  include the neighbours of isub in the list of subdomains to
            !  consider
            !-------------------------------------------------------------------
            DO i=1,nneigh(isub)
               jsub = ineigh(i,isub)
               !----------------------------------------------------------------
               !  check if ineigh(i,isub) is already in the list
               !----------------------------------------------------------------
               IF (not_assigned(jsub).AND.not_listed(jsub)) THEN
                  nlist            = nlist + 1
                  list(nlist)      = jsub
                  not_listed(jsub) = .FALSE.
               ENDIF 
            ENDDO

            !-------------------------------------------------------------------
            !  If there is no new elements evailable, it can mean two things:
            !  1) either we are done (have assigned all subs) or 
            !  2) we have stranded on an island and need to get a ferry
            !-------------------------------------------------------------------
            IF (ilist+1.GT.nlist) THEN
               !----------------------------------------------------------------
               !  Find a subdomain that has not yet been assigned
               !  We do that by looping over all the subs and checking if they
               !  have been assigned
               !----------------------------------------------------------------
#ifdef __VECTOR
               !----------------------------------------------------------------
               !  The vector version is a simple do loop: the not assigned with
               !  the higest id will win
               !----------------------------------------------------------------
               isub = 0
               DO i=1,nsubs,1
                  IF (not_assigned(i)) THEN
                     isub = i 
                  ENDIF
               ENDDO
#else
               !----------------------------------------------------------------
               !  The scalar version starts the loop at nsubs to obtain the
               !  save solution as the vector version. Both return isub = 0
               !  if all subs have been assigned - which cannot be and will
               !  be treated as a bug !
               !----------------------------------------------------------------
               isub = nsubs
               DO WHILE (.NOT.not_assigned(isub)) 
                  isub = isub - 1
                  IF (isub.LT.1) EXIT
               ENDDO
#endif
               !----------------------------------------------------------------
               !  if isub is zero all subs have been assigned and we exit the
               !  loop
               !----------------------------------------------------------------
               IF (isub.LT.1) THEN
                  EXIT
               ELSE
                  IF (ppm_debug .GT. 0) THEN
                     WRITE(mesg,'(A,I6)')   &
     &                   'found an island. Chosen sub as new land: ',isub
                     CALL ppm_write(ppm_rank,'ppm_topo_subs2proc',mesg,info)
                  ENDIF
               ENDIF 

               !----------------------------------------------------------------
               !  or we did reach an island and push the new land onto the stack
               !----------------------------------------------------------------
               nlist       = nlist + 1
               list(nlist) = isub
            ENDIF 
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Assign last subs to the last processor
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
         IF (not_assigned(isub)) THEN
            sub2proc(isub)     = ppm_nproc - 1
            not_assigned(isub) = .FALSE.
         ENDIF    
      ENDDO

      !-------------------------------------------------------------------------
      !  Debugging: check that all subdomains have been assigned
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         i = 0
         DO isub=1,nsubs
            IF (not_assigned(isub)) THEN
               i = i + 1 
            ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  If we found any, write an error message and bail out !
         !----------------------------------------------------------------------
         IF (i.NE.0) THEN
            info = ppm_error_error
            WRITE (mesg,'(A,I10,A)') 'missed: ',i,' sub domains'
            CALL ppm_error(ppm_err_subs_map,'ppm_topo_subs2proc',     &
     &                   mesg,__LINE__,info)
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  Do some statistics (only root)
         !----------------------------------------------------------------------
         IF (ppm_rank.EQ.0) THEN
            DO rank=0,ppm_nproc-1
               costsum = 0.0_MK
               DO isub=1,nsubs
                  IF (sub2proc(isub).EQ.rank) THEN
                     costsum = costsum + cost(isub)
                  ENDIF
               ENDDO
               WRITE(mesg,'(A,I5,A,E12.4)') 'cost for processor: ',rank, &
     &                                      ' is ',costsum
               CALL ppm_write(ppm_rank,'ppm_topo_subs2proc',mesg,info) 
            ENDDO
         ENDIF
         !----------------------------------------------------------------------
         !  End of Debugging: check that all subdomains have been assigned
         !----------------------------------------------------------------------
      ENDIF 

      !-------------------------------------------------------------------------
      !  Allocate the local list of subdomains assigned to the local processor
      !-------------------------------------------------------------------------
      ldc(1) = nsubs
      iopt   = ppm_param_alloc_fit
      CALL ppm_alloc(isublist,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)') 'isublist with dim ',nsubs,' failed'
          CALL ppm_error(ppm_err_alloc,'ppm_topo_subs2proc',     &
     &                   mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store a local list of subdomains assigned to the local processor
      !-------------------------------------------------------------------------
      nsublist = 0
      DO i=1,nsubs
         IF (sub2proc(i).EQ.ppm_rank) THEN
            nsublist           = nsublist + 1
            isublist(nsublist) = i
         ENDIF 
      ENDDO

      !-------------------------------------------------------------------------
      !  Free the memory again
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(list,ldc,iopt,info)
      CALL ppm_alloc(not_assigned,ldc,iopt,info)
      CALL ppm_alloc(not_listed  ,ldc,iopt,info)

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_subs2proc',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_subs2proc_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_subs2proc_d
#endif
