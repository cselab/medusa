!-------------------------------------------------------------------------------
!  Subroutine   :                  pse_check_abort
!-------------------------------------------------------------------------------
!
!  Purpose      : Rank 0 checks if a particular files exists on
!                 the local directory. If so, all processors halt.
!
!  Input        :
!
!  Input/output : 
!
!  Output       : lExist     (L) .TRUE. if the file exists, .FALSE.
!                                otherwise
!                 info       (I) return status (0 if no error)
!
!  Routines     :
!
!  Remarks      : Could check validity of file name (length etc). The
!                 Abort name conflicts with the C routine on some 
!                 systems, thus we rename it to ChkAbort.
!
!  References   :
!
!  Revisions    :
!-------------------------------------------------------------------------------
!  $Log: wvic_check_abort.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:46  menahel
!  initial import
!
!  Revision 1.1  2005/09/28 11:40:27  michaebe
!  Fork from ppm_pvc
!
!-------------------------------------------------------------------------------
!  Ivo F. Sbalzarini
!  Institute of Computational Science
!  ETH Zentrum, Hirschengraben 84
!  CH-8092 Zurich, Switzerland
!-------------------------------------------------------------------------------

SUBROUTINE wvic_check_abort(lExist,info)
  
  !-----------------------------------------------------------------------------
  !  Modules
  !-----------------------------------------------------------------------------
  USE module_wvic
  USE MPI

  !-----------------------------------------------------------------------------
  !  Arguments     
  !-----------------------------------------------------------------------------
  LOGICAL, INTENT(  OUT)          :: lExist
  INTEGER, INTENT(  OUT)          :: info
  !-----------------------------------------------------------------------------
  !  Local variables 
  !-----------------------------------------------------------------------------
  INTEGER          :: ios
  !-----------------------------------------------------------------------------
  !  Externals 
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !  Initialise 
  !-----------------------------------------------------------------------------
  lExist = .FALSE.
  
  !-----------------------------------------------------------------------------
  !  Check if file exists
  !-----------------------------------------------------------------------------
  IF (rank.EQ.0) THEN
     INQUIRE(FILE='ABORT',EXIST=lExist,IOSTAT=ios)
  ENDIF
  !-----------------------------------------------------------------------------
  !  Broad cast the result
  !-----------------------------------------------------------------------------
  CALL MPI_BCast(ios   ,1,MPI_INTEGER,0,comm,info)
  CALL MPI_BCast(lExist,1,MPI_LOGICAL,0,comm,info)

  
  !-----------------------------------------------------------------------------
  !  Check ios
  !-----------------------------------------------------------------------------
  IF (ios.NE.0) THEN
     info = -1
     GOTO 9999
  ENDIF
  
  !-----------------------------------------------------------------------------
  !  Return 
  !-----------------------------------------------------------------------------
9999 CONTINUE
  RETURN
END SUBROUTINE wvic_check_abort
