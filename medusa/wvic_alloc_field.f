!-------------------------------------------------------------------------------
!* filename: wvic_alloc_field                                                *!
!* project : ppm                                                              *!
!* purpose : well, allocate a field                                           *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 18:17:05 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
! $Log: wvic_alloc_field.F,v $
! Revision 1.1.1.1  2006/07/25 15:13:46  menahel
! initial import
!
! Revision 1.2  2005/11/11 17:09:52  michaebe
! gave alloc an lda argument
!
! Revision 1.1  2005/09/28 11:40:26  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
!  Allocate a field
!  [ you need an interface for this ]
!-------------------------------------------------------------------------------
SUBROUTINE wvic_alloc_field(vfield_up, ilda, info)

  USE module_wvic
  IMPLICIT NONE

  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: vfield_up
  INTEGER                   , INTENT(out) :: info
  INTEGER                   , INTENT(in ) :: ilda


  INTEGER, DIMENSION(3)                   :: vndata, ldl, ldu
  
  IF(ASSOCIATED(vfield_up)) THEN
     DEALLOCATE(vfield_up,stat=info)
  END IF

  IF(info.NE.0) GOTO 9999

  vndata(1) = MAXVAL(ndata(1,isublist(1:nsublist)))
  vndata(2) = MAXVAL(ndata(2,isublist(1:nsublist)))
  vndata(3) = MAXVAL(ndata(3,isublist(1:nsublist)))

  ldl = 1 - ghostsize
  ldu = vndata + ghostsize
  
  ALLOCATE(vfield_up(ilda,ldl(1):ldu(1),ldl(2):ldu(2),ldl(3):ldu(3),nsublist),&
       & stat = info)
9999 CONTINUE

END SUBROUTINE wvic_alloc_field

  
SUBROUTINE wvic_alloc_field_s(vfield_up, info)

  USE module_wvic
  IMPLICIT NONE

  REAL(mk), DIMENSION(:,:,:,:), POINTER :: vfield_up
  INTEGER                   , INTENT(out) :: info


  INTEGER, DIMENSION(3)                   :: vndata, ldl, ldu
  
  IF(ASSOCIATED(vfield_up)) THEN
     DEALLOCATE(vfield_up,stat=info)
  END IF

  IF(info.NE.0) GOTO 9999


  vndata(1) = MAXVAL(ndata(1,isublist(1:nsublist)))
  vndata(2) = MAXVAL(ndata(2,isublist(1:nsublist)))
  vndata(3) = MAXVAL(ndata(3,isublist(1:nsublist)))

  ldl = 1 - ghostsize
  ldu = vndata + ghostsize
 
  ALLOCATE(vfield_up(ldl(1):ldu(1),ldl(2):ldu(2),ldl(3):ldu(3),nsublist),&
       & stat = info)

9999 CONTINUE

END SUBROUTINE wvic_alloc_field_s

  
