!-------------------------------------------------------------------------------
!* filename: wvic_dumpfield_plot3d                                           *!
!* project : ppm                                                              *!
!* purpose : well                                                             *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Wed Aug 11 17:43:00 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
! $Log: wvic_dumpfield_plot3d.F,v $
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.1  2005/09/28 11:40:34  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------



SUBROUTINE wvic_dumpfield_plot3d (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  INTEGER , INTENT(inout) :: info
  
  CHARACTER(len=256) :: xfile, qfile
  INTEGER            :: ixfile, iqfile
  INTEGER :: i,j,k,m,isub
  INTEGER :: nnx, nny, nnz
  CHARACTER(len=256) :: msg
  !-----------------------------------------------------------------------------
  !  construct filename
  !-----------------------------------------------------------------------------
  if(verbose) CALL ppm_write(rank,'wvic_dumpfield_*','entered',info)
  
  WRITE(qfile,'(A,A,I2.2,A,I4.4,A)') &
       &runtag(1:iruntag),'R',rank,'I',itime,'.q'
  WRITE(xfile,'(A,A,I2.2,A,I4.4,A)') &
       &runtag(1:iruntag),'R',rank,'I',itime,'.x'

  iqfile = 20 + 2*rank
  ixfile = 20 + 2*rank + 1
  
  OPEN(iqfile,file=qfile,form='unformatted')
  !  WARNING:  problem for ncpu>=20 !!
  OPEN(ixfile,file=xfile,form='unformatted')


  !-----------------------------------------------------------------------------
  !  write number of grids
  !-----------------------------------------------------------------------------
  WRITE(iqfile) nsublist
  WRITE(ixfile) nsublist

  !-----------------------------------------------------------------------------
  !  write the grid sizes
  !-----------------------------------------------------------------------------
  WRITE(iqfile) (ndata(1,isublist(isub)),&
       &         ndata(2,isublist(isub)),&
       &         ndata(3,isublist(isub)),isub=1,nsublist)
  WRITE(ixfile) (ndata(1,isublist(isub)),&
       &         ndata(2,isublist(isub)),&
       &         ndata(3,isublist(isub)),isub=1,nsublist)
  DO isub=1,nsublist

     nnx = ndata(1,isublist(isub))
     nny = ndata(2,isublist(isub))
     nnz = ndata(3,isublist(isub))
     
     WRITE(ixfile) (((REAL(k-1,mk)*dx+min_sub(1,isublist(isub)),&
          & k=1,nnx),j=1,nny),i=1,nnz),&
          &        (((REAL(j-1,mk)*dy+min_sub(2,isublist(isub)),&
          & k=1,nnx),j=1,nny),i=1,nnz),&
          &        (((REAL(i-1,mk)*dz+min_sub(3,isublist(isub)),&
          & k=1,nnx),j=1,nny),i=1,nnz)

  END DO

  CLOSE(ixfile)

  !-----------------------------------------------------------------------------
  !  write the values
  !-----------------------------------------------------------------------------
  DO isub=1,nsublist
     nnx = ndata(1,isublist(isub))
     nny = ndata(2,isublist(isub))
     nnz = ndata(3,isublist(isub))
     

     WRITE(iqfile) 0.0_mk, 0.0_mk, 0.0_mk, time

     WRITE(iqfile) ((((SQRT(SUM(field_wp(1:3,k,j,i&
          &,isub)**2)),&
          & k=1,nnx),j=1,nny),i=1,nnz),m=1,1),&
          &((((field_up(m,k,j,i,isub),&
          & k=1,nnx),j=1,nny),i=1,nnz),m=1,3),&
          &(((((field_wp(4,k,j,i&
          &,isub)),&
          & k=1,nnx),j=1,nny),i=1,nnz),m=1,1)

  END DO

  CLOSE(iqfile)
  if(verbose) CALL ppm_write(rank,'wvic_dumpfield_*','exit',info)
  
END SUBROUTINE wvic_dumpfield_plot3d





#ifdef __AGM_TESTRUN__
SUBROUTINE wvic_dumpfield_plot3d_agm (info)

  USE module_wvic
  USE ppm_module_write
  IMPLICIT NONE

  INTEGER , INTENT(inout) :: info
  
  CHARACTER(len=256) :: xfile, qfile
  INTEGER            :: ixfile, iqfile
  INTEGER :: i,j,k,m,isub
  INTEGER :: nnx, nny, nnz
  CHARACTER(len=256) :: msg
  !-----------------------------------------------------------------------------
  !  construct filename
  !-----------------------------------------------------------------------------
  if(verbose) CALL ppm_write(rank,'wvic_dumpfield_*','entered',info)
  
  WRITE(qfile,'(A,A,I2.2,A,I4.4,A)') &
       &runtag(1:iruntag),'R',rank,'I',itime,'.agmq'
  WRITE(xfile,'(A,A,I2.2,A,I4.4,A)') &
       &runtag(1:iruntag),'R',rank,'I',itime,'.amgx'

  iqfile = 20 + 2*rank
  ixfile = 20 + 2*rank + 1
  
  OPEN(iqfile,file=qfile,form='unformatted')
  !  WARNING:  problem for ncpu>=20 !!
  OPEN(ixfile,file=xfile,form='unformatted')


  !-----------------------------------------------------------------------------
  !  write number of grids
  !-----------------------------------------------------------------------------
  WRITE(iqfile) nsublist
  WRITE(ixfile) nsublist

  !-----------------------------------------------------------------------------
  !  write the grid sizes
  !-----------------------------------------------------------------------------
  WRITE(iqfile) (ndata(1,isublist(isub)),&
       &         ndata(2,isublist(isub)),&
       &         ndata(3,isublist(isub)),isub=1,nsublist)
  WRITE(ixfile) (ndata(1,isublist(isub)),&
       &         ndata(2,isublist(isub)),&
       &         ndata(3,isublist(isub)),isub=1,nsublist)
  DO isub=1,nsublist

     nnx = ndata(1,isublist(isub))
     nny = ndata(2,isublist(isub))
     nnz = ndata(3,isublist(isub))
     
     WRITE(ixfile) (((field_xp(3,i,j,k,isub),i=1,nnx),j=1,nny),k=1,nnz),&
          &(((field_xp(2,i,j,k,isub),i=1,nnx),j=1,nny),k=1,nnz),&
          &(((field_xp(1,i,j,k,isub),i=1,nnx),j=1,nny),k=1,nnz)

     

  END DO

  CLOSE(ixfile)

  !-----------------------------------------------------------------------------
  !  write the values
  !-----------------------------------------------------------------------------
  DO isub=1,nsublist
     nnx = ndata(1,isublist(isub))
     nny = ndata(2,isublist(isub))
     nnz = ndata(3,isublist(isub))
     

     WRITE(iqfile) 0.0_mk, 0.0_mk, 0.0_mk, time

     WRITE(iqfile) ((((SQRT(SUM(field_wp(1:3,k,j,i&
          &,isub)**2)),&
          & k=1,nnx),j=1,nny),i=1,nnz),m=1,1),&
          &((((field_up(m,k,j,i,isub),&
          & k=1,nnx),j=1,nny),i=1,nnz),m=1,3),&
          &(((((field_wp(4,k,j,i&
          &,isub)),&
          & k=1,nnx),j=1,nny),i=1,nnz),m=1,1)

  END DO

  CLOSE(iqfile)
  if(verbose) CALL ppm_write(rank,'wvic_dumpfield_*','exit',info)
  
END SUBROUTINE wvic_dumpfield_plot3d_agm
#endif
