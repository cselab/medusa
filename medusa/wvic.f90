!------------------------------------------------------------------------------!
!* filename: ppm_vortex_client                                                *!
!* project : ppm                                                              *!
!* purpose : vortex method client for ppm                                     *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 10 14:12:01 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic.F,v $
!  Revision 1.1.1.1  2006/07/25 15:13:46  menahel
!  initial import
!
!  Revision 1.3  2005/11/11 10:04:06  michaebe
!  adapted call to init
!
!  Revision 1.2  2005/10/14 09:05:03  michaebe
!  cosmetics
!
!  Revision 1.1  2005/09/28 11:40:25  michaebe
!  Fork from ppm_pvc
!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!  Main driver for wall-enabled vortex client
!------------------------------------------------------------------------------!
PROGRAM wvic

  USE module_wvic
  USE ppm_module_finalize
  USE MPI

  INTEGER :: info, taskid
  INTEGER :: steps
  CHARACTER(len=256) :: ctrlfile
#ifdef __WITH_MEMORY_MEASUREMENT__
  !-----------------------------------------------------
  !  memory measurements
  !-----------------------------------------------------
  INTEGER    :: eye
  INTEGER*4  :: fragments
  INTEGER*8  :: total_free, largest_free, total_used
  INTEGER    :: heap_info
#endif  

  steps = 2400000
  WRITE(ctrlfile,'(A4)') 'Ctrl'
#ifdef __WITH_CRAYPAT__
  CALL PAT_HWPC_INIT( taskID, "wvic main")
#endif
  !----------------------------------------------------------------------------!
  ! = init =
  ! double shear layer
  !----------------------------------------------------------------------------!
!  CALL wvic_init(ctrlfile,info)
  CALL wvic_init_cart(ctrlfile,info)
#ifdef __WITH_MEMORY_MEASUREMENT__
  eye = heap_info(fragments, total_free, largest_free, total_used)
  WRITE(6,*) '[after wvic_init] heap_info fragments =',fragments,' total_free =', &
       & total_free/1024/1024,' largest_free = ',largest_free/1024/1024,' TOTAL USED = '&
       &,total_used/1024/1024, 'MB i=',eye
#endif

  !----------------------------------------------------------------------------!
  ! runnit
  !----------------------------------------------------------------------------!
#ifdef __WITH_CRAYPAT__
  CALL PAT_HWPC_BEGIN(10, "williamson x 2")
#endif
  CALL wvic_tvdrk3(steps,info)
#ifdef __WITH_MEMORY_MEASUREMENT__
  eye = heap_info(fragments, total_free, largest_free, total_used)
  WRITE(6,*) '[after wvic_tvdrk3] heap_info fragments =',fragments,' total_free =', &
       & total_free/1024/1024,' largest_free = ',largest_free/1024/1024,' TOTAL USED = '&
       &,total_used/1024/1024, 'MB i=',eye
#endif
#ifdef __WITH_CRAYPAT__
  CALL PAT_HWPC_END(10)
  CALL PAT_HWPC_FINALIZE()
#endif
  !----------------------------------------------------------------------------!
  ! only mpi finalize because ppm_finalize dont work on xlf
  !----------------------------------------------------------------------------!
  CALL ppm_finalize(info)
  CALL MPI_Finalize(info)
  
END PROGRAM wvic
  


