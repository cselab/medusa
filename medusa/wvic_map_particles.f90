!-------------------------------------------------------------------------------
!* filename: wvic_map_particles                                              *!
!* project : ppm                                                              *!
!* purpose : maps the particles                                               *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 10 20:06:20 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!-------------------------------------------------------------------------------
! $Log: wvic_map_particles.F,v $
! Revision 1.1.1.1  2006/07/25 15:13:47  menahel
! initial import
!
! Revision 1.1  2005/09/28 11:40:38  michaebe
! Fork from ppm_pvc
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Map the particles
!-------------------------------------------------------------------------------
SUBROUTINE wvic_map_particles (pmaptype,info)

  USE module_wvic
  USE ppm_module_map
  USE ppm_module_data
  USE ppm_module_write
  USE MPI

  !-----------------------------------------------------------------------------
  !  arguments
  !-----------------------------------------------------------------------------
  INTEGER, INTENT(in) :: pmaptype  ! primary maptype, shoudl be either
                                   ! partial or global
  INTEGER, INTENT(inout) :: info

  
  !-----------------------------------------------------------------------------
  !  locals
  !-----------------------------------------------------------------------------
  INTEGER             :: maptype, mpart
  CHARACTER(len=256)  :: msg
  REAL(mk)            :: cpu1, cpu0 
  !-----------------------------------------------------------------------------
  !  make sure the particles are all in the domain
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------
  !  memory measurements
  !-----------------------------------------------------
  INTEGER    :: eye
  INTEGER*4  :: fragments
  INTEGER*8  :: total_free, largest_free, total_used
  INTEGER    :: heap_info


  CALL wvic_pbc(info)
  


#ifndef __SMALL_FOOTPRINT__
     maptype = ppm_param_map_partial
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     ! ===   push   ====
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(wp0,lda,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(dwp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(wp,lda,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#else
     maptype = ppm_param_map_partial
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(up,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(wp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_push
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_send
     CALL ppm_map_part(xp,dime,np,mpart,topo_id,maptype,info)
     maptype = ppm_param_map_pop
     CALL ppm_map_part(dwp,dime,np,mpart,topo_id,maptype,info)
     np = mpart
#endif

  !-----------------------------------------------------------------------------
  !  done
  !-----------------------------------------------------------------------------

END SUBROUTINE wvic_map_particles
  
