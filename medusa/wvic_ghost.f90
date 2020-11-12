!-------------------------------------------------------------------------------
!* filename: wvic_ghost                                                      *!
!* project : ppm                                                              *!
!* purpose : get ghost values                                                 *!
!*         :                                                                  *!
!* author  : Michael Bergdorf                                                 *!
!*         : Computational Science and Engineering Lab (CSE-Lab)              *!
!*         : ICOS, ETH Zurich                                                 *!
!*         :                                                                  *!
!* date    : Tue Aug 17 14:12:03 2004                                         *!
!* please return to <bergdorf@inf.ethz.ch> currently at the ['si-sE 'lab]     *!
!
!  $Log: wvic_ghost.F,v $
!  Revision 1.4  2006/09/27 09:30:21  pchatela
!  Fixes, spectra calculation,
!  most importantly: moved the u_infty out, so it does not kill the dgammadt
!
!  Revision 1.3  2006/09/16 00:22:03  pchatela
!  Implemented the kinetic energy spectrum, dumped into an ascii file.
!
!  Revision 1.2  2006/09/11 14:57:27  pchatela
!  Fixed velocity computation with odd/even symmetries
!  Added flag for adaptive time step
!
!  Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!  initial import
!
!  Revision 1.2  2005/11/11 14:04:24  michaebe
!  clean up, additions, comments
!
!  Revision 1.1  2005/09/28 11:40:35  michaebe
!  Fork from ppm_pvc
!
!-------------------------------------------------------------------------------




SUBROUTINE wvic_ghost(which, info)

  USE module_wvic
  USE ppm_module_write
  USE ppm_module_data
  USE ppm_module_map
  

  INTEGER , INTENT(in)    :: which
  INTEGER , INTENT(inout)   :: info
  
  INTEGER :: maptype,dlda
  INTEGER            :: isub, isubl
  INTEGER            :: i,j,k,kl,comp
  LOGICAL, DIMENSION(3) :: evenlower, evenupper
  LOGICAL, SAVE :: firsttime = .TRUE.
  REAL(mk), DIMENSION(:,:,:,:,:), POINTER :: who
  REAL(mk), DIMENSION(3)                  :: offsets
  CHARACTER(len=256)      :: msg 

  IF(verbose) WRITE(msg,*) 'is ghosting'
  IF(Verbose) CALL ppm_write(rank,'wvic_ghost',msg,info)
  
  SELECT CASE(which)

  CASE(wvic_prm_vorticity)
     who => field_wp
     dlda = lda
  CASE(wvic_prm_dgammadt)
     who => field_dwp
     dlda = lda
  CASE(wvic_prm_stream)
     who => field_wps
     dlda = dime
  CASE(wvic_prm_velocity)
     who => field_up
     dlda = dime
  END SELECT

  IF(firsttime) THEN
     firsttime = .FALSE.
     maptype = ppm_param_map_init
     
     CALL ppm_map_field_ghost(who,dlda,topo_id,mesh_id,&
          & ghostsize,maptype,info)
     
  END IF
  
  maptype = ppm_param_map_ghost_get
  CALL ppm_map_field_ghost(who,dlda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  
  maptype = ppm_param_map_push
  CALL ppm_map_field_ghost(who,dlda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  
  maptype = ppm_param_map_send
  CALL ppm_map_field_ghost(who,dlda,topo_id,mesh_id, &
       & ghostsize,maptype,info)
  
  maptype = ppm_param_map_pop
  CALL ppm_map_field_ghost(who,dlda,topo_id,mesh_id, &
       & ghostsize,maptype,info)

  IF (trailvortex) THEN
      SELECT CASE(which)
      CASE(wvic_prm_vorticity)
      DO isub=1,nsublist
         isubl = isublist(isub)
         IF (min_sub(3,isubl).LT.(min_physg(3)+0.5*dz)) THEN
            DO k=1-ghostsize(3),0
               DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                  DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                     who(1,i,j,k,isub) = who(1,i,j,2-k,isub)
                     who(2,i,j,k,isub) = who(2,i,j,2-k,isub)
                     who(3,i,j,k,isub) = who(3,i,j,2-k,isub)
                     !who(1,i,j,k,isub) = who(1,i,j,1,isub) - &
                     !& ( who(1,i,j,2-k,isub) - who(1,i,j,1,isub) )
                     !who(2,i,j,k,isub) = who(2,i,j,1,isub) - &
                     !& ( who(2,i,j,2-k,isub) - who(2,i,j,1,isub) )
                     !who(3,i,j,k,isub) = who(3,i,j,1,isub) - &
                     !& ( who(3,i,j,2-k,isub) - who(3,i,j,1,isub) )
                  END DO
               END DO
            END DO
         END IF
         IF (max_sub(3,isubl).GT.(max_physg(3)-0.5*dz)) THEN
            DO kl=1,ghostsize(3)
               k = ndata(3,isubl)+kl
               DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                  DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                     who(1,i,j,k,isub) = who(1,i,j,ndata(3,isubl),isub)
                     who(2,i,j,k,isub) = who(2,i,j,ndata(3,isubl),isub)
                     who(3,i,j,k,isub) = who(3,i,j,ndata(3,isubl),isub)
                  END DO
               END DO
            END DO
         END IF
     END DO
     CASE(wvic_prm_dgammadt)
     DO isub=1,nsublist
         isubl = isublist(isub)
         IF (min_sub(3,isubl).LT.(min_physg(3)+0.5*dz)) THEN
            DO k=1-ghostsize(3),0
               DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                  DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                     who(1,i,j,k,isub) = who(1,i,j,1,isub)
                     who(2,i,j,k,isub) = who(2,i,j,1,isub)
                     who(3,i,j,k,isub) = who(3,i,j,1,isub)
                  END DO
               END DO
            END DO
         END IF
         IF (max_sub(3,isubl).GT.(max_physg(3)-0.5*dz)) THEN
            DO kl=1,ghostsize(3)
               k = ndata(3,isubl)+kl
               DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                  DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                     who(1,i,j,k,isub) = 0.0_mk
                     who(2,i,j,k,isub) = 0.0_mk
                     who(3,i,j,k,isub) = 0.0_mk
                  END DO
               END DO
            END DO
         END IF
     END DO
     CASE(wvic_prm_velocity)
     DO isub=1,nsublist
         isubl = isublist(isub)
         IF (min_sub(3,isubl).LT.(min_physg(3)+0.5*dz)) THEN
            DO k=1-ghostsize(3),0
               DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                  DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                     who(1,i,j,k,isub) = who(1,i,j,2-k,isub)
                     who(2,i,j,k,isub) = who(2,i,j,2-k,isub)
                     who(3,i,j,k,isub) = who(3,i,j,2-k,isub)
                     !who(1,i,j,k,isub) = who(1,i,j,1,isub) - &
                     !& ( who(1,i,j,2-k,isub) - who(1,i,j,1,isub) )
                     !who(2,i,j,k,isub) = who(2,i,j,1,isub) - &
                     !& ( who(2,i,j,2-k,isub) - who(2,i,j,1,isub) )
                     !who(3,i,j,k,isub) = who(3,i,j,1,isub) - &
                     !& ( who(3,i,j,2-k,isub) - who(3,i,j,1,isub) )
                  END DO
               END DO
            END DO
         END IF
         IF (max_sub(3,isubl).GT.(max_physg(3)-0.5*dz)) THEN
            DO kl=1,ghostsize(3)
               k = ndata(3,isubl)+kl
               DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                  DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                     who(1,i,j,k,isub) = who(1,i,j,ndata(3,isubl),isub)
                     who(2,i,j,k,isub) = who(2,i,j,ndata(3,isubl),isub)
                     who(3,i,j,k,isub) = who(3,i,j,ndata(3,isubl),isub)
                  END DO
               END DO
            END DO
         END IF
     END DO
     END SELECT
  END IF
  
#ifdef __ODDEVENGHOSTS
  ! Now we extend odd/even fields...
  IF (trailvortex) THEN
      SELECT CASE(which)
      CASE(wvic_prm_vorticity)
      evenlower(1) = .FALSE.
      evenlower(2) = .FALSE.
      evenlower(3) = .TRUE.
      evenupper(1) = .FALSE.
      evenupper(2) = .FALSE.
      evenupper(3) = .TRUE.
      offsets = 0.0_mk
      
      CASE(wvic_prm_dgammadt)
      evenlower(1) = .FALSE.
      evenlower(2) = .FALSE.
      evenlower(3) = .TRUE.
      evenupper(1) = .FALSE.
      evenupper(2) = .FALSE.
      evenupper(3) = .TRUE.
      offsets = 0.0_mk
      
      CASE(wvic_prm_stream)
      evenlower(1) = .FALSE.
      evenlower(2) = .FALSE.
      evenlower(3) = .TRUE.
      evenupper(1) = .FALSE.
      evenupper(2) = .FALSE.
      evenupper(3) = .TRUE.
      offsets = 0.0_mk
      
      CASE(wvic_prm_velocity)
      evenlower(1) = .TRUE.
      evenlower(2) = .TRUE.
      evenlower(3) = .FALSE.
      evenupper(1) = .TRUE.
      evenupper(2) = .TRUE.
      evenupper(3) = .FALSE.
      offsets(1) = u_infty(1)
      offsets(2) = u_infty(2)
      offsets(3) = u_infty(3)
      
      END SELECT

      DO isub=1,nsublist
         isubl = isublist(isub)
         IF (min_sub(3,isubl).LT.(min_physg(3)+0.5*dz)) THEN
            DO comp=1,3
               IF (evenlower(comp)) THEN
                  DO k=1-ghostsize(3),0
                      DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                          DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                             who(comp,i,j,k,isub) = who(comp,i,j,2-k,isub)
                          END DO
                      END DO
                   END DO
               ELSE
                  DO k=1-ghostsize(3),0
                     DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                        DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                            who(comp,i,j,k,isub) = 2.0_mk*offsets(comp) - who(comp,i,j,2-k,isub)
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END IF
         IF (max_sub(3,isubl).GT.(max_physg(3)-0.5*dz)) THEN
            DO comp=1,3
               IF (evenupper(comp)) THEN
                  DO kl=1,ghostsize(3)
                     k = ndata(3,isubl)+kl
                     DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                        DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                           who(comp,i,j,k,isub) = who(comp,i,j,ndata(3,isubl)-kl,isub)
                        END DO
                     END DO
                  END DO
               ELSE
                  DO kl=1,ghostsize(3)
                     k = ndata(3,isubl)+kl
                     DO j=1-ghostsize(2),ndata(2,isubl)+ghostsize(2)
                        DO i=1-ghostsize(1),ndata(1,isubl)+ghostsize(1)
                           who(comp,i,j,k,isub) = 2.0_mk*offsets(comp) - who(comp,i,j,ndata(3,isubl)-kl,isub)
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END IF
     END DO
  END IF
#endif
  
  IF(verbose) WRITE(msg,*) 'done'
  IF(Verbose) CALL ppm_write(rank,'wvic_ghost',msg,info)
  

END SUBROUTINE wvic_ghost
