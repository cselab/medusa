!-----------------------------------------------------------------------
!  Subroutine   :            ppm_mg_res 
!-----------------------------------------------------------------------
!  Purpose      : In this routine we compute the residula in each level
!            
!                  
!  
!  Input       :

!  Input/output :
! 
!  Output       : info        (I) return status. 0 upon success
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!-------------------------------------------------------------------------
!  $Log: ppm_mg_res_coarse.f,v $
!  Revision 1.1.1.1  2006/07/25 15:18:20  menahel
!  initial import
!
!  Revision 1.6  2006/02/08 19:56:24  kotsalie
!  fixed multiple domains
!
!  Revision 1.5  2005/12/08 12:44:45  kotsalie
!  commiting dirichlet
!
!  Revision 1.4  2005/03/14 13:21:51  kotsalie
!  COMMITED THE VECTOR CASE. IT IS FOR LDA=3
!
!  Revision 1.3  2004/09/28 14:07:40  kotsalie
!  Changes concerning 4th order
!
!  Revision 1.2  2004/09/23 12:16:49  kotsalie
!  Added USE statement
!
!  Revision 1.1  2004/09/22 18:47:32  kotsalie
!  MG new version
!
!
!------------------------------------------------------------------------  
!  Parallel Particle Mesh Library (PPM)
!  Institute of Computational Science
!  ETH Zentrum, Hirschengraben 84
!  CH-8092 Zurich, Switzerland
!------------------------------------------------------------------------- 

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_2D_sca_s(mlev,c1,c2,c3,c4,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_2D_sca_d(mlev,c1,c2,c3,c4,E,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_3D_sca_s(mlev,c1,c2,c3,c4,c5,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_3D_sca_d(mlev,c1,c2,c3,c4,c5,E,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_2D_vec_s(mlev,c1,c2,c3,c4,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_2D_vec_d(mlev,c1,c2,c3,c4,E,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_3D_vec_s(mlev,c1,c2,c3,c4,c5,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_coarse_3D_vec_d(mlev,c1,c2,c3,c4,c5,E,info)
#endif
#endif
#endif

        !---------------------------------------------------------------------- 
        !  Includes
        !----------------------------------------------------------------------
#include "ppm_define.h"

        !-------------------------------------------------------------------    
        !  Modules 
        !--------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_write
        USE ppm_module_data_mg
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_data_mesh



#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !-------------------------------------------------------------------    
        !  Arguments     
        !-------------------------------------------------------------------
        INTEGER,                   INTENT(IN)      ::  mlev
        REAL(MK),                  INTENT(OUT)     ::  E
#if  __MESH_DIM == __2D
        REAL(MK),                  INTENT(IN)      ::  c1,c2,c3,c4 
#elif __MESH_DIM == __3D
        REAL(MK),                  INTENT(IN)      ::  c1,c2,c3,c4,c5 
#endif
        INTEGER,                   INTENT(INOUT)   ::  info
        !---------------------------------------------------------------------  
        !  Local variables 
        !---------------------------------------------------------------------
        CHARACTER(LEN=256) :: cbuf
        INTEGER                                    ::  i,j,isub,color
        INTEGER                                    ::  ilda,isweep,count
        INTEGER                                    ::  aa,bb,cc,dd,ee,gg
        REAL(MK)                                   ::  c11,c22,c33,c44,c55 
        INTEGER                                    ::  k,idom
        REAL(MK)                                   ::  x,y
        REAL(MK)                                   ::  res
#if __MESH_DIM == __2D
        INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
        INTEGER,DIMENSION(3)                       ::  ldl3,ldu3
#endif
#if __MESH_DIM == __3D
        INTEGER,DIMENSION(5)                       ::  ldl5,ldu5
        INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
#endif
        INTEGER                                    ::  iopt,iface,topoid
        REAL(MK)                                   ::  t0
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_2d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_2d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_3d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_3d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_2d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_2d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        TYPE(mg_field_3d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
        TYPE(mg_field_3d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#endif

        !-----------------------------------------------------------------------
        !Externals
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !Initialize
        !-----------------------------------------------------------------------

        CALL substart('ppm_mg_res',t0,info)
        IF (l_print) THEN
         WRITE(cbuf,*) 'RESIDUAL in LEVEL:',mlev
         CALL PPM_WRITE(ppm_rank,'mg_res_coarse',cbuf,info)
        ENDIF

        !-----------------------------------------------------------------------
        !  Check arguments
        !-----------------------------------------------------------------------
        IF (ppm_debug .GT. 0) THEN
          IF (c1.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth',  &
     &            'Factor c1 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c2.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth',  &
     &            'Factor c2 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c3.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth',  &
     &            'Factor c3 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#if __MESH_DIM == __3D
          IF (c4.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_smooth',  &
     &            'Factor c4 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
        ENDIF
        !-----------------------------------------------------------------------
        !Definition of necessary variables and allocation of arrays
        !-----------------------------------------------------------------------
        topoid=ppm_field_topoid


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        mgfield=>mgfield_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
        mgfield=>mgfield_3d_vec_d
#endif
#endif
#endif


#if  __DIM == __SFIELD
#if  __MESH_DIM == __2D

        !-----------------------------------------------------------------------
        !Implementation
        !----------------------------------------------------------------------- 
        IF (order.EQ.ppm_param_order_2) THEN
                 DO isub=1,nsubs
                   aa=0
                   bb=0
                   cc=0
                   dd=0

                IF (.NOT.lperiodic) THEN
                 DO iface=1,4
                    IF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                   !DO NOTHING
                    ELSEIF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                        aa=1
                      ELSEIF (iface.EQ.2) THEN
                        bb=1
                      ELSEIF (iface.EQ.3) THEN
                       cc=1
                      ELSEIF (iface.EQ.4) THEN
                       dd=1
		             ENDIF
                     ENDIF 
                    ENDDO !iface
	              endif !periodic
             ENDDO
         E=-HUGE(E)
          DO isub=1,nsubs
            DO j=start(2,isub,mlev)+cc,stop(2,isub,mlev)-dd
               DO i=start(1,isub,mlev)+aa,stop(1,isub,mlev)-bb
                     res =(mgfield(isub,mlev)%uc(i-1,j)+&
     &                     mgfield(isub,mlev)%uc(i+1,j))*c2 + &
     &                    (mgfield(isub,mlev)%uc(i,j-1)+    &
     &                     mgfield(isub,mlev)%uc(i,j+1))*c3 - &
     &                     mgfield(isub,mlev)%uc(i,j)*c4 - &
     &                     mgfield(isub,mlev)%fc(i,j)
                     E=MAX(ABS(res),E)
                     mgfield(isub,mlev)%err(i,j)=-res
               ENDDO
            ENDDO
          ENDDO

        ELSEIF (order.EQ.ppm_param_order_4) THEN  
  

        c22=c2/12.0_MK
        c33=c3/12.0_MK
        c44=c4*1.25_MK

         E=-HUGE(E)
          DO isub=1,nsubs
            DO j=start(2,isub,mlev),stop(2,isub,mlev)
               DO i=start(1,isub,mlev),stop(1,isub,mlev)
                     res =(16.0_MK*mgfield(isub,mlev)%uc(i-1,j)+&
     &                     16.0_MK*mgfield(isub,mlev)%uc(i+1,j)-&
     &                     mgfield(isub,mlev)%uc(i-2,j)-&
     &                     mgfield(isub,mlev)%uc(i+2,j))*c22 + &
     &                    (16.0_MK*mgfield(isub,mlev)%uc(i,j-1)+    &
     &                     16.0_MK*mgfield(isub,mlev)%uc(i,j+1)-&
     &                     mgfield(isub,mlev)%uc(i,j-2)-&
     &                     mgfield(isub,mlev)%uc(i,j+2))*c33 - &
     &                     mgfield(isub,mlev)%uc(i,j)*c44 - &
     &                     mgfield(isub,mlev)%fc(i,j)
                     E=MAX(ABS(res),E)
                     mgfield(isub,mlev)%err(i,j)=-res
               ENDDO
            ENDDO
          ENDDO



        ENDIF 

#elif __MESH_DIM == __3D

                 DO isub=1,nsubs
                   aa=0
                   bb=0
                   cc=0
                   dd=0
                   ee=0
                   gg=0

                IF (.NOT.lperiodic) THEN
                 DO iface=1,6
                    IF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                   !DO NOTHING
                    ELSEIF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                        aa=1
                      ELSEIF (iface.EQ.2) THEN
                        bb=1
                      ELSEIF (iface.EQ.3) THEN
                       cc=1
                      ELSEIF (iface.EQ.4) THEN
                       dd=1
                      ELSEIF (iface.EQ.5) Then
                       ee=1
                      ELSEIF (iface.EQ.6) Then
                       gg=1
		             ENDIF
                     ENDIF 
                    ENDDO !iface
	              endif !periodic
             ENDDO
        !-----------------------------------------------------------------------
        !Implementation
        !----------------------------------------------------------------------- 
        E=-HUGE(E)
          DO isub=1,nsubs
           DO k=start(3,isub,mlev)+ee,stop(3,isub,mlev)-gg
              DO j=start(2,isub,mlev)+cc,stop(2,isub,mlev)-dd
                 DO i=start(1,isub,mlev)+aa,stop(1,isub,mlev)-bb
                       res =(mgfield(isub,mlev)%uc(i-1,j,k)+&
     &                       mgfield(isub,mlev)%uc(i+1,j,k))*c2 + &
     &                      (mgfield(isub,mlev)%uc(i,j-1,k)+    &
     &                       mgfield(isub,mlev)%uc(i,j+1,k))*c3 +&
     &                      (mgfield(isub,mlev)%uc(i,j,k-1)+    &
     &                       mgfield(isub,mlev)%uc(i,j,k+1))*c4 -&
     &                       mgfield(isub,mlev)%uc(i,j,k)*c5 - &
     &                       mgfield(isub,mlev)%fc(i,j,k)
                       E=MAX(ABS(res),E)
                       mgfield(isub,mlev)%err(i,j,k)=-res
                 ENDDO
              ENDDO
           ENDDO
        ENDDO




#endif
#elif __DIM == __VFIELD
#if  __MESH_DIM == __2D

                 DO isub=1,nsubs
                   DO ilda=1,vecdim
                   aa=0
                   bb=0
                   cc=0
                   dd=0

                IF (.NOT.lperiodic) THEN
                 DO iface=1,4
                    IF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                   !DO NOTHING
                    ELSEIF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                        aa=1
                      ELSEIF (iface.EQ.2) THEN
                        bb=1
                      ELSEIF (iface.EQ.3) THEN
                       cc=1
                      ELSEIF (iface.EQ.4) THEN
                       dd=1
		             ENDIF
                     ENDIF 
                    ENDDO !iface
	              endif !periodic
                 ENDDO
             ENDDO
        !-----------------------------------------------------------------------
        !Implementation
        !----------------------------------------------------------------------- 

                  E=-HUGE(E)
        DO isub=1,nsubs
           DO j=start(2,isub,mlev)+cc,stop(2,isub,mlev)-dd
              DO i=start(1,isub,mlev)+aa,stop(1,isub,mlev)-bb
               DO ilda=1,vecdim
                    res =(mgfield(isub,mlev)%uc(ilda,i-1,j)+&
     &                    mgfield(isub,mlev)%uc(ilda,i+1,j))*c2 + &
     &                   (mgfield(isub,mlev)%uc(ilda,i,j-1)+    &
     &                    mgfield(isub,mlev)%uc(ilda,i,j+1))*c3 - &
     &                    mgfield(isub,mlev)%uc(ilda,i,j)*c4 - &
     &                    mgfield(isub,mlev)%fc(ilda,i,j)
                    E=MAX(ABS(res),E)
                    mgfield(isub,mlev)%err(ilda,i,j)=-res
               ENDDO
              ENDDO
           ENDDO
        ENDDO



#elif __MESH_DIM == __3D

        !-----------------------------------------------------------------------
        !Implementation
        !----------------------------------------------------------------------- 
         
     IF (order.EQ.ppm_param_order_2) THEN

                 DO isub=1,nsubs
                   aa=0
                   bb=0
                   cc=0
                   dd=0
                   ee=0
                   gg=0
                   DO ilda=1,vecdim

                IF (.NOT.lperiodic) THEN
                 DO iface=1,6
                    IF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                   !DO NOTHING
                    ELSEIF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                        aa=1
                      ELSEIF (iface.EQ.2) THEN
                        bb=1
                      ELSEIF (iface.EQ.3) THEN
                       cc=1
                      ELSEIF (iface.EQ.4) THEN
                       dd=1
                      ELSEIF (iface.EQ.5) Then
                       ee=1
                      ELSEIF (iface.EQ.6) Then
                       gg=1
		             ENDIF
                     ENDIF 
                    ENDDO !iface
	              endif !periodic
                 ENDDO
             !ENDDO
 
            E=-HUGE(E)
            !DO isub=1,nsubs
             DO k=start(3,isub,mlev)+ee,stop(3,isub,mlev)-gg
               DO j=start(2,isub,mlev)+cc,stop(2,isub,mlev)-dd
                  DO i=start(1,isub,mlev)+aa,stop(1,isub,mlev)-bb
#ifdef __VECTOR
                        res =(mgfield(isub,mlev)%uc(1,i-1,j,k)+&
     &                        mgfield(isub,mlev)%uc(1,i+1,j,k))*c2 + &
     &                       (mgfield(isub,mlev)%uc(1,i,j-1,k)+    &
     &                        mgfield(isub,mlev)%uc(1,i,j+1,k))*c3 +&
     &                       (mgfield(isub,mlev)%uc(1,i,j,k-1)+    &
     &                        mgfield(isub,mlev)%uc(1,i,j,k+1))*c4 -&
     &                        mgfield(isub,mlev)%uc(1,i,j,k)*c5 - &
     &                        mgfield(isub,mlev)%fc(1,i,j,k)
                        E=MAX(ABS(res),E)
                        mgfield(isub,mlev)%err(1,i,j,k)=-res

                        res =(mgfield(isub,mlev)%uc(2,i-1,j,k)+&
     &                        mgfield(isub,mlev)%uc(2,i+1,j,k))*c2 + &
     &                       (mgfield(isub,mlev)%uc(2,i,j-1,k)+    &
     &                        mgfield(isub,mlev)%uc(2,i,j+1,k))*c3 +&
     &                       (mgfield(isub,mlev)%uc(2,i,j,k-1)+    &
     &                        mgfield(isub,mlev)%uc(2,i,j,k+1))*c4 -&
     &                        mgfield(isub,mlev)%uc(2,i,j,k)*c5 - &
     &                        mgfield(isub,mlev)%fc(2,i,j,k)
                        E=MAX(ABS(res),E)
                        mgfield(isub,mlev)%err(2,i,j,k)=-res

                        res =(mgfield(isub,mlev)%uc(3,i-1,j,k)+&
     &                        mgfield(isub,mlev)%uc(3,i+1,j,k))*c2 + &
     &                       (mgfield(isub,mlev)%uc(3,i,j-1,k)+    &
     &                        mgfield(isub,mlev)%uc(3,i,j+1,k))*c3 +&
     &                       (mgfield(isub,mlev)%uc(3,i,j,k-1)+    &
     &                        mgfield(isub,mlev)%uc(3,i,j,k+1))*c4 -&
     &                        mgfield(isub,mlev)%uc(3,i,j,k)*c5 - &
     &                        mgfield(isub,mlev)%fc(3,i,j,k)
                        E=MAX(ABS(res),E)
                        mgfield(isub,mlev)%err(3,i,j,k)=-res

#else
                   DO ilda=1,vecdim
                        res =(mgfield(isub,mlev)%uc(ilda,i-1,j,k)+&
     &                        mgfield(isub,mlev)%uc(ilda,i+1,j,k))*c2 + &
     &                       (mgfield(isub,mlev)%uc(ilda,i,j-1,k)+    &
     &                        mgfield(isub,mlev)%uc(ilda,i,j+1,k))*c3 +&
     &                       (mgfield(isub,mlev)%uc(ilda,i,j,k-1)+    &
     &                        mgfield(isub,mlev)%uc(ilda,i,j,k+1))*c4 -&
     &                        mgfield(isub,mlev)%uc(ilda,i,j,k)*c5 - &
     &                        mgfield(isub,mlev)%fc(ilda,i,j,k)
                        E=MAX(ABS(res),E)
                        mgfield(isub,mlev)%err(ilda,i,j,k)=-res
                   ENDDO
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
     ELSEIF (order.EQ.ppm_param_order_4) THEN


     ENDIF

#endif
#endif


        !---------------------------------------------------------------------- 
        !  Return 
        !----------------------------------------------------------------------
9999    CONTINUE
        CALL substop('ppm_mg_res',t0,info)
        RETURN
#if __DIM == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_2D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_2D_sca_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_3D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_3D_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_2D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_2D_vec_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_3D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_coarse_3D_vec_d
#endif
#endif
#endif




