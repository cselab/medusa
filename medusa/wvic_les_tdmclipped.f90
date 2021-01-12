!------------------------------------------------------------------------------!
!
!
!
!
!
!
!
!
!------------------------------------------------------------------------------!
!  project:  WVIC
!  purpose:  Clipped tensor diffusivity model (Cottet:1996)
!   author:  Michael Bergdorf,
!         :  ETHZ Computational Science & Engineering Laboratory
!    email:  bergdorf@inf.ethz.ch
!------------------------------------------------------------------------------!
!    $Log: wvic_les_tdmclipped.F,v $
!    Revision 1.2  2006/10/23 08:19:11  pchatela
!    LES models
!    Bugfixes in KE spectra and factor for Parseval identity
!    Removed the reset of noise in init_tvphysics_0 and_1
!
!    Revision 1.1.1.1  2006/07/25 15:13:47  menahel
!    initial import
!
!    Revision 1.1  2005/12/10 21:10:14  michaebe
!    implementation initiale
!
!------------------------------------------------------------------------------!

SUBROUTINE wvic_les_tdm

  USE module_wvic

  !----------------------------------------------------------------------------!
  !  Argumentation
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER                :: i,j,k,isub,isubl,info
  REAL(mk)               :: un,us,ve,vw,wt,wb
  REAL(mk)               :: uc,vc,wc
  INTEGER                :: lx,ly,lz,i2,j2,k2
  REAL(mk), DIMENSION(3) :: on,os,oe,ow,ot,ob,oc
  REAL(mk), DIMENSION(3) :: les, dw, dv
  REAL(mk)               :: ales, coef2
  REAL(mk)               :: dxi
    
  dxi = 1.0_mk/dx
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              dw = 0.0_mk
              DO lz=-1,1
                 k2=k+lz
                 DO ly=-1,1
                    j2 = j+ly
                    DO lx=-1,1
                       i2 = i+lx
                       dv(1) = (field_up(1,i2,j2,k2,isub)-field_up(1,i,j,k,isub))*REAL(lx,mk)*dxi
                       dv(2) = (field_up(2,i2,j2,k2,isub)-field_up(2,i,j,k,isub))*REAL(ly,mk)*dxi
                       dv(3) = (field_up(3,i2,j2,k2,isub)-field_up(3,i,j,k,isub))*REAL(lz,mk)*dxi
                       ales=ABS(dv(1)+dv(2)+dv(3))/9.
                       IF (dv(1)+dv(2)+dv(3).LT.0.) THEN
                          coef2 = les_tdm_C_comp
                       ELSE
                          coef2 = les_tdm_C_dil
                       END IF
                       dw(1) = dw(1) + (field_wp(1,i2,j2,k2,isub)-field_wp(1,i,j,k,isub))*ales*coef2
                       dw(2) = dw(2) + (field_wp(2,i2,j2,k2,isub)-field_wp(2,i,j,k,isub))*ales*coef2
                       dw(3) = dw(3) + (field_wp(3,i2,j2,k2,isub)-field_wp(3,i,j,k,isub))*ales*coef2
                    END DO
                 END DO
              END DO
              field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub) + dw(1)
              field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub) + dw(2)
              field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub) + dw(3)
              
           END DO

        END DO

     END DO

  END DO

END SUBROUTINE wvic_les_tdm

SUBROUTINE wvic_les_tdmclipped
  !----------------------------------------------------------------------------!
  !  Clipped LES from GSW
  !  Keep the positive part of 2 \nu + (u_p-u_q)\cdot (x_p-x_q)
  !----------------------------------------------------------------------------!
  USE module_wvic

  !----------------------------------------------------------------------------!
  !  Argumentation
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER                :: i,j,k,isub,isubl,info
  REAL(mk)               :: un,us,ve,vw,wt,wb
  REAL(mk)               :: uc,vc,wc
  INTEGER                :: lx,ly,lz,i2,j2,k2
  REAL(mk), DIMENSION(3) :: on,os,oe,ow,ot,ob,oc
  REAL(mk), DIMENSION(3) :: les, dw, dv
  REAL(mk)               :: ales, coef2
  REAL(mk)               :: dxi
  
  dxi = 1.0_mk/dx
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)
        DO j=1,ndata(2,isubl)
           DO i=1,ndata(1,isubl)
              dw = 0.0_mk
              DO lz=-1,1
                 k2=k+lz
                 DO ly=-1,1
                    j2 = j+ly
                    DO lx=-1,1
                       i2 = i+lx
                       dv(1) = (field_up(1,i2,j2,k2,isub)-field_up(1,i,j,k,isub))*REAL(lx,mk)*dx
                       dv(2) = (field_up(2,i2,j2,k2,isub)-field_up(2,i,j,k,isub))*REAL(ly,mk)*dx
                       dv(3) = (field_up(3,i2,j2,k2,isub)-field_up(3,i,j,k,isub))*REAL(lz,mk)*dx
                       IF (2*nu+dv(1)+dv(2)+dv(3).GT.0.0_mk) THEN
                          ales = dxi**2 * (dv(1)+dv(2)+dv(3)) * les_tdmclipped_C
                          dw(1) = dw(1) + (field_wp(1,i2,j2,k2,isub)-field_wp(1,i,j,k,isub))*ales
                          dw(2) = dw(2) + (field_wp(2,i2,j2,k2,isub)-field_wp(2,i,j,k,isub))*ales
                          dw(3) = dw(3) + (field_wp(3,i2,j2,k2,isub)-field_wp(3,i,j,k,isub))*ales
                       END IF
                    END DO
                 END DO
              END DO
              field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub) + dw(1)
              field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub) + dw(2)
              field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub) + dw(3)
              
           END DO

        END DO

     END DO

  END DO

END SUBROUTINE wvic_les_tdmclipped


SUBROUTINE wvic_les_tdmclipped_old

  USE module_wvic

  !----------------------------------------------------------------------------!
  !  Argumentation
  !----------------------------------------------------------------------------!
  
  !----------------------------------------------------------------------------!
  !  Localities
  !----------------------------------------------------------------------------!
  INTEGER                :: i,j,k,isub,isubl,info
  REAL(mk)               :: un,us,ve,vw,wt,wb
  REAL(mk)               :: uc,vc,wc
  REAL(mk), DIMENSION(3) :: on,os,oe,ow,ot,ob,oc
  REAL(mk), DIMENSION(3) :: les
  REAL(mk)               :: dxi

  dxi = 1.0_mk/dx
  
  DO isub=1,nsublist
     isubl = isublist(isub)
     DO k=1,ndata(3,isubl)

        DO j=1,ndata(2,isubl)

           DO i=1,ndata(1,isubl)

              !----------------------------------------------------------------!
              ! compute the velocity differeces
              !----------------------------------------------------------------!
              uc = field_up(1,i,j,k,isub)
              vc = field_up(2,i,j,k,isub)
              wc = field_up(3,i,j,k,isub)
              oc = field_wp(1:3,i,j,k,isub)
              
              un = uc - field_up(1,i+1,j,k,isub)
              us = uc - field_up(1,i-1,j,k,isub)

              ve = vc - field_up(2,i,j+1,k,isub)
              vw = vc - field_up(2,i,j-1,k,isub)

              wt = wc - field_up(3,i,j,k+1,isub)
              wb = wc - field_up(3,i,j,k-1,isub)

              !----------------------------------------------------------------!
              ! compute fluxes
              !----------------------------------------------------------------!
              on(1) = oc(1) - field_wp(1,i+1,j,k,isub)
              os(1) = oc(1) - field_wp(1,i-1,j,k,isub)
              oe(1) = oc(1) - field_wp(1,i,j+1,k,isub)
              ow(1) = oc(1) - field_wp(1,i,j-1,k,isub)
              ot(1) = oc(1) - field_wp(1,i,j,k+1,isub)
              ob(1) = oc(1) - field_wp(1,i,j,k-1,isub)

              on(2) = oc(2) - field_wp(2,i+1,j,k,isub)
              os(2) = oc(2) - field_wp(2,i-1,j,k,isub)
              oe(2) = oc(2) - field_wp(2,i,j+1,k,isub)
              ow(2) = oc(2) - field_wp(2,i,j-1,k,isub)
              ot(2) = oc(2) - field_wp(2,i,j,k+1,isub)
              ob(2) = oc(2) - field_wp(2,i,j,k-1,isub)

              on(3) = oc(3) - field_wp(3,i+1,j,k,isub)
              os(3) = oc(3) - field_wp(3,i-1,j,k,isub)
              oe(3) = oc(3) - field_wp(3,i,j+1,k,isub)
              ow(3) = oc(3) - field_wp(3,i,j-1,k,isub)
              ot(3) = oc(3) - field_wp(3,i,j,k+1,isub)
              ob(3) = oc(3) - field_wp(3,i,j,k-1,isub)

              !----------------------------------------------------------------!
              ! compute les term
              !----------------------------------------------------------------!
              les = 0.0_mk
              les(1) = les(1) + MAX( un,0.0_mk) * on(1)
              les(1) = les(1) + MAX(-us,0.0_mk) * os(1)
              les(1) = les(1) + MAX( ve,0.0_mk) * oe(1)
              les(1) = les(1) + MAX(-vw,0.0_mk) * ow(1)
              les(1) = les(1) + MAX( wt,0.0_mk) * ot(1)
              les(1) = les(1) + MAX(-wb,0.0_mk) * ob(1)

              les(2) = les(2) + MAX( un,0.0_mk) * on(2)
              les(2) = les(2) + MAX(-us,0.0_mk) * os(2)
              les(2) = les(2) + MAX( ve,0.0_mk) * oe(2)
              les(2) = les(2) + MAX(-vw,0.0_mk) * ow(2)
              les(2) = les(2) + MAX( wt,0.0_mk) * ot(2)
              les(2) = les(2) + MAX(-wb,0.0_mk) * ob(2)

              les(3) = les(3) + MAX( un,0.0_mk) * on(3)
              les(3) = les(3) + MAX(-us,0.0_mk) * os(3)
              les(3) = les(3) + MAX( ve,0.0_mk) * oe(3)
              les(3) = les(3) + MAX(-vw,0.0_mk) * ow(3)
              les(3) = les(3) + MAX( wt,0.0_mk) * ot(3)
              les(3) = les(3) + MAX(-wb,0.0_mk) * ob(3)

              field_dwp(1,i,j,k,isub) = field_dwp(1,i,j,k,isub) &
                   & + les_tdmclipped_C * les(1) * dxi
              field_dwp(2,i,j,k,isub) = field_dwp(2,i,j,k,isub) &
                   & + les_tdmclipped_C * les(2) * dxi
              field_dwp(3,i,j,k,isub) = field_dwp(3,i,j,k,isub) &
                   & + les_tdmclipped_C * les(3) * dxi
              
           END DO

        END DO

     END DO

  END DO

END SUBROUTINE wvic_les_tdmclipped_old


