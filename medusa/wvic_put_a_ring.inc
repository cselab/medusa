
  rpos = RING_CENTER(3)
  DO isub=1,nsublist
     isubl = isublist(isub)
     
     DO k=1,ndata(3,isubl)
        
        DO j=1,ndata(2,isubl)
           
           DO i=1,ndata(1,isubl)
              twp = 0.0_mk
              DO imirr=-1,1
                 DO jmirr=-1,1
                    DO kmirr=-1,1
!              imirr=0;jmirr=0;kmirr=0
                          ! locations of grid points
              tx = min_sub(1,isubl) + REAL(i-1,mk)*dx + &
                   &REAL(imirr,mk)*length(1)
              ty = min_sub(2,isubl) + REAL(j-1,mk)*dy + &
                   &REAL(jmirr,mk)*length(2)
              tz = min_sub(3,isubl) + REAL(k-1,mk)*dz + &
                   &REAL(kmirr,mk)*length(3)
              !-----------------------------------------------------------------
              ! calculate angle1, and radius1
              !-----------------------------------------------------------------
              theta1 = ATAN2((tx-RING_CENTER(1)),(ty-RING_CENTER(2)))
              !-----------------------------------------------------
              !  compute noise
              !-----------------------------------------------------
              omega_x = theta1
              noise   = 0.0_mk
              DO kk=1,nk
                 noise = noise &
                      & + SIN(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*aky(kk))) &
                      & + COS(REAL(kk,mk)*(omega_x+2.0_mk*M_PI*bky(kk)))
              END DO
              
              rad1r  = SQRT((tx-RING_CENTER(1))**2 + (ty-RING_CENTER(2))**2) -&
                   & (vrrad+noise*amp)
              rad1t  = tz-rpos
              rad1sq = rad1r**2 + rad1t**2
              radstrength = -EXP(-rad1sq*radstr)*radstr*gamma/M_PI
              twp(1) = twp(1) - radstrength*COS(theta1)
              twp(2) = twp(2) + radstrength*SIN(theta1)
              twp(3) = twp(3) + 0.0_mk
              !-----------------------------------------------------------------
              ! substract mirrored in order to get zero vorticity in
              ! the center
              !-----------------------------------------------------------------
              rad1r  = SQRT((tx-RING_CENTER(1))**2 + (ty-RING_CENTER(2))**2) +&
                   & (vrrad+noise*amp) 
              rad1t  = tz-rpos
              rad1sq = rad1r**2 + rad1t**2
              rad1sqTILDA=rad1sq*radstr
              radstrength =  EXP(-rad1sqTILDA)*radstr*gamma/M_PI
              !-----------------------------------------------------
              !  lets try without this
              !  Thu Jun 15 13:28:20 CEST 2006
              twp(1) = twp(1) - radstrength*COS(theta1)
              twp(2) = twp(2) + radstrength*SIN(theta1)
              
              END DO; END DO; END DO 
              field_wp(1,i,j,k,isub) = field_wp(1,i,j,k,isub) + twp(1)
              field_wp(2,i,j,k,isub) = field_wp(2,i,j,k,isub) + twp(2)
              field_wp(3,i,j,k,isub) = field_wp(3,i,j,k,isub) + twp(3)
              
           END DO

        END DO

     END DO

  END DO
