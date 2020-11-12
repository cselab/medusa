FUNCTION sisec(s)
  USE module_wvic
  REAL(mk) :: s
  REAL(mk) :: phi
  REAL(mk) :: T
  REAL(mk) :: R
  REAL(mk) :: tx
  REAL(mk) :: ty
  REAL(mk) :: tz
  REAL(mk) :: sisec
  COMMON / varsisec / tx, ty, tz, phi, T, R
  sisec = - 2_mk * R * SIN((2_mk * M_PI * s) / T + phi ) &
       * M_PI * tx / T + &
       2_mk * R * COS((2_mk * M_PI * s) / T + phi ) &
       * M_PI * ty / T - s + tz;
END FUNCTION sisec
  
FUNCTION sisec2(s)
  USE module_wvic
  USE module_wvic
  REAL(mk) :: s
  REAL(mk) :: phi
  REAL(mk) :: T
  REAL(mk) :: R
  REAL(mk) :: tx
  REAL(mk) :: ty
  REAL(mk) :: tz
  REAL(mk) :: sisec2
  COMMON / varsisec / tx, ty, tz, phi, T, R
  sisec2 = -((4*R*COS((2*M_PI*s)/T+phi)*M_PI**2*tx+4*R*SIN((2*M_PI*s)/T+phi)*M_PI**2*ty)/T**2 + 1)
END FUNCTION sisec2

  
