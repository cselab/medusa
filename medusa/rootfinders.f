  FUNCTION rtbis(func,x1,x2,xacc)
!  COMMON / varsisec1 / tx, ty, tz, phi, T, R, mk
  USE module_wvic
  IMPLICIT NONE
  REAL(mk), INTENT(IN) :: x1,x2,xacc
  REAL(mk) :: rtbis
  INTERFACE
      FUNCTION func(x)
      USE module_wvic
      IMPLICIT NONE
      REAL(mk), INTENT(IN) :: x
      REAL(mk)                :: func
      END FUNCTION func
  END INTERFACE
  INTEGER, PARAMETER :: MAXIT=100
  INTEGER  :: j
  REAL(mk) :: deltax,f,fmid,xmid
  fmid=func(x2)
  f=func(x1)
  if (f*fmid >= 0.0) then
    stop 'rtbis: root must be bracketed'
  end if
  if (f < 0.0) then                 !Orient the search so that f>0 lies at x+delta.
      rtbis=x1
      deltax=x2-x1
  else
      rtbis=x2
      deltax=x1-x2
  end if
  do j=1,MAXIT                      !Bisection loop.
      deltax=deltax*0.5_mk
      xmid=rtbis+deltax
      fmid=func(xmid)
      if (fmid <= 0.0) rtbis=xmid
      if (ABS(deltax) < xacc .or. fmid == 0.0) RETURN
  end do
  stop 'rtbis: too many bisections'
  END FUNCTION rtbis

FUNCTION rtflsp(func,x1,x2,xacc)
  USE module_wvic
  IMPLICIT NONE
  REAL(mk), INTENT(IN) :: x1,x2,xacc
  REAL(mk) :: rtflsp
  INTERFACE
    FUNCTION func(x)
      USE module_wvic
      IMPLICIT NONE
      REAL(mk), INTENT(IN) :: x
      REAL(mk) :: func
    END FUNCTION func
  END INTERFACE
  INTEGER, PARAMETER :: MAXIT=100
  INTEGER :: j
  REAL(mk) :: del,deltax,f,fh,fl,xh,xl,temp
  fl=func(x1)
  fh=func(x2)
  if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) then
    stop 'rtflsp: root must be bracketed between arguments'
  end if
  if (fl < 0.0) then
    xl=x1 
    xh=x2
  else
    xl=x2
    xh=x1
    temp=fl
    fl=fh
    fh=temp
  end if
  deltax=xh-xl
  do j=1,MAXIT
    rtflsp=xl+deltax*fl/(fl-fh)
    f=func(rtflsp)
    if (f < 0.0) then
    del=xl-rtflsp
    xl=rtflsp
    fl=f
  else
    del=xh-rtflsp
    xh=rtflsp
    fh=f
  end if
  deltax=xh-xl
  if (abs(del) < xacc .or. f == 0.0) RETURN
    end do
    stop 'rtflsp exceed maximum iterations'
END FUNCTION rtflsp


FUNCTION rtsec(func,x1,x2,xacc)
  USE module_wvic
  IMPLICIT NONE
  REAL(mk), INTENT(IN) :: x1,x2,xacc
  REAL(mk) :: rtsec
  INTERFACE
    FUNCTION func(x)
      USE module_wvic
      IMPLICIT NONE
      REAL(mk), INTENT(IN) :: x
      REAL(mk) :: func
    END FUNCTION func
  END INTERFACE
  INTEGER, PARAMETER :: MAXIT=100
  INTEGER :: j
  REAL(mk) :: deltax,f,fl,xl,temp
  fl=func(x1)
  f=func(x2)
  if (abs(fl) < abs(f)) then
    rtsec=x1 
    xl=x2
    temp=f
    f=fl
    fl=temp
  else
    xl=x1
    rtsec=x2
  end if
  do j=1,MAXIT
    deltax=(xl-rtsec)*f/(f-fl)
    xl=rtsec
    fl=f
    rtsec=rtsec+deltax
    f=func(rtsec)
    if (abs(deltax) < xacc .or. f == 0.0) RETURN
  end do
stop 'rtsec: exceed maximum iterations'
END FUNCTION rtsec


FUNCTION zriddr(func,x1,x2,xacc,test)
  USE module_wvic  
  IMPLICIT NONE
  REAL(mk), INTENT(IN) :: x1,x2,xacc
  INTEGER, INTENT(IN) :: test
  REAL(mk) :: zriddr
  INTERFACE
    FUNCTION func(x)
      USE module_wvic
      IMPLICIT NONE
      REAL(mk), INTENT(IN) :: x
      REAL(mk) :: func
    END FUNCTION func
  END INTERFACE
  INTEGER, PARAMETER :: MAXIT=100
  REAL(mk), PARAMETER :: UNUSED=-1.11e30_mk
  INTEGER :: j
  REAL(mk) :: fh,fl,fm,fnew,s,xh,xl,xm,xnew
  fl=func(x1)
  fh=func(x2)
  if ((fl > 0.0 .and. fh < 0.0) .or. (fl < 0.0 .and. fh > 0.0)) then
    xl=x1
    xh=x2
    zriddr=UNUSED
    do j=1,MAXIT
      xm=0.5_mk*(xl+xh)
      fm=func(xm)
      s=sqrt(fm**2-fl*fh)
      if (s == 0.0) RETURN
      xnew=xm+(xm-xl)*(sign(1.0_mk,fl-fh)*fm/s)
      if (abs(xnew-zriddr) <= xacc) RETURN
      zriddr=xnew
      fnew=func(zriddr)
      if (fnew == 0.0) RETURN
      if (sign(fm,fnew) /= fm) then
        xl=xm 
        fl=fm
        xh=zriddr
        fh=fnew
      else if (sign(fl,fnew) /= fl) then
        xh=zriddr
        fh=fnew
      else if (sign(fh,fnew) /= fh) then
        xl=zriddr
        fl=fnew
      else
        stop 'zriddr: never get here'
      end if
      if (abs(xh-xl) <= xacc) RETURN
    end do
    stop 'zriddr: exceeded maximum iterations'
  else if (fl == 0.0) then
    zriddr=x1
  else if (fh == 0.0) then
    zriddr=x2
  else
    if (test== 1) then
      stop 'zriddr: root must be bracketed, sisec1'
    else if (test == 2) then
      stop 'zriddr: root must be bracketed, sisec2a'
    else if (test == 3) then
      stop 'zriddr: root must be bracketed, sisec2b'
    else
      stop 'zriddr: root must be bracketed, ??'
    end if
  end if
END FUNCTION zriddr


FUNCTION zbrent(func,x1,x2,tol)
  USE module_wvic
  IMPLICIT NONE
  REAL(mk), INTENT(IN) :: x1,x2,tol
  REAL(mk) :: zbrent
  INTERFACE
    FUNCTION func(x)
      USE module_wvic
      IMPLICIT NONE
      REAL(mk), INTENT(IN) :: x
      REAL(mk) :: func
    END FUNCTION func
  END INTERFACE
  INTEGER, PARAMETER :: ITMAX=100
  REAL(mk), PARAMETER :: EPS=epsilon(x1)
  INTEGER :: iter
  REAL(mk) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  a=x1
  b=x2
  fa=func(a)
  fb=func(b)
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
    stop 'root must be bracketed for zbrent'
  end if
  c=b
  fc=fb
  do iter=1,ITMAX
    if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
      c=a 
      fc=fa 
      d=b-a
      e=d
    end if
    if (abs(fc) < abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    end if
    tol1=2.0_mk*EPS*abs(b)+0.5_mk*tol
    xm=0.5_mk*(c-b)
    if (abs(xm) <= tol1 .or. fb == 0.0) then
      zbrent=b
      RETURN
    end if
    if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
      s=fb/fa 
      if (a == c) then
        p=2.0_mk*xm*s
        q=1.0_mk-s
      else
        q=fa/fc
        r=fb/fc
        p=s*(2.0_mk*xm*q*(q-r)-(b-a)*(r-1.0_mk))
        q=(q-1.0_mk)*(r-1.0_mk)*(s-1.0_mk)
      end if
      if (p > 0.0) q=-q
      p=abs(p)
      if (2.0_mk*p < min(3.0_mk*xm*q-abs(tol1*q),abs(e*q))) then
        e=d
        d=p/q
      else
        d=xm
        e=d
      end if
    else
      d=xm
      e=d
    end if
    a=b
    fa=fb
    b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
    fb=func(b)
  end do
  stop 'zbrent: exceeded maximum iterations'
  zbrent=b
END FUNCTION zbrent










