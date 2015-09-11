!
!  program to test the laplace inversion of a diffusion equation
!
PROGRAM fracd

  IMPLICIT NONE

  INCLUDE "fftw3.f"  ! FFTW constants

  INTEGER :: n,m,kx
  INTEGER, PARAMETER :: nt = 100, nx = 128, nkx = nx/2 + 1
  REAL, PARAMETER :: pi = 4.d0*ATAN(1.d0)
  COMPLEX(kind=8), PARAMETER :: czero = CMPLX(0.d0,0.d0)
  REAL(kind=8) :: std,beta,tmin,tmax,dt,xmax,xmin,dx,x0,theta,ytnorm
  REAL(kind=8), DIMENSION(nt) :: t
  REAL(kind=8), DIMENSION(nx) :: x
  REAL(kind=8), DIMENSION(nx) :: ya, yl, yt
  REAL(kind=8), DIMENSION(nx,1) :: yt2d
  REAL(kind=8), DIMENSION(nx,nt) :: yanal, ylap, ytemp, rerr, aerr
  COMPLEX(kind=8) :: xphase
  COMPLEX(kind=8), DIMENSION(0:nkx-1) :: ylc, ytc
  COMPLEX(kind=8), DIMENSION(0:nkx-1,1) :: ytc2d

  INTEGER(kind=8) :: planb
  REAL(kind=8), EXTERNAL :: gaussian
  COMPLEX(kind=8), EXTERNAL :: talbot, cmult

  ! initialize variables
  tmin = 0.0; tmax = 1.0
  xmin = -5.0; xmax = 5.0
  x0 = 0.5d0

  std = 1.0
  beta = 1.0

  dt = (tmax-tmin)/REAL(nt)    ! avoid t = 0.0 (tmin,tmax]
  dx = (xmax-xmin)/REAL(nx-1)  ! inclusive [xmin,xmax]

  t = (/ (tmin + REAL(n)*dt, n = 1, nt) /)
  x = (/ (xmin + REAL(n-1)*dx, n = 1, nx) /)

  ! create the analytic solution to classical diffusion
  OPEN(42,file="ya.dat",status="replace",form="formatted")
  DO m = 1, nt
     DO n = 1, nx
        ya(n) = gaussian(x(n),t(m),std,beta)
     END DO
     yanal(:,m) = ya
     WRITE(42,*) (ya(n), n = 1, nx)

  END DO
  CLOSE(42)

  ! test FFTW3 - transform a Gaussian function from k-space with exponential decay in time
  ytc = czero; yt = 0.0
  ytc2d = czero; yt2d = czero
  ytnorm = (2.d0*pi)
  CALL dfftw_plan_dft_c2r_1d(planb,nx,ytc,yt,FFTW_ESTIMATE)
  OPEN(44,file="yt.dat",status="replace",form="formatted")
  DO m = 1, nt
     DO kx = 0, nkx-1
        xphase = CMPLX(COS(2.0d0*pi*REAL(kx)*x0),-1.0d0*SIN(2.0d0*pi*REAL(kx)*x0))
        ytc(kx) = xphase*EXP(-REAL(kx*kx)*t(m)/ytnorm)! generate test data
     END DO
     ytc = ytc/(4.d0*pi)
     CALL dfftw_execute_dft_c2r(planb,ytc,yt)
     WRITE(44,*) (yt(n), n = 1, nx)
     ytemp(:,m) = yt
  END DO
  CLOSE(44)
  CALL dfftw_destroy_plan(planb)



  ! calculate the transformed data
  ! begin FFT section
  ylc = czero; yl = 0.0

  CALL dfftw_plan_dft_c2r_1d(planb,nx,ylc,yl,FFTW_ESTIMATE)

  OPEN(43,file="yl.dat",status="replace",form="formatted")
  DO m = 1, nt
     DO kx = 0, nkx-1
        theta = 2.0d0*pi*REAL(kx)*x0
        xphase = CMPLX(COS(theta),-SIN(theta))
        ylc(kx) = talbot(t(m),20,0.d0,REAL(kx,kind=8),beta)  ! calculate the laplace transform in time
        ylc(kx) = cmult(xphase,ylc(kx))
     END DO

     CALL dfftw_execute_dft_c2r(planb,ylc,yl)
!     yl = yl/REAL(nkx)
     yl = yl/(4.d0*pi)
     ylap(:,m) = yl
     WRITE(43,*) (yl(n), n = 1, nx)

  END DO
  CLOSE(43)

  CALL dfftw_destroy_plan(planb)


  ! find the absolute error
  OPEN(45,file="aerr.dat",status="replace",form="formatted")
  OPEN(46,file="rerr.dat",status="replace",form="formatted")
  OPEN(47,file="maxaerr.dat",status="replace",form="formatted")
  DO m = 1, nt
     DO n = 1, nx
        CALL find_err(yanal(n,m),ylap(n,m),aerr(n,m),rerr(n,m))
     END DO
     WRITE(45,*) (aerr(n,m), n = 1, nx)
     WRITE(46,*) (rerr(n,m), n = 1, nx)
     WRITE(47,*) t(m), MAXVAL(aerr(:,m)), MINVAL(aerr(:,m))
  END DO
  CLOSE(45); CLOSE(46); CLOSE(47)

END PROGRAM fracd





SUBROUTINE find_err(ya,ye,aerr,rerr)
!
!  calculates the absolute relative error weighted by actual value "ya"
!
  IMPLICIT NONE

  REAL(kind=8), INTENT(IN) :: ya,ye
  REAL(kind=8), INTENT(OUT) :: aerr,rerr
  REAL(kind=8) :: yt

  IF (ya .EQ. 0.0) THEN
     yt = 1.d-5
  ELSE
     yt = ya
  END IF

  aerr = ABS(ya - ye)
  rerr = aerr/yt
  

END SUBROUTINE find_err


REAL(kind=8) FUNCTION gaussian(x,t,std,beta) RESULT(y)
!
!  function for rescaled Gaussian
!
  IMPLICIT NONE

  REAL(kind=8), PARAMETER :: pi = 4.d0*ATAN(1.d0)

  REAL(kind=8), INTENT(IN) :: x,t,std,beta
  REAL(kind=8) :: norm

  norm = 1.d0/SQRT(2.d0*pi*std**2*t**beta)
  y = norm*EXP(-x*x/(2.d0*std*std*t))

END FUNCTION gaussian


COMPLEX(kind=8) FUNCTION talbot(t,ns,shift,k,beta) RESULT(fsval)
!
! port from Python script - Talbot quadratures with optimized coefficients for contour integration
!
  IMPLICIT NONE

  REAL(kind=8), PARAMETER :: pi = 4.d0*ATAN(1.d0), eps0 = 1.d-8
  REAL(kind=8), PARAMETER ::  & ! constants for contour integration
       c1 = 0.5017d0, &
       c2 = 0.6407d0, &
       c3 = 0.6122d0, &
       c4 = 0.2645d0
  
  COMPLEX(kind=8), PARAMETER :: ic = CMPLX(0.d0,1.d0), czero = CMPLX(0.d0,0.d0)

  INTEGER :: n
  INTEGER, INTENT(IN) :: ns
  REAL(kind=8), INTENT(IN) :: t, shift  ! shift on the REAL frequency axis from poles
  REAL(kind=8), INTENT(IN) :: k, beta
  REAL(kind=8) :: h, theta
  COMPLEX(kind=8) :: z,dz,tval,cfsval,texp

  COMPLEX(kind=8), EXTERNAL :: sfunc, cmult, cdiv

  h = 2.d0*pi/REAL(ns)  ! step-size

  ! check for t = 0.0
  IF (t <= 0.0) THEN
     WRITE(*,*) "Cannot invert at t <= 0."
     RETURN
  END IF

  fsval = 0.0
  cfsval = czero
  DO n = 0, ns
     theta = -pi + (REAL(n) + 0.5d0)*h
     IF ( MOD(theta,pi) == 0.0) theta = eps0  ! shift away from singularity, since there is a sin(c2*theta) in the denominator
     z = shift + REAL(ns)/t*(c1*theta/TAN(c2*theta) - c3 + c4*ic*theta)
     dz = REAL(ns)/t*( -c1*c2*theta/(SIN(c2*theta)*SIN(c2*theta)) + c1/TAN(c2*theta) + c4*ic)
     tval = cmult(sfunc(z,k,beta),dz)  ! function evaluation in Laplace space
     tval = cmult(EXP(z*t),tval)
     cfsval = cfsval + tval
  END DO

  fsval = REAL( h/(2.d0*pi)*cmult(-ic,cfsval) )
!  fsval = h/(2.d0*pi)*cmult(-ic,cfsval)

END FUNCTION talbot



COMPLEX(kind=8) FUNCTION sfunc(s,k,beta) RESULT(fval)
!
! defines the laplace transformed function
!
  IMPLICIT NONE
  REAL(kind=8), PARAMETER :: pi = 4.d0*ATAN(1.d0)
  REAL(kind=8), INTENT(in) :: k,beta
  REAL(kind=8) :: kr
  COMPLEX(kind=8), INTENT(IN) :: s

  COMPLEX(kind=8), EXTERNAL :: cdiv,cpow

  kr = k/SQRT(2.d0*pi)
!  kr = k

  ! f(t) = gaussian(x,t)
  fval = cdiv( cpow(s,beta-1.0), cpow(s,beta) + ABS(kr)**2.d0 )

END FUNCTION sfunc



COMPLEX(kind=8) FUNCTION cmult(a,b) RESULT(c)
!
! multiplication of two complex numbers a*b
!
  IMPLICIT NONE
  COMPLEX(kind=8), INTENT(IN) :: a,b

  c = CMPLX( REAL(a)*REAL(b) - AIMAG(a)*AIMAG(b), REAL(a)*AIMAG(b) + AIMAG(a)*REAL(b) )

END FUNCTION cmult



COMPLEX(kind=8) FUNCTION cdiv(a,b) RESULT(c)
!
!  division of two complex numbers a/b
!
  IMPLICIT NONE
  REAL(kind=8) :: a2, b2
  COMPLEX(kind=8), INTENT(IN) :: a,b

  COMPLEX(kind=8), EXTERNAL :: cmult

  a2 = ABS(a)*ABS(a)
  b2 = ABS(b)*ABS(b)

  c = cmult(a,CONJG(b))/b2

END FUNCTION cdiv


COMPLEX(kind=8) FUNCTION cpow(a,b) RESULT(c)
!
!  exponentiate a complex number a to power b; a^b
!  -- uses Euler notation z = r*exp(ic*theta)
!
  IMPLICIT NONE

  REAL(kind=8), INTENT(IN) :: b
  REAL(kind=8) :: r,theta
  COMPLEX(kind=8), INTENT(IN) :: a
  
  IF (b .EQ. 0.0) THEN
     c = 1.d0  ! catch the trivial case
  ELSE
     r = ABS(a)
     theta = ATAN2(AIMAG(a),REAL(a))

     c = r**b*CMPLX( (COS(b*theta)), SIN(b*theta) )

  END IF

!  ELSE
!     WRITE(*,*) "Incorrect power in CPOW. b = ", b
!     STOP
!  END IF

END FUNCTION cpow
