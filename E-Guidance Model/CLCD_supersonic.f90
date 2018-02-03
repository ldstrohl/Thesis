	SUBROUTINE CLCD(Mach, alpha, CL, CD)
	!*************************************************************************
	!! Calculate CL and CD of a mid L/D vehicle; 
    !! Curve fitting with Chebyshev Rational Order 3/4
	!! Inputs:
    !!	    Mach  - Mach number (0.2<=Mach<=4.0)
    !!      alpha - angle of attack in rad (in the range of 0-90 deg )
	!! Outputs:
    !!      CL    - Aerodynamic lift coefficient
	!!      CD    - Aerodynamic drag coefficient
	!!
	!************************************************************************** 
    IMPLICIT NONE

    REAL(8), INTENT(IN) :: Mach
    REAL(8), INTENT(IN) :: alpha
    REAL(8), INTENT(OUT) :: CL
    REAL(8), INTENT(OUT) :: CD

    REAL(8) :: x,y
    REAL*8 a(14+1), c(14+1)
    ! Coefficients for CL
    DATA a(1)/-0.2424228870950600D0/
    DATA a(2)/0.05921206835755609D0/
    DATA a(3)/-0.3380659491734871D0/
    DATA a(4)/-0.07281123439928469D0/
    DATA a(5)/-0.05001198031877806D0/
    DATA a(6)/0.06248622392818820D0/
    DATA a(7)/0.3516747384983234D0/
    DATA a(8)/-0.02679704023052718D0/
    DATA a(9)/-0.5999160805871470D0/
    DATA a(10)/-0.09724607581034851D0/
    DATA a(11)/-0.05261262605549739D0/
    DATA a(12)/0.03250128613098667D0/
    DATA a(13)/0.04122927675359531D0/
    DATA a(14)/0.03160999711896428D0/
    DATA a(15)/0.005282611155395087D0/
    ! Coefficients for CD
    DATA c(1)/-0.4414108246546230D0/
    DATA c(2)/-0.1171017381396072D0/
    DATA c(3)/-0.8076156872155501D0/
    DATA c(4)/0.08214522584113917D0/
    DATA c(5)/1.163405207629440D0/
    DATA c(6)/0.02796962306484996D0/
    DATA c(7)/0.4756624586285432D0/
    DATA c(8)/-0.1565306453086145D0/
    DATA c(9)/-0.5009766764136040D0/
    DATA c(10)/-0.05157386885167327D0/
    DATA c(11)/-0.1078811591642589D0/
    DATA c(12)/0.08834231956268198D0/
    DATA c(13)/0.1247092097041742D0/
    DATA c(14)/0.06858882310706151D0/
    DATA c(15)/0.02671615282031569D0/
    
    x   = Mach
    y   = alpha
    
    IF(x <0.2D0)    x = 0.2D0
    IF (x > 4.0D0)  x = 4.0D0
    IF(y<0.0)       y =0.0
    IF(y>1.5708D0)  y=1.5708D0
    CALL EVALCRATL(14,0,0,x,y,a,CL)
    CALL EVALCRATL(14,0,0,x,y,c,CD)
    
    RETURN
  END
    
    
!----------------------------------------------------------
SUBROUTINE EVALCRATL(order, logx, logy, x, y, c, z)
!----------------------------------------------------------
INTEGER order,logx,logy
REAL*8 x,y,c(*),z
INTEGER tcnt,j,m
REAL*8 tx(12),ty(12),num,den
IF(logx.NE.1) THEN
  x=(x-(2.100000000000000D0))/(1.900000000000000D0)
ELSE
  x=(DLOG(x)-(-0.1115717756571049D0))/(1.497866136776995D0)
END IF
IF(logy.NE.1) THEN
  y=(y-(50.00000000000000D0))/(50.00000000000000D0)
ELSE
  y=(DLOG(y)-(0.000000000000000D0))/(0.000000000000000D0)
END IF
SELECT CASE (order)
  CASE (6)
    tcnt=3
  CASE (10)
    tcnt=4
  CASE (14)
    tcnt=5
  CASE (18)
    tcnt=6
  CASE (22)
    tcnt=7
  CASE (26)
    tcnt=8
  CASE (30)
    tcnt=9
  CASE (34)
    tcnt=10
  CASE (38)
    tcnt=11
  CASE (42)
    tcnt=12
  CASE DEFAULT
    z=0.0
    RETURN
END SELECT
IF(tcnt.GT.7) THEN
  IF(x.LT.-1.D0) x=-1.D0
  IF(x.GT. 1.D0) x= 1.D0
  IF(y.LT.-1.D0) y=-1.D0
  IF(y.GT. 1.D0) y= 1.D0
END IF
tx(1)=1.D0
ty(1)=1.D0
tx(2)=x
ty(2)=y
DO 10 j=3,tcnt
  tx(j)=2*x*tx(j-1)-tx(j-2)
  ty(j)=2*y*ty(j-1)-ty(j-2)
10 CONTINUE
m=2
num=c(1)
den=1.0+c(2)*tx(m)+c(3)*ty(m)
DO 20 j=4,order,4
  num=num+c(j)*tx(m)
  num=num+c(j+1)*ty(m)
  m=m+1
  den=den+c(j+2)*tx(m)
  den=den+c(j+3)*ty(m)
20 CONTINUE
IF(den.EQ.0.0) THEN
  z=0.0
ELSE
  z= (num/den)*(1.858055000000000D0)+(0.9183550000000000D0)
END IF
RETURN
END

