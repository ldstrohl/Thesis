   SUBROUTINE density(ALT, Rho, Rho_a)
	!*************************************************************************
	!! Calculate Mars atmosphere density and density gradient
	!!  (Based on Mars GRAM2010 data)				
	!!  
	!! RHO=DENSITY IN KG/M^3    
	!! RHO_A= d(RHO)/dalt
    !!  
	!! Input:
    !!	    ALT				geodetic altitude (KM)
	!! Outputs:
    !!      Rho			    atmosphere density (KG/M^3)
	!!      Rho_a           dRHO/dALT(ALT is km)
	!!
    !! June 17, 2015
	!**************************************************************************
		IMPLICIT NONE

		REAL(8), INTENT(IN) :: ALT
		REAL(8), INTENT(OUT) :: Rho
		REAL(8), INTENT(OUT)::  Rho_a
		REAL(8) :: a,b,c,d,e,f,g,h,i,x,y,num,den, rho2
!-----------------------------------------------------------
!   X= alt km)
!   y= log(rho) (rho in km/m^3)
!   y=(a+cx+ex^2+gx^3+ix^4)/(1+bx+dx^2+fx^3+hx^4) [NL]
         a= -4.422345201370313
         b= -0.04399163896948651
         c= 0.1130629761317862
         d= 0.0006555968430113561
         e= 0.0003003890333711271
         f= -3.709114148189494E-06
         g= -2.148583579589203E-05
         h= 8.137552734673534E-09
         i= 7.041319756467658E-08
!-----------------------------------------------------------
        x = ALT
        IF(ALT>125.0D0) THEN 
            x = 125.0D0
        ELSE IF(ALT<0.2D0) THEN
            x = 0.2D0
        END IF
        num = -4.422345201370313D0+x*(0.1130629761317862D0+&
            &x*(0.0003003890333711271D0+x*(-2.148583579589203D-05+&
            &x*(7.041319756467658D-08))))
        den = 1.0+x*(-0.04399163896948651D0+x*(0.0006555968430113561D0+&
            &x*(-3.709114148189494D-06+x*(8.137552734673534D-09))))
        
        y   = num/den
        Rho = DEXP(y)
        Rho_a = Rho*(den*(c+x*(2.0*e+x*(3.0*g+4.0*i*x))) &
            & - num*(b+x*(2.*d+x*(3.0*f+4.0*h*x))))/den**2
        
  RETURN
END
