function Vs = Vsound(ALT)
% 	!**************************************************************************
% 	!! @brief		Calculate the Speed of Sound at Mars
% 	!!
% 	!! Calculate the speed of sound
% 	!!
% 	!! @param[in]	ALT				planetodetic altitude (KM)
% 	!! @param[in]	Vs				speed of sound at Alt (m/s)
% 	!!
% 	!<
% 	!!**************************************************************************
% 	    USE GLOBAL_CONSTANTS
%
% 		IMPLICIT NONE
% 		REAL(8), INTENT(IN) :: ALT
% 		REAL(8), INTENT(OUT) :: Vs
% 		REAL(8) :: HIGHT,A,B,C,D,E,F,G,C_H,C_I,C_J,C_K,X,DUM1,DUM2
%         REAL(8) :: KT2,KT3,GAM_MARS,GAS_CONST,T, ALT_M, H,I,ALT_MAX
%         !INTEGER :: Planet    !(1=Earth, 2 = Mars)


ALT_M = ALT * 1000.0; % ! ALTITUDE IN METERS

%             !COEFFICIENTS FOR CURVE FIT FOR 8TH ORDER POLY
A = 3.90812820e+02;
B = -5.07169986e+03;
C = 6.00818111e+04;
D = -3.71909017e+05;
E = 1.27052434e+06;
F = -2.47999976e+06;
G = 2.74880412e+06;
H = -1.60759656e+06;
I = 3.84982694e+05;

if (ALT_M > 130000.0)
    ALT_M = 130000.0;
elseif (ALT_M < 10000.0)
    ALT_M = 10000.0;
end

ALT_MAX = 130081.410992;
%             ! SCALING AND CENTERING FOR POLYNOMIAL FIT
X = ALT_M/ALT_MAX;
Vs = A + X*(B + X*(C + X*(D + X*(E + X*(F + X*(G + X*(H + I*X))))))) ;


end


