    !**************************************************************************
	!> @ingroup		Guidance_Lib
	!! @brief		Calculate the commanded body axes; All vectores are in the same coordinate frame
	!!
	!! @param[in]   R       Radius vector       
	!! @param[in]   V       Velocity vector  
    !! @param[in]   PV      A vector pointing in the direction of the thrust
    !! @param[out]  ONE_B   Unit vector for body-longitudinal axis   
    !! @param[out]  ONE_N   Unit vector for body-normal axis       
    !! @param[out]  ONE_Y   Unit vector for body-Y axis 
    !! @param[out]  ALPHA   Angle of attack (rad) 
 	!!
	!! @date		June 28, 2016
	!<
	!**************************************************************************
	SUBROUTINE VACUUM_COMMAND(R, V, PV, ALPHA, ONE_B, ONE_N, ONE_Y)

		IMPLICIT NONE
		REAL(8), DIMENSION(3), INTENT(IN)  :: R, V, PV
		REAL(8), DIMENSION(3), INTENT(OUT) :: ONE_B, ONE_N, ONE_Y
		REAL(8), INTENT(OUT) :: ALPHA
        REAL(8) :: P_norm
        
        P_norm=DSQRT(DOT_PRODUCT(PV, PV))
        IF (P_norm < 5.0D-16) THEN
            P_norm = 5.0D-16
        END IF
        
        ONE_B = PV/P_norm
		CALL CROSS(ONE_B,R,ONE_Y)           ! GET Body Y-axis

		P_norm=DSQRT(DOT_PRODUCT(ONE_Y, ONE_Y))
        IF (P_norm < 5.0D-16) THEN
            P_norm = 5.0D-16
        END IF
        ONE_Y = ONE_Y/P_norm
        		
		CALL CROSS(ONE_Y,ONE_B,ONE_N)       ! Body-Normal axis
		
        ALPHA = DATAN2(-DOT_PRODUCT(V, ONE_N), DOT_PRODUCT(V, ONE_B))
        
		RETURN
    END SUBROUTINE VACUUM_COMMAND 
    
    
    !**************************************************************************	  
	!> @ingroup		Guidance_Lib
	!! @brief		COMPUTE YAW, PITCH and ROLL EULER ANGLE COMMANDS with respect to NED FRAME
    !!              in 3-2-1 rotation sequence
    !!
	!! @param[in]       ONEB			unit vector of body x-axis in guidance frame
	!! @param[in]       ONEY			unit vector of body y-axis in guidance frame
	!! @param[in]       ONEN			unit vector of body normal axis (the negative of z-axis)
	!! @param[out]      YAW             YAW angle of body frame with respect to NED frame (rad)			
	!! @param[out]      PITCH			PITCH angle of body frame with respect to NED frame (rad)	
	!! @param[out]      ROLL			ROLL angle of body frame with respect to NED frame (rad)	
	!! @param[out]      PITCH_UP		Flag to indicate if PITCH is near +/- 90 deg (Yes 1; No 0)	
	!!
	!! @date	    June 10, 2016
	!! @see			GLOBAL_CONSTANT
	!<
	!**************************************************************************	    
	SUBROUTINE NED_EULER_GUID(ONEB, ONEN, ONEY,YAW,PITCH, ROLL,PITCH_UP)
		USE GLOBAL_CONSTANT
		USE GUIDANCE_VARIABLE
		
		IMPLICIT NONE
		REAL(8), DIMENSION(3), INTENT(IN) :: ONEB, ONEN, ONEY
		REAL(8), INTENT(OUT) :: YAW, PITCH, ROLL
		INTEGER, INTENT(OUT) :: PITCH_UP
		REAL(8), DIMENSION(3) :: ONE_B, ONE_N, ONE_Y, ONE_Z
		REAL(8), DIMENSION(3,3) :: T_GL
        

	    PITCH_UP = 0
	       
        T_GL    = G_T_NG                        ! NED-to-Guidance coordinate transformation
        ONE_B= MATMUL(TRANSPOSE(T_GL), ONEB)    ! Convert body axes into NED frame
    	ONE_Y= MATMUL(TRANSPOSE(T_GL), ONEY)
    	ONE_N= MATMUL(TRANSPOSE(T_GL), ONEN)
    	
		ONE_Z   =-ONE_N
		PITCH   =-DASIN(ONE_B(3))
		! Note that YAW and Roll are not well defined if Pitch = +/- 90 deg
		   
        YAW     = DATAN2(ONE_B(2),ONE_B(1))     ! Yaw angle from inputs
        ROLL    = DATAN2(ONE_Y(3),ONE_Z(3))     ! Roll angle from inputs
        
	
        RETURN
	END SUBROUTINE NED_EULER_GUID
	
    
  