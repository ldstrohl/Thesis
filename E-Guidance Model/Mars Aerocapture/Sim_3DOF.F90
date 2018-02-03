	PROGRAM MAIN
    !	*******************************************************************
    !   Main program to simulate 3DOF aerobraking trajectories
    !   Closed-Loop Simulation
    !
    !   For different vehicle model, initial conditions, and targeted apoapsis 
    !   the user simply changes the data at the top part of this code, 
    !   and replaces the aerodynamic model CLCD(mach, alpha)and
    !   trim_alpha(mach) used in this code with the new models.
    !  
    !   This version uses the Orion MPCV model
    !   10/4/2013
    !*********************************************************************
        USE PC_GUID_PARAM
        USE INITIAL_GUID_PARAM
        USE GLOBAL_CONSTANTS
        USE MarsGRAM2010Interface
        IMPLICIT NONE                    

        REAL(8) :: S_ref, M,  PHIBKRATEMAXDEG, PHIBKACCELMAXDEG
        REAL(8) :: Amax, Qmax, qbar_max, Damax, DQmax
        REAL(8) :: inclination, Omega
        REAL(8) :: altF, VF, r_at, r_pt, dt, t0, bank0
        REAL(8), DIMENSION(6) :: x0, xx, xnew
        INTEGER :: Mission_Case, iter, Nstate, STOP_END
        REAL(8) :: r,lon, lat, v, gamma, psi, lonsite, latsite, dt_guid
        REAL(8) :: e, ef,  bankNAV, range, Az, Psi0,range0, time
        REAL(8) :: alt, vel, Vs, rho, rho_a, CL, CD, L, D
        REAL(8) :: VV,VJ, h, mach, sigma_CMD, delta_t, dummy
        REAL(8) :: t1, t2,rf, alt_F, Lat_geocentric, radius
        REAL(8) :: lat_geod, lat_geoc, height, r_initial
        REAL(8) :: Qdot,a_g, qbar, V315
        INTEGER :: Guid_init, i
        INTEGER ::  UP =0
        REAL(8) :: L_factor, D_factor
        REAL(8) :: rdot,r_f, r_exit, V_exit, gamma_exit, r_a, V_a
        REAL(8) :: r_p, V_at, V_pt
        INTEGER :: MODE
        INTEGER, SAVE :: flag
        REAL(8), SAVE :: sigma
        REAL(8), SAVE :: alpha
        
           !! bankCMD:    Commanded bank angle                        (deg)
    !! alpha:      angle of attack                             (deg)
    !! flag        an integer flag;   flag=-1: no trajectory found to meet the

		REAL(8) :: trim_alpha		
        REAL(8) :: temperature, pressure	
        REAL(8), DIMENSION(3) :: r_I, V_I, h_vec
        REAL(8) :: mag_h, incl, Vi
        REAL(8) :: M_ref, Guid_cycle, t_guid
        
        REAL(8) :: Alt_table(500), Rho_table(500)
        
        ! MC variables
        CHARACTER(20) :: str
        INTEGER :: numruns, counter, disp_flag, writer, write_count
        REAL(8) :: mass_mult, CL_mult, CD_mult, rdisp1, rdisp2, rdisp3, vdisp1, vdisp2, vdisp3
        REAL(8) :: dlat, dlon, dv, dgamma, dpsi, r_disp, V_disp, mass_disp, CL_disp, CD_disp, month_rand, fpa_disp
        REAL(8), dimension(6,1) :: DEOstate
        REAL(8), dimension(5,1) :: DEIstate
        REAL(8), dimension(5,6) :: stmEI
        
        ! Mars GRAM variables
        type(gram_state) :: state
        type(gram_settings) :: settings
        INTEGER :: error
        
        EXTERNAL EOM_3DOF	
        
        OPEN(UNIT = 16, FILE = "../Data/Mars_density.DAT", STATUS = 'UNKNOWN')
        OPEN(UNIT = 17, FILE = "TRAJ1.DAT", STATUS = 'UNKNOWN')
        OPEN(UNIT = 18, FILE = "GRAM_density.DAT", STATUS = 'UNKNOWN')
        OPEN(UNIT = 26, FILE = "Final_states.DAT", STATUS = 'UNKNOWN')
        
        DO i=1,186
            READ(16,100)    Alt_table(i), Rho_table(i)
        END DO
100     FORMAT(2E18.7)    
        CLose (16)
    
        dt        = 0.2                 ! simulation step size, sec  
        Guid_cycle= 1.0                 ! Guidance cycle time
        t0        = 0.0                 ! initial time, sec
        t_guid    = Guid_cycle
      
    
        PHIBKRATEMAXDEG   = 10.D0       ! max bank rate, deg/sec
        PHIBKACCELMAXDEG  = 5.0D0       ! max bank acceleration, deg/sec^2
    
        !   Set the  initial conditions
        
        Mission_Case  = 2
        Numruns = 3
        Counter = 0
        write_count = 0
        writer = 1
        disp_flag = 1                   ! Flag to disperse ICs, mass, CL, CD, and density
        if (disp_flag == 0) then
            Numruns = 1
        end if 
        !   If running Monte Carlo, set dispersion percentages. The loaded data files are
        !   normal with 0 mean and 3-sigma values of +-1, so the data files will be 
        !   multiplied by these values to achieve the desired dispersions. I.E. for a 
        !   normal mass dispersion with 3-sigma +- 10 percent, set mass_disp = 0.10.
        r_disp = 0.0        ! +- 0 m
        V_disp = 0.0        ! +- 0 m/s
        mass_disp = 0.0     ! +- 0-percent dispersion
        CL_disp = 0.10      ! +- 10-percent dispersion
        CD_disp = 0.10      ! +- 10-percent dispersion
        
        !   Open the necessary dispersion seeds
        IF (disp_flag == 1) THEN
            OPEN(UNIT = 19, FILE = "r_seeds.dat", STATUS = 'UNKNOWN')          ! 3,000x3 normal position dispersion seeds, 0 mean, 3-sigma = +- 1
            OPEN(UNIT = 20, FILE = "V_seeds.dat", STATUS = 'UNKNOWN')          ! 3,000x3 normal velocity dispersion seeds, 0 mean, 3-sigma = +- 1
            OPEN(UNIT = 21, FILE = "mass_seeds.dat", STATUS = 'UNKNOWN')       ! 3,000x1 normal mass dispersion seeds, 0 mean, 3-sigma = +- 1
            OPEN(UNIT = 25, FILE = "CL_seeds.dat", STATUS = 'UNKNOWN')         ! 3,000x1 normal CL dispersion seeds, 0 mean, 3-sigma = +- 1
            OPEN(UNIT = 23, FILE = "CD_seeds.dat", STATUS = 'UNKNOWN')         ! 3,000x1 normal CD dispersion seeds, 0 mean, 3-sigma = +- 1
            OPEN(UNIT = 30, FILE = "fpa_seeds.dat", STATUS = 'UNKNOWN')        ! 3,000x1 normal fpa dispersion seeds, 0 mean, 3-sigma = +- 1
        END IF
        
        
        !   Set the GRAM settings
        settings%datadir = '../../Release1.0_Nov10/binFiles/'   ! directory for COSPAR data and topographic height data, must give a local path
        settings%lendatadir = LEN(settings%datadir)
        settings%gcmdir = '../../Release1.0_Nov10/binFiles/'    ! Directory for GCM binary data files, must give a local path
        settings%lengcmdir = LEN(settings%gcmdir)
        !settings%myear = 2003     !
        !settings%month = 4        !
        !settings%mday = 10           ! initial date and time, myear can be 2-digit between 1970-2069, otherwise must be 4-digit
        !settings%ihr = 0           !
        !settings%imin = 0         !
        !settings%sec = 0         !
        settings%deltatex = 0.0     ! default is 0.0, adjustment for exospheric temperature
        settings% rpscale = 1.0     ! default is 1.0, density scale
        settings%mapyear = 0        ! default is 1, 1 or 2 for TES mapping year 1,2 GCM input data or 0 for Mars-GRAM 2001 GCM input data sets
        settings%rwscale = 1.0      ! default is 1.0, random wind perturbation scale factor (>=0)
        settings%wlscale = 1.0      ! default is 1.0, scale factor for perturbation wavelengths (0.1-10) 
        settings%blwinfac = 1.0     ! default is 1.0, scale factor for boundary layer slope winds (0 = none)
        settings%wmscale = 1.0      ! default is 1.0, scale factor for mean winds
        settings%corlmin = 0.0      ! default is 0.0, Minimum relative step size for perturbations (0.0 - 1.0)
        settings%dusttau = 0.45      ! default is 0.3, Optical depth of background dust level (no time-developing dust storm, just uniformly mixed dust), 0.1 to 3.0, or use 0 for assumed seasonal variation of background dust  
        settings%dustmin = 0.3      ! default is 0.3, Minimum seasonal dust tau if input Dusttau=0 (>=0.1)
        settings%dustmax = 1.0      ! default is 1.0, Maximum seasonal dust tau if input Dusttau=0 (<=1.0)
        settings%dustnu = 0.003     ! default is 0.003, Parameter for vertical distribution of dust density (Haberle et al., J. Geophys. Res., 104, 8957, 1999)
        settings%dustdiam = 5.0     ! default is 5.0, Dust particle diameter (micrometers, assumed monodisperse)
        settings%dustdens = 3000.0  ! default is 3000.0, Dust particle density (kg/m**3)
        settings%als0 = 0.0         ! default is 0.0, starting Ls value (degrees) for dust storm (0 = none)
        settings%intens = 0.0       ! default is 0.0, dust storm intensity (0.0 - 3.0). Storm intensity (>0) is added to Dusttau. 
        settings%radmax = 0.0       ! default is 0.0, max. radius (km) of dust storm (0 or >10000 = global)
        settings%dustlat = 0.0      ! default is 0.0, Latitude (degrees) for center of dust storm
        settings%dustlon = 0.0      ! default is 0.0, Longitude (degrees) (West positive if LonEW = 0 ,or East positive if LonEW = 1) for center of dust storm
        settings%f107 = 68.0        ! default is 68.0, 10.7 cm solar flux (10**-22 W/cm**2 at 1 AU)
        settings%nr1 = 1001         ! default is 1001, starting random number (0 < NR1 < 30000)
        settings%zoffset = 3.25     ! default is 3.25, constant height offset (km) for MTGCM data or constant part of Ls-dependent (Bougher) height offset (0.0 means no constant offset). Positive offset increases density, negative offset decreases density.
        settings%ibougher = 2
        
        
    !   Loop containing Monte Carlo runs
    DO WHILE (Counter < Numruns)
        
        !   Re-initialize step sizes and initial time
        dt        = 0.2                 ! simulation step size, sec  
        Guid_cycle= 1.0                 ! Guidance cycle time
        t0        = 0.0                 ! initial time, sec
        t_guid    = Guid_cycle
        
        Counter = Counter + 1
        write_count = write_count + 1
        print*, 'run number: ', Counter
        
        IF (write_count > 200) THEN
            write_count = 0
            writer = writer+1
            write(str, '(i10)')writer
            CLOSE(17)
            OPEN(UNIT = 17, FILE = "TRAJ"//trim(adjustl(str))//".DAT", STATUS = 'UNKNOWN')
        END IF
        
        !   Generate dispersions for density, CL/CD, and mass, and month
        IF (disp_flag == 1) THEN
            settings%nr1 = settings%nr1 + 5 ! GRAM random number seed
            read(21,*) mass_mult        ! mass multiplier
            read(25,*) CL_mult          ! CL multiplier
            read(23,*) CD_mult          ! CD multiplier
            read(30,*) fpa_disp         ! fpa dispersion
            
            !CALL RANDOM_NUMBER(month_rand)
            !config%mn = CEILING(12*month_rand)
            !print*,'month:',config%mn
            
        ELSE IF (disp_flag == 0) THEN
            settings%rpscale = 0.0
            settings%rwscale = 0.0
        END IF
        
        IF (Mission_Case == 1) THEN     ! CobraMRV from Dan Matz (1/27/2016)
            
            settings%myear = 2003     !
            settings%month = 4        !
            settings%mday = 10           ! initial date and time, myear can be 2-digit between 1970-2069, otherwise must be 4-digit
            settings%ihr = 0           !
            settings%imin = 0         !
            settings%sec = 0         !
            
            S_ref   = 62.21             ! reference area of the vehicle (m^2)
            M_ref   = 58700.0           ! nominal mass of the vehicle (kg)
            M       = M_ref             ! dispersed mass      
            
      	    x0(1)	= 125.0             ! plantocentric altitude above reference ellipsoid (km)
            x0(2)   =-176.4016724       ! plantocentric longitude (deg)
            if (x0(2)<0) then
                x0(2) 	= (360.0+x0(2))
            end if
            x0(3)   =-22.4067783        ! plantodetic latitude (deg)
            x0(4)   = 6204.31           ! velocity (m/s), relative
            x0(5)   =-10.9923           ! flight path angle(relative),deg
            x0(6)   =-2.1753            ! heading angle (deg)
      
            bank0   = 0.0               ! initial bank angle (deg)
            
        ELSE IF (Mission_Case == 2) THEN  ! From Dam Matz (2/9/2016)
            settings%myear = 2017       !
            settings%month = 5          !
            settings%mday  = 5          ! initial date and time, myear can be 2-digit between 1970-2069, otherwise must be 4-digit
            settings%ihr   = 8          !
            settings%imin  = 46         !
            settings%sec   = 17.9       !
            
            S_ref   = 62.21             ! reference area of the vehicle (m^2)
            M_ref   = 60000.0           ! nominal mass of the vehicle (kg)
            M       = M_ref             ! dispersed mass      
            
      	    x0(1)	= 125.0             ! plantocentric altitude above reference ellipsoid (km)
            x0(2)   =-176.4016724       ! plantocentric longitude (deg)
            if (x0(2)<0) then
                x0(2) 	= (360.0+x0(2))
            end if
            x0(3)   =-22.4067783        ! plantodetic latitude (deg)
            x0(4)   = 6204.31           ! velocity (m/s), relative
            x0(5)   =-10.6925           ! flight path angle(relative),deg
            x0(6)   =-2.17053            ! heading angle (deg)
      
            bank0   = 0.0               ! initial bank angle (deg)           
        END IF
   
        !   Disperse mass
        IF (disp_flag == 1) THEN
            M = M*(1 + mass_mult*mass_disp)
        END IF

        !   Disperse FPA for Powell Comparison
        IF (disp_flag == 1) THEN
            x0(5) = x0(5) + 0.1*fpa_disp
        END IF
                
    ! The following ususally need not be changed    
    ! -------------------------------------------------------------------
        
        height   = x0(1)                ! initial geodetic altitude in km
        ! Get initial geocentric radius r in km and geocentric latitude in rad
        CALL radius_from_geodetic(x0(2)/r_to_d,x0(3)/r_to_d,height,r_initial,lat_geoc)
        
        r       = r_initial
        lon     = x0(2)
        lat     = lat_geoc*r_to_d
        v       = x0(4)
        gamma   = x0(5)
        psi     = x0(6)
        
    ! Normalization
        x0(1)   = 1000.0*r/R0
        if (x0(2)<0.0) then
            x0(2) 	= (360.0+x0(2))
        end if
        x0(2)   = x0(2)/r_to_d
        !x0(3)   = x0(3)/r_to_d
        
        x0(3)   = lat_geoc
        
        x0(4)   = x0(4)/Vscale
        x0(5)   = x0(5)/r_to_d
        x0(6)   = x0(6)/r_to_d
        
        iter    = 0
        t1      = 0
        dt_guid = dt
        dt      = dt/tscale
        Nstate  = 6
        STOP_END= 0
     
        xx      = x0
        sigma   = bank0/r_to_d
        bankNAV = sigma*r_to_d
        !-----------------------------------------------------------------
        
        ! Disperse entry interface conditions
            IF (disp_flag == 1) THEN
                stmEI = reshape((/ 0.000431877, 0.000474938, 0.000523995, 1.98722e-05, -0.000282583, &
                                 & 0.000845483, 0.000947229, 0.000741702, 2.84545e-05, -0.000554336, &
                                 & -9.49675e-05, -9.86429e-05, 8.68356e-05, 2.80315e-06, 6.17977e-05, &
                                 & -0.43473, -0.481676, -0.351638, -0.021562, 0.278215, &
                                 & 0.562317, 0.624282, 0.248497, 0.0163075, -0.364653, &
                                 & -0.830862, -0.932346, -0.62266, -0.0338761, 0.550808 /), (/5, 6/))
                read(19,*) rdisp1,rdisp2,rdisp3        ! initial position vector dispersion
                read(20,*) vdisp1,vdisp2,vdisp3        ! initial velocity vecotr dispersion
                DEOstate = reshape((/ rdisp1*r_disp, rdisp2*r_disp, rdisp3*r_disp, &
                    & vdisp1*v_disp, vdisp2*v_disp, vdisp3*v_disp /), (/6,1/))
                DEIstate = MATMUL(stmEI,DEOstate)
                dlat = DEIstate(1,1)/r_to_d
                dlon = DEIstate(2,1)/r_to_d
                dV = DEIstate(3,1)/Vscale
                dgamma = DEIstate(4,1)/r_to_d
                dpsi = DEIstate(5,1)/r_to_d
                
                xx(2) = xx(2) + dlon
                xx(3) = xx(3) + dlat
                xx(4) = xx(4) + dV
                xx(5) = xx(5) + dgamma
                xx(6) = xx(6) + dpsi
            END IF
        
            
        ! Disperse FPA for Powell Comparison
        !read(21,*) mass_mult  
        !xx(5) = xx(5)*(1 + 0.1*mass_mult)
        
        ! Start the closed-loop simulation
     
        time    = 0

        !=========== Closed-Loop Simulation Loop ================================    
        DO While (STOP_END == 0)             ! see the stopping criteria later
            iter    = iter+1
            t2      = t1+dt
            IF (iter == 1) THEN
                Guid_init = 0
            ELSE
                Guid_init = 1
            END IF  
            radius  = xx(1)*R0
            lon     = xx(2)*r_to_d
            lat     = xx(3)*r_to_d
            vel     = xx(4)*Vscale
            gamma   = xx(5)*r_to_d
            psi     = xx(6)*r_to_d
            bankNAV = sigma*r_to_d

            IF (iter == 1) THEN
                CALL Altitude_above_Ellipsoid(xx(2),xx(3),xx(1)*R0/1000, alt,lat_geod)
                CALL Vsound(alt,Vs)                     ! speed of sound, m/s
                mach        = vel/Vs                    ! Mach number
                alpha       = trim_alpha(mach)
                !CALL  CLCD_sim(mach,alpha, CL, CD) 
                
                CALL  CLCD(mach,alpha, CL, CD) 
                IF (disp_flag == 1) THEN
                    CL = CL*(1 + CL_mult*CL_disp)
                    CD = CD*(1 + CD_mult*CD_disp)
                END IF
                
                !   Setup Mars GRAM and get the first density output
                CALL gram_initialize(t1, alt, lat, lon, settings, state, error)
                rho = state%density
                
              	!CALL density_sim(alt, rho,rho_a)            ! air density (kg/m^3)
                !CALL density(alt,rho, rho_a)        ! air density (kg/m^3)
                !CALL density_Mars(Alt_table, Rho_table, alt, rho,rho_a)            ! air density (kg/m^3)
                
                L   = 0.5*rho*vel**2*S_ref*CL/M         ! Lift acceleration (m/s^2)
                D   = 0.5*rho*vel**2*S_ref*CD/M         ! Drag acceleration (m/s^2)
            ELSE
                L   = L*g0
                D   = D*g0
            END IF

           IF(t_guid >= Guid_cycle) THEN
            ! Call the entry guidance interface to get the guidance commands
                CALL Guidance_Interface(radius,lon,lat,vel, gamma,psi,bankNAV,      &
                    & M_ref,S_ref,Guid_cycle, L,D,time,Guid_init, sigma, alpha, flag)
                t_guid = 0.0
                sigma  = sigma/r_to_d
                alpha  = alpha/r_to_d
           
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
        !   Apply rate and acceleration limiter for bank angle     
                if (iter == 1) then   !initialization for the bank angle limter function
                    sigma_CMD   = bankNAV/r_to_d 
                    delta_t     = 0
                    CALL shortest_bank_limiter(sigma_CMD,PHIBKRATEMAXDEG,&
                     & PHIBKACCELMAXDEG,delta_t,dummy)
                end if
              
                sigma_CMD   = sigma 
                delta_t     = Guid_cycle
            
                CALL shortest_bank_limiter(sigma_CMD,PHIBKRATEMAXDEG,&
                     & PHIBKACCELMAXDEG,delta_t,sigma)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
           END IF   
           t_guid = t_guid+dt_guid
            
        ! Integrate the 3DOF equations of motion one time step ahead
     
            CALL RK4_3DOF(EOM_3DOF, t1, t2, xx, Nstate,M,S_ref,sigma,L/g0,D/g0,xnew)  
            xx 	= xnew
            t1  = t2
            time= t2*tscale

            IF(xx(5) > 0.0) THEN
                UP = 1
            END IF
            
            
            ! Check the conditions for stopping the simulation; the user
            !  can add or modify these stopping criteria
            CALL Altitude_above_Ellipsoid(xx(2),xx(3),xx(1)*R0/1000,alt,lat_geod);

            rdot    = xx(4)*sin(xx(5))
            r_f     = alt*1000/R0 + rdot*dt
            r_exit  = 125000.0/R0
            IF(r_f > r_exit  .AND. UP == 1) THEN
                IF(abs(rdot)>1.0D-8) THEN
                    dt = (r_exit-alt*1000/R0)/rdot
                    IF(abs(dt)>1./tscale) THEN
                        dt = 1.*dt/(abs(dt)*tscale)
                    END IF
                ELSE
                    STOP_END = 3
                END IF
            END IF
            IF(abs(alt*1000/R0-r_exit) <= 100./R0 .AND. UP == 1) THEN
                STOP_END = 1
            END IF
            IF(time > 1700) THEN
                STOP_END = 2
            END IF

            ! print*, 'h=', alt
            ! The following records the data for plotting and output
            ! The user can modify the following as desired
 
            VV          = xx(4)
            CALL Vsound(alt,Vs)                 ! speed of sound, m/s            
            VJ   		= Vscale*VV             ! m/s
            mach 		= VJ/Vs                 ! Mach number
            
            ! Update Mars GRAM and get new density
            CALL gram_update(t1, alt, xx(3)*r_to_d, xx(2)*r_to_d, settings, state, error)
            rho = state%density
 
            IF (UP == 0) THEN
             !   Write(18, 112) (xx(1)*R0 - R0)/1000, rho
                Write(18, 112) alt, DLOG(rho)
112             FORMAT(2F18.12)
            END IF     
            
            !CALL density_sim (alt,rho, rho_a)        ! air density (kg/m^3)
            !CALL density(alt,rho, rho_a)        ! air density (kg/m^3)
            !CALL density_Mars(Alt_table, Rho_table, alt, rho,rho_a)            ! air density (kg/m^3)

            !CALL CLCD_sim(mach,alpha,CL,CD)
            CALL CLCD(mach,alpha,CL,CD)
            IF (disp_flag == 1) THEN
                CL = CL*(1 + CL_mult*CL_disp)
                CD = CD*(1 + CD_mult*CD_disp)
            END IF
            
            ! Add uncertainty
            CL = 1.0*CL
            CD = 1.0*CD
            
            L   = 0.5*rho*VV**2*R0*S_ref*CL/M   ! Lift acceleration (g)
            D   = 0.5*rho*VV**2*R0*S_ref*CD/M   ! Drag acceleration (g)
   
            A_g = sqrt(D**2+L**2)                !load factor (g)
            A_g  = A_g*g0/9.798287622535D0          ! in Earth g
                        
            V315 = VJ**3.15
            Qdot = 5.21e-5*sqrt(rho/0.3048)*V315/11348.93 ! heat rate BTU/ft^2-sec
            qbar	= 0.5*rho*VJ*VJ/47.88           ! dynamic pressure, psf
            
            CALL RELATIVE_TO_INERTIAL(xx, r_I, V_I)  
            Vi = sqrt(DOT_PRODUCT(V_I,V_I))*Vscale
 
            Write(17, 111) TIME, ALT, xx(2)*r_to_d, xx(3)*r_to_d,VV*Vscale,xx(5)*r_to_d,&
            &sigma*r_to_d,a_g,alpha*r_to_d,mach,Qdot, qbar, Vi !, CL, CD, rho

 111    	    FORMAT(13F14.3)
            
            !print*, 'time =', TIME
        
        END DO   

        !end of integration for closed-loop simulation; 
        ! Retrieve pperation mode 
        MODE = GUID_MODE                ! 1 = Optimal agogee-targeting mode
                                        ! 2 = Optimal DeltaV-minimization mode
                                        
        r_pt = GUID_r_pt +R0
        r_at = GUID_r_at +R0

        CALL Orbit(xx, r_a,r_p, e, V_a, V_exit,gamma_exit)
        CALL RELATIVE_TO_INERTIAL(xx, r_I, V_I)
        CALL CROSS(r_I, V_I, h_vec)             ! compute angular momentum vector
        mag_h=sqrt(DOT_PRODUCT(h_vec, h_vec))
        incl = acos(h_vec(3)/mag_h)*r_to_d      ! final orbital inclination 
        V_at = sqrt(2.0*(-1.0/(r_at+r_pt)+1.0/r_at)*R0) ! targeted apoapsis velocity
        V_pt = sqrt(2.0*(-1.0/(r_at+r_pt)+1.0/r_pt)*R0) ! targeted apoapsis velocity
        print*, 'predicted apogee altitude (km)        =', (r_a-1)*R0/1000
        print*, 'predicted apogee velocity (km/s)      =', V_a*Vscale/1000
        print*, 'exit inertial velocity (km/s)         =', V_exit*Vscale/1000
        print*, 'exit inertial flight-path angle (deg) =', gamma_exit*r_to_d
        print*, 'final orbital inclination (deg)       =', incl
      IF(MODE == 1 .OR. MODE == 3) THEN
        print*, 'Delta V (m/s)                         =', Vscale*V_at-V_a*Vscale        
      ELSE ! MODE = 2 or 4
        print*, 'Delta V1 (m/s)                        =', PC_DV1*Vscale
        print*, 'Delta V2 (m/s)                        =', PC_DV2*Vscale
        print*, 'Total Delta V (m/s)                   =', Vscale*(PC_DV1+PC_DV2)       
      END IF 
       
      ! Write data for final apogee altitude, velocity, intertial velocity, inertial fpa,
      ! final orbital inclination, and Delta V
      IF (MODE == 1 .OR. MODE == 3) THEN    
        Write(26, 111) (r_a-1)*R0/1000, V_a*Vscale/1000, V_exit*Vscale/1000, gamma_exit*r_to_d,&
        & incl, abs(Vscale*V_at-V_a*Vscale), PC_DV1*Vscale, PC_DV2*Vscale, (r_p-1)*R0/1000

        
      ELSE IF (MODE == 2 .OR. MODE == 4) THEN
        Write(26, 111) (r_a-1)*R0/1000, V_a*Vscale/1000, V_exit*Vscale/1000, gamma_exit*r_to_d,&
        & incl, PC_DV1*Vscale, PC_DV2*Vscale, (r_p-1)*R0/1000

      END IF
            
      ! Close Mars GRAM data files so they can be reopened for the next MC run
      CLOSE(9)      ! molatoph.bin
      CLOSE(10)     ! COSPAR2.DAT
      CLOSE(11)     ! hgtoffst.dat
      CLOSE(12)     ! albedo1.bin
    END DO
    
    CLOSE(19)       ! r_seeds.dat
    CLOSE(20)       ! V_seeds.dat
    CLOSE(21)       ! mass_seeds.dat
    CLOSE(25)       ! CL_seeds.dat
    CLOSE(23)       ! CD_seeds.dat
    CLOSE(17)		! TRAJ.DAT
    CLOSE(18)       ! Density log
    CLOSE(26)       ! Final orbital data
    CLOSE(30)       ! fpa_seeds.dat
	END PROGRAM  MAIN
	
