function params = getFixedParameters()
    % Fixed parameters as derived and defined in the model documentation 
    addpath('/Users/johannes/Dropbox/Mein Mac (Johanness MacBook Pro (2))/Documents/MATLAB/Master_Project/MDO_MOO_Model_v100724/');

    % Standard Environmental Settings
    params.g = 9.81; % Acceleration due to gravity in m/s^2
   
    % Model Flight Height Settings
    params.h_tofldg = 0; % [m], take-off and landing height AMSL, (above mean sea level)
    params.h_hover = 15.24; % [m] 30.5m = 100ft, 15.24m = 50ft, Hover height AGL, (above ground level)
    params.h_cruise = 1219.2; % [m], 304.7m = 1000ft, Cruise height AGL (above ground level), 1219.2m = 4000ft
    
    % Atmospheric Model 
    rho_msl = 1.225; % [kg/m^3], air density at main sea level (MSL)
    rho_tofldg = rho_msl*(1-22.558*10^(-6)*params.h_tofldg)^4.2559; % [kg/m^3], air density at take-off and landing area
    rho_cruise = rho_msl*(1-22.558*10^(-6)*(params.h_cruise + params.h_tofldg))^4.2559; % [kg/m^3], air density at cruise height
    rho_climb = rho_msl*(1-22.558*10^(-6)*(params.h_cruise/2 + params.h_tofldg))^4.2559; % [kg/m^3], assumed air density during climb 
    params.rho_hover = rho_tofldg; 
    params.rho = rho_cruise; 
    params.rho_climb = rho_climb; 
    params.rho_msl = rho_msl;
    
    % Aerodynamic Paramaters 
    params.alpha_deg_cruise = 3; % wing angle of attack in cruise [-]
    params.alpha_deg_climb = 8; % wing angle of attack in climb [-]
    params.alpha_deg_max = 15; % max wing angle of attack pre stall [-]

    % Design: Number of Props in eVTOL Configuration 
    params.num_props_cruise = 1; % number of props in cruise (horizontal flight) configuration -
    params.num_props_hover = 8; % number of prop in hover (vertical flight) configuration -
    params.d_tip = 0.025; % [m] % spacing between rotors and fuselage
    params.w_fuselage = 1.5; % [m] width of fuselage 
    params.l_fus_m = 6; % assuming 9m for fuselage length (= aircraft length) [see Tecnam P2006T]
    params.r_fus_m = 0.75; % assuming 1.5m for cabin width [see Tecnam P2006T]

    % Mass assumptions 
    params.M_pax = 82; % [kg], Mass per passenger 
    params.M_lug = 16; % [kg], Mass per luggage 
    params.m_crew = 82.5 + 14; % [kg], EASA standard crew weight
    params.m_other = 668; % [kg], fuselage weight, and others (seats etc.). from tecnam P2006T 
    params.m_motor_initial = 500; % [kg], initial guess for motor mass
    params.m_bat_initial = 200; % [kg], initial guess for battery mass
    
    % Battery assumptions
    params.rho_bat = 400; % [Wh/kg] fixed assumption on battery density (400Wh/kg acc. Uber Eleveate requirement)
    params.C_charge =1; % [1/h] C-rate (discharge rate of battery)


    % Economic assumptions
    params.N_s = 4; % Number of seats  
    params.P_s_empty = 1436.5; % [€/kg of M_empty], weight specific aircraft acquisition price  
    params.N_wd = 260; % [days],  Number of working days per year, based on uber
    params.N_df = 6; % [-], Number of daily flights 
    params.P_e = 0.096668; % [€/kWh], energy price 
    params.P_bat_s = 115; % [€/kWh], energy capacity specific battery replacement price 
    params.t_res = 60 * 30; % [seconds], 30min reserve for VFR acc. FAA 
    params.t_hover = 60; % [seconds], hover time in seconds
    params.T_D = 8*60*60; % assuming a fixed operation window of 8hours daily, based on uber 
    params.S_P = 45400; % [€] annual salary for one pilot, based on Uber 
    params.N_AC = 1; % [-] number of aircraft controlled by one pilot (>1, bunker based operation)
    params.U_pilot = 2000*60*60; %[sec] annual pilot utilization in seconds (2000hrs)
    params.GWP_battery = 124.5, % [kg CO2e/kWh] LIB NMC battery life cycle global warming impact
    params.GWP_energy = 0.37896; %[kg CO2/kWh] electrcity generation GWP for battery charging
    params.LF = 0.67; % load factor, based on uber 
    params.fare_km = 1.98; % euro per km per passenger, based on london taxi day fare



    % Efficiencies 
    params.eta_e = 0.9; % electric efficiency 
    params.eta_p = 0.85; % horizontal propulsive efficiency 
    params.eta_hp = 0.7; % hover propulsion efficiency
    params.eta_h = params.eta_hp * params.eta_e; % Hover efficiency
    params.eta_c = params.eta_p * params.eta_e; % horizontal flight efficiency 

    % Transportation Mode Comparison Figure of Merit Weighting 
    params.time_weight = 0.333;
    params.co2_weight = 0.333;
    params.energy_weight = 0;
    params.costs_weight = 0.333;


end
