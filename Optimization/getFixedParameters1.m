function params = getFixedParameters()
    % Fixed parameters as derived and defined in the model documentation 
    addpath('/Users/johannes/Dropbox/Mein Mac (Johanness MacBook Pro (2))/Documents/MATLAB/Master_Project/MDO_MOO_Model_v100724/');

    % Standard Environmental Settings
    params.g = 9.81; % Acceleration due to gravity in m/s^2
   
    % Model Flight Height Settings
    params.h_tofldg = 0; % [m], take-off and landing height AMSL, (above mean sea level)
    params.h_hover = 30.5; % [m] Hover height AGL, (above ground level)
    params.h_cruise = 304.7; % [m], Cruise height AGL (above ground level)
    
    % Atmospheric Model 
    rho_msl = 1.225; % [kg/m^3], air density at main sea level (MSL)
    rho_tofldg = rho_msl*(1-22.558*10^(-6)*params.h_tofldg)^4.2559; % [kg/m^3], air density at take-off and landing area
    rho_cruise = rho_msl*(1-22.558*10^(-6)*(params.h_cruise + params.h_tofldg))^4.2559; % [kg/m^3], air density at cruise height
    rho_climb = rho_msl*(1-22.558*10^(-6)*(params.h_cruise/2 + params.h_tofldg))^4.2559; % [kg/m^3], assumed air density during climb 
    params.rho_hover = rho_tofldg; 
    params.rho = rho_cruise; 
    params.rho_climb = rho_climb; 
    
    % Aerodynamic Paramaters 
    params.alpha_deg_cruise = 3; % wing angle of attack in cruise [-]
    params.alpha_deg_climb = 8; % wing angle of attack in climb [-]
    params.alpha_deg_max = 15; % max wing angle of attack pre stall [-]

    % Design: Number of Props in eVTOL Configuration 
    params.num_props_cruise = 1; % number of props in cruise (horizontal flight) configuration -
    params.num_props_hover = 8; % number of prop in hover (vertical flight) configuration -
    params.d_tip = 0.4; % [m] % spacing between rotors and fuselage
    params.w_fuselage = 1.5; % [m] width of fuselage 

    % Mass assumptions 
    params.M_pax = 82; % [kg], Mass per passenger 
    params.M_lug = 16; % [kg], Mass per luggage 
    params.m_crew = 82.5 + 14; % [kg], EASA standard crew weight
    params.m_other = 668; % [kg], fuselage weight, and others (seats etc.). from tecnam P2006T 
    params.m_motor_initial = 500; % [kg], initial guess for motor mass
    params.m_bat_initial = 200; % [kg], initial guess for battery mass
    
    % Economic assumptions
    params.N_s = 4; % Number of seats  
    params.P_s_empty = 1436.5; % [€/kg of M_empty], weight specific aircraft acquisition price  
    params.N_wd = 260; % [days],  Number of working days per year 
    params.N_df = 6; % [-], Number of daily flights 
    params.P_e = 0.09; % [€/kWh], energy price 
    params.P_bat_s = 115; % [€/kWh], energy capacity specific battery replacement price 
    params.t_res = 60 * 30; % [seconds], 30min reserve for VFR acc. FAA 
    params.t_hover = 60; % [seconds], hover time in seconds

    % Efficiencies 
    params.eta_e = 0.9; % electric efficiency 
    params.eta_p = 0.85; % horizontal propulsive efficiency 
    params.eta_hp = 0.7; % hover propulsion efficiency
    params.eta_h = params.eta_hp * params.eta_e; % Hover efficiency
    params.eta_c = params.eta_p * params.eta_e; % horizontal flight efficiency 

    % Transportation Mode Comparison Figure of Merit Weighting 
    params.time_weight = 0.25;
    params.co2_weight = 0.25;
    params.energy_weight = 0.25;
    params.costs_weight = 0.25;


end
