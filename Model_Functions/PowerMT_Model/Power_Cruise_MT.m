function Power_Cruise_BEM = Power_Cruise_Calculation(MTOM, rho, D, V_cruise, R_prop, num_props, eta_c, alpha_deg_cruise)
    g = 9.81; % earth acceleration in m/s^2
    alpha_rad_cruise = deg2rad(alpha_deg_cruise);
    T_required = D + 0* sin(alpha_rad_cruise) * MTOM * g; % thrust required in cruise in N
    A_prop = pi * R_prop^2; % propeller disk area in m^2
    T_prop = T_required / num_props; % thrust per propeller in N
    sigma = T_prop / A_prop; % Disk Loading in N/m^2
    vi = sqrt(T_prop / (2 * rho * A_prop)); % induced velocity in m/s
    Power_Cruise_BEM = (T_required * V_cruise) / eta_c + (T_prop * vi) / eta_c * num_props; % power in cruise condition in W 
end