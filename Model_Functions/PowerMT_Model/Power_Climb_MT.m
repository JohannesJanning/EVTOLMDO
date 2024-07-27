function Power_Climb_BEM = Power_Climb_Calculation(MTOM, rho, c, b, V_climb, ROC, C_D, R_prop, num_props, eta_c, alpha_deg_climb)
    g = 9.81; % gravitational acceleration m/s^2
    S = c * b; % wing area m^2 
    W = MTOM * g; % aircraft weight in N
    V_climb_calc = sqrt(V_climb^2 + ROC^2); % total climb velocity m/s
    alpha_rad_climb = deg2rad(alpha_deg_climb);
    D = rho / 2 * V_climb_calc^2 * S * C_D; 
    T_required = D + W * sin(alpha_rad_climb); % thrust required to overcome Drag and climb N
    A_prop = pi * R_prop^2; % propeller disk area in m^2 
    T_prop = T_required / num_props; % thrust required per propeller in N
    sigma = T_prop / A_prop; % Disk Loading in N/m^2
    vi = sqrt(T_prop / (2 * rho * A_prop)); % induced velocity at rotor m/s
    Power_Climb_BEM = (T_required * V_climb_calc) / eta_c + (T_prop * vi) / eta_c * num_props; 
end