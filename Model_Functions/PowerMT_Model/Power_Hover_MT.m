function Power_Hover_BEM = Power_Hover_BEM_Calculation(MTOM, rho, R_prop, num_rotors, eta_h)
    g = 9.81; % earth acceleration in m/s^2
    T_total = MTOM * g; % total thrust requiref for hover in N
    T_rotor = T_total / num_rotors; % thrust per rotor in N
    A_rotor = pi * R_prop^2; % single rotor disk area in m^2 
    A_total = num_rotors * A_rotor; % total rotor disk area in m^2 
    sigma = T_total / A_total; % total disk loading in N/m^2 
    vi = sqrt(T_rotor / (2 * rho * A_rotor)); % induced velocity in m/s
    P_i_BEM = (T_total * vi) / eta_h; % induced power in W
    Power_Hover_BEM = (MTOM * g / eta_h) * sqrt(sigma / (2 * rho)); % power equation in W 
end