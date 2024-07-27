function Power_Descent_BEM = Power_Descent_BEM_calculation(MTOM, rho, c, b, V_descent, ROD, CD, R_prop, num_props, eta_c)
    P_min = 5000; % Minimum power requirement for stabilization in W
    g = 9.81; % earth acceleration in m/s^2
    W = MTOM * g; % aircraft weight in N
    S = c * b; % wing area in m^2 
    V_descent_calc = sqrt(V_descent^2 + ROD^2); % total descent velocity in m/s
    D = 0.5 * rho * V_descent_calc^2 * S * CD; % drag in descent N
    T_required = D - (W * ROD / V_descent_calc); % thrust required in descent in N
    % Ensure T_required is non-negative
    if T_required < 0
     T_required = 0;
    end
    A_prop = pi * R_prop^2; % propeller disk area in m^2 
    T_prop = T_required / num_props; %thrust per propeller in N
    sigma = T_prop / A_prop; % disk loading in N/m^2 
    vi = sqrt(T_prop / (2 * rho * A_prop)); % induced velocity in m/s
    Power_Descent = (T_required * V_descent_calc) / eta_c + (T_prop * vi) / eta_c * num_props; % Power required for descent (P_descent) including induced power
    Power_Descent_BEM = max(Power_Descent, P_min); % total power in descent, considering total minimum power required
end