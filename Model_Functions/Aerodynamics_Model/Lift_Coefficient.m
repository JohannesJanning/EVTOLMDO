function C_L = C_L_calculation(alpha_deg,c, b)
    alpha_rad = deg2rad(alpha_deg);
    AR = b/c;
    e = 0.8;
    a_airfoil = 2 * pi;
    a_wing = a_airfoil / (1 + (a_airfoil / (pi * AR * e))); % wing lift-curve slope
    cl_0 = 0.236; % NACA 2412 profile
    C_L = a_wing * alpha_rad + cl_0; % total aircraft lift coefficient
end