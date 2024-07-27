function C_D = C_D_calculation(c_l, c, b)
    C_D_min = 0.0397; 
    AR = b/c;
    e = 0.8;
    C_D_i = (c_l.^2) / (pi * AR * e); % wing lift-curve slope
    C_D = C_D_min + C_D_i; % total aircraft lift coefficient
end