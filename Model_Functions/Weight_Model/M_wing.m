function M_wing = Mass_Wing(c, b);
    S = c * b; % wing area in m^2
    M_wing = -0.0802 + 2.2854 * S; % wing weight, where S is the wing area in m^2 
end