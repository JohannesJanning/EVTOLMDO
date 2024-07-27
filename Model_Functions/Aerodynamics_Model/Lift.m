function Lift = Lift_calculation(rho, v, cl, c, b);
    S = c * b; % wing area in m^2 
    Lift = 0.5 * rho * v^2 * cl * S;
end