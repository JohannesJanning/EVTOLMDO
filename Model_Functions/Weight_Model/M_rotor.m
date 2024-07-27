function Mass_rotor = Mass_rotor_calculation(NP_hover, NP_cruise, R_hover, R_cruise);
    % NP = number of props , R = radius of the rotors
    Mass_rotor = NP_hover * (0.7484 * R_hover^2 - 0.0403 * R_hover) + NP_cruise * (0.7484 * R_cruise^2 - 0.0403 * R_cruise); % rotor weight, where R is the rotor radius in m 
end