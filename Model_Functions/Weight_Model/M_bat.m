% Battery Mass in kg 

function M = M_bat(rho_bat, E_trip, E_res)
    % Calculate the total battery mass
    M = (E_trip + E_res) / (0.64 * rho_bat);
end
