% MTOM in kg 

function M = M_TOM(M_pay, M_bat, omega_empty)
    % Calculate the Maximum Takeoff Mass (MTOM)
    M = (M_pay + M_bat) / (1 - omega_empty);
end
