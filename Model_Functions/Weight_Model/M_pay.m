% Payload Mass in kg 

function M = M_pay(N_s, M_pax, M_lug)
    % Calculate the total payload mass
    M = (M_pax + M_lug) * N_s;
end
