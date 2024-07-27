%% Function for Reserve Energy (Wh)
function E = E_res(P_cruise, t_res)
    E = P_cruise * t_res / 3600;
end