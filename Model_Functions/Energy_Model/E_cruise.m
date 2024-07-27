%% Function for Cruise Energy (Wh)
function E = E_cruise(P_cruise, t_cruise)
    E = P_cruise * t_cruise / 3600;
end