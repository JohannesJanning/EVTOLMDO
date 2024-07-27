%% Function for Climb Energy (Wh) 
function e_climb = E_climb(P_climb, t_climb)
    e_climb = P_climb * t_climb / 3600;
end
