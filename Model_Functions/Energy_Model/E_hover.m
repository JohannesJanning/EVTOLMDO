%% Function for Hover Energy (Wh)
function E = E_hover(P_hover, t_hover)
    E = P_hover * t_hover / 3600;
end