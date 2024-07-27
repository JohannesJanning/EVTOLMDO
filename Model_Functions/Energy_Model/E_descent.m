%% Function for Descent Energy (Wh)
function E = E_descent(P_descent, t_descent)
    E = P_descent * t_descent / 3600;
end