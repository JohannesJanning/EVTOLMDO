%% Function for max Trip Energy (Wh)
function E = E_trip_max(rho_bat, M_bat, E_res)
    E = 0.64 * rho_bat * M_bat - E_res;
end