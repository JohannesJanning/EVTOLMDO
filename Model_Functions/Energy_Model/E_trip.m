%% Function for Trip Energy (Wh)
function E_trip = E_trip_calc(E_hover, E_climb, E_cruise, E_descent)
    E_trip = E_hover + E_climb + E_cruise + E_descent;
end