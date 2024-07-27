% Battery replacement costs in € per trip


function C = C_MB(P_bat_s, E_trip, E_tot, E_trip_max)
    C = 0.0024 * P_bat_s * E_tot / 1000 * E_trip / E_trip_max;
end
