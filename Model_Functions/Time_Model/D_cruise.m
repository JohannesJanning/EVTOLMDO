%% Function to calculate the distance for cruise (m)
function d_cr = D_cruise(D_trip_ground_in, D_climb, D_descent)
  d_cr = (D_trip_ground_in*1000-D_climb-D_descent);
end