%% Function to calculate the time for the trip (s)
function t_tr = t_trip(t_cruise, t_climb, t_descent, t_hover)
  t_tr = (t_cruise+t_climb+t_descent+t_hover);
 end