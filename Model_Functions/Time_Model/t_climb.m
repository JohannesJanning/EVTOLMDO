%% Function to calculate the time in climb (s)
function t_cl = t_climb(h_cruise, h_hover, ROC)
  t_cl = (h_cruise-h_hover)/ROC;
end 