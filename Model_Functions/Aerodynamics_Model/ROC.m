function ROC = ROC_calculation(alpha_deg_climb, V_climb)
    alpha_rad_climb = deg2rad(alpha_deg_climb);
    ROC = tan(alpha_rad_climb) * V_climb;
end