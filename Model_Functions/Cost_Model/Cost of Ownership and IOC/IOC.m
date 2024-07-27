% Indirect Operating Costs as described in of the model manual
function IOC_value = IOC(COC, omega_empty, MTOM, P_s_empty, N_wd, N_df)
    IOC_value = 0.233 * COC + 0.0175 * omega_empty * MTOM * P_s_empty / (N_wd * N_df);
end
