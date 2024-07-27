% Cost of Ownership per trip in € 


function COO_value = COO(COC, omega_empty, M_TOM, P_s_empty, N_wd, N_df)
    COO_value = 0.06 * COC + 0.0796 * omega_empty * M_TOM * P_s_empty / (N_wd * N_df);
end
