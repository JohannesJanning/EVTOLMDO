function validationResults = recheckModel(rho_bat, MTOM, c, b, V_cruise, V_climb, R_prop_cruise, R_prop_hover, set_trip_distance, params)

    % Fixed parameters
    alpha_deg_cruise = params.alpha_deg_cruise;
    alpha_deg_climb = params.alpha_deg_climb;
    alpha_deg_max = params.alpha_deg_max;
    num_props_cruise = params.num_props_cruise;
    num_props_hover = params.num_props_hover;
    g = params.g;
    rho = params.rho;
    rho_hover = params.rho_hover;
    rho_climb = params.rho_climb;
    h_hover = params.h_hover;
    h_cruise = params.h_cruise;
    N_s = params.N_s;
    M_pax = params.M_pax;
    M_lug = params.M_lug;
    P_s_empty = params.P_s_empty;
    N_wd = params.N_wd;
    N_df = params.N_df;
    P_e = params.P_e;
    P_bat_s = params.P_bat_s;
    t_res = params.t_res;
    t_hover = params.t_hover;
    D_trip = set_trip_distance;
    eta_e = params.eta_e;
    eta_p = params.eta_p;
    eta_hp = params.eta_hp;
    eta_h = params.eta_h;
    eta_c = params.eta_c;



    c_l_cruise = Lift_Coefficient(alpha_deg_cruise, c, b);
    c_l_climb = Lift_Coefficient(alpha_deg_climb, c, b);

    c_d_cruise = Drag_Coefficient(c_l_cruise, c, b);
    c_d_climb = Drag_Coefficient(c_l_climb, c, b);

    drag_cruise = Drag(rho, V_cruise, c, b, c_d_cruise);


    roc = ROC(alpha_deg_climb, V_climb);

    % Power Processing (via Momentum Theory)
    Power_hover = Power_Hover_MT(MTOM, rho_hover, R_prop_hover, num_props_hover, eta_h);
    Power_climb = Power_Climb_MT(MTOM, rho_climb, c, b, V_climb, roc, c_d_climb, R_prop_cruise, num_props_cruise, eta_c, alpha_deg_climb);
    Power_cruise = Power_Cruise_MT(MTOM, rho, drag_cruise, V_cruise, R_prop_cruise, num_props_cruise, eta_c, 0);


    % Time and distance processing
    t_cl = t_climb(h_cruise, h_hover, roc);
    d_cl = D_climb(V_climb, t_cl);
    d_cru = D_cruise(D_trip, d_cl, 0);
    t_cr = t_cruise(d_cru, V_cruise);
    t_dst = 0;
    t_tot = (t_cl + t_hover + 0 + t_cr);

    % Energy processing
    e_hvr = E_hover(Power_hover, t_hover);
    e_cl = E_climb(Power_climb, t_cl);
    e_cru = E_cruise(Power_cruise, t_cr);
    e_dst = E_descent(0, t_dst);
    e_trip = E_trip(e_hvr, e_cl, e_cru, 0);
    e_res = E_res(Power_cruise, t_res);

    % Mass processing
    m_wing = M_wing(c, b);
    m_rotor = M_rotor(num_props_hover, num_props_cruise, R_prop_hover, R_prop_cruise);
    m_motor = M_motor(Power_hover, Power_climb);
    m_crew = 82.5 + 14;
    m_other = 668;
    m_empty = m_wing + m_rotor + m_motor + m_crew + m_other;

    m_bat = M_bat(rho_bat, e_trip, e_res);
    m_pay = M_pay(N_s, M_pax, M_lug);

    omega_empty = m_empty / MTOM;

    % Cost processing
    c_e = C_E(e_trip, P_e);
    c_n = C_N(D_trip);
    c_cob = C_Cob(t_tot);
    c_mwr = C_Mwr(t_tot);

    e_tot = E_tot(rho_bat, m_bat);
    e_tripm = E_trip_max(rho_bat, m_bat, e_res);

    c_mb = C_MB(P_bat_s, e_trip, e_tot, e_tripm);
    c_m = C_M(c_mwr, c_mb);
    coc = COC(c_e, c_m, c_n, c_cob);

    coo = COO(coc, omega_empty, MTOM, P_s_empty, N_wd, N_df);
    ioc = IOC(coc, omega_empty, MTOM, P_s_empty, N_wd, N_df);

    toc = TOC(coc, coo, ioc);
    toc_s = TOC_s(toc, N_s, D_trip);

    % Return results in a structured format
    validationResults = struct('Trip_Energy', e_trip, 'Trip_Time', t_tot, 'V_cruise', V_cruise, 'V_climb', V_climb, 'R_cruise', R_prop_cruise, 'R_hover', R_prop_hover, 'chord', c, 'wingspan_b', b, 'Power_Hover', Power_hover, 'Power_Climb', Power_climb, 'Power_Cruise', Power_cruise, 'MTOM', MTOM, 'm_empty', m_empty, 'm_pay', m_pay, 'm_bat', m_bat, 'm_motor', m_motor, 'm_rotor', m_rotor, 'm_wing', m_wing, 'c_l_cruise', c_l_cruise, 'c_l_climb' , c_l_climb ,'toc_s', toc_s, 'toc', toc, 'coc', coc, 'c_e', c_e, 'c_m', c_m, 'c_mwr', c_mwr, 'c_mb', c_mb, 'c_n', c_n, 'c_cob', c_cob, 'coo', coo, 'ioc', ioc, 'd_cru', d_cru, 'd_cl', d_cl);


end
