function results = processOptimizationResults(x, fval, params, set_trip_distance)
    results = [];
    for i = 1:size(x, 1)
        MTOM = x(i, 1);
        R_prop_cruise = x(i, 2);
        R_prop_hover = x(i, 3);
        c = x(i, 4);
        b = x(i, 5);
        % rho_bat = x(i, 6);

        % Use the fixed parameters
        alpha_deg_cruise = params.alpha_deg_cruise;
        alpha_deg_climb = params.alpha_deg_climb;
        num_props_cruise = params.num_props_cruise;
        num_props_hover = params.num_props_hover;
        g = params.g;
        rho_hover = params.rho_hover;
        rho = params.rho;
        rho_msl = params.rho_msl;
        rho_climb = params.rho_climb;
        h_hover = params.h_hover;
        h_cruise = params.h_cruise;
        N_s = params.N_s;
        M_pax = params.M_pax;
        M_lug = params.M_lug;
        P_s_empty = params.P_s_empty;
        N_wd = params.N_wd;
        T_D = params.T_D;
        N_df = params.N_df;
        P_e = params.P_e;
        P_bat_s = params.P_bat_s;
        t_res = params.t_res;
        t_hover = params.t_hover;
        D_trip = set_trip_distance;
        eta_h = params.eta_h;
        eta_c = params.eta_c;
        m_crew = params.m_crew;
        time_weight = params.time_weight;
        co2_weight = params.co2_weight;
        energy_weight = params.energy_weight;
        costs_weight = params.costs_weight;
        l_fus_m = params.l_fus_m;
        r_fus_m = params.r_fus_m;
        rho_bat = params.rho_bat;
        C_charge = params.C_charge;
        U_pilot = params.U_pilot;
        N_AC = params.N_AC;
        T_D = params.T_D; 
        S_P = params.S_P;
        GWP_battery = params.GWP_battery;
        GWP_energy = params.GWP_energy;
        LF = params.LF;
        fare_km = params.fare_km;

        % Disciplinary Analysis
    
        % Aerodynamic Processing 
        c_l_cruise = Lift_Coefficient(alpha_deg_cruise, c, b);
        c_l_climb = Lift_Coefficient(alpha_deg_climb, c, b);
        c_d_cruise = Drag_Coefficient(c_l_cruise, c, b);
        c_d_climb = Drag_Coefficient(c_l_climb, c, b);
        V_cruise = sqrt(2*MTOM*g/(c_l_cruise*c*b*rho));
        V_climb = sqrt(2*MTOM*g/(c_l_climb*c*b*rho));
        drag_cruise = Drag(rho, V_cruise, c, b, c_d_cruise);
        roc = ROC(alpha_deg_climb, V_climb); 

        AR = b / c; % aspect ratio 
        WL = MTOM / (b * c); % wing loading in kg / m^2
        DL = MTOM / (num_props_hover * pi()*R_prop_hover^2); % disk loading in kg/m^2
        LDcruise = c_l_cruise / c_d_cruise; % LD ratio in cruise 
        LDclimb = c_l_climb / c_d_climb; % LD ratio in climb 

        % Power Processing (via Momentum Theory)
        Power_hover = Power_Hover_MT(MTOM, rho_hover, R_prop_hover, num_props_hover, eta_h);
        Power_climb = Power_Climb_MT(MTOM, rho_climb, c, b, V_climb, roc, c_d_climb, R_prop_cruise, num_props_cruise, eta_c, alpha_deg_climb);
        Power_cruise = Power_Cruise_MT(MTOM, rho, drag_cruise, V_cruise, R_prop_cruise, num_props_cruise, eta_c, 0);

        % Time and distance processing 
        t_cl = t_climb(h_cruise, h_hover, roc);
        d_cl = D_climb(V_climb, t_cl);
        d_cru = D_cruise(D_trip, d_cl, 0);
        t_cr = t_cruise(d_cru, V_cruise);
        t_tot = t_cl + t_hover + 0 + t_cr; 

        % Energy processing 
        e_hvr = E_hover(Power_hover, t_hover);
        e_cl = E_climb(Power_climb, t_cl);
        e_cru = E_cruise(Power_cruise, t_cr);
        e_trip = E_trip(e_hvr, e_cl, e_cru, 0);
        e_res = E_res(Power_cruise, t_res);        


        % Mass processing
        m_furnish = M_furnish(MTOM, V_cruise, rho);
        m_gear = M_gear(MTOM, R_prop_cruise, r_fus_m);
        m_fuselage = M_fuselage(l_fus_m, r_fus_m, MTOM, rho, rho_msl, V_cruise);
        m_system = M_system(l_fus_m, b, MTOM);
        m_bat = M_bat(rho_bat, e_trip, e_res);
        m_motor = M_motor(R_prop_hover, Power_hover, R_prop_cruise, Power_climb);
        m_wing = M_wing(c, b, V_cruise, rho, MTOM);
        m_pay = M_pay (N_s, M_pax, M_lug);
        m_rotor = M_rotor(num_props_hover, num_props_cruise, R_prop_hover, R_prop_cruise);
        m_empty = m_wing + m_motor + m_rotor + m_crew + m_furnish + m_fuselage + m_system + m_gear;
        omega_empty = m_empty / MTOM;

        e_tot = E_tot(rho_bat, m_bat); % total battery energy 
        e_tripm = E_trip_max(rho_bat, m_bat, e_res); % max trip energy


        % Battery and Operations Design Module

        E_battery = rho_bat * m_bat; % total energy capcity of the battery in Wh
        C_rate_hover = Power_hover / E_battery;
        C_rate_climb = Power_climb / E_battery;
        C_rate_cruise = Power_cruise / E_battery;
        C_rate_average = (C_rate_hover * t_hover + C_rate_climb * t_cl + C_rate_cruise * t_cr) / t_tot; % c-rate discharge 
        DOD = e_trip / E_battery; % depth of discharge 
     
        T_T = 1 / C_charge * DOD * 3600;  % [seconds] turnaround time 
        DH = T_T / t_tot + 1; % deadhead ratio 

        n_cycles = N_cycles(DOD, C_rate_average, C_charge); % battery usage cycles, before replacement (degradation to 80% BOL)
        n_bat_req = 1/n_cycles * N_wd * T_D / (t_tot*DH) * e_trip / e_tripm; % number of batteries required to fullfill number of flight ; 416.67
        
        FC_d = T_D / (t_tot * DH); % daly flight cycles 
        FC_a = N_wd * FC_d; % flight cycles flown per year
        FH_d = FC_d * t_tot / 3600; % daily flight hours 
        FH_a = FC_a * t_tot / 3600; % annual flight hours 
        P_bat_single = P_bat_s * e_tot/1000; % single battery price 
        P_bat_annual = n_bat_req * P_bat_single; % yearly battery price        
        


        % --> Cost processing 

        % Cash Operating Costs
        c_e = C_E(e_trip, P_e);
        c_n = C_N(MTOM ,D_trip);
        c_cob = C_Cob(N_wd, T_D, U_pilot, N_AC, t_tot, DH, S_P);
        c_mwr = C_Mwr(t_tot);
        %supplementary energy processing

        c_mb = C_MB(P_bat_s, e_tot, n_bat_req, t_tot, DH, T_D, N_wd);
        c_m = C_M(c_mwr, c_mb);
        coc = COC(c_e, c_m, c_n, c_cob);
           
        % Cost of Wwnership & Indirect Operating Costs 
        coo = COO(coc, omega_empty, MTOM, P_s_empty, N_wd, T_D, t_tot, DH);
        ioc = IOC(coc, omega_empty, MTOM, P_s_empty, N_wd, T_D, t_tot, DH);

        % Total Operating Costs
        toc = TOC(coc, coo, ioc);
        toc_s = TOC_s(toc, N_s, D_trip);


        % Economics & Environment

        GWP_flight = GWP_cycle(E_battery, n_bat_req, FC_a, e_trip, GWP_battery, GWP_energy);
        GWP_annual = GWP_flight * FC_a;
        GWP_energy_fraction = (e_trip/1000 * GWP_energy) / GWP_flight;  
        


        revenue = fare_km * D_trip * N_s * LF; % estimated revenue per trip in €
        ticket_price_pax = revenue / (N_s*LF); % estimated ticket price in €
        Profit_flight = (revenue-toc); % profit per flight in €
        Profit_annual = (revenue-toc) * FC_a; % annual profit in €
    
    
        %%% Efficiency Metrices
        DOC = coc + coo; % € - direct operating cost per flight 
        
        ECE = DOC / (e_trip/1000); % €/kWh - Energy Cost Efficiency
        BCE = c_mb / (e_trip/1000); % €/kWh - Battery Cost Efficienty 
        RpK = revenue / D_trip; % revenue per km ratio
        PR = Profit_flight / toc; % profit ratio 
        
        EpK = e_trip/1000 / D_trip; % kWh/km
        GWPpK = GWP_flight / D_trip; % kg CO2e / km 
        CEE = GWP_flight / (e_trip/1000); % Carbon Efficiency of Energy (CEE)
   
        AC_cost = m_empty * P_s_empty; % aircraft acquisition cost €

        CSI = DOC * GWP_flight / D_trip; % Combined Sustainability Index
        
        COI = (toc * GWP_flight) / (D_trip * Profit_annual/1000); % combined operation index


        % -> Transportation Mode Comparison processing

        % Calculate eVTOL specific parameters
        time_requirement = t_tot/60; % trip time in min 
        eVTOL_CO2_kg_skm = e_trip / 1000 / 4 / D_trip * 0.378; % kg CO2 per seat km
        eVTOL_energy_consumption = e_trip / 4 / D_trip; % Wh/km
        eVTOL_costs_eur_skm = toc_s; % EUR/seat km

        % Transportation mode comparison data 
        modes = {'Airplane', 'Gasonline Vehicle (1)', 'Diesel Vehicle (1)', ...
                 'Electric Vehicle (1)', 'Gasonline Vehicle (5)', 'Diesel Vehicle (5)', ...
                 'Public Bus (100%)', 'Electric Vehicle (5)', 'Train (100%)', ...
                 'Bicycle', 'Airplane (79.6%)', 'Diesel Vehicle (1.3)', ...
                 'Electric Vehicle (1.3)', 'Gasonline Vehicle (1.3)', ...
                 'Public Bus (60%)', 'Train (50%)', 'eVTOL'};

        CO2_kg_skm = [0.198, 0.157, 0.128, 0.065, 0.031, 0.026, 0.013, 0.013, 0.007, 0.000, ...
                      0.249, 0.099, 0.050, 0.120, 0.022, 0.012, eVTOL_CO2_kg_skm];
        energy_consumption = [216.06, 632.40, 480.00, 172.70, 126.48, 96.00, 49.42, 34.54, 57.02, -111.11, ...
                              340.99, 369.23, 132.85, 486.46, 82.84, 114.04, eVTOL_energy_consumption];
        costs_eur_skm = [0.46, 0.117, 0.083, 0.105, 0.023, 0.017, 0.060, 0.021, 0.200, -0.491, ...
                         0.579, 0.064, 0.081, 0.090, 0.104, 0.402, eVTOL_costs_eur_skm];

        velocities = struct('Cars_up_to_60km', 60.00, ...
                            'Cars_above_60km', 85.00, ...
                            'Airplane_below_400km', 74.00, ...
                            'Airplane_above_400km', 151.00, ...
                            'Bicycle', 18.80, ...
                            'Public_Bus_up_to_60km', 39.70, ...
                            'Public_Bus_above_60km', 64.00, ...
                            'Train_up_to_60km', 49.10, ...
                            'Train_above_60km', 99.00);

        circuity_ratios = struct('Cars_up_to_180km', 1.30, ...
                                 'Cars_above_180km', 1.20, ...
                                 'Train_all', 1.20, ...
                                 'Airplane_all', 1.05, ...
                                 'Bicycle_all', 1.28, ...
                                 'Bus_up_to_100km', 1.60, ...
                                 'Bus_above_100km', 1.25);

        adjusted_trip_distance = zeros(size(modes));
        for j = 1:length(modes)
            mode = modes{j};
            if contains(mode, 'Gasonline Vehicle') || contains(mode, 'Diesel Vehicle') || contains(mode, 'Electric Vehicle')
                if D_trip <= 180
                    adjusted_trip_distance(j) = D_trip * circuity_ratios.Cars_up_to_180km;
                else
                    adjusted_trip_distance(j) = D_trip * circuity_ratios.Cars_above_180km;
                end
            elseif contains(mode, 'Airplane')
                adjusted_trip_distance(j) = D_trip * circuity_ratios.Airplane_all;
            elseif strcmp(mode, 'Bicycle')
                adjusted_trip_distance(j) = D_trip * circuity_ratios.Bicycle_all;
            elseif contains(mode, 'Public Bus')
                if D_trip <= 100
                    adjusted_trip_distance(j) = D_trip * circuity_ratios.Bus_up_to_100km;
                else
                    adjusted_trip_distance(j) = D_trip * circuity_ratios.Bus_above_100km;
                end
            elseif contains(mode, 'Train')
                adjusted_trip_distance(j) = D_trip * circuity_ratios.Train_all;
            elseif strcmp(mode, 'eVTOL')
                % eVTOL does not have a circuity ratio adjustment
                adjusted_trip_distance(j) = D_trip;
            end
        end

        CO2_total = CO2_kg_skm .* adjusted_trip_distance;
        energy_total = energy_consumption .* adjusted_trip_distance;
        costs_total = costs_eur_skm .* adjusted_trip_distance;

        CO2_total(end) = eVTOL_CO2_kg_skm * adjusted_trip_distance(end);
        energy_total(end) = eVTOL_energy_consumption * adjusted_trip_distance(end);
        costs_total(end) = eVTOL_costs_eur_skm * adjusted_trip_distance(end);

        % Determine time demand based on average velocity and adjusted trip distance
        time_demand = zeros(size(modes));
        for j = 1:length(modes)
            mode = modes{j};
            if contains(mode, 'Gasonline Vehicle') || contains(mode, 'Diesel Vehicle') || contains(mode, 'Electric Vehicle')
                if D_trip <= 60
                    time_demand(j) = adjusted_trip_distance(j) / velocities.Cars_up_to_60km * 60; % converting hours to minutes
                else
                    time_demand(j) = adjusted_trip_distance(j) / velocities.Cars_above_60km * 60; % converting hours to minutes
                end
            elseif contains(mode, 'Airplane')
                if D_trip <= 400
                    time_demand(j) = (adjusted_trip_distance(j) / velocities.Airplane_below_400km * 60) + 120; % converting hours to minutes and adding 2h airport time
                else
                    time_demand(j) = (adjusted_trip_distance(j) / velocities.Airplane_above_400km * 60) + 120; % converting hours to minutes and adding 2h airport time
                end
            elseif strcmp(mode, 'Bicycle')
                time_demand(j) = adjusted_trip_distance(j) / velocities.Bicycle * 60; % converting hours to minutes
            elseif contains(mode, 'Public Bus')
                if D_trip <= 60
                    time_demand(j) = adjusted_trip_distance(j) / velocities.Public_Bus_up_to_60km * 60; % converting hours to minutes
                else
                    time_demand(j) = adjusted_trip_distance(j) / velocities.Public_Bus_above_60km * 60; % converting hours to minutes
                end
            elseif contains(mode, 'Train')
                if D_trip <= 60
                    time_demand(j) = adjusted_trip_distance(j) / velocities.Train_up_to_60km * 60; % converting hours to minutes
                else
                    time_demand(j) = adjusted_trip_distance(j) / velocities.Train_above_60km * 60; % converting hours to minutes
                end
            elseif strcmp(mode, 'eVTOL')
                % Using eVTOL velocity
                time_demand(j) = time_requirement; % converting hours to minutes
            end
        end

        % Calculate ratings between 1 and 10 for each mode
        calculate_rating = @(x, X_min, X_max) (1 * (x - X_min) - 10 * (x - X_max)) / (X_max - X_min);

        time_rating = calculate_rating(time_demand, min(time_demand), max(time_demand));
        co2_rating = calculate_rating(CO2_total, min(CO2_total), max(CO2_total));
        energy_rating = calculate_rating(energy_total, min(energy_total), max(energy_total));
        costs_rating = calculate_rating(costs_total, min(costs_total), max(costs_total));

        % Calculate the Figure of Merit (FoM) for each mode
        FoM = (time_weight * time_rating') + ...
        (co2_weight * co2_rating') + ...
        (energy_weight * energy_rating') + ...
        (costs_weight * costs_rating');

        % Store the FoM for the eVTOL
        eVTOL_FoM = FoM(end);

        % Append the results
        results = [results; fval(i, :), x(i, :), rho_bat, V_cruise, V_climb, Power_hover, Power_climb, Power_cruise, E_battery/1000, EpK, m_empty, m_pay, m_bat, m_motor, m_rotor, m_wing, m_gear, m_furnish, m_fuselage, m_system, c_l_cruise, c_l_climb, c_d_cruise, c_d_climb, LDcruise, LDclimb, roc*196.85, AR, WL, DL, toc_s, toc, DOC, coc, c_e, c_m, c_mwr, c_mb, c_n, c_cob, coo, ioc, AC_cost, ECE, BCE, n_bat_req, n_cycles, DOD, C_rate_hover, C_rate_climb, C_rate_cruise, C_rate_average, C_charge, T_T, DH, P_bat_single, P_bat_annual, GWP_flight, GWP_annual, GWPpK, CEE, GWP_energy_fraction, CSI, COI, FC_a, FC_d, FH_a, FH_d, Profit_annual, Profit_flight, PR, RpK, ticket_price_pax, d_cru, d_cl, eVTOL_FoM];
    end
end
