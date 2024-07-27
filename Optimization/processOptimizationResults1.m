function results = processOptimizationResults(x, fval, params, set_trip_distance)
    results = [];
    for i = 1:size(x, 1)
        MTOM = x(i, 1);
        R_prop_cruise = x(i, 2);
        R_prop_hover = x(i, 3);
        c = x(i, 4);
        b = x(i, 5);
        rho_bat = x(i, 6);

        % Use the fixed parameters
        alpha_deg_cruise = params.alpha_deg_cruise;
        alpha_deg_climb = params.alpha_deg_climb;
        num_props_cruise = params.num_props_cruise;
        num_props_hover = params.num_props_hover;
        g = params.g;
        rho_hover = params.rho_hover;
        rho = params.rho;
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
        eta_h = params.eta_h;
        eta_c = params.eta_c;
        m_crew = params.m_crew;
        m_other = params.m_other;
        time_weight = params.time_weight;
        co2_weight = params.co2_weight;
        energy_weight = params.energy_weight;
        costs_weight = params.costs_weight;

      

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
        m_bat = M_bat(rho_bat, e_trip, e_res);
        m_motor = M_motor(Power_hover, Power_climb);
        m_wing = M_wing(c, b);
        m_pay = M_pay (N_s, M_pax, M_lug);
        m_rotor = M_rotor(num_props_hover, num_props_cruise, R_prop_hover, R_prop_cruise);
        m_empty = m_wing + m_motor + m_rotor + m_crew + m_other;
        omega_empty = m_empty / MTOM;


        % --> Cost processing 

        % Cash Operating Costs
        c_e = C_E(e_trip, P_e);
        c_n = C_N(D_trip);
        c_cob = C_Cob(t_tot);
        c_mwr = C_Mwr(t_tot);
        %supplementary energy processing
        e_tot = E_tot(rho_bat, m_bat); % total battery energy 
        e_tripm = E_trip_max(rho_bat, m_bat, e_res); % max trip energy
        c_mb = C_MB(P_bat_s, e_trip, e_tot, e_tripm);
        c_m = C_M(c_mwr, c_mb);
        coc = COC(c_e, c_m, c_n, c_cob);

        % Cost of Wwnership & Indirect Operating Costs 
        coo = COO(coc, omega_empty, MTOM, P_s_empty, N_wd, N_df);
        ioc = IOC(coc, omega_empty, MTOM, P_s_empty, N_wd, N_df);

        % Total Operating Costs
        toc = TOC(coc, coo, ioc);
        toc_s = TOC_s(toc, N_s, D_trip);


        
        % -> Trnapsortation Mode Comparison processing

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
        results = [results; fval(i, :), x(i, :), V_cruise, V_climb, Power_hover, Power_climb, Power_cruise, m_empty, m_pay, m_bat, m_motor, m_rotor, m_wing, c_l_cruise, c_l_climb, toc_s, toc, coc, c_e, c_m, c_mwr, c_mb, c_n, c_cob, coo, ioc, d_cru, d_cl, eVTOL_FoM];
    end
end
