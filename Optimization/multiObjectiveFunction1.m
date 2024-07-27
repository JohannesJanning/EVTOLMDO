function f = multiObjectiveFunction(x, params, set_trip_distance)
    % Extracting design variables
    MTOM = x(1);
    R_prop_cruise = x(2);
    R_prop_hover = x(3);
    c = x(4);
    b = x(5);
    rho_bat = x(6);


    % Using the fixed parameters
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
    t_hover = params.t_hover;
    D_trip = set_trip_distance;
    eta_h = params.eta_h;
    eta_c = params.eta_c;
    
    
    % Aerodynamic & Performance Processing 
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

    % Time and distance Processing 
    t_cl = t_climb(h_cruise, h_hover, roc);
    d_cl = D_climb(V_climb, t_cl);
    d_cru = D_cruise(D_trip, d_cl, 0);
    t_cr = t_cruise(d_cru, V_cruise);
    t_tot = (t_cl + t_hover + t_cr); 

    % Energy processing 
    e_hvr = E_hover(Power_hover, t_hover);
    e_cl = E_climb(Power_climb, t_cl);
    e_cru = E_cruise(Power_cruise, t_cr);
    e_trip = E_trip(e_hvr, e_cl, e_cru, 0);

    % Define objective functions
    f = [e_trip, t_tot]; % Minimizing energy trip and trip time
end