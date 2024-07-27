function [c, ceq] = nonlinearConstraints(x, params, set_trip_distance)
    % Extract design variables
    MTOM = x(1);
    R_prop_cruise = x(2);
    R_prop_hover = x(3);
    c = x(4);
    b = x(5);
    rho_bat = x(6);

    % Use the fixed parameters
    alpha_deg_cruise = params.alpha_deg_cruise;
    alpha_deg_climb = params.alpha_deg_climb;
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
    t_res = params.t_res;
    t_hover = params.t_hover;
    D_trip = set_trip_distance;
    eta_h = params.eta_h;
    eta_c = params.eta_c;
    m_crew = params.m_crew;
    m_other = params.m_other; 
    d_tip = params.d_tip;
    w_fuselage = params.w_fuselage;


    % --> Disciplinary Analysis:
    
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


    % -> Constraint Definition:

    % not in use: lift_cruise_constraint = Lift(rho, 0.95*V_cruise, c_l_cruise, c, b) - 1 * MTOM * g; % equlibrium constraint horizontal flight
    % not in use: lift_climb_constraint = Lift(rho_climb, 0.95*V_climb, c_l_climb, c, b) - 1 * MTOM * g * cos(atan(roc/V_climb)); % equilibrium constraint climb flight 
    geometry_constraint = 2 * (num_props_hover / 4 * 2 * R_prop_hover * 3/4 + 2 * d_tip + w_fuselage / 2) - b; % constraint on lifting rotor and fuselage distance
    vertiport_constraint1 = 2 * (2 * d_tip + 4 * R_prop_hover + w_fuselage / 2) - 18 * 1.1; % 18 m vertiport constraint (according to EASA PTS VTP) with 10% margin
    vertiport_constraint2 = b - 18 * 1.1; %18 m vertiport constraint (according to EASA PTS VTP) with 10% margin
    % not in use: MTOM_constraint1 = MTOM - 5700; % SC-VTOL EASA max. certifying take-off mass (kg)
    % not in use: cruise_constraint = sqrt( MTOM * g * 2 / (rho * c_l_cruise * c * b)) -  1.05 * V_cruise; % below 1.3 factor, the computation consumes long time and the results are unrealsitic
    % not in use: climb_constraint = sqrt(cosd(alpha_deg_climb) * MTOM * g *2 /(rho_climb*c_l_climb*c*b)) -  1.05 * V_climb;
    MTOM_constraint2 = m_empty + m_bat + m_pay - 1.04*MTOM;

    ceq = [];
    c = [
         geometry_constraint;
         vertiport_constraint1;
         vertiport_constraint2;
         MTOM_constraint2];
end