function Weight_motorESC = Weight_motorESC_calculation(p_hover, p_climb);

    Weight_motorESC =  6.1e-4 * p_hover +  6.1e-4 * p_climb; % motor weight, where P_line is the motor power rating 
end