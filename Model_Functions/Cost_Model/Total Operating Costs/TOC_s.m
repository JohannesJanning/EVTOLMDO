% Specific Total Operating Costs per seat per km in €/skm
function TOC_s_value = TOC_s(TOC_value, N_s, D_trip_ground_in)
    TOC_s_value = TOC_value / (N_s * D_trip_ground_in);
end
