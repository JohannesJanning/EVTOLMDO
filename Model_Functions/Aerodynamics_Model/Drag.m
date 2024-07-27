function D = Drag(rho, V, c, b, C_D)
    S = c * b;
    D = rho / 2 * V^2 * S * C_D;
end