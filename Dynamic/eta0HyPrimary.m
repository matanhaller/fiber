% Primary magnetic field (times free space wave impedance) in the y
% direction of a point charge moving at a constant velocity in the x
% direction

function eta0Hy = eta0HyPrimary(x, y, z, t, x0, y0, z0, beta)
    gamma = (1 - beta^2) ^ (-0.5);
    eta0Hy = gamma * beta * (z - z0) ./ (gamma^2 * (x - x0 - beta*t).^2 + (y - y0).^2 + (z - z0).^2).^(1.5);
end