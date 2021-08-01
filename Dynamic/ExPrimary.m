% Primary electric field in the x direction of a point charge moving at a
% constant velocity in the x direction

function Ex = ExPrimary(x, y, z, t, x0, y0, z0, beta)
    gamma = (1 - beta^2) ^ (-0.5);
    Ex = -gamma * (x - x0 - beta*t) ./ (gamma^2 * (x - x0 - beta*t).^2 + (y - y0).^2 + (z - z0).^2).^(1.5);
end