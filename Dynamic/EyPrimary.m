% Primary electric field in the y direction of a point charge moving at a
% constant velocity in the x direction

function Ey = EyPrimary(x, y, z, t, x0, y0, z0, beta)
    gamma = (1 - beta^2) ^ (-0.5);
    Ey = -gamma * (y - y0) ./ (gamma^2 * (x - x0 - beta*t).^2 + (y - y0).^2 + (z - z0).^2).^(1.5);
end