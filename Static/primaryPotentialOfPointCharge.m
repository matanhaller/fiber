% Primary potential of point charge (basically Coulomb's Law)

function phi = primaryPotentialOfPointCharge(eta, x, y, z)
    R = 1;
    phi = -1 ./ sqrt(x.^2 + (y - eta).^2 + z.^2) .* heaviside(sqrt(x.^2 + y.^2) - R);
end

