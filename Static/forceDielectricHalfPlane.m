% Force due to a dielectric half-plane. Proportional to
% (eta-1)^-2 and is the asymptotic limit of the exact force due
% to the fiber as eta-->1.

function F0 = forceDielectricHalfPlane(epsilonR, eta)
    F0 = 0.25 * (epsilonR - 1) / (epsilonR + 1) * (eta - 1) .^ (-2);
end