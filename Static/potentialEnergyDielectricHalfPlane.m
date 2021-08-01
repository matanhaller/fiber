% Potential energy due to a dielectric half-plane. Proportional to
% (eta-1)^-1 and is the asymptotic limit of the exact potential energy due
% to the fiber as eta-->1.

function W0 = potentialEnergyDielectricHalfPlane(epsilonR, eta)
    W0 = 0.5 * (epsilonR - 1) / (epsilonR + 1) * (eta - 1).^(-1);
end