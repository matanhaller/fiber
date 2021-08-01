% Coefficients of the potential energy on the point charge in the limit of
% a conducting cylinder.

function res = potentialEnergyOnPointChargeConductingCoeff(n, k, epsilonR, eta)
    res = (besseli(n,k) .* besselk(n,k*eta)) .* (besselk(n,k*eta) ./ besselk(n,k));
end