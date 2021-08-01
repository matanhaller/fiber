% Coefficients of the potential energy on the point charge

function res = potentialEnergyOnPointChargeCoeff(n, k, epsilonR, eta)
    res = (epsilonR - 1) * ((besselk(n,k*eta) .* besselip(n,k)) .* (besselk(n,k*eta) .* besseli(n,k))) ...
        ./ (epsilonR * besselip(n,k) .* besselk(n,k) - besseli(n,k) .* besselkp(n,k));
end