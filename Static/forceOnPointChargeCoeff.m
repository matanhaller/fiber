% Coefficients of the force on the point charge

function res = forceOnPointChargeCoeff(n, k, epsilonR, eta)
    res = (epsilonR - 1) * k * ((besselkp(n,k*eta) .* besselip(n,k)) .* (besselk(n,k*eta) .* besseli(n,k))) ...
        ./ (epsilonR .* besselip(n,k) .* besselk(n,k) - besseli(n,k) .* besselkp(n,k));
end