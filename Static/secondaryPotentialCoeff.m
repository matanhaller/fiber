% Coefficients of the secondary potential induced by the fiber

function phi = secondaryPotentialCoeff(n, k, epsilonR, eta, x, y, z)
    R = 1;
    [theta, r, z] = cart2pol(x, y, z);
    phi = cos(n * (theta - 0.5*pi)) .* cos(k*z) .* (...
        (1/epsilonR) .* AI(n,k,r,eta,epsilonR) .* (1 - heaviside(r - R)) + ...
        BK(n,k,r,eta,epsilonR) .* heaviside(r - R));
end

function res = AI(n, k, r, eta, epsilonR)
    res = epsilonR * (besselk(n,k*eta).*besseli(n,k) + BK(n,k,1,eta,epsilonR)) .* (besseli(n,k*r) ./ (besseli(n,k)));
end

function res = BK(n, k, r, eta, epsilonR)
    res = -(epsilonR - 1) * ((besselk(n,k*eta) .* besselip(n,k)) .* (besselk(n,k*r) .* besseli(n,k))) ...
        ./ (epsilonR .* besselip(n,k) .* besselk(n,k) - besseli(n,k) .* besselkp(n,k));
end