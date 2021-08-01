% Coefficients of the Intensity of the Radial Part of the E
% field

function res = intensityCoeffRadial(n, k, epsilonR, eta, x, y, z)
    [theta, r, z] = cart2pol(x,y,z);
    res = cos(n * (theta - 0.5*pi)) .* cos(k*z) .* (...
         AI(n,k,r,eta,epsilonR));
end

function res = AI(n, k, r, eta, epsilonR)
    res = besselk(n, k*eta) .* ((besselip(n,k) .* besselk(n,k) - besseli(n,k) .* besselkp(n,k))...
        ./(epsilonR .* besselip(n,k) .* besselk(n,k) - besseli(n,k) .* besselkp(n,k))) .* k .* besselip(n, k * r);
end