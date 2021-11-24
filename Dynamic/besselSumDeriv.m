% Derivative of sum of product of Bessel functions.

function bsDeriv = besselSumDeriv(rho, n, kz, omega, beta, nuMax, calcSum)
    gamma = (1 - beta^2) ^ (-0.5);
    omegaNorm = omega / (gamma*beta);
    hypot = sqrt(kz.^2 + omegaNorm.^2);
    if calcSum == 1
        nuVec = -nuMax:nuMax; nuVec = reshape(nuVec, 1, [], numel(nuVec));
        nuMat = ones(size(omega));
        nuMat = repmat(nuMat, 1, 1, numel(nuVec)).* nuVec;
        W = repmat(omega, 1, 1, numel(nuVec));
        H = repmat(hypot, 1, 1, numel(nuVec));
        
        bsDeriv = ((W/beta) .* besseljp(nuMat, (W/beta).* rho) .* besseli(-n-nuMat, H.*rho) ...
                 + H .* besselj(nuMat, (W/beta).* rho) .* besselip(-n-nuMat, H.*rho));
        bsDeriv = sum(bsDeriv, 3);
    else
        phiVec = linspace(-pi, pi, nuMax); phiVec = reshape(phiVec, 1, [], numel(phiVec));
        phiMat = ones(size(omega));
        phiMat = repmat(phiMat, 1, 1, numel(phiVec)) .* phiVec;
        W = repmat(omega, 1, 1, numel(phiVec));
        H = repmat(hypot, 1, 1, numel(phiVec));

        bsDerivCoeff = exp(1j*n*phiMat) .* (-1j*W/beta.*cos(phiMat) + hypot.*sin(phiMat)) ...
            .* exp(-1j*W/beta*rho.*cos(phiMat)) .* exp(hypot*rho.*sin(phiMat));
        bsDeriv = 1 / (2*pi) * exp(-1j*(pi/2)*n) * trapz(phiVec(:), bsDerivCoeff, 3);
    end
end