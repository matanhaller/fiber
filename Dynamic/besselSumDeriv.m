% Derivative of sum of product of Bessel functions.

function bsDeriv = besselSumDeriv(rho, n, kz, omega, beta, nuMax)
    gamma = (1 - beta^2) ^ (-0.5);
    omegaNorm = omega / (gamma*beta);
    hypot = sqrt(kz.^2 + omegaNorm.^2);
    
    nuVec = -nuMax:nuMax; nuVec = reshape(nuVec, 1, [], numel(nuVec));
    nuMat = ones(size(omega));
    nuMat = repmat(nuMat, 1, 1, numel(nuVec)).* nuVec;
    W = repmat(omega, 1, 1, numel(nuVec));
    H = repmat(hypot, 1, 1, numel(nuVec));
    
    bsDeriv = ((W/beta) .* besseljp(nuMat, (W/beta).* rho) .* besseli(-n-nuMat, H.*rho) ...
             + H .* besselj(nuMat, (W/beta).* rho) .* besselip(-n-nuMat, H.*rho));
    bsDeriv = sum(bsDeriv, 3);
    
end