% Sum of product of Bessel functions, used for calculation of the Fourier
% transform of the primary fields.

function bs = besselSum(rho, n, kz, omega, beta, nuMax)
    gamma = (1 - beta^2) ^ (-0.5);
    omegaNorm = omega / (gamma*beta);
    hypot = sqrt(kz.^2 + omegaNorm.^2);
        
    nuVec = -nuMax:nuMax; nuVec = reshape(nuVec, 1, [], numel(nuVec));
    nuMat = ones(size(omega));
    nuMat = repmat(nuMat, 1, 1, numel(nuVec)).* nuVec;
    W = repmat(omega, 1, 1, numel(nuVec));
    H = repmat(hypot, 1, 1, numel(nuVec));
    
    bs = besselj(nuMat, (W/beta) .* rho) .* besseli(-n-nuMat, H.*rho);
    bs = sum(bs, 3);

end