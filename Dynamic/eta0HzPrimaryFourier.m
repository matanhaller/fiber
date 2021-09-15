% Fourier trasform of the z component of the primary magnetic field
% (times free space wave impedance).

function eta0HzFourier = eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta)
    gamma = (1 - beta^2) ^ (-0.5);
    omegaNorm = omega / (gamma*beta);
    hypot = sqrt(kz.^2 + omegaNorm.^2);
    
    C = 1 / (2*pi) .* exp(1j*(pi/2)*n) .* exp(1j*kz*z0 + 1j*(omega/beta) * x0 - hypot*abs(y0));
    
    nuVec = -40:40; nuVec = reshape(nuVec, 1, [], numel(nuVec));
    nuMat = ones(size(omega));
    nuMat = repmat(nuMat, 1, 1, numel(nuVec)).* nuVec;
    W = repmat(omega, 1, 1, numel(nuVec));
    H = repmat(hypot, 1, 1, numel(nuVec));
    
    besselSum = besselj(nuMat, (W/beta) .* rho) .* besseli(-n-nuMat, H.*rho);
    besselSum = sum(besselSum, 3);
    
%     besselSum = 0;
%     nuMax = 40;
%     for nu=(-nuMax:nuMax)
%         besselSum = besselSum + besselj(nu, (omega/beta) .* rho) .* besseli(-n-nu, hypot.*rho);
%     end
    
    eta0HzFourier = C .* besselSum;
end