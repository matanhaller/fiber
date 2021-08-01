% Fourier trasform of the z component of the primary electric field.

function EzFourier = EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta)
    gamma = (1 - beta^2) ^ (-0.5);
    omegaNorm = omega / (gamma*beta);
    hypot = sqrt(kz.^2 + omegaNorm.^2);
    
    C = -1j*kz / (2*pi) .* exp(1j*(pi/2)*n)  .* (beta * hypot).^(-1) .* ...
        exp(1j*kz*z0 + 1j*(omega/beta)*x0 - hypot*abs(y0));
    
    besselSum = 0;
    nuMax = 40;
    for nu=(-nuMax:nuMax)
        besselSum = besselSum + besselj(nu, (omega/beta) .* rho) .* besseli(-n-nu, hypot.*rho);
    end
    
    EzFourier = C .* besselSum;
end