% Fourier trasform of the z component of the primary electric field.

function EzFourier = EzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs)
    gamma = (1 - beta^2) ^ (-0.5);
    omegaNorm = omega / (gamma*beta);
    hypot = sqrt(kz.^2 + omegaNorm.^2);
    
    C = -1j*kz / (2*pi) .* exp(1j*(pi/2)*n)  .* (beta * hypot).^(-1) .* ...
        exp(1j*kz*z0 + 1j*(omega/beta)*x0 - hypot*abs(y0));
      
    EzFourier = C .* bs;
end