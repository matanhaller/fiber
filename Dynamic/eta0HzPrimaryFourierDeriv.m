% Derivative of Fourier transform of the z component of the primary magnetic field
% (times free space wave impedance).

function eta0HzFourierDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma)
    gamma = (1 - beta^2) ^ (-0.5);
    omegaNorm = omega / (gamma*beta);
    hypot = sqrt(kz.^2 + omegaNorm.^2);
    
    C = 1 / (2*pi) .* exp(1j*(pi/2)*n) .* exp(1j*kz*z0 + 1j*(omega/beta) * x0 - hypot*abs(y0) ...
        + 0.5*sigma^2*((-(omega/beta).^2 - kz.^2)));

    eta0HzFourierDeriv = C .* bsDeriv;
end