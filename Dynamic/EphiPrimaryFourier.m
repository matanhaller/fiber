% Fourier transform of the phi component of the primary electric field.

function EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta)
    kVacSquared = kz.^2 - omega.^2;
    
    % Using Maxwell's equations
    EphiFourier = (n.*kz./rho .* EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta) ...
        - 1j*omega.*eta0HzPrimaryFourierDeriv(rho, n, kz, omega, x0, y0, z0, beta)) ./ kVacSquared;
end