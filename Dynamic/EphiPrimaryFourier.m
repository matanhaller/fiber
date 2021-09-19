% Fourier transform of the phi component of the primary electric field.

function EphiFourier = EphiPrimaryFourier(rho, n, kz, omega, Ez, eta0HzDeriv)
    kVacSquared = kz.^2 - omega.^2;
    
    % Using Maxwell's equations
    EphiFourier = (n.*kz./rho .* Ez - 1j*omega.*eta0HzDeriv) ./ kVacSquared;
end