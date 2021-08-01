% Fourier transform of the rho component of the primary electric field.

function ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta)
    kVacSquared = kz.^2 - omega.^2;
    
    % Using Maxwell's equations
    ErhoFourier = ((1j * kz .* EzPrimaryFourierDeriv(rho, n, kz, omega, x0, y0, z0, beta)) ...
        +(omega.*n./rho).*eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta)) ./ kVacSquared;
end