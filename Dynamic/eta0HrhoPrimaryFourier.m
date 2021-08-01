% Fourier transform of the rho component of the primary magnetic field.

function eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta)
    kVacSquared = kz.^2 - omega.^2;
    
    % Using Maxwell's equations
    eta0HrhoFourier = ((-omega.*n./rho) .* EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta) ...
        +1j*kz.*eta0HzPrimaryFourierDeriv(rho, n, kz, omega, x0, y0, z0, beta)) ./ kVacSquared;
end