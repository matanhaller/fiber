% Fourier transform of the rho component of the primary magnetic field.

function eta0HrhoFourier = eta0HrhoPrimaryFourier(rho, n, kz, omega, Ez, eta0HzDeriv)
    kVacSquared = kz.^2 - omega.^2;
    
    % Using Maxwell's equations
    eta0HrhoFourier = ((-omega.*n./rho) .* Ez + 1j*kz.*eta0HzDeriv) ./ kVacSquared;
end