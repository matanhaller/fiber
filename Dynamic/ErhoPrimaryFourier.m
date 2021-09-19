% Fourier transform of the rho component of the primary electric field.

function ErhoFourier = ErhoPrimaryFourier(rho, n, kz, omega, EzDeriv, eta0Hz)
    kVacSquared = kz.^2 - omega.^2;
    
    % Using Maxwell's equations
    ErhoFourier = ((1j * kz .* EzDeriv) + (omega.*n./rho).*eta0Hz) ./ kVacSquared;
end