% Fourier transform of the phi component of the primary magnetic field
% (times free space wave impedance).

function eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, EzDeriv, eta0Hz, eps)
    kVacSquared = kz.^2 - omega.^2 + eps^2;
    
    % Using Maxwell's equations
    eta0HphiFourier = (n.*kz./rho .* eta0Hz + 1j * omega .* EzDeriv) ./ kVacSquared;
end