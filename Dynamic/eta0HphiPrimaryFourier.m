% Fourier transform of the phi component of the primary magnetic field
% (times free space wave impedance).

function eta0HphiFourier = eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta)
    kVacSquared = kz.^2 - omega.^2;
    
    % Using Maxwell's equations
    eta0HphiFourier = (n.*kz./rho .* eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta) ...
        +1j * omega .* EzPrimaryFourierDeriv(rho, n, kz, omega, x0, y0, z0, beta)) ./ kVacSquared;
end