% Effective phi component of the primary magnetic field (times free-space
% wave impedance).

function eta0HphiFourierEff = eta0HphiPrimaryFourierEff(rho, n, kz, omega, x0, y0, z0, beta, epsilon)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    
    eta0HphiFourierEff = - n.*kz./(kCyl.^2).*eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta) ...
        - 1j*epsilon*(omega./kCyl).*(besselip(n,kCyl)./besseli(n,kCyl)).*EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta) ...
        + eta0HphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
end