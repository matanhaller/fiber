% Effective phi component of the primary electric field.

function EphiFourierEff = EphiPrimaryFourierEff(rho, n, kz, omega, x0, y0, z0, beta, epsilon)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    
    EphiFourierEff = - n.*kz./(kCyl.^2).*EzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta) ...
        + 1j*(omega./kCyl).*(besselip(n,kCyl)./besseli(n,kCyl)).*eta0HzPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta) ...
        + EphiPrimaryFourier(rho, n, kz, omega, x0, y0, z0, beta);
end