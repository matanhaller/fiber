% Effective phi component of the primary electric field.

function EphiFourierEff = EphiPrimaryFourierEff(n, kz, omega, epsilon, Ez, Ephi, eta0Hz)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    
    EphiFourierEff = - n.*kz./(kCyl.^2).*Ez ...
        + 1j*(omega./kCyl).*(besselip(n,kCyl)./besseli(n,kCyl)).*eta0Hz ...
        + Ephi;
end