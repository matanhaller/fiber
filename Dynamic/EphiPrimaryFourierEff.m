% Effective phi component of the primary electric field.

function EphiFourierEff = EphiPrimaryFourierEff(n, kz, omega, epsilon, Ez, Ephi, eta0Hz, eps)
    kVac = sqrt(kz.^2 - omega.^2) + eps; kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    EphiFourierEff = Ephi ...
        - n.*kz./(kVac.^2).*Ez ...
        + 1j*(omega./kVac).*(besselkp(n,kVac)./besselk(n,kVac)).*eta0Hz;
end