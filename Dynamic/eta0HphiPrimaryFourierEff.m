% Effective phi component of the primary magnetic field (times free-space
% wave impedance).

function eta0HphiFourierEff = eta0HphiPrimaryFourierEff(n, kz, omega, epsilon, Ez, eta0Hz, eta0Hphi, eps)
    kVac = sqrt(kz.^2 - omega.^2) + eps; kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    eta0HphiFourierEff = eta0Hphi ...
        - n.*kz./(kVac.^2).*eta0Hz ...
        - 1j*(omega./kVac).*(besselkp(n,kVac)./besselk(n,kVac)).*Ez;
end