% Effective phi component of the primary magnetic field (times free-space
% wave impedance).

function eta0HphiFourierEff = eta0HphiPrimaryFourierEff(n, kz, omega, epsilon, Ez, eta0Hz, eta0Hphi)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    
    eta0HphiFourierEff = - n.*kz./(kCyl.^2).*eta0Hz ...
        - 1j*epsilon*(omega./kCyl).*(besselip(n,kCyl)./besseli(n,kCyl)).*Ez ...
        + eta0Hphi;
end