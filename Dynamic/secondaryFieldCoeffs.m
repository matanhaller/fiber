% Coefficients Ank,..,Dnk of the secondary field.

function [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon)
    % Dispersion relation wavevectors
    kVac = sqrt(kz.^2 - omega.^2); kVac = real(kVac) + 1j*sign(omega).*imag(kVac);    
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    
    % Primary fields
    Ep = EzPrimaryFourier(1, n, kz, omega, x0, y0, z0, beta);
    eta0Hp = eta0HzPrimaryFourier(1, n, kz, omega, x0, y0, z0, beta);
    
    % Effective primary fields
    EpEff = EphiPrimaryFourierEff(1, n, kz, omega, x0, y0, z0, beta, epsilon);
    eta0HpEff = eta0HphiPrimaryFourierEff(1, n, kz, omega, x0, y0, z0, beta, epsilon);
    
    % Matrix elements
    M11 = 1j*omega.*(besselkp(n,kVac)./kVac - epsilon.*besselip(n,kCyl)./besseli(n,kCyl).*besselk(n,kVac)./kCyl);
    M12 = n.*kz.*(1./(kVac.^2) - 1./(kCyl.^2)) .* besselk(n,kVac);
    M21 = M12;
    M22 = -1j*omega.*(besselkp(n,kVac)./kVac - besselip(n,kCyl)./besseli(n,kCyl).*besselk(n,kVac)./kCyl);
    
    % Matrix determinant
    Delta = M11 .* M22 - M12 .* M21;
    
    Bnk = 1./Delta .* (M22.*(-eta0HpEff) - M12.*(-EpEff));
    eta0Dnk = 1./Delta .* (-M21.*(-eta0HpEff) + M11.*(-EpEff));
    
    Ank = (1./besseli(n,kCyl)) .* (Ep + Bnk.*besselk(n,kVac));
    eta0Cnk = (1./besseli(n,kCyl)) .* (eta0Hp + eta0Dnk.*besselk(n,kVac));
    
end