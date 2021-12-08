% Coefficients Ank,..,Dnk of the secondary field, times the determinant.
% (Relevant only at the poles)

function [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffsTimesDelta(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax, sigma, calcSum)
    % Dispersion relation wavevectors
    kVac = sqrt(kz.^2 - omega.^2); kVac = real(kVac) + 1j*sign(omega).*imag(kVac);    
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    
    % Bessel function sums
    bs = besselSum(1, n, kz, omega, beta, nuMax, calcSum);
    bsDeriv = besselSumDeriv(1, n, kz, omega, beta, nuMax, calcSum); 
    
    % Primary fields
    Ezp = EzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
    eta0Hzp = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
    EzpDeriv = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
    eta0HzpDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
    Ephip = EphiPrimaryFourier(1, n, kz, omega, Ezp, eta0HzpDeriv);
    eta0Hphip = eta0HphiPrimaryFourier(1, n, kz, omega, EzpDeriv, eta0Hzp);
    
    % Effective primary fields
    EpEff = EphiPrimaryFourierEff(n, kz, omega, epsilon, Ezp, Ephip, eta0Hzp);
    eta0HpEff = eta0HphiPrimaryFourierEff(n, kz, omega, epsilon, Ezp, eta0Hzp, eta0Hphip);
    
    % Matrix elements
    M11 = 1j*omega.*(besselkp(n,kVac)./kVac - epsilon.*besselip(n,kCyl)./besseli(n,kCyl).*besselk(n,kVac)./kCyl);
    M12 = n.*kz.*(1./(kVac.^2) - 1./(kCyl.^2)) .* besselk(n,kVac);
    M21 = M12;
    M22 = -1j*omega.*(besselkp(n,kVac)./kVac - besselip(n,kCyl)./besseli(n,kCyl).*besselk(n,kVac)./kCyl);
    
    Bnk = M22.*(-eta0HpEff) - M12.*(-EpEff);
    eta0Dnk = -M21.*(-eta0HpEff) + M11.*(-EpEff);
    
    Ank = Bnk.*besselk(n,kVac)./besseli(n,kCyl);
    eta0Cnk = eta0Dnk.*besselk(n,kVac)./besseli(n,kCyl);
    
end