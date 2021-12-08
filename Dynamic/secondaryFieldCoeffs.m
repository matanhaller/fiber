% Coefficients Ank,..,Dnk of the secondary field.

function [Ank, Bnk, eta0Cnk, eta0Dnk] = secondaryFieldCoeffs(n, kz, omega, x0, y0, z0, beta, epsilon, nuMax, sigma, calcSum, eps)
    % Dispersion relation wavevectors
    kVac = sqrt(kz.^2 - omega.^2) + eps; kVac = real(kVac) + 1j*sign(omega).*imag(kVac);    
    kCyl = sqrt(kz.^2 - epsilon * omega.^2) + eps; kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    
    % Bessel function sums
    bs = besselSum(1, n, kz, omega, beta, nuMax, calcSum);
    bsDeriv = besselSumDeriv(1, n, kz, omega, beta, nuMax, calcSum); 
    
    % Primary fields
    Ezp = EzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
    eta0Hzp = eta0HzPrimaryFourier(n, kz, omega, x0, y0, z0, beta, bs, sigma);
    EzpDeriv = EzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
    eta0HzpDeriv = eta0HzPrimaryFourierDeriv(n, kz, omega, x0, y0, z0, beta, bsDeriv, sigma);
    Ephip = EphiPrimaryFourier(1, n, kz, omega, Ezp, eta0HzpDeriv, eps);
    eta0Hphip = eta0HphiPrimaryFourier(1, n, kz, omega, EzpDeriv, eta0Hzp, eps);

    % Effective primary fields
    EpEff = EphiPrimaryFourierEff(n, kz, omega, epsilon, Ezp, Ephip, eta0Hzp, eps);
    eta0HpEff = eta0HphiPrimaryFourierEff(n, kz, omega, epsilon, Ezp, eta0Hzp, eta0Hphip, eps);
    
    % Modified Bessel function values on the cylinder's boundary
    I = besseli(n, kCyl);
    Ip = besselip(n, kCyl);
    K = besselk(n, kVac);
    Kp = besselkp(n, kVac);

    % Matrix elements
    M11 = 1j*omega.*(epsilon.*Ip./kCyl - (Kp./K).*I./kVac);
    M12 = n.*kz.*(1./(kCyl.^2) - 1./(kVac.^2)) .* I;
    M21 = M12;
    M22 = -1j*omega.*(Ip./kCyl - (Kp./K).*I./kVac);
    
    % Matrix determinant
    Delta = M11 .* M22 - M12 .* M21;
    
    % Solving equation systems
    Ank = (M22.*eta0HpEff - M12.*EpEff) ./ Delta;
    eta0Cnk = (- M21.*eta0HpEff + M11.*EpEff) ./ Delta;

    Bnk = (Ank.*I - Ezp) ./ K;
    eta0Dnk = (eta0Cnk.*I - eta0Hzp) ./ K;

end