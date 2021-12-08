% Azimuthal Fourier component of the secondary electric field.

function EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2) + eps; kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    kVac = sqrt(kz.^2 - omega.^2) + eps; kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    if rho < 1
        EphiFourier = ((n.*kz./rho).*Ank.*besseli(n,kCyl.*rho) ...
            - 1j*omega.*kCyl.*eta0Cnk.*besselip(n,kCyl.*rho)) ./ (kCyl.^2);
    else
        EphiFourier = ((n.*kz./rho).*Bnk.*besselk(n,kVac.*rho) ...
            - 1j*omega.*kVac.*eta0Dnk.*besselkp(n,kVac.*rho)) ./ (kVac.^2);
    end
   
end