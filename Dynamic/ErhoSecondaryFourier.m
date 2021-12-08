% Radial Fourier component of the secondary electric field.

function ErhoFourier = ErhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2) + eps; kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    kVac = sqrt(kz.^2 - omega.^2) + eps; kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    if rho < 1
        ErhoFourier = (1j*kz.*kCyl.*Ank.*besselip(n,kCyl.*rho) ...
            + (omega.*n./rho).*eta0Cnk.*besseli(n,kCyl.*rho)) ./ (kCyl.^2);
    else
        ErhoFourier = (1j*kz.*kVac.*Bnk.*besselkp(n,kVac.*rho) ...
            + (omega.*n./rho).*eta0Dnk.*besselk(n,kVac.*rho)) ./ (kVac.^2);
    end
   
end