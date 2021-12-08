% Radial Fourier component of the secondary magnetic field.

function eta0HrhoFourier = eta0HrhoSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon, eps)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2) + eps; kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    kVac = sqrt(kz.^2 - omega.^2) + eps; kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    if rho < 1
        eta0HrhoFourier = (-(epsilon*omega.*n./rho).*Ank.*besseli(n,kCyl.*rho) ...
            + 1j.*kz.*kCyl.*eta0Cnk.*besselip(n,kCyl.*rho)) ./ (kCyl.^2);
    else
        eta0HrhoFourier = (-(omega.*n./rho).*Bnk.*besselk(n,kVac.*rho) ...
            + 1j.*kz.*kVac.*eta0Dnk.*besselkp(n,kVac.*rho)) ./ (kVac.^2);
    end
   
end