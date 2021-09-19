% Azimuthal Fourier component of the secondary magnetic field.

function eta0HphiFourier = eta0HphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    kVac = sqrt(kz.^2 - omega.^2); kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    if rho < 1
        eta0HphiFourier = (1j*omega.*epsilon.*kCyl.*Ank.*besselip(n,kCyl.*rho) ...
            + (n.*kz./rho).*eta0Cnk.*besseli(n,kCyl.*rho)) ./ (kCyl.^2);
    else
        eta0HphiFourier = (1j*omega.*kVac.*Bnk.*besselkp(n,kVac.*rho) ...
            + (n.*kz./rho).*eta0Dnk.*besselk(n,kVac.*rho)) ./ (kVac.^2);
    end
   
end