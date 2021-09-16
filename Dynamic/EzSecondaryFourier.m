% Longitudinal Fourier component of the secondary electric field.

function EzFourier = EzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, beta, epsilon)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    kVac = sqrt(kz.^2 - omega.^2); kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    if rho < 1
        EzFourier = Ank.*besseli(n,kCyl.*rho);
    else
        EzFourier = Bnk.*besselk(n,kVac.*rho);
    end
   
end