% Longitudinal Fourier component of the secondary magnetic field.

function eta0HzFourier = eta0HzSecondaryFourier(eta0Cnk, eta0Dnk, rho, n, kz, omega, epsilon)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);
    kVac = sqrt(kz.^2 - omega.^2); kVac = real(kVac) + 1j*sign(omega).*imag(kVac);
    
    if rho < 1
        eta0HzFourier = eta0Cnk.*besseli(n,kCyl.*rho);
    else
        eta0HzFourier = eta0Dnk.*besselk(n,kVac.*rho);
    end
   
end