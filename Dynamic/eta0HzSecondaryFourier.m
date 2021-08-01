% Longitudinal Fourier component of the secondary magnetic field.

function eta0HzFourier = eta0HzSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, beta, epsilon)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2);
    kVac = sqrt(kz.^2 - omega.^2);
    
    if rho < 1
        eta0HzFourier = eta0Cnk.*besseli(n,kCyl.*rho);
    else
        eta0HzFourier = eta0Dnk.*besselk(n,kVac.*rho);
    end
   
end