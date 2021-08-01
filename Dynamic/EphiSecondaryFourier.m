% Azimuthal Fourier component of the secondary electric field.

function EphiFourier = EphiSecondaryFourier(Ank, Bnk, eta0Cnk, eta0Dnk, rho, n, kz, omega, x0, y0, z0, beta, epsilon)
    kCyl = sqrt(kz.^2 - epsilon * omega.^2);
    kVac = sqrt(kz.^2 - omega.^2);
    
    if rho < 1
        EphiFourier = ((n.*kz./rho).*Ank.*besseli(n,kCyl.*rho) ...
            - 1j*omega.*kCyl.*eta0Cnk.*besselip(n,kCyl.*rho)) ./ (kCyl.^2);
    else
        EphiFourier = ((n.*kz./rho).*Bnk.*besselk(n,kVac.*rho) ...
            - 1j*omega.*kVac.*eta0Dnk.*besselkp(n,kVac.*rho)) ./ (kVac.^2);
    end
   
end