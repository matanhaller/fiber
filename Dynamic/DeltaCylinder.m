% Determinant of the dielectric cylinder's eigenvalue equation.

function d = DeltaCylinder(n, kz, omega, epsilon)
    kVac = sqrt(kz.^2 - omega.^2); kVac = real(kVac) + 1j*sign(omega).*imag(kVac);    
    kCyl = sqrt(kz.^2 - epsilon * omega.^2); kCyl = real(kCyl) + 1j*sign(omega).*imag(kCyl);

    I = besseli(n,kCyl);
    Ip = besselip(n,kCyl);
    K = besselk(n,kVac);
    Kp = besselkp(n,kVac);
    
    M11 = 1j*omega.*(epsilon.*Ip./kCyl - (Kp./K).*I./kVac);
    M12 = n.*kz.*(1./(kCyl.^2) - 1./(kVac.^2)) .* I;
    M21 = M12;
    M22 = -1j*omega.*(Ip./kCyl - (Kp./K).*I./kVac);
    
    d = M11 .* M22 - M12 .* M21;
end