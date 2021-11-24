% Determinant of the dielectric cylinder's eigenvalue equation.

function d = DeltaCylinder(n, kz, omega, epsilon)
    kVac = sqrt(kz.^2 - omega.^2);    
    kCyl = sqrt(kz.^2 - epsilon * omega.^2);

    I = besseli(n,kCyl);
    Ip = besselip(n,kCyl);
    K = besselk(n,kVac);
    Kp = besselkp(n,kVac);
    
    M11 = 1j*omega.*(Kp./kVac - epsilon.*Ip./I.*K./kCyl);
    M12 = n.*kz.*(1./(kVac.^2) - 1./(kCyl.^2)) .* K;
    M21 = M12;
    M22 = -1j*omega.*(Kp./kVac - Ip./I.*K./kCyl);
    
    d = M11 .* M22 - M12 .* M21;
end