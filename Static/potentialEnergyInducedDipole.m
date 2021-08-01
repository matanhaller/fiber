% Potential energy due to a cylindrical dipole with ~1/(eta-1) dipole
% moment. Proportional to (eta-1)^-3 and is the asymptotic limit of the
% exact potential energy due to the fiber as eta-->inf.

function Winf = potentialEnergyInducedDipole(epsilonR, eta)
    Winf = (epsilonR - 1) / (epsilonR + 1) * (eta - 1).^(-3);
end