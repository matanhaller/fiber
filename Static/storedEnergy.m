% Calculating stored Energy for some z value

function energy = storedEnergy(intensity, R, dr, dtheta)
    S = (dr)^2 * sin(dtheta) / 2;
    vec = S * reshape(1:2:(2 * size(intensity,1)-1), [], 1);
    weight_mat = zeros(size(intensity));
    weight_mat = repmat(vec,1,size(intensity,2));
    energy = sum(sum(intensity.*weight_mat));
end 