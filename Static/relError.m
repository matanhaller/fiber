% Relative error between two values.

function err = relError(x, y)
    err = abs(x - y) ./ abs(x);
end