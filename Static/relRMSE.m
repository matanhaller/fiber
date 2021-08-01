% Calculating the relative Root-Mean-Square-Error (RMSE) between two functions

function err = relRMSE(F, G)
    N = numel(F);
    rmse = sqrt((1/N) * sum((F - G) .^ 2));
    Fnorm = sqrt((1/N) * sum(F .^ 2));
    err = rmse / Fnorm;
end