function d = mahalFS(Y,MU,SIGMA)
%mahalFS computes Mahalanobis distances (in squared units) for each row of matrix Y
%
%   d = mahalFS(Y,MU,SIGMA) returns the Mahalanobis distance (in squared
%   units) of each observation in Y using centroid MU and covariance SIGMA.

%% Beginning of code

Ytilde = Y - MU;

% Cholesky factorization exploits SPD structure of covariance matrices.
% R'*R = SIGMA, then inv(SIGMA) = inv(R)*inv(R'), so
% Ytilde*inv(SIGMA).*Ytilde summed = sum((Ytilde/R).^2, 2)
[R, p] = chol(SIGMA);
if p == 0
    Z = Ytilde / R;
    d = sum(Z .* Z, 2);
else
    d = sum((Ytilde / SIGMA) .* Ytilde, 2);
end

end
%FScategory:UTISTAT
