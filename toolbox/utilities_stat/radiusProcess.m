function C = radiusProcess(d,v,nu)

% ellipsoids of decreasing radius define increasing levels of trimming


if nargin < 3 || nu <= 0
    nu = [];
end
d = d(:);
n = size(d(:),1);

% sort Mahalanobis distances
d = sort(d,'ascend');

W = zeros(n,1);
for i = 1:n
    alpha = (n-i+1) / (n+1);
    % radius
    r    = radiusQuantile(alpha , v , nu); 
    % radius density
    den  = radiusDensity(r , v , nu);
    W(i) = sqrt(n) * den  * (d(n-i+1) -  r);
end

absW = abs(W);
C    = max(absW);
end
