% WNChygepdf.m - evaluates a Wallenius Probability Density.
%   See "Calculation Methods for Wallenius' Noncentral Hypergeometric
%   Distribution", A. Fog, 10/26/2004.
%
%  Created by Jim Huntley,  8/03/06

function [pdf] = WNChygepdf(n,nn,m,N,omega)

n2  = nn-n;
m2  = N-m;
d   = omega*(m-n) + m2 - n2;
bc1 = bc(m,n);   
bc2 = bc(m2,n2); 
sum1 = 0;
for jn = 1:n+1
    sum1 = sum1 + (-1)^(jn-1) * bc(n,jn-1) * beta(n2+1,d+omega*(jn-1));
end
pdf = bc1 * bc2 * d * sum1;

return
