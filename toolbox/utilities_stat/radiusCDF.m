function F_R = radiusCDF(x,v,nu)
% Distribution function of the radius for the Mahalanobis Squared Distance
%
%<a href="matlab: docsearchFS('radiusCDF')">Link to the help function</a>
%
% This is $F_{R}$. 
%

%  Rbeta(x^2/(nu-2+x^2), p/2, nu/2)

if nargin < 3 || isempty(nu) || nu <= 0
    F_R = chi2cdf(x^2,v);
else
    F_R = betacdf(x^2/(nu-2+x^2), v/2, nu/2);
end
end

%FScategory:UTISTAT
