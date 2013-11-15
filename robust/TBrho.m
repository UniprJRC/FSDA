function rhoTB = TBrho(u,c)
%TBrho computes (rho) function for Tukey biweight
%
%<a href="matlab: docsearch('Tbrho')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator 
%               (beta in regression or mu in the location case ...) 
%
% Function TBrho transforms vector u as follows 
% TBrho(u)=
% 1-[1-(u/c)^2]^3    if |u/c|<=1 
% 1                  if |u/c|>1
% See equation (2.37) p. 29 of Maronna et al. (2006)
%
%
%
% References:
%
% ``Robust Statistics, Theory and Methods'' by Maronna, Martin and Yohai;
% Wiley 2006.
%
%
% Copyright 2008-2013.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti
%
%
%<a href="matlab: docsearch('Tbrho')">Link to the help page for this function</a>
% Last modified 02-May-2013
%
% Examples:

%{

x=-6:0.01:6;
rhoTB=TBrho(x,2);
plot(x,rhoTB)
xlabel('x','Interpreter','Latex')
ylabel('$\rho (x)$','Interpreter','Latex')

%}

%% Beginning of code

w = (abs(u)<=c);
rhoTB = (u.^2/(2).*(1-(u.^2/(c^2))+(u.^4/(3*c^4)))).*w +(1-w)*(c^2/6);

end
