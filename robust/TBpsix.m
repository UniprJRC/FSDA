function psix=TBpsix(u,c)
%TBpsix computes psi function (derivative of rho function) times x for Tukey's biweight  
%
%<a href="matlab: docsearchFS('tbpsix')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator 
%               (beta in regression or mu in the location case ...) 
%
% Function psibi transforms vector x as follows 
% TBpsix(x)=
% x^2[1-(x/c)^2]^2    if |x/c|<=1 
% 0                 if |x/c|>1
% See equation (2.38) p. 29 of Maronna et al. (2006)
%
% Remark: Tukey's biweight  psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('tbpsix')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%{

    x=-6:0.01:6;
    psixTB=TBpsix(x,2);
    plot(x,psixTB)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x)$','Interpreter','Latex')

%}

%% Beginning of code

psix = (abs(u) < c) .* (u.^2) .* ( 1 - (u./c).^2 ).^2 ;
end
