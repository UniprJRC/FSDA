function psider=OPTpsider(x,c)
%OPTpsider computes derivative of psi function (second derivative of rho function) for optimal weight function
%
%<a href="matlab: docsearchFS('optpsider')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    x:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...)
%
% Function OPTpsider transforms vector x as follows
%
%               |  (1/3.25*c^2)                                                          |x|<=2c
%               |   
%\psider(x,c) = |  (1/3.25) * (-1.944 / c^2 + 1.728 * 3*x.^2 / c^4 - 0.312 * 5* x.^4 / c^6 + 0.016 * 7* x.^6 / c^8)     2c<=|x|<=3c
%               |
%               |   0                                                                      |x|>3c                              
%
%
% Remark: Optimal psi-function is almost linear around u = 0 in accordance with
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
%<a href="matlab: docsearchFS('optpsider')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

x=-6:0.01:6;
psiOPTder=OPTpsider(x,1.2);
plot(x,psiOPTder)
xlabel('x','Interpreter','Latex')
ylabel('$\psi''(x)$','Interpreter','Latex')
title('Optimal') 
%}

%% Beginning of code

% Computes Standardized optimal psi function (first derivative of rho function)
% \rho'(x)

psider = zeros(length(x),1);
absx=abs(x);

% 1 /(3.25c^2) if |r| <=2*c
inds1 = absx <= 2*c;
psider(inds1) = 1 / (3.25*c^2);

% 1/(3.25) * ( -1.944*2* (r/c) .... +8*0.002 (r/c)^8 )    if    2c< |r| <3c
inds1=(absx > 2*c)&(absx <= 3*c);
x1 = x(inds1);
psider(inds1) = (-1.944/ c^2 + 5.184 * x1.^2 / c^4 - 1.56 * x1.^4 / c^6 + 0.112 * x1.^6 / c^8) / 3.25;

% 0 if r >3*c

end
