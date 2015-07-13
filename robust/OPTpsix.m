function psi=OPTpsix(x,c)
%OPTpsix computes psi function (derivative of rho function) times x
%
%<a href="matlab: docsearchFS('optpsix')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    x:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...)
%
% Function OPTpsix transforms vector x as follows
%
%               |  (1/3.25*c^2) x^2                                                         |x|<=2
%               |   
%  \psix(x,c) = |  (1/3.25) * (-1.944 * x^2 / c^2 + 1.728 * x.^4 / c^4 - 0.312 * x.^6 / c^6 + 0.016 * x.^8 / c^8)     2c<=|x|<=3c
%               |
%               |   0                                                                      |x|>3c                              
%
%
% Remark: function \psi(x) *x has been precalculated by it is used to find
% the weights in case of Tau estimators
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
%<a href="matlab: docsearchFS('OPTpsix')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{
    x=-6:0.01:6;
    psixOPT=OPTpsix(x,1.2);
    plot(x,psixOPT)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x) \times x$','Interpreter','Latex')

%}

%% Beginning of code

% Computes Standardized optimal psi function times x(first derivative of rho function)
% \rho'(x) *x

psi = zeros(length(x),1);
absx=abs(x);

% r^2 /(3.25c^2) if r <=2*c
inds1 = (absx <= 2*c);
psi(inds1) = x(inds1).^2 / (3.25*c^2);

% 1/(3.25) * ( -1.944* (x/c)^2 .... +8*0.002 (x/c)^8 )    if    2c< |x| <3c
inds1 = (absx > 2*c)&(absx <= 3*c);
x1 = x(inds1);
psi(inds1) = (-1.944 * x1.^2 / c^2 + 1.728 * x1.^4 / c^4 - 0.312 * x1.^6 / c^6 + 0.016 * x1.^8 / c^8) / 3.25;

% 0 if |x| >3*c

end
