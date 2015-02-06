function w = OPTwei(x,c)
%OPTwei computes weight function psi(u)/u for optimal weight function
%
%<a href="matlab: docsearchFS('optwei')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    x:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...)
%
%  Output:
%
%    w :         n x 1 vector contains the optimal weights associated to the residuals or
%                Mahalanobis distances for the n units of the sample
%
%
% Function OPTwei transforms vector u as follows
%
%               |  1/(3.25*c^2)                                                          |x|<=2
%               |   
%  \psix(x,c) = |  (1/3.25) * (-1.944 * 1 / c^2 + 1.728 * x.^2 / c^4 - 0.312 * x.^2 / c^6 + 0.016 * x.^6 / c^8)     2c<=|x|<=3c
%               |
%               |   0 
%
%
% References:
%
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('optwei')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

x=-6:0.1:6;
weiOPT=OPTwei(x,2);
plot(x,weiOPT)
xlabel('x','Interpreter','Latex')
ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')

%}

%% Beginning of code


% weights are = 0 if r >3*c
w = zeros(length(x),1);
absx=abs(x);

% weights are 1 /(3.25c^2) if |x| <=2*c
w(absx <= 2*c) = 1 / (3.25*c^2);

% weights are 1/(3.25) * (-1.944* (1/c^2)+1.728 * (x^2/c^4) .... +8*0.002 * (x^6/c^8) )    if    2c< |x| <3c
inds = (absx > 2*c)&(absx <= 3*c);
x1 = x(inds);
w(inds) = (-1.944  / c^2 + 1.728 * x1.^2 / c^4 - 0.312 * x1.^4 / c^6 + 0.016 * x1.^6 / c^8) / 3.25;

% 0 for |x| >3c
end