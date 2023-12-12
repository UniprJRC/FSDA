function psix=OPTpsix(u,c)
%OPTpsix computes psi function (derivative of rho function) times x
%
%<a href="matlab: docsearchFS('OPTpsix')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. 
%               vector of  length n containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        tuning parameters. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...) 
%
%
%  Optional input arguments:
%
%  Output:
%
%
%   psix  :    vector of length n which contains the values of the derivative of the optimal psi
%                function multiplied by u, associated to the residuals or Mahalanobis
%                distances for the n units of the sample.
%
% More About:
%
%
% Function OPTpsix transforms vector u as follows
%
% \[
% \psi(x)x =x\rho' (x) = \begin{cases}  \frac{2.7692 x^2}{c^2}   \qquad |x| \leq \frac{2}{3} c \\
%  -\frac{5.3834 x^2}{c^2} +\frac{43.0672 x^4}{c^4}   -\frac{69.9840 x^6}{c^6}  +\frac{32.3 x^8}{c^8} 
%   \qquad \frac{2}{3} c < |x|  \leq c   \\
%   0 & \; \vert x \vert > c. \end{cases}
% \]
% 
%
% Remark: function $\psi(x) *x$ has been precalculated because it is used to find
% the weights in case of Tau estimators
%
% See also: HYPpsix, HApsix, TBpsix
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('OPTpsix')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Plot of psi(x) function (derivative of rho function) times x.
    x=-6:0.01:6;
    psixOPT=OPTpsix(x,1.2);
    plot(x,psixOPT)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi (x) \times x$','Interpreter','Latex')

%}

%% Beginning of code

% Computes Standardized optimal psi function times x(first derivative of rho function)
% \rho'(x) *x

% psix = zeros(size(u));
% absx=abs(u);
% 
% % r^2 /(3.25c^2) if r <=2*c
% inds1 = (absx <= 2*c);
% psix(inds1) = u(inds1).^2 / (3.25*c^2);
% 
% % 1/(3.25) * ( -1.944* (x/c)^2 .... +8*0.002 (x/c)^8 )    if    2c< |x| <3c
% inds1 = (absx > 2*c)&(absx <= 3*c);
% x1 = u(inds1);
% psix(inds1) = (-1.944 * x1.^2 / c^2 + 1.728 * x1.^4 / c^4 - 0.312 * x1.^6 / c^6 + 0.016 * x1.^8 / c^8) / 3.25;

% 0 if |x| >3*c



% psiOPT = zeros(length(u),1);
psix = zeros(size(u));

absx=abs(u);

%  if r <=(2/3)*c
inds1 = absx <= (2/3)*c;
psix(inds1) = 2.7692*u(inds1).^2 / c^2;

%     if    (2/3)*c< |x| <c
inds1 = (absx > (2/3)*c)&(absx <= c);
x1 = u(inds1);
psix(inds1) = -5.3834 * x1.^2 / c^2 + 43.0672 * x1.^4 / c^4 - 69.9840 * x1.^6 / c^6 + 32.3 * x1.^8 / c^8 ;

% 0 if r >c

%% Old implementation in terms of 0<|u|<2c, 2c<|u|<3c, |u|>3c
% \[
%   OPTpsix(u,c) = \left\{
% \begin{array}{cc}                                                                              
%                 \frac{1}{3.25*c^2} u^2                           &                         |u| \leq 2c   \\         
%                
%    (1/3.25) * (-1.944 * \frac{u^2}{c^2} + 1.728 \frac{u^4}{c^4} - 0.312 \frac{x^6}{c^6} + 0.016 \frac{x^8}{c^8}  &   \qquad 2c \leq |u| \leq 3c \\
%                
%                  0                                                                         &                      |u|>3c \\    
% \end{array}
%    \right.
%  \]
% 
% psix = zeros(size(u));
% absx=abs(u);
% 
% % r^2 /(3.25c^2) if r <=2*c
% inds1 = (absx <= 2*c);
% psix(inds1) = u(inds1).^2 / (3.25*c^2);
% 
% % 1/(3.25) * ( -1.944* (x/c)^2 .... +8*0.002 (x/c)^8 )    if    2c< |x| <3c
% inds1 = (absx > 2*c)&(absx <= 3*c);
% x1 = u(inds1);
% psix(inds1) = (-1.944 * x1.^2 / c^2 + 1.728 * x1.^4 / c^4 - 0.312 * x1.^6 / c^6 + 0.016 * x1.^8 / c^8) / 3.25;
% 0 if |x| >3*c
end
%FScategory:UTISTAT