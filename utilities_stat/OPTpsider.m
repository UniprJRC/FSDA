function psider=OPTpsider(u,c)
%OPTpsider computes derivative of psi function (second derivative of rho function) for optimal weight function
%
%<a href="matlab: docsearchFS('OPTpsider')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. 
%               vector containing residuals or Mahalanobis distances
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
%   psider :     vector which contains the values of the derivative of the optimal psi
%                function associated to the residuals or Mahalanobis
%                distances for the n units of the sample.
%
% More About:
%
% Function OPTpsider transforms vector x as follows
%
% \[
% \psi'(x) = \begin{cases}  \frac{2.7692}{c^2}   \qquad |x| \leq \frac{2}{3} c \\
%  -\frac{5.3834}{c^2} +\frac{3*43.0672 x^2}{c^4}   -\frac{5*69.9840 x^4}{c^6}  +\frac{7*32.3 x^6}{c^8} 
%   \qquad \frac{2}{3} c < |x|  \leq c   \\
%   0 & \; \vert x \vert > c. \end{cases}
% \]

%
% Remark: Optimal psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
% See also: HYPpsider, HApsider, TBpsider
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('OPTpsider')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Plot of derivative of psi function. 
    x=-6:0.01:6;
    psiOPTder=OPTpsider(x,1.2);
    plot(x,psiOPTder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')
    title('Optimal') 
%}

%% Beginning of code

% Computes derivative of  optimal psi function (psi function = first derivative of rho function)

psider = zeros(size(u));

absx=abs(u);

%  if r <=(2/3)*c
inds1 = absx <= (2/3)*c;
psider(inds1) = 2.7692 / c^2;

%     if    (2/3)*c< |x| <c
inds1 = (absx > (2/3)*c)&(absx <= c);
x1 = u(inds1);
psider(inds1) = -5.3834 / c^2 + 3*43.0672 * x1.^2 / c^4 - 5*69.9840 * x1.^4 / c^6 + 7*32.3 * x1.^6 / c^8 ;

%% Old implementation in terms of 0<|u|<2c, 2c<|u|<3c, |u|>3c
% \[
%   OPTpsider(u,c) = \left\{
% \begin{array}{cc}
% (1/3.25*c^2)           &           |u| \leq 2c   \\                                                     
% (1/3.25) \left( -1.944  \frac{1}{c^2} + 1.728  \frac{3u^2}{c^4} - 0.312 \frac{5u^4}{c^6} + 0.016 \frac{7u^6}{c^8}   \right)  &  \qquad 2c \leq |u| \leq 3c \\
%    0            &                      |u|>3c \\       
% \end{array}
%    \right.
%  \]                                                                 
%
% Computes Standardized optimal psi function (first derivative of rho function)
% \rho'(x)
% 
% psider = zeros(size(u));
% absx=abs(u);
% 
% % 1 /(3.25c^2) if |r| <=2*c
% inds1 = absx <= 2*c;
% psider(inds1) = 1 / (3.25*c^2);
% 
% % 1/(3.25) * ( -1.944*2* (r/c) .... +8*0.002 (r/c)^8 )    if    2c< |r| <3c
% inds1=(absx > 2*c)&(absx <= 3*c);
% x1 = u(inds1);
% psider(inds1) = (-1.944/ c^2 + 5.184 * x1.^2 / c^4 - 1.56 * x1.^4 / c^6 + 0.112 * x1.^6 / c^8) / 3.25;
% 
% % 0 if r >3*c
end
%FScategory:UTISTAT