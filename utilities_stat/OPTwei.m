function w = OPTwei(u,c)
%OPTwei computes weight function psi(u)/u for optimal weight function
%
%<a href="matlab: docsearchFS('OPTwei')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        tuning parameters. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...) 
%
%  Optional input arguments:
%
%  Output:
%
%    w :         n x 1 vector which contains the optimal weights
%                associated to the scaled residuals or Mahalanobis
%                distances for the n units of the sample.
%
%
% More About:
%
% Function OPTwei transforms vector u as follows
%
%  \[
%   OPTwei(u,c) = \left\{
%    \begin{array}{cc}
%                  1/(3.25*c^2)                      &           |u| \leq  2 \\
%                
%    (1/3.25) \left( -1.944 * 1 / c^2 + 1.728  \frac{u^2}{c^4} - 0.312\frac{u^2}{c^6} + 0.016 \frac{u.^6}{c^8} \right) &  \qquad  2c\leq  |u|\leq 3c \\
%                
%                   0
% \end{array}
%    \right.
%  \]
%
%  Remark: Yohai and Zamar (1997)  showed that the optimal $\rho$ function 
%  is optimal in the following highly desirable sense: the final M estimate
%  has a breakdown point of one-half and minimizes the maximum bias under
%  contamination distributions (locally for small fraction of
%  contamination), subject to achieving a desidered nominal asymptotic
%  efficiency when the data are Gaussian.
%
% See also: HYPwei, HAwei, TBwei
%
% References:
%
% YOHAI V.J., ZAMAR R.H. (1997) Optimal locally robust M-estimates of
% regression. J Plan Stat Inference 64, pp. 309-323
%
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('OPTwei')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Plot of weight function.
    x=-6:0.1:6;
    weiOPT=OPTwei(x,2);
    plot(x,weiOPT)
    xlabel('x','Interpreter','Latex')
    ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')

%}

%% Beginning of code


% weights are = 0 if r >3*c
w = zeros(length(u),1);
absx=abs(u);

% weights are 1 /(3.25c^2) if |x| <=2*c
w(absx <= 2*c) = 1 / (3.25*c^2);

% weights are 1/(3.25) * (-1.944* (1/c^2)+1.728 * (x^2/c^4) .... +8*0.002 * (x^6/c^8) )    if    2c< |x| <3c
inds = (absx > 2*c)&(absx <= 3*c);
x1 = u(inds);
w(inds) = (-1.944  / c^2 + 1.728 * x1.^2 / c^4 - 0.312 * x1.^4 / c^6 + 0.016 * x1.^6 / c^8) / 3.25;

% 0 for |x| >3c
end
%FScategory:UTISTAT