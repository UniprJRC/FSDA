function psider=HUpsider(u,c)
%HUpsider computes derivative of psi function (second derivative of rho function) for Tukey's biweight  
%
%<a href="matlab: docsearchFS('HUpsider')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator 
%               (beta in regression or mu in the location case ...) 
%
% More About:
%
% \[
% HUpsider(u)= \left\{
%    \begin{array}{cc}
%  1                            &  \mbox{if  }  |u/c| \leq 1 \\
%  0                            &  |u/c|>1   \\
% \end{array}
%    \right.
%  \]
%
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
%<a href="matlab: docsearchFS('HUpsider')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

    x=-6:0.01:6;
    psiHUder=HUpsider(x,2);
    plot(x,psiHUder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')

%}

%% Beginning of code
psider=zeros(length(u),1);
psider(abs(u) < c)=1;
end
