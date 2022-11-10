function psider=PDpsider(u,alpha)
%PDpsider computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator
%
%<a href="matlab: docsearchFS('PDpsider')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    u:         scaled residuals or Mahalanobis distances. Vector. n x 1
%               vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    alpha :    tuning parameter. Scalar. Scalar in the interval (0,1] which
%               controls the robustness/efficiency of the estimator
%               (beta in regression or mu in the location case ...). The
%               greater is alpha the greater is the bdp and smaller is the
%               efficiency.
%
%
%  Optional input arguments:
%
%  Output:
%
%
%   psider :     derivative of psi function. Vector. 
%                n x 1 vector which contains the values of the derivative
%                of the psi function of minimum power divergence estimator
%                associated to the residuals or Mahalanobis distances for
%                the n units of the sample.
%
% More About:
%
% Function PDpsider transforms vector u as follows 
% \[
% PDpsider(u)= \alpha (1- \alpha u^2) \exp{-\alpha u^2/2} 
% \] 
%
%
% See also TBpsider, HUpsider, HYPpsider, HApsider, OPTpsider
%
% References:
%
%  Riani, M. Atkinson, A.C., Corbellini A. and Perrotta A. (2020), Robust
%  Regression with Density Power Divergence: Theory, Comparisons and Data
%  Analysis, submitted.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('PDpsider')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Plot the second derivative of minimum power divergence estimator.
    x=-6:0.01:6;
    alpha=1;
    c=2;
    psiPDder=PDpsider(x,c);
    plot(x,psiPDder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')
%}

%% Beginning of code

psider = alpha *(1- alpha *u.^2).* exp(-alpha * u.^2/2); 

end
%FScategory:UTISTAT