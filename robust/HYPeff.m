function [c,A,B,d] = HYPeff(eff, ~,k,traceiter)
%HYPeff finds constant c which is associated to the requested efficiency
%for hyperbolic estimator
%
%<a href="matlab: docsearchFS('hypeff')">Link to the help page for this function</a>
%
%  Required input arguments:
%
%    eff:       scalar which contains the required efficiency (of location
%               of scale estimator)
%               Generally eff=0.85, 0.9 or 0.95
%    p :        scalar, number of response variables (UP TO NOW p=1 JUST
%                   REGRESSION, TODO FOR MUULTIVARIATE ANALYSIS)
%
%  Optional input arguments:
%
%   k        : supremum of the change of variance curve
%              supCVC(psi,x) x \in R
%              Default value is k=4.5
%  traceiter : scalar. If traceiter = 1 it is possible to monitor
%              how the value of the objective function B^2/A
%              gets closer to the target (eff) during the iterations
%
% Output:
%
%  c,A,B,d = scalars associated to the nominal requested efficiency to be
%  inserted in the hyperbolic tangent estimator
%  For example, inside function psi (derivative of rho)  we have
%
% HYPpsi(u) = 	{ u,			                               |u| <= d,
%               {
%		        { \sqrt(A * (k - 1)) * tanh(sqrt((k - 1) * B^2/A)*(c -|u|)/2) .* sign(u)
%		        { 	                 d <= |u| <  c,
%               {
%		        { 0,			                         |u| >= c.
%
%	It is necessary to have 0 < A < B < 2 *normcdf(c)-1- 2*c*normpdf(c) <1
%
%
%
%
% References:
%
%
% Frank R. Hampel, Peter J. Rousseeuw and Elvezio Ronchetti (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% Journal of the American Statistical Association , Vol. 76, No. 375,
% pp. 643-648 (HRR)
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('hypeff')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
%
% Examples:

%{

    % Find value of c, A, B, for a nominal efficiency of 0.8427
    % when k=4.5
 
    ktuning=4.5;
    [c,A,B,d]=HYPeff(0.8427,1,ktuning);
    % In this case
    % c = 3.000130564905703
    % A = 0.604298601602487
    % B = 0.713612241773758
    % d= 1.304379168746527
    % See also Table 2 of HRR p. 645
%}

%% Beginning of code


if (nargin >2),
    if (k < 0) ,
        error('FSDA:HYPeff:WrongK',[' Illegal choice of parameters in hyperbolic tangent estimator: ' ...
            num2str(k) ]')
    end
else
    k=4.5;
end


if nargin <4
    traceiter =0;
end

% c=4 starting value of the iteration
c=4;
step=2;
empeff=0;

eps=1e-8;
iter=0;
while abs(empeff-eff)>eps
    
    % Find parameters A, B and d using routine HYPck
    % disp(c)
    [A,B,d]=HYPck(c,k);
    
    iter=iter+1;
    
    if iter==100
        disp(['Effective tolerance in routine HYPeff=' num2str(abs(empeff-eff))])
        break
    end
    
    if traceiter==1
        disp([iter c empeff eff])
    end
    
    step=step/2;
    empeff=B^2/A;
    if empeff<eff
        c=c+step;
    elseif empeff>eff
        c=max(c-step,0.1);
    else
    end
end

end
