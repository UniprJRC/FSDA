function x = corrinv(p, rho, n)
%computes the quantiles of the sampling distribution of the correlation coefficient
%
% Inverse of the cumulative distribution function of the sampling
% distribution of the correlation coefficient under the null hypothesis
% that the n bivariate observations come from a bivariate normal
% distribution with correlation parameter rho
%
%<a href="matlab: docsearchFS('corrinv')">Link to the help function</a>
%
% Required input arguments:
%
%    p:         Probability at which the inverse of the cdf must be evaluated
%               $0 \leq p \leq 1$.
%               Scalar, vector or matrix 3D array of the same size of x and b.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               Data Types - single | double
%    rho :      Value of the correlation coefficient in the population.
%               Scalar, vector or matrix or 3D array. If rho is not a
%               scalar all the 3 input arguments (p,rho and n) must have
%               the same size or just the numel of one of the 3 input
%               arguments must be greater than 1
%               Data Types - single | double
%    n :        Sample size.
%               Scalar, vector or matrix or 3D array. If n is not a
%               scalar all the 3 input arguments (p,rho and n) must have
%               the same size or just the numel of one of the 3 input
%               arguments must be greater than 1
%               Data Types - single | double
%
%
% Optional input arguments:
%
%
%  Output:
%
%    x:         Inverse cdf value. Scalar, vector or matrix or 3D array of the same size
%               of input arguments p, rho and n. $p=\int_{-1}^x f_{r}(r | \rho,n) dr$ is the
%               inverse of the sample correlation coefficient cdf, given
%               $\rho$, the population correlation coefficient and the
%               sample size $n$, for the corresponding
%               probabilities in p.
%
%
% See also: corrpdf, corrcdf
%
%
% References:
%
% Das Gupta, S. (1980). Distribution of the Correlation Coefficient,
% in: Fienberg, S.E., Hinkley, D.V. (eds) R.A. Fisher: An Appreciation,
% Lecture Notes in Statistics, vol 1. Springer, New York, NY.
% https://doi.org/10.1007/978-1-4612-6079-0_3
%
% Acknowledgements:
%
% For additional information see
% https://mathworld.wolfram.com/CorrelationCoefficientBivariateNormalDistribution.html
% This function follows the lines of MATLAB code developed by Xu Cui,
% and the file exchange submission Joshua Carmichael (2022), sample
% correlation distribution function
% https://www.mathworks.com/matlabcentral/fileexchange/45785-sample-correlation-distribution-function/
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('corrinv')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%

%{
    %% An example where p, rho and n are all scalars.
    % Find x030 such that Pr(r<x030|rho=0.1, n=12)=0.3.
    p=0.3;
    rho=0.1;
    n=12;
    rhos=num2str(rho);
    ns=num2str(n);
    x030=corrinv(p,rho,n);
    disp(['Quantile x030 in f(r |rho=' rhos, ', n=' ns ')= '  num2str(x030)])
    disp('In other words, the probability of obtaining values')
    disp(['of the sample correlation coefficient smaller than' num2str(x030)])
    disp(['is equal to 0.3 when rho='  rhos, ' and n=' ns ])
%}

%{
    %% p is not scalar.
    p=(0.1:0.1:0.9)';
    rho=0;
    n=12;
    rhos=num2str(rho);
    ns=num2str(n);
    x=corrinv(p,rho,n);
    nam={['rho=' rhos, ', n=' ns]};
    Xt=array2table(x,"RowNames","x"+p,"VariableNames",nam);
    disp(Xt)
%}

%{
    %% An example where rho is not scalar.
    p=0.3;
    rho=(0:0.1:0.8)';
    n=12;
        x=corrinv(p,rho,n);
    nam={['x' num2str(p) ' when, n=' num2str(n)]};
    Xt=array2table(x,"RowNames","rho="+rho,"VariableNames",nam);
    disp(Xt)
%}

%% Beginning of code
arguments
    p {mustBeNumeric, mustBeReal, mustBeInRange(p, 0, 1)}
    rho {mustBeNumeric, mustBeReal, mustBeInRange(rho, -1, 1,'exclusive')}
    n {mustBeNumeric, mustBeInteger, mustBePositive}
end

if numel(p)>1 && isscalar(n) &&  isscalar(rho)
    n=repmat(n,size(p));
    rho=repmat(rho,size(p));
    x=zeros(size(rho));
elseif isscalar(p) && numel(n)>1 &&  isscalar(rho)
    p=repmat(p,size(n));
    rho=repmat(rho,size(n));
    x=zeros(size(n));
elseif isscalar(p) && isscalar(n) &&  numel(rho)>1
    p=repmat(p,size(rho));
    n=repmat(n,size(rho));
    x=zeros(size(rho));
else
    % In this case rho, p and n all have the same size
    assert(isequal(size(rho),size(p),size(n)),'rho, p and n have different sizes')
    x=zeros(size(rho));
end

% Define function handle to cdf (assumed available)
% corrcdf(r, n) should return cdf evaluated at r
for i = 1:numel(p)
    if p(i)==0
        x(i)=-Inf;
    elseif p(i)==1
        x(i)=Inf;
    else
        target_p = p(i);

        % Define root function: cdf(r) - p = 0
        f = @(r) corrcdf(r, rho(i), n(i)) - target_p;

        % The correlation coefficient is bounded in [-1, 1]
        % but it make sense to compute the quantiles when rho
        % is strictly inside [-1, 1].
        % Use fzero with a good initial guess
        x(i) = fzero(f, [ -0.999999, 0.999999 ]);
    end
end
end


%FScategory:ProbDist