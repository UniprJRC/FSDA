function y = corrcdf(x, rho, n, varargin)
%corrpdf correlation coefficient probability density function
%
% Cumulative distribution function of the correlation coefficient under the
% null hypothesis that the n bivariate observations come from a bivariate
% normal distribution with correlation parameter rho
%
%<a href="matlab: docsearchFS('corrcdf')">Link to the help function</a>
%
% Required input arguments:
%
%    x:         Value at which the cdf must be evaluated.
%               Scalar, vector or matrix 3D array. If x is not a
%               scalar all the 3 input arguments (r,rho and n) must have
%               the same size or just the numel of one of the 3 input
%               arguments must be greater than 1
%               Data Types - single | double
%    rho :      value of the correlation coefficient in the population.
%               Scalar, vector or matrix or 3D array. If rho is not a
%               scalar all the 3 input arguments (r,rho and n) must have
%               the same size or just the numel of one of the 3 input
%               arguments must be greater than 1
%               Data Types - single | double
%    n :        sample size.
%               Scalar, vector or matrix or 3D array. If n is not a
%               scalar all the 3 input arguments (r,rho and n) must have
%               the same size or just the numel of one of the 3 input
%               arguments must be greater than 1
%               Data Types - single | double
%
%
% Optional input arguments:
%
%   upper:      upper or lower tail. Scalar character. 
%               if nargin>3  normcdf(...,'upper') computes the upper tail probability of the 
%                normal distribution.  
%               Example - 'upper'
%               Data Types - char
%
%  Output:
%
%    y:         CDF value. Scalar, vector or matrix or 3D array of the same size
%               of input arguments x, rho and n. $y=\int_-1^x f_{t}(x |
%               \rho,n) dt$ is the value of the cdf of the distribution of
%               the correlation coefficient evaluated at x.
%
%
% See also: corrpdf
%
%
% https://mathworld.wolfram.com/CorrelationCoefficientBivariateNormalDistribution.html,
%
% Acknowledgements:
%
% This function follows the lines of MATLAB code developed by Xu Cui,
% https://www.alivelearn.net/?p=709 Stanford University and the file
% exchange submission Joshua Carmichael (2022). sample correlation
% distribution function
% (https://www.mathworks.com/matlabcentral/fileexchange/45785-sample-correlation-distribution-function)
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('corrcdf')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%


%{
    %% An example where x, rho and n are all scalars.
    % Find Pr(r<x|rho=0.1|n=12)
    x=0;
    rho=0.1;
    n=12;
    xs=num2str(x);
    rhos=num2str(rho);
    ns=num2str(n);
    prob=corrcdf(x,rho,n);
    disp(['Pr(r<' xs '|rho=' rhos, ', n=' ns ')'])
    disp(prob);
%}

%{
    %% x is not scalar
    x=-1:0.01:1;
    rho=0;
    n=12;
    rhos=num2str(rho);
    ns=num2str(n);
    prob=corrcdf(x,rho,n);
    % disp(['Pr(r<x|rho=' rhos, ', n=' ns ')'])
    plot(x,prob)
    xlabel("x")
    ylabel(['Pr(r<x|rho=' rhos, ', n=' ns ')'])
%}

%{
    %% Check quality of the approximation when rho=0.
    % Note that if rho=0 the test statistic
    % sqrt(n-2)*r./sqrt(1-r.^2) is distributed as a Student T with (n-2)
    % degrees of freedom
    x=-0.9:0.01:0.9;
    rho=0;
    n=12;
    rhos=num2str(rho);
    ns=num2str(n);
    prob=corrcdf(x,rho,n);
    plot(x,prob)
    xlabel("x")
    ylabel(['Pr(r<x|rho=' rhos, ', n=' ns ')'])
    hold('on')
    testt=sqrt(n-2)*x./sqrt(1-x.^2);
    probt=tcdf(testt,n-2);
    plot(x,probt)
    legend({'cdf using function corrcdf' 'cdf based on Student t'})
%}

%{
    %% An example where rho is not scalar.
    x=0.3;
    rho=(0:0.1:0.8)';
    n=12;
    xs=string(x);
    rhos=string(rho);
    ns=string(n);
    Prob=corrcdf(x,rho,n);
    nameRows="Pr(r<"+xs+"|rho="+ rhos+ ", n="+ ns+ ")=";
    nameRowsT=array2table(Prob,"RowNames",nameRows);
    disp(nameRowsT)
%}

%{
    %%  An example where n is not scalar.
    x=0.3;
    rho=0';
    n=(5:5:50)';
    xs=string(x);
    rhos=string(rho);
    ns=string(n);
    Prob=corrcdf(x,rho,n);
    nameRows="Pr(r<"+xs+"|rho="+ rhos+ ", n="+ ns+ ")=";
    nameRowsT=array2table(Prob,"RowNames",nameRows);
    disp(nameRowsT)
%}

%{
    % An example where n, r and rho have the same dimension.
    x=(0.2:0.1:0.6)';
    rho=(0.1:0.1:0.5)';
    n= (10:10:50)';
    xs=string(x);
    rhos=string(rho);
    ns=string(n);
    Prob=corrcdf(x,rho,n);
    nameRows="Pr(r<"+xs+"|rho="+ rhos+ ", n="+ ns+ ")=";
    nameRowsT=array2table(Prob,"RowNames",nameRows);
    disp(nameRowsT)
%}

%% Beginning of code

if nargin>3 && strcmpi(varargin{end},'upper')
    % Compute upper tail
    uflag=true;
elseif nargin>3 && ischar(varargin{end})&& ~strcmpi(varargin{end},'upper')
    error(message('FSDA:corrcdf:UpperTailProblem'));
else
    uflag=false;
end

if numel(x)>1 && numel(n)==1 &&  numel(rho)==1
    n=repmat(n,size(x));
    rho=repmat(rho,size(x));
    y=zeros(size(rho));
elseif numel(x)==1 && numel(n)>1 &&  numel(rho)==1
    x=repmat(x,size(n));
    rho=repmat(rho,size(n));
    y=zeros(size(n));
elseif numel(x)==1 && numel(n)==1 &&  numel(rho)>1
    x=repmat(x,size(rho));
    n=repmat(n,size(rho));
    y=zeros(size(rho));
else
    % In this case rho, x and n all have the same size
    y=zeros(size(rho));
end

for i=1:numel(x)
    y(i)=integral(@(r)corrpdf(r,rho(i),n(i)),-1,x(i),'AbsTol',1e-12);
end

if uflag==true
    y=1-y;
end

end
%FScategory:UTISTAT