function [MDRenv] = FSRenvmdr(n,p,varargin)
%FSRenvmdr computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search
%
%<a href="matlab: docsearchFS('FSRenvmdr')">Link to the help function</a>
%
%  Required input arguments:
%
%    n : scalar, number of obseravtions
%    p : scalar, number of explanatory variables (including the intercept if present)
%
%  Optional input arguments:
%
%   init:       Search initialization. Scalar.
%               Scalar which specifies the initial subset size to monitor
%               minimum deletion residual, if init is not specified it will
%               be set equal to
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%               Example - 'init',100 starts monitoring from step m=100 
%               Data Types - double
%  prob:    quantiles for which envelopes have
%               to be computed. Vector.
%               1 x k vector containing quantiles for which envelopes have
%               to be computed. The default is to produce 1%, 50% and 99%
%               envelopes.
%               Example - 'prob',[0.01 0.99] 
%               Data Types - double
%  exact:    Method for the calculation of the quantiles. Scalar.
%                If it is equal to 1 (default) the calculation of the quantiles
%               of the T and F distribution is based on functions finv and
%               tinv from the Matlab statistics toolbox, otherwise the
%               calculations of the former quantiles is based on functions
%               finvFS and tinvFS. The solution has a tolerance of 1e-8
%               (change variable tol in files finvFS.m and tinvFS.m if
%               required.
%               Example - 'exact',0
%               Data Types - double
%               Remark: the use of functions tinv and finv is more precise
%               but requires more time.
%
%  Output:
%
%  MDRenv:      matrix with n-m0+1 rows and length(prob)+1 columns
%               1st col = fwd search index from m0 to n-1
%               2nd col = envelope for quantile prob(1)
%               3rd col = envelope for quantile prob(2)
%               ...
%               (k+1) col = envelope for quantile prob(k)
%
% Subfunctions: tinvFS, finvFS, tcdfFS, fpdfFS, fcdfFS
%
% Other function dependencies: none.
%
% See also LXS.m, FSREDA.m
%
% References:
%
%   Atkinson, A.C. and Riani, M. (2006). Distribution theory and
%   simulations for tests of outliers in regression. Journal of
%   Computational and Graphical Statistics, Vol. 15, pp. 460–476
%   Riani, M. and Atkinson, A.C. (2007). Fast calibrations of the forward
%   search for testing multiple outliers in regression, Advances in Data
%   Analysis and Classification, Vol. 1, pp. 123–141.
%
% Copyright 2008-2015
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRenvmdr')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{  
    % FSRenvmdr with all default options.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent 
    % envelopes based on 1000observations and 5 explanatory variables using 
    % exact method.
    MDRenv=FSRenvmdr(10000,5);
%}

%{  
    % FSRenvmdr with optional argument.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent 
    % envelopes based on 1000 observations and 5 explanatory variables 
    % using approximate method.
    MDRenv1=FSRenvmdr(10000,5,'exact',0);
%}

%{
    % Example with plot of the envelopes.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent
    % envelopes based on 100observations and 5 explanatory variables using 
    % exact method.
    Menv=FSRenvmdr(100,5,'exact',1);
    plot(Menv(:,1),Menv(:,2:4));
%}

%{
    %% Checking the accurary of the envelopes.
    n=100;
    p=5;
    state=1000;
    randn('state', 1000);
    X=randn(n,5);

    nsimul=1000;
    for j=1:nsimul
        y=randn(n,1);
        [out]=LXS(y,X);
        mdr = FSRmdr(y,X,out.bs,'init',10);
        mdrStore(:,j)=mdr(:,2);
    end
%}



%% Input parameters checks

if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSRenvmdr:missingInputs','n must be scalar non empty and non missing!!');
end

if ~isscalar(p) || isempty(n) || isnan(p)
    error('FSDA:FSRenvmdr:missingInputs','p must be scalar non empty and non missing!!!');
end

if n<40
    inisearch=p+1;
else
    inisearch=min(3*p+1,floor(0.5*(n+p+1)));
end

% Notice that prob must be a row vector
prob=[0.01 0.5 0.99];
options=struct('init',inisearch,'prob',prob,'exact',1);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRenvmdr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

if nargin>2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
m0=options.init;
prob=options.prob;
exact=options.exact;

% Check that the initial subset size is not greater than n-1
if m0>n-1
    error('FSDA:FSRenvmdr:TooLargen',['Initial starting point of the search (m0=' num2str(m0) ') is greater than n-1(n-1=' num2str(n-1) ')']);
end

%% Envelopes generation

% Make sure that prob is a row vector.
if size(prob,1)>1;
    prob=prob';
end

lp = length(prob);
prob = 1-prob;


% m = column vector which contains fwd search index.
m =(m0:n-1)';
lm=length(m);

% mm = fwd search index replicated lp times.
mm = repmat(m,1,lp);

if exact;
    % finv finds the inverse of the F distribution.
    quant=finv(repmat(prob,lm,1),2*(n-mm),2*(mm+1));
else
    % finvFS finds the inverse of the F distribution.
    quant=finvFS(repmat(prob,lm,1),2*(n-mm),2*(mm+1));
end

% from the equivalence with Incomplete beta distribution.
q=(mm+1)./(mm+1+(n-mm).*(quant));

% Minsca = matrix of the scaled MDR envelopes in each step of the search.
if exact;
    MinSca= abs(tinv(0.5*(1+q), mm-p));
else
    MinSca= abs(tinvFS(0.5*(1+q), mm-p));
end


% Compute variance of the truncated normal distribution.
% mm/n is the percentage of observations inside subset.
a=norminv(0.5*(1+mm/n));
%corr=1-(2*a.*normpdf(a))./(2*normcdf(a)-1);
corr=1-2*(n./mm).*a.*normpdf(a);

MDRenv=[m MinSca./sqrt(corr)];

end


%% Subfunction: tinvFS
%
function [x0] = tinvFS(p,v)
%  tinvFS - Student t Inverse Cumulative Distribution function.
%            Used to compute the quantiles of t distribution.
%
% Remark: It is equivalent to MATLAB function tinv
%
%  Usage:  x0 = tinvFS(p,v)
%
%  Input:  p - matrix of percentage points ( 0 < p < 1 )
%          v - matrix of degrees of freedom for t distribution
%
%  Output: x0 - matrix of critical T values st Pr(X < x0) = p and x ~ T(v)
%
%  REMARK:  Uses initial Normal approximation and then iterates for accuracy
%           The default tolerance if 1e-8
%

% Use normal approx to start
t = sqrt(-2.*log(abs((p>0.5) - p))) ;
z = 2.515517 + t.*(0.802853 + t.*0.010328) ;
d = 1 + t.*(1.432788 + t.*(0.189269 + t.*0.001308)) ;
z = t - (z./d) ;
d = sqrt(0.25.*(v<=2) + (1 - 2./v).*(v>2)) ;
x0 = ((p>0.5).*z - (p<=0.5).*z)./d ;

tol = 1e-8; % Change the tolerance if necessary
p = 1 - p;

converge=0;
k=0;
while (converge==0 && k<=500)
    f0 = 1-tcdfFS(x0,v);
    df0 = tpdf(x0,v);
    
    x1 = x0 - (p-f0)./df0 ;
    converge = max(max(abs(x0 - x1))) < tol ;
    x0 = x1;
    k = k + 1;
end;
if not(converge);
    disp('Warning: tinvFS has not converged');
end;

end

%% Subfunction: finvFS
%
function [x0] = finvFS(p,v1,v2)
%  finvFS - Inverse of the F Cumulative Distribution function
%            with v1,v2 degrees of freedom
%
% Remark: It is equivalent to MATLAB function finv
%
%  Usage: x0 = finvFS(p,v1,v2)
%
%  Input:  p  - matrix of percentage points ( 0 < p < 1 )
%          v1 - matrix of numerator df (conformable with p)
%          v2 - matrix of denominator df (conformable with p)
%
%  Output: X  - matrix of critical values st Pr(X < x0) = p and X ~ F(v1,v2)
%
% REMARK: the default tolerance is 1e-8

tol = 1e-8 ;
tol2 = tol^2 ;
% Paulson normal approximation as starting value
t  = sqrt(-2*log(abs((p>0.5) - p)));
z  = 2.515517 + t.*(0.802853 + t.*0.010328);
d  = 1 + t.*(1.432788 + t.*(0.189269 + t.*0.001308));
z  = t - (z./d);
z  = (p>0.5).*z - (p<=0.5).*z;
c  = 2./(9.*v2);
d  = 2./(9.*v1);
a  = 1 - c;
b  = 1 - d;
d  = (a.^2).*d + (b.^2).*c - c.*d.*(z.^2);
c  = a.^2 - c.*(z.^2);
x0 = abs((a.*b + z.*sqrt(d + (tol-d).*(d<0)))./(c + (1-c).*(c<0.3))).^3;
x0 = x0 + (d<=0 |  c< 0.3).*(0.5.*(p<0.5) + 2.0.*(p>=0.5) - x0);
p  = 1 - p ;
converge=0;
k=0;
while (converge==0 && k<=50)
    f0  = 1-fcdfFS(x0,v1,v2) ;
    df0 = fpdfFS(x0,v1,v2);
    
    % Routines from MATLAB programmers
    %  f0  = 1-fcdf(x0,v1,v2) ;
    % df0 = fpdf(x0,v1,v2);
    
    x1 = x0 - (p-f0)./df0 ;
    negative = sum(sum(not((x1 > tol2)))) ;
    if negative>0;
        x1 = x1 + (x1<=tol2).*(x0.*(0.5 + 1.5.*(p< f0)) - x1) ;
    end;
    converge = max(max(abs(x0 - x1))) < tol & not(negative);
    x0 = x1;
    k = k + 1;
end;
if not(converge);
    disp('Warning: finvFS has not converged, exaxt routine finv is used');
    x0=finv(1-p,v1,v2);
end;

end
%% Subfunction: tcdfFS
%
function F = tcdfFS(x,v)
% tcdfFS computes the cdf of the student T distribution
% It is equivalent to MATLAB function tcdf
%
%         F = tcdfFS(x,v)
%

if any(any(v<=0))
    error('FSDA:FSRenvmdr:WrongDf','Degrees of freedom must be positive')
end

v = min(v,1000000); % make it converge and also accept Inf.

neg = x<0;
F = fcdfFS(x.^2,ones(size(v)),v);
F = 1-(1-F)./2;
F = F + (1-2*F).*neg;

end

%% Subfunction: fpdfFS
%
function f = fpdfFS(x,v1,v2)
%   fpdfFS computes the density function for the F distribution
% It is equivalent to MATLAB function fpdf
%         f = fpdfFS(x,df1,df2)
%

c = v2./v1;
xx = x./(x+c);
f = betapdfFS(xx,v1/2,v2/2);
f = f./(x+c).^2.*c;

%% Nested function: betapdfFS
    function d = betapdfFS(x,v1,v2)
        % betapdfFS computes the density function of the beta distribution
        % It is equivalent to MATLAB function betapdf
        % f = betapdfFS(x,v1,v2)
        %
        
        if any(any((v1<=0)|(v2<=0)))
            error('FSDA:FSRenvmdr:WrongV1OrV2','Parameter v1 or v2 is nonpositive')
        end
        
        I = find((x<0)|(x>1));
        
        d = x.^(v1-1) .* (1-x).^(v2-1) ./ beta(v1,v2);
        d(I) = 0*I;
    end
end

%% Subfunction: fcdfFS
%
function F = fcdfFS(x,v1,v2)
% fcdfFS computes the cdf of the F distribution
% It is equivalent to MATLAB function fcdf
%
%         F = fcdfFS(x,v1,v2)
%

x = x./(x+v2./v1);
F = pbeta(x,v1./2,v2./2);

%% Nested function: pbeta
    function F = pbeta(x,v1,v2)
        % pbeta  computes the cdf of the beta distribution
        % It is equivalent to MATLAB function betacdf
        %         F = pbeta(x,v1,v2)
        %
        
        if any(any((v1<=0)|(v2<=0)))
            error('FSDA:FSRenvmdr:WrongV1OrV2','Parameter v1 or v2 is nonpositive')
        end
        
        % Il = find(x<=0);
        % Iu = find(x>=1);
        Ii = find(x>0 & x<1);
        
        F = 0*(x+v1+v2); % Stupid allocation trick
        % F(Il) = 0*Il; % zeros for all values x<=0
        F(x>=1) = 1;  % 0*Iu + 1; % one for all values x>=1
        % Computation of the cdf for all values 0<x<1
        if ~isempty(Ii)
            F(Ii) = betainc(x(Ii),v1(Ii),v2(Ii));
        end
        
    end

end