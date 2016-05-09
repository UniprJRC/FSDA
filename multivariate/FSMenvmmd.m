function [MMDenv] = FSMenvmmd(n,v,varargin)
%FSMenvmmd computes the theoretical envelopes of Minimum MD outside subset during the search
%
%<a href="matlab: docsearchFS('FSMenvmmd')">Link to the help function</a>
%
% Required input arguments:
%
% n :           Number of observations. Scalar. Number of observations.
%               Data Types - single | double
% v :           Number of variables. Scalar. Number of variables.
%               Data Types - single | double
%
% Optional input arguments:
%
% init :       Point where to start monitoring required diagnostics. Scalar. 
%              Note that if bsb is supplied, init>=length(bsb). If init is not
%              specified it will be set equal to floor(n*0.6).
%                 Example - 'init',50 
%                 Data Types - double
% prob:        quantiles for which envelopes have
%               to be computed. Vector. Vector containing 1 x k elements .
%               The default is to produce 1 per cent, 50 per cent and 99 per cent envelopes.
%                 Example - 'prob',[0.05 0.95] 
%                 Data Types - double
% exact:      It indicates how to calculate the quantiles of F
%               distribution. Scalar. If it is equal to 1 (default)  is based on function
%               finv and from the Matlab statistics toolbox, otherwise the
%               calculations of the former quantiles is based on functions
%               invcdff. The solution has a tolerance of 1e-8 (change
%               variable tol in files invcdff.m)
%               Remark. the use of function finv is more precise
%               but requires more time.
%                 Example - 'exact',0 
%                 Data Types - double
%   scaled:  It indicates how to compute the envelopes. Scalar. 
%               If scaled=1 the envelopes are produced for
%               scaled Mahalanobis distances (no consistency factor is
%               applied) else the traditional consistency factor is applies
%               (this is the default)
%                 Example - 'scaled',0 
%                 Data Types - double
%
% Subfunctions.
%   invcdff.
%
% Other function dependencies:
%   none.
%
% Output:
%
%  MMDenv=      n-m0+1 x length(prob)+1 columns containing the envelopes
%               for the requested quantiles.
%               1st col = fwd search index from m0 to n-1; 
%               2nd col = envelope for quantile prob[1]; 
%               3rd col = envelope for quantile prob[2]; 
%               ...; 
%               (k+1) col = envelope for quantile prob[k].
%
% See also FSMenvmmd.m, FSM.m
%
% References:
%
%       Riani, M., Atkinson A.C., Cerioli A. (2009). Finding an unknown
%       number of multivariate outliers. Journal of the Royal Statistical
%       Society Series B, Vol. 71, pp. 201–221.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRenvmmd')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%
%{
    %% FSMenvmmd with all default options.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent envelopes based on 10000
    % observations and 5 explanatory variables using exact method.
    MMDenv=FSMenvmmd(10000,5);
    plot(MMDenv(:,1),MMDenv(:,2:end))
%}

%{
    %% FSMenvmmd with otpional arguments.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent envelopes based on 10000
    % observations and 5 explanatory variables. The envelopes are produced for
    % scaled Mahalanobis distances (no consistency factor is applied)
    MMDenv=FSMenvmmd(10000,5,'scaled',1);
    plot(MMDenv(:,1),MMDenv(:,2:end))
%}

%{
    %% Order statistics and simulations envelopes .
    % In this example we compare the accuracy of the envelopes computed with 
    % order statistics with those which come from simulations. 

    % Fix a seed 
    state=1000;

    mtstream = RandStream('shr3cong','Seed',state);
    RandStream.setGlobalStream(mtstream);
    defaultStream = RandStream.getGlobalStream();
    reset(defaultStream)

    % If you run this example in a version older than 7.9 replace the previous
    % four lines with 
    % randn('state', 1000);
    n=200;
    p=3;


    init=25;
    nsimul=1000;
    mmdStore=zeros(n-init,nsimul);

    for j=1:nsimul
        Y=randn(n,p);
        [fre]=unibiv(Y);
        %create an initial subset with the 20 observations with the lowest
        %Mahalanobis Distance
        fre=sortrows(fre,4);
        bs=fre(1:25,1);
        mmd = FSMmmd(Y,bs,'init',init);
        mmdStore(:,j)=mmd(:,2);
    end

    % Sort rows of matrix mmdStore
    mmdStore=sort(mmdStore,2);

    % Create figure which compares empirical and theoretical forward envelopes
    % for minimum deletion residual
    figure;
    hold('on');
    quant=[0.01 0.5 0.99];
    sel=round(nsimul*quant);
    % Plot lines of empirical quantiles
    line(mmd(:,1),mmdStore(:,sel),'LineStyle','--','Color','g');
    % Plots lines of theoretical quantiles using order statistics
    mmdT=FSMenvmmd(n,p,'exact',1,'init',init);
    line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
    xlabel('Subset size m');
%}

%% Input parameters checks

if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSMenvmmd:Wrongn','n must be scalar non empty and non missing!!');
end

if ~isscalar(v) || isempty(n) || isnan(v)
    error('FSDA:FSMenvmmd:Wrongv','v must be scalar non empty and non missing!!!');
end

inisearch=floor(n*0.6);

% Note that prob must be a row vector
prob=[0.01 0.5 0.99];
options=struct('init',inisearch,'prob',prob,'exact',1,'scaled',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSMenvmmd:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
scaled=options.scaled;

% Check that the initial subset size is not greater than n-1
if m0>n-1
    error('FSDA:FSMenvmmd:WrongM0',['Initial starting point of the search (m0=' num2str(m0) ') is greater than n-1(n-1=' num2str(n-1) ')']);
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
    % invcdff finds the inverse of the F distribution.
    quant=invcdff(repmat(prob,lm,1),2*(n-mm),2*(mm+1));
end


% from the equivalence with the incomplete beta distribution.
q=(mm+1)./(mm+1+(n-mm).*(quant));


% cor=(n/(n-1))*v*(mm-1)./(mm-v);
cor=v*((mm+1)./mm).*(mm-1)./(mm-v);
%cor=v*((mm.^2-1)./mm)./(mm-v);

% Minsca = matrix of the scaled minMD envelopes in each step of the search.
if exact;
    MinSca= sqrt(cor.*finv(q,v, mm-v));
else
    MinSca= sqrt(cor.*invcdff(q,repmat(v,size(mm)), mm-v));
end


% Compute Tallis correction factor based on the chi^2 distribution
% mm/n is the percentage of observations inside subset if scaled is not equal to 1.
if scaled==1;
    corr=1;
else
    a=chi2inv(mm/n,v);
    corr=(n./mm).*(chi2cdf(a,v+2));
end

MMDenv=[m MinSca./sqrt(corr)];

end


%% Subfunction: invcdff
%
function [x0] = invcdff(p,n1,n2)
%  INVCDFF - Inverse of the F Cumulative Distribution function
%            with n1,n2 degrees of freedom
%
%  Usage: x0 = invcdff(p,n1,n2)
%
%  Input:  p  - matrix of percentage points ( 0 < p < 1 )
%          n1 - matrix of numerator df (conformable with p)
%          n2 - matrix of denominator df (conformable with p)
%
%  Output: X  - matrix of critical values st Pr(x < X) = P and x ~ F(n1,n2)
%
% Notes: the default tolerance is 1e-8

tol = 1e-8 ;
tol2 = tol^2 ;
% Use Paulson normal approx to start
t  = sqrt(-2*log(abs((p>0.5) - p)));
z  = 2.515517 + t.*(0.802853 + t.*0.010328) ;
d  = 1 + t.*(1.432788 + t.*(0.189269 + t.*0.001308)) ;
z  = t - (z./d) ;
z  = (p>0.5).*z - (p<=0.5).*z ;
c  = 2./(9.*n2) ;
d  = 2./(9.*n1) ;
a  = 1 - c ;
b  = 1 - d ;
d  = (a.^2).*d + (b.^2).*c - c.*d.*(z.^2) ;
c  = a.^2 - c.*(z.^2) ;
x0 = abs((a.*b + z.*sqrt(d + (tol-d).*(d<0)))./(c + (1-c).*(c<0.3))).^3 ;
x0 = x0 + (d<=0 |  c< 0.3).*(0.5.*(p<0.5) + 2.0.*(p>=0.5) - x0) ;
p  = 1 - p ;
converge=0;
k=0;
while (converge==0 && k<=50)
    f0  = 1-pf(x0,n1,n2) ;
    df0 = df(x0,n1,n2);
    
    % Routines from MATLAB programmers
    %  f0  = 1-fcdf(x0,n1,n2) ;
    % df0 = fpdf(x0,n1,n2);
    
    x1 = x0 - (p-f0)./df0 ;
    negative = sum(sum(not((x1 > tol2)))) ;
    if negative>0;
        x1 = x1 + (x1<=tol2).*(x0.*(0.5 + 1.5.*(p< f0)) - x1) ;
    end;
    converge = max(max(abs(x0 - x1))) < tol & not(negative);
    x0 = x1;
    k = k + 1;
end
if not(converge);
    disp('Warning: INVCDFF has not converged, exact routine finv is used');
    x0=finv(1-p,n1,n2);
end

end

%% Subfunction: df
%
function f = df(x,a,b)
%   DF The F density function
%         f = df(x,df1,df2)
%
c = b./a;
xx = x./(x+c);
f = dbeta(xx,a/2,b/2);
f = f./(x+c).^2.*c;

%% Nested function: dbeta
    function d = dbeta(x,a,b)
        %DBETA    The beta density function f = dbeta(x,a,b)
        %
        
        if any(any((a<=0)|(b<=0)))
            error('FSDA:FSMenvmmd:WrongAorB','Parameter a or b is nonpositive')
        end
        
        I = find((x<0)|(x>1));
        
        d = x.^(a-1) .* (1-x).^(b-1) ./ beta(a,b);
        d(I) = 0*I;
    end
end

%% Subfunction: pf
%
function F = pf(x,a,b)
%PF       The F distribution function
%
%         F = pf(x,df1,df2)
%

x = x./(x+b./a);
F = pbeta(x,a./2,b./2);

%% Nested function: pbeta
    function F = pbeta(x,a,b)
        %PBETA    The beta distribution function
        %
        %         F = pbeta(x,a,b)
        %
        
        if any(any((a<=0)|(b<=0)))
            error('FSDA:FSMenvmmd:WrongAorB','Parameter a or b is nonpositive')
        end
        
        % Il = find(x<=0);
        % Iu = find(x>=1);
        Ii = find(x>0 & x<1);
        
        F = 0*(x+a+b); % Stupid allocation trick
        % F(Il) = 0*Il; % zeros for all values x<=0
        F(x>=1) = 1;  % 0*Iu + 1; % one for all values x>=1
        % Computation of the cdf for all values 0<x<1
        if ~isempty(Ii)
            F(Ii) = betainc(x(Ii),a(Ii),b(Ii));
        end
        
    end

end
%FScategory:MULT-Multivariate

