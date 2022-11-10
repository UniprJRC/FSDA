function [MDRenv] = FSRenvmdr(n,p,varargin)
%FSRenvmdr computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search
%
%<a href="matlab: docsearchFS('FSRenvmdr')">Link to the help function</a>
%
%  Required input arguments:
%
%    n : number of observations. Scalar. Number of observations on which
%       the envelopes are based.
%    p : number of explanatory variables (including the intercept if
%    present). Scalar. Number of expl. variables on which
%       the envelopes are based.
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
%
%  prob:    quantiles for which envelopes have
%               to be computed. Vector.
%               1 x k vector containing quantiles for which envelopes have
%               to be computed. The default is to produce 1%, 50% and 99%
%               envelopes.
%               Example - 'prob',[0.01 0.99]
%               Data Types - double
%
%  Output:
%
%  MDRenv:      forward envelopes of mdr. Matrix. Matrix with n-m0+1 rows
%               and length(prob)+1 columns.
%               1st col = fwd search index from m0 to n-1;
%               2nd col = envelope for quantile prob(1);
%               3rd col = envelope for quantile prob(2)
%               ...
%               (k+1) col = envelope for quantile prob(k).
%
%
% See also: LXS.m, FSReda.m
%
% References:
%
% Atkinson, A.C. and Riani, M. (2006), Distribution theory and
% simulations for tests of outliers in regression, "Journal of
% Computational and Graphical Statistics", Vol. 15, pp. 460-476.
% Riani, M. and Atkinson, A.C. (2007), Fast calibrations of the forward
% search for testing multiple outliers in regression, "Advances in Data
% Analysis and Classification", Vol. 1, pp. 123-141.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRenvmdr')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % FSRenvmdr with all default options.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent
    % envelopes based on 1000 observations and 5 explanatory variables
    MDRenv=FSRenvmdr(10000,5);
%}

%{
    % Example with plot of the envelopes.
    % Example of creation of 1%, 50% and 99%
    % envelopes based on 100 observations and 5 explanatory variables
    Menv=FSRenvmdr(100,5);
    plot(Menv(:,1),Menv(:,2:4));
%}

%{
    %%Comparing the accuracy of the envelopes computed with order statistics with the simulated ones.
    % Fix a seed
    state=1000;
    mtstream = RandStream('shr3cong','Seed',state);
    %RandStream.setDefaultStream(mtstream);
    RandStream.setGlobalStream(mtstream);
    %defaultStream = RandStream.getDefaultStream();
    defaultStream = RandStream.getGlobalStream();
    reset(defaultStream)

    % If you run this example in a version older than 7.9 replace the previous four lines with
    % randn('state', 1000);

    n=200;
    p=3;
    X=randn(n,p);

    init=20;
    nsimul=1000;
    mdrStore=zeros(n-init,nsimul);

    for j=1:nsimul
        y=randn(n,1);
        [out]=LXS(y,X,'nsamp',1000','msg',0);
        mdr = FSRmdr(y,X,out.bs,'init',init);
        mdrStore(:,j)=mdr(:,2);
    end

    % Sort rows of matrix mdrStore
    mdrStore=sort(mdrStore,2);

    % Create figure which compares empirical and theoretical forward envelopes
    % for minimum deletion residual
    figure;
    hold('on');
    quant=[0.01 0.5 0.99];
    sel=round(nsimul*quant);
    % Plot lines of empirical quantiles
    line(mdr(:,1),mdrStore(:,sel),'LineStyle','--','Color','g');
    % Plots lines of theoretical quantiles using order statistics
    mdrT=FSRenvmdr(n,p+1,'init',init);
    line(mdrT(:,1),mdrT(:,2:4),'LineStyle','-','Color','r');
    xlabel('Subset size m');
%}

%% Beginning of code

% Input parameters checks

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


if coder.target('MATLAB')
    % Notice that prob must be a row vector
    prob=[0.01 0.5 0.99];
    
    options=struct('init',inisearch,'prob',prob);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSRenvmdr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
end

if nargin>2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end
m0=options.init;
prob=options.prob;
% Check that the initial subset size is not greater than n-1
if m0>n-1
    error('FSDA:FSRenvmdr:TooLargen','Initial starting point of the search (m0=%d) is greater than n-1 (n-1=%d)', m0, n-1);
end

%% Envelopes generation

% Make sure that prob is a row vector.
if size(prob,1)>1
    probf=prob';
else
    probf=prob;
end

lp = length(probf);
probf = 1-probf;


% m = column vector which contains fwd search index.
m =(m0:n-1)';
lm=length(m);

% mm = fwd search index replicated lp times.
mm = repmat(m,1,lp);

% finv finds the inverse of the F distribution.
quant=finv(repmat(probf,lm,1),2*(n-mm),2*(mm+1));


% from the equivalence with Incomplete beta distribution.
q=(mm+1)./(mm+1+(n-mm).*(quant));

% Minsca = matrix of the scaled MDR envelopes in each step of the search.
MinSca= abs(tinv(0.5*(1+q), mm-p));


% Compute variance of the truncated normal distribution.
% mm/n is the percentage of observations inside subset.
a=norminv(0.5*(1+mm/n));
%corr=1-(2*a.*normpdf(a))./(2*normcdf(a)-1);
corr=1-2*(n./mm).*a.*normpdf(a);

MDRenv=[m MinSca./sqrt(corr)];

end
%FScategory:REG-Regression