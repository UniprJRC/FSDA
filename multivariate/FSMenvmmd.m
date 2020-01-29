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
%
% prob:        quantiles for which envelopes have
%               to be computed. Vector. Vector containing 1 x k elements .
%               The default is to produce 1 per cent, 50 per cent and 99 per cent envelopes.
%                 Example - 'prob',[0.05 0.95] 
%                 Data Types - double
%
%   scaled:  It indicates how to compute the envelopes. Scalar. 
%               If scaled>0 the envelopes are produced for
%               scaled Mahalanobis distances (no consistency factor is
%               applied) else the traditional consistency factor is applied
%               (this is the default)
%                 Example - 'scaled',0 
%                 Data Types - double
%
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
% Riani, M., Atkinson, A.C. and Cerioli, A. (2009), Finding an unknown
% number of multivariate outliers, "Journal of the Royal Statistical
% Society Series B", Vol. 71, pp. 201-221.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSMenvmmd')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%
%{
    %% FSMenvmmd with all default options.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent envelopes based on 10000
    % observations and 5 explanatory variables.
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
    mmdT=FSMenvmmd(n,p,'init',init);
    line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
    xlabel('Subset size m');
%}

%% Beginning of code

% Input parameters checks

if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSMenvmmd:Wrongn','n must be scalar non empty and non missing!!');
end

if ~isscalar(v) || isempty(n) || isnan(v)
    error('FSDA:FSMenvmmd:Wrongv','v must be scalar non empty and non missing!!!');
end

inisearch=floor(n*0.6);

% Note that prob must be a row vector
prob=[0.01 0.5 0.99];
options=struct('init',inisearch,'prob',prob,'scaled',0);

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
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

m0=options.init;
prob=options.prob;
scaled=options.scaled;

% Check that the initial subset size is not greater than n-1
if m0>n-1
    error('FSDA:FSMenvmmd:WrongM0',['Initial starting point of the search (m0=' num2str(m0) ') is greater than n-1(n-1=' num2str(n-1) ')']);
end

%% Envelopes generation

% Make sure that prob is a row vector.
if size(prob,1)>1
    prob=prob';
end

lp = length(prob);
prob = 1-prob;


% m = column vector which contains fwd search index.
m =(m0:n-1)';
lm=length(m);

% mm = fwd search index replicated lp times.
mm = repmat(m,1,lp);

    % finv finds the inverse of the F distribution.
    quant=finv(repmat(prob,lm,1),2*(n-mm),2*(mm+1));


% from the equivalence with the incomplete beta distribution.
q=(mm+1)./(mm+1+(n-mm).*(quant));


% cor=(n/(n-1))*v*(mm-1)./(mm-v);
cor=v*((mm+1)./mm).*(mm-1)./(mm-v);
%cor=v*((mm.^2-1)./mm)./(mm-v);

% Minsca = matrix of the scaled minMD envelopes in each step of the search.
    MinSca= sqrt(cor.*finv(q,v, mm-v));


% Compute Tallis correction factor based on the chi^2 distribution
% mm/n is the percentage of observations inside subset if scaled is not equal to 1.
if scaled > 0
    corr=1;
else
    a=chi2inv(mm/n,v);
    corr=(n./mm).*(chi2cdf(a,v+2));
end

MMDenv=[m MinSca./sqrt(corr)];

end

%FScategory:MULT-Multivariate

