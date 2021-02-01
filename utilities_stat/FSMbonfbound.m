function Bbound = FSMbonfbound(n,p,varargin)
%FSMbonfbound computes Bonferroni bounds for each step of the  search (in mult analysis)
%
%
%<a href="matlab: docsearchFS('FSMbonfbound')">Link to the help function</a>
%
%  Required input arguments:
%
%    n : number of observations. Scalar. Number of observations on which
%       the envelopes are based.
%               Data Types - single | double
%    p : number of variables. Scalar. Number of variables on which
%       the envelopes are based.
%               Data Types - single | double
%
%  Optional input arguments:
%
% init :       Point where to start monitoring required diagnostics. Scalar. 
%              Note that if bsb is supplied, init>=length(bsb). If init is not
%              specified it will be set equal to floor(0.5*(n+p+1))+1.
%                 Example - 'init',50 
%                 Data Types - double
%
% prob:        quantiles for which envelopes have
%               to be computed. Vector. Vector containing 1 x k elements .
%               The default is to produce 1 per cent, 50 per cent and 99 per cent envelopes.
%                 Example - 'prob',[0.05 0.95] 
%                 Data Types - double
%
% distrib:      Reference distribution to use. Character.
%               The statistical distribution used to compute the
%               approximated Bonferroni bounds. Distributions implemented
%               are 'chi2' and 'F' (default).
%                 Example - 'distrib','chi2'
%                 Data Types - char
%
%  Output:
%
%  Bbound:      Bonferroni forward envelopes of mmd. Matrix.
%               Matrix with n-m0+1 rows and length(prob)+1 columns:
%               1st col = fwd search index from m0 to n-1,
%               2nd col = bound for quantile prob[1],
%               3rd col = bound for quantile prob[2],
%               ...,
%               (k+1) col = bound for quantile prob[k].
%
% See also: FSMenvmmd, FSRbonfbound
%
% References:
%
% Atkinson, A.C. and Riani, M. (2006), Distribution theory and
% simulations for tests of outliers in regression, "Journal of
% Computational and Graphical Statistics", Vol. 15, pp. 460-476.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSMbonfbound')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:
%{
    %% Example using default options.
    n=1000;
    p=5;
    init=floor(0.5*(n+p+1))+1; 
    MMDenv = FSMenvmmd(n,p,'init',init);
    Bbound = FSMbonfbound(n,p,'init',init);
    figure;
    plot(MMDenv(:,1),MMDenv(:,2:end),'r',Bbound(:,1),Bbound(:,2:end),'b');
%}
%{
    % Options init and prob.
    % Example using option, init=10 and prob=[0.01 0.05 0.99 0.999]
    n=2000;
    p=15;
    init=100;
    prob=[0.95 0.99 0.999];
    MMDenv = FSMenvmmd(n,p,'init',init,'prob',prob);
    Bbound = FSMbonfbound(n,p,'init',init,'prob',prob);
    figure;
    plot(MMDenv(:,1),MMDenv(:,2:end),'r',Bbound(:,1),Bbound(:,2:end),'b');
%}

%{
      % Comparison between chi2 and F distributions.
      % Example plotting distrib=chi2 and F, init=100 and prob=[0.999].
      n=2000;
      p=10;
      init=100;
      prob=[0.99];
      MMDenv = FSMenvmmd(n,p,'init',init,'prob',prob);
      distrib='chi2';
      BboundC = FSMbonfbound(n,p,'init',init,'prob',prob,'distrib',distrib);
      distrib='F';
      BboundF = FSMbonfbound(n,p,'init',init,'prob',prob,'distrib',distrib);
      figure;
      plot(MMDenv(:,1),MMDenv(:,2:end),BboundC(:,1),BboundC(:,2:end),BboundF(:,1),BboundF(:,2:end));
      legend('Order statistic envelope','Bonferroni Chi2 bound','Bonferroni F bound','Location','best');
%}

%% Beginning of code

% Input parameters checks
if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSMbonfbound:Wrongn','n must be scalar non empty and non missing!!');
end

if ~isscalar(p) || isempty(n) || isnan(p)
    error('FSDA:FSMbonfbound:Wrongp','p must be scalar non empty and non missing!!!');
end

% The default starting point to monitor mdr is equal to the integer part of
% floor(0.5*(n+p+1))
inisearch=floor(0.5*(n+p+1))+1;

% Notice that prob must be a row vector
prob=[0.01 0.5 0.99];

distrib='F';

options=struct('init',inisearch,'prob',prob,'distrib',distrib);


if nargin>2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end
m0=options.init;
prob=options.prob;
distrib=options.distrib;
% Check that the initial subset size is not greater than n-1
if m0>n-1
    error('FSDA:FSMbonfbound:WrongM0',['Initial starting point of the search (m0=' num2str(m0) ') is greater than n-1(n-1=' num2str(n-1) ')']);
end

%% Bonferroni bound generation

% Make sure that prob is a row vector.
if size(prob,1)>1
    prob=prob';
end

% m = column vector which contains fwd search index.
m =(m0:n-1)';

% mm = fwd search index replicated lp times.
lp = length(prob);
mm = repmat(m,1,lp);
probm = repmat(prob,length(m),1);

if strcmp(distrib,'chi2')
    MinBonf = sqrt((chi2inv(1-((1-probm)./(mm+1)),p)));
else
    MinBonf = sqrt((n/(n-1)*p).*((mm-1)./(mm-p)).*finv(1-((1-probm)./(mm+1)),p,(mm-p)));
    %MinBonf = sqrt(((m-1).^2./m).*betainv(1-((1-probm)./(mm+1)),p/2,(mm-p-1)/2));
end          
Bbound = [m MinBonf];
end

%FScategory:UTISTAT