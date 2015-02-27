function Bbound = FSMbonfbound(n,p,varargin)
%FSMbonfbound computes Bonferroni bounds for each step of the  search (in mult analysis)
%
%
%<a href="matlab: docsearchFS('fsmbonfbound')">Link to the help function</a>
%
%  Required input arguments:
%
%    n : scalar, number of observations
%    p : number of variables
%
%  Optional input arguments:
%
%       init:   scalar which specifies the initial subset size to compute
%               Bonferroni bound. If init is not specified it
%               will be set equal to floor(0.5*(n+p+1))+1
%  prob:        1 x k vector containing quantiles for which envelopes have
%               to be computed. The default is to produce 1%, 50% and 99%
%               envelopes.
% distrib:      the statistical distribution used to compute the
%               approximated Bonferroni bounds. Distributions implemented
%               are 'chi2' and 'F' (default).
%
%  Output:
%
%  MBenv:       matrix with n-m0+1 rows and length(prob)+1 columns
%               1st col = fwd search index from m0 to n-1
%               2nd col = bound for quantile prob[1]
%               3rd col = bound for quantile prob[2]
%               ...
%               (k+1) col = bound for quantile prob[k]
%
% Subfunctions: 
%
% Other function dependencies: none.
%
% See also FSMenvmdr and FSRbonfbound
%
% References:
%
%
% Copyright 2008-2015
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('fsmbonfbound')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:
%{
% Example using default options
  n=1000;
  p=5;
  init=floor(0.5*(n+p+1))+1; 
  MMDenv = FSMenvmmd(n,p,'exact',1,'init',init);
  Bbound = FSMbonfbound(n,p,'init',init);
  figure;
  plot(MMDenv(:,1),MMDenv(:,2:end),'r',Bbound(:,1),Bbound(:,2:end),'b');
%}
%{
% Example using default option, init=10 and prob=[0.01 0.05 0.99 0.999]
  n=2000;
  p=15;
  init=100;
  prob=[0.95 0.99 0.999];
  MMDenv = FSMenvmmd(n,p,'exact',1,'init',init,'prob',prob);
  Bbound = FSMbonfbound(n,p,'init',init,'prob',prob);
  figure;
  plot(MMDenv(:,1),MMDenv(:,2:end),'r',Bbound(:,1),Bbound(:,2:end),'b');
%}
%{
% Example plotting distrib=chi2 and F, init=100 and prob=[0.999]
  n=2000;
  p=10;
  init=100;
  prob=[0.99];
  MMDenv = FSMenvmmd(n,p,'exact',1,'init',init,'prob',prob);
  distrib='chi2';
  BboundC = FSMbonfbound(n,p,'init',init,'prob',prob,'distrib',distrib);
  distrib='F';
  BboundF = FSMbonfbound(n,p,'init',init,'prob',prob,'distrib',distrib);
  figure;
  plot(MMDenv(:,1),MMDenv(:,2:end),BboundC(:,1),BboundC(:,2:end),BboundF(:,1),BboundF(:,2:end));
  legend('Order statistic envelope','Bonferroni Chi2 bound','Bonferroni F bound','Location','best');
%}
%% Input parameters checks

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
    for i=1:2:length(varargin);
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
if size(prob,1)>1;
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