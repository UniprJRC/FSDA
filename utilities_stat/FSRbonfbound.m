function Bbound = FSRbonfbound(n,p,varargin)
%FSRbonfbound computes Bonferroni bounds for each step of the search (in linear regression)
%
%<a href="matlab: docsearchFS('FSRbonfbound')">Link to the help function</a>
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
%        init :      Search initialization. Scalar.
%                      It specifies the point where to initialize the search
%                       and start monitoring minimum deletion residual. if init is not
%                       specified it will be set equal to :
%                       p+1, if the sample size is smaller than 40;
%                       min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%                       Example - 'init',100 starts monitoring from step m=100
%                       Data Types - double
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
%  Bbound:      Bonferroni forward envelopes of mdr. Matrix.
%               Matrix with n-m0+1 rows and length(prob)+1 columns:
%               1st col = fwd search index from m0 to n-1,
%               2nd col = bound for quantile prob[1],
%               3rd col = bound for quantile prob[2],
%               ...,
%               (k+1) col = bound for quantile prob[k].
%
%
% See also: FSRenvmdr
%
% References:
%
%   Atkinson, A.C. and Riani, M. (2006). Distribution theory and
%   simulations for tests of outliers in regression. Journal of
%   Computational and Graphical Statistics, Vol. 15, pp. 460-476 
%   Riani, M. and Atkinson, A.C. (2007). Fast calibrations of the forward
%   search for testing multiple outliers in regression, Advances in Data
%   Analysis and Classification, Vol. 1, pp. 123-141.
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRbonfbound')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% mdr with Bonferroni envelopes.
    % Example of creation of 1, 50 and 99 per cent envelopes based on 1000
    % observations and 5 explanatory variables using exact method
    MDRenv = FSRenvmdr(1000,5,'init',10);
    Bbound = FSRbonfbound(1000,5,'init',10);
    plot(MDRenv(:,1),MDRenv(:,2:end),Bbound(:,1),Bbound(:,2:end));
%}

%{
    % mdr with personalized Bonferroni envelopes.
    % Example of creation of 1, 50, 99  and 99.9 per cent envelopes based on 1000
    % observations and 5 explanatory variables using exact method
    MDRenv = FSRenvmdr(1000,5,'init',10,'prob',[0.01 0.5 0.99 0.999]);
    Bbound = FSRbonfbound(1000,5,'init',10,'prob',[0.01 0.5 0.99 0.999]);
    plot(MDRenv(:,1),MDRenv(:,2:5),Bbound(:,1),Bbound(:,2:5));
%}
%% Input parameters checks

if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSRbonfbound:Wrongn','n must be scalar non empty and non missing!!');
end

if ~isscalar(p) || isempty(n) || isnan(p)
    error('FSDA:FSRbonfbound:Wrongp','p must be scalar non empty and non missing!!!');
end


% The default starting point to monitor mdr is equal to the integer part of
% floor(0.5*(n+p+1))
inisearch=floor(0.5*(n+p+1))+1;

% Notice that prob must be a row vector
prob=[0.01 0.5 0.99];
options=struct('init',inisearch,'prob',prob);

if nargin>2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end
m0=options.init;
prob=options.prob;

% Check that the specified starting point is not greater than n-1
if m0>n-1
    error('FSDA:FSRbonfbound:InvalidArg',['Initial starting point of the search (m0=' num2str(m0) ') is greater than n-1(n-1=' num2str(n-1) ')']);
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
% Set to 1 all m-p not greater than 0
mmminusp=mm-p;
mm(mmminusp<=0)=1;
mmminusp(mmminusp<=0)=1;

MinBonf = abs(tinv((1-probm)./(mm+1), mmminusp));

Bbound = [m MinBonf];

end

%FScategory:UTISTAT
