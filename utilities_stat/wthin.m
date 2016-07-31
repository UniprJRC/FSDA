function [Wt,pretain,varargout] = wthin(X,varargin)
%WTHIN applies thinning to a uni/bi-dimensional dataset
%
%<a href="matlab: docsearchFS('wthin')">Link to the help page for this function</a>
% Last modified 06-Feb-2016
%
%   Computes retention probabilities and bernoulli (0/1) weights on the
%   basis of data density estimate.
%
%  Required input arguments:
%
%   X :          Input data. Vector or 2-column matrix. The structure
%                contains the uni/bi-variate data to be thinned on the
%                basis of a probability density estimate.
%
%
%  Optional input arguments:
%
%   bandwidth :  bandwidth value. Scalar. The bandwidth used to estimate
%                the density. It can be estimated from the data using
%                function bwe.
%                Data Types - scalar
%                Example - bandwidth,0.35
%
%   retainby  :  retention method. String. The function used to retain the
%                observations. It can be:
%                - 'comp2one', i.e. 1 - pdfe/max(pdfe))
%                - 'inverse' (default),  i.e. (1 ./ pdfe) / max((1 ./ pdfe)))
%                Data Types - char
%                Example - 'method','comp2one'
%
%
%  Output:
%
%   Wt :        vector of Bernoulli weights. Vector. Contains 1 for retained
%               units and 0 for thinned units.
%               Data Types - single | double.
%
%   pretain :   vector of retention probabilities. Vector. These are the
%               probabilities that each point in X will be retained,
%               estimated using a gaussian kernel using function ksdensity.
%               Data Types - single | double.
%
%  Optional Output:
%
%   Xt :        vector of retained units. Vector. It is X(Wt,:).
%               Data Types - single | double.
%
%  See also ksdensity, mvksdensity, bwe
%
%
% References:
%
%   A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
%      Techniques for Data Analysis," Oxford University Press.
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('wthin')">Link to the help page for this function</a>
% Last modified 06-Feb-2016
%
%
% Examples:
%
%{
    % univariate thinning.
    % The dataset is bi-dimensional and contain two collinear groups with
    % regression structure. One group is dense, with 1000 units; the second
    % has 100 units. Thinning in done according to the density of the values
    % predicted by the OLS fit.
    x1 = randn(1000,1);
    x2 = 8 + randn(100,1);
    x = [x1 ; x2];
    y = 5*x + 0.9*randn(1100,1);
    b = [ones(1100,1) , x] \ y;
    yhat = [ones(1100,1) , x] * b;
    plot(x,y,'.',x,yhat,'--');

    % thinning over the predicted values
    [Wt,pretain] = wthin(yhat, 'retainby','comp2one');

    plot(x(Wt,:),y(Wt,:),'k.',x(~Wt,:),y(~Wt,:),'r.');
    axis manual
    clickableMultiLegend(['Retained: ' num2str(sum(Wt))],['Thinned:   ' num2str(sum(~Wt))]);

%}

%{
    % Bi-dimensional thinning.
    % Same dataset, but thinning is done on the original bi-variate data.
   
    plot(x,y,'.');

    % thinning over the original bi-variate data
    [Wt2,pretain2] = wthin([x,y]);

    plot(x(Wt2,:),y(Wt2,:),'k.',x(~Wt2,:),y(~Wt2,:),'r.');
    axis manual
    clickableMultiLegend(['Retained: ' num2str(sum(Wt2))],['Thinned:   ' num2str(sum(~Wt2))]);
%}

%{
    % Use of 'retainby' option.
    % Since the thinning on the original bi-variate data with the default
    % retention method ('inverse') removes too many units, let's try with
    % the less conservative 'comp2one' option.
   
    plot(x,y,'.');

    % thinning over the original bi-variate data
    [Wt2,pretain2] = wthin([x,y], 'retainby','comp2one');

    plot(x(Wt2,:),y(Wt2,:),'k.',x(~Wt2,:),y(~Wt2,:),'r.');
    axis manual
    clickableMultiLegend(['Retained: ' num2str(sum(Wt2))],['Thinned:   ' num2str(sum(~Wt2))]);
    title('"comp2one" thinning over the original bi-variate data');
    
%}

%{
    % Optional output Xt.
    % Same dataset, the retained data are also returned using varagout option.
   
    % thinning over the original bi-variate data
    [Wt2,pretain2,RetUnits] = wthin([x,y]);
    RetUnits
%}


%% options

% for reasons of performance options are checked only if necessary
if nargin > 1
    
    options     = struct('retainby','inverse','bandwidth',0);
    UserOptions = varargin(1:2:length(varargin));
    if ~isempty(UserOptions) && (length(varargin) ~= 2*length(UserOptions))
        error('FSDA:kdebiv:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    if nargin>1
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    
    % retention method: it can be
    % 'comp2one' (i.e. 1 - pdfe/max(pdfe))
    % 'inverse'  (i.e. (1 ./ pdfe) / max((1 ./ pdfe)))
    retainby    = options.retainby;
    
    % the bandwidth used to estimate the density
    bandwidth   = options.bandwidth;
    
    if ~isempty(UserOptions)
        retainby_types  = {'inverse' , 'comp2one'};
        if  isempty(retainby) || ~(ischar(retainby) && max(strcmp(retainby,retainby_types)))
            retainby = 'inverse';
        end
        if  ~isscalar(bandwidth)
            bandwidth = 0;
        end
    end
    
else
    bandwidth = 0;
    retainby  = 0;
end

%% Compute the density along the predicted values

[~,d] = size(X);

if d > 1
    support = 'unbounded';
else
    minX = min(X); maxX = max(X); e = (maxX-minX)/10^(6);
    support = [ (min(X)-e) , (max(X)+e) ];
end

% for some reason, ksdensity is faster if bandwidth is not provided. The
% 'if' statement is only for performance reasons.
if bandwidth == 0
    % bandwidth selection: Remark: ksdensity uses by default Scott's rule.
    [pdfe,xout,u]  = ksdensity(X,X,'Support',support);
else
    [pdfe,xout,u]  = ksdensity(X,X,'Support',support,'bandwidth',bandwidth);
end

varargout{2} = u;
varargout{3} = xout;

% substitute the zero sampling probability with a very small value
pdfe(pdfe<=0)=0.00001;

% convert the density values into the vector of retention probabilities;
% sampling probability should be inversely proportional to the density, but
% different functions are possible.
if retainby == 0
    % in this case the user has not provided optional arguments and accept
    % all defaults. For performance reasons, the 'switch' statement is
    % skipped and the default 'inverse' function is applied.
    pretain = (1 ./ pdfe) / max((1 ./ pdfe));
else
    switch retainby
        case 'comp2one'
            % complement to 1
            pretain = 1 - pdfe/max(pdfe);
        case 'inverse'
            % inverse
            pretain = (1 ./ pdfe) / max((1 ./ pdfe));
    end
end

% Thinning: Xt is the retained vector; Wt are the indices of the retained
% points in the original data X.
[Xt , Wt] =  rthin(X , pretain);
varargout{1} = Xt;

end

