function [Wt,pretain,varargout] = wthin(X,varargin)
%WTHIN Estimates thinning weights on the basis of kernel density estimate
%
%<a href="matlab: docsearchFS('wthin')">Link to the help page for this function</a>
% Last modified 06-Feb-2016
%
%   Computes a probability density estimate of the sample in the N-by-D
%   matrix X, at the values in X or, optionally, in XI given as varargin.
%
%  Required input arguments:
%
%   X: vector or 2-column matrix with the uni/bi-variate data sample on which a
%      probability density estimate is computed. Matrix.
%      Data Types - single | double.
%
%  Optional input arguments:
%
%   bandwidth :  bandwidth value. Scalar. The bandwidth used to estimate
%                the density. It can be one of the following: 'scott' ,
%                'normal', 'robust'.
%               Data Types - char
%               Example - 'method','robust'
%
%   retainby  :  retention method. String. The function used to retain the
%                observations. It can be 
%                'comp2one',  i.e. 1 - pdfe/max(pdfe))
%                'inverse',    i.e. (1 ./ pdfe) / max((1 ./ pdfe)))
%               Data Types - char
%               Example - 'method','comp2one'
%
%
%  Output:
%
%   Wt :        vector of Bernoulli weights. Vector. It is 1 for retained
%               units and 0 for thinned units.
%               Data Types - single | double.
%
%   pretain :   vector of retention probabilities. Vector. The retention
%               probabilities are estimated using a gaussian kernel using
%               function ksdensity.
%               Data Types - single | double.
%
%  Optional Output:
%
%   Xt :        vector of retained units. Vector. It is X(Wt,:).
%               Data Types - single | double.
%
%  See also ksdensity, mvksdensity
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
    % Uni-dimensional thinning.
    % Dataset of regression structures formed by a dense group of 1000 units
    % and another of 100 units. Thinning in along the direction of the
    % predicted values.
    x1 = randn(1000,1);
    %x1 = x1.^2;
    x2 = 8 + randn(100,1);
    %x2 = x2.^2;
    x = [x1 ; x2];
    y = 5*x + 0.9*randn(1100,1);
    plot(x,y,'.');

    % thinning over the predicted values
    [Wt,pretain] = wthin(y);

    plot(x(Wt,:),y(Wt,:),'k.',x(~Wt,:),y(~Wt,:),'r.');
    axis manual
    clickableMultiLegend(['Retained: ' num2str(sum(Wt))],['Thinned:   ' num2str(sum(~Wt))]);
%}

%{
    % Bi-dimensional thinning.
    % Dataset as before, but thinning is done on the bi-variate space.
   
    plot(x,y,'.');

    % thinning over the predicted values
    [Wt2,pretain2] = wthin([x,y]);

    plot(x(Wt2,:),y(Wt2,:),'k.',x(~Wt2,:),y(~Wt2,:),'r.');
    axis manual
    clickableMultiLegend(['Retained: ' num2str(sum(Wt2))],['Thinned:   ' num2str(sum(~Wt2))]);
%}

%% options

options     = struct('retainby','inverse','bandwidth','scott');
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
        bandwidth_types = {'scott' , 'normal', 'robust'};
        
        if  isempty(retainby) || ~(ischar(retainby) && max(strcmp(retainby,retainby_types)))
            retainby = 'inverse';
        end
        if  isempty(bandwidth) || ~(ischar(bandwidth) && max(strcmp(bandwidth,bandwidth_types)))
            bandwidth = 'scott';
        end

end

%% bandwidth selection: Remark: ksdensity uses by default Scott's rule
[n,d] = size(X);
switch bandwidth
    case 'scott'
        sig = mad(X,1,1) / 0.6745;           % Robust estimate of sigma
        bw  = sig * (4/((d+2)*n))^(1/(d+4)); % Scott's rule, optimal for normal distribution
    case 'normal'
        bw = 1.06  * std(X) * n^(-1/5);  % normal reference rule
    case 'robust'
        bw = 0.786 * iqr(X) * n^{-1/5};  % robustified normal reference rule
    otherwise
        bw = 1.06  * std(X) * n^(-1/5);  % normal reference rule
end

%% Compute the density along the predicted values

if d > 1
    support = 'unbounded';
else
    %support = 'positive';
    minX = min(X); maxX = max(X); e = (maxX-minX)/10^(6);
    support = [ (min(X)-e) , (max(X)+e) ];
end

[pdfe,xout,u]  = ksdensity(X,X,'Function','pdf','Support',support,'bandwidth',bw);
varargout{2} = u;
varargout{3} = xout;

% substitute the zero sampling probability with a very small value
pdfe(pdfe<=0)=0.00001;

% convert the density values into the vector of retention probabilities;
% sampling probability should be inversely proportional to the density, but
% different functions are possible
switch retainby
    case 'comp2one'
        % complement to 1
        pretain = 1 - pdfe/max(pdfe);
    case 'inverse'
        % inverse
        pretain = (1 ./ pdfe) / max((1 ./ pdfe));
end

% Thinning: Xt is the retained vector; Wt are the indices of the retained
% points in the original data X.
[Xt , Wt] =  rthin(X , pretain);
varargout{1} = Xt;

end

