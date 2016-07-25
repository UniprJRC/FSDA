function [Wt,pretain,varargout] = wthin(X,varargin)
%WTHIN Estimates thinning weights on the basis of kernel density estimate
%
%<a href="matlab: docsearchFS('wthin')">Link to the help page for this function</a>
% Last modified 06-Feb-2016
%
%   Computes a probability density estimate of the sample in the N-by-D
%   matrix X, at the values in X or, optionally, in XI given as varargin.
%   The estimation is based on a product Gaussian kernel function using a
%   bandwidth estimated along ... .
%
%  Required input arguments:
%
%
%
%  Optional input arguments:
%
%
%
%  Output:
%
%
%  Optional Output:
%
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
    % xxxx
    % Compute ...
    x1 = randn(1000,1);
    %x1 = x1.^2;
    x2 = 8 + randn(100,1);
    %x2 = x2.^2;
    x = [x1 ; x2];
    y = 5*x + 0.2*randn(1100,1);

    plot(x,y,'.');

    [Wt,pretain,varargout] = wthin(y);

    plot(x(Wt,:),y(Wt,:),'k.',x(~Wt,:),y(~Wt,:),'r.');
    clickableMultiLegend('Retained','Thinned');

%}

n = size(X,1);

% retention method: 
% retain_method = 'comp2one' (i.e. 1 - pdfe/max(pdfe))
% retain_method = 'inverse'  (i.e. (1 ./ pdfe) / max((1 ./ pdfe)))
retain_method = 'inverse';

% bandwidth selection
bandwidth = 'normal';
switch bandwidth
    case 'normal'
        bw = 1.06  * std(X) * n^(-1/5);  % normal reference rule
    case 'robust'
        bw = 0.786 * iqr(X) * n^{-1/5};  % normal reference rule
    otherwise
        bw = 1.06  * std(X) * n^(-1/5);  % normal reference rule
end

%Compute the density along the predicted values
minX = min(X); maxX = max(X); e = (maxX-minX)/10^(6);
support = [ (min(X)-e) , (max(X)+e) ];

%support = 'positive';
[pdfe,xout,u]  = ksdensity(X,X,'Function','pdf','Support',support,'bandwidth',bw);
varargout{2} = u;
varargout{3} = xout;

% substitute the zero sampling probability with a very small value
pdfe(pdfe<=0)=0.00001;

% convert the density values into the vector of retention probabilities;
% sampling probability should be inversely proportional to the density, but
% different functions are possible
switch retain_method
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

