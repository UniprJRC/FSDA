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

%}

n = size(X,1);
%Compute the density along the predicted values
support = [ (min(X)) , (max(X)+0.1) ];
%support = 'positive';
bandwidth_normal = 1.06  * std(X) * n^{-1/5};  % normal reference rule
%bandwidth_robust = 0.786 * iqr(X) * n^{-1/5}; % normal reference rule
[pdfe,xout,u]  = ksdensity(X,X,'Function','pdf','Support',support,'bandwidth',bandwidth_normal);
varargout{2} = u;
varargout{3} = xout;

% substitute the zero sampling probability with a very small value
pdfe(pdfe<=0)=0.00001;

% Vector of retention probabilities: convert the density values into the
% inverse values, in order to apply a sampling probability that is
% inversely proportional to the density
pretain = 1 - pdfe/max(pdfe);
%pretain = (1 ./ pdfe) / max((1 ./ pdfe));

% Thinning: Xt is the retained vector; Wt are the indices of the retained
% points in the original data X.
[Xt , Wt] =  rthin(X , pretain);
varargout{1} = Xt;

end

