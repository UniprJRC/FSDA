function d = mahalFS(Y,MU,SIGMA)
%mahalFS computes Mahalanobis distances (in squared units) for each row of matrix Y 
%
%   d = mahalFS(Y,MU,SIGMA) returns the Mahalanobis distance (in squared units) of
%   each observation (point) in Y using centroid MU and covariance matrix SIGMA
%
%      d(i) = (Y(i,:)-MU) * SIGMA^(-1) * (Y(i,:)-MU)',
%
%
%<a href="matlab: docsearchFS('mahalFS')">Link to the help function</a>
%
%  Required input arguments:
%
% Y :           Input data. Matrix. 
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%        MU :   Centroid. Vector.  1 x v vector containing centroid to use
%      SIGMA:   Covariance matrix. Matrix. v x v matrix containing covariance matrix which must be used
%       
%
%  Optional input arguments:
%
%  Output:
%
%    d :         Mahalanobis distances. Vector.
%                n x 1 vector which contains the squared Mahalanobid distances.
%   \[
%      d(i) = (y_i-\mu)^T \times  \Sigma^{-1} \times (y_i-\mu), \qquad 
%      i=1, 2, \ldots, n
%   \]
%
% See also: mahal
%
%
% References:
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mahalFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

%{
    % Example of computation of MD.
    Y=randn(10,2);
    MU=median(Y); 
    SIGMA=[0.3 0.4; 0.4 1];
    % Compute MD using as centroid the medians and shape matrix SIGMA
    d=mahalFS(Y,MU,SIGMA);
%}
%% Beginning of code

Ytilde = bsxfun(@minus,Y, MU);
d=sum((Ytilde/SIGMA).*Ytilde,2);

end
%FScategory:UTISTAT