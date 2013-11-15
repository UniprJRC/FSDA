function d = mahalFS(Y,MU,SIGMA)
%mahalFS computes Mahalanobis distances (in squared units) for each row of matrix Y 
%
%   d = mahalFS(Y,MU,SIGMA) returns the Mahalanobis distance (in squared units) of
%   each observation (point) in Y using centroid MU and covariance matrix SIGMA
%
%      d(I) = (Y(I,:)-MU) * SIGMA^(-1) * (Y(I,:)-MU)',
%
%
%<a href="matlab: docsearch('mahalFS')">Link to the help function</a>
%
%  Required input arguments:
%
%         Y :   n x v data matrix; n observations
%               and v variables
%               Rows of Y represent observations, and columns represent
%               variables.
%        MU :   v x 1 vector containing centroid to use
%      SIGMA:   v x v matrix containing covariance matrix which must be used
%       
%
% Copyright 2008-2013.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti 
%            and Vytis Kopustinskas (2009-2010)
%
%<a href="matlab: docsearch('mahalFS')">Link to the help function</a>
% Last modified 02-May-2013


% Examples:

%{
    Y=randn(10,2);
    MU=median(Y); 
    SIGMA=[0.3 0.4; 0.4 1];
    % Compute MD using as centroid the medians and shape matrix SIGMA
    d=mahalFS(Y,MU,SIGMA);
%}

Ytilde = bsxfun(@minus,Y, MU);
d=sum((Ytilde/SIGMA).*Ytilde,2);

end
