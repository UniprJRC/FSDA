function d = mahalCorAna(Y,MU)
%mahalFS computes Mahalanobis distances (in squared units) for each row of matrix Y 
%
%   d = mahalCorAna(Y,MU,SIGMA) returns the Mahalanobis distance (in squared units) of
%   each observation (point) in Y using centroid MU and covariance matrix
%   (diag(MU))^-1
%
%      d(i) = (Y(i,:)-MU) * diag(MU)^(-1) * (Y(i,:)-MU)',
%
%   Generarally Y is the matrix of profile rows and MU is the centroid of
%   profile rows c' (row vector of column masses) of the contingencey table
%
%<a href="matlab: docsearchFS('mahalCorAna')">Link to the help function</a>
%
%  Required input arguments:
%
% Y :           Input data. Matrix. 
%               I x J Profile rows matrix; n observations and v variables. 
%               Rows of
%               Y represent observations, and columns represent variables.
%                Data Types - single|double
%        MU :   Centroid. Vector. 1 x J vector containing centroid and
%               covariance matrix to use
%       
%
%  Optional input arguments:
%
%  Output:
%
%    d :         Mahalanobis distances. Vector.
%                n x 1 vector which contains the squared Mahalanobid distances.
%   \[
%      d(i) = (y_i-\mu)^T \times  diag(\mu)^{-1} \times (y_i-\mu), \qquad 
%      i=1, 2, \ldots, n
%   \]
%
% See also: mahal
%
%
% References:
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mahalCorAna')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

%{
    % Example of computation of MD.
    % Generate a contingency table.
    nrow=30;
    ncol=5;
    % Fix the marginals of the two rows
    nrowt=200*ones(1,nrow);
    % Fix the marginals of the three columns
    ncolt=1200*ones(ncol,1);
    % Generate the contingency table
    out=rcontFS(nrow,ncol,nrowt,ncolt);
    N=out.m144;
    MU=sum(N,1)/sum(nrowt); 
    % Compute MD ;
    n=sum(N,'all');
    P=N/n;
    ProfileRows=P./sum(P,1);
    d=mahalCorAna(ProfileRows,MU);
%}

%{
    % Find total inertial of contingency table
    % Generate a contingency table.
    nrow=30;
    ncol=5;
    % Fix the marginals of the two rows
    nrowt=200*ones(1,nrow);
    % Fix the marginals of the three columns
    ncolt=1200*ones(ncol,1);
    % Generate the contingency table
    out=rcontFS(nrow,ncol,nrowt,ncolt);
    N=out.m144;
    n=sum(N,'all');
    P=N/n;
    r=sum(P,2); % row masses
    c=sum(P,1); % centroid of row masses
    ProfileRows=P./r;
    d2=r.*mahalCorAna(ProfileRows,c);
    % d2 is the total inertia
    disp('Total inertia')
    disp(sum(d2));
%}

%% Beginning of code
% For versions of MATLAB before 2015b it is necessary to use bxsfun
%Ytilde = bsxfun(@minus,Y, MU);
% d=sum((Ytilde./MU).*Ytilde,2);

d=sum(((Y-MU).^2./MU),2);

end
%FScategory:UTISTAT