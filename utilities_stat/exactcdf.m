function p=exactcdf(x,empdist)
%exactcdf finds exact p-values
%
%<a href="matlab: docsearchFS('exactcdf')">Link to the help page for this function</a>
%
%
% Function for finding the exact cdf of each element in the vector x with
% respect to the empirical distribution, represented by the vector empdist, i.e.
% the generic element i of the output vector p is the result of:
% \[
%  \frac{ \displaystyle \sum_{j=1}^K I_{empdist(j) \leq x_i}}{K}
% \] 
% where $I$ is the indicator function and $K$ is the length of vector
% $empdist$
%
% Required input arguments:
%
%       x : empirical replicates of a test.
%           Vector. Vetor of length k containing the empirical realization
%           of a generic test
%                 Data Types - double
%
%  Optional input arguments:
%
%    empdist: empirical distribution of the same test. Vector.
%             Vector of length $K$ generally with $K \geq k$ containig the
%             empirical simulated distribution of the  test. If this
%             optional argument is not supplied the empirical distribution
%             is taken from input vector x.
%                 Example - randn(K,1)
%                 Data Types - double
%
%  Output:
%
%     p    : empirical cdf. Vector. Vector with the same length of
%           input vector x containing the empirical cdf of each element
%           of input vector x. More precisely:  $p(i)$ is computed as
%            \[
%             \frac{ \displaystyle \sum_{j=1}^K I_{empdist(j) \leq x_i}}{K}
%            \] 
%
% See also: normcdf.m
%
% References:
%
% Athey, S., Eckles, D., & Imbens, G. W. (2018). Exact p-values for network
% interference, "Journal of the American Statistical Association", Vol.
% 113, pp. 230-240.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('exactcdf')">Link to the help function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit

% Examples:

%{
    % exactcdf with just one input argument.
    k=1000;
    x=randn(k,1);
    p=exactcdf(x);
%}

%{
    %% exactcdf with two input arguments.
    k=10;
    x=randn(k,1);
    K=100000;
    empdist=randn(K,1);
    % Compute empirical cdf for each element of vector x.
    p=exactcdf(x,empdist);

    % Compute theoretical cdf based on normcdf
    pTheo=normcdf(x);

    % Compare empirical cdf with theoretical cdf 
    plot(p,pTheo,'o')
    xlabel('Empirical cdf')
    ylabel('Theoretical cdf')
%}

%{
    % Using exactcdf for calculating exact p-values.
    k=10;
    x=randn(k,1);
    K=100000;
    empdist=randn(K,1);
    % Compute empirical cdf for each of element of vector x.
    p=exactcdf(x,empdist);

    % Compute exact p-values for an unilateral right-tailed test
    pval_rt=1-p;

    % Compute exact p-values for an unilateral left-tailed test
    pval_lt=p;

%}

%% Beginning of code

if nargin<2
    empdist=x;
end

transpose=0;
[nx,k]=size(x);
if k>1
    x=x';
    [nx,~]=size(x);
    transpose=1; 
end

[ny,k]=size(empdist);
if k>1
    empdist=empdist';
    [ny,~]=size(empdist);
end
X=[x zeros(nx,1) (1:nx)'];
Y=[empdist ones(ny,1) zeros(ny,1)];
XY=[X;Y];
[~,ord]=sort(XY(:,1));
XY=XY(ord,:);
F=cumsum(XY(:,2))/ny;
pos=XY(:,3);
tmp=[pos(pos>0) F(pos>0)];
[~,ord]=sort(tmp(:,1));
p=tmp(ord,2);
if transpose==1
    p=p';
end
end
%FScategory:UTISTAT
