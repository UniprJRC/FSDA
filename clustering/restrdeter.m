function [out]  = restrdeter(eigenvalues, niini, restr, tol, userepmat)
%restrdeter computes determinant restriction 
%
%<a href="matlab: docsearchFS('restrdeter')">Link to the help function</a>
%
%   restrdeter restricts the determinant according to the constraint
%   specified in scalar restr. This function is called in every
%   concentration step of function tclust in case determinant restriction
%   is needed
%
%  Required input arguments:
%
%eigenvalues: Eigenvalues. Matrix. v x k matrix containing the eigenvalues
%             of the covariance matrices of the k groups.
%             v is the number of variables of the dataset which has to be
%             clustered.
%     niini: Cluster size. Column vector. k x 1 vector containing the size
%             of the k clusters
%     restr: Restriction factor. Scalar. Scalar containing the restr parameter in tclust program.
%            More in detail, parameter restr defines the cluster's shape
%            restrictions, which are applied on all clusters during each
%            iteration.
%            Setting restr to 1, yields the strongest restriction,
%            forcing all eigenvalues/determinants to be equal and so the
%            method looks for similarly scattered (respectively spherical)
%            clusters.
%
%  Optional input arguments:
%
%      tol : tolerance. Scalar defining the tolerance of the procedure.
%            The default value is 1e-8
%               Example - 'tol',[1e-18] 
%               Data Types - double
% userepmat : use builtin repmat. Scalar. If userepmat is true function repmat is used instead
%             of bsxfun inside the procedure. Remark: repmat is built in
%             from MATLAB 2013b so it is faster to use repmat if the
%             current version of MATLAB is >2013a
%               Example - 'userepmat',1 
%               Data Types - double
%
%  Output:
%
%
%            out      : Restricted eigenvalues which satisfy the
%                       determinant constraint. Matrix. v-by-k matrix
%                       containing restricted eigenvalues. 
%                       The ratio between the determinants (that is the
%                       product of the columns of matrix out) is not
%                       greater than restr
%
% See also tclust, restreigen, tclustreg
%
% References:
%
% This function implements the algorithm described in 
% Fritz H. Garcia-Escudero, L.A. and Mayo-Iscar, A. (2012), A fast
% algorithm for robust constrained clustering. Available at
% http://www.eio.uva.es/infor/personas/tclust_algorithm.pdf
%
% Copyright 2008-2017.
% Written by FSDA team
%
% DETAILS. This algorithm solves the minimization problem with constraints
% without resorting to the Dykstra algorithm. This implementation is based
% on the paper  "A fast algorithm for robust constrained clustering" by
% Fritz H., Garcia Escudero L.A. and Mayo-Iscar A. (2012). (FGM2012 in the
% code below)
%
%
%
%<a href="matlab: docsearchFS('restrdeter')">Link to the help function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit

% Examples:
%
%{
   %% Example using all default options.
   % Suppose v=3 and k=4 so the matrix containing the eigenvalues is 3-by-4
   % First column of matrix eigenvalues contains the eigenvalues of the first group
   % Second column of matrix eigenvalues contains the eigenvalues of the second group
   % Thrid column of matrix eigenvalues contains the eigenvalues of the third group
   % Fourth column of matrix eigenvalues contains the eigenvalues of the fourth group
   rng(10,'twister')
   eigenvalues=abs(10*randn(3,4));
   % niini is the column vector containing the sizes of the 4 groups
   niini=[30;40;20;10];
   out=restrdeter(eigenvalues,niini,1.1)
   disp('Input matrix of unrestricted eigenvalues')
   disp(eigenvalues)
   disp('Output matrix of restricted eigenvalues which satisfy determinant constraint')
   disp(out)
   disp('Ratio between largest and smallest determinant')
   disp(max(prod(eigenvalues))/min(prod(eigenvalues)))
   disp('Ratio between largest and smallest restricted determinants')
   disp(max(prod(out))/min(prod(out)))
%}


%{
    %% Determinant restriction when an eigenvalue is 0.
    % Suppose 5 variables and six groups
    av=abs(randn(5,6));
    % The third eigenvalue of the second groups is set to 0
    av(3,2)=0;
    % Maximum ratio among determinants must be equal to 1.6.
    restr=1.6;
    % Group sizes
    niini=[30;40;20;10;50;100];
    disp('Original values of the determinants')
    disp(prod(av))
    % Apply the restriction
    a=restrdeter(av,niini,restr);
    disp('Restricted eigenvalues which satisfy determinant constraint')
    disp(a)
    disp('Values of restricted determinants')
    disp(prod(a))
    disp('Maximum value of ratio among determinants')
    disp(max(prod(a))/min(prod(a)))
%}

%{
    % An example using option arguments tol and repmat.
    % Suppose 3 variables and six groups
    av=abs(randn(3,6));
    % Maximum ratio among determinants must be equal to 1.6.
    restr=1.6;
    % Group sizes
    niini=[30;40;20;10;50;100];
    % Apply the restriction using a tolerance of 1e-12 and use MATLAB
    % function repmat for the computations
    tol=1e-12;
    repm=1;
    a=restrdeter(av,niini,restr,tol,repm);
%}


%{
    % Determinant restriction when all eigenvalues of a group are 0.
    % Two variables and five groups.
    av=abs(randn(2,5));
    niini=[30;40;20;10;50];
    av(:,2)=0;
    a=restrdeter(av,niini,restr);
    disp('Maximum value of ratio among determinants')
    disp(max(prod(a))/min(prod(a)))
%}

%{
    % Determinant restriction when all eigenvalues of two groups are 0.
    av=abs(randn(2,5));
    av(:,2:3)=0;
    a=restrdeter(av,niini,restr);
    disp('Maximum value of ratio among determinants')
    disp(max(prod(a))/min(prod(a)))
%}

%% Beginning of code

if nargin<4
    tol =1e-8;
end

% userepmat specifies if it is necessary to use function repmat or bsxfun
% Remark: repmat has become built in from Release 2013b so it is faster to
% use it
if nargin<5
    userepmat=0;
end


% Get number of variables (v) and number of clusters (k)
[v,k]=size(eigenvalues);

eigenvalues(eigenvalues<1e-15)=0;
% Eigenvalue restriction using a restriction facto 10^10 is initially
% applied separately to each group
for j=1:k
    eigenvalues(:,j) = restreigen(eigenvalues(:,j), 1, 10^10,tol,userepmat);
end
% product of the elements of all columns
es=prod(eigenvalues,1);
es(es==0)=1;
gm=bsxfun(@rdivide,eigenvalues,es.^(1/v));

d=sum(eigenvalues./gm,1)/v;
d(isnan(d))=0;

dfin=restreigen(d,niini,restr^(1/v),tol,userepmat);
gm(gm==0)=1;
out=repmat(dfin,v,1).*gm;

end
%FScategory:CLUS-RobClaMULT