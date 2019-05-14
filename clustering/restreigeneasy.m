function [out]  = restreigeneasy(eigenvalues, niini, restr, tol)
%restreigen computes eigenvalues restriction (without Dykstra algorithm)
%
%<a href="matlab: docsearchFS('restreigeneasy')">Link to the help function</a>
%
%   restreigeneasy restricts the eigenvalues according to the constraint
%   specified in scalar restr. This function is a readable and understandble
%   (but slow)  version of routine restreigen. Routine restreigen is called
%   in any concentration step of function tclust and can also be used
%   inside function MixSim to generate groups with a prespecified level of
%   overlapping
%
%  Required input arguments:
%
%eigenvalues: Eigenvalues. Matrix. v x k matrix containing the eigenvalues
%             of the covariance matrices of the k groups.
%             v is the number of variables of the dataset which has to be
%             clustered.
%     niini: Cluster size. Vector. k x 1 vector containing the size of the k clusters
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
%
%  Output:
%
%
%            out      : Restricted eigenvalues. Matrix. v-by-k matrix
%                       containing restricted eigenvalues.
%                       The ratio between two possible elements in matrix
%                       out is not greater than restr
%
% See also restreigen, tclust, restrdeter, tclustreg
%
% References:
%
% This function implements the algorithm described in
% Fritz H., Garcia-Escudero, L.A. and Mayo-Iscar, A. (2013), A fast
% algorithm for robust constrained clustering,
%"Computational Satistics and Data Analysis", Vol. 61, pp. 124-136.
% [Available at
% http://www.eio.uva.es/infor/personas/tclust_algorithm.pdf]
%
% Copyright 2008-2019.
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
%<a href="matlab: docsearchFS('restreigeneasy')">Link to the help function</a>
%
%$LastChangedDate:: 2018-01-24 13:44:34 #$: Date of the last commit

% Examples:
%
%{
   % Example using all default options.
   % Suppose v=3 and k=4 so the matrix containing the eigenvalues is 3-by-4
   % First column of matrix eigenvalues contains the eigenvalues of the first group
   % Second column of matrix eigenvalues contains the eigenvalues of the second group
   % Thrid column of matrix eigenvalues contains the eigenvalues of the third group
   % Fourth column of matrix eigenvalues contains the eigenvalues of the fourth group
   rng(10,'twister')
   eigenvalues=abs(10*randn(3,4));
   % niini is the vector containing the sizes of the 4 groups
   niini=[30;40;20;10];
   out=restreigeneasy(eigenvalues,niini,1.1)
   disp('Input matrix of unrestricted eigenvalues')
   disp(eigenvalues)
   disp('Output matrix of restricted eigenvalues')
   disp(out)
   disp('Ratio between largest and smallest unrestricted eigenvalues')
   disp(max(max(eigenvalues))/min(min(eigenvalues)))
   disp('Ratio between largest and smallest restricted eigenvalues')
   disp(max(max(out))/min(min(out)))
%}


%{
    % Compare speed.
    % We compare the speed of restreigneasy with that of restreigen. We use
    % userepmat=2 if the current MATLAB version if >=R2017a or userepmat =1
    % if MATLAB version is >=R2013a but <R2017a else we use userepmat =0
    v=10;
    k=8;
    tol=1e-8;

    if verLessThanFS(9.2)== false
        % If MATLAB version is at least 2017a
        userepmat=2;
    elseif verLessThanFS(8.1) == false
        % if MATLAB version is at least R2013b  
        userepmat=1;
    else
        userepmat=0;
    end

    oldroutinetime=0;
    newroutinetime=0;
    rng(1)
    for j=1:10000
        eigenvalues=100*abs(randn(v,k));
        % niini is the vector containing the sizes of the 4 groups
        niini=randi(100,[k,1]);
        tic;
        outold=restreigeneasy(eigenvalues,niini,1.1);
        % Uncomment the line below if you want 
        % outold=restreigen(eigenvalues,niini,1.1,tol,1);
        oldroutinetime=oldroutinetime+toc;
        tic;
        outnew=restreigen(eigenvalues,niini,1.1,tol,userepmat);
        newroutinetime=newroutinetime+toc;
        if max(max(abs(outold-outnew)))>1e-5
            error('The two routines are different')
        end
    end
    disp(['Computing time of restreigeneasy: ' num2str(oldroutinetime)])
    disp(['Computing time of restreigen: ' num2str(newroutinetime)])
%}



%% Beginning of code

if nargin<4
    tol =1e-8;
end

% Initializations
c=restr;
d=eigenvalues';

% Get number of variables (v) and number of clusters (k)
[v,k]=size(eigenvalues);


% Get number of units of the dataset
n=sum(niini);

% We assume that niini is a column vector
% nis is a matrix which replicates in the columns the sizes of the goups
% First row of nis = size of first group repated v times
% Second row of nis = size of second group repated v times
% ....
% kth row of nis = size of kth group repated v times

nis=repmat(niini,1,v);

% Below is the alternative inefficient code
% nis = repmat(niini,1,v);


% niini=niini';

% dsor is 2*k*v the ordered set of values in which the restriction objective
% function change the definition The elements in dsor correspond to  the
% frontiers for the intervals in which this objective function has the same
% definition
% In other words
% dsor=(d_{11}, ........, d_{kv},d_{11}/restr, ........, d_{kv}/restr)

% OLD was
% dsor=sort([eigenvalues(:);eigenvalues(:)/c])';
eigvector=eigenvalues(:);
dsor=sort([eigvector/c;eigvector])';
kv=k*v;
% dimsor=length(dsor);
dimsor=kv*2;


% d1 is like dsor but contains an additional element which is larger than the largest element of dsor
d1=dsor;
d1(dimsor+1)=d1(dimsor)*2;

% d2 is like dsor but contains an additional element which smaller than the smallest element of dsor
d2=[0,dsor];


% ed is a set with the middle points of these intervals
ed=(d1+d2)/2;
dimsor=dimsor+1;

% the only relevant eigenvalues are those belong to a clusters with sample
% size greater than 0. eigenvalues corresponding to a cluster with 0
% elements has no influence in the objective function
dnis=d(nis>0);

maxdnis=max(dnis);


if maxdnis <= tol
    % if all the eigenvalues are 0 this means all points are concentrated
    % in k groups and there is a perfect fit
    % no further changes on the eigenvalues required, so return them
    % immediately and stop the procedure!
    out = eigenvalues;
else
    % we check if the  eigenvalues verify the restrictions
    
    % abs here is just for computational purposes
    if abs(maxdnis/min(dnis))<=c
        % If all eigenvalues satisy the constraint
        % no further changes on the eigenvalues required, so return them immediately!
        % Simply replace the 0 eigenvalues with the mean of the eigenvalues
        % which are greater than zero
        if min(nis)==0
            d(nis==0)=mean(dnis);
        end
        out=d';
    else
        % ----------------------------------------------------------------------
        % t, s and r are matrices of size k-times-2*k*v+1
        % Each column of matrices r s and t is associated to a candidate
        % value of m
        % Each row of matrices r s and t is associated to a group
        % The final candidate value of m will be obtained as a weighted sum
        % of the rows (the weights are given by the group sizes)
        % Each element of r, s and t respectively contains
        % r = \sum_{l=1}^v (d_l <m) + \sum_{l=1}^v (d_l >cm)
        % s = \sum_{l=1}^v (d_l <m)
        % t = \sum_{l=1}^v (d_l >c m)
        t=zeros(k,dimsor);
        s=t;
        r=t;
        
        % sol and obj are two column vectors of size 2*k*v+1
        % sol contains all the candidates solutions of m
        % Vector sol contains the
        % critical values of the interval functions which define the m
        % objective function we use the centers of the interval to get a
        % definition for the function in each interval this set with the
        % critical values (in the array sol) contains the optimum m value
        sol=zeros(dimsor,1);
        
        % obj is the vector which contains the function which must be
        % minized
        % There is one value of the objective function for each candidate
        % solution m
        % The optimum value of obj is the smallest
        % The optmimum m is the one associated to the smallest value of
        % function obj
        
        obj=sol;
        % The candidates values of mp which minimize the objective function are
        dimsor=2*k*v+1;
        for mp =1:dimsor
            edmp=ed(mp);
            edmpc=edmp*c;
            
            % Computation of r s and t
            % Note that the sum goes from 1 to v
            r(:,mp)=sum(d<edmp,2)+sum(d>edmpc,2);
            s(:,mp)=sum(d.*(d<edmp),2);
            t(:,mp)=sum(d.*(d>edmpc),2);
            
            % The natural loop is clearly slower
            %             for i=1:k
            %                 di=d(i,:);
            %                 r(i,mp)=sum(di<edmp)+sum(di>edmpc);
            %                 s(i,mp)=sum(di.*(di<edmp));
            %                 t(i,mp)=sum(di.*(di>edmpc));
            %             end
            
            
            % Note that here the sum goes from 1 to k
            % sol(mp) \sum_{j=1}^k n_j (s_j +t_j/c)  / \sum_{j=1}^k n_j r_j
            solmp=sum(niini/n.*(s(:,mp)+t(:,mp)/c))/(sum(niini/n.*(r(:,mp))));
            sol(mp)=solmp;
            
            % solnum(mp)=sum(niini/n.*(s(:,mp)+t(:,mp)/c));
            % solden(mp)=(sum(niini/n.*(r(:,mp))));
            
            e = solmp*(d<solmp)+d.*(d>=solmp).*(d<=c*solmp)+(c*solmp)*(d>c*solmp);
            o=nis/n.*(log(e)+d./e);
            
            % sol (mp) contains \sum_{j=1}^k \sum_{l=1}^v ( log d_{jl} + d_{jl}^m / d_{jl}
            % equation (3.4)
            obj(mp)=sum(sum(o));
        end
        
        [~,indmax]=min(obj);
        % m is the optimum value for the eigenvalues procedure
        m=sol(indmax);
        out= ((m*(d<m)+d.*(d>=m).*(d<=c*m)+(c*m)*(d>c*m)))';
        
    end
end
%FScategory:CLUS-RobClaMULT