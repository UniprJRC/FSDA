function [AR,RI,MI,HI]=RandIndexFS(c1,c2, noisecluster)
%RandIndexFS calculates Rand type Indices to compare two partitions
%
%<a href="matlab: docsearchFS('RandIndexFS')">Link to the help function</a>
%
% Suppose we want to compare two partitions summarized by the contingency
% table $T=[n_{ij}]$ where $i=1, 2, ..., r$ and $j=1,...,c$ and $n_{ij}$
% denotes the number of data points which are in cluster i in the first
% partition and in cluster j in the second partition. Let A denote the
% number of all pairs of data points which are either put into the same
% cluster by both partitions or put into different clusters by both
% partitions. Conversely, let D denote the number of all pairs of data
% points that are put into one cluster in one partition, but into different
% clusters by the other partition. The partitions disagree for all pairs D
% and agree for all pairs A. A+D=totcomp= total number of comparisons.
% We can measure the agreement by the Rand index A/(A+D)=A/(totcomp) which
% is invariant with respect to permutations of the columns or rows of T.
% The index has to be corrected for agreement by chance if the sizes of the
% clusters are not uniform (which is usually the case). Since the Rand
% index lies between 0 and 1, the expected value of the Rand index
% (although not a constant value) must be greater than or equal to 0. On
% the other hand, the expected value of the adjusted Rand index has value
% zero and the maximum value of the adjusted Rand index is also 1. Hence,
% there is a wider range of values that the adjusted Rand index can take
% on, thus increasing the sensitivity of the index. The formula of the
% adjusted Rand index (AR) is given below
% \[
%  AR= \frac{\mbox{RI- Expected value of RI}}{\mbox{Max Index  - Expected value of RI}}
% \]
%
%
% Required input arguments:
%
%    c1:   labels of first partition or contingency table. A numeric or
%          character vector containining the class labels of
%          the first partition or a 2-dimensional numeric matrix which
%          contains the cross-tabulation of cluster assignments.
%          Data Types: single | double | char | logical
%    c2:   labels of second partition. A numeric or character vector
%          containining the class labels of
%          the second partition. The length of vector c2 must be equal to
%          the length of vector c1. This second input is required just if
%          c1 is not a 2-dimensional numeric matrix.
%          Data Types: single | double | char | logical
%
%
%Optional input arguments:
%
%   noisecluster: label or number associated to the 'noise class' or 'noise level'.
%                 Scalar, numeric or character.
%                 Number or character label which
%                 denotes the points which do not belong to any cluster.
%                 These points are not takern into account for the
%                 computation of the Rand type indexes. The default is to
%                 consider all points in order to compute the ARI index.
%                 Example - 0 (in this case the units which in of the
%                 two partitions have 0 class are not taken into account in the
%                 index calculations)
%                 Data Types - double or character
%
% Output:
%
%  AR:          Adjusted Rand index. Scalar. A number between -1 and 1.
%               The adjusted Rand index is the corrected-for-chance version
%               of the Rand index.
%  RI:          Rand index (unadjusted). Scalar. A number between 0 and 1.
%               Rand index computes the fraction of pairs of objects for
%               which both classification methods agree.
%               RI ranges from 0 (no pair classified in the same way under
%               both clusterings) to 1 (identical clusterings).
%  MI:          Mirkin's index. Scalar. A number between 0 and 1.
%               Mirkin's index computes the percentage of pairs of objects for
%               which both classification methods disagree. MI=1-RI.
%  HI :         Hubert index. Scalar.  A number between -1 and 1.
%               HI index is equal to the fraction of pairs of objects for
%               which both classification methods agree minus the
%               fraction of pairs of objects for which both
%               classification methods disagree. HI= RI-MI.
%
%
% See also: crosstab, tclust
%
% References:
%
% Hubert L. and Arabie P. (1985), Comparing Partitions, "Journal of
% Classification", Vol. 2, pp. 193-218.
%
% Acknowledgements:
%
% This function follows the lines of MATLAB code developed by
% David Corney (2000) 	D.Corney@cs.ucl.ac.uk
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('RandIndexFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % RandindexFS with the contingency table as input.
    T=[1 1 0;
    1 2 1;
    0 0 4];
    ARI=RandIndexFS(T);
%}

%{

    % RandindexFS with the two vectors as input.
     c=[1 1;
        1 2
        2 1;
        2 2 ;
        2 2;
        2 3;
        3 3;
        3 3;
        3 3;
        3 3];
    % c1= numeric vector containing the labels of the first partition
    c1=c(:,1);
    % c1= numeric vector containing the labels of the second partition
    c2=c(:,2);
    ARI=RandIndexFS(c1,c2);
%}

%{
   % Computation of ARI, RI, MI and HI.
   [ARI,RI,MI,HI]=RandIndexFS(c1,c2);
    disp('Adjusted Rand index')
    disp(ARI)
    disp('Rand index (RI)')
    disp(RI)
    disp('Mirkin index = 1-RI')
    disp(MI)
    disp('Hubert index = RI-MI ')
    disp(HI)
%}

%{
        % Compare ARI for iris data (true classification against tclust classification).
        load fisheriris
        % first partition c1 is the true partition
        c1=species;
        % second partition c2 is the output of tclust clustering procedure
        out=tclust(meas,3,0,100,'msg',0);
        c2=out.idx;
        [ARI,RI,MI,HI]=RandIndexFS(c1,c2);
%}

%{
        % Compare ARI for iris data (exclude unassigned units from tclust).
        load fisheriris
        % first partition c1 is the true partition
        c1=species;
        % second partition c2 is the output of tclust clustering procedure
        out=tclust(meas,3,0.1,100,'msg',0);
        c2=out.idx;
        % Units inside c2 which contain number 0 are referred to trimmed observations
        noisecluster=0;
        [ARI,RI,MI,HI]=RandIndexFS(c1,c2,noisecluster);
%}

%% Beginning of code

if nargin < 2 || isempty(c2)
    if size(c1,2)<2
        error('FSDA:RandIndexFS:WrongInput','RandIndex requires a contingency table with at least two columns')
    end
    % Supplied input is the contingency table
    C=c1;
else
    
    if min(size(c1)) > 1 || min(size(c2)) > 1
        error('FSDA:RandIndexFS:WrongInput','RandIndex: Requires two vector arguments')
    end
    
    if nargin>2
        if ischar(c1) || iscell(c1)
            boo1=strcmp(c1,noisecluster);
        else
            boo1=c1==noisecluster;
        end
        if ischar(c2) || iscell(c2)
            boo2=strcmp(c2,noisecluster);
        else
            boo2=c2==noisecluster;
        end
        
        boo=~(boo1 |boo2);
        c1=c1(boo);
        c2=c2(boo);
    end
    C=crosstab(c1,c2);	% form contingency matrix
    
end

% From the contingency table it is necessary to remove the noise
% clusterclass


n=sum(sum(C));
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

% totcomp= total number of pairs of entities = total number of comparisons
% totcomp =n(n-1)/2
totcomp=nchoosek(n,2);
t2=sum(sum(C.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

% Expected value of the index (for adjustment)
% Hubert and Arabie proposed an adjustment which assumes a generalized
% hypergeometric distribution as null hypothesis: the two clusterings are
% drawn randomly with a fixed number of clusters and a fixed number of
% elements in each cluster (the number of clusters in the two clusterings
% need not be the same). Then the adjusted Rand Index is the (normalized)
% difference of the Rand Index and its expected value under the null
% hypothesis.
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));


% A = total number of agreements
A=totcomp+t2-t3;		%no. agreements
% D = total number of disagreements
D=  -t2+t3;		%no. disagreements
% Remark: totcomp=A+D = total number of comparisons

if totcomp==nc
    AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
    AR=(A-nc)/(totcomp-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/totcomp;			%Rand 1971		%Probability of agreement
MI=D/totcomp;			%Mirkin 1970	%p(disagreement)

%Hubert 1977
% prob of agreement - prob of disagreement
% p(agree)-p(disagree)
HI=(A-D)/totcomp;
end

%FScategory:UTISTAT


