function [ABk,Bk,EBk,VarBk] = FowlkesMallowsIndex(c1,c2)
%FowlkesMallowsIndex computes the Fowlkes and Mallows index.
%
%<a href="matlab: docsearchFS('FowlkesMallowsIndex')">Link to the help function</a>
%
% Fowlkes-Mallows index (see references) is an external evaluation method
% that is used to determine the similarity between two clusterings
% (clusters obtained after a clustering algorithm). This measure of
% similarity could be either between two hierarchical clusterings or a
% clustering and a benchmark classification. A higher the value for the
% Fowlkes-Mallows index indicates a greater similarity between the clusters
% and the benchmark classifications.
% This index can be used to compare either two cluster label sets or a
% cluster label set with a true label set. The formula of the
% adjusted Fowlkes-Mallows index (ABk) is given below
% \[
%  ABk= \frac{\mbox{Bk- Expected value of Bk}}{\mbox{Max Index  - Expected value of Bk}}
% \]
%
%
% Required input arguments:
%
%    c1:   labels of first partition or contingency table. 
%          Numeric or character vector. 
%          A numeric or character vector containining the class labels of
%          the first partition or a 2-dimensional numeric matrix which
%          contains the cross-tabulation of cluster assignments.
%          Data Types: single | double | char | logical
%    c2:   labels of second partition. 
%          Numeric or character vector. 
%          A numeric or character vector containining the class labels of
%          the second partition. The length of vector c2 must be equal to
%          the length of vector c1. This second input is required just if
%          c1 is not a 2-dimensional numeric matrix.
%          Data Types: single | double | char | logical
%
%Optional input arguments:
%
%   noisecluster: label or number associated to the 'noise class' or 'noise level'.
%                 Scalar, numeric or character.
%                 Number or character label which
%                 denotes the points which do not belong to any cluster.
%                 These points are not takern into account for the
%                 computation of the Fowlkes and Mallows index
%                 Example - 0 (in this case the units which in of the
%                 two partitions have 0 class are not taken into account in the
%                 index calculations)
%                 Data Types - double or character
%
% Output:
%
%  ABk  :       Adjusted Fowlkes and Mallows index. Scalar. A number between -1 and 1.
%               The adjusted Fowlkes and Mallows index is the
%               corrected-for-chance version of the Fowlkes and Mallows
%               index.
%  Bk:          Value of the Fowlkes and Mallows index. Scalar. A number between 0 and 1.            
%  EBk:         Expectation of the Fowlkes and Mallows index. Scalar.
%               Expected value of the index computed under the null
%               hypothesis of no-relation.
%  VarBk:        Variance of the Fowlkes and Mallows index. Scalar.
%               Variance of the index computed under the null
%               hypothesis of no-relation.
%
%
% See also: RandIndexFS
%
% References:
%
% Fowlkes, E.B. and Mallows, C.L. (1983), A Method for Comparing Two
% Hierarchical Clusterings, "Journal of the American
% Statistical Association", Vol. 78, pp. 553-569.
% [http://en.wikipedia.org/wiki/Fowlkes-Mallows_index]
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FowlkesMallowsIndex')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{

    % FowlkesMallowsIndex (adjusted) with the two vectors as input.
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
    FM=FowlkesMallowsIndex(c1,c2);
%}

%{
    % FM index (adjusted) with the contingency table as input.
    T=[1 1 0;
    1 2 1;
    0 0 4];
    FM=FowlkesMallowsIndex(T);
%}

%{
        % Compare FM (unadjusted) for iris data (true classification against tclust classification).
        load fisheriris
        % first partition c1 is the true partition
        c1=species;
        % second partition c2 is the output of tclust clustering procedure
        out=tclust(meas,3,0,100,'msg',0);
        c2=out.idx;
        [~,FM,EFM,VARFM]=FowlkesMallowsIndex(c1,c2);
%}

%{
        % Compare FM index (unadjusted) for iris data (exclude unassigned units from tclust).
        load fisheriris
        % first partition c1 is the true partition
        c1=species;
        % second partition c2 is the output of tclust clustering procedure
        out=tclust(meas,3,0.1,100,'msg',0);
        c2=out.idx;
        % Units inside c2 which contain number 0 are referred to trimmed observations
        noisecluster=0;
        [~,FM,EFM,VARFM]=RandIndexFS(c1,c2,noisecluster);
%}

%{
    % FM index (unadjusted) for iris data with 3 groups coming from single linkage.
    % FM index between true and empirical classification
    load fisheriris
    d = pdist(meas);
    Z = linkage(d);
    C = cluster(Z,'maxclust',3);
    [AFM,FM,FMexp,FMvar]=FowlkesMallowsIndex(C,species);
    disp('FM index is equal to')
    disp(FM)
    disp('Expectation of FM index is')
    disp(FMexp)
    disp('Variance of FM index is')
    disp(FMvar)
    disp('Adjsuted FM index is equal to')
    disp(AFM)
%}

%{
    % Monitoring of (adjusted) FM index for iris data using true classification as benchmark.
    load fisheriris
    d = pdist(meas);
    Z = linkage(d);
    kk=1:15;
    % Produce agglomerative hierarchical cluster tree
    C = cluster(Z,'maxclust',kk);
    FM =zeros(length(kk)-1,1);
    for j=kk
        FM(j)=FowlkesMallowsIndex(C(:,j),species);
    end
    plot(kk,FM)
    xlabel('Number of groups')
    ylabel('Fowlkes and Mallows Index')
%}

%% Beginning of code

if nargin < 2 || isempty(c2)
    if size(c1,2)<2
        error('FSDA:FowlkesMallowsIndex:InvalidArg','FowlkesMallowsIndex: Requires a contingency table with at least two columns')
    end
    % Supplied input is the contingency table
    M=c1;
else
    
    if min(size(c1)) > 1 || min(size(c2)) > 1
        error('FSDA:FowlkesMallowsIndex:InvalidArg','FowlkesMallowsIndex: Requires two vector arguments')
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
    
    if iscell(c1) || iscell(c2)
        M=crosstab(c1,c2);
    else
        class_c1 = sort(unique(c1), 'ascend');
        class_c2 = sort(unique(c2), 'ascend');
        
        Rows = numel(class_c1);
        Cols = numel(class_c2);
        
        % M is the contingency table.
        % The loop below is faster than the use of crosstab
        M = zeros(Rows, Cols);
        for i = 1:Rows
            for j = 1:Cols
                % Function nnz is the number o nonzero matrix elements
                M(i,j) = nnz(c1(:) == class_c1(i) & c2(:) == class_c2(j));
            end
        end
    end
end

%midot = vector containing in position i (i=1, ...,R) m_{i.}=\sum_{j=1^C} m_{ij}
midot = sum(M,2);
%mdotj = vector containing in position j j=1, ...,C m_{.j}=\sum_{i=1^R} m_{ij}  ,
mdotj = sum(M,1);

n=length(c1);
Tk = sum(M(:).^2) - n;

Pk = sum(midot.^2) - n;
Qk = sum(mdotj.^2) - n;
Bk=Tk/sqrt(Pk*Qk);

% Ebk = Expected value
EBk=sqrt(Pk*Qk)/(n*(n-1));

% VArbk = variance
Pk2 = sum(midot.*(midot - 1).*(midot - 2));
Qk2=  sum(mdotj.* (mdotj - 1).*(mdotj - 2));
VarBk = 2/(n * (n - 1)) + 4 * Pk2 * Qk2/((n * (n - 1) * (n -2)) * Pk * Qk)...
    + (Pk - 2 - 4 * Pk2/Pk) * (Qk - 2 - 4 *Qk2/Qk)/((n * (n - 1) * (n - 2) * (n - 3)))...
    - Pk * Qk/(n^2 * (n-1)^2);

% Compute adjusted version of the index
ABk = (Bk - EBk)/(1- EBk);
end

%FScategory:UTISTAT
