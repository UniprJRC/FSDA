function [C,nselected] = subsets(nsamp,n,p,ncomb,msg)
%subsets creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements
%
%<a href="matlab: docsearch('subsets')">Link to the help function</a>
%
%  Required input arguments:
%
%       nsamp : Number of subsamples which have to be extracted
%               Remark: if nsamp=0 all subsets will be extracted.
%               They will be (n choose p).
%         n   : scalar, number of observations of the dataset
%         p   : size of the subsets. In regression with p explanatory
%               variable the size of the elmental subsets is p. In
%               multivariate analysis in presente of v variables, the size
%               of the elemental subsets is v+1.
%
%  Optional input arguments:
%
%       ncomb : scalar (n choose p). If the user has already computed this
%               value it can supply it directly, otherwise the program will
%               calculate it automatically
%        msg  : scalar which controls whether to display or not messages
%               on the screen. If msg=1 (default), messages are displayed
%               on the screen about estimated time
%
%  Output:
%
%
%           C = matrix with nselected rows and p columns (stored in int16 format)
%               which contains the index of the subsets which
%               need to be extracted.
%   nselected = scalar, number of rows of matrix C.
%
%
% See also randsampleFS.m, lexunrank.m, bc.m
%
% References: see references in randsampleFS.m, lexunrank.m and bc.m
%
%
% Copyright 2008-2014.
% Written by FSDA team
%
%
%<a href="matlab: docsearch('subsets')">Link to the help function</a>
% Last modified 08-Dec-2013
%
% Examples:
%
%
%{
       % Create a matrix which contains the indexes of 20 subsets
       % when n=100, p=3
       n=100;
       p=3;
       C=subsets(20,100,3)
%}

%{
        % Extract 80000 unique subsets
        C=subsets(80000,100,5);
        size(unique(C,'rows'))
%}


%% Beginning of code

seq=1:n;

if nargin<4;
    ncomb=bc(n,p);
    
end

if ncomb<nsamp
    disp(['Warning: number of subsets which have been chosen (' num2str(ncomb) ') are greater than possibile number (' num2str(nsamp) ')']);
    disp('All subsets are extracted')
    nsamp=0;
end


if nargin<5
    msg=1;
end

%% Combinatorial part to extract the subsamples
% Key combinatorial variables used:
% C = matrix containing the indexes of the subsets (subsamples) of size p
% nselected = size(C,1), the number of all selected subsamples.
% nselected = number of combinations on which the procedure is run.
% rndsi = vector of nselected indexes randomly chosen between 1 e ncomb.
Tcomb = 5e+7; T2comb = 1e+8;

if nsamp==0 || ncomb <= Tcomb
    if nsamp==0
        if ncomb > 100000 && msg==1
            disp(['Warning: you have specified to extract all subsets (ncomb=' num2str(ncomb) ')']);
            disp('The problem is combinatorial and the run may be very lengthy');
        end
        nselected = ncomb;
    else
        nselected = nsamp;
    end
    
    % If nsamp = 0 matrix C contains the indexes of all possible subsets
    C=combsFS(seq,p);
    
    
    % If nsamp is > 0 just select randomly ncomb rows from matrix C
    if nsamp>0
        % Extract without replacement nsamp elements from ncomb
        rndsi=randsampleFS(ncomb,nsamp,2);
        C = C(rndsi,:);
    end
    
else
    if nsamp > 100000 && msg==1
        disp('Warning: you have specified to extract more than 100000 subsets');
        disp('To iterate so many times may be lengthy');
    end
    nselected = nsamp;
    
    usepascal=0;
    
    if ncomb>Tcomb && ncomb<T2comb;
        
        % Extract without replacement nsamp elements from ncomb
        rndsi=randsampleFS(ncomb,nsamp,2);
        
        if ispc;
            
            [~,sys]=memory;
            bytesavailable=sys.PhysicalMemory.Available;
            if bytesavailable > 2*8*n^2;
                pascalM=pascal(n);
            else
                pascalM=pascal(n);
            end
            usepascal=1;
            
        end
        
    end
    
    % Create matrix C which will contain in each row the indexes forming the
    % subset which is extracted at step i, where i=1....number of selected
    % subsamples (nselected)
    if n < 2^15
        C=zeros(nselected,p,'int16');
    else
        C=zeros(nselected,p,'int32');
    end
    
    for i=1:nselected
        
        if ncomb>Tcomb && ncomb<T2comb
            
            if usepascal
                s=lexunrank(n,p,rndsi(i),pascalM);
            else
                s=lexunrank(n,p,rndsi(i));
            end
        else
            s=randsampleFS(n,p);
        end
        C(i,:)=s;
        
    end
end

end




