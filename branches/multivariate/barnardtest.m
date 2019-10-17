function pval=barnardtest(N,varargin)
% Barnard's unconditional test.
%
%<a href="matlab: docsearchFS('barnardtest')">Link to the help function</a>
%
%
% This function computes the Barnard's unconditional test.
% The Barnard test is a powerful alternative
% of Fisher's exact test for 2x2 contingency tables.
%
%
%
%  Required input arguments:
%
%       N    :    Contingency table (default) or n-by-2 input dataset.
%                 Matrix or Table.
%                 Matrix or table which contains the input contingency
%                 table (say of size I-by-J) or the original data matrix.
%                 In this last case N=crosstab(N(:,1),N(:,2)). As default
%                 procedure assumes that the input is a contingency table.
%
%  Optional input arguments:
%
%   resolution: The resolution of the search space for the nuisance
%               parameter.
%               Scalar. Small number which defines the resolution. See the
%               More About section for more details.
%               The default value of the resolution is 0.001.
%               Example - 'resolution',0.01
%               Data Types - single | double
% datamatrix :  Data matrix or contingency table. Boolean.
%               If datamatrix is true the first input argument N is forced
%               to be interpreted as a data matrix, else if the input
%               argument is false N is treated as a contingency table. The
%               default value of datamatrix is false, that is the procedure
%               automatically considers N as a contingency table.
%               Example - 'datamatrix',true
%               Data Types - logical
%
%
%  Output:
%
%         pval:  p-value of the test. Scalar.
%                pval is the p-value, i.e. the probability of
%                observing the given result, or one more extreme, by
%                chance if the null hypothesis of independence between
%                rows and columns is true. Small values of pval cast doubt
%                on the validity of the null hypothesis.
%
% More About:
%
% For a 2x2 contingency table, such as $N=[n_{11},n_{12};n_{21},n_{22}]$, the normalized
% difference in proportions between the two categories, given in each
% column, can be written with pooled variance (Score statistic) as
% \[
% T(X)=\frac{\hat{p}_2-\hat{p}_1}{\sqrt{\hat{p}(1-\hat{p})(\frac{1}{c_1}+\frac{1}{c_2})}},
% \]
%
% where $\hat{p}=(n_{11}+n_{21})/n$ , $\hat{p}_2=n_{12}/(n_{12}+n_{22})$,
% $\hat{p}_1=n_{11}/(n_{11}+n_{21})$, $c_1=n_{11}+n_{21}$ and $c_2=n_{12}+n_{22}$.
%
% The probability of observing $N$ (the input contingency table) is
%
% \[
% P(N)=\frac{c_1!c_2!}{n_{11}!n_{12}!n_{21}!n_{22}!} p^{n_{11}+n_{12}}(1-p)^{n_{21}+n_{22}},
% \]
%
% where $p$ is the unknown nuisance parameter.
%
% Barnard's test considers all tables with category sizes $c_1$ and $c_2$ for a
% given $p$. The p-value is the sum of probabilities of the tables having a
% score in the rejection region, e.g. having significantly large difference
% in proportions for a two-sided test. The p-value of the test is the
% maximum p-value calculated over all $p$ between 0 and 1. The input
% resolution parameter controls the resolution to search for.
%
%
% See also crosstab, fishertest, CressieRead, corrNominal
%
% References:
%
% Barnard, G.A. (1945), A new test for 2x2 tables, "Nature", pp. 156-177.
% Barnard, G.A. (1947), Significance tests for 2x2 tables, "Biometrika",
% Vol. 34, pp. 123-138.
% Suissa, S. and Shuster, J.J. (1985), Exact Unconditional Sample Sizes
% for the 2x2 Binomial Trial, "Journal of the Royal Statistical Society",
% Ser. A, Vol. 148, pp. 317-327.
% Lin, C.Y., Yang, M.C. (2009), Improved p-value tests for comparing two
% independent binomial proportions, "Communications in Statistics-Simulation
% and Computation", Vol. 38, pp.78-91.
%
%
%
% Acknowledgements:
%
% This file was inspired by Trujillo-Ortiz, A., R. Hernandez-Walls, A. Castro-Perez, L.
% Rodriguez-Cardozo N.A. Ramos-Delgado and R. Garcia-Sanchez. (2004).
% Barnardextest:Barnard's Exact Probability Test.
% https://www.mathworks.com/matlabcentral/fileexchange/6198-barnardextest .
% and  by Cardillo G. (2009) MyBarnard: a very compact routine for Barnard's exact
% test on 2x2 matrix.
% http://www.mathworks.com/matlabcentral/fileexchange/25760 .
% A comparison with the the current implementation is given in the last
% example.
% The FSDA team wishes to thank Dr. Ivano Azzini for the current implementation
% of the Barnard test.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('barnardtest')">Link to the help function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit
%
%

% Examples:

%{
    %%  Barnard test with all the default options.
    % N= 2x2 Input contingency table
    N=[8,14; 1,3];
    pval=barnardtest(N);
    disp(['The p-value of the test is: ' num2str(pval)])
%}

%{
    %%  Resolution option.
    % N= 2x2 Input contingency table
    N=[20,14; 10,13];
    % pvalue with the default resolution (0.001)
    pval001=barnardtest(N);
    % p value with a resolution of 0.01
    pval01=barnardtest(N,'resolution',0.01);
    disp(['The p-value with a resolution 0.01 is: ' num2str(pval01)])
    disp(['The p-value a resolution 0.001 is: ' num2str(pval001)])
%}


%{
    % An example when the input is a MATLAB table.
    rownam={'OutcomeI', 'OutcomeII'};
    colnam={'TreatmentI' 'TreatmentII'};
    N=[40,14;10,30];
    if verLessThan('matlab','8.2.0') ==0
        Ntable=array2table(N,'RowNames',rownam,'VariableNames',colnam);
    else
        Ntable=N;
    end
    pval=barnardtest(Ntable);
%}

%{
    % An example when the input is a datamatrix.
    N=[40,14;10,30];
    % Recreate the orginal data matrix
    X=crosstab2datamatrix(N);
    % barnardtest when input is a datamatrix
    pval=barnardtest(X,'datamatrix',true);
%}

%{
    % Comparison with other existing implementations.
    % Using the example below
    N=[16 40; 1 2];
    pval=barnardtest(N);
    % our p-value is 0.456054
    % This value coincides with the R implementation (package barnard)
    % based on a C routine. On the other hand the vectorized implementation
    % of Barnard test http://www.mathworks.com/matlabcentral/fileexchange/25760
    % called using mybarnard(N,1000) gives a p-value of   0.456051 
%}


%% Beginning of code

% Check MATLAB version. If it is not smaller than 2013b than output is
% shown in table format
verMatlab=verLessThan('matlab','8.2.0');

% Check whether N is a contingency table or a n-by-p input dataset (in this
% last case the contingency table is built using the first two columns of the
% input dataset).
if ~isempty(varargin)
    UserOptions=varargin(1:2:length(varargin));
    checkdatamatrix = strcmp(UserOptions,'datamatrix');
    if sum(checkdatamatrix)
        datamatrix = varargin{2*find(checkdatamatrix)};
    else
        datamatrix=false;
    end
else
    datamatrix=false;
end

% If input is a datamatrix it is necessary to construct the contingency
% table
if datamatrix == true
    N =crosstab(N(:,1),N(:,2));
end

% 1. The deafault resolution of the search space.
resolution=0.001;

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    options=struct('datamatrix',false,...
        'resolution',resolution);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:barnardtest:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    resolution=options.resolution;
end

% transform table in array
if verMatlab ==0 && istable(N)
    N=table2array(N);
end

a=N(1,1);
b=N(1,2);
c=N(2,1);
d=N(2,2);


% Contingency Table columns sums.
c1 = a+c;
c2 = b+d;

% Some useful f(c1,c2).
C=c1+c2;
IC=(1/c1)+(1/c2);
PC=c1*c2;

% The initial normalized difference in proportions
% between the two categories, for each column (TX0).

pao = a/c1;
pbo = b/c2;
pxo = (a+b)/C;

TX0 = abs(pao-pbo)/sqrt(pxo*(1-pxo)*IC);
if isnan(TX0)
    TX0=0;
end


% PREPARE data for the FOR cycles.
P = zeros(ceil(1/resolution),1);
iP=1;

% 2. The required factorials.
%    ... Transformed in logarithms.
cMax=max(c1,c2);
F=zeros(cMax,1);
for i=1:cMax+1
    F(i)=sum(log(2:i-1));
end
F1=F(1:c1+1);
F2=F(1:c2+1);
F1end=F1(end);
F2end=F2(end);

% 3. The required integers.
cInt=0:C;
pa=cInt(1:c1+1)/c1;
pb=cInt(1:c2+1)/c2;
px=cInt/C;

% 4. To store the #'p' S and TX.
S = zeros(PC,1);
TX =  S;

c2seq=(0:c2)';
c2seq1=c2seq+1;

% FOR cycles.
for p=0:resolution:1
    %  rI=1;
    
    logp=log(p);
    log1mp=log(1-p);
    
    for i = 0:c1
        
        % The loop below is replaced by a series of vectorized instructions
        %         for j = 0:c2
        %             if p>0
        %                 S(rI) = exp(F1(end)-F1(end-i)-F1(i+1)+ F2(end)-F2(end-j)-F2(j+1) +(i+j)*log(p)+(C-(i+j))*log(1-p));
        %             else
        %                 S(rI)=0;
        %             end
        %
        %             % Check possible division by 0, etc. (Locally are irrelevant
        %             % for the final result see #1).
        %             if isnan(S(rI)) || (S(rI)==Inf)
        %                 S(rI)=0;
        %             end
        %
        %             TX(rI)=(pa(i+1)-pb(j+1))/sqrt(px(i+j+1)*(1-px(i+j+1))*IC);
        %
        %             rI=rI+1;
        %         end
        
        c2seqi=c2seq+i;
        q2store=exp(F1end-F1(end-i)-F1(i+1)+ F2end-F2(end-c2seq)-F2(c2seq1) +(c2seqi)*logp+(C-(c2seqi))*log1mp);
        % Check possible division by 0, etc. (Locally are irrelevant
        % for the final result see #1).
        % It seems that ~isfinite is slightly faster than
        % q2store(isnan(q2store) | isinf(q2store))=0;
        q2store(~isfinite(q2store))=0;
        
        indexes=i*(c2+1)+1:+(i+1)*(c2+1);
        S(indexes) = q2store;
        
        % Schk(i*(c2+1)+1:+(i+1)*(c2+1)) = exp(F1(end)-F1(end-i)-F1(i+1)+ F2(end)-F2(end-c2seq)-F2(c2seq+1) +(i+c2seq)*log(p)+(C-(i+c2seq))*log(1-p));
        pxic2=px(i+c2seq1);
        TX(indexes)= (pa(i+1)-pb(c2seq1))./sqrt(pxic2.*(1-pxic2)*IC);
        
    end
    P(iP) = sum(S(TX>=TX0));
    iP=iP+1;
end

pval=max(P);

end
%FScategory:MULT-Categorical