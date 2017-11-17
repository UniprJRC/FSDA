function out = SparseTableTest(N,varargin)
%SparseTableTest computes independence test for large and sparse contingency tables
%
%<a href="matlab: docsearchFS('SparseTableTest')">Link to the help function</a>
%
% This function implements a new test of indipendence between row variables
% distribution ('outcomes') and columns ('treatments') which is expecially
% suited for the analysis of large and sparse $I$-by-$J$ contingency
% tables. The procedure is based on the collapsing of the original table
% into a set of 2-by-2 tables for each cell of the original table which has
% no less than a small number of counts (set in the optional input
% parameter 'threshold') and testing each of the resulting collapsed tables
% for independence by any test (Fisher exact test (default), Barnard test
% or those belonging to the power divergence family of Cressie and Read).
% Because of the Bonferroni inequality, a sufficient condition for
% attaining a significance level $\alpha$ for this test (i.e., the
% probability of detecting a positive association between two levels of the
% response variables when in fact there is no such association) is that
% each test done for each cell of the $I$-by-$J$ table rejects with
% significance level equal to $\alpha$ divided by the number of comparisons
% done. An additional bonus of the procedures is that it enables to
% highlight the most important contribution to the association of each
% single entry of the original I-by-J-table two way table. The original
% idea of this test is due to Spyros Arsenis (Joint Research Centre of the
% European Commission) and has been successfully applied to the analysis of
% contingency table coming from international trade data.
%
%  Required input arguments:
%
%       N    :    Contingency table (default) or n-by-2 input dataset.
%                 Matrix or Table. Matrix or table which contains the input
%                 contingency table (say of size I-by-J) or the original
%                 data matrix. In this last case N=crosstab(N(:,1),N(:,2)).
%                 By default the procedure assumes that the input is a
%                 contingency table.
%
%  Optional input arguments:
%
% threshold  : Threshold to select collapsed contigencey tables. Scalar.
%              Scalar which specifies above which value collapsed
%              contingency tables have to be produced. The default value of
%              threshold is 2.
%              Example - 'threshold',3
%              Data Types - single | double | int32 | int64
%     alpha  : Significance level. Scalar value in the range (0,1).
%              Significance level of the hypothesis test, specified as the
%              comma-separated pair consisting of 'alpha' and a scalar
%              value in the range (0,1). The default value of alpha is
%              0.01.
%              Example - 'alpha',0.05
%              Data Types - single | double
%  testname :  Test to use on collapsed 2-by-2 tables. Char or double.
%              If testname is a number, it identifies the value of $\lambda$
%              to use of the power divergence family. See function
%              CressieRead.m for further details. If testname is a
%              character, possible values are 'Fisher' (to use the Fisher
%              exact test, see function fishertest) or 'Barnard' (to use
%              Barnard exact test, see function barnardtest). The default
%              value of testname is 1, that is $\chi^2$ test is used. Note
%              also that fishertest has been introduced in MATLAB in
%              release 2014b.
%              Example - 'testname',1
%              Data Types - single | double | char
% datamatrix : Data matrix or contingency table. Boolean. If datamatrix
%              is true the first input argument N is forced to be
%              interpreted as a data matrix, else if the input argument is
%              false N is treated as a contingency table. The default value
%              of datamatrix is false, that is the procedure automatically
%              considers N as a contingency table
%              Example - 'datamatrix',true
%              Data Types - logical
%  Output:
%
% out:   structure which contains the following fields:
%
% out.TestResults = p-values based on collapsed contingency tables. 
%                   I-by-J matrix.
%                   The $(i,j)$-th entry of the TestResults matrix is the p-value
%                   of the Fisher exact test based on the collapsed $(i,j)$-th
%                   table. If the $(i,j)$-th entry of input matrix UserData is smaller
%                   or equal than the input parameter threshold, the test is not
%                   performed and the corresponding $(i,j)$-th entry of matrix
%                   TestResults is equal to Inf.
% out.RejectedBonf = Results of the tests based on Bonferrroni threshold. 
%                   Boolean matrix.
%                   The $(i,j)$-th entry of the RejectedBonf matrix is true 
%                   if the corresponding test based on the collapsed
%                   $(i,j)$-th table is significant. Bonferroni threshold is used.
% out.RejectedSidak = Results of the tests based on Sidak threshold. 
%                   Boolean matrix.
%                   The $(i,j)$-th entry of the RejectedSidak matrix is true if the
%                   corresponding test based on the collapsed
%                   $(i,j)$-th table is significant. Sidak threshold is used.
%
%
%
% More About:
%
%  $N$ = $I$-by-$J$-contingency table. The $(i,j)$-th element is equal to
%  $n_{ij}$,  $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$. The sum of the
%  elements of N is $n$ (the grand total). The sum of the elements of the
%  $i$-th row of the contingency table is denoted with $n_{i.}$ (n_idot in
%  the code). The sum of the elements of the $j$-th column of the
%  contingency table is denoted with $n_{.j}$ (n_dotj in the code).
%  $P$=$I$-by-$J$-table containing correspondence matrix
%  (proportions). The $(i,j)$-th element is equal to
%  $n_{ij}/n$, $i=1, 2, \ldots, I$ and $j=1, 2,
%  \ldots, J$.  The sum of the elements of $P$ is 1.
%  $P^*$=$I$-by-$J$-table containing correspondence matrix (proportions)
%  under the hypothesis of independence. The $(i,j)$-th element is equal to
%  $p_{ij}^*=p_{i.}p_{.j}$, $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$.  
%  The sum of the elements of $P^*$ is 1.
%
% See also: fishertest, barnardtest, CressieRead, rcontFS
%
%
% References:
%
% Arsenin, S., Riani M. (2018), Data mining large contingency tables
% standard approaches and a new method, submitted
%
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('SparseTableTest')">Link to the help function</a>
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%
%{
    % SparseTableTest with all default options.
    % For reproducibility
    rng default;  
    x1 = unidrnd(3,50,1);
    x2 = unidrnd(3,50,1);
    % Cross-tabulate x1 and x2.
    InputTable = crosstab(x1,x2);
    out = SparseTableTest(InputTable);
%}

%{
    % Cressie and Read test on collapsed contingency table.
    la=2/3;
    % Input is a matrix.
    % T = Contingency Table for Car Accident Type (rows) by
    % Accident Severity (columns)
    T=[2365 944 412; 249 585 276];
    out=SparseTableTest(T,'testname',la);
%}
 
%{
    %% Chi-squared test on collapsed contingency table.
    % Input is a data matrix and contingency table has to be built
    load smoke
    % X = original data matrix
    X=smoke.data;
    % Chi-squared test is used on collapsed 2-by-2 tables. 
    % Cells which have a frequency smaller or equal than 15 are ignored. 
    out=SparseTableTest(X,'datamatrix',true,'threshold',15,'testname',1);
    % show the output obtained
    RejectedBonf = out.RejectedBonf
    RejectedSidak = out.RejectedSidak
    TestResults = out.TestResults
%}
    
%% Beginning of code

% Check MATLAB version. If it is not smaller than 2013b than output is
% shown in table format
verMatlab=verLessThan('matlab','8.2.0');


% Check whether N is a contingency table or a n-by-p input dataset (in this
% last case the contigency table is built using the first tow columns of the
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
    N = crosstab(N(:,1),N(:,2));
    [I,J]=size(N);
else
    [I,J]=size(N);
end

if verMatlab ==0 && istable(N)
    N=table2array(N);
end

if any(any(N < 0))
    error('FSDA:SparseTableTest:WrongInput','Test expects counts that are nonnegative values');
end

if I < 2 || J < 2
    error('FSDA:SparseTableTest:WrongInput','Input contingencey table must at least be of size 2-by-2');
end

threshold = 2;
alpha     = 0.01;
testname  = 1;

if nargin > 1
    options=struct('testname',testname,'threshold',threshold,'alpha',alpha,'datamatrix',false);
    
    % UserOptions=varargin(1:2:length(varargin));
    if ~isempty(varargin)
        UserOptions=varargin(1:2:length(varargin));
        if ~isempty(UserOptions)
            % Check if number of supplied options is valid
            if length(varargin) ~= 2*length(UserOptions)
                error('FSDA:SparseTableTest:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
        
        threshold=options.threshold;
        alpha=options.alpha;
        testname=options.testname;
    end
end

if ischar(testname) && strcmp(testname,'fisher')
    lambda=100;
elseif ischar(testname) && strcmp(testname,'barnard')
    lambda=200;
else
    lambda=testname;
end

% Number of rows and number of columns of input data
[I,J] = size(N);

% The sample size
n = sum(sum(N));

% TestResults = matrix to store results. The i-th and j-th element of
% TestResults is equal to 1 if the collapsed test based on ij-th entry is
% significant
TestResults=zeros(I,J);

for j=1:J
    
    n_dotj = sum(N(:,j));
    
    for i = 1:I
        n_ij=N(i,j);
        
        % apply test just if n_ij is greater than threshold
        if (n_ij<=threshold)
            % UserData(i,j)=0;
            TestResults(i,j)=inf;
        else
            % Build the 2x2 Contingency Matrix
            
            n_idot=sum(N(i,:));
            
            % Collapsed 2-by-2 contingency table
            N_2by2=[n_ij n_idot-n_ij; n_dotj-n_ij n-n_idot-n_dotj+n_ij];
            
            % test on collapsed contingency table
            if lambda==100
                % Fisher exact test
                [~,pval] = fishertest(N_2by2);
            elseif lambda==200
                % Barnard test
            else
                % Power divergence family
                [~,pval]=CressieRead(N_2by2,'la',lambda);
            end
            
            % Store p-value of test
            TestResults(i,j)=pval;
        end
    end
end

% Select how many comparisons have to be made
nCompAfterT = sum(sum(N>threshold));

% Use Sidak or Bonferroni threshold
RejectedSidak = TestResults < (1-(1-alpha)^(1/nCompAfterT));
RejectedBonf  = TestResults < alpha/nCompAfterT;

out=struct;
out.RejectedBonf=RejectedBonf;
out.RejectedSidak=RejectedSidak;
out.TestResults=TestResults;

% An alternative approach is the one based on threshold which controls
% false discrovery rate
% k=(1:1:nCompAfterT)';
%pvalsor=sort(TestResults(:));
% Hoch=sum(pvalsor<alpha./(nCompAfterT+1-k));

end
%FScategory:MULT-Categorical
