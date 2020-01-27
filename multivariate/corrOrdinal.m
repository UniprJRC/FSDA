function out=corrOrdinal(N, varargin)
%corrOrdinal measures strength of association between two ordered categorical variables.
%
%<a href="matlab: docsearchFS('corrOrdinal')">Link to the help function</a>
%
% corrOrdinal computes Goodman-Kruskal's $\gamma$, $\tau_a$,
% $\tau_b$, $\tau_c$ of Kendall and $d_{y|x}$ of Somers.
% All these indexes measure the correlation among two ordered qualitative
% variables and go between -1 and 1. The sign of the coefficient indicates
% the direction of the relationship, and its absolute value indicates the
% strength, with larger absolute values indicating stronger
% relationships. Values close to an absolute value of 1 indicate a strong
% relationship between the two variables. Values close to 0 indicate little
% or no relationship. More in detail:
% $\gamma$ is a symmetric measure of association.
% Kendall's $\tau_a$ is a symmetric measure of association that does not
% take ties into account. Ties happen when both members of the data pair
% have the same value.
% Kendall's $\tau_b$ is a symmetric measure of association which takes ties
% into account. Even if $\tau_b$ ranges from -1 to 1, a value of -1 or
% +1 can be obtained only from square tables.
% $\tau_c$ (also called Stuart-Kendall $\tau_c$) is a symmetric measure of
% association which makes an adjustment for table size in addition to a
% correction for ties. Even if $\tau_c$ ranges from -1 to 1, a value of -1
% or +1 can be obtained only from square tables.
% Somers' $d$ is an asymmetric extension of $\tau_b$ in that it uses a
% correction only for pairs that are tied on the independent variable
% (which in this implementation it is assumed to be on the rows of the
% contingency table).
% Additional details about these indexes can be found in the "More About"
% section of this document.
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
%   NoStandardErrors:  Just indexes without standard errors and p-values.
%               Boolean.
%               if NoStandardErrors is true just the indexes are computed
%               without standard errors and p-values. That is no
%               inferential measure is given. The default value of
%               NoStandardErrors is false.
%               Example - 'NoStandardErrors',true
%               Data Types - Boolean
%
% dispresults : Display results on the screen. Boolean.
%               If dispresults is true (default) it is possible to see on the
%               screen all the summary results of the analysis.
%               Example - 'dispresults',false
%               Data Types - Boolean
%
%       Lr   :  Vector of row labels. Cell.
%               Cell containing the labels of the rows of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table. because in this case  Lr=N.Properties.RowNames;
%               Example - 'Lr',{'a' 'b' 'c'}
%               Data Types - cell array of strings
%
%       Lc   :  Vector of column labels. Cell.
%               Cell containing the labels of the columns of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table because in this case Lc=N.Properties.VariableNames;
%               Example - 'Lc',{'c1' c2' 'c3' 'c4'}
%               Data Types - cell array of strings
%
% datamatrix :  Data matrix or contingency table. Boolean.
%               If datamatrix is true the first input argument N is forced
%               to be interpreted as a data matrix, else if the input
%               argument is false N is treated as a contingency table. The
%               default value of datamatrix is false, that is the procedure
%               automatically considers N as a contingency table
%               Example - 'datamatrix',true
%               Data Types - logical
%
%  conflev:     Confidence levels to be used to
%               compute confidence intervals. Scalar.
%               The default value of conflev is 0.95, that
%               is 95 per cent confidence intervals
%               are computed for all the indexes (note that this option is
%               ignored if NoStandardErrors=true).
%               Example - 'conflev',0.99
%               Data Types - double
%
%  Output:
%
%         out:   structure which contains the following fields:
%
% 		out.N    =  $I$-by-$J$-array containing contingency table
%                   referred to active rows (i.e. referred to the rows which
%                   participated to the fit).
%                   The $(i,j)$-th element is equal to $n_{ij}$,
%                   $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$. The
%                   sum of the elements of out.N is $n$ (the grand
%                   total).
% 	out.Ntable   =  Same as out.N but in table format (with row and
%                   column names).
%                   This output is present just if your MATLAB
%                   version is not<2013b.
%      out.gam   =  1 x 4 vector which contains Goodman and Kruskall gamma index,
%                   standard error, test and p-value.
%      out.taua  =  1 x 4 vector which contains index $\tau_a$,
%                   standard error, test and p-value.
%      out.taub  =  1 x 4 vector which contains index $\tau_b$,
%                   standard error, test and p-value.
%       out.tauc =  1 x 4 vector which contains index $\tau_c$,
%                   standard error, test and p-value.
%       out.som  =  1 x 4 vector which contains Somers index $d_{y|x}$,
%                   standard error, test and p-value.
% out.TestInd    =  5-by-4 matrix containing index values (first column),
%                   standard errors (second column), zscores (third column),
%                   p-values (fourth column). Note that the
%                   standard errors in this matrix are computed assuming
%                   the null hypothesis of independence.
% out.TestIndtable  = 5-by-4 table containing index values (first column),
%                   standard errors (second column), zscores (third column),
%                   p-values (fourth column). Note that the
%                   standard errors in this table are computed assuming
%                   the null hypothesis of independence.
% out.ConfLim    =  5-by-4 matrix containing index values (first column),
%                   standard errors (second column), lower confidence limit
%                   (third column), upper confidence limit (fourth column).
%                   Note that the standard errors in this matrix are computed not
%                   assuming the null hypothesis of independence.
% out.ConfLimtable  = 5-by-4 table containing index values (first column),
%                   standard errors (second column), lower confidence limit
%                   (third column), upper confidence limit (fourth column).
%                   Note that the standard errors in this table are computed not
%                   assuming the null hypothesis of independence.
%
%
% More About:
%
%
% All these indexes are based on concordant and discordant pairs.
% A pair of observations is concordant if the subject who is higher on one
% variable also is higher on the other variable, and a pair of observations
% is discordant if the subject who is higher on one variable is lower on
% the other variable.
% Let $C$ be the total number of concordant pairs (concordances) and $D$
% the total number of discordant pairs (discordances) . If $C > D$ the
% variables have a positive association, but if $C < D$ then the variables
% have a negative association.
%
% In symbols, given an $I \times J$ contingency table the concordant pairs
% with cell $i,j$ are
% \[
% a_{ij} = \sum_{k<i} \sum_{l<j} n_{kl} + \sum_{k>i} \sum_{l>j}   n_{kl}
% \]
%
% the number of discordant pairs is
%
% \[
% b_{ij} = \sum_{k>i} \sum_{l<j} n_{kl} + \sum_{k<i} \sum_{l>j} n_{kl}
% \]
%
%
% Twice the number of concordances, $C$ is given by:
% \[
% 2 \times C = \sum_{i=1}^I \sum_{j=1}^J n_{ij} a_{ij}
% \]
% Twice the number of discordances, $D$ is given by:
% \[
% 2 \times D = \sum_{i=1}^I \sum_{j=1}^J n_{ij} b_{ij}
% \]
%
% Goodman-Kruskal's $\gamma$ statistic is equal to the ratio:
%
% \[
% \gamma= \frac{C-D}{C+D}
% \]
%
%
% $\tau_a$ is equal to concordant minus discordant pairs, divided by a
% factor which takes into account the total number of pairs.
%
% \[
%  \tau_a= \frac{C-D}{0.5 n(n-1)}
% \]
%
% $\tau_b$ is equal to concordant minus discordant
% pairs divided by a term representing the geometric mean between the
% number of pairs not tied on x and the number not tied on y.
% More precisely:
% \[
%  \tau_b= \frac{C-D}{\sqrt{ (0.5 n(n-1)-T_x)(0.5 n(n-1)-T_y)}}
% \]
% where $T_x= \sum_{i=1}^I 0.5 n_{i.}(n_{i.}-1)$ and
% $T_y=\sum_{j=1}^J 0.5 n_{.j}(n_{.j}-1)$
% Note that $\tau_b \leq \gamma$.
%
% $\tau_c$ is equal to concordant minus discordant pairs multiplied by a factor that adjusts
% for table size.
% \[
%  \tau_c= \frac{C-D}{ n^2(m-1)/(2m)}
% \]
% where $m= min(I,J)$;
%
% Somers' $d_{y|x}$ is an
% asymmetric extension of $\gamma$ that differs only in the inclusion of the
% number of pairs not tied on the independent variable. More precisely
%
% \[
%  d_{y|x} = \frac{C-D}{0.5 n(n-1)-T_x}
% \]
%
% Null hypothesis:
% corresponding index = 0. Alternative hypothesis (one-sided) index < 0 or
% index > 0.
%
%
%
% In order to compute confidence intervals and test hypotheses, this
% routine computes the standard error of the various indexes.
% Note that the expression of the standard errors which is used to compute
% the confidence intervals is different from the expression which is used
% to test the null hypothesis of no association (no relationship or independence)
% between the two variables.
%
%As concerns the Goodman-Kruskal's $\gamma$ index we have that:
% \[
%  var(\gamma) =   \frac{4}{(C + D)^4}  \sum_{i=1}^I \sum_{j=1}^J
%  n_{ij} (D a_{ij} - C b_{ij} )^2
% \]
% where
% \[
% d_{ij}=a_{ij}- b_{ij}
% \]
%
% The variance of $\gamma$  assuming the independence hypothesis is:
% \[
% var_0(\gamma) =\frac{1}{(C + D)^2} \left( \sum_{i=1}^I \sum_{j=1}^J
%  n_{ij} d_{ij}^2  -4(C-D)^2/n  \right)
% \]
%
% As concerns $\tau_a$ we have that:
% \[
% var(\tau_a)= \frac{2}{n(n-1)} \left\{ \frac{2(n-2)}{n(n-1)^2}  \sum_{i=1}^I \sum_{j=1}^J (d_{ij} - \overline d)^2 + 1 - \tau_a^2 \right\}
% \qquad \mbox{with $i,j$ such that $N(i,j)>0$}
% \]
% where
%
% \[
% \overline d = \sum_{i=1}^I \sum_{j=1}^J d_{ij} /n  \qquad \mbox{with $i,j$ such that $N(i,j)>0$}
% \]
%
% The variance of $\tau_a$  assuming the independence hypothesis is:
% \[
% var_0(\tau_a) =\frac{2 (2n+5)}{9n(n-1) }
% \]
%
%
% As concerns $\tau_b$ we have that:
% \[
% var(\tau_b)= \frac{n}{w^4} \left\{ n \sum_{i=1}^I \sum_{j=1}^J n_{ij} \tau_{ij}^2 - \left( \sum_{i=1}^I \sum_{j=1}^J  n_{ij}\tau_{ij}\right)^2  \right\}
% \]
% where
%
% \[
%  \tau_{ij} = 2n d_{ij} +2(C-D) n_{.j} w /n^3+2(C-D) (n_{i.}/n) \sqrt{
%  w_c/w_r} \qquad \mbox{and} \qquad w= \sqrt{w_rw_c}
% \]
%
% The variance of $\tau_b$  assuming the independence hypothesis is:
% \[
% var_0(\tau_b) =\frac{4}{w_r w_c}  \left\{ \sum_{i=1}^I \sum_{j=1}^J n_{ij} d_{ij} ^2 -4(C-D)^2/n \right\}
% \]
%
%
% As concerns Stuart's $\tau_c$ we have that:
% \[
% var(\tau_c)= \frac{4m^2}{(m-1)^2 n^4} \left\{ \sum_{i=1}^I \sum_{j=1}^J n_{ij} d_{ij} ^2 -4(C-D)^2/n \right\}
% \]
%
%
%
% The variance of $\tau_c$  assuming the independence hypothesis is:
% \[
% var_0(\tau_c) =var(\tau_c)
% \]
%
%
% As concerns $d_{y|x}$ we have that:
% \[
% var( d_{y|x})= \frac{4}{w_r^4} \left\{ \sum_{i=1}^I \sum_{j=1}^J n_{ij}
% (w_r d_{ij} -2(C-D) (n-n_{i.}) \right\}^2
% \]
% where
% \[
% w_r= n^2- \sum_{i=1}^I n_{i.}^2
% \]
%
% The variance of $d_{y|x}$  assuming the independence hypothesis is:
% \[
% var_0(d_{y|x}) = \frac{4}{w_r^2} \left\{ \sum_{i=1}^I \sum_{j=1}^J n_{ij} d_{ij} ^2 -4(C-D)^2/n \right\}
% \]
%
% From the theoretical point of view, Simon (1978) showed that all sample
% measures having the same numerator $(C-D)$ have the same efficacy and hence
% the same local power, for testing independence.
%
% See also crosstab, rcontFS, CressieRead, corr, corrNominal
%
% References:
%
% Agresti, A. (2002), "Categorical Data Analysis", John Wiley & Sons. [pp.
% 57-59]
% Agresti, A. (2010), "Analysis of Ordinal Categorical Data", Second Edition,
% Wiley, New York, pp. 194-195.
% Hollander, M, Wolfe, D.A., Chicken, E. (2014), "Nonparametric Statistical
% Methods", Third edition, Wiley,
% Goktas, A. and Oznur, I. (2011), A comparision of the most commonly used
% measures of association for doubly ordered square contingency tables via
% simulation, "Metodoloski zvezki", Vol. 8, pp. 17-37, [available at:
% www.stat-d.si/mz/mz8.1/goktas.pdf]
% Goodman, L.A. and Kruskal, W.H. (1954), Measures of association for
% cross classifications, "Journal of the American Statistical
% Association", Vol. 49, pp. 732-764.
% Goodman, L.A. and Kruskal, W.H. (1959), Measures of association for
% cross classifications II: Further Discussion and References,
% "Journal of the American Statistical Association", Vol. 54, pp. 123-163.
% Goodman, L.A. and Kruskal, W.H. (1963), Measures of association for
% cross classifications III: Approximate Sampling Theory,
% "Journal of the American Statistical Association", Vol. 58, pp. 310-364.
% Goodman, L.A. and Kruskal, W.H. (1972), Measures of association for
% cross classifications IV: Simplification of Asymptotic
% Variances, "Journal of the American Statistical Association", Vol. 67,
% pp. 415-421.
% Liebetrau, A.M. (1983), "Measures of Association", Sage University Papers
% Series on Quantitative Applications in the Social Sciences, 07-004,
% Newbury Park, CA: Sage. [pp. 49-56]
% SAS documentation (2009), See
% http://support.sas.com/documentation/cdl/en/statugfreq/63124/PDF/default/statugfreq.pdf,
% pp. 1738-1740.
% Morton, B.B. and Benedetti, J.K. (1977), Sampling Behavior of Tests for
% Correlation in Two-Way Contingency Tables, "Journal of the American
% Statistical Association", Vol. 72, pp. 309-315.
% Simon, G. (1978), Alternative analysis for the singly ordered contingency
% table, "Journal of the American Statistical Association", Vol. 69, pp.
% 971-976.
%
% Acknowledgements:
%
% This file was inspired by Trujillo-Ortiz, A. and R. Hernandez-Walls.
% gkgammatst: Goodman-Kruskal's gamma test. URL address
% http://www.mathworks.com/matlabcentral/fileexchange/42645-gkgammatst
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('corrOrdinal')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%

% Examples:

%{
    %%  corrOrdinal with all the default options.
    % Rows of N indicate the results of a written test with levels:
    % 'Sufficient' 'Good' Very good'
    % Columns of N indicate the results of an oral test with levels:
    % 'Sufficient' 'Good' Very good'
    N=[20    40    20;
       10    45    45;
        0     5    15];
    out=corrOrdinal(N);
    % Because the asymptotic 95 per cent confidence limits do not contain
    % zero, this indicates a strong positive association between the
    % written and the oral examination.
%}

%{
    %% Compare calculation of tau-b with that which comes from
    % Matlab function corr.
    % Starting from a contingency table, create the original data matrix to
    % te able to call corr.
    N=[20    23    20;
       21    25    22;
       18     18    19];
    n11=N(1,1); n12=N(1,2); n13=N(1,3);
    n21=N(2,1); n22=N(2,2); n23=N(2,3);
    n31=N(3,1); n32=N(3,2); n33=N(3,3);
    x11=[1*ones(n11,1) 1*ones(n11,1)];
    x12=[1*ones(n12,1) 2*ones(n12,1)];
    x13=[1*ones(n13,1) 3*ones(n13,1)];
    x21=[2*ones(n21,1) 1*ones(n21,1)];
    x22=[2*ones(n22,1) 2*ones(n22,1)];
    x23=[2*ones(n23,1) 3*ones(n23,1)];
    x31=[3*ones(n31,1) 1*ones(n31,1)];
    x32=[3*ones(n32,1) 2*ones(n32,1)];
    x33=[3*ones(n33,1) 3*ones(n33,1)];
    % X original data matrix
    X=[x11; x12; x13; x21; x22; x23; x31; x32; x33];
    % Find taub and pvalue of taub using MATLAB routine corr
    [RHO,pval]=corr(X,'type','Kendall');
    % Compute tau-b using FSDA corrOrdinal routine.
    out=corrOrdinal(X,'datamatrix',true,'dispresults',false);
    disp(['tau-b from MATLAB routine corr=' num2str(RHO(1,2))])
    disp(['tau-b from FSDA routine corrOrdinal=' num2str(out.taub(1))])
    % Remark the p-values are slightly different
    disp(['pvalue of H0:taub=0 from MATLAB routine corr=' num2str(pval(1,2))])
    disp(['pvalue of H0:taub=0 from FSDA routine corrOrdinal=' num2str(out.taub(4))])
%}

%{
    %  corrOrdinal with option conflev.
    N=[26 26 23 18  9;
       6  7  9 14 23];
    out=corrOrdinal(N,'conflev',0.999);
%}

%{
    % corrOrdinal with with option NoStandardErrors.
    N=[26 26 23 18  9;
       6  7  9 14 23];
    out=corrOrdinal(N,'NoStandardErrors',true);
%}

%{
    % Income and job satisfaction.
    % Relationship between the income (with levels '< 5000' '5000-25000' and
    % '>25000') and  job satisfaction (with levels 'Dissatisfied' 'Moderately satisfied'
    % and 'Very satisfied') for a sample of 300 persons
    % Input data is matlab table Ntable:
    N = [24 23 30;19 43 57;13 33 58];
    rownam={'Less_than_5000',  'Between_5000_and_25000' 'Greater_than_25000'};
    colnam= {'Dissatisfied' 'Moderately_satisfied' 'Very_satisfied'};
    if verLessThan('matlab','8.2.0') ==0
        Ntable=array2table(N,'RowNames',matlab.lang.makeValidName(rownam),'VariableNames',matlab.lang.makeValidName(colnam));
        %  Check relationship
        out=corrOrdinal(Ntable);
    else
        out=corrOrdinal(N);
    end
%}

%{
    % Input is the contingency table in matrix format, labels for rows and
    % columns are supplied.
    N=[20    40    20;
       10    45    45;
        0     5    15];
    % labels for rows and columns
    labels_rows= {'Sufficient' 'Good' 'Very_good'};
    labels_columns= {'Sufficient' 'Good' 'Very_good'};
    out=corrOrdinal(N,'Lr',labels_rows,'Lc',labels_columns,'dispresults',false);
    if verLessThan('matlab','8.2.0') ==0
        % out.Ntable uses labels for rows and columns which are supplied
        disp(out.Ntable)
    end
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
    [N,~,~,labels] =crosstab(N(:,1),N(:,2));
    [I,J]=size(N);
    % default labels for rows of contingency table
    Lr=labels(1:I,1);
    % default labels for columns of contingency table
    Lc=labels(1:J,2);
else
    [I,J]=size(N);
    % Size of N
    % default labels for rows of contingency table
    Lr=cellstr(strcat('r',num2str((1:I)')));
    % default labels for columns of contingency table
    Lc=cellstr(strcat('c',num2str((1:J)')));
end

dispresults=true;
NoStandardErrors=false;
conflev=0.95;

options=struct('Lr',{Lr},'Lc',{Lc},'datamatrix',false,...
    'dispresults',dispresults,'NoStandardErrors',NoStandardErrors,'conflev',conflev);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:CorrOrdinal:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    dispresults=options.dispresults;
    NoStandardErrors=options.NoStandardErrors;
    conflev=options.conflev;
    Lr  = options.Lr;
    Lc  = options.Lc;
    
end

% Extract labels for rows and columns
if verMatlab ==0 && istable(N)
    Ntable=N;
    N=table2array(N);
else
    if isempty(Lr)
        Lr=cellstr(num2str((1:I)'));
    else
        % Check that the length of Lr is equal to I
        if length(Lr)~=I
            error('Wrong length of row labels');
        end
    end
    
    if isempty(Lc)
        Lc=cellstr(num2str((1:J)'));
    else
        % Check that the length of Lc is equal to J
        if length(Lc)~=J
            error('Wrong length of column labels');
        end
    end
    if verMatlab ==0
        Ntable=array2table(N,'RowNames',matlab.lang.makeValidName(Lr),'VariableNames',matlab.lang.makeValidName(Lc));
    end
end


[I,J] = size(N);

n = sum(sum(N)); %sample size

if NoStandardErrors == true
    % Fast way to compute number of concordances and discordances if the
    % standard errors are not required
    C=0;
    D=0;
    for i=1:I
        for j=1:J
            if i<I && j<J
                xsel=N(i+1:I,j+1:J);
                C=C+N(i,j)*sum(xsel(:));
            end
            if j>1 && i<I
                xsel=N(i+1:I,1:j-1);
                D=D+N(i,j)*sum(xsel(:));
            end
        end
    end
    
else
    % con = nr \times nc matrix which will contain concordant pairs
    con = zeros(I,J);
    % dis = nr \times nc matrix which will contain discordace pairs
    dis = con;
    
    for i = 1:I
        for j = 1:J
            % xgyg = xgreater and y greater
            xgyg=N(i+1:I,j+1:J);
            % xsys = xsmaller and y smaller
            xsys=N(1:i-1,1:j-1);
            % store concordant pairs with i,j
            con(i,j) =   sum(xgyg(:))+ sum(xsys(:));
            % xgys = xgreater and y smaller
            xgys=N(i+1:I,1:j-1);
            % xsyg = xsmaller and y greater
            xsyg=N(1:i-1,j+1:J);
            % store discordant pairs with i,j
            dis(i,j) =   sum(xgys(:))+ sum(xsyg(:));
        end
    end
    
    C = sum(sum(N.*con))/2; %number of concordances
    D = sum(sum(N.*dis))/2; %number of discordances
    
end

out=struct;
out.N=N;
if verMatlab ==0
    out.Ntable=Ntable;
end

% totpairs= total number of pairs
totpairs=0.5*n*(n-1);
% nidot = row sums
nidot=sum(N,2);
% ndotj = columns sums
ndotj=sum(N,1);

% Gamma index
gam = (C - D)/(C + D);
% Tau-a index
taua = (C - D)/totpairs;
% Taub-b index
% Tx = number of pairs tied on X
Tx=sum(nidot.*(nidot-1)/2);
% Ty = number of pairs tied on Y
Ty=sum(ndotj.*(ndotj-1)/2);
taub=(C - D)/sqrt( (totpairs-Tx)*(totpairs-Ty) );
% Tau-c index
m=min(I,J);
tauc=m*2*(C-D)/( n^2 *(m-1));
% Somers index
wr=n^2-sum(nidot.^2);
som = 2*(C-D)/wr;

if NoStandardErrors
    segamH0=NaN;  zgam=NaN;  pvalgam=NaN;
    setauaH0=NaN; ztaua=NaN; pvaltaua=NaN;
    setaubH0=NaN; ztaub=NaN; pvaltaub=NaN;
    setaucH0=NaN; ztauc=NaN; pvaltauc=NaN;
    sesomH0=NaN; zsom=NaN; pvalsom=NaN;
    segam=NaN;
    setaua=NaN;
    setaub=NaN;
    setauc=NaN;
    sesom=NaN;
else
    
    % Compute required elements to find standard errors of the various indexes
    wc=n^2-sum(ndotj.^2);
    
    d=con-dis;
    sumdij2nij=sum( N(:) .* (d(:).^2) );
    
    %%  Goodman-Kruskal's gamma statistic
    % Find the standard error of the Goodman-Kruskal's gamma statistic
    s2gam=  (4/ (C+D)^4 )*sum(N(:).* ((D*con(:)-C*dis(:)).^2) );
    
    % The lines below contain an alternative formula to compute the
    % variance
    %  psi = 2*(D*con-C*dis)/(C+D)^2;
    % s2 = Goodman-Kruskal's gamma variance
    % s2gamCHK = sum(sum(N.*(psi.^2))) - sum(sum((N.*psi)))^2;
    
    % segam = Goodman-Kruskal's asymtotic standard error
    segam = sqrt(s2gam);
    % Standard error used to find the value of the test under the independence
    % hypothesis.
    segamH0=sqrt((1/(C+D)^2)*( sumdij2nij  - 4*(C-D)^2/n  ));
    
    zgam = gam/segamH0; % z-score
    pvalgam = 2*(1 - normcdf(abs(zgam))); %p-value (two-sided)
    
    %% tau-a statistic
    % Find standard error of tau-a
    Ci=con-dis;
    Ci=Ci(N~=0);
    Cbar=sum(Ci)/n;
    % setaua = standard error used to compute the confidence interval
    setaua= sqrt( 2/(n*(n-1)) * ((2*(n-2))/(n*(n-1)^2)* sum((Ci(:) - Cbar).^2) + 1 - taua^2));
    
    % setauaH0 = standard error used to find the value of the test under the independence
    % hypothesis.
    setauaH0 = sqrt(2*(2*n+5)/(9*n*(n-1)));
    ztaua = taua/setauaH0; % z-score
    pvaltaua = 2*(1 - normcdf(abs(ztaua))); %p-value (two-sided)
    
    %% tau-b statistic
    % For computationl purposes it is better to use relative frequencies
    % rather than absolute frequencies
    Pi=N/n; % matrix of relative frequencies
    pdiff=(con-dis)/n;
    Pdiff=2*(C-D)/n^2;
    delta1= sqrt(1 - sum((nidot/n).^2));
    delta2=sqrt(1 - sum((ndotj/n).^2));
    tauij=(2 * pdiff + Pdiff * repmat(ndotj/n,I,1) ) * delta2 * delta1 + ...
        (Pdiff * repmat(nidot/n,1,J) * delta2)/delta1;
    % setaub = standard errot used to compute the confidence interval
    setaub= sqrt(( (  sum(Pi(:) .* tauij(:).^2) - sum(Pi(:) .* tauij(:)).^2)/(delta1 * delta2)^4)/n);
    
    % The formula written in the help section (which uses the absolute
    % frequencies rather than the relative frequencies would be
    %     w=sqrt(wr*wc);
    %     tauijCHK=(2*n*d+2*(C-D)*repmat(ndotj/n,I,1) ) * w/n^2 + ...
    %         2*(C-D)* repmat(nidot/n,1,J) * sqrt(wc/wr);
    %     setaubCHK= sqrt((n/w^4)*  (   n*sum(N(:) .* tauijCHK(:).^2) -    sum(N(:) .* tauijCHK(:)).^2   ));
    
    % Formula below from the computational point of view is highly
    % inefficient
    %v=nidot*ones(1,J)+wc*ones(I,1)*ndotj;
    %setaubnew=sqrt((1/w^4)*( sum( N(:).*(  (2*w.*d(:)+taub*v(:)).^2 ) ) -n^3 *taub^2*((wr+wc)^2)));
    
    % Standard error used to find the value of the test under the independence
    % hypothesis.
    setaubH0=sqrt( (4/(wr*wc))*(sumdij2nij - (2*(C-D))^2/n ));
    ztaub = taub/setaubH0; % z-score
    pvaltaub = 2*(1 - normcdf(abs(ztaub))); %p-value (two-sided)
    
    %% (Stuart's) tau-c statistic
    % setauc = standard errot used to compute the confidence interval
    setauc= sqrt( 4*m^2/((m-1)^2*n^4) *  (sumdij2nij- 4*(C-D)^2/n ));
    
    % Standard error used to find the value of the test under the independence
    % hypothesis.
    setaucH0=setauc;
    
    ztauc = tauc/setaucH0; % z-score
    pvaltauc = 2*(1 - normcdf(abs(ztauc))); %p-value (two-sided)
    
    %% Somers' D statistic
    % Find standard error of Somers D stat
    nidotrep=repmat(nidot,1,J);
    
    % sesom = standard errot used to compute the confidence interval
    sesom=sqrt((4/(wr.^4)) *  sum( N(:) .* ((wr*(con(:)-dis(:))) - (2*(C-D))*(n-nidotrep(:)) ).^2));
    
    % sesomH0 = standard error used to test the independence hypothesis
    sesomH0= sqrt( 4/(wr^2) * (sumdij2nij   - (2*(C-D))^2/n ));
    
    zsom = som/sesomH0; % z-score
    pvalsom = 2*(1 - normcdf(abs(zsom))); %p-value (two-sided)
end

% Store results for Goodman-Kruskal's gamma statistic
out.gam=[gam segamH0 zgam pvalgam];

% Store results for tau-a statistic
out.taua=[taua setauaH0 ztaua pvaltaua];

% Store results for tau-b statistic
out.taub=[taub setaubH0 ztaub pvaltaub];

% Store results for tau-c statistic
out.tauc=[tauc setaucH0 ztauc pvaltauc];

% Store results for Somers  statistic
out.som=[som sesomH0 zsom pvalsom];

% Store results to test independence hypothesis
TestInd=[gam segamH0 zgam pvalgam;
    taua setauaH0 ztaua pvaltaua;
    taub setaubH0 ztaub pvaltaub;
    tauc setaucH0 ztauc pvaltauc;
    som sesomH0 zsom pvalsom];

out.TestInd=TestInd;
rownam={'gamma' 'taua' 'taub' 'tauc' 'dyx'};
colnamTestInd={'Coeff' 'se' 'zscore' 'pval'};
if verMatlab ==0
    TestIndtable=array2table(TestInd,'RowNames',rownam,'VariableNames',colnamTestInd);
    out.TestIndtable=TestIndtable;
end

% Store confidence intervals
talpha=-norminv((1-conflev)/2);
gamconflim=[gam segam gam-talpha*segam gam+talpha*segam];
tauaconflim=[taua setaua taua-talpha*setaua taua+talpha*setaua];
taubconflim=[taub setaub taub-talpha*setaub taub+talpha*setaub];
taucconflim=[tauc setauc tauc-talpha*setauc tauc+talpha*setauc];
somconflim=[som sesom som-talpha*sesom som+talpha*sesom];
ConfLim=[gamconflim; tauaconflim; taubconflim; taucconflim; somconflim];
out.ConfLim=ConfLim;
colnamConfInt={'Value' 'StandardError' 'ConflimL' 'ConflimU'};

if verMatlab ==0
    ConfLimtable=array2table(ConfLim,'RowNames',rownam,'VariableNames',colnamConfInt);
    out.ConfLimtable=ConfLimtable;
end

if dispresults == true
    if NoStandardErrors == false
        if verMatlab ==0
            
            % Test H_0
            % Test of independence
            disp('Test of H_0: independence between rows and columns')
            disp('The standard errors are computed under H_0')
            disp(TestIndtable);
            disp('-----------------------------------------')
            disp(['Indexes and ' num2str(conflev*100) '% confidence limits'])
            disp('The standard error are computed under H_1')
            disp(ConfLimtable);
        else
            % Test H_0
            % Test of independence
            disp('Test of H_0: independence between rows and columns')
            disp('The standard errors are computed under H_0')
            disp(colnamTestInd)
            disp(TestInd);
            disp('-----------------------------------------')
            disp(['Indexes and ' num2str(conflev*100) '% confidence limits'])
            disp('The standard error are computed under H_1')
            disp(colnamConfInt)
            disp(ConfLim);
        end
        
    else
        disp('-----------------------------------------')
        if verMatlab ==0
            disp(TestIndtable(:,1));
        else
            disp(TestInd(:,1));
        end
    end
end

end
%FScategory:MULT-Categorical