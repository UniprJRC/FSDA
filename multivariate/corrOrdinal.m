function out=corrOrdinal(N, varargin)
%corrOrdinal measures strength of association between two ordered qualitative variables.
%
%<a href="matlab: docsearchFS('corrOrdinal')">Link to the help function</a>
%
% corrOrdinal computes Goodman-Kruskal's $\gamma$, $\tau_a$,
% $\tau_b$, $\tau_c$, and $d_{y|x}$ of Somers.
% All these indexes measure the correlation among two ordered qualitative
% variables and go between -1 and 1. The sign of the coefficient indicates
% the direction of the relationship, and its absolute value indicates the
% strength, with larger absolute values indicating stronger relationships.
% Values close to an absolute value of 1 indicate a strong relationship
% between the two variables. Values close to 0 indicate little or no
% relationship. More in detail:
% $\gamma$ is a symmetric measure of association.
% $\tau_a$ is a symmetric measure of association that does not take ties
% into account. Ties happen when both members of the data pair have the
% same value.
% Kendall's $\tau_b$ is a symmetric measure of association which takes ties
% into account. Even if $\tau_b$ ranges from -1 to 1, a value of -1 or
% +1 can be obtained only from square tables.
% $\tau_c$ (also called Stuart-Kendall Tau-c), differs from $\tau_b$ as in
% being more suitable for rectangular tables than for square tables.
% $\tau_c$ is a symmetric measure of association which makes an
% adjustment for table size in addition to a correction for ties. Even if
% $\tau_c$ ranges from -1 to 1, a value of -1 or +1 can be obtained only
% from square tables.
% Somers' $d$ is an asymmetric extension of $\tau_b$ in that it uses a
% correction only for pairs that are tied on the independent variable
% (which in this implementation it is assumed to be on the rows of the
% contingency table).
% Additional details about these indexes can be found in the More About
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
%   NoStandardErrors:  Just index without standard errors and p-values.
%               Boolean.
%               if NoStandardErrors is true just the indexes are computed
%               without standard errors and p-values. That is no
%               inferential measure is given. The default value of
%               NoStandardErrors is false.
%                 Example - 'NoStandardErrors',true
%                 Data Types - Boolean
%  dispresults :  display results on the screen. Boolean.
%                 If dispresults is true (default) it is possible to see on the
%                 screen all the summary results of the analysis.
%                 Example - 'dispresults',false
%                 Data Types - Boolean
%
%  Output:
%
%         out:   structure which contains the following fields:
%
%                out.taua= 1 x 4 vector which contains index $\tau-a$,
%                standard error, test and p-value.
%                out.taub= 1 x 4 vector which contains index $\tau-b$,
%                standard error, test and p-value.
%                out.tauc= 1 x 4 vector which contains index $\tau-c$,
%                standard error, test and p-value.
%                out.dyx= 1 x 4 vector which contains Somers index $d_y|x$,
%                standard error, test and p-value.
%                out.gk= 1 x 4 vector which contains Goodman and Kruskall tau index,
%                standard error, test and p-value.
%
% More About:
%
% Null hypothesis:
% corresponding index = 0. Alternative hypothesis (one-sided) index < 0 or
% index > 0.
%
% All these indexes are based on concordant and discordant pairs.
% A pair of observations is concordant if the subject who is higher on one
% variable also is higher on the other variable, and a pair of observations
% is discordant if the subject who is higher on one variable is lower on
% the other variable.
% If C > D the variables have a positive association, but if C < D then the
% variables have a negative association. C and D are, respectivelly, the
% total number of concordances and discordances.
%
%
% In symbols, given an $I \times J$ contingency table
% \[
%            a_{ij} = \sum_{k<i} \sum_{l<j} n_{kl} + \sum_{k>i} \sum_{l>j}   n_{kl}
% \]
%
% \[
%            b_{ij} = \sum_{k>i} \sum_{l<j} n_{kl} + \sum_{k<i} \sum_{l>j} n_{kl}
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
% $\tau_a$ is equal to concordant minus discordant pairs, divided by a factor to
% account for total number of pairs (sample size).
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
% for table size).
% \[
%  \tau_c= \frac{C-D}{ n^2(m-1)/(2m)}
% \]
% where $m= min(I,J)$;
%
% Somers' $d_{y|x}$. is an
% asymmetric extension of $\gamma$ that differs only in the inclusion of the
% number of pairs not tied on the independent variable. More precisely
%
% \[
%  d_{y|x} = \frac{C-D}{(0.5 n(n-1)-T_x)}
% \]
%
%
% See
% http://support.sas.com/documentation/cdl/en/statugfreq/63124/PDF/default/statugfreq.pdf,
% pp. 1739 for the estimation of the asymptotic variance.
%
% Standard error of $\tau-a$ is computed as follows:
% \[
% s.e. \tau-a =\frac{1}{3} \sqrt{n(n-1)(2n+5)/2}
% \]
%
% In order to compute the confidence intervals and test hypothesis this
% routine computes the standard error of the various indexes.
% Note that the expression of the standard errors which is used to compute
% the confidence intervals is different from the expression which is used
% to test the null hypothesis of no association (no relationship or independence)
% between the two variables.
%
% The asymtotic Goodman-Kruskal's gamma variance is calculated as,
%
%   varg = 4/(C + D)^2*(Sum_i Sum_j n_ij*((A_ij - B_ij)^2 - (C - D)^2/n)

% Example. We take the example given in the Lecture 14 of Simon Jackman's
% Political Science 151B course, Stanford University. Where there is the
% interest to found and make a statistical inference on the relationship
% between the job satisfaction and income for a sample of 104 black
% americans [URL address http://jackman.stanford.edu/classes/151B/06/
% Lecture14.pdf]. As you can see, the categorical (qualitative) independent
% variables (Job satisfaction and Income) are ordered. We use an one tailed
% test and alpha-value of 0.05. Data are:
%
%        ------------------------------------------------------
%                                  Job satisfaction
%        ------------------------------------------------------
%                                     Moderately     Very
%          Income        Dissatisfied  satisfied   satisfied
%        ------------------------------------------------------
%          < 5000             6           13           3
%         5000-25000          9           37          12
%          > 25000            3           13           8
%        ------------------------------------------------------
%
% Input data:
% x = [6 13 3;9 37 12;3 13 8];
%
% x=[20 40 20; 10 45 45; 0 5 15];
% x=[20 40 30; 10 45 55; 15 25 35];

% Calling on Matlab the function:
%                gkgammatst(x,0.05,1)
%
% Answer is:
% ---------------------------------------------------------------------------------------
% Sample size: 104
% Contingency table: 3 x 3
% Goodman-Kruskal's gamma statistic: 0.2873
% Goodman-Kruskal's asymtotic standard error: 0.1506
% z-value: 1.9081
% Sample size: 104
% P-value associated to the Goodman-Kruskal's gamma statistic: 0.0282
% In a one-sided test:
% With a given significance = 0.050
% There is a significant positive relationship between the ordered qualitative variables.
% ---------------------------------------------------------------------------------------
%
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A. and R. Hernandez-Walls. (2013). gkgammatst:
%    Goodman-Kruskal's gamma test. URL address
%    http://www.mathworks.com/matlabcentral/fileexchange/42645-gkgammatst
%
% See also crosstab, rcontFS, CressieRead
%
% References:
%
% Goktas, A. and Oznur, I. (2011). A comparision of the most commonly used
%              measures of association for doubly ordered square
%              contingency tables via simulation. Metodoloski zvezki 8(1):
%              17-37. (URL address: www.stat-d.si/mz/mz8.1/goktas.pdf)
% Goodman, L. A. and Kruskal, W. H. (1954). Measures of association for
%              cross classifications. Journal of the American Statistical
%              Association, 49:732-764.
% Goodman, L. A. and Kruskal, W. H. (1959). Measures of association for
%              cross classifications II:Further Discussion and References.
%              Journal of the American Statistical Association, 54:123-163.
% Goodman, L. A. and Kruskal, W. H. (1963). Measures of association for
%              cross classifications III: Approximate Sampling Theory,
%              Journal of the American Statistical Association, 58:310-364.
% Goodman, L. A. and Kruskal, W. H. (1972). Measures of association for
%              cross classifications IV: Simplification of Asymptotic
%              Variances. Journal of the American Statistical Association,
%              67:415-421.
%Agresti, A. (2002) Categorical Data Analysis. John Wiley & Sons, pp. 57–59.
% Hollander, M, Wolfe, D. A., Chicken, E. (2014) Nonparametric Statistical
% Methods, Third edition, Wiley,
%Liebetrau, A. M. (1983) Measures of Association, Sage University Papers
%Series on Quantitative Applications in the Social Sciences, 07-004,
%Newbury Park, CA: Sage, pp. 49-56

%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('corrOrdinal')">Link to the help function</a>
% Last modified 05-06-2016
%
%

% Examples:

%{
    %%  corrOrdinal with all the default options.
    % Rows of N indicate the results of an oral test with levels:
    % 'Sufficient' 'Good' Very good'
    % Columns of N indicate the results of an oral test with levels:
    % 'Sufficient' 'Good' Very good'
        N=[20    40    20;
            10    45    45;
             0     5    15];
    out=corrOrdinal(N);
    % Because the asymptotic 95 per cent confidence limits do not contain
    % zero, this indicates a strong positive association between the
    % written and oral examination.
%}

%{
    %% test
    x11=[1*ones(20,1) 1*ones(20,1)];
    x12=[1*ones(40,1) 2*ones(40,1)];
    x13=[1*ones(20,1) 3*ones(20,1)];
    x21=[2*ones(10,1) 1*ones(10,1)];
    x22=[2*ones(45,1) 2*ones(45,1)];
    x23=[2*ones(45,1) 3*ones(45,1)];
    x31=[3*ones(0,1) 1*ones(0,1)];
    x32=[3*ones(5,1) 2*ones(5,1)];
    x33=[3*ones(15,1) 3*ones(15,1)];
    X=[x11; x12; x13; x21; x22; x23; x31; x32; x33];

    [RHO,pval]=corr(X,'type','Kendall');

    N=[26 26 23 18  9;
     6  7  9 14 23];

%}

%% Beginning of code

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

options=struct('datamatrix',false,...
    'dispresults',dispresults,'NoStandardErrors',NoStandardErrors,'conflev',conflev);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:CorAna:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
end

% Extract labels for rows and columns
if istable(N)
    Lc=N.Properties.VariableNames;
    Lr=N.Properties.RowNames;
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
Ntable=array2table(N,'RowNames',matlab.lang.makeValidName(Lr),'VariableNames',matlab.lang.makeValidName(Lc));
end


[I,J] = size(N);

n = sum(sum(N)); %sample size

if NoStandardErrors == true
    % Fast way to compute number of concordances and discordances if the
    % standard errors are nor required
    C=0;
    D=0;
    for i=1:I
        for j=1:J
            if i<I && j<J
                xsel=N(i+1:I,j+1:J);
                C=C+sum(xsel(:));
            end
            if j>1 && i<I
                xsel=N(i+1:I,1:j-1);
                D=D+sum(xsel(:));
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
out.Ntable=Ntable;

% Compute required elements to find standard errors

% totpairs= total number of pairs
totpairs=0.5*n*(n-1);

% nidot = row sums
nidot=sum(N,2);
% ndotj = columns sums
ndotj=sum(N,1);

wr=n^2-sum(nidot.^2);
wc=n^2-sum(ndotj.^2);
%w=sqrt(wr*wc);
%v=nidot*ones(1,J)+wc*ones(I,1)*ndotj;
d=con-dis;
sumdij=sum( N(:) .* (d(:).^2) );

%%  Goodman-Kruskal's gamma statistic
gam = (C - D)/(C + D);
% Find the standard error of the Goodman-Kruskal's gamma statistic
psi = 2*(D*con-C*dis)/(C+D)^2;
% s2 = Goodman-Kruskal's gamma variance
s2gam = sum(sum(N.*(psi.^2))) - sum(sum((N.*psi)))^2;
% segam = Goodman-Kruskal's asymtotic standard error
segam = sqrt(s2gam);
% Standard error used to find the value of the test under the independence
% hypothesis.
segamH0=sqrt((1/(C+D)^2)*( sumdij  - (2*(C-D))^2/n  ));

zgam = gam/segamH0; %inference via z-score
pvalgam = 2*(1 - normcdf(abs(zgam))); %p-value (two-sided)
% Store results for Goodman-Kruskal's gamma statistic
out.gam=[gam segamH0 zgam pvalgam];

%% tau-a statistic
taua = (C - D)/totpairs;
% Find standard error of tau-a
Ci=con-dis;
Ci=Ci(N~=0);
Cbar=sum(Ci)/n;
% setaua = standard error used to compute the confidence interval
setaua= sqrt( 2/(n*(n-1)) * ((2*(n-2))/(n*(n-1)^2)* sum((Ci(:) - Cbar).^2) + 1 - taua^2));

% setauaH0 = standard error used to find the value of the test under the independence
% hypothesis.
setauaH0 = sqrt(2*(2*n+5)/(9*n*(n-1)));
ztaua = taua/setauaH0; %inference via z-score
pvaltaua = 2*(1 - normcdf(abs(ztaua))); %p-value (two-sided)
% Store results for tau-a statistic
out.taua=[taua setauaH0 ztaua pvaltaua];

%% tau-b statistic
% Tx = number of pairs tied on X
Tx=sum(nidot.*(nidot-1)/2);
% Ty = number of pairs tied on Y
Ty=sum(ndotj.*(ndotj-1)/2);

taub=(C - D)/sqrt( (totpairs-Tx)*(totpairs-Ty) );

Pi=N/n; % matrix of relative frequencies
pdiff=(con-dis)/n;
Pdiff=2*(C-D)/n^2;
delta1= sqrt(1 - sum((nidot/n).^2));
delta2=sqrt(1 - sum((ndotj/n).^2));
tauphi=(2 * pdiff + Pdiff * repmat(ndotj/n,I,1) ) * delta2 * delta1 + ...
    (Pdiff * repmat(nidot/n,1,J) * delta2)/delta1;
% setaub = standard errot used to compute the confidence interval
setaub= sqrt(((sum(Pi(:) .* tauphi(:).^2) - sum(Pi(:) .* tauphi(:)).^2)/(delta1 * delta2)^4)/n);
% Notation from Agresti, A. (2002) Categorical Data Analysis. John Wiley & Sons, pp. 57-59.

% Formula to check
%setaub=sqrt((1/w^4)*( sum( N(:).*(  (2*w.*d(:)+taub*v(:)).^2 ) ) -n^3 *taub^2*((wr+wc)^2)));

% Standard error used to find the value of the test under the independence
% hypothesis.
setaubH0=sqrt( (4/(wr*wc))*(sumdij - (2*(C-D))^2/n ));
ztaub = taub/setaubH0; %inference via z-score
pvaltaub = 2*(1 - normcdf(abs(ztaub))); %p-value (two-sided)
% Store results for tau-b statistic
out.taub=[taub setaubH0 ztaub pvaltaub];

%% (Stuart's) tau-c statistic
m=min(I,J);
tauc=m*2*(C-D)/( n^2 *(m-1));

% setauc = standard errot used to compute the confidence interval
setauc= sqrt( 4*m^2/((m-1)^2*n^4) *  (sumdij- (2*(C-D))^2/n ));

% Standard error used to find the value of the test under the independence
% hypothesis.
setaucH0=setauc;

ztauc = tauc/setaucH0; %inference via z-score
pvaltauc = 2*(1 - normcdf(abs(ztauc))); %p-value (two-sided)
% Store results for tau-c statistic
out.tauc=[tauc setaucH0 ztauc pvaltauc];

%% Somers' D statistic
som = 2*(C-D)/wr;
% Find standard error of Somers D stat
nidotrep=repmat(nidot,1,J);

% sesom = standard errot used to compute the confidence interval
sesom=sqrt((4/(wr.^4)) *  sum( N(:) .* ((wr*(con(:)-dis(:))) - (2*(C-D))*(n-nidotrep(:)) ).^2));

% sesomH0 = standard error used to test the independence hypothesis
sesomH0= sqrt( 4/(wr^2) * (sumdij   - (2*(C-D))^2/n ));

zsom = som/sesomH0; %inference via z-score
pvalsom = 2*(1 - normcdf(abs(zsom))); %p-value (two-sided)
% Store results for Somers  statistic
out.som=[som sesomH0 zsom pvalsom];

if dispresults == true
    
    % Test H_0
    % Test of independence
    disp('Test of H_0: independence between rows and columns')
    disp('The standard errors are computed under H_0')
    TestInd=[gam segamH0 zgam pvalgam;
        taua setauaH0 ztaua pvaltaua;
        taub setaubH0 ztaub pvaltaub;
        tauc setaucH0 ztauc pvaltauc;
        som sesomH0 zsom pvalsom];
    rownam={'gamma' 'taua' 'taub' 'tauc' 'dyx'};
    colnam={'Coeff' 'se' 'zscore' 'pval'};
    TestIndtable=array2table(TestInd,'RowNames',rownam,'VariableNames',colnam);
    disp(TestIndtable);
    
    disp(['Indexes and ' num2str(conflev*100) '% confidence limits'])
    disp('The standard error are computed under H_1')
    talpha=-norminv((1-conflev)/2);
    gamconflim=[gam segam gam-talpha*segam gam+talpha*segam];
    tauaconflim=[taua setaua taua-talpha*setaua taua+talpha*setaua];
    taubconflim=[taub setaub taub-talpha*setaub taub+talpha*setaub];
    taucconflim=[tauc setauc tauc-talpha*setauc tauc+talpha*setauc];
    somconflim=[som sesom som-talpha*sesom som+talpha*sesom];
    conflim=[gamconflim; tauaconflim; taubconflim; taucconflim; somconflim];
    out.conflim=conflim;
    colnam={'Value' 'StandardError' 'ConflimL' 'ConflimU'};
    ConfLimtable=array2table(conflim,'RowNames',rownam,'VariableNames',colnam);
    disp(ConfLimtable);
end

end
