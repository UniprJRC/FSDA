function out=corrNominal(N, varargin)
%corrNominal measures strength of association between two unordered (nominal) categorical variables.
%
%<a href="matlab: docsearchFS('corrNominal')">Link to the help function</a>
%
% corrNominal computes $\chi2$, $\Phi$, Cramer's $V$, Goodman-Kruskal's
% $\lambda_{y|x}$, Goodman-Kruskal's  $\tau_{y|x}$, and Theil's $H_{y|x}$
% (uncertainty coefficient).
% All these indexes measure the association among two unordered qualitative
% variables.
% Additional details about these indexes can be found in the "More About"
% section or in the "Output section" of this document.
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
%       Lr   :  Vector of row labels. Cell.
%               Cell containing the labels of the rows of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table. because in this case  Lr=N.Properties.RowNames;
%               Example - 'Lr',{'a' 'b' 'c'}
%               Data Types - cell array of strings
%       Lc   :  Vector of column labels. Cell.
%               Cell containing the labels of the columns of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table because in this case Lc=N.Properties.VariableNames;
%               Example - 'Lc',{'c1' c2' 'c3' 'c4'}
%               Data Types - cell array of strings
% datamatrix :  data matrix or contingency table. Boolean.
%               If datamatrix is true the first input argument N is forced
%               to be interpreted as a data matrix, else if the input
%               argument is false N is treated as a contingency table. The
%               default value of datamatrix is false, that is the procedure
%               automatically considers N as a contingency table
%               Example - 'datamatrix',true
%               Data Types - logical
%  conflev:     confidence levels to be used to
%               compute confidence intervals. Scalar.
%               The default value of conflev is 0.95  that
%               is 95 per cent confidence intervals
%               are computed for all the indexes (note that this option is
%               ignored if NoStandardErrors=true
%               Example - 'conflev',0.99
%               Data Types - double
%
%  Output:
%
%         out:   structure which contains the following fields:
%
% 		out.N         =   $I$-by-$J$-array containing contingency table
%                         referred to active rows (i.e. referred to the rows which
%                         participated to the fit).
%                         The $(i,j)$-th element is equal to $n_{ij}$,
%                         $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$. The
%                         sum of the elements of out. is $n$ (the grand
%                         total).
% 		out.Ntable   =   same as out.N but in table format (with row and
%                         column names).
%        out.Chi2      = 1 x 2 vector which contains $\chi^2$ index,
%                        and p-value.
%       out.Phi    = 1 x 2 vector which contains index $\Phi$ index,
%                        and p-value.  Phi is a chi-square-based measure of
%                        association that involves dividing the chi-square
%                        statistic by the sample size and taking the square
%                        root of the result. More precisely
%                        \[
%                        \Phi= \sqrt{ \frac{\chi^2}{n} }
%                        \]
%                       This index lies in the interval $[0 \sqrt{\min[(I-1),(J-1)]}$.
%      out.CramerV     = 1 x 4 vector which contains Cramer's V index,
%                       standard error, z test, and p-value. Cramer'V index is index $\Phi$
%                       divided by its maximum. More precisely
%                        \[
%                       V= \sqrt{\frac{\Phi}{\min[(I-1),(J-1)]}}=\sqrt{\frac{\chi^2}{n \min[(I-1),(J-1)]}}
%                        \]
%                       The range of Cramer index is [0 1]. A Cramer's V in
%                       the range of [0, 0.3] is considered as weak,
%                       [0.3,0.7] as medium and > 0.7 as strong.
%                       In order to compute the confidence interval for
%                       this index we firt find a confidence interval for
%                       the non centrality parameter $\Delta$ of the
%                       $\chi^2$ distribution with $df=(I-1)(J-1)$ degrees of
%                       freedom. (see Smithson (2003) pp. 39 41) $[\Delta_L
%                       \Delta_U]$. A confidence interval for $\Delta$ is
%                       transformed into one for $V$ by the following
%                       transformation.
%                       \[
%                       V_L=\sqrt{\frac{\Delta_L+ df }{n \min[(I-1),(J-1)]}}
%                       \]
%                       and 
%                       \[
%                       V_U=\sqrt{\frac{\Delta_U+ df }{n \min[(I-1),(J-1)]}}
%                       \]
%       out.GKlambdayx  = 1 x 4 vector which contains index Goodman and Kruskall $\lambda_{y|x}$,
%                       standard error, z test, and p-value.
%                       \[
%                       \lambda_{y|x} = \sum_{i=1}^I \frac{r_i- r}{n-r}
%                       \]
%                       \[
%                       r_i =\max(n_{ij})
%                       \]
%                       \[
%                       r =\max(n_{.j})
%                       \]
%       out.tauyx      = 1 x 4 vector which contains tau index $\tau_{y|x}$,
%                       standard error, ztest and p-value.
%                       \[
%                       \tau_{y|x}= \frac{\sum_{i=1}^I \sum_{j=1}^J f_{ij}^2/f_{i.} -\sum_{j=1}^J f_{.j}^2 }{1-\sum_{j=1}^J f_{.j}^2 }
%                       \]
%          out.Hyx     = 1 x 4 vector which contains the uncertainty
%                       coefficient index (proposed by Theil) $H_{y|x}$,
%                       standard error, ztest and p-value.
%                       \[
%                       H_{y|x}= \frac{\sum_{i=1}^I \sum_{j=1}^J f_{ij} \log( f_{ij}/ (f_{i.}f_{.j}))}{\sum_{j=1}^J f_{.j} \log  f_{.j} }
%                       \]
% out.TestIndtable  = 4-by-4 table containing index values (first column),
%                   standard errors (second column), zscores (third column),
%                   p-values (fourth column).
% out.ConfLimtable  = 4-by-4 table containing index values (first column),
%                   standard errors (second column), lower confidence limit
%                   (third column), upper confidence limit (fourth column).
%
% More About:
%
%                       $\lambda_{y|x}$ is a measure of association that
%                       reflects the proportional reduction in error when
%                       values of the independent variable (variable in the
%                       rows of the contingency table) are used to predict
%                       values of the dependent variable (variable in the
%                       columns of the contingency table). The range of
%                       $\lambda_{y|x}$ is [0 1].  A value of 1
%                       means that the independent variable perfectly
%                       predicts the dependent variable. On the other hand,
%                       a value of 0 means that the independent variable is
%                       no help in predicting the dependent variable.
%                       More generally, let $V(y)$ a measure of variation
%                       for the marginal distribution $(f_{.1}=n_{.1}/n,
%                       ..., f_{.J}=n_{.J}/n)$ of the response $y$ and let
%                       $V(y|i)$ denote the same measure computed for the
%                       conditional distribution  $(f_{1|i}=n_{1|i}/n, ...,
%                       f_{J|i}=n_{J|i}/n)$ of $y$ at the $i$-th setting of
%                       the the explanatory variable $x$. A proportional
%                       reduction in variation measure has the form.
%                       \[
%                         \frac{V(y) - E[V(y|x)]}{V(y|x)}
%                       \]
%                       where  $E[V(y|x)]$ is the expectation of the
%                       conditional variation taken with respect to the
%                       distribution of $x$. When $x$ is a categorical
%                       variable having marginal distribution,
%                       $(f_{1.}, \ldots, f_{I.})$, 
%                       \[
%                       E[V(y|x)]= \sum_{i=1}^I (n_{i.}/n) V(y|i) =  \sum_{i=1}^I f_{i.} V(y|i)
%                       \]
%                       If we take as measure of variation $V(y)$ the Gini coefficient
%                       \[
%                       V(y)=1 -\sum_{j=1}^J f_{.j} \qquad V(y|i)=1 -\sum_{j=1}^J f_{j|i}
%                       \]
%                       we obtain index the proportionl reduction in
%                       variation $\tau_{y|x}$ of Goodman and Kruskall.
%                       \[
%                       \tau_{y|x}= \frac{\sum_{i=1}^I \sum_{j=1}^J f_{ij}^2/f_{i.} -\sum_{j=1}^J f_{.j}^2 }{1-\sum_{j=1}^J f_{.j}^2 }
%                       \]
%                       If, on the other hand, we take as measure of
%                       variation $V(y)$ the entropy index
%                       \[
%                       V(y)=-\sum_{j=1}^J f_{.j} \log f_{.j}  \qquad V(y|i) -\sum_{j=1}^J f_{j|i} \log f_{j|i}
%                       \]
%                       we obtain the index $H_{y|x}$, the uncertainty
%                       coefficient of Theil.
%                       \[
%                       H_{y|x}= \frac{\sum_{i=1}^I \sum_{j=1}^J f_{ij} \log( f_{ij}/ (f_{i.}f_{.j}))}{\sum_{j=1}^J f_{.j} \log  f_{.j} }
%                       \]
%                       The range of  $\tau_{y|x}$ and $H_{y|x}$ is [0 1]. 
%                       A large value of
%                       of the index represents a strong association, in
%                       the sense that we can guess $y$ much better when we
%                       know x than when we do not.
%                       In other words, $\tau_{y|x}=H_{y|x} =1$ is equivalent to no
%                       conditional variation in the sense that for each
%                       $i$, $n_{j|i}=1$. For example, a value of:
%                       $\tau_{y|x}=0.85$ indicates that knowledge of x
%                       reduces error in predicting values of y by 85 per
%                       cent (when the variation measure which is used is
%                       the Gini's index). 
%                       $H_{y|x}=0.85$ indicates that
%                       knowledge of x reduces error in predicting values
%                       of y by 85 per cent (when variation measure which
%                       is used is the entropy index)
%
% See also crosstab, rcontFS, CressieRead, corr, corrOrdinal
%
% References:
%
% Agresti, A. (2002). Categorical Data Analysis. John Wiley & Sons, pp. 23-26.
% Goodman, L. A. and Kruskal, W. H. (1959). Measures of association for
% cross classifications II: Further Discussion and References,
% Journal of the American Statistical Association, 54, pp. 123-163.
% Goodman, L. A. and Kruskal, W. H. (1963). Measures of association for
% cross classifications III: Approximate Sampling Theory,
% Journal of the American Statistical Association, 58, pp. 310-364.
% Goodman, L. A. and Kruskal, W. H. (1972). Measures of association for
% cross classifications IV: Simplification of Asymptotic
% Variances. Journal of the American Statistical Association, 67, pp.
% 415-421.
% Liebetrau, A. M. (1983). Measures of Association, Sage University Papers
% Series on Quantitative Applications in the Social Sciences, 07-004,
% Newbury Park, CA: Sage, pp. 49-56.
% Smithson, M.J. (2003) Confidence Intervals, Quantitative Applications in
% the Social Sciences Series, No. 140. Thousand Oaks, CA: Sage. pp. 39-41.
%
% Acknowledgements:
% In order to find the confidence interval for the non centrality parameter
% of the Chi squared distribution we use routine ncpci from the Effect Size Toolbox 
% Code by Harald Hentschke (University of Tübingen) and 
% Maik Stüttgen (University of Bochum)
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('corrNominal')">Link to the help function</a>
% Last modified 24-10-2017
%
%

% Examples:

%{
    %%  corrNominal with all the default options.
    % Rows of N indicate type of Bachelor degree:
    % 'Economics' 'Law' 'Literature'
    % Columns of N indicate employment type:
    % 'Private_firm' 'Public_firm' 'Freelance' 'Unemployed'
          N=[150	80	20	50
            80	250	30	140
            30	50	0	120];
    out=corrNominal(N);
%}

%{
    %% Example of option conflev.
    %  Use data from Goodman Kruskal (1954).
        N=[1768   807    189 47
           946   1387    746 53 
           115    438    288 16];
    out=corrNominal(N,'conflev',0.99);
%}


%{
    % corrNominal with option dispresults.
        N=[ 6 14 17 9;
            30 32 17 3];
    out=corrNominal(N,'dispresults',false);
%}

%{
    % Example which starts from the original data matrix.
    N=[26    26    23 18 9;
       6      7     9 14 23];
    % From the contingency table reconstruct the original data matrix.
    n11=N(1,1); n12=N(1,2); n13=N(1,3); n14=N(1,4); n15=N(1,5);
    n21=N(2,1); n22=N(2,2); n23=N(2,3); n24=N(2,4); n25=N(2,5);
    x11=[1*ones(n11,1) 1*ones(n11,1)];
    x12=[1*ones(n12,1) 2*ones(n12,1)];
    x13=[1*ones(n13,1) 3*ones(n13,1)];
    x14=[1*ones(n14,1) 4*ones(n14,1)];
    x15=[1*ones(n15,1) 5*ones(n15,1)];
    x21=[2*ones(n21,1) 1*ones(n21,1)];
    x22=[2*ones(n22,1) 2*ones(n22,1)];
    x23=[2*ones(n23,1) 3*ones(n23,1)];
    x24=[2*ones(n24,1) 4*ones(n24,1)];
    x25=[2*ones(n25,1) 5*ones(n25,1)];
    % X original data matrix
    X=[x11; x12; x13; x14; x15; x21; x22; x23; x24; x25];
    out=corrNominal(X,'datamatrix',true)
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
options.Lr=Lr;
options.Lc=Lc;

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:CorrNominal:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
if istable(N)
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

% sample size
n = sum(N(:));

% ndotj= column sums
ndotj=sum(N,1);

% nidot = row sums
nidot=sum(N,2);

% Matrix of theoretical frequencies
Ntheo=(nidot*ndotj/n);

% Chi2 index
Chi2=sum(((N(:)-Ntheo(:)).^2)./Ntheo(:));

% Phi index
Phi=sqrt(Chi2/n);

% Cramer index
CramerV=Phi/sqrt(min([I-1 J-1]));

% Goodman and Kruskall lambda
ndotjmax=max(ndotj);
nidotmax=max(nidot);

GKlambdayx=(sum(max(N,[],2))-ndotjmax)/(n-ndotjmax);

nidotmat=repmat(nidot,1,J);
ndotjmat=repmat(ndotj,I,1);

sumdotj2=sum(ndotj.^2);
boo=(N(:)~=0);

% tau index (proportional reduction in variation of Goodman and Kruskall),
% when variation is mesured using Gini's coefficient
tauyx= (  n*sum(N(:).^2./nidotmat(:))-sumdotj2)/(n^2-sumdotj2);

% see equation (4.4.3) GK JASA 1963 or p. 120 of GK Measures of association
% NumA=sum( ((N(:)./nidotmat(:)).*(1-sum( (ndotjmat(:)/n).^2 )) - (ndotjmat(:)/n).*(1-(1/n)*(N(:).^2 ./ nidotmat(:))) ).^2 .*N(:)/n)  ;
% NumB=sum((1/n)* N(:).^2 ./ nidotmat(:)) - sum( (ndotj(:)/n).^2 );
%
% vartauyx=(4/n)*(NumA- NumB^2)/  (1- sum( (ndotj/n).^2 ) ).^4;
% setauyx=sqrt(vartauyx);

% H index = proportional reduction in variation when variation is measured
% using the entropy
% H index (uncertainty coefficient) of Theil
Hyx= -sum(N(boo).*log(N(boo)./Ntheo(boo))  )/sum(ndotj.*log(ndotj/n));

Hx=-sum( (nidot/n).*log(nidot/n));
Hy=-sum( (ndotj/n).*log(ndotj/n));
Hxy=-sum( (N(:)/n).* log(N(:)/n) );
% Hyxchk=(Hx+Hy-Hxy)/Hy;


if NoStandardErrors
    seCramerV=NaN; zCramerV=NaN; pvalCramerV=NaN;
    seGKlambdayx=NaN;  zGKlambdayx=NaN;  pvalGKlambdayx=NaN;
    setauyx=NaN;     ztauyx=NaN;     pvaltauyx=NaN;
    seHyx=NaN;     zHyx=NaN;    pvalHyx=NaN;
else
    
    %  n * sum(N(i, )^2/sum.row[i])
    nerrunconditional=n^2- n*sum(N(:).^2 ./ nidotmat(:));
    
    nerrconditional= n^2- sum(ndotj.^2);
    errunconditional = nerrunconditional/(n^2);
    errconditional = nerrconditional/(n^2);
    f = errconditional * (errunconditional + 1) - 2 * errunconditional;
    % Unefficient implementation using loops
    %            vartauCR=0;
    %             for i=1:I
    %                 for j=1:J
    %                     vartauCR =vartauCR + N(i,j)*(-2*errunconditional*ndotj(j)/n +errconditional*(2*N(i,j)/nidot(i)-sum( (N(i,:)/nidot(i)).^2)) - f)^2/(n^2 * errconditional^4);
    %                 end
    %             end
    
    Ndivnidot=repmat(sum((N./nidotmat).^2,2),1,J);
    vartauyx=sum( N(:).* (-2*errunconditional*ndotjmat(:)/n +errconditional*(2*N(:)./nidotmat(:)-Ndivnidot(:) ) - f).^2 )/(n^2 * errconditional^4);
    setauyx=sqrt(vartauyx);
    
    % Find standard error for Cramer V
    % use external routine ncci to find confidence interval for non
    % centrality parameter of the chi2 distribution
    df=(I-1)*(J-1);
    k=min(I,J);
    
    % use external routine ncpci to find confidence interval
    ncpConfInt=ncpci(Chi2,'X2',df,'confLevel',conflev);
    ConfIntCramerV=sqrt((ncpConfInt+df)/(n*(k-1)));
    
    % Store confidence intervals
    talpha=-norminv((1-conflev)/2);
    seCramerV=(CramerV-ConfIntCramerV(1))/talpha;
    
    
    % The asymptotic variance of Gk index is given by
    %
    % \[
    % \lambda_{y|x} = \frac{n- \sum_{i=1}^I r_i }{(n-r)^3} \left( \sum_{i=1}^I
    % r_i +r -2 \sum_{i=1}^I (r_i|l_i=l) \right)
    % \]
    %
    % \[
    % r_i =\max(n_{ij})
    % \]
    % \[
    % r =\max(n_{.j})
    % \]
    %
    seqJ=1:J;
    seqI=1:I;
    
    % column index associated to maximal column frequency
    %rmax = max_j n_ij
    nijmax=max(N,[],2);
    % Note that Lcolmax is always a scalar (also in the case of ties)
    Lcolmax=min(seqJ(ndotj==ndotjmax));
    Lcol=zeros(I,1);
    for i=1:I
        testi= N(i, intersect(seqJ(N(i,:) ==ndotjmax), seqJ(N(i,:) == nidotmax)));
        if testi==n
            Lcol(i) = min(testi==n);
        elseif N(i, Lcolmax) == ndotjmax
            Lcol(i) = Lcolmax;
        else
            Lcol(i) = min(seqJ(N(i,:) == nijmax(i)));
        end
    end
    
    varGKlambdayx= (n - sum(nijmax)) * (sum(nijmax) + ndotjmax -2*(sum(nijmax(seqI(Lcol == Lcolmax)))))/(n-ndotjmax)^3;
    seGKlambdayx=sqrt(varGKlambdayx);
    
    % variance of uncertainty coefficient of Theil
    varHyx=sum(N(:).*(  Hy*log(N(:)./nidotmat(:)) +(Hx-Hxy)*log(ndotjmat(:)/n) ).^2)/(n^2*Hy^4);
    seHyx=sqrt(varHyx);
    
    % Compute zscores and p-values
    zCramerV = CramerV/seCramerV; % z-score
    pvalCramerV = 2*(1 - normcdf(abs(zCramerV))); %p-value (two-sided)
    zGKlambdayx = GKlambdayx/seGKlambdayx; % z-score
    pvalGKlambdayx = 2*(1 - normcdf(abs(zGKlambdayx))); %p-value (two-sided)
    ztauyx = tauyx/setauyx; % z-score
    pvaltauyx = 2*(1 - normcdf(abs(ztauyx))); %p-value (two-sided)
    zHyx = Hyx/seHyx; % z-score
    pvalHyx = 2*(1 - normcdf(abs(zHyx))); %p-value (two-sided)
    
end

% Store results in output structure out
out=struct;
out.N=N;
out.Ntable=Ntable;

out.Chi2=Chi2;
out.Phi=Phi;

out.CramerV=[CramerV seCramerV  zCramerV  pvalCramerV];
out.GKlambdayx=[GKlambdayx seGKlambdayx zGKlambdayx pvalGKlambdayx];
out.tauyx=[tauyx setauyx ztauyx pvaltauyx];
out.Hyx=[Hyx seHyx zHyx pvalHyx];

rownam={'CramerV' 'GKlambdayx' 'tauyx' 'Hyx'};

CramerVconflim=[CramerV seCramerV ConfIntCramerV];
GKlambdayxconflim=[GKlambdayx seGKlambdayx GKlambdayx-talpha*seGKlambdayx GKlambdayx+talpha*seGKlambdayx];
tauyxconflim=[tauyx setauyx tauyx-talpha*setauyx tauyx+talpha*setauyx];
Hyxconflim=[Hyx seHyx Hyx-talpha*seHyx Hyx+talpha*seHyx];


ConfLim=[CramerVconflim; GKlambdayxconflim; tauyxconflim;  Hyxconflim];
out.ConfLimtable=ConfLim;
colnam={'Value' 'StandardError' 'ConflimL' 'ConflimU'};
ConfLimtable=array2table(ConfLim,'RowNames',rownam,'VariableNames',colnam);
out.ConfLimtable=ConfLimtable;

% Store results to test independence hypothesis
TestInd=[CramerV seCramerV  zCramerV  pvalCramerV;
    GKlambdayx seGKlambdayx zGKlambdayx pvalGKlambdayx;
    tauyx setauyx ztauyx pvaltauyx;
    Hyx seHyx zHyx pvalHyx];

rownam={'CramerV' 'GKlambdayx' 'tauyx' 'Hyx'};
colnam={'Coeff' 'se' 'zscore' 'pval'};
TestIndtable=array2table(TestInd,'RowNames',rownam,'VariableNames',colnam);
out.TestIndtable=TestIndtable;

% if dispresults == true
%
%     if NoStandardErrors == false
%         disp('Goodmand and Kruskall lambda')
%         disp(GKlambdayx)
%         disp('tau index')
%         disp(tauyx)
%         disp('Uncertainty coefficient of Theil')
%         disp(Hyx)
%
%
%     else
%
%
%         % TODO
%
%     end
%
% else
%     % TODO
% end

if dispresults == true
    
    disp('Chi2 index')
    disp(Chi2)
    disp('Phi index')
    disp(Phi)
    disp('Cramer''s V ')
    disp(CramerV)
    
    if NoStandardErrors == false
        
        % Test H_0
        % Test of independence
        disp('Test of H_0: independence between rows and columns')
        disp(TestIndtable);
        disp('-----------------------------------------')
        disp(['Indexes and ' num2str(conflev*100) '% confidence limits'])
        disp(ConfLimtable);
    else
        disp('-----------------------------------------')
        disp(TestIndtable(:,1));
    end
end

end
