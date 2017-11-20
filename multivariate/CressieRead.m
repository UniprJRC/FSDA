function [PD , pval] = CressieRead(N,varargin)
%CressieRead computes the power divergence family
%
%<a href="matlab: docsearchFS('CressieRead')">Link to the help function</a>
%
%  Required input arguments:
%
%       N    :    Contingency table (default) or n-by-2 input dataset.
%                 Matrix or Table. Matrix or table which contains the input
%                 contingency table (say of size I-by-J) or the original
%                 data matrix. In this last case N=crosstab(N(:,1),N(:,2)).
%                 As default procedure assumes that the input is a
%                 contingency table.
%
%  Optional input arguments:
%
%       la   :  Parameter $\lambda$ of the family. Scalar. Scalar which 
%               contains the power in the Cressie-Read power divergence
%               statistics. The default value of la is 2/3.
%               If $\lambda=1$ we obtain Pearson's chi-squared statistic, 
%               see http://en.wikipedia.org/wiki/Chi-squared_test.
%               If $\lambda=0$ we obtain the Log-likelihood ratio (G, or G^2
%               test), see http://en.wikipedia.org/wiki/G-test.
%               If $\lambda=-0.5$ we obtain the Freeman-Tukey statistic, or
%               minumum Matusita distance (Hellinger distance).
%               If $\lambda=-1$ we obtain the modified log-likelihood ratio.
%               If $\lambda=-2$ we obtain the Neyman's statistic (or modified
%               chi squared estimate).
%               $\lambda=2/3$ is the value suggested by Cressie and Read (2004).
%               Example - 'la',0
%               Data Types - double
% datamatrix  : Data matrix or contingency table. Boolean. If datamatrix
%               is true the first input argument N is forced to be
%               interpreted as a data matrix, else if the input argument is
%               false N is treated as a contingency table. The default
%               value of datamatrix is false, that is the procedure
%               automatically considers N as a contingency table
%               Example - 'datamatrix',true
%               Data Types - logical
%
%  Output:
%
%       PD    : Cressie-Read power divergence test statistic. Scalar.
%               Scalar which measures the discrepancy/distance between
%               observed and expected frequencies, under the null hypothesis
%               that there is no difference in the row variable distribution
%               ('outcomes') between the columns ('treatments').
%
%     pval    : p-value of the test. Scalar. Value in the range [0,1] which 
%               represents the p-value of the test. The p-value is the
%               probability, under the null hypothesis, of observing a
%               value as extreme or more extreme of the power divergence
%               test statistic.
%
% More About:
%
%  $N$ = $I$-by-$J$-contingency table. The $(i,j)$-th element is equal to
%  $n_{ij}$,  $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$. The sum of the
%  elements of N is $n$ (the grand total).
%  $P$=$I$-by-$J$-table containing correspondence matrix
%  (proportions). The $(i,j)$-th element is equal to
%  $n_{ij}/n$, $i=1, 2, \ldots, I$ and $j=1, 2,
%  \ldots, J$.  The sum of the elements of $P$ is 1.
%  $P^*$=$I$-by-$J$-table containing correspondence matrix (proportions)
%  under the hypothesis of independence. The $(i,j)$-th element is equal to
%  $p_{ij}^*=p_{i.}p_{.j}$, $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$.
%  The sum of the elements of $P^*$ is 1.
%
% The power divergence family is defined:
%   
% \[
%  2n I^\lambda(P,P^*,\lambda) = \frac{2}{\lambda(\lambda+1)} 
%  n \sum_{i=1}^{I} \sum_{j=1}^{J} p_{ij} 
%  \left[ \left( \frac{p_{ij}}{p^*_{ij}} \right)^\lambda -1 \right]    
% \]
% where $\lambda$ is the family parameter.
% The term power divergence describes the fact that the statistic $2n
% I^\lambda(P,P^*,\lambda)$ measures the divergence of $P$ from $P^*$
% through a (weighted) sum of powers of the terms $\frac{p_{ij}}{p^*_{ij}}$
% for $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$.
% The reference distribution (independently of the value of $\lambda$) is
% $\chi^2$ with $(I-1) \times (J-1)$ degrees of freedom if (a) the hypothesis $H_0$ of
% independence is true; (b) the sample size $n$ is large; (c) the number of
% cells $I \times J$  is small relative to $n$ (and that all the expected cell
% frequencies $n^*_{ij}$ are large); (d) unknown parameters are estimated with BAN
% estimates; and (e) the models satisfy the regularity conditions of Birch
% (1964) (Cressie and Read, 1984; p. 63).
% If some of the expected frequencies are very small while others are
% greater than 1, there are no firm recommendations regarding the best
% approximation to use for the critical value. In some of these cases the
% normal approximation may be appropriate, however the convergence of
% this normal approximation is slow (Cressie and Read, 1984; p. 63).
% The test we have just performed is intrinsically non-directional, that is 
% the significant value says nothing at all about the
% particular direction of the difference.
% Suggestions that X^2 approximates a chi-squared random variable more
% closely than does G2 (for various multinomial and contingency table models)
% have been made in enumeration and simulation studies.
% Cressie and Read p. 85: extreme values of $\lambda$ are most useful for
% detecting deviations from an expected cell frequency in a single cell
% (i.e., bump or dip alternatives). Use of the maximum cell frequency to
% test for "spikes" is discussed by Levin (1983); he also derives the
% appropriate distributional results for hypothesis testing.
%
% See also CorAna, crosstab, rcontFS, corrNominal, corrOrdinal
%
% References:
%
% Rudas, T. (1986), A Monte Carlo Comparision of Small Sample Behaviour of
% The Pearson, the Likelihood Ratio and the Cressie-Read Statistics,
% Journal Statistcal Computation and Simulation, vol 24, pp 107-120.
% Read, TRC and Cressie, NAC (1988), Goodness of Fit Statistics for
% Discrete Multivariate Data. Springer Verlag.
% Ewens, WJ and Grant, GR (2001), Statistical Methods in Bioinformatics,
% Springer Verlag.
%
% Copyright 2008-2017.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('CressieRead')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Compute Cressie and Read test using lambda=2/3 (input is a matrix).
    % Input is a matrix.
    % T = Contingency Table for Car Accident Type (rows) by
    % Accident Severity (columns)
    T=[2365 944 412; 249 585 276];
    [PD,pval]=CressieRead(T);
%}

%{
    % Compute Cressie and Read test using lambda=2/3 (input is a table).
    % T = Contingency Table for Car Accident Type (rows) by
    % Accident Severity (columns)
    % See Table 3.1 of Cressie and Read (1984)
    T=[2365 944 412; 249 585 276];
    AccidentType={'Rollover' 'NotRollover'};
    AccidentSeverity={'NotSevere' 'ModeratelySevere' 'Severe'};
    if verLessThan('matlab','8.2.0') ==0
        Ttable=array2table(T,'RowNames',AccidentType,'VariableNames',AccidentSeverity);
        [PD,pval]=CressieRead(Ttable);
    else
        [PD,pval]=CressieRead(T);
    end
%}

%{
    %% Compute Cressie and Read statistic for a series of values of lambda.
    % T = Contingency Table for Car Accident Type (rows) by
    % Accident Severity (columns)
    T=[2365 944 412; 249 585 276];
    la=-2:0.1:3;
    PD=zeros(length(la),1);
    for i=1:length(la)
        PD(i)=CressieRead(T,'la',la(i));
    end
    plot(la,PD) 
    xlabel('\lambda')
    ylabel('Cressie and Read test statistic')
%}

%{
    % Compute Cressie and Read test using a set of values of lambda (input is a table).
    % T = Contingency Table for Dumping Severity by Operation
    % See Table 3.4 of Cressie and Read (1984)
    T=[61 28 7; 68 23 13; 58 40 12; 53 38 16];
    Operation={'A1' 'A2' 'A3' 'A4'};
    DumpingSeverity={'None' 'Slight' 'Moderate'};
    if verLessThan('matlab','8.2.0') ==0
        T=array2table(T,'RowNames',Operation,'VariableNames',DumpingSeverity);
    end
    la=[-5 -2 -1 -0.5 0 0.5 2/3 1 2 5];
    PD=[la' zeros(length(la),2)];
    for i=1:length(la)
        [PD(i,2),PD(i,3)]=CressieRead(T,'la',la(i));
    end
    disp(PD)
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
    N =crosstab(N(:,1),N(:,2));
    [I,J]=size(N); 
else
    [I,J]=size(N);
end

if verMatlab ==0 && istable(N) 
    N=table2array(N);
end

if any(any(N < 0))
    error('Test expects counts that are nonnegative values');
end

if I < 2 || J < 2
    error('Matrix of observation must at least be of size 2-by-2');
end

% default value of lambda
la=2/3;
options=struct('la',la,'datamatrix',false);

% UserOptions=varargin(1:2:length(varargin));
if ~isempty(varargin)
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:CressieRead:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    la=options.la;
end

n=sum(sum(N));
% Tstar = matrix containing expected frequencies under the hypothesis of
% independence
Tstar=sum(N,2)*sum(N,1)/n;
sel=N(:)>0;

tol=1e-13;
if abs(la)<tol
    PD=2*sum(N(sel).*log((N(sel)./Tstar(sel))));
elseif abs(la+1)<tol
    PD=2*sum(Tstar(sel).*log((Tstar(sel)./N(sel))));
else
    PD=(2/(la*(la+1)))*sum(N(sel).*((N(sel)./Tstar(sel)).^la-1));
end

df=(I-1)*(J-1);  % degrees of freedom
% p-value of the test
pval = 1-chi2cdf(PD,df);
end
%FScategory:MULT-Categorical