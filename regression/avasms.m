function [BestSol,corMatrix]=avasms(y,X,varargin)
%avasms computes avas using a series of alternative options
%
%<a href="matlab: docsearchFS('avasms')">Link to the help page for this function</a>
%
%   This function applies avas with a series of options and produces the
%   augmented star plot.
%
% Below we list the five options giving an informative text description
% followed by the optional input arguments. We also give in parenthesis, the
% abbreviations used in the augmented star plot, where these are different.
%
% - Robustness (rob).
% - Remove effect of order of explanatory variables by regression (scail).
% - Trapezoidal or rectangular rule for numerical integration in function
%   ctsub (trapezoid).
% - Ordering inclusion of the variables in the backfitting algorithm using
%   $R^2$ values (PredictorOrderR2).
% - Initial robust transformation of the response (tyinitial).
%
% We reduce the number of analyses for investigation by removing all those
% for which the residuals fail the Durbin-Watson and Jarque-Bera tests, at
% the $5\%$ level (two-sided for Durbin-Watson). Unlike the Durbin-Watson
% test, JarqueBera test uses a combination of estimated
% skewness and kurtosis to test for distributional shape of the residuals.
% We order the analyses by the Durbin-Watson significance level, from the
% largest value downwards. The rays in individual plots are of equal length
% for those features used in an analysis. Robustness is always plotted to
% the East and all rays are in identical places in each plot. Information
% in the plot is augmented by making the length of the rays for each
% analysis reflect the properties of the analysis; they are proportional to
% $p_{DW}$, the significance level of the Durbin-Watson test.
%
%
% Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :          Predictor variables. Matrix. Matrix of explanatory
%               variables (also called 'regressors') of dimension n x (p-1)
%               where p denotes the number of explanatory variables
%               including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless you explicitly remove it using input option
%               intercept, so do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%
%  Optional input arguments:
%
%   delrsq : termination threshold. Scalar. Iteration (in the outer loop)
%            stops when rsq changes less than delrsq in nterm. The default
%            value of delrsq is 0.01.
%           Example - 'delrsq',0.001
%           Data Types - double
%
%       l :  type of transformation. Vector. Vector of length p which
%           specifies the type of transformation for the explanatory
%           variables.
%           l(j)=1 => j-th variable assumes orderable values and any
%                   transformation (also non monotone) is allowed.
%           l(j)=2 => j-th variable assumes circular (periodic) values
%                 in the range (0.0,1.0) with period 1.0.
%           l(j)=3 => j-th variable transformation is to be monotone.
%           l(j)=4 => j-th variable transformation is to be linear.
%           l(j)=5 => j-th variable assumes categorical (unorderable) values.
%           j =1, 2,..., p+1.
%           The default value of l is a vector of ones of length p,
%           that is the procedure assumes that all the explanatory
%           variables have orderable values (and therefore can be
%           transformed non monothonically).
%           Note that in avas procedure the response is always transformed
%           (in a monothonic way).
%           Example - 'l',[3 3 1]
%           Data Types - double
%
%    maxit : maximum number of iterations for the outer loop. Scalar. The
%            default maximum number of iterations before exiting the outer
%            loop is 20.
%           Example - 'maxit',30
%           Data Types - double
%
%   nterm  : minimum number of consecutive iteration below the threshold
%           to terminate the outer loop. Positive scalar. This value
%           specifies how many consecutive iterations below the threshold
%           it is necessary to have to declare convergence in the outer
%           loop. The default value of nterm is 3.
%           Example - 'nterm',5
%           Data Types - double
%
%
%       w  : weights for the observations. Vector. Row or column vector of
%           length n containing the weights associated to each
%           observations. If w is not specified we assum $w=1$ for $i=1,
%           2, \ldots, n$.
%           Example - 'w',1:n
%           Data Types - double
%
%    plots : plots on the screen. Boolean. If plots is true it is possible
%           to visualize on the screen the augmented star plot and heatmap
%           associated to the correlation among the admissible solutions.
%           The default value of plots is false, that is no plot is shown
%           on the screen.
%           Example - 'plots',true
%           Data Types - Logical
%
%
% Output:
%
%      BestSol :  Best solutions. table.
%                 A table containing the details of the admissible
%                 solutions which have been found. We define a solution as
%                 admissible if the residuals pass the Durbin-Watson and
%                 Jarque-Bera tests, at the 5 per cent level. The rows of
%                 BestSol are ordered in a non increasing way using the
%                 p-value of the Durbin Watson test (rescaled by the value
%                 of R2 and the number of units not declared as outliers).
%                 Colums 1-5 contain boolean information about the usage of
%                 options PredictorOrderR2, scail, trapezoid, rob,
%                 tyinitial. These 5 columns are ordered in non decreasing
%                 way depending on the frequency in which a particular
%                 option has been used in the set of admissiable solutions.
%                 For example, if option PredictorOrderR2 is the only one
%                 which is always used in the set of admissible solutions,
%                 then column 1 is associated to PredictorOrderR2.
%                 Sixth column contains the value of R2.
%                 Seventh column contains the p-value of Durbin Watson
%                 test.
%                 Eight column contains the p-value of Jarque Bera test.
%                 Ninth column contains the number of units which have not
%                 been declared as outliers.
%                 10th column is a cell which contains the residuals for
%                 the associated solution.
%                 11th column contains the struct out which is the
%                 output of the call to avas.
%                 12th column contains the numbers obtained by the product
%                 of p-value of Durbin Watson test, the values of R2 and
%                 the number of units which have not been declared as
%                 outliers.
%
%   corMatrix  :  Correlation matrix. Correlation matrix among the
%                 residuals associated to the admissible solutions which
%                 have been found. The size of CorMatrix is
%                 size(BestSol,1)-times-size(BestSol,1).
%
%
% More About:
%
% Routine avas, (additivity and variance stabilization) provides
% non-parametric alternatives to the widely used Box-Cox parametric
% transformation of the response and the Box-Tidwell family of power
% transformations of explanatory variables. In routine avas it is possible:
%
% [1] To use option rob in order to replace standard OLS regression with
% robust regression in all steps of the iterative procedure. More precisely
% using robust regression we identify a subset $S_m$ of clean observations
% to use in the backfitting algorithm to calculate the transformations
% $f(X)$ at each iteration, the $n-m$ outliers being ignored. We also use
% the subset in function ctsub in the calculation of the numerical variance
% stabilising transformation. Since different response transformations can
% indicate different observations as outliers, the identification of
% outliers occurs repeatedly during our robust algorithm, once per
% iteration.
%
% [2] To use option scail in order to linearly transform the explanatory
% variables in order to let each column of matrix X have the same weight
% when predicting g(y).
%
% [3]To use option PredictorOrderR2 to completely eliminate dependence on
% the order of the variables. That is, in each iteration of the backfitting
% algorithm we impose an ordering which is based on the variable which
% produces the highest increment of $R^2$.
%
% [4] To use option tyinitial which calculates a robust parametric
% transformation of the response before starting the main backfitting loop.
%
% [5] To use option trapeziod in order to specify how to compute the
% integral to obtain trasnformed values for the units fall outside the
% range of integration determined by the units not declared as outliers.
%
% The purpose of this routine is to apply avas with and without the five
% options listed above. There are therefore 32 combinations of options that
% could be chosen. It is not obvious that all will be necessary when
% analysing all sets of data. This routine returns the list of options
% which produced the best fit and provides the necessary ingredients for a
% graphical method, the augmented star plot, to indicate which options are
% important in a particular data analysis.
%
% See also: avas.m, ace.m, aceplot.m, smothr.m, supsmu.m, ctsub.m
%
% References:
%
%
% Riani M., Atkinson A.C. and Corbellini A. (2021), Robust Transformations
% for Multiple Regression via Additivity and Variance Stabilization,
% submitted.
% Tibshirani R. (1987), Estimating optimal transformations for regression,
% "Journal of the American Statistical Association", Vol. 83, 394-405.
% Wang D.  and Murphy M. (2005), Identifying nonlinear relationships
% regression using the ACE algorithm, "Journal of Applied Statistics",
% Vol. 32, pp. 243-258.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('avasms')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:


%{
    % Example 2 in RAC (2021).
    % There are four explanatory variables and 151 observations with $x_1 $
    % equally spaced from 0:0.1:15. The linear model is:
    % 
    % \[
    % z = \text{sin}(x) + \sum_{j=2}^4 (j-1)x_j + 0.5\mathcal{U}(-0.5,0.5), 
    % \]
    %  where $x_{ij} (j=2, \ldots, 4),$ are independently $\mathcal{N}(0,1).$
    %  Eight outliers of value one replace the values of $z$ at $x = 8.9, 9.9,
    %  10.4, 10.5  ... 10.9$. The response $y = \exp z$. There are thus four
    %  explanatory variables, only one of which requires transformation, as
    %  does the response to log$y$.
    rng(100)
    x1 = (0:0.1:15)';
    n=length(x1);
    X3=randn(n,3);
    Xb=X3*(1:3)';
    z = sin(x1) + 0.5*(rand(size(x1))-0.5)+Xb;
    z([90,100,105:110]) = 1;
    X=[x1 X3];
    y=exp(z);
    % Automatic model selection
    [BestSol,corMatrix]=avasms(y,X); 
%}


%{
    %% Wang and Murphy data.
    rng('default')
    seed=100;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    y([121 80 34 188 137 110 79 86 1])=1.9+randn(9,1)*0.01;
    [BestSol,corMatrix]=avasms(y,X);
%}


%{
   %% Example of the use of avasms using the Fish data.
    load("fish.mat");
    Y=fish;
    % pike is removed. We use just the first 3 explanatory variables.
    sel=categorical(Y{:,1})~='Pike';
    y=Y{sel,2};
    X=Y{sel,3:5};
    [BestSol,corMatrix]=avasms(y,X,'l',3*ones(size(X,2),1));
%}

%% Beginning of code

if nargin <2
    error('FSDA:avasms:missingInputs','A required input argument is missing.')
end

[n,p]=size(X);

% l specifies how to transform the variables
% The first p values of l refer to the p explanatory variables.
% The last refers to the response
% 4 = linear transformation
% 3 = monotonic transformation
% ........
ll=ones(p+1,1);
%  termination threshold for outer loop
delrsq=0.01;
% maxit = max. number of iterations for the outer loop
maxit=20;
% nterm : number of consecutive iterations for which
%         rsq must change less than delrsq for convergence.
nterm=3;
w=ones(n,1);
plots=false;

UserOptions=varargin(1:2:length(varargin));

if ~isempty(UserOptions)

    options=struct('l',ll,'delrsq',delrsq,'nterm',nterm,...
        'w',w,'maxit',maxit,'plots',plots);

    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:avasms:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)

    % We now overwrite inside structure options the default values with
    % those chosen by the user
    % Notice that in order to do this we use dynamic field names
    for j=1:2:length(varargin)
        options.(varargin{j})=varargin{j+1};
    end

    ll=options.l;
    delrsq=options.delrsq;
    w=options.w;
    nterm=options.nterm;
    maxit=options.maxit;

end

tyini=struct;
tyini.la=[0 0.1 0.2 0.3 0.4 1/2];
% tyini.la=[-1:0.1:1];
tyini.la=-1:0.2:1; %   -0.5 0.5];

rob=[true false];
scail=rob;
iniFanPlot=rob;
trapezoid=rob;

if size(X,2)==1
    PredictorOrderR2=false;
else
    PredictorOrderR2=rob;
end
VAL=zeros(32,9);
Res=cell(32,1);
Out=Res;
ijkl=1;
for i=1:length(PredictorOrderR2)
    for j=1:2
        for k=1:2
            for l=1:2
                if iniFanPlot(l)==true
                    tyinitial=tyini;
                else
                    tyinitial=false;
                end
                for m=1:2
                    out=avas(y,X,'l',ll,'delrsq',delrsq,'nterm',nterm,'maxit',maxit','w',w,...
                        'PredictorOrderR2',PredictorOrderR2(i),'rob',rob(j),'scail',scail(k),'tyinitial',tyinitial,'Trapezoid',trapezoid(m));
                    if rob(j)==true && ~isempty(out.outliers)
                        nused=n-length(out.outliers);
                    else
                        nused=n;
                    end
                    valused=[PredictorOrderR2(i) rob(j) scail(k) iniFanPlot(l) trapezoid(m) out.rsq out.pvaldw out.pvaljb nused];
                    VAL(ijkl,:)=valused;
                    Res{ijkl}=(out.ty-sum(out.tX,2))';
                    Out{ijkl}=out;
                    ijkl=ijkl+1;
                end
            end
        end
    end
end
%% Build table
names=["PredictorOrderR2" "rob" "scail" "iniFanPlot" "trapezoid" "R2" "pvalDW" "pvalJB" "nused"];
VALt=array2table(VAL,'VariableNames',names);
% Add to the column the cell containing the residuals for each solution
VALt.res=Res;
VALt.Out=Out;

% plot(VALt.R2,VALt.pvalDW,'o')
% xlabel('R2')
% ylabel('DW test')


boopval=VALt.pvalDW>0.10 & VALt.pvalJB>0.10;
VALtsel=VALt(boopval,:);
if isempty(VALtsel)
    boopval=VALt.pvalDW>0.05 & VALt.pvalJB>0.05;
    VALtsel=VALt(boopval,:);
end

if isempty(VALtsel)
    warning('FSDA:avasms:noGoodModel','no model satisfies the criteria of pvalues')
end

% Reorder the first five columns
[~,sortind]=sort(sum(VALtsel{:,1:5},1),'desc');
VALtfin=VALtsel;
VALtfin=[VALtfin(:,sortind) VALtsel(:,6:end)];

ord=VALtfin{:,"pvalDW"}.*VALtfin{:,"R2"}.*VALtfin{:,"nused"}; % .*VALtfin{:,"pvalJB"}
VALtfin.ord=ord;
[~,indsor]=sort(ord,'descend');
BestSol=VALtfin(indsor,:);

lab="R2="+string(num2str(BestSol{:,"R2"},3))...
    +" n="+string(num2str(BestSol{:,"nused"}))...
    + newline ...
    +" dw="+string(num2str(BestSol{:,"pvalDW"},2))...
    +" jb="+string(num2str(BestSol{:,"pvalJB"},2));

Rescell=BestSol{:,"res"};
% n x numbsol
Resarray=cell2mat(Rescell)';
% corMatrix correlation matrix among the residuals of the solutions which
% have been found.
corMatrix=corr(Resarray);

if plots==true
    figure
    varlabs=BestSol.Properties.VariableNames(1:5);
    VALtadj=BestSol{:,1:5}.*BestSol{:,"ord"};


    d = logical(eye(size(corMatrix,2)));

    corMatrix(d)=0;

    if size(VALtadj,1)<=8
        glyphplotFS(VALtadj,'obslabels',lab,'varlabels',varlabs,'Standardize','matrix')
        hold off
        figure
        heatmap(corMatrix)

    else
        glyphplotFS(VALtadj(1:8,:),'obslabels',lab(1:8),'varlabels',varlabs,'Standardize','matrix')
        hold off
        figure
        heatmap(corMatrix(1:8,1:8))
    end

end

end

%FScategory:REG-Transformations