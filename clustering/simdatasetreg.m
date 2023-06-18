function [y,X,id]=simdatasetreg(n, Pi, Beta, S, Xdistrib, varargin)
%simdatasetReg simulates a regression dataset given the parameters of a mixture regression model
%
%
%<a href="matlab: docsearchFS('simdatasetreg')">Link to the help function</a>
%
%   [y,X,id]=simdatasetreg(n, Pi, Beta, S) generates a regression dataset
%   of size n from a mixture model with parameters 'Pi' (mixing
%   proportions), 'Beta' (matrix of regression coefficients), and 'S'
%   (vector of variances of the distributions of the points around each
%   regression hyperplane). Component sample sizes are produced as a
%   realization from a multinomial distribution with probabilities given by
%   mixing proportions. For example, if n=200, k=4 and Pi=(0.25, 0.25,
%   0.25, 0.25) function Nk1=mnrnd( n-k, Pi) is used to generate k integer
%   numbers (whose sum is n-k) from the multinominal distribution with
%   parameters n-k and Pi. The size of the groups is given by Nk1+1. The
%   first Nk1(1)+1  observations are generated using vector of regression
%   coefficients Beta(:,1) and variance S(1), ..., and the X simulated as
%   specified in structure Xdistrib, the last Nk1(k)+1 observations are
%   generated using using vector of regression coefficients Beta(:,k),
%   variance S(k) and the X simulated as specified in structure Xdistrib
%
%
%  Required input arguments:
%
%         n   : sample size  of the dataset. Scalar.
%        Pi   : vector of length k defining mixing proportions. Vector.
%               $\sum_{j=1}^k \pi=1$.
%      Beta   : p-by-k matrix containing (in the columns) regression
%               coefficients for the k groups. Matrix.
%         S   : vector of length k containing the variances of the k
%               regression hyperplanes. Vector.
%    Xdistrib : information about how to generate
%               each explanatory variable inside each group. Structure.
%               Structure which contains the following fields:
%                   Xdistrib.intercept = scalar equal to 1 if intercept is
%                       present. The default value of Xdistrib.intercept is 1.
%                   Xdistrib.type=type of distribution. Possible values for
%                       are 'normal', 'halfnormal', 'uniform', 'User'.
%                   Xdistrib.mu = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters mu
%                       for each explanatory variable and each group. The
%                       default value of Xdistrib.mu is zeros(p-1, k). Note
%                       that  Xdistrib.mu is used just if Xdistrib.type is
%                   Xdistrib.sigma = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters
%                       sigma for each explanatory variable and each group.
%                       The default value of Xdistrib.sigma is ones(p-1,k). 
%                       Notethat  Xdistrib.sigma is used just if
%                       Xdistrib.type is 'normal' or 'halfnormal'.
%                   Xdistrib.a = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters a
%                       for each explanatory variable and each group. The
%                       default value of Xdistrib.a is zeros(p-1, k).
%                       Notethat  Xdistrib.a is used just if Xdistrib.type
%                       is 'uniform', that is if we have U(a, b).
%                   Xdistrib.b = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters b
%                       for each explanatory variable and each group. The
%                       default value of Xdistrib.a is zeros(p-1, k).
%                       Note that  Xdistrib.b is used just if Xdistrib.type
%                       is 'uniform'; that is we have U(a, b).
%                   Xdistrib.X = matrix with at least n rows and p-1 (if
%                       intercept is present) or p (if intercept is not
%                       present) columns containing the values of the
%                       explanatory variables for the k groups.
%                       Notethat  Xdistrib.X is used just if Xdistrib.type
%                       is 'User'.
%               Data Types - struct
%
%  Optional input arguments:
%
%   noiseunits : number of type of outlying observations. Scalar or
%                structure. Missing value, scalar or structure.
%                This input parameter specifies the number
%                and type of outlying observations. The default value of
%                noiseunits is 0.
%                - If noiseunits is a scalar t different from 0, then t
%                  units from the uniform distribution in the interval
%                  min([X y]) max([X y]) are generated in such a way that their
%                  squared distance from the fitted value (squared residual) of each
%                  existing group is larger then the quantile 1-0.999 of
%                  the Chi^2 distribution with 1 degree of freedom. In
%                  order to generate these units the maximum number
%                  of attempts is equal to 10000.
%                - If noiseunits is a structure it may contain the following
%                  fields:
%                  number = scalar, or vector of length f. The sum of the
%                       elements of vector 'number' is equal to the total
%                       number of outliers which are simulated.
%                  alpha = scalar or vector of legth f containing the
%                       level(s) of simulated outliers. The default value
%                       of alpha is 0.001.
%                  maxiter = maximum number of trials to simulate outliers.
%                       The default value of maxiter is 10000.
%                  interval= missing value or vector of length 2, or matrix
%                         of size 2-by-2 or matrix of size 2-by-(p+1) which
%                         controls the min and max of the generated
%                         outliers for each dimension.
%                         * If interval is a vector of length 2 each outlier
%                         has a value for each column of X and y which lies
%                         inside interval(1) and interval(2).
%                         * If interval is a matrix of size 2-by-2 each
%                         outlier has a value for each column of X which
%                         lies inside interval(1,1) and interval(2,1) and a
%                         value of y which lies inside interval(1,2) and
%                         interval(2,2).
%                         * If interval is a 2-by-(p+1) matrix outliers are
%                         simulated in:
%                         interval(1,1) interval (2,1) for expl variable 1
%                         ...
%                         interval(1,p) interval (2,p) for expl variable p
%                         interval(1,p+1) interval (2,p+1) for response y.
%                         If interval is empty (default), the outliers are
%                         simulated in the interval min(X) max(X) and
%                         min(y) max (y).
%                  typeout = list of length f containing the type of
%                       outliers which must be simulated. Possible values
%                       for typeout are:
%                       * unif (or uniform), if the outliers must be
%                         generated using the uniform distribution;
%                       * norm (or normal), if the outliers must be
%                         generated using the normal distribution;
%                       * Chisquarez, if the outliers must be generated
%                         using the Chi2 distribution with z degrees of
%                         freedom;
%                       * Tz or tz, if the outliers must be generated using
%                         the Student T distribution with z degrees of
%                         freedom;
%                       * pointmass, if the outliers are concentrated on a
%                         particular point;
%                       * by_comp, if the outliers are distributed along a
%                         linear component. The option was introduced to add
%                         dense area in one linear component.
%                       * componentwise, if the outliers must have the same
%                         coordinates of the existing rows of matrix X apart
%                         from the single coordinate of y (which will be the
%                         min or max of y or to the
%                         min or max specified in interval).
%                For example, the code:
%                   noiseunits=struct;
%                   noiseunits.number=[100 100];
%                   noiseunits.typeout={'uniform' 'componentwise'};
%                   noiseunits.interval=[-2 2];
%                adds 200 outliers, the first 100 generated using
%                a uniform distribution and the last 100 using
%                componentwise scheme. Outliers are generated in the
%                interval [-2 2] for each variable.
%               Example - 'noiseunits', 10
%               Data Types - double
%
%    noisevars : Type of noise explanatory variables. Scalar or structure.
%                Empty value, scalar or structure.
%                - If noisevars is not specified or is an empty value
%                  (default) no noise variable is added to the matrix of
%                  simulated data.
%                - If noisevars is a scalar equal to r, then r new noise
%                  explnatory variables are added to the matrix of
%                  simulated data using the uniform distribution in the
%                  range [min(X) max(X)].
%                - If noisevars is a structure it may contain the following
%                  fields:
%                  noisevars.number= a scalar or a vector of length f. The sum of
%                       elements of vector 'number' is equal to the total
%                       number of noise variables to be addded.
%                  noisevars.distribution= string or cell array of strings of length
%                       f which specifies the distribution to be used to
%                       simulate the noise variables.
%                       If field distribution is not present then the
%                       uniform distribution is used to simulate the noise
%                       variables.
%                       String 'distribution' can be one of the following
%                       values:
%                       * uniform = uniform distribution
%                       * normal  = normal distribution
%                       * t or T followed by a number which controls the
%                         degrees of freedom. For example, t6 specifies to
%                         generate the data according to a Student T with 6
%                         degrees of freedom.
%                       * chisquare followed by a number which controls the
%                         degreess of freedom. For example, chisquare8
%                         specifies to generate the data according to a Chi
%                         square distribution with 8 degrees of freedom.
%                  noisevars.interval= string or vector of length 2 or matrix of size
%                         2-by-f (where f is the number of noise variables)
%                         which controls for each element of vector
%                         'number' or each element of cell 'distribution',
%                         the min and max of the noise variables. For
%                         example, interval(1,3) and interval(2,3) are
%                         respectively the minimum and maximum values of
%                         simulated the data for the third noise variable
%                         If interval is empty (default), the noise
%                         variables are simulated uniformly between the
%                         smallest and the largest coordinates of the
%                         simulated data matrix X.
%                For example, the code:
%                   noisevars=struct;
%                   noisevars.number=[3 2];
%                   noisevars.distribution={'Chisquare5' 'T3'};
%                adds 5 noise explaantory variables, the first 3 generated using
%                the Chi2 with 5 degrees of freedom and the last two
%                using the Student t with 3 degrees of freedom. Noise
%                variables are generated in the interval min(X) max(X).
%               Example - 'noisevars', 5
%               Data Types - double
%
%       lambda : Transformation coefficient. Scalar. Scalar containing
%                inverse Box-Cox transformation coefficient to apply to y.
%                The value false (default) implies that no transformation
%                is applied to response variable.
%               Example - 'lambda',2;
%               Data Types - double
%
%
%  Output:
%
%           y  : Response variable. Vector.
%                Vector of dimension (n+nout)-by-1  containing the values
%                of the responses for the k groups.
%
%           X  : Explanatory variables. Matrix. Matrix of size 
%                (n + nout)-by-(p + nnoise) containinng the values of the
%                explanatory variables for the k groups. Noise coordinates
%                are provided in the last nnoise columns.
%
%           id : classification vector. Vector. Classification vector of
%                length n + nout; 0 represents an outlier.
%
%            REMARK: If nout outliers could not be generated a warning is
%                produced. In this case matrix X and vector id will have
%                just n rows.
%
%   More About:
%
% To make a dataset more challenging for clustering, a user might want to
% simulate noise variables or outliers. Parameter 'nnoise' specifies the
% desired number of noise variables. If an interval 'int' is specified,
% noise will be simulated from a Uniform distribution on the interval given
% by 'int'. Otherwise, noise will be simulated uniformly between the
% smallest and largest coordinates of mean vectors. 'nout' specifies the
% number of observations outside (1 - 'alpha') ellipsoidal contours for the
% weighted component distributions. Outliers are simulated on a hypercube
% specified by the interval 'int'. A user can apply an inverse Box-Cox
% transformation of y providing a coefficient 'lambda'. The value 1
% implies that no transformation is needed for the response.
%
%
%
% See also: simdataset, MixSimreg
%
% References:
%
% Maitra, R. and Melnykov, V. (2010), Simulating data to study performance
% of finite mixture modeling and clustering algorithms, "The Journal of
% Computational and Graphical Statistics", Vol. 19, pp. 354-376. [to refer to
% this publication we will use "MM2010 JCGS"]
%
% Melnykov, V., Chen, W.-C. and Maitra, R. (2012), MixSim: An R Package
% for Simulating Data to Study Performance of Clustering Algorithms,
% "Journal of Statistical Software", Vol. 51, pp. 1-25.
%
% Davies, R. (1980), The distribution of a linear combination of
% chi-square random variables, "Applied Statistics", Vol. 29, pp. 323-333.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('simdatasetreg')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:
%
%{
    %% Generate mixture of regression. 
    % Use an average overlapping at centroids = 0.01 and all default options:
    % 1) Beta is generated according to random normal for each group with
    % mu=0 and sigma=1;
    % 2) X in each dimension and each group is generated according to U(0, 1);
    % 3) regression hyperplanes contain intercepts.
    % The value of p includes the intercept
    p=5;  
    k=3;
    Q=MixSimreg(k,p,'BarOmega',0.01);
    n=200;
    % Q.Xdistrib.BarX in this case has dimension 5-by-3 and is equal to
    % 1.0000    1.0000    1.0000
    % 0.5000    0.5000    0.5000
    % 0.5000    0.5000    0.5000
    % 0.5000    0.5000    0.5000
    % 0.5000    0.5000    0.5000
    % Probabilities of overlapping are evaluated at
    % Q.Beta(:,1)'*Q.Xdistrib.BarX(:,1) ... Q.Beta(:,3)'*Q.Xdistrib.BarX(:,3)
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    % spmplot([y X(:,2:end)],id)
    yXplot(y,X,'group',id);
%}

%{
    %% Generate 2 groups in 4 dimensions and add outliers from uniform distribution.
    rng('default')
    rng(100)
    p=4; % p includes the intercept
    k=2;
    out=MixSimreg(k,p,'BarOmega',0.01);
    n=300;
    noisevars=0;
    noiseunits=300;
    [y,X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib,'noisevars',noisevars,'noiseunits',noiseunits);
    yXplot(y,X,'group',id);
    suplabel('2 regression lines with outliers from uniform','t')
%}

%{
    %% Generate 4 groups in 4 dimensions and add outliers from uniform distribution.
    clear all
    close all
    rng('default')
    rng(10000)
    p=2;  % p includes the intercept
    k=4;
    out=MixSimreg(k,p,'BarOmega',0.01);
    n=300;
    noisevars=0;
    noiseunits=3000;
    [y,X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib,'noisevars',noisevars,'noiseunits',noiseunits);
    yXplot(y,X,'group',id);
     suplabel('2 regression lines with outliers from uniform','t')
%}

%{
    %% Add outliers generated from Chi2 with 5 degrees of freedom.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add asymmetric very concentrated noise
    noiseunits.typeout={'Chisquare5'};
    [y,X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib,'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    title(BigAx,'2 groups with outliers from $\chi^2_5$','Interpreter','Latex')
%}

%{
    %% Add outliers generated from Chi2 with 40 degrees of freedom.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add asymmetric concentrated noise
    noiseunits.typeout={'Chisquare40'};
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib,'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    title(BigAx,'4 groups with outliers from $\chi^2_{40}$','Interpreter','Latex')
%}

%{
    %% Add outliers generated from normal distribution.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add normal noise
    noiseunits.typeout={'normal'};
    [y,X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S,out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    title(BigAx,'4 groups with outliers from normal distribution','Interpreter','Latex')
%}

%{
    %% Add outliers generated from Student T with 5 degrees of freedom.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=3000;
    % Add outliers from T5
    noiseunits.typeout={'T5'};
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S,out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    suplabel('4 groups with outliers from Student T with 5 degrees if freedom','t')
%}

%{
    %% Add componentwise contamination.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars='';
    noiseunits=struct;
    noiseunits.number=3000;
    % Add asymmetric concentrated noise
    noiseunits.typeout={'componentwise'};
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S,out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    yXplot(y,X,'group',id);
    suplabel('4 groups with component wise outliers','t')
%}

%{
    %% Add outliers generated from Chisquare and T distribution.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=5000*ones(2,1);
    noiseunits.typeout={'Chisquare3','T20'};
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    title(BigAx,'4 groups with outliers from $\chi^2_{3}$ and $T_{20}$','Interpreter','Latex')
%}

%{
    %% Add outliers from Chisquare and T distribution and use a personalized value of alpha.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars=0;
    noiseunits=struct;
    noiseunits.number=5000*ones(2,1);
    noiseunits.typeout={'Chisquare3','T20'};
    noiseunits.alpha=0.002;
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    title(BigAx,'4 groups with outliers from $\chi^2_{3}$ and $T_{20}$','Interpreter','Latex')
%}

%{
    %% Add outliers from Chi2 and point mass contamination and add one noise variable.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noisevars=struct;
    noisevars.number=1;
    noiseunits=struct;
    noiseunits.number=[100 100];
    noiseunits.typeout={'pointmass' 'Chisquare5'};
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    title(BigAx,'4 groups with outliers from $\chi^2_{5}$ and point mass $+1$ noise var','Interpreter','Latex')
%}

%{
    %% Example of the use of personalized interval to generate outliers.
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noiseunits=struct;
    noiseunits.number=1000;
    noiseunits.typeout={'uniform'};
    % Generate outliers in the interval [-1 1] for the first variable and
    % interval [1 2] for the second variable
    noiseunits.interval=[-1 1;
                         1 2];
    % Finally add a noise variable
    noisevars=struct;
    noisevars.number=1;
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    [H,AX,BigAx]=yXplot(y,X,'group',id);
    title(BigAx,'4 groups with outliers from uniform using a personalized interval $+1$ noise var','Interpreter','Latex')
%}

%{
    % Example of the use of personalized interval to generate outliers (1).
    % Generate 1000 outliers from uniform in the interval [-2 3] and
    % 1000 units using componentwise contamination in the interval [-2 3]
    n=300;
    k=4;
    p=2;  % p includes the intercept
    out=MixSimreg(k,p,'BarOmega',0.01);
    noiseunits=struct;
    noiseunits.number=[1000 1000];
    noiseunits.typeout={'uniform' 'componentwise'};
    noiseunits.interval=[-2; 3];
    % Finally add a noise variable
    noisevars=struct;
    noisevars.number=1;
    [y, X,id]=simdatasetreg(n, out.Pi, out.Beta, out.S, out.Xdistrib, 'noisevars',noisevars,'noiseunits',noiseunits);
    yXplot(y,X,'group',id);
    suplabel('4 groups with outliers componentwise and from uniform in the interval [-2 3]','t')
%}

%{
    %% Example with user defined explanatory variables values (1).
    % Here the X distribution is the same for each component.
    clear all
    close all

    rng(1234,'twister');

    % mixture parameters
    intercept = 0; % 1/0 = intercept yes/no
    p=1+intercept;
    k=2;
    n=200;

    % beta distributed as halfnormal
    betadistrib=struct;
    betadistrib.type='HalfNormal';
    betadistrib.sigma=3;

    % explanatory variables distribution chosen by the User from a beta
    XdistribB=struct;
    XdistribB.intercept=intercept;
    XdistribB.type='User';

    X1=random('beta',1,5,n,1);             % data generation: user distribution is a beta
    XdistribB.BarX = ones(1,k)*mean(X1);   % mean of the generated data: one per group

    % overlap level baromega: chosen at random here, in a given range
    mino = 0.01; maxo = 0.1;
    baromega = mino + (maxo-mino).*rand(1,1);

    % estimated mixsim parameters
    Q=MixSimreg(k,p,'BarOmega',baromega,'Xdistrib',XdistribB,'betadistrib',betadistrib);

    % Simulate the data from the mixim parameters and the user values for X
    if intercept
        Q.Xdistrib.X = [ones(n,1) X1];
    else
        Q.Xdistrib.X = X1;
    end

    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    yXplot(y,X,'group',id,'tag','X_beta');
    set(gcf,'Name','X Beta distributed');
    title('User-defined distribution for X');
%}

%{
    %% Example with user defined explanatory variables values (2).
    % Here the X distribution is specific for each component.

    clear all
    close all

    rng(12345,'twister');

    % mixture parameters
    intercept = 0;      % 1/0 = intercept yes/no
    n=200;
    p=1+intercept;
    k=2;                %do not change k: it would not work (see below to generalise)

    % beta distributed as halfnormal
    betadistrib=struct;
    betadistrib.type='HalfNormal';
    betadistrib.sigma=3;

    % explanatory variables distribution chosen by the User from a beta
    XdistribB=struct;
    XdistribB.intercept=intercept;
    XdistribB.type='User';

    %for i=1:10
    % X beta distributed
    X2=random('beta',0.5,1,n,1);
    muBeta2 = mean(X2);

    X1=random('beta',1,0.5,n,1);
    muBeta1 = mean(X1);
    % data generation: user distribution is a beta
    XdistribB.BarX = [muBeta1 muBeta2]; % mean of the generated data: one per group

    % overlap level baromega: chosen at random here, in a given range
    mino = 0.01; maxo = 0.05;
    maxomega = mino + (maxo-mino).*rand(1,1);

    % estimated mixsim parameters
    Q=MixSimreg(k,p,'hom',true,'MaxOmega',maxomega,'Xdistrib',XdistribB,'betadistrib',betadistrib);

    % Simulate the data from the mixim parameters and the user values for X
    if intercept
        Q.Xdistrib.X = [ones(k*n,1) , [X1 ; X2]];
    else
        Q.Xdistrib.X = [X1 ; X2];
    end

    [y,X,id]=simdatasetreg(k*n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    yXplot(y,X,'group',id,'tag','X_beta');
    set(gcf,'Name','X Beta distributed');
    title('User-defined distribution for X');

%}

%% Beginning of code


if (n < 1)
    error('FSDA:simdatasetreg:Wrongn','Wrong sample size n...')
end

[p,k]=size(Beta);


% If the user selects a X distribution with more units than n, than we
% change the number of units to generate.
if strcmp(Xdistrib.type,'User') && size(Xdistrib.X,1) > n 
     n = size(Xdistrib.X,1);
     disp('Warning: size of Xdistrib.X is greater than n: ');
     disp(['         We generate size(Xdistrib.X,1) = ' num2str(n) ' units']);
end      
        
if sum(Pi <= 0)~=0 || sum(Pi > 1) ~= 0
    error('FSDA:simdatasetreg:WrongPi','Wrong vector of mixing proportions Pi: the values must be in the interval (0 1)')
end

if length(Pi) ~= k
    error('FSDA:simdatasetreg:WrongLengthPi', ['Vector of mixing proportions Pi must have a length equal to ' num2str(k)] )
end


noiseunitsdef   = '';
noisevarsdef    = '';
lambdadef='';

options=struct('noisevars',noisevarsdef,'noiseunits',noiseunitsdef,...
    'lambda',lambdadef);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:simdatasetreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:simdatasetreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for ii=1:2:length(varargin)
        options.(varargin{ii})=varargin{ii+1};
    end
    
end

lambda=options.lambda;
noiseunits = options.noiseunits;
noisevars  = options.noisevars;


if (n >= k) 
    mrr = mnrnd( n-k, Pi);
    
    % Nk contains the sizes of the clusters
    Nk = ones(1,k)+(mrr);
else
    error('FSDA:simdatasetreg:Wrongn','Sample size (n) cannot be less than the number of clusters (k)')
end


y=zeros(n,1);
id=zeros(n,1);

% Check if intercept is present
intercept=Xdistrib.intercept;
if intercept==1
    X=[ones(n,1) zeros(n,p-1)];
else
    X=zeros(n,p);
end

for j=1:k
    aa=sum(Nk(1:j-1))+1;
    bb=sum(Nk(1:j));
    id(aa:bb)=j*ones(Nk(j),1);
    
    % The dimension of a and b is p-1-by-k (if intercept is 1)
    % or p-by-k (if intercept is 0).
    if find(strcmp('Uniform',Xdistrib.type))
        if intercept==1
            Xab=bsxfun(@times, rand(Nk(j),p-1),(Xdistrib.b(:,j)-Xdistrib.a(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.a(:,j)');
            X(aa:bb,2:end)=Xab;
        else
            Xab=bsxfun(@times, rand(Nk(j),p),(Xdistrib.b(:,j)-Xdistrib.a(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.a(:,j)');
            X(aa:bb,:)=Xab;
        end
        
    elseif find(strcmp('Normal',Xdistrib.type))
        
        if intercept==1
            Xab=bsxfun(@times, randn(Nk(j),p-1),(Xdistrib.sigma(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.mu(:,j)');
            X(aa:bb,2:end)=Xab;
        else
            Xab=bsxfun(@times, randn(Nk(j),p),(Xdistrib.sigma(:,j))');
            Xab=bsxfun(@plus,Xab,Xdistrib.mu(:,j)');
            X(aa:bb,:)=Xab;
        end
        
    elseif find(strcmp('HalfNormal', Xdistrib.type))
        if intercept==1
            Xab=bsxfun(@times, abs(randn(Nk(j),p-1)),(Xdistrib.sigma(:,j))');
            X(aa:bb,2:end)=Xab;
        else
            Xab=bsxfun(@times, abs(randn(Nk(j),p)),(Xdistrib.sigma(:,j))');
            X(aa:bb,:)=Xab;
        end
        
    elseif find(strcmp('User',Xdistrib.type))
        d=find(strcmp(fieldnames(Xdistrib),'X'),1);
        if ~isempty(d)
            X = Xdistrib.X;
            % Lines below added just in case the user forgets to put the
            % constant column to account for the intercept
            constcols = find(max(X,[],1)-min(X,[],1) == 0);
            if intercept==1 && numel(constcols)==0
                X = [ones(n,1) , Xdistrib.X];
            end
        else
            error('FSDA:simdatasetreg:MissingField','If string Xdistrib = ''User'' then the user must provide input matrix X')
        end
    else
        error('FSDA:simdatasetreg:Wrongbetadistrib','Possible values for option betadistrib are ''Normal'' ''Uniform'' ''HalfNormal'' and ''User'' ')
    end
    
    % generation of the data
    y(aa:bb)=X(aa:bb,:)*Beta(:,j)+sqrt(S(j))*randn(Nk(j),1);
    
    %    X(a:b,:) = rand(Nk(j),p)*BarX(j)*2;
    %     y(a:b)=X(a:b,:)   +sqrt(S(j))*randn(Nk(j),1);
    %        X(a:b,:) = randn(Nk(j),p) +Mu(j);
    
end

% Now, we add contamination
if isstruct(noiseunits) || ~isempty(noiseunits)
    
    % Set all default values for outlier generation
    % One single outliers is generated from the uniform distribution.
    % Maximum 10000 iterations are used and this point must have a minimum
    % (squared) Mahalanobis distance from each centroid of existing groups greater
    % than 0.999 confidence level
    fnoiseunitsdef=struct;
    fnoiseunitsdef.number=1;
    fnoiseunitsdef.typeout={'uniform'};
    fnoiseunitsdef.alpha=0.001;
    fnoiseunitsdef.maxiter=10000;
    % interval to simulate outliers (if it is empty than min(X) and max(X)
    % will be used
    fnoiseunitsdef.interval='';
    
    if isstruct(noiseunits)
        
        fnoiseunits=fieldnames(noiseunits);
        
        % Check if user options inside options.fnoiseunits are valid options
        chkoptions(fnoiseunitsdef,fnoiseunits)
        
        % labeladd option
        d=find(strcmp('number',fnoiseunits));
        if d>0
            number=noiseunits.number;
            if (min(number) < 0)
                error('FSDA:simdataset:Wrongnnoise','Wrong value of outliers: it cannot be smaller than 0')
            end
        else
            number=fnoiseunitsdef.number;
        end
        
        % labeladd option
        d=find(strcmp('interval',fnoiseunits));
        if d>0
            intervalout=noiseunits.interval;
            if (isvector(intervalout) && (length(intervalout)~=2))...
                    || ~(numel(intervalout)==4 || numel(intervalout)==2*(p+1) || numel(intervalout)==2)
                error('FSDA:simdatasetreg:Wrongintervalout',...
                    ['Wrong input interval to generate outliers '...
                    '(it can be or a vector of length 2 or a matrix of size 2-by-2'...
                    ' or a matrix of size 2-by-(p+1))'])
            end
        else
            intervalout='';
        end
        
        d=find(strcmp('typeout',fnoiseunits));
        if d>0
            typeout=noiseunits.typeout;
        else
            typeout=fnoiseunitsdef.typeout;
        end
        
        d=find(strcmp('alpha',fnoiseunits));
        if d>0
            alpha=noiseunits.alpha;
            % Just in case alpha is supplied as a scalar and length(number)
            % is greater than 1, then it is necessary to resize alpha to
            % make it have length equal to length(number)
            if isscalar(alpha) &&  length(number)>1
                alpha = alpha*ones(length(number),1);
            end
        else
            alpha = fnoiseunitsdef.alpha*ones(length(number),1);
        end
        
        d=find(strcmp('maxiter',fnoiseunits));
        if d>0
            maxiter=noiseunits.maxiter;
        else
            maxiter=fnoiseunitsdef.maxiter;
        end
        
        if (number < 0)
            error('FSDA:simdatasetreg:Wrongnout','Wrong value of number: it cannot be smaller than 0')
        end
        
        if ((max(alpha) >= 1) || (min(alpha) <= 0))
            error('FSDA:simdatasetreg:WrongAlpha','Wrong value of alpha: it must be in the interval (0 1)')
        end
        
        if (maxiter < 1)
            error('FSDA:simdatasetreg:WrongMaxIter','Wrong value for maximum number of iterations: it cannot be <1')
        end
        
        % nout = total number of outliers  which has to be
        % simulated
        noiseunits=sum(number);
    else % in this case noiseunits is a scalar different from missing
        number=noiseunits;
        noiseunits=number;
        typeout=fnoiseunitsdef.typeout;
        alpha=fnoiseunitsdef.alpha;
        maxiter=fnoiseunitsdef.maxiter;
        intervalout='';
    end
    
    Xout=zeros(noiseunits, p);
    yout=zeros(noiseunits, 1);
    ni=0;
    idtmp=zeros(sum(number),1);
    for ii=1:length(number)
        typeouti=typeout{ii};
        
        [Xouti, youti,  ~] = getOutliersreg(number(ii), Beta, S, alpha(ii), maxiter, typeouti, intervalout,ii);
        Xout(ni+1:ni+size(Xouti,1),:)=Xouti;
        yout(ni+1:ni+size(Xouti,1),:)=youti;
        
        idtmp(ni+1:ni+size(Xouti,1))=-ii;
        ni=ni+size(Xouti,1);
    end
    
    if ni<noiseunits
        warning('FSDA:simdataset:Modifiedn',['Output matrix X will have just ' num2str(n+ni) ...
            ' rows and not ' num2str(n+noiseunits)])
        noiseunits=ni;
        
    end
    X =[X;Xout(1:ni,:)];
    y =[y;yout(1:ni)];
    
    id =[id;idtmp(1:ni)];
else
    noiseunits=0;
end

if isstruct(noisevars) || ~isempty(noisevars)
    % Set all default values for noise variable generation
    % One single noise variable is generated from the uniform distribution.
    % The values of the noise variable range from min(min(Mu)) and max(max(Mu))
    noisevarsdef=struct;
    noisevarsdef.number=1;
    noisevarsdef.distribution={'uniform'};
    noisevarsdef.interval='';
    
    if isstruct(noisevars)
        fnoisevars=fieldnames(noisevars);
        
        % Check if user options inside options.fnoisevars are valid options
        chkoptions(noisevarsdef,fnoisevars)
        
        % number option
        d=find(strcmp('number',fnoisevars));
        if d>0
            number=noisevars.number;
            if (min(number) < 0)
                error('FSDA:simdatasetreg:Wrongnnoise','Wrong value of number of noisevars: it cannot be smaller than 0')
            end
        else
            number=noisevarsdef.number;
        end
        
        % distribution
        d=find(strcmp('distribution',fnoisevars));
        if d>0
            distribution=noisevars.distribution;
        else
            distribution=noisevarsdef.distribution;
        end
        
        % interval
        d=find(strcmp('interval',fnoisevars));
        if d>0
            interval=noisevars.interval;
        else
            interval=noisevarsdef.interval;
        end
        
        % nvars = total number of noise variables which has to be
        % simulated
        nvars=sum(number);
        
        %  if noisevars ~= 0
        if isempty(interval)
            L = min(min(X));
            U = max(max(X));
            L = L* ones(1,nvars);
            U = U* ones(1,nvars);
        else
            L = interval(1,:);
            U = interval(2,:);
        end
    else % in this case noisevars is a scalar different from missing
        L = min(min(X));
        U = max(max(X));
        number=noisevars;
        nvars=number;
        distribution={'uniform'};
    end
    
    
    rrr=zeros(n + noiseunits, nvars);
    for ii=1:length(number)
        distributioni=distribution{ii};
        
        if strcmp(distributioni(1),'T') || strcmp(distributioni(1),'t')
            nu=str2double(distributioni(2:end));
            rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescaleFS(trnd(nu, n + noiseunits,  number(ii)));
        elseif strcmp(distributioni,'norm') || strcmp(distributioni,'normal')
            % data generated from the normal distribution rescaled in the
            % interval [0 1]
            rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescaleFS(randn(n + noiseunits,  number(ii)));
        elseif strcmp(distributioni,'uniform')
            rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rand(n + noiseunits,  number(ii));
        elseif strcmp(distributioni(1:9),'Chisquare')
            nu=str2double(distributioni(10:end));
            rrr(:,sum(number(1:ii-1))+1:sum(number(1:ii))) = rescaleFS(chi2rnd(nu, n + noiseunits,  number(ii)));
        else
            error('FSDA:simdatasetreg:WrongDistrib','Variable distribution type not supported')
        end
    end
    
    % Values of noise variables were constrained to lie in the interval [0 1]
    % Now we rescale them to the interval [L U]
    Xnoise=bsxfun(@times,rrr,U-L);
    Xnoise=bsxfun(@plus,Xnoise,L);
    X = [X, Xnoise];
end

if ~isempty(lambda)
    if ~isscalar(lambda)
        error('FSDA:simdatasetreg:WrongLambda','lambda must be  a scalar (just y is transformed)')
    end
    
    y = (lambda *y + 1).^(1/lambda) - 1;
    if (sum(isnan(y)) ~= 0)
        warning('FSDA:simdatasetreg:NaNs','NaNs were produced during transformation')
    end
end


%% Inner functions
% Xout (yout) with nout rows which contains the outliers.
% fail = scalar. If fail =1 than it was not possible to generate the
% outliers in the interval specified by input option int in maxiter trials
% else fail = 0
    function [Xout,yout, fail] = getOutliersreg(nout, Beta, S, alpha, maxiter, typeout, intervalout,current_group)
        fail = 0;
        % maxiter = maximum number of iterations to generate outliers
        critval =chi2inv(1-alpha,1);
        
        Xout = zeros(nout,p);
        yout = zeros(nout,1);
        
        % If intervalout is not specified, just simulate outliers in the
        % interval min(X) and max(X) else use L and U supplied by the ser
        if isempty(intervalout)
            Lout = min(X);
            Uout = max(X);
            Louty = min(y);
            Uouty = max(y);
            
        else
            % Just in case a two elements vector is supplied instead of a
            % matrix of size 2-by-2 or 2-by-(p+1)
            if isvector(intervalout)
                Lout=intervalout(1)*ones(1,p);
                Uout=intervalout(2)*ones(1,p);
                Louty = Lout(1,1);
                Uouty = Uout(1,1);
            elseif  numel(intervalout)==4
                Lout=intervalout(1,1)*ones(1,p);
                Uout=intervalout(2,1)*ones(1,p);
                Louty = intervalout(1,2);
                Uouty = intervalout(2,2);
            else
                Lout=intervalout(1,p);
                Uout=intervalout(2,p);
                Louty=intervalout(1,p+1);
                Uouty=intervalout(2,p+1);
            end
        end
        
        i = 1;
        
        % Remark: maxiter1 must be much greater than nout
        if nout<2000
            maxiter1=20000;
        else
            maxiter1=nout*10;
            maxiter=maxiter1;
        end
        
        % generate y values rrally
        if strcmp(typeout(1),'T') || strcmp(typeout(1),'t')
            nuT=str2double(typeout(2:end));
            rrally = rescaleFS(trnd(nuT, maxiter1, 1));
        elseif strcmp(typeout(1:4),'norm')
            rrally = rescaleFS(randn(maxiter1,1));
        elseif strcmp(typeout(1:4),'unif')
            rrally = rand(maxiter1,1);
        elseif   strcmp(typeout,'by_comp')
            %rrally = rand(maxiter1,1);
        elseif strcmp(typeout(1:9),'Chisquare')
            nuC=str2double(typeout(10:end));
            rrally = rescaleFS(chi2rnd(nuC,maxiter1,1));
        elseif   strcmp(typeout,'pointmass')
            rrally=rand(maxiter1,1);
        elseif   strcmp(typeout,'componentwise')
            % component wise contamination
        else
            error('FSDA:simdatasetreg:WrongDistrib','Outlier distribution type not supported')
        end
        
        % X values are generated from uniform distribution
        rrallX = rand(maxiter1, p);
        
        
        iter=0;
        while (i <= nout  &&  iter<maxiter)
            iter=iter+1;
            
            if strcmp(typeout,'componentwise')
                % extract one unit among those already extracted and
                % contaminate just a single random coordinate
                
                Xout(i,:)=X(randsample(1:n,1),:);
                
                if rand(1,1)>0.5
                    yout(i)=Uouty;
                else
                    yout(i)=Louty;
                end
            else
                
                % extract one unit from rrallX and rrally
                rindex=randsample(maxiter1,1);
                rrX=rrallX(rindex,:);
                % Rescale the unit in the interval Lout and Uout for X
                Xout(i,:) = (Uout-Lout).*rrX+Lout;
                if strcmp(typeout,'by_comp')
                    yout(i) = Xout(i,:) * Beta(:,current_group) + sqrt(S(current_group))* rand/max(Xout(i,:),0.5);
                else
                    rry=rrally(rindex);
                    % Rescale the unit in the interval Louty and Uouty for y
                    yout(i) = (Uouty-Louty).*rry+Louty;
                end
                
            end
            
            % With the intercep, the first column of Xout is always 1 (constant)
            if intercept==1
                Xout(i,1)=1;
            end
            
            % calculate the distance of each potential outlier from the k
            % linear regression components. Each distance must be greater
            % than critval
            ij=0;
            
            if ~strcmp(typeout,'by_comp')
                for jj=1:k
                    yhatij=Xout(i,:)*Beta(:,jj);
                    if ((yout(i)-yhatij)^2)/(S(jj)) <critval
                        
                        ij=1;
                        break
                    end
                end
            end
            % disp(jj)
            
            if ij==0
                i = i + 1;
            end
            
            if ij==0 && strcmp(typeout,'pointmass')
                % Row of point mass contamination which has been found is
                % replicated nout times and stored inside matrix Xout
                % In other words, Xout will have nout equal rows
                Xalleq=repmat(Xout(1,:),nout,1);
                yalleq=repmat(yout(1),nout,1);
                Xout=Xalleq;
                yout=yalleq;
                break
            end
        end
        
        % If iter = maxiter then it was not possible  to generate nout
        % outliers in maxiter simulations.
        if iter== maxiter
            disp(['Warning: it was not possible to generate ' num2str(nout) ' outliers'])
            disp(['in ' num2str(maxiter) ' replicates in the interval [' num2str(Lout(1)) ...
                '--' num2str(Uout(1)) ']'])
            disp(['Number of values which was possible to generate is equal to ' num2str(i)])
            disp('Please modify the type of outliers using option ''typeout'' ')
            disp('or increase input option ''alpha''')
            disp(['The value of alpha now is ' num2str(alpha)]);
            disp(['Outliers have been generated according to ' typeout])
            % If max number of iteration has been reached fail is 1
            fail=1;
            Xout=Xout(1:i,:);
            yout=yout(1:i);
        end
    end

end
%FScategory:CLUS-MixSim
