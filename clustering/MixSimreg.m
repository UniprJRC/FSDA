function [out] = MixSimreg(k,p,varargin)
%MixSimreg generates k regression hyperplanes in p dimensions with given overlap
%
%<a href="matlab: docsearchFS('MixSimreg')">Link to the help function</a>
%
%  MixSimreg(k,p) generates k groups in p dimensions. It is possible to
%  control the average and maximum or standard deviation of overlapping.
%
%  Notation and background.
%
%  Given two generic clusters $i$ and $j$ with $i \ne j=1,...,k$, indexed by
%  $\phi(x,\mu_i,\sigma_i^2)$ and $\phi(x,\mu_j, \sigma_j^2)$ with probabilities of
%  occurrence $\pi_i$ and $\pi_j$, the misclassification probability with
%  respect to cluster $i$ (denoted with $w_{j|i}$) is defined as 
%  \[ 
%   Pr[\pi_i \phi(x,\mu_i,\sigma_i^2) < \pi_j \phi(x,\mu_j,\sigma_j^2)] 
%   \]
%  where, in the regression context, $\mu_i={\overline x}_i' \beta_i$ and
%  $\mu_j= \overline x_j' \beta_j$. We assume that the length of vectors $x_i$,
%  $x_j$, $\beta_i$, and $\beta_j$ is $p$ (number of explanatory variables
%  including or excluding intercept). In our implementation, the
%  distribution of the elements of vectors $\beta_i$ ($\beta_j$) can be 'Normal'
%  (with parameters $\mu$ and $\sigma$), 'HalfNormal' (with parameter $\sigma$) or
%  uniform (with parameters $a$ and $b$). Same thing for the distribution of
%  the elements of $x_i$ ($x_j$). However, while the parameters of the
%  distributions are the same for all elements of $\beta$ in all groups, the
%  parameters of the distribution of the elements of vectors $x_i$ ($x_j$) can
%  vary for each group and each explanatory variable. In other words, it is
%  possible to specify (say) that the distribution of the second
%  explanatory variable in the first group is $U(2, 3)$ while the distribution
%  of the third explanatory variable in the second group is $U(2, 10)$.
%
%  The matrix containing the misclassification probabilities $w_{j|i}$ is
%  called OmegaMap.
%  The probability of overlapping between groups i and j is given by
%      \[      w_{j|i} + w_{i|j}    \qquad       i,j=1,2, ..., k      \]
%  The diagonal elements of OmegaMap are equal to 1.
%  The average overlap (BarOmega, in the code) is defined as the sum of the
%  off diagonal elements of OmegaMap (containing the misclassification
%  probabilities) divided by $k*(k-1)/2$.
%  The maximum overlap (MaxOmega, in the code) is defined as:
%     \[       \max (w_{j|i} + w_{i|j})  \qquad i \ne j=1,2, ..., k   \]
%
%  The probability of overlapping $w_{j|i}$ is nothing but the cdf of a linear
%  combination of non central $\chi^2$ distributions with 1 degree of freedom
%  plus a linear combination of $N(0,1)$ evaluated in a point $c$.
%  The coefficients of the linear combinations of non central $\chi^2$ and
%  $N(0,1)$ and point c depend on 
%  $\sigma_{j|i}^2 = \sigma_i^2/\sigma_j^2$.
%  This probability is computed using routine ncx2mixtcdf
%
%  Required input arguments:
%
%            k: Number of groups (components). Scalar.
%               Desired number of groups.
%               Data Types - int16|int32|int64|single|double
%            p: Number of explanatory variables for each regression
%               hyperplane (including intercept). Scalar.
%               Desired number of variables.
%               Data Types - int16|int32|int64|single|double
%
%  Optional input arguments:
%
%    BarOmega : Requested average overlap. Scalar. Value of desired average
%               overlap. The default value is ''
%               Example - 'BarOmega',0.05 
%               Data Types - double
%
%    MaxOmega : Requested maximum overlap. Scalar. Value of desired maximum
%               overlap. If BarOmega is empty the default value of MaxOmega
%               is 0.15.
%               Example - 'MaxOmega',0.05 
%               Data Types - double
%
%    StdOmega : Requested std of overlap. Scalar. Value of desired standard
%               deviation of overlap.
%               Remark 1: The probability of overlapping between two
%               clusters $i$ and $j$ ($i \ne j=1, 2, ..., k$), called $p_{ij}$, is
%               defined as the sum of the two misclassification
%               probabilities $p_{ij}=w_{j|i} + w_{i|j}$.
%
%               Remark 2: it is possible to specify up to two values among
%               BarOmega MaxOmega and StdOmega.
%               Example - 'StdOmega',0.05 
%               Data Types - double
%
%         hom : Equal Sigmas. Scalar boolean. 
%               Scalar boolean which specifies if the desired clusters have
%               to be heterogeneous or homogeneous:
%               hom=false (default) ==> heterogeneous,
%               hom=true            ==> homogeneous \Sigma_1 = ... = \Sigma_k
%               Example - 'hom',false 
%               Data Types - boolean
%
%  restrfactor: restriction factor. Scalar. 
%               Scalar in the interval $[1, \infty]$ which specifies the
%               maximum ratio to allow between the largest $\sigma^2$ and
%               the smallest $\sigma^2$ which are generated. If, for example,
%               restrfactor=10, after generating the mixtures of regression
%               lines (hyperplanes) we check that the ratio
%                     $\sigma^2_i/\sigma^2_j$, $i \ne j=1, ..., k$,
%               is not larger than restrfactor. In order to apply this
%               restriction, which is typical of tclust.m, we call routine
%               restreigen.m.
%               Example - 'restrfactor',8 
%               Data Types - double
%
%       PiLow : Smallest miximg proportion. Scalar. 
%               Value of the smallest mixing proportion (if 'PiLow'
%               is not reachable with respect to k, equal proportions are
%               taken; PiLow = 1.0 implies equal proportions by default).
%               PiLow must be a number in the interval (0 1]
%               Example - 'PiLow',0.1 
%               Data Types - double
%
%    Xdistrib : distribution to use for each explanatory variable. Scalar
%               or structure. It specifies the distribution to use
%               for each explanatory variable and each group. Once chosen,
%               the distribution is fixed for each explanatory variable and
%               each group; however, the parameters of the chosen
%               distribution may vary across variables and groups. For
%               example, once decided that the distibution of the X is
%               uniform, the second variable of the first group can be
%               defined in [a21 b21] while the third variable of the second
%               group can be defined in [a32 b32].
%               - If Xdistrib = 1 the default is to assume that the explanatory
%                 variables come from U(0, 1) and that the first explanatory
%                 variable is a constant term (the intercept).
%               - If Xdistrib is a structure, it may contain information about
%                 the distribution (in the fieldname 'type') and the
%                 parameters of the distribution. The following options are
%                 admitted for Xdistrib:
%                   Xdistrib.intercept = scalar equal to 1 if intercept is
%                   present. The default value of Xdistrib.intercept is 1.
%                   The other fields of Xdistrib depend on the distribution
%                   which is chosen.
%                   Xdistrib.type = string which identifies the kind of distribution. 
%                   Possibile values are:
%                   'Normal'; NORMAL DISTRIBUTION N(mu, sigma); In this
%                   case the use must supply mu and sigma.
%                   'Uniform'; UNIFORM DISTRIBUTION U(a, b).
%                   'HalfNormal'; HALF NORMAL DISTRIBUTION Half(sigma)= |N(0 sigma)|.
%                   'User'.  OTHER DISTRIBUTION. In this case the user must directly provide
%                   means of the p explanatory variables for each group.
%                   Xdistrib.mu   = matrix of size (p-1)-by-k if
%                        Xdistrib.intercept=1 or p-by-k if
%                       Xdistrib.intercept=0 containing the parameters mu
%                       for each explanatory variable and each group. The
%                       default value of Xdistrib.mu is 0.5*ones(p-1, k).
%                       Xdistrib.mu is used just if X.distrib.type is normal.
%                   Xdistrib.sigma = matrix of size (p-1)-by-k if
%                       (Xdistrib.intercept=1) or p-by-k if
%                       (Xdistrib.intercept=0) containing the parameters
%                       sigma for each explanatory variable and each group.
%                       The default value of Xdistrib.sigma is ones(p-1,k).
%                       Xdistrib.sigma is used just if X.distrib.type is
%                       'Normal' or if Xdistrib.type is 'HalfNormal'.
%                   Xdistrib.a    = matrix of size (p-1)-by-k 
%                       if Xdistrib.intercept=1 or p-by-k if
%                       Xdistrib.intercept=0 containing the parameters
%                       a (the lower limits) for each explanatory variable
%                       and each group. The default value of Xdistrib.a is
%                       zeros(p-1, k).  
%                       Xdistrib.a is used just if Xdistrib.type is 'Uniform'
%                   Xdistrib.b = matrix of size (p-1)-by-k if
%                       Xdistrib.intercept=1 or p-by-k if
%                       Xdistrib.intercept=0 containing the parameters b
%                       (the upper limits) for each explanatory variable
%                       and each group. The default value of Xdistrib.b is
%                       ones(p-1, k). 
%                       Xdistrib.b is used just if Xdistrib.type is 'Uniform'
%                   Xdistrib.BarX= (p-1)-by k matrix if intercept is present
%                   or p-by-k matrix if intercept is not present containing the
%                   means of the p explanatory variables for each group.
%                       Xdistrib.BarX is used just if Xdistrib.type is 'User'
%                 Example - 'Xdistrib',1 
%                 Data Types - double
%
% betadistrib : distribution to use for regression coefficients. 
%               Scalar or structure. It specifies the distribution to use
%               for each element of the vectors of regression coefficients.
%               Scalar or structure.
%               Once chosen, the distribution together with its parameters
%               is fixed for each element of beta, across each group.
%               - If betadistrib = 1 the default is to assume that the vector
%                 of regression coefficients come from N(0, 1).
%               - If betadistrib is a structure it may contain information
%                 about the distribution (in the fieldname type) and the
%                 parameters of the distribution.
%                 The following options are admitted for betadistrib:
%                  betadistrib.type = string which identifies the kind of distribution. 
%                   Possibile values are:
%                   'Normal'; NORMAL DISTRIBUTION N(mu, sigma); In this
%                   case the user must supply mu and sigma.
%                   'Uniform'; UNIFORM DISTRIBUTION U(a, b).
%                   'HalfNormal'; HALF NORMAL DISTRIBUTION Half(sigma)= |N(0 sigma)|.
%                   'User'.  OTHER DISTRIBUTION. 
%                   betadistrib.mu = scalar, containing parameter mu for the
%                       distribution of each element of beta across each
%                       group. The default value of betadistrib.mu is 0.
%                       betadistrib.mu is used just if betadistrib.type is normal.
%                   betadistrib.sigma = scalar, containing parameter sigma for
%                       the distribution of each element of beta across
%                       each group. The default value of betadistrib.sigma
%                       is 1.
%                       betadistrib.sigma is used just if betadistrib.type
%                       is Normal or if betadistrib is HalfNormal.
%                   betadistrib.a = scalar, containing parameter a for the
%                     distribution of each element of beta across each
%                     group. The default value of betadistrib.a is 0.
%                       betadistrib.a is used just if betadistrib.type
%                       is 'Uniform'.
%                    betadistrib.b = scalar, containing parameter b for
%                     the distribution of each element of beta across
%                     each group. The default value of betadistrib.b is 1.
%                       betadistrib.a is used just if betadistrib.type
%                       is 'Uniform'.
%                     betadistrib.Beta = matrix of size (p-1)-by k
%                     (if intercept is present) or p-by-k (if intercept is
%                     not present) containing the vectors of regression
%                     coefficients for the k groups.
%                       betadistrib.Beta is used just if betadistrib.type
%                       is 'User'.
%                 Example - 'betadistrib',1 
%                 Data Types - double
%
%        resN : maximum number of attempts. Scalar integer.
%               Maximum number of mixture re-simulations to find a
%               simulation setting with prespecified level of overlapping.
%               The default value of resN is 100.
%                 Example - 'resN',3 
%                 Data Types - double
%
%         tol : tolerance. Vector of length 2. 
%               - tol(1) (which will be called tolmap) specifies
%                 the tolerance between the requested and empirical
%                 misclassification probabilities (default is 1e-06).
%               - tol(2) (which will be called tolnxc2) specifies the
%                 tolerance to use in routine ncx2mixtcdf (which computes
%                 the cdf of linear combinations of non central chi2
%                 distributions). The default value of tol(2) 1e-06.
%                 Example - 'tol',[0.01 0.02] 
%                 Data Types - double
%
%         lim : maximum number of integration terms to use inside routine
%               ncx2mixtcdf. Integer. Default is 1e06.
%                 Example - 'lim',0.001 
%                 Data Types - double
%               REMARK: Parameters tolncx2=tol(2) and lim are used by
%               function ncx2mixtcdf.m which computes the cdf of a linear
%               combination of non central chi2 r.v.. This is the
%               probability of misclassification.
%
%     Display : Level of display. Logical.
%               - 'off' displays no output.
%               - 'notify' (default) displays output if requested
%                  overlap cannot be reached in a particular simulation.
%               - 'iter' displays output at each iteration of each
%                 simulation.
%                 Example - 'Display','off' 
%                 Data Types - char
%
%       Remark: The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%       Remark: If 'BarOmega' is not specified, the function generates a
%               mixture solely based on 'MaxOmega';
%               if 'MaxOmega' is not specified, the function generates a
%               mixture solely based on 'BarOmega'.
%               If both BarOmega and MaxOmega are not specified the
%               function generates a mixture using MaxOmega=0.15.
%               If both BarOmega and MaxOmega are empty values
%               (e.g. out=MixSimreg(3,4,'MaxOmega','','BarOmega','')
%               the following message appears on the screen
%               Error. At least one overlap characteristic between BarOmega
%               and MaxOmega should be specified...
%  Output:
%
%         out:   structure which contains the following fields
%
%       out.OmegaMap = matrix of misclassification probabilities (k-by-k);
%                      OmegaMap(i,j) = $w_{j|i}$ is the probability that X,
%                      coming from the $i$-th component (group), is classified
%                      to the $j$-th component.
%
%       out.BarOmega = scalar. Value of average overlap. BarOmega is computed
%                      as (sum(sum(OmegaMap))-k)/(0.5*k(k-1))
%
%       out.MaxOmega = scalar. Value of maximum overlap. MaxOmega is the
%                      maximum of OmegaMap(i,j)+OmegaMap(j,i)
%                      (i ~= j)=1, 2, ..., k. In other words, MaxOmega=
%                      OmegaMap(rcMax(1),rcMax(2))+OmegaMap(rcMax(2),rcMax(1))
%
%       out.StdOmega = scalar. Value of standard deviation (std) of overlap.
%                      StdOmega is the standard deviation of the k*(k-1)/2
%                      probabilities of overlapping.
%
%         out.rcMax  = vector of length 2. It containes the row and column
%                      numbers associated with the pair of components
%                      producing maximum overlap 'MaxOmega'
%
%              fail  = scalar, flag value. 0 indicates a successful mixture
%                      generation, 1 represents failure.
%
%            out.Pi  = vector of length k containing the mixing proportions.
%                      Clearly, sum(out.Pi)=1.
%
%          out.Beta = p-by-k matrix containing (in each column) the
%                      regression coefficients for each group.
%
%            out.Mu  = vector of length k, consisting of components' mean vectors
%                      for each regression hyperplane.
%                      out.Mu(1)=BarX'Beta(:,1) ... out.Mu(p)=BarX'Beta(:,k)
%
%             out.S =  k-by-1 vector containing the variances for the k
%                      groups.

%
% See also tkmeans, tclust, tclustreg, lga, rlga, ncx2mixtcdf, restreigen
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
% Parlett, B.N. and Reinsch, C. (1969), Balancing a matrix for calculation
% of eigenvalues and eigenvectors, "Numerische Mathematik", Vol. 13,
% pp. 293-304.
% Parlett, B.N. and Reinsch, C. (1971), Balancing a Matrix for Calculation of
% Eigenvalues and Eigenvectors, in Bauer, F.L. Eds, "Handbook for Automatic
% Computation", Vol. 2, pp. 315-326, Springer.
%   Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A. (2008), 
%   A General Trimming Approach to Robust Cluster Analysis. Annals
%   of Statistics, Vol. 36, 1324-1345. [Technical Report available at:
%   http://www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf].
% Torti F., Perrotta D., Riani, M. and Cerioli A. (2018). Assessing Robust
% Methodologies for Clustering Linear Regression Data, "Advances in Data
% Analysis and Classification". [doi
% https://doi.org/10.1007/s11634-018-0331-4].

% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('MixSimreg')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples:
%
%{
    %% Example 1: Mixture of regression with prefixed average overlap.
    % Generate mixture of regression using an average overlapping at
    % centroids =0.01. Use all default options
    % 1) Beta is generated according to random normal for each group with
    % mu=0 and sigma=1
    % 2) X in each dimension and each group is generated according to U(0, 1)
    % 3) regression hyperplanes contain intercepts
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
    yXplot(y,X(:,2:end),'group',id);
%}

%{
    % Example 2: Mixture of regression with prefixed average overlap.
    % Generate mixture of regression hyperplanes using an average overlapping at
    % centroids =0.01.
    % 1) we use all the default options for Beta (random normal for each group with
    % mu=0.5 and sigma=1)
    % 2) X in the second dimension for the third group is generated according to U(1, 3)
    rng(10,'twister')
    % Specify the distribution of the explanatory variables
    Xdistrib=struct;
    Xdistrib.type='Uniform';
    Xdistrib.a=zeros(p-1,k);
    Xdistrib.a(2,3)=1;
    Xdistrib.b=ones(p-1,k);
    Xdistrib.b(2,3)=3;
    % 3) regression hyperplanes contain intercepts
    Q=MixSimreg(k,p,'BarOmega',0.01,'Xdistrib',Xdistrib);
    n=200;
    % Q.Xdistrib.BarX in this case has dimension 5-by-3 and is equal to
    %     1.0000    1.0000    1.0000
    %     0.5000    0.5000    0.5000
    %     0.5000    0.5000    2.0000
    %     0.5000    0.5000    0.5000
    %     0.5000    0.5000    0.5000
    % Probabilitties of overlapping are evaluated at
    % Q.Beta(:,1)'*Q.Xdistrib.BarX(:,1) ... Q.Beta(:,3)'*Q.Xdistrib.BarX(:,3)
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    yXplot(y,X(:,2:end),'group',id);
%}

%{
    % Example 3: Mixture of regression with prefixed average overlap.
    % Exactly as before but now the distribution of beta is N(0 6)
    rng(10,'twister')
    p=5;
    k=3;
    % Specify the distribution of the explanatory variables
    Xdistrib=struct;
    Xdistrib.type='Uniform';
    Xdistrib.a=zeros(p-1,k);
    Xdistrib.a(2,3)=1;
    Xdistrib.b=ones(p-1,k);
    Xdistrib.b(2,3)=3;
    % Specify the distribution of the beta coefficients
    betadistrib=struct;
    betadistrib.type='Normal';
    betadistrib.sigma=6;
    Q=MixSimreg(k,p,'BarOmega',0.01,'Xdistrib',Xdistrib,'betadistrib',betadistrib);
    n=200;
    % Probabilitties of overlapping are evaluated at
    % Q.Beta(:,1)'*Q.Xdistrib.BarX(:,1) ... Q.Beta(:,3)'*Q.Xdistrib.BarX(:,3)
    % Q.betadistrib in this case is equal to
    %      type: 'Normal'
    %     sigma: 6
    %        mu: 0.5000
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    yXplot(y,X(:,2:end),'group',id)
%}

%{
    % Example 4: Internation trade data example.
    % All slopes are positive (beta generated using half normal) p=1 and there
    % is no intercept
    rng(10,'twister')
    p=1;
    k=5;
    Xdistrib=struct;
    Xdistrib.intercept=0;
    Xdistrib.type='Uniform';
    Xdistrib.a=zeros(p,k);
    Xdistrib.b=10*ones(p,k);

    betadistrib=struct;
    betadistrib.type='HalfNormal';
    betadistrib.sigma=6;
    Q=MixSimreg(k,p,'BarOmega',0.01,'Xdistrib',Xdistrib,'betadistrib',betadistrib);
    n=200;

    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    yXplot(y,X,'group',id)
%}


%{
    % Example 5:  Another international trade data example.
    % Here the strips of certain groups are limited up to certain values
    % There is no intercept.
    % In this example we compare high and low overlap among regression hyperplanes
    p=1;
    k=4;
    Xdistrib=struct;
    Xdistrib.intercept=0;
    Xdistrib.type='Uniform';
    Xdistrib.a=zeros(p,k);
    Xdistrib.b=[4 2 10 5];

    betadistrib=struct;
    betadistrib.type='HalfNormal';
    betadistrib.sigma=6;
    n=200;

    % Strong overlap BarOmega=0.2
    close all
    rng(10,'twister')
    Q=MixSimreg(k,p,'BarOmega',0.2,'Xdistrib',Xdistrib,'betadistrib',betadistrib);
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    yXplot(y,X,'group',id,'tag','Strong_Overlap')
    set(gcf,'Name','Strong overlap')

    % Small overlap BarOmega=0.01
    rng(10,'twister')
    Q=MixSimreg(k,p,'BarOmega',0.01,'Xdistrib',Xdistrib,'betadistrib',betadistrib);
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    yXplot(y,X,'group',id,'tag','Small_Overlap')
    set(gcf,'Name','Small overlap')
%}

%{ 
    % Example 6: Betadistrib with a specific parameter for each group.
    rng(10,'twister')
    clear all
    close all
    intercept = 1;
    n=200;
    p=1;
    k=2;

    Xdistrib=struct;
    Xdistrib.intercept=intercept;
    Xdistrib.type='Uniform';
    Xdistrib.a=zeros(p,k);
    Xdistrib.b=10*ones(p,k);

    % beta distributed as HalfNormal
    betadistrib=struct;
    betadistrib.type='HalfNormal';
    betadistrib.sigma=6;

    % beta distributed as normal
    betadistrib1=struct;
    betadistrib1.type='Normal';
    betadistrib2 = betadistrib1;
    betadistrib1.mu=  -1;
    betadistrib1.sigma=1;
    betadistrib2.mu=[-1 , 1];
    betadistrib2.sigma=[0.01 , 1];

    Q1=MixSimreg(k,p+intercept,'BarOmega',0.01,'Xdistrib',Xdistrib,'betadistrib',betadistrib1);
    Q2=MixSimreg(k,p+intercept,'BarOmega',0.01,'Xdistrib',Xdistrib,'betadistrib',betadistrib2);

    [y1,X1,id1]=simdatasetreg(n,Q1.Pi,Q1.Beta,Q1.S,Q1.Xdistrib);
    [y2,X2,id2]=simdatasetreg(n,Q2.Pi,Q2.Beta,Q2.S,Q2.Xdistrib);
    yXplot(y1,X1,'group',id1,'tag','scalar')
    title('Betadistrib is a scalar: same parameters for all betas')
    yXplot(y2,X2,'group',id2,'tag','vector')
    title('Betadistrib is a vector: a parameter for each beta')
    cascade;

%}

%% Beginning of code 

% User options

% Default
if nargin<2
    error('FSDA:MixSimreg:Missingp','k=number of components and p = number of explanatory variables (includind intercept) must be specified');
end

if (p < 1)
    error('FSDA:MixSimreg:Wrongp','Wrong number of explanatory variables p')
end

if k<=1
    error('FSDA:MixSimreg:Wrongk','Wrong number of mixture components k')
end

BarOmegadef = '';
MaxOmegadef = 0.15;
StdOmegadef = '';
PiLowdef    = 0;
resNdef     = 100;
toldef      = [1e-06; 1e-06];
limdef      = 1e06;
restrfactordef='';

options=struct( 'BarOmega',BarOmegadef,'MaxOmega',MaxOmegadef,...
    'StdOmega',StdOmegadef, 'Display', 'notify', ...
    'hom',false,'PiLow',PiLowdef,'betadistrib', 1,'Xdistrib',1, ...
    'resN',resNdef,'tol',toldef,'lim',limdef,'restrfactor',restrfactordef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:MixSimreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present in
    % structure options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:MixSimreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin > 2
    
    % If the user inside user options has only specified BarOmega but not
    % MaxOmega then MaxOmega is initialied with an empty value
    checkBarOmega = strcmp(UserOptions,'BarOmega');
    checkMaxOmega = strcmp(UserOptions,'MaxOmega');
    
    if sum(checkBarOmega) && ~sum(checkMaxOmega)
        options.MaxOmega='';
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
end

% Default values for the optional parameters are set inside structure
% 'options'

rcMax  = [0 0];


MaxOmega = options.MaxOmega;
BarOmega = options.BarOmega;
StdOmega = options.StdOmega;
tol=options.tol;
% tolmap = tolerance between the requested and empirical misclassification
% probabilities
tolmap    = options.tol(1);
% tolncx2 = tolerance to use in routine ncx2mixtcdf (which computes cdf of
% linear combinations of non central chi2 distributions)
tolncx2  = options.tol(2);
hom      = options.hom;
PiLow    = options.PiLow;
resN     = options.resN;
lim      = options.lim;
restrfactor= options.restrfactor;
Display  = options.Display;
betadistrib =options.betadistrib;
Xdistrib = options.Xdistrib;

% If Xdistrib is a structure then we control whether field Xdistrib.type is
% present. If the answer is yes and it is uniform, normal or halfnormal we
% provide all the default options (that is uniform in [0 1] for all the
% explanatory variables and intercept present). If Xdistrib is not a
% structure then Xdistrib is redefined as a structure with
% Xdistrib.type='Uniform'
% Xdistrib.a=zeros(p,1)
% Xdistrib.a=ones(p,1)
% Xdistrib.intercept=1

if isstruct(Xdistrib)
    fld=fieldnames(Xdistrib);
    if find(strcmp('intercept',fld))
        intercept=Xdistrib.intercept;
    else
        intercept=1;
    end
    
    if ~isempty(find(strcmp('type',fld),1))
        if strcmp(Xdistrib.type,'Uniform')
            if intercept==1
                Xdistribdef.a=zeros(p-1,k);
                Xdistribdef.b=ones(p-1,k);
                Xdistribdef.intercept=1;
            else
                Xdistribdef.a=zeros(p,k);
                Xdistribdef.b=ones(p,k);
                Xdistribdef.intercept=0;
            end
            Xdistribdef.type='Uniform';
            
        elseif strcmp(Xdistrib.type,'Normal')
            if intercept==1
                Xdistribdef.mu=0.5+zeros(p-1,k);
                Xdistribdef.sigma=ones(p-1,k);
                Xdistribdef.intercept=1;
            else
                Xdistribdef.mu=0.5+zeros(p,k);
                Xdistribdef.sigma=ones(p,k);
                Xdistribdef.intercept=0;
            end
            Xdistribdef.type='Normal';
            
        elseif strcmp(Xdistrib.type,'HalfNormal')
            if intercept==1
                Xdistribdef.sigma=ones(p-1,k);
                Xdistribdef.intercept=1;
            else
                Xdistribdef.sigma=ones(p,k);
                Xdistribdef.intercept=0;
            end
            Xdistribdef.type='HalfNormal';
            
        elseif strcmp(Xdistrib.type,'User')
            % If Xdistrib.type is User the user must directly supply matrix
            % of the means of the p explnatory variabels for the k gruops
            % so we check whether field BarX is present
            
            if ~isempty(find(strcmp('BarX',fld),1))
                if intercept ==1
                    BarX= [ones(1,k); Xdistrib.BarX];
                else
                    BarX=Xdistrib.BarX;
                end
            else
                str = sprintf(['\n'...
                    'Distribution is ''User'' but the matrix BarX containing the means\n'...
                    'of the p explanatory has not been given']);
                sprintf('%s',str);
                error('FSDA:MixSimreg:MissingField','Please also supply inside Xdistrib field BarX')
            end
            Xdistribdef.type='User';
            Xdistribdef.BarX=BarX;
            Xdistribdef.intercept=intercept;
            
        else
            error('FSDA:MixSimreg:WrongXdistrib','Possible values for option Xdistrib are ''Normal'' ''Uniform'' ''HalfNormal'' and ''User'' ')
        end
    else
        Xdistribdef.type='Uniform';
        Xdistribdef.a=zeros(p,k);
        Xdistribdef.b=ones(p,k);
        Xdistribdef.intercept=1;
    end
    
    % Check if user options inside options.fground are valid options
    chkoptions(Xdistribdef,fld)
    for i=1:length(fld)
        Xdistribdef.(fld{i})=Xdistrib.(fld{i});
    end
    
    % For the options not set by the user use their default value
    Xdistrib=Xdistribdef;
else
    Xdistrib=struct;
    Xdistrib.type='Uniform';
    Xdistrib.a=zeros(p-1,k);
    Xdistrib.b=ones(p-1,k);
    Xdistrib.intercept=1;
end

% Now construct BarX  (the X design point in which overlap will be computed)
intercept=Xdistrib.intercept;

if find(strcmp('Uniform',Xdistrib.type))
    
    if intercept==1
        BarX=[ones(1,k); (Xdistrib.a+Xdistrib.b)/2];
    else
        BarX= (Xdistrib.a+Xdistrib.b)/2;
    end
    
elseif find(strcmp('Normal',Xdistrib.type))
    
    if intercept==1
        BarX=[ones(1,k); Xdistrib.mu];
    else
        BarX=Xdistrib.mu;
    end
    
elseif find(strcmp('HalfNormal',Xdistrib.type))
    
    if intercept==1
        BarX=[ones(1,k);  Xdistrib.sigma*sqrt(2/pi)];
    else
        BarX=Xdistrib.sigma*sqrt(2/pi);
    end
    
elseif find(strcmp('User',Xdistrib.type))
    % In this case there is no additional computation to do because the user
    % has already supplied matrix BarX
else
    error('FSDA:MixSimreg:WrongXdistrib','Possible values for option Xdistrib are ''Normal'' ''Uniform'' ''HalfNormal'' and ''User'' ')
end

Xdistrib.BarX=BarX;

if isstruct(betadistrib)
    fbetadistrib=fieldnames(betadistrib);
    
    if strcmp('Uniform',betadistrib.type)
        a=find(strcmp('a',fbetadistrib),1);
        if isempty(a)
            betadistrib.a=0;
        end
        
        b=find(strcmp('b',fbetadistrib),1);
        if isempty(b)
            betadistrib.b=0;
        end
        
    elseif strcmp('Normal',betadistrib.type)
        mu=find(strcmp('mu',fbetadistrib), 1);
        if isempty(mu)
            betadistrib.mu=0.5;
        end
        
        sigma=find(strcmp('sigma',fbetadistrib), 1);
        if isempty(sigma)
            betadistrib.sigma=1;
        end
        
    elseif strcmp('HalfNormal',betadistrib.type)
        sigma=find(strcmp('sigma',fbetadistrib), 1);
        if isempty(sigma)
            betadistrib.sigma=sqrt(2/pi) ;
        end
        
    elseif strcmp('User',betadistrib.type)
        d=find(strcmp('Beta',fbetadistrib),1);
        if ~isempty(d)
            Beta= betadistrib.Beta;
        else
            error('FSDA:MixSimreg:MissingField','If betadistrib =''User'' then the user must provide input matrix Beta')
        end
    else
        error('FSDA:MixSimreg:Wrongbetadistrib','Possible values for option betadistrib are ''Normal'' ''Uniform'' ''HalfNormal'' and ''User'' ')
    end
else
    betadistrib=struct;
    betadistrib.type='Normal';
    betadistrib.mu=0;
    betadistrib.sigma=1;
end

if ~islogical(hom)
    error('FSDA:MixSimreg:Wronghom','option hom must be a logical value')
end

if PiLow < 0 || PiLow > 1
    error('FSDA:MixSimreg:WrongPiLow','Option PiLow must be in interval [0 1]')
end


if resN < 1
    error('FSDA:MixSimreg:WrongresN','Number of resimulations cannot be smaller than 1')
end

if (min(options.tol) <= 0)
    error('FSDA:MixSimreg:Wrongtol','Wrong value of tolerance, it must a scalar stricty greater than 0')
end

if lim < 1
    error('FSDA:MixSimreg:Wronglim','Wrong value of lim, it cannot be smaller than 1')
end

if isempty(MaxOmega) && isempty(StdOmega)  && ~isempty(BarOmega)
    % method =0 ==> just BarOmega has been specified
    method = 0;
    Omega = BarOmega;
elseif isempty(BarOmega) && isempty(StdOmega)
    % method =1 ==> just MaxOmega has been specified
    method = 1;
    if isempty(MaxOmega)
        Omega=0.15;
    else
        Omega = MaxOmega;
    end
elseif isempty(BarOmega) && isempty(MaxOmega) && ~isempty(StdOmega)
    % method =1.5 ==> Just StdOmega has been specified
    method = 1.5;
    Omega=StdOmega;
elseif ~isempty(BarOmega) && ~isempty(MaxOmega) && isempty(StdOmega)
    % method =2 ==> both BarOmega and MaxOmega have been specified
    method = 2;
elseif  ~isempty(BarOmega) && ~isempty(StdOmega) && isempty(MaxOmega)
    % method =3 ==> both BarOmega and StdOmega have been specified
    method = 3;
elseif isempty(BarOmega) && isempty(MaxOmega)
    % method =-1 ==> both BarOmega and MaxOmega have not been specified
    method = -1;
elseif ~isempty(BarOmega) && ~isempty(StdOmega) && ~isempty(MaxOmega)
    % method =4 ==> both BarOmega MaxOmega and StdOmega have been specified
    method = 4;
    
else
    method= 2;
end

% Get in vector indabovediag the linear indices of the elements above
% diagonal in a matrix of size k-by-k. This will be necessary to compute the
% standard deviation of overlapping
indabovediag=triu2vec(k,1);


if method == 0 || method == 1 || method == 1.5
    
    Q = OmegaClustReg(Omega, method, k, PiLow, BarX, betadistrib, ...
        tol, lim, resN, hom, restrfactor, Display);
    
elseif method == 2
    
    % In this case both BarOmega and MaxOmega have been specified
    Q = OmegaBarOmegaMaxReg(k, PiLow, BarX, betadistrib, ...
        tol, lim, resN, hom, BarOmega, MaxOmega, restrfactor, Display);
    
elseif method ==3
    % In this case both BarOmega and StdOmega have been specified
    StdOmega=options.StdOmega;
    Q=OmegaBarOmegaStdReg(k, PiLow, BarX, betadistrib, ...
        tol, lim, resN, hom, BarOmega, StdOmega, restrfactor, Display);
elseif method ==4
    % In this case both MaxOmega, BarOmega and StdOmega have been specified
    error('FSDA:MixSimreg:TooManyConstr','It is not possible to specify both MaxOmega, BarOmega and StdOmega at the same time')
    
elseif method~=-1
    error('FSDA:MixSimreg:WrongMethod','Should never enter here')
else
    % isempty(BarOmega) && isempty(MaxOmega)
    error('FSDA:MixSimreg:ConstrRequired','At least one overlap characteristic between MaxOmega and BarOmega should be specified')
end

% Add details about the distribution of X which has been used
Q.Xdistrib=Xdistrib;
% Add details about the distribution of beta which has been used
Q.betadistrib=betadistrib;
% Given that the variances are scalars, it is unnecessary to store a
% 1-by-1-by-k 3D vector
Q.S=squeeze(Q.S);
out = Q;


%% Beginning of inner functions

    function  Q = OmegaClustReg(Omega, method, k, PiLow, BarX, betadistrib,  ...
            tol, lim, resN, hom, restrfactor, Display)
        % OmegaClustReg procedure when average or maximum overlap or Std of
        % overlap is specified (not more than one overlapping measure)
        %
        %  INPUT parameters
        %
        % Omega     : scalar containing requested overlap value
        % method    : scalar which specifies whether average or maximum
        %             overlap is requested. If method == 0 average overlap
        %             is requested elseif method == 1 max overlap is
        %             required elseif method == 1.5 std of overlap is
        %             requested.
        % k         : scalar, number of components (groups)
        % PiLow     : smallest mixing proportion allowed
        % BarX      : p-by-k matrix containing the p dimensional centroids
        %             of the k groups
        % betadistrib : structure which specifies the distribution to use
        %               for each element of the vectors of regression coefficients.
        %               Once chosen, the distribution together with its parameters
        %               is fixed for each element of beta, across each group.
        %                 The following options are admitted for betadistrib:
        %               NORMAL DISTRIBUTION N(mu, sigma)
        %                   betadistrib.type='Normal';
        %                   betadistrib.mu = scalar, containing parameter mu for the
        %                       distribution of each element of beta across each
        %                       group. The default value of betadistrib.mu is 0
        %                   betadistrib.sigma = scalar, containing parameter sigma for
        %                       the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   UNIFORM DISTRIBUTION U(a, b)
        %                   > betadistrib.type='Uniform';
        %                   > betadistrib.a = scalar, containing parameter a for the
        %                     distribution of each element of beta across each
        %                     group. The default value of betadistrib.a is 0
        %                   > betadistrib.b = scalar, containing parameter b for
        %                     the distribution of each element of beta across
        %                     each group. The default value of betadistrib.b is 1.
        %                   HALF NORMAL DISTRIBUTION Half(sigma)= |N(0 sigma)|
        %                   > betadistrib.type='HalfNormal';
        %                   > betadistrib.sigma = scalar containing parameter sigma
        %                       for the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   USER DISTRIBUTION
        %                   > betadistribtion.type='User';
        %                     If betadistribtion.type='User' the user must directly
        %                     provide the values of the beta coefficients.
        %                     So, if betadistribtion.type is 'User', we expect there
        %                     is a field called Beta.
        %                   > betadistribution.Beta = matrix of size (p-1)-by k
        %                     (if intercept is present) or p-by-k (if intercept is
        %                     not present) containing the vectors of regression
        %                     coefficients for the k groups.
        % tol, lim  : parameters for ncx2mixtcdf.m which computes the probability
        %             of overlapping
        % resN      : scalar, number of resamplings allowed
        % hom       : scalar. If hom =1 equal variances are
        %             generated
        %
        % Optional input parameters
        %
        %restrfactor: scalar in the interval [1 \infty] which specifies the
        %             maximum ratio to allow between the largest eigenvalue and
        %             the smallest eigenvalue of the k covariance matrices which
        %             are generated.
        %    Display: Level of display.
        %             'off' displays no output;
        %             'notify' (default) displays output if requested
        %             overlap cannot be reached in a particular simulation
        %             'iter' displays output at each iteration of each
        %             simulation
        %
        %
        %  OUTPUT parameters
        %
        %   A structure Q containing the following fields
        %           Pi   : vector of length k containing mixing proportions
        %           Mu   : matrix of size k-by-v contains mean vectors of
        %                  the k groups
        %           S    : array of size v-by-v-by-k containing covariance matrices
        %       OmegaMap : k-by-k matrix containing misclassification probabilities
        %       BarOmega : scalar, average overlap of the groups which have been
        %                  generated
        %       MaxOmega : scalar, maximum overlap of the groups which have been
        %                  generated
        %       StdOmega : scalar, standard deviation of overlap of the groups which have been
        %                  generated
        %        rcMax   : vector of length 2 containing the pair of
        %                  components producing the highest overlap
        %           fail : flag indicating if the process failed (1). If everything went
        %                  well fail=0
        %
        
        
        if nargin< 13
            prnt = 1;
        else
            switch Display
                case {'none','off'}
                    prnt = 0;
                case {'notify','notify-detailed'}
                    prnt = 1;
                case {'iter','iter-detailed'}
                    prnt=2;
                otherwise
                    prnt = 1;
            end
        end
        
        fixcl=zeros(k,1);
        
        for isamp=1:resN
            
            fail=0;
            
            % Generate parameters
            % procedure genPi generates (mixture proportions) k numbers
            % whose sum is equal to 1 and the smallest value is not smaller
            % than PiLow
            Pigen=genPi(k,PiLow);
            % procedure genMu generates random centroids
            [Mugen, Beta]=genMu(k, BarX, betadistrib);
            % Mugen = k-by-1
            % Beta = p-by-k
            
            % Generate the covariance matrices
            if hom == 0
                Sgen= rand(k,1);
            else % homogeneous clusters
                Sgen= rand(1,1)*ones(k,1);
            end
            
            
            if nargin>11 && ~isempty(restrfactor)
                %S05  = Sgen;
                %Sinv = Sgen;
                detS      = zeros(k,1);
                Lambda_vk = Sgen';
                
                % Apply the restrictions to matrix Lambda_vk
                autovalues=restreigen(Lambda_vk,Pigen,restrfactor);
                
                Sgen  = zeros(1,1,k);
                Sinv  = Sgen;
                S05   = Sgen;
                for j=1:k
                    %disp(j)
                    Sgen(1,1,j) = autovalues(j);
                    S05(1,1,j)  = sqrt(Sgen(1,1,j));
                    Sinv(1,1,j) = Sgen(1,1,j)^(-1);
                    detS(j)     = autovalues(j);       
                end
                Sgen = reshape(Sgen,1,1,k);
                [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
            else
                Sgen=reshape(Sgen,1,1,k);
                [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen);
            end
            
            
            % check if desired overlap is reachable
            asympt = 1;
            c = 0.0;
            
            [OmegaMap, Balpha, Malpha, rcMax] = ...
                GetOmegaMap(c, 1, k, li, di, const1, fixcl, tolncx2, lim, asympt);
            
            if method == 0
                diff = Balpha - Omega;
                if prnt>1
                    disp(['Average empirical overlap - Average requested overlap=' num2str(diff)])
                end
            else
                diff = Malpha - Omega;
                if prnt>1
                    disp(['Max empirical overlap - Max requested overlap=' num2str(diff)])
                end
            end
            
            
            
            if (diff < -tolmap) % Prefixed overlapping is not reachable
                if prnt>=1
                    disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                end
                fail = 1;
            else
                lower=0;
                upper=4;
                % c is a constant which is used to scale the covariance
                % matrices of the groups in order to obtain a prespecified
                % level of average or maximum overlap
                c=0;
                while c==0
                    
                    [c, OmegaMap, Balpha, Malpha] = FindC(lower, upper, Omega, ...
                        method, k, li, di, const1, fixcl, tol, lim);
                    lower =upper;
                    upper=upper^2;
                    if upper>1e+10 % TOCKECK 100000
                        if prnt>=1
                            disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                        end
                        fail=1;
                        break
                    end
                end
            end
            
            if fail==0
                Sgen=c*Sgen;
                break  % this break enables to get out from the resampling loop
            end
        end
        
        if method == 0
            if prnt>1
                disp(['Average empirical overlap - Average requested overlap=' num2str(Balpha-Omega)])
            end
        else
            if prnt>1
                disp(['Max empirical overlap - Max requested overlap=' num2str(Malpha-Omega)])
            end
        end
        
        
        if isamp == resN
            warning('FSDA:MixSimreg:OverlapNotReached',['The desired overlap has not been reached in ' num2str(resN) ' simulations']);
            warning('FSDA:MixSimreg:NsimulTooSmall','Please increase the number of simulations allowed (option resN) or change the value of overlap');
            
            fail = 1;
        end
        
        % Compute standard deviation of overlap
        cand=triu(OmegaMap,1)+(tril(OmegaMap,-1))';
        overlapv=cand(:);
        % Recompute rcMax (the indices of the two groups producing the
        % highest overlap)
        [~,indmaxoverlapv]=max(overlapv);
        [rcMax(1), rcMax(2)]=ind2sub([k k],indmaxoverlapv);
        
        %         overlapv=overlapv(overlapv>0);
        %         if length(overlapv)<0.5*k*(k-1);
        %             overlapc=[overlapv; zeros(0.5*k*(k-1)-length(overlapv),1)];
        %         else
        %             overlapc=overlapv;
        %         end
        %         % Compute standard deviation of overlap for current
        %         % solution
        %         stdoverlap=std(overlapc);
        stdoverlap=std(overlapv(indabovediag));
        
        
        Q = struct;
        Q.OmegaMap=OmegaMap;
        Q.BarOmega=Balpha;
        Q.MaxOmega=Malpha;
        Q.StdOmega=stdoverlap;
        Q.rcMax=rcMax;
        Q.fail=fail;
        Q.Pi=Pigen;
        Q.Mu=Mugen;
        Q.S=Sgen;
        Q.Beta=Beta;
    end


% In this case both BarOmega and MaxOmega have been specified
    function Q = OmegaBarOmegaMaxReg(k, PiLow, BarX, betadistrib, ...
            tol, lim, resN, hom, BarOmega, MaxOmega, restrfactor, Display)
        % OmegaBarOmegaMaxReg procedure when average and maximum overlap are both specified
        %
        %
        %  INPUT parameters
        %
        %       p  : scalar, dimensionality
        %    PiLow : scalar, smallest mixing proportion allowed
        % BarX     : p-by-k matrix containing the p dimensional centroids
        %             of the k groups
        % betadistrib : structure which specifies the distribution to use
        %               for each element of the vectors of regression coefficients.
        %               Once chosen, the distribution together with its parameters
        %               is fixed for each element of beta, across each group.
        %                 The following options are admitted for betadistrib:
        %               NORMAL DISTRIBUTION N(mu, sigma)
        %                   betadistrib.type='Normal';
        %                   betadistrib.mu = scalar, containing parameter mu for the
        %                       distribution of each element of beta across each
        %                       group. The default value of betadistrib.mu is 0
        %                   betadistrib.sigma = scalar, containing parameter sigma for
        %                       the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   UNIFORM DISTRIBUTION U(a, b)
        %                   > betadistrib.type='Uniform';
        %                   > betadistrib.a = scalar, containing parameter a for the
        %                     distribution of each element of beta across each
        %                     group. The default value of betadistrib.a is 0
        %                   > betadistrib.b = scalar, containing parameter b for
        %                     the distribution of each element of beta across
        %                     each group. The default value of betadistrib.b is 1.
        %                   HALF NORMAL DISTRIBUTION Half(sigma)= |N(0 sigma)|
        %                   > betadistrib.type='HalfNormal';
        %                   > betadistrib.sigma = scalar containing parameter sigma
        %                       for the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   USER DISTRIBUTION
        %                   > betadistribtion.type='User';
        %                     If betadistribtion.type='User' the user must directly
        %                     provide the values of the beta coefficients.
        %                     So, if betadistribtion.type is 'User', we expect there
        %                     is a field called Beta.
        %                   > betadistribution.Beta = matrix of size (p-1)-by k
        %                     (if intercept is present) or p-by-k (if intercept is
        %                     not present) containing the vectors of regression
        %                     coefficients for the k groups.
        % tol, lim : parameters for ncx2mixtcdf function which computes the
        %            probabilities of overlapping
        %     resN : scalar, number of resamplings allowed
        %      hom : scalar. If hom =1 equal variances are
        %            generated.
        % BarOmega : scalar, required average overlap
        % MaxOmega : scalar, required maximum overlap
        %
        % Optional input parameters
        %
        %restrfactor: scalar in the interval [1 \infty] which specifies the
        %             maximum ratio to allow between the largest eigenvalue and
        %             the smallest eigenvalue of the k covariance matrices which
        %             are generated
        %    Display: Level of display.
        %             'off' displays no output;
        %             'notify' (default) displays output if requested
        %             overlap cannot be reached in a particular simulation
        %             'iter' displays output at each iteration of each
        %             simulation
        %
        %  OUTPUT parameters
        %
        %   A structure Q containing the following fields
        %
        %             Pi : mixing proportions
        %           Mu   : matrix of size k-by-v containing mean vectors of
        %                  the k groups
        %           S    : array of size v-by-v-by-k containing covariance matrices
        %       OmegaMap : k-by-k matrix containing misclassification probabilities
        %       BarOmega : scalar, average overlap of the groups which have been
        %                  generated
        %       MaxOmega : scalar, maximum overlap of the groups which have been
        %                  generated
        %       StdOmega : scalar, standard deviation of overlap of the groups
        %                  which have been generated
        %        rcMax   : vector of length 2 containing the pair of
        %                  components producing the highest overlap
        %           fail : flag indicating if the process failed (1). If
        %                  everything went well fail=0
        Balpha=BarOmega;
        Malpha=MaxOmega;
        if Malpha<Balpha || Malpha>Balpha*k*(k-1)/2
            disp('Both conditions should hold:')
            disp('1. MaxOverlap > AverOverlap')
            disp('2.  MaxOverlap < AverOverlap * K (K - 1) / 2')
            error('FSDA:MixSimreg:WrongOverlapSupplied','Incorrect values of average and maximum overlaps...');
            
        else
            
            if nargin< 13
                prnt = 1;
            else
                switch Display
                    case {'none','off'}
                        prnt = 0;
                    case {'notify','notify-detailed'}
                        prnt = 1;
                    case {'iter','iter-detailed'}
                        prnt=2;
                    otherwise
                        prnt = 1;
                end
            end
            li2=zeros(2, 2, 1);
            di2=li2;
            const12=zeros(2);
            
            fix2=zeros(2,1);
            
            for isamp=1:resN
                
                % generate parameters
                % procedure genPi generates (mixture proportions) k numbers
                % whose sum is equal to 1 and the smallest value is not
                % smaller than PiLow
                Pigen=genPi(k,PiLow);
                
                
                % procedure genMu generates random centroids
                [Mugen, Beta]=genMu(k, BarX, betadistrib);
                % Mugen = k-by-1
                % Beta = p-by-k
                
                % Generate the covariance matrices
                if hom == 0
                    Sgen= abs(randn(1))*rand(k,1);
                else % homogeneous clusters
                    Sgen= rand(1,1)*ones(k,1);
                end
                
                
                if nargin>11 && ~isempty(restrfactor)
                    S05=Sgen;
                    Sinv=Sgen;
                    detS=zeros(k,1);
                    Lambda_vk=Sgen';
                    
                    
                    % Apply the restrictions to matrix Lambda_vk
                    autovalues=restreigen(Lambda_vk,Pigen,restrfactor);
                    
                    for j=1:k
                        %disp(j)
                        Sgen = autovalues';
                        S05 = sqrt(Sgen);
                        Sinv(:,:,j) = Sgen.^(-1);
                        detS(j)=autovalues(j);
                        
                    end
                    Sgen=reshape(Sgen,1,1,k);
                    [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
                else
                    Sgen=reshape(Sgen,1,1,k);
                    [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen);
                end
                
                
                
                % Check if maximum overlap is reachable.  Maximum
                % overlapping takes place when c \rightarrow \infty that is
                % when asympt =1
                asympt = 1;
                % If asympt=1 c does not matter so any value of c will do
                c = 0;
                
                % Initialize fail. A flag which indicates
                % whether the required overlapped has been
                % reached. fail=1 means required overlap
                % not reached
                fail=1;
                
                % fixcl = vector which specifies what are the clusters which
                % participate to the process of inflation. In this stage
                % all clusters partecipate to the process of inflation, so
                % fixcl is a vector of zeroes
                fixcl=zeros(k,1);
                
                [OmegaMap, Balpha, Malpha, rcMax]=GetOmegaMap(c, 1, k, li, di, const1, fixcl, tolncx2, lim, asympt);
                
                % Malpha is the maximum level of overlapping which can be
                % reached with asympt=1
                diff = Malpha - MaxOmega;
                
                
                if diff >= -tolmap %  maximum level of overlapping is reachable
                    
                    lower = 0.0;
                    upper = 2^10;
                    
                    while fail ~=0
                        
                        % find constant c for two clusters which show the
                        % highest ovelapping
                        
                        % Extract parameters for the two clusters with the
                        % highest overlapping
                        li2(1,2,:) = li(rcMax(1),rcMax(2),:);
                        di2(1,2,:) = di(rcMax(1),rcMax(2),:);
                        const12(1,2)=const1(rcMax(1),rcMax(2));
                        li2(2,1,:) = li(rcMax(2),rcMax(1),:);
                        di2(2,1,:) = di(rcMax(2),rcMax(1),:);
                        const12(2,1)=const1(rcMax(2),rcMax(1));
                        
                        Malpha = MaxOmega;
                        
                        % the fourth input element of FindC is method (in
                        % this case 1 is supplied because maximum overlap
                        % is requested)
                        % The fifth input element of findC is k. In this
                        % case k=2 (the two groups which show the highest
                        % overlap)
                        c=FindC(lower, upper, Malpha, 1, 2, li2, di2, const12, fix2, tol, lim);
                        
                        if c == 0 % abnormal termination
                            if prnt>=1
                                disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                            end
                            fail = 1;
                            break
                        end
                        
                        asympt = 0;
                        % Compute map of misclassification probabilities
                        % (OmegaMap)
                        % using the value of c which comes out from procedure
                        % findC
                        % Balpha = Average overlap using c from findC
                        % Malpha = Maximum overlap using c from findC
                        [OmegaMap, Balpha, Malpha, rcMax]=GetOmegaMap(c, 1, k, li, di, const1, fixcl, tolncx2, lim, asympt);
                        upper = c;
                        
                        % We hope that Balpha (overall average overlapping)
                        % obtained using c associated with the two clusters
                        % which showed the highest overlapping is greater
                        % than BarOmega (average requested overlapping).
                        diff = Balpha - BarOmega;
                        
                        % If diff<tolmap the desired average overlap characteristic is
                        % possibly unattainable using the candidate \mu
                        % \Sigma and \pi and a new candidate is needed
                        if (diff < -tolmap) % BarOmega is not reachable
                            if prnt>=1
                                disp(['Warning: the desired average overlap cannot be reached in simulation '  num2str(isamp)]);
                            end
                            fail = 1;
                            break
                        end
                        
                        % Now we make sure that none of pairwise overlaps
                        % (that is  make sure that the maximum pairwise
                        % overlap (which is Malpha) does not exceed MaxOmega
                        % (the maximum requested overlap). If this is the
                        % case  do another iteration of the loop (while
                        % fail ~=0) using rcMax which has just been found.
                        
                        
                        
                        diff = Malpha - MaxOmega;
                        if prnt>1
                            disp(['Average empirical overlap - Average requested overlap=' num2str(Balpha - BarOmega)])
                            disp(['Max empirical overlap - Max requested overlap=' num2str(diff)])
                        end
                        
                        if (diff < tolmap) %  MaxOmega has been reached
                            fail = 0;
                            break
                        end
                        
                    end
                    
                end
                
                if fail == 0
                    %  OmegaMax is reached and OmegaBar is reachable
                    %  correct covariances by multiplier c
                    
                    if nargin>11 && ~isempty(restrfactor)
                        S05=(c^0.5)*S05;
                        Sinv=(1/c)*Sinv;
                        detS=c*detS;
                        [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
                    else
                        Sgen=c*Sgen;
                        [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen);
                    end
                    
                    % The two clusters which enabled to obtain the highest
                    % overlap are kept unchanged all the way through the
                    % termination of the algorithm
                    fixcl(rcMax(1)) = 1;
                    fixcl(rcMax(2)) = 1;
                    upper = 1;
                    
                    % Now find the c which guarrantees the average
                    % requested ovelapping BarOmega
                    % method =0 because average overlapping is requested
                    method = 0;
                    Malphain=Malpha;
                    [c, OmegaMap, Balpha, Malpha, rcMax]=FindC(lower, upper, BarOmega, method, k, li, di, const1, fixcl, tol, lim);
                    
                    % If c =0 max number of iterations has been reached
                    % inside findc therefore another simulation is
                    % requested
                    if c==0 || abs(Malphain-Malpha)>tolmap
                        if prnt>=1
                            disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                        end
                        fail=1;
                    else
                        % correct covariances by multiplier c
                        for jj=1:k
                            if fixcl(jj) == 0
                                Sgen(:,:,jj) = c * Sgen(:,:,jj);
                            end
                        end
                        
                        
                        % Inside if  fail==0 OmegaMax is reached and OmegaBar is reachable
                        break
                    end
                    
                end
                
            end
            
            if prnt>1
                disp(['Average empirical overlap - Average requested overlap=' num2str(Balpha - BarOmega)])
                disp(['Max empirical overlap - Max requested overlap=' num2str(Malpha - MaxOmega)])
            end
            
            if isamp == resN && resN>1
                warning('FSDA:MixSimreg:OverlapNotReached',['The desired overlap has not been reached in ' num2str(resN) ' simulations']);
                warning('FSDA:MixSimreg:NsimulTooSmall','Increase the number of simulations allowed (option resN) or change the value of overlap');
                
                fail = 1;
                
            end
            
            % Compute standard deviation of overlap
            cand=triu(OmegaMap,1)+(tril(OmegaMap,-1))';
            overlapv=cand(:);
            %             overlapv=overlapv(overlapv>0);
            %             if length(overlapv)<0.5*k*(k-1);
            %                 overlapc=[overlapv; zeros(0.5*k*(k-1)-length(overlapv),1)];
            %             else
            %                 overlapc=overlapv;
            %             end
            %             % Compute standard deviation of overlap for current
            %             % solution
            %             stdoverlap=std(overlapc);
            %
            stdoverlap=std(overlapv(indabovediag));
            
            Q=struct;
            Q.OmegaMap=OmegaMap;
            Q.BarOmega=Balpha;
            Q.MaxOmega=Malpha;
            Q.StdOmega=stdoverlap;
            Q.rcMax=rcMax;
            Q.fail=fail;
            Q.Pi=Pigen;
            Q.Beta=Beta;
            Q.Mu=Mugen;
            Q.S=Sgen;
            
        end
    end


    function Q = OmegaBarOmegaStdReg(k, PiLow, BarX, betadistrib, ...
            tol, lim, resN, hom, BarOmega, StdOmega, restrfactor, Display)
        % OmegaBarOmegaStdReg procedure when average and std of overlap are both specified
        %
        %
        %  INPUT parameters
        %
        %       k  : scalar, number of components
        %    PiLow : scalar, smallest mixing proportion allowed
        % BarX     : p-by-k matrix containing the p dimensional centroids
        %             of the k groups
        % betadistrib : structure which specifies the distribution to use
        %               for each element of the vectors of regression coefficients.
        %               Once chosen, the distribution together with its parameters
        %               is fixed for each element of beta, across each group.
        %                 The following options are admitted for betadistrib:
        %               NORMAL DISTRIBUTION N(mu, sigma)
        %                   betadistrib.type='Normal';
        %                   betadistrib.mu = scalar, containing parameter mu for the
        %                       distribution of each element of beta across each
        %                       group. The default value of betadistrib.mu is 0
        %                   betadistrib.sigma = scalar, containing parameter sigma for
        %                       the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   UNIFORM DISTRIBUTION U(a, b)
        %                   > betadistrib.type='Uniform';
        %                   > betadistrib.a = scalar, containing parameter a for the
        %                     distribution of each element of beta across each
        %                     group. The default value of betadistrib.a is 0
        %                   > betadistrib.b = scalar, containing parameter b for
        %                     the distribution of each element of beta across
        %                     each group. The default value of betadistrib.b is 1.
        %                   HALF NORMAL DISTRIBUTION Half(sigma)= |N(0 sigma)|
        %                   > betadistrib.type='HalfNormal';
        %                   > betadistrib.sigma = scalar containing parameter sigma
        %                       for the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   USER DISTRIBUTION
        %                   > betadistribtion.type='User';
        %                     If betadistribtion.type='User' the user must directly
        %                     provide the values of the beta coefficients.
        %                     So, if betadistribtion.type is 'User', we expect there
        %                     is a field called Beta.
        %                   > betadistribution.Beta = matrix of size (p-1)-by k
        %                     (if intercept is present) or p-by-k (if intercept is
        %                     not present) containing the vectors of regression
        %                     coefficients for the k groups.
        % tol, lim : parameters for ncx2mixtcdf function which computes the
        %            probabilities of overlapping
        %     resN : scalar, number of resamplings allowed
        %      hom : scalar. If hom =1 equal variances are
        %            generated.
        % BarOmega : scalar, required average overlap
        % StdOmega : scalar, required standard deviation of overlap
        %
        % Optional input parameters
        %
        %restrfactor: scalar in the interval [1 \infty] which specifies the
        %             maximum ratio to allow between the largest eigenvalue and
        %             the smallest eigenvalue of the k covariance matrices which
        %             are generated
        %    Display: Level of display.
        %             'off' displays no output;
        %             'notify' (default) displays output if requested
        %             overlap cannot be reached in a particular simulation
        %             'iter' displays output at each iteration of each
        %             simulation
        %
        %  OUTPUT parameters
        %
        %   A structure Q containing the following fields
        %
        %       OmegaMap : k-by-k matrix containing misclassification probabilities
        %       BarOmega : scalar, average overlap of the groups which have been
        %                  generated
        %       MaxOmega : scalar, maximum overlap of the groups which have been
        %                  generated
        %        rcMax   : vector of length 2 containing the pair of
        %                  components producing the highest overlap
        %           fail : flag indicating if the process failed (1). If
        %                  everything went well fail=0
        %             Pi : mixing proportions
        %           Mu   : matrix of size k-by-v containing mean vectors of
        %                  the k groups
        %           S    : array of size v-by-v-by-k containing covariance matrices
        
        if k<=2
            error('FSDA:MixSimreg:WrongOverlapSupplied','Average and std of overlap can be both set when k>2');
        end
        
        if nargin< 13
            prnt = 1;
        else
            switch Display
                case {'none','off'}
                    prnt = 0;
                case {'notify','notify-detailed'}
                    prnt = 1;
                case {'iter','iter-detailed'}
                    prnt=2;
                otherwise
                    prnt = 1;
            end
        end
        
        li2=zeros(2, 2, 1);
        di2=li2;
        const12=zeros(2);
        
        fix2=zeros(2,1);
        
        MaxOmegaloopini=BarOmega*1.1;
        
        eps=tolmap*10;
        stdoverlap=NaN;
        
        for isamp=1:resN
            
            % generate parameters
            % procedure genPi generates (mixture proportions) k numbers
            % whose sum is equal to 1 and the smallest value is not
            % smaller than PiLow
            Pigen=genPi(k,PiLow);
            
            % procedure genMu generates random centroids
            [Mugen, Beta]=genMu(k, BarX, betadistrib);
            % Mugen = k-by-1
            % Beta = p-by-k
            
            % Generate the covariance matrices
            if hom == 0
                Sgen= abs(randn(1))*rand(k,1);
            else % homogeneous clusters
                Sgen= rand(1,1)*ones(k,1);
            end
            
            
            if nargin>11 && ~isempty(restrfactor)
                S05=Sgen;
                Sinv=Sgen;
                detS=zeros(k,1);
                Lambda_vk=Sgen';
                
                
                % Apply the restrictions to matrix Lambda_vk
                autovalues=restreigen(Lambda_vk,Pigen,restrfactor);
                
                for j=1:k
                    %disp(j)
                    Sgen = autovalues';
                    S05 = sqrt(Sgen);
                    Sinv(:,:,j) = Sgen.^(-1);
                    detS(j)=autovalues(j);
                    
                end
                Sgen=reshape(Sgen,1,1,k);
                [liini, diini, const1ini]=ComputePars(1, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
            else
                Sgen=reshape(Sgen,1,1,k);
                [liini, diini, const1ini]=ComputePars(1, k, Pigen, Mugen, Sgen);
            end
            
            
            % Check if maximum overlap is reachable.  Maximum
            % overlapping takes place when c \rightarrow \infty that is
            % when asympt =1
            asympt = 1;
            % If asympt=1 c does not matter so any value of c will do
            c = 0;
            
            % fixcl = vector which specifies what are the clusters which
            % participate to the process of inflation. In this stage
            % all clusters partecipate to the process of inflation, so
            % fixcl is a vector of zeroes
            fixclini=zeros(k,1);
            
            [OmegaMap, Balphaini, Malphaini, rcMaxini]=GetOmegaMap(c, 1, k, liini, diini, const1ini, fixclini, tolncx2, lim, asympt);
            
            % Malpha is the maximum level of overlapping which can be
            % reached with asympt=1
            
            
            % Check if sigmamax is reachable:
            % sigmamax= sqrt((M-xmin)*(xmax-M))
            % Note that here xmin=0
            sigmamax=sqrt(BarOmega*(Malphaini-BarOmega));
            
            
            diff = sigmamax - StdOmega;
            
            % Initialize fail
            fail=1;
            
            
            if diff < -tolmap % Requested StdOmega of overlapping is reachable
                disp('Requested sigma of ovelap must be smaller than')
                disp('BarOmega*(MaxOmega(achievable) - BarOmega)')
                disp(['In simulation ' num2str(isamp)])
                disp(['MaxOmega(achievable)=' num2str(Malphaini) ' and'])
                disp(['Maximum achievable sigma is=' num2str(sigmamax)]);
                Balpha=Balphaini;
                Malpha=Malphaini;
                
            else
                
                % Now we loop for different values of MaxOmega in order to
                % find the one which guarrantees the required Std of
                % overlap
                if nargin>11 && ~isempty(restrfactor)
                    S05ini=S05;
                    Sinvini=Sinv;
                    detSini=detS;
                end
                Sgenini=Sgen;
                
                Erho1=10;
                step=0.03;
                step05=0;
                
                MaxOmegaloop=MaxOmegaloopini;
                iter=0;
                while abs(Erho1-1)>eps
                    
                    if step<1e-15
                        break
                    end
                    % if the value of MaxOmegaloop is greater than than the max
                    % overlap achivable (which is Malphaini) than requested
                    % standard deviation is too large and it is necessary to
                    % decrease it
                    if MaxOmegaloop> Malphaini
                        fail=1;
                        if prnt>=1
                            disp('Please decrease requested std of overlap')
                        end
                        break
                    end
                    
                    fail=1;
                    Sgen=Sgenini;
                    if nargin>13 && ~isempty(restrfactor)
                        S05=S05ini;
                        Sinv=Sinvini;
                        detS=detSini;
                    end
                    
                    Malpha=Malphaini;
                    Balpha=Balphaini;
                    % rcMaxini=[1;2];
                    rcMax=rcMaxini;
                    li=liini;
                    di=diini;
                    const1=const1ini;
                    fixcl=fixclini;
                    
                    lower = 0.0;
                    upper = 2^10;
                    
                    while fail ~=0
                        
                        % Extract parameters for the two clusters with the
                        % highest overlapping
                        li2(1,2,:) = li(rcMax(1),rcMax(2),:);
                        di2(1,2,:) = di(rcMax(1),rcMax(2),:);
                        const12(1,2)=const1(rcMax(1),rcMax(2));
                        li2(2,1,:) = li(rcMax(2),rcMax(1),:);
                        di2(2,1,:) = di(rcMax(2),rcMax(1),:);
                        const12(2,1)=const1(rcMax(2),rcMax(1));
                        
                        Malpha = MaxOmegaloop;
                        
                        % The fourth input element of FindC is method (in
                        % this case 1 is supplied because maximum overlap
                        % is requested)
                        % The fith input element of findC is k. In this
                        % case k=2 (the two groups which show the highest
                        % overlap)
                        c=FindC(lower, upper, Malpha, 1, 2, li2, di2, const12, fix2, tol, lim);
                        
                        if c == 0 && prnt >=1 % abnormal termination
                            disp(['Warning: the desired overlap cannot be reached in simulation '  num2str(isamp)]);
                            fail = 1;
                            break
                        end
                        
                        asympt = 0;
                        % Compute map of misclassification probabilities
                        % (OmegaMap)
                        % using the value of c which comes out from procedure
                        % findC
                        % Balpha = Average overlap using c from findC
                        % Malpha = Maximum overlap using c from findC
                        [OmegaMap, Balpha, Malpha, rcMax]=GetOmegaMap(c, 1, k, li, di, const1, fixcl, tolncx2, lim, asympt);
                        upper = c;
                        
                        % We hope that Balpha (overall average overlapping)
                        % obtained using c associated with the two clusters
                        % which showed the highest overlapping is greater
                        % than BarOmega (average requested overlapping).
                        diff = Balpha - BarOmega;
                        % If diff<tolmap the desired average overlap characteristic is
                        % possibly unattainable using the candidate \mu
                        % \Sigma and \pi and a new candidate is needed
                        if (diff < -tolmap) % BarOmega is not reachable
                            if Erho1<1 && prnt>=1
                                disp(['Warning: sigma is too small in simulation '  num2str(isamp)]);
                            end
                            fail = 1;
                            break
                        end
                        
                        % Now we make sure that none of pairwise overlaps
                        % (that is  make sure that the maximum pairwise
                        % overlap (which is Malpha) does not exceed MaxOmega
                        % (the maximum requested overlap). If this is the
                        % case  do another iteration of the loop (while
                        % fail ~=0) using rcMax which has just been found.
                        
                        % TO CHECK REMOVED
                        % diff = Malpha - MaxOmegaloop;
                        %if (diff < tolmap) %  MaxOmega has been reached
                        fail = 0;
                        %  break
                        % end
                        
                    end
                    
                    
                    
                    if fail == 0
                        %  OmegaMax is reached and OmegaBar is reachable
                        %  correct covariances by multiplier C
                        
                        if nargin>11 && ~isempty(restrfactor)
                            S05=(c^0.5)*S05;
                            Sinv=(1/c)*Sinv;
                            detS=c*detS;
                            [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen, S05, Sinv, detS);
                        else
                            Sgen=c*Sgen;
                            [li, di, const1]=ComputePars(1, k, Pigen, Mugen, Sgen);
                        end
                        
                        % The two clusters which enabled to obtain the highest
                        % overlap are kept unchanged all the way through the
                        % termination of the algorithm
                        fixcl(rcMax(1)) = 1;
                        fixcl(rcMax(2)) = 1;
                        upper = 1;
                        
                        % Now find the c which guarrantees the average
                        % requested ovelapping BarOmega
                        % method =0 because average overlapping is requested
                        method = 0;
                        % Malphain=Malpha;
                        % FindC(lower, upper, Balpha, method, p, K, li, di, const1, fix, pars, lim, &c, OmegaMap, &Balpha, &Malpha, rcMax);
                        [c, OmegaMap, Balpha, Malpha, rcMax]=FindC(lower, upper, BarOmega, method, k, li, di, const1, fixcl, tol, lim);
                        
                        % If c =0 max number of iterations has been reached
                        % inside Findc therefore another simulation is
                        % requested
                        if c==0 % || abs(Malphain-Malpha)>1*tolmap
                            if Erho1<10 &&  prnt >=1
                                disp(['Warning: sigma is too large in simulation '  num2str(isamp)]);
                            end
                            fail=1;
                            break
                        else
                            % correct covariances by multiplier c
                            for jj=1:k
                                if fixcl(jj) == 0
                                    Sgen(:,:,jj) = c * Sgen(:,:,jj);
                                end
                            end
                            
                        end
                        
                    end
                    
                    % Compute standard deviation of overlap
                    cand=triu(OmegaMap,1)+(tril(OmegaMap,-1))';
                    overlapv=cand(:);
                    %                     overlapv=overlapv(overlapv>0);
                    %                     if length(overlapv)<combk2;
                    %                         overlapc=[overlapv; zeros(combk2-length(overlapv),1)];
                    %                     else
                    %                         overlapc=overlapv;
                    %                     end
                    %                     % Compute standard deviation of overlap for current
                    %                     % solution
                    %                     stdoverlap=std(overlapc);
                    stdoverlap=std(overlapv(indabovediag));
                    
                    Erho1old=Erho1;
                    Erho1=StdOmega/stdoverlap;
                    if prnt == 2
                        disp(['Iteration ' num2str(iter) ' in simulation ' num2str(isamp)])
                        disp(['Average overlap =' num2str(Balpha)])
                        disp('Ratio between std of required overlap and std of empirical overlap')
                        disp(Erho1)
                        % disp(step)
                    end
                    
                    iter=iter+1;
                    if step05==1 || (Erho1old>1 && Erho1<1)
                        step=step*0.5;
                        step05=1;
                    else
                        step=step*0.95;
                    end
                    
                    if Erho1>1
                        MaxOmegaloop=MaxOmegaloop+step;
                    else
                        
                        if MaxOmegaloop-step<=BarOmega && iter >=5
                            if prnt >=1
                                disp(['Warning: sigma is too small in simulation '  num2str(isamp)]);
                            end
                            fail=1;
                            break
                        end
                        %MaxOmegaloop=max(MaxOmegaloop-step,MaxOmegaloopini);
                        MaxOmegaloop=max(MaxOmegaloop-step,BarOmega);
                        
                    end
                end
                if fail==0
                    break
                end
            end
        end
        if isamp == resN && resN>1
            warning('FSDA:MixSimreg:OverlapNotReached',['The desired overlap has not been reached in ' num2str(resN) ' simulations']);
            warning('FSDA:MixSimreg:NsimulTooSmall','Increase the number of simulations allowed (option resN) or change the value of overlap');
            
            fail = 1;
            
        end
        
        
        
        Q=struct;
        Q.OmegaMap=OmegaMap;
        Q.BarOmega=Balpha;
        Q.MaxOmega=Malpha;
        Q.StdOmega=stdoverlap;
        Q.rcMax=rcMax;
        Q.fail=fail;
        Q.Pi=Pigen;
        Q.Mu=Mugen;
        Q.Beta=Beta;
        Q.S=Sgen;
        
    end

    function Pigen=genPi(k,PiLow)
        % genPi generates vector of mixing proportions
        %
        %  Required input arguments:
        % 		k       : scalar, number of components
        % 		PiLow   : scalar, smallest possible mixing proportion
        %
        %  OUTPUT parameters
        %
        % 		Pigen : vector of length k containing mixing proportions
        %               Vector Pigen satisfies the following constraints:
        %               sum(Pigen)=1
        %               min(Pigen)>=PiLow
        %
        %  For example Pigen=genPi(4,0.24)
        %  produces Pigen=
        %    0.2440
        %    0.2517
        %    0.2565
        %    0.2478
        %
        flag = 0;
        
        if PiLow >= 1 || PiLow <= 0
            if PiLow < 0 || PiLow >= 1
                disp('Warning: PiLow is out of range... generated equal mixing proportions...');
            end
            Pigen=zeros(k,1);
            Pigen=Pigen+1/k;
            
        else
            
            Pigen = randg(1,k,1);
            
            s=sum(Pigen);
            
            for j=1:k
                Pigen(j) = PiLow + Pigen(j) / s * (1 - k * PiLow);
                if (Pigen(j) < PiLow)
                    flag = 1;
                    break
                end
            end
            if (flag == 1)
                warning('FSDA:MixSimreg:WrongPiLow','PiLow is too high... generated equal mixing proportions...');
                Pigen=zeros(k,1)+1/k;
            end
        end
    end



    function [Mugen, Beta]=genMu(k, BarX, betadistrib)
        % genMu generates vector of length k containing the k centroids or matrix of means of size k-by-p
        %
        %  Required input arguments:
        %           k : number of components
        %        BarX : p-by-k matrix containing the means of the
        %               p explanatory variables for the k regression
        %               hyperplanes
        % betadistrib : structure which specifies the distribution to use
        %               for each element of the vectors of regression coefficients.
        %               Once chosen, the distribution together with its parameters
        %               is fixed for each element of beta, across each group.
        %                 The following options are admitted for betadistrib:
        %               NORMAL DISTRIBUTION N(mu, sigma)
        %                   betadistrib.type='Normal';
        %                   betadistrib.mu = scalar, containing parameter mu for the
        %                       distribution of each element of beta across each
        %                       group. The default value of betadistrib.mu is 0
        %                   betadistrib.sigma = scalar, containing parameter sigma for
        %                       the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   UNIFORM DISTRIBUTION U(a, b)
        %                   > betadistrib.type='Uniform';
        %                   > betadistrib.a = scalar, containing parameter a for the
        %                     distribution of each element of beta across each
        %                     group. The default value of betadistrib.a is 0
        %                   > betadistrib.b = scalar, containing parameter b for
        %                     the distribution of each element of beta across
        %                     each group. The default value of betadistrib.b is 1.
        %                   HALF NORMAL DISTRIBUTION Half(sigma)= |N(0 sigma)|
        %                   > betadistrib.type='HalfNormal';
        %                   > betadistrib.sigma = scalar containing parameter sigma
        %                       for the distribution of each element of beta across
        %                       each group. The default value of betadistrib.sigma
        %                       is 1
        %                   USER DISTRIBUTION
        %                   > betadistribtion.type='User';
        %                     If betadistribtion.type='User' the user must directly
        %                     provide the values of the beta coefficients.
        %                     So, if betadistribtion.type is 'User', we expect there
        %                     is a field called Beta.
        %                   > betadistribution.Beta = matrix of size (p-1)-by k
        %                     (if intercept is present) or p-by-k (if intercept is
        %                     not present) containing the vectors of regression
        %                     coefficients for the k groups.  %
        %  OUTPUT parameters
        % 		Mugen : vector of length k containing
        %               \overline x_1'*\beta_1, \overline x_2'*\beta_2, ..., \overline x_k'*\beta_k
        %       Beta  : matrix of size p-by-k containing regression
        %               coefficients associated to the k hyperplanes
        
        
        if find(strcmp('Uniform',betadistrib.type))
            %Beta=betadistrib.a+rand(p,k)*(betadistrib.b-betadistrib.a);
            Beta=betadistrib.a+bsxfun(@times,rand(p,k),(betadistrib.b-betadistrib.a));

        elseif find(strcmp('Normal',betadistrib.type))
            %Beta=betadistrib.mu+randn(p,k)*betadistrib.sigma;
            if isscalar(betadistrib.mu)
                betamu = betadistrib.mu;
                %betamu = repmat(betadistrib.mu,p,k); %not necessary because
                % betadistrib.mu is a scalar and is automatically expanded by MATLAB
            elseif size(betadistrib.mu,2)==k
                betamu = repmat(betadistrib.mu,p,1);
            else
                error('FSDA:MixSimreg:Wrongbetadistrib','betadistrib.mu can only be a scalar or a 1 x k vector');
            end
            if isscalar(betadistrib.sigma)
                betasigma = betadistrib.sigma;
                % betasigma = repmat(betadistrib.sigma,p,k); %not necessary because
                % betadistrib.sigma is a scalar and is automatically expanded by MATLAB
            elseif size(betadistrib.sigma,2)==k
                betasigma = repmat(betadistrib.sigma,p,1);
            else
                error('FSDA:MixSimreg:Wrongbetadistrib','betadistrib.sigma can only be a scalar or a 1 x k vector');
            end
            % Beta=betamu+bsxfun(@times,randn(p,k),betasigma); %bsxfun not necessary anymore
            Beta=betamu+randn(p,k).*betasigma;
            
        elseif find(strcmp('HalfNormal',betadistrib.type))
            %Beta=abs(randn(p,k))*betadistrib.sigma;
            Beta=bsxfun(@times,abs(randn(p,k)),betadistrib.sigma);
            
        elseif find(strcmp('User',betadistrib.type))
            Beta=betadistrib.Beta;
        else
            error('FSDA:MixSimreg:Wrongbetadistrib','Possible values for option betadistrib are ''Normal'' ''Uniform'' ''Halfnormal'' and ''User'' ')
        end
        
        Mugen=   sum(Beta.*BarX,1)';
        % Mugen is vector of length k which contains
        % Mugen = \overline x_1'*\beta_1, \overline x_2'*\beta_2, ..., \overline x_k'*\beta_k
    end


    function  [c, OmegaMap2, BarOmega2, MaxOmega2, rcMax]=FindC(lower, upper, Omega, method, k, li, di, const1, fix, tol, lim)
        %find multiplier c to be applied to the variances in the
        %interval [lower upper] in order to reach the required average or
        %maximum overlap
        %
        %  Required input arguments:
        %
        % lower : scalar - lower bound of the interval
        % upper : scalar - upper bound of the interval
        % Omega : scalar, associated with maximum or average overlapping requested
        % method : scalar which specifies whether average (method=0) or maximum
        %          overlap is requested
        %     k  : number of components
        % li, di, const1 : parameters needed for computing overlap,
        %          precalculated using routine ComputePars
        %    fix : vector of length k containing zeros or ones
        %          if fix(j) =1 cluster j does not participate to inflation
        %          or deflation. If fix=zeros(k,1) all clusters participate
        %          in inflation/deflation This parameter is used just if
        %          heterogeneous clusters are used
        %    tol : vector of length 2.
        %          tol(1) (which will be called tolmap) specifies
        %          the tolerance between the requested and empirical
        %          misclassification probabilities (default is 1e-06)
        %          tol(2) (which will be called tolnxc2) specifies the
        %          tolerance to use in routine ncx2mixtcdf (which computes cdf
        %          of linear combinations of non central chi2 distributions).
        %          The default value of tol(2) 1e-06
        %    lim : maximum number of integration terms default is 1e06.
        %          REMARK: Optional parameters tol and lim will be used by
        %          function ncx2mixtcdf.m which computes the cdf of a linear
        %          combination of non central chi2 r.v.. This is the
        %          probability of overlapping
        %
        %  OUTPUT parameters
        %
        %        c    : scalar inflation parameter
        %               c is a constant which is used to scale the covariance
        %               matrices of the groups in order to obtain a
        %               prespecified level of average or maximum overlap
        %   OmegaMap2 : k-by-k matrix containing map of misclassification
        %               probabilities
        %   BarOmega2 : scalar. Average overlap found using c
        %   MaxOmega2 : scalar. Maximum overlap found using c
        %      rcMax  : vector of length 2 containing the indexes associated
        %              to the pair of components producing the highest overlap
        %
        
        diff = Inf;
        stopIter = 200; % 500
        tolmap=tol(1);
        tolncx2=tol(2);
        
        sch = 0;
        asympt = 0;
        
        % Intervals which contain positive or negative powers of 2 are
        % considered. For example first interval is [0 1024] then if
        % MaxOmega2 (maximum overlap which has been found using c=512) is <
        % Omega (maximum required overlap), then the new interval becomes
        % [512 1024]  (c has to be increased and the new candidate c is
        % 0.5*(512+1024)) else the new interval becomes [0 512] (c has to
        % be decreased and the new candidate c is 0.5*(0+512)=256
        while abs(diff) > tolmap
            
            c = (lower + upper) / 2.0;
            
            [OmegaMap2, BarOmega2, MaxOmega2, rcMax]=GetOmegaMap(c, 1, k, li, di, const1, fix, tolncx2, lim, asympt);
            
            if method == 0 % in this case average overlap is requested
                
                if BarOmega2 < Omega
                    % clusters are too far
                    lower = c;
                else
                    upper = c;
                end
                
                diff = BarOmega2 - Omega;
                
            else % in this case maximum overlap is requested
                
                if MaxOmega2 < Omega
                    % clusters are too far
                    lower = c;
                else
                    upper = c;
                end
                
                diff = MaxOmega2 - Omega;
                
            end
            
            sch = sch + 1;
            
            if sch == stopIter
                c = 0.0;
                disp(['Warning: required overlap was not reached in routine findC after ' num2str(stopIter) ' iterations...'])
                break
            end
        end
    end
end
%FScategory:CLUS-MixSim
