function out = mdMAARtest(Y, varargin)
%mdMAARtest Diagnostic tests for Missing Always At Random (MAAR).
%
%<a href="matlab: docsearchFS('mdMAARtest')">Link to the help function</a>
%
%  mdMAARtest implements the three diagnostic procedures proposed by
%  Bojinov, Pillai and Rubin (2020) for investigating whether partially
%  observed variables are associated with missingness indicators after
%  conditioning on the fully observed variables.
%
%  The available procedures are:
%
%  1) 'ccm': comparison of conditional means;
%  2) 'dtmm': direct test of a postulated missingness mechanism;
%  3) 'cop': semiparametric Gaussian-copula diagnostic.
%
%  The null hypothesis is that the restrictions implied by MAAR, row
%  exchangeability and conditional independence of the missingness
%  indicators are compatible with the observed data. Rejection indicates
%  that at least one of these assumptions is violated. The procedures are
%  diagnostic tools and do not establish that a nonrejected mechanism is
%  MAAR.
%
%  Required input arguments:
%
%    Y : Input data. Matrix, table or timetable. n x p data matrix possibly
%        containing NaN values. Rows of Y represent observations and columns
%        represent variables.
%        Data Types - single | double | table | timetable
%
%  Optional input arguments:
%
%  Options common to all methods:
%
%       method : Diagnostic method. Character or string.
%                Possible values are:
%                'ccm'  = comparison of conditional means;
%                'dtmm' = direct test of the missingness mechanism;
%                'cop'  = Gaussian-copula diagnostic.
%                The default value is 'ccm'.
%                Example - 'method','cop'
%                Data Types - char | string
%
%        alpha : Significance level. Scalar in the interval (0,1).
%                Bonferroni corrections are applied internally to the
%                component tests. The default value is 0.05.
%                Example - 'alpha',0.01
%                Data Types - single | double
%
%        ridge : Nonnegative numerical regularization used only when a
%                covariance or information matrix is nearly singular. The
%                default value is 1e-8.
%                Example - 'ridge',1e-6
%                Data Types - single | double
%
%          msg : Display a summary of the results. Boolean. The default value
%                is true.
%                Example - 'msg',false
%                Data Types - logical | single | double
%
%        plots : Produce a method-specific diagnostic plot. Boolean.
%                The default value is false.
%                See More About for details. The default is false.
%                Example - 'plots',true
%                Data Types - logical | single | double
%
%  Options used only by method 'dtmm':
%
%      imputed : Completed data sets used by method 'dtmm'. Cell array or
%                n x p x M numeric array. Each completed data set must have
%                the same dimensions as Y and must contain no missing values.
%                If empty, mdEM and mdImputeStochastic are used to generate
%                the imputations. The default value is [].
%                Example - 'imputed',{Yimp1,Yimp2,Yimp3}
%                Data Types - cell | single | double
%
% nimputations : Number of stochastic imputations for method 'dtmm' when
%                option imputed is empty. Positive integer greater than 1.
%                The default value is 10.
%                Example - 'nimputations',20
%                Data Types - single | double
%
%      maxiter : Maximum number of iterations passed to mdEM when imputations
%                must be generated and to the logistic-regression fitting
%                algorithm . The default value is 100. This option is used
%                just by method 'dtmm'.
%                Example - 'maxiter',200
%                Data Types - single | double
%
%          tol : Convergence tolerance passed to mdEM when imputations must
%                be generated and to the logistic-regression fitting
%                algorithm. The default value is 1e-5.  This option is used
%                just by method 'dtmm'.
%                 Example - 'tol',1e-8
%                Data Types - single | double
%
%  Options used only by method 'cop':
%
%        niterMCMC : Number of MCMC iterations for method 'cop'. The default
%                value is 2000.
%                Example - 'niterMCMC',5000
%                Data Types - single | double
%
%        nburnMCMC : Number of initial MCMC iterations discarded for method
%                'cop'. The default value is 500.
%                Example - 'nburnMCMC',1000
%                Data Types - single | double
%
%         thinMCMC : MCMC thinning interval for method 'cop'. The default value
%                is 1.
%                Example - 'thinMCMC',5
%                Data Types - single | double
%
% pluginthreshold : Threshold used by the Gaussian-copula sampler. Margins
%                having more than pluginthreshold distinct observed values
%                use fixed normal scores, as in sbgcop. The default value is
%                100 (this option is used just by method 'cop').
%                Example - 'pluginthreshold',50
%                Data Types - single | double
%
%  savesamples : Save the posterior conditional-covariance draws produced by
%                method 'cop'. Boolean. The default value is false (this
%                option is used just by method 'cop').
%                Example - 'savesamples',true
%                Data Types - logical | single | double
%
%  Output:
%
%    out : Structure containing the following fields:
%          out.method        = selected diagnostic method.
%          out.reject        = true if at least one component test rejects.
%          out.whichReject   = logical vector associated with the variables
%                              containing missing values. Element k is true
%                              when the missingness indicator of that variable
%                              is implicated by the diagnostic.
%          out.pvalue        = component significance measures.
%
%                              For methods 'ccm' and 'cop', out.pvalue is a
%                              q-by-q matrix, where q is the number of
%                              partially observed variables. Rows correspond
%                              to missingness indicators $R_i$ and columns to
%                              partially observed variables $Y_j$.
%
%                              For method 'ccm', out.pvalue(i,j) is the
%                              p-value of the nested-model F test comparing
%                              the conditional mean of $Y_j$ across the groups
%                              defined by $R_i$, after conditioning on the
%                              fully observed variables. A small value
%                              indicates that knowing whether variable i is
%                              missing provides additional information about
%                              the conditional mean of $Y_j$. Diagonal elements
%                              are not tested.
%
%                              For method 'cop', out.pvalue(i,j) is the
%                              two-sided posterior tail probability for the
%                              latent conditional covariance between
%                              missingness indicator $R_i$ and partially
%                              observed variable $Y_j$, given the fully
%                              observed variables. A small value indicates
%                              that most posterior mass lies on one side of
%                              zero. These quantities are posterior tail
%                              probabilities, not classical frequentist
%                              p-values. Diagonal elements are not tested.
%
%                              For method 'dtmm', out.pvalue is a q-by-1
%                              vector. Element out.pvalue(i) is the pooled
%                              Meng-Rubin $D_L$ p-value for the missingness
%                              mechanism of variable i. It compares a reduced
%                              logistic model for $R_i$ containing only fully
%                              observed variables with a full model that also
%                              contains all imputed partially observed
%                              variables. A small value indicates that at
%                              least one partially observed variable provides
%                              additional information about whether variable
%                              i is missing.
%          out.stat          = method-specific test statistics.
%          out.alpha         = nominal significance level.
%          out.alphaAdjusted = Bonferroni-adjusted component level.
%          out.missingVariables = indices of variables containing NaNs.
%          out.fullyObservedVariables = indices of fully observed variables.
%          out.variableNames = variable names.
%          out.missingIndicator = n x q matrix; 1 denotes a missing value.
%          out.n             = number of observations.
%          out.p             = number of variables.
%          out.nvarmiss      = number q of partially observed variables.
%          out.interpretation = concise interpretation of the result.
%
%          Additional fields for method 'ccm':
%          out.df1, out.df2  = degrees of freedom of the nested-model F tests.
%          out.rejectPairs   = q x q logical matrix of rejected comparisons.
%
%          Additional fields for method 'dtmm':
%          out.df1, out.df2  = D_L reference-distribution degrees of freedom.
%          out.relativeIncrease = relative increase in variance due to
%                              nonresponse in the pooled likelihood-ratio test.
%          out.DLtable       = table with one row for each partially observed
%                              variable and columns dL, rL, DL and pvalue.
%                              Row names are the corresponding variable names.
%          out.deviance      = detailed deviances used in the D_L calculation.
%          out.imputationInfo = information about the completed data sets.
%
%          Additional fields for method 'cop':
%          out.CI            = q x q x 2 array of simultaneous credible limits.
%          out.posteriorMean = posterior means of conditional covariances.
%          out.posteriorMedian = posterior medians.
%          out.rejectPairs   = q x q logical matrix of rejected associations.
%          out.nSaved        = number of retained MCMC draws.
%          out.samples       = posterior draws when savesamples is true.
%
%  More About:
%
% This function implements the three diagnostics from Bojinov, Pillai and
% Rubin:
%
% $\mathbf{ccm}$: comparison of conditional means using nested Gaussian linear models
%       and Bonferroni correction.
%
% $\mathbf{dtmm}$: direct testing of the missingness mechanism through logistic
%       models and multiple imputation. The likelihood-ratio statistics are
%       combined using the Meng-Rubin $D_L$ procedure, often referred to as
%       the $D_3$ procedure in the multiple-imputation literature.
% ​
% $\mathbf{cop}$: gaussian-copula diagnostic based on posterior conditional
%      covariances between partially observed variables and missingness
%       indicators.
%
% ---------------------------------------------------.
%
% $\mathbf{ccm}$ method.
%
%  Let $J_m$ denote the set of variables that contain missing values and
%  $J_f$ the set of fully observed variables, and let
%
%  \[
%      q=|J_m|
%  \]
%
%  denote the number of partially observed variables. Under the assumptions
%  considered by Bojinov, Pillai and Rubin, a missingness indicator $R_k$
%  should not depend on another partially observed variable $Y_j$ after
%  conditioning on $Y_{J_f}$.
%
%  The comparison-of-conditional-means procedure tests this implication by
%  comparing, for every $j$ different from $k$, the nested linear models
%
%    \[
%     Y_j \sim Y_{J_f}
%    \]
%
%  and
%
%    \[
%     Y_j \sim Y_{J_f} + R_k + Y_{J_f}:R_k.
%    \]
%
% $\mathbf{dtmm}$ method.
%
%
%  The direct procedure (dtmm) tests the postulated missingness mechanism
%  separately for each missingness indicator $R_k$. For a given $R_k$, the
%  reduced logistic model contains only an intercept and the fully observed
%  variables $Y_{J_f}$:
%
%  \[
%      \mathrm{logit}\{
%      \Pr(R_k=1\mid Y_{J_f})
%      \}
%      =
%      \alpha_k+Y_{J_f}^{\mathsf T}\beta_k.
%  \]
%
%  The full model additionally contains all the partially observed
%  variables after imputation. Thus, in completed data set $m$,
%
%  \[
%      \mathrm{logit}\{
%      \Pr(R_k=1\mid Y_{J_f},Y_{J_m}^{(m)})
%      \}
%      =
%      \alpha_k^{(m)}
%      +Y_{J_f}^{\mathsf T}\beta_k^{(m)}
%      +(Y_{J_m}^{(m)})^{\mathsf T}\gamma_k^{(m)}.
%  \]
%
%  The response $R_k$ and the fully observed variables do not change across
%  imputations. Consequently, the reduced model is identical in all
%  completed data sets and is fitted only once. The full model, on the
%  other hand, depends on the imputed values and is fitted separately in
%  each of the $M$ completed data sets.
%
%  Let $\widehat\theta_k^{(m)}$ denote the vector of full-model coefficient
%  estimates obtained from completed data set $m$. The Meng-Rubin $D_L$
%  procedure first forms the common pooled coefficient vector
%
%  \[
%      \overline\theta_k
%      =
%      \frac{1}{M}
%      \sum_{m=1}^M
%      \widehat\theta_k^{(m)}.
%  \]
%
%  Each completed data set is then reevaluated at this same pooled
%  coefficient vector, without refitting the model. If
%  $D_{F,m}(\widehat\theta_k^{(m)})$ is the full-model deviance evaluated
%  at its own maximum-likelihood estimate and
%  $D_{F,m}(\overline\theta_k)$ is the deviance evaluated at the pooled
%  coefficients, their average difference
%
%  \[
%      \frac{1}{M}
%      \sum_{m=1}^M
%      \left\{
%      D_{F,m}(\overline\theta_k)
%      -
%      D_{F,m}(\widehat\theta_k^{(m)})
%      \right\}
%  \]
%
%  measures the loss of fit caused by forcing all imputations to use a
%  common full-model parameter vector. It therefore measures the
%  disagreement among the imputation-specific estimates and is used to
%  estimate the additional uncertainty due to missing information.
%
%  Since the full model adds the $q$ partially observed variables, the
%  estimated relative increase in variance is
%
%  \[
%      r_L
%      =
%      \frac{M+1}{q(M-1)}
%      \frac{1}{M}
%      \sum_{m=1}^M
%      \left\{
%      D_{F,m}(\overline\theta_k)
%      -
%      D_{F,m}(\widehat\theta_k^{(m)})
%      \right\}.
%  \]
%
%  Let $D_R$ denote the deviance of the reduced model. The numerator used
%  by this function is
%
%  \[
%      d_L
%      =
%      \max\left\{0,
%      D_R
%      -
%      \frac{1}{M}
%      \sum_{m=1}^M
%      D_{F,m}(\overline\theta_k)
%      \right\}.
%  \]
%
%  The maximum with zero is only a numerical safeguard against tiny negative
%  values caused by finite-precision calculations. The pooled Meng-Rubin
%  statistic used by this function is then
%
%  \[
%      D_L
%      =
%      \frac{d_L}{q(1+r_L)}.
%  \]
%
%
%  Thus $d_L$ measures the improvement of the full model over the reduced
%  model after imposing common pooled coefficients, whereas $r_L$ quantifies
%  the additional uncertainty arising when the imputation-specific
%  full-model estimates disagree substantially.
%
%  The relative increase $r_L$ can also be expressed through the
%  corresponding fraction of missing information,
%
%  \[
%      \lambda_L=\frac{r_L}{1+r_L}.
%  \]
%
%  This transformation maps $r_L\in[0,\infty)$ into
%  $\lambda_L\in[0,1)$. Values close to zero correspond to little
%  missing-information uncertainty, whereas values approaching one
%  correspond to a large missing-information component. The DTMM
%  diagnostic plot therefore uses $\lambda_L$ on the vertical axis, while
%  the original $r_L$ is reported in out.DLtable.
% 
%  The resulting statistic is compared with the
%  finite-$M$ $F$ reference distribution of Meng and Rubin (1992).
%  A small $p$-value indicates that, after conditioning on the fully
%  observed variables, at least one partially observed variable provides
%  additional information about the missingness indicator $R_k$.
%
%
% $\mathbf{cop}$ method.
%
%  The Gaussian-copula procedure uses the extended rank-likelihood sampler of
%  Hoff (2007). For every retained posterior correlation matrix, the function
%  computes the covariance between the partially observed outcomes and their
%  missingness indicators conditional on the fully observed variables. A
%  component is rejected when its Bonferroni-adjusted credible interval does
%  not contain zero.
%
%  The original R implementation uses mice for method 'dtmm'. When option
%  imputed is empty, this MATLAB implementation instead uses the joint-normal
%  EM and stochastic-imputation functions already available in FSDA.
%
%   $\large{\text{PLOTS IN OUTPUT}}$.
%
%   If plots is true, the graphical summary depends on the
%   selected diagnostic method:
%
%   $\mathbf{ccm}$
%       Produces a $q$-by-$q$ heatmap of the pairwise F-test
%       p-values. Rows correspond to missingness indicators
%       $R_i$ and columns to partially observed variables
%       $Y_j$. Diagonal cells are not tested. Small p-values
%       are shown using reddish colours, and cells significant
%       after the Bonferroni correction are displayed in red.
%
%   $\mathbf{dtmm}$
%       Produces a bubble scatter plot that separates the two
%       components entering the Meng-Rubin $D_L$ statistic,
%       \[
%           D_L=\frac{d_L}{q(1+r_L)}.
%       \]
%       The horizontal coordinate $d_L/q$ measures the
%       likelihood-ratio improvement in fit per tested
%       restriction obtained by augmenting the reduced
%       missingness model with the partially observed variables.
%       More precisely, $d_L$ is the reduction in deviance
%       obtained by moving from the reduced model to the full
%       model when the latter is evaluated using the common
%       coefficient vector pooled across imputations.
%       The vertical coordinate is
%       \[
%           \lambda_L=\frac{r_L}{1+r_L},
%       \]
%       a bounded representation of the missing-information
%       component. There is one bubble for each missingness
%       indicator, labelled with the corresponding variable
%       name. Bubble size and colour increase as the associated
%       $p$-value decreases; a black ring identifies components
%       rejected at the Bonferroni-adjusted significance level.
%       The vertical reference line is the Bonferroni-adjusted
%       threshold for $d_L/q$ under $r_L=0$. The horizontal
%       reference line $\lambda_L=0.5$, equivalently $r_L=1$,
%       is an interpretative benchmark for the
%       missing-information component. These reference lines
%       are graphical guides; the formal decision is based on
%       the finite-$M$ Meng-Rubin $p$-value.
%       For graphical readability, numerator values exceeding
%       six times the reference threshold are displayed at that
%       upper bound and marked by a right-pointing symbol. Their
%       exact values remain available in out.DLtable.
%    
%   $\mathbf{cop}$
%       Produces a $q$-by-$q$ heatmap of the two-sided posterior
%       tail probabilities for the latent conditional
%       covariances between missingness indicators $R_i$ and
%       partially observed variables $Y_j$, conditional on the
%       fully observed variables. Rows correspond to $R_i$ and
%       columns to $Y_j$. Diagonal cells are not tested. Small
%       posterior tail probabilities are shown using reddish
%       colours, and rejected components are displayed in red.
%
%   The plot is created in a new MATLAB figure. Setting plots
%   to false suppresses all graphical output but does not
%   affect the numerical results returned in out.
%
%
%
%  See also: mdEM, mdImputeStochastic, mdMCARtest, mdLittleTest, mdJJtest
%
%  References:
%
%  Bojinov, I., Pillai, N. S. and Rubin, D. B. (2020), "Diagnosing missing
%  always at random in multivariate data", Biometrika, Vol. 107, pp. 246-253.
%
%  Hoff, P. D. (2007), "Extending the rank likelihood for semiparametric
%  copula estimation", Annals of Applied Statistics, Vol. 1, pp. 265-283.
%
%  Meng, X.-L. and Rubin, D. B. (1992), "Performing likelihood ratio tests
%  with multiply-imputed data sets", Biometrika, Vol. 79, pp. 103-111.
%
%  Bojinov, I. (2018). diagMAAR: Diagnostic tests for missing always at
%  random [Repository GitHub], https://github.com/bojinov/diagMAAR/
%
%  Copyright 2008-2026.
%  Written by FSDA team
%
%<a href="matlab: docsearchFS('mdMAARtest')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
   %% Example 1: MAAR data and comparison of conditional means.
   rng(1)
   n = 500;
   p = 5;
   Yfull = randn(n,p);
   eta = -1 + Yfull(:,4) - Yfull(:,5);
   prob = 1./(1+exp(-eta));
   Y = Yfull;
   for j = 1:3
       Y(rand(n,1)<prob,j) = NaN;
   end
   out = mdMAARtest(Y,'method','ccm','plots',true);
   disp(out.pvalue)
%}
%
%{
   %% Example 2: Violation of MAAR detected by conditional means.
   rng(2)
   n = 800;
   Yfull = randn(n,5);
   Y = Yfull;
   p1 = 1./(1+exp(-(-1+Yfull(:,4)-Yfull(:,5))));
   R1 = rand(n,1)<p1;
   p2 = 1./(1+exp(-(-1+1.5*Yfull(:,1)+Yfull(:,4)-Yfull(:,5))));
   R2 = rand(n,1)<p2;
   p3 = 1./(1+exp(-(-1+Yfull(:,1)+Yfull(:,2)+Yfull(:,4)-Yfull(:,5))));
   R3 = rand(n,1)<p3;
   Y(R1,1)=NaN; Y(R2,2)=NaN; Y(R3,3)=NaN;
   out = mdMAARtest(Y,'method','ccm','plots',true);
   disp(out.whichReject)
%}
%
%{
   %% Example 3: Direct test with multiple stochastic imputations.
   rng(3)
   n = 500;
   Yfull = randn(n,5);
   Y = Yfull;
   pr = 1./(1+exp(-(-1+Yfull(:,4)-Yfull(:,5))));
   for j=1:3
       Y(rand(n,1)<pr,j)=NaN;
   end
   out = mdMAARtest(Y,'method','dtmm','nimputations',10);
   disp(out.DLtable)
%}

%{
   %% Example 4: Gaussian-copula diagnostic.
   % A small number of iterations is used only to keep the example fast.
   rng(4)
   n = 300;
   Yfull = randn(n,5);
   Y = Yfull;
   pr = 1./(1+exp(-(-1+Yfull(:,4)-Yfull(:,5))));
   for j=1:3
       Y(rand(n,1)<pr,j)=NaN;
   end
   out = mdMAARtest(Y,'method','cop','niterMCMC',600,'nburnMCMC',200, ...
       'plots',true);
   disp(out.posteriorMedian)
%}

%{
    %% Complex DTMM example: variables in all four diagnostic quadrants.
    rng(24680)
    
    n = 2500;
    q = 8;
    M = 4;
    
    % Two fully observed variables.
    X1 = randn(n,1);
    X2 = randn(n,1);
    
    % Eight variables that will contain missing values.
    Z = randn(n,q);
    
    % Complete data before introducing missing values.
    Ycomplete = [Z X1 X2];
    
    % Eight different missingness mechanisms.
    %
    % Variables 1-2:
    %   Missingness depends only on fully observed variables.
    %   Expected: low numerator, low missing information.
    %
    % Variables 3-4:
    %   Missingness additionally depends on another partially observed
    %   variable, but imputations will be very stable.
    %   Expected: high numerator, low missing information.
    %
    % Variables 5-6:
    %   No systematic dependence on partially observed variables, but
    %   imputations will deliberately disagree across completed data sets.
    %   Expected: low numerator, high missing information.
    %
    % Variables 7-8:
    %   Missingness depends on another partially observed variable and
    %   imputations will deliberately disagree.
    %   Expected: high numerator, high missing information.
    %
    % The simulation is designed to populate the four regions of the DTMM
    % diagnostic plot.
    %
    % LL1 and LL2 have a small likelihood-ratio numerator and stable
    % imputations. They are therefore expected in the lower-left region.
    %
    % HL1 and HL2 have a systematic dependence of missingness on partially
    % observed variables, producing a large numerator, while their imputations
    % are stable. They are therefore expected in the lower-right region.
    %
    % LH1 and LH2 have little systematic evidence in the numerator but their
    % imputation-specific fits disagree strongly. They are therefore expected
    % in the upper-left region.
    %
    % HH1 and HH2 combine a systematic missingness signal with substantial
    % disagreement among imputations. They are therefore expected in the
    % upper-right region.
    %
    % Importantly, the position of a point relative to the vertical reference
    % line is not by itself the final test decision. The final D_L statistic
    % also accounts for the missing-information component r_L.
    %
    % The comparison between HL and HH variables illustrates the role of r_L.
    % Both groups can have large values of d_L/q, but the HH variables have
    % much larger missing-information fractions. Consequently, their final
    % D_L statistics can be appreciably smaller and their p-values larger.
    %
    % Conversely, the LH variables show that a large missing-information
    % fraction alone does not produce rejection when the numerator contains
    % little evidence against the reduced missingness model.
    eta = zeros(n,q);
    
    % Low numerator / low missing information.
    eta(:,1) = -1.5 + 1.1*X1;
    eta(:,2) = -1.5 - 1.1*X2;
    
    % High numerator / low missing information.
    eta(:,3) = -1.0 + 0.35*Z(:,1) + 0.2*X1;
    eta(:,4) = -1.0 + 0.35*Z(:,2) - 0.2*X2;
    
    % Low numerator / high missing information.
    % These mechanisms are MCAR relative to all variables.
    eta(:,5) = 0;
    eta(:,6) = 0;
    
    % High numerator / high missing information.
    eta(:,7) = -0.5 + 0.35*Z(:,3) + 0.2*X1;
    eta(:,8) = -0.5 + 0.35*Z(:,4) - 0.2*X2;
    
    prob = 1./(1+exp(-eta));
    R = rand(n,q) < prob;
    
    % Introduce missing values.
    Y = Ycomplete;
    for j = 1:q
        Y(R(:,j),j) = NaN;
    end
    
    % Add informative variable names.
    varNames = {'LL1','LL2', ...
        'HL1','HL2', ...
        'LH1','LH2', ...
        'HH1','HH2', ...
        'X1','X2'};
    
    Y = array2table(Y,'VariableNames',varNames);
    
    % Construct four completed data sets.
    
    Imputations = cell(M,1);
    
    % Alternating signs are used to generate conflicting imputations for
    % variables 5-8.
    signImp = [-1 1 -1 1];
    
    % Small perturbation for stable imputations.
    stableNoise = 0.03;
    
    % Parameters controlling disagreement among unstable imputations.
    shift = 0.15;
    unstableNoise = 0.35;
    
    for m = 1:M
    
        % Start from the known complete simulated data.
        Ym = Ycomplete;
    
        % Variables 1-4: very stable imputations.
        for j = 1:4
            idx = R(:,j);
    
            Ym(idx,j) = Z(idx,j) + ...
                stableNoise*randn(sum(idx),1);
        end
    
        % Variables 5-8: deliberately conflicting imputations.
        %
        % The sign of the shift reverses across imputations. Consequently,
        % the imputation-specific logistic coefficients can disagree strongly,
        % increasing r_L.
        for j = 5:8
            idx = R(:,j);
    
            Ym(idx,j) = Z(idx,j) + ...
                signImp(m)*shift + ...
                unstableNoise*randn(sum(idx),1);
        end
    
        Imputations{m} = Ym;
    end
    
    %Run DTMM and produce the diagnostic plot.
    out = mdMAARtest(Y,'method','dtmm', ...
        'imputed',Imputations, ...
        'plots',true, ...
        'msg',true);
%}

%% Beginning of code

if nargin < 1
    error('FSDA:mdMAARtest:TooFewInputs', ...
        'At least one input argument is required.');
end

% Preserve variable names before converting tables to arrays.
if istimetable(Y)
    Y = timetable2table(Y,'ConvertRowTimes',false);
end
if istable(Y)
    variableNames = string(Y.Properties.VariableNames);
    Y = table2array(Y);
else
    variableNames = "Y" + string((1:size(Y,2))');
end

if ~ismatrix(Y) || ~isnumeric(Y)
    error('FSDA:mdMAARtest:WrongInput', ...
        'Input argument Y must be a numeric matrix, table or timetable.');
end

Y = double(Y);
[n,p] = size(Y);
if n < 5 || p < 2
    error('FSDA:mdMAARtest:SmallData', ...
        'At least five observations and two variables are required.');
end
if any(all(isnan(Y),1))
    error('FSDA:mdMAARtest:AllMissingVariable', ...
        'At least one variable contains only missing values.');
end

% Default options.
method = 'ccm';
alpha = 0.05;
imputed = [];
nimputations = 10;
maxiter = 100;
tol = 1e-5;
niterMCMC = 2000;
nburnMCMC = 500;
thinMCMC = 1;
pluginthreshold = 100;
ridge = 1e-8;
msg = true;
plots = false;
savesamples = false;

% Names of the options explicitly supplied by the user. This list is also
% used below to identify options which are irrelevant to the selected
% diagnostic method.
UserOptions = {};

if ~isempty(varargin)
    options = struct('method',method,'alpha',alpha,'imputed',imputed, ...
        'nimputations',nimputations,'maxiter',maxiter,'tol',tol, ...
        'niterMCMC',niterMCMC,'nburnMCMC',nburnMCMC,'thinMCMC',thinMCMC, ...
        'pluginthreshold',pluginthreshold,'ridge',ridge, ...
        'msg',msg,'plots',plots,'savesamples',savesamples);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mdMAARtest:WrongInputOpt', ...
            'Number of supplied options is invalid. Values may be missing.');
    end

    if ~isempty(UserOptions)
        % Match option names without regard to letter case and replace each
        % recognized name by the canonical field name used in options.
        %
        % This step is necessary for mixed-case names such as niterMCMC,
        % nburnMCMC and thinMCMC. For example, 'nitermcmc', 'NITERMCMC' and
        % 'niterMCMC' are all converted to the canonical name 'niterMCMC'
        % before validation and assignment.
        canonicalOptionNames = fieldnames(options);
        for iOpt = 1:numel(UserOptions)
            match = find(strcmpi(UserOptions{iOpt}, ...
                canonicalOptionNames),1);
            if ~isempty(match)
                UserOptions{iOpt} = canonicalOptionNames{match};
                varargin{2*iOpt-1} = canonicalOptionNames{match};
            end
        end

        % Unknown option names, which have not been canonicalized above,
        % are detected here.
        aux.chkoptions(options,UserOptions)
    end

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end

    method = options.method;
    alpha = options.alpha;
    imputed = options.imputed;
    nimputations = options.nimputations;
    maxiter = options.maxiter;
    tol = options.tol;
    niterMCMC = options.niterMCMC;
    nburnMCMC = options.nburnMCMC;
    thinMCMC = options.thinMCMC;
    pluginthreshold = options.pluginthreshold;
    ridge = options.ridge;
    msg = options.msg;
    plots = options.plots;
    savesamples = options.savesamples;
end

method = lower(char(method));
if ~ismember(method,{'ccm','dtmm','cop'})
    error('FSDA:mdMAARtest:WrongMethod', ...
        'Option method must be ''ccm'', ''dtmm'' or ''cop''.');
end

% Warn when the user supplies options which are not used by the selected
% diagnostic method. The warning is informative; the values are still
% checked by the general input-validation code below.
commonOptions = {'method','alpha','ridge','msg','plots'};
dtmmOptions = {'imputed','nimputations','maxiter','tol'};
copOptions = {'niterMCMC','nburnMCMC','thinMCMC','pluginthreshold','savesamples'};

switch method
    case 'ccm'
        validMethodOptions = commonOptions;
    case 'dtmm'
        validMethodOptions = [commonOptions dtmmOptions];
    case 'cop'
        validMethodOptions = [commonOptions copOptions];
end

if ~isempty(UserOptions)
    % Compare option names without regard to case. This is important for
    % canonical mixed-case names such as niterMCMC, nburnMCMC and thinMCMC.
    unusedMask = ~cellfun(@(name) ...
        any(strcmpi(name,validMethodOptions)),UserOptions);
    unusedOptions = UserOptions(unusedMask);

    if ~isempty(unusedOptions)
        warning('FSDA:mdMAARtest:UnusedOptions', ...
            ['The following option(s) are not used by method ''%s'': ' ...
            '%s.'], ...
            method,strjoin(unusedOptions,', '));
    end
end
if ~isscalar(alpha) || ~isnumeric(alpha) || alpha <= 0 || alpha >= 1
    error('FSDA:mdMAARtest:WrongAlpha', ...
        'Option alpha must be a scalar in the interval (0,1).');
end
if ~isscalar(nimputations) || nimputations <= 1 || ...
        nimputations ~= floor(nimputations)
    error('FSDA:mdMAARtest:WrongNImputations', ...
        'Option nimputations must be an integer greater than 1.');
end
if ~isscalar(maxiter) || maxiter <= 0 || maxiter ~= floor(maxiter)
    error('FSDA:mdMAARtest:WrongMaxiter', ...
        'Option maxiter must be a positive integer.');
end
if ~isscalar(tol) || tol <= 0
    error('FSDA:mdMAARtest:WrongTol', ...
        'Option tol must be a positive scalar.');
end
if ~isscalar(niterMCMC) || niterMCMC <= 1 || niterMCMC ~= floor(niterMCMC)
    error('FSDA:mdMAARtest:WrongniterMCMC', ...
        'Option niterMCMC must be an integer greater than 1.');
end
if ~isscalar(nburnMCMC) || nburnMCMC < 0 || nburnMCMC ~= floor(nburnMCMC) || nburnMCMC >= niterMCMC
    error('FSDA:mdMAARtest:WrongnburnMCMC', ...
        'Option nburnMCMC must be a nonnegative integer smaller than niterMCMC.');
end
if ~isscalar(thinMCMC) || thinMCMC < 1 || thinMCMC ~= floor(thinMCMC)
    error('FSDA:mdMAARtest:WrongThin', ...
        'Option thinMCMC must be a positive integer.');
end
if ~isscalar(pluginthreshold) || pluginthreshold < 2 || ...
        pluginthreshold ~= floor(pluginthreshold)
    error('FSDA:mdMAARtest:WrongPluginThreshold', ...
        'Option pluginthreshold must be an integer greater than 1.');
end

if ~isscalar(ridge) || ~isnumeric(ridge) || ridge < 0
    error('FSDA:mdMAARtest:WrongRidge', ...
        'Option ridge must be a nonnegative scalar.');
end
if ~(isscalar(msg) && (islogical(msg) || isnumeric(msg)))
    error('FSDA:mdMAARtest:WrongMsg', ...
        'Option msg must be a logical or numeric scalar.');
end
if ~(isscalar(plots) && (islogical(plots) || isnumeric(plots)))
    error('FSDA:mdMAARtest:WrongPlots', ...
        'Option plots must be a logical or numeric scalar.');
end
if ~(isscalar(savesamples) && (islogical(savesamples) || isnumeric(savesamples)))
    error('FSDA:mdMAARtest:WrongSaveSamples', ...
        'Option savesamples must be a logical or numeric scalar.');
end
msg = logical(msg);
plots = logical(plots);
savesamples = logical(savesamples);


% Identify the variables involved in the MAAR diagnostics.
% missingMask(i,j) is true when observation i is missing on variable j.
% locmis contains the column indices of the partially observed variables.
% locfull contains the column indices of the fully observed variables.
% q is the number of partially observed variables.
missingMask = isnan(Y);
locmis = find(any(missingMask,1));
locfull = find(~any(missingMask,1));
q = numel(locmis);
if q == 0
    error('FSDA:mdMAARtest:NoMissingValues', ...
        'Input data do not contain missing values.');
end

% Constant variables make all three diagnostics ill-defined.
for j = 1:p
    yj = Y(~isnan(Y(:,j)),j);
    if numel(unique(yj)) < 2
        error('FSDA:mdMAARtest:ConstantVariable', ...
            'Variable %d is constant on its observed values.',j);
    end
end

% Construct the missingness-indicator matrix for the partially observed
% variables. R(i,k)=1 if variable locmis(k) is missing for observation i
% and R(i,k)=0 otherwise.
R = double(missingMask(:,locmis));

% Common output fields.
out = struct;
out.method = method;
out.alpha = alpha;
out.n = n;
out.p = p;
out.nvarmiss = q;
out.variableNames = variableNames;
out.missingVariables = locmis(:);
out.fullyObservedVariables = locfull(:);
out.missingIndicator = R;
out.reject = false;
out.whichReject = false(q,1);
out.pvalue = [];
out.stat = [];
out.alphaAdjusted = alpha;
out.interpretation = '';

switch method
    case 'ccm'
        % Inputs passed to ccmTest:
        % Y       = original n-by-p data matrix containing NaNs;
        % R       = n-by-q missingness-indicator matrix;
        % locmis  = indices of the q partially observed variables;
        %           For example locmis=[2 3 5] means that the variables
        %           which contain missing are 2nd 3rd and 5th
        % locfull = indices of the fully observed variables;
        %           For example locfull=[1 4] means that the variables
        %           which contain missing are 1st and 4th
        % alpha   = nominal significance level;
        % ridge   = numerical regularization for nearly singular systems.
        res = ccmTest(Y,R,locmis,locfull,alpha,ridge); % comparison of conditional means test
        out = copyFields(out,res);

    case 'dtmm'
        % Create or validate the completed data sets used by the direct test.
        % imputed is either empty, a cell array, or an n-by-p-by-M array;
        % nimputations is used only when imputations must be generated.
        [Imputations,impInfo] = prepareImputations(Y,imputed,nimputations, ...
            maxiter,tol);

        % Inputs passed to dtmmTest:
        % Imputations = cell of length nimputations containing M completed n-by-p data matrices;
        % R           = n-by-q missingness-indicator matrix;
        % locmis      = indices of partially observed variables;
        % locfull     = indices of fully observed variables;
        % alpha       = nominal significance level;
        % ridge       = numerical regularization;
        % maxiter,tol = controls for the logistic-regression iterations.
        res = dtmmTest(Imputations,R,locmis,locfull,alpha,ridge,maxiter,tol);
        out = copyFields(out,res);

        % Create a compact diagnostic table for the Meng-Rubin D_L test.
        % Rows correspond to the variables containing missing values.
        rowNames = cellstr(out.variableNames(out.missingVariables));
        out.DLtable = table(out.dL,out.relativeIncrease,out.stat,out.pvalue, ...
            'VariableNames',{'dL','rL','DL','pvalue'}, ...
            'RowNames',rowNames);

        % dL is returned by dtmmTest only to construct DLtable; the table is
        % the public summary output, while out.deviance retains full details.
        out = rmfield(out,'dL');
        out.imputationInfo = impInfo;

    case 'cop'
        % Inputs passed to copulaTest:
        % Y,R,locmis,locfull identify the data and missingness structure;
        % alpha is the nominal significance level;
        % niterMCMC, nburnMCMC and thinMCMC control the MCMC simulation;
        % pluginthreshold selects margins treated through fixed normal scores;
        % ridge stabilizes matrix operations; savesamples controls storage;
        % msg controls progress messages.
        res = copulaTest(Y,R,locmis,locfull,alpha,niterMCMC,nburnMCMC,thinMCMC, ...
            pluginthreshold,ridge,savesamples,msg);
        out = copyFields(out,res);
end

out.interpretation = interpretMAAR(out);

if msg
    printMAAR(out)
end
if plots
    plotMAAR(out)
end

end

% -------------------------------------------------------------------------
function out = copyFields(out,res)
%copyFields copies all fields of a method-specific result into out.
%
% Inputs:
%   out : structure containing the fields common to all MAAR methods.
%   res : structure returned by ccmTest, dtmmTest or copulaTest.
%
% Output:
%   out : input structure augmented or overwritten with the fields of res.
fields = fieldnames(res);
for i = 1:numel(fields)
    out.(fields{i}) = res.(fields{i});
end
end

% -------------------------------------------------------------------------
function res = ccmTest(Y,R,locmis,locfull,alpha,ridge)
%ccmTest performs the comparison-of-conditional-means diagnostic.
%
% Inputs:
%   Y       : n-by-p data matrix containing NaNs in missing positions.
%   R       : n-by-q missingness-indicator matrix for the variables in
%             locmis. R(i,k)=1 means that Y(i,locmis(k)) is missing.
%   locmis  : 1-by-q vector containing the column indices of variables with
%             at least one missing value.
%   locfull : vector containing the column indices of variables observed for
%             every unit. It can be empty.
%   alpha   : nominal significance level used before Bonferroni correction.
%   ridge   : nonnegative regularization used only when a least-squares
%             system is nearly singular.
%
% Output:
%   res     : structure containing pairwise p-values, F statistics, degrees
%             of freedom, adjusted alpha and rejection indicators.
%
% For each pair (k,j), with k and j = 1,2, ..., q, with k different from j,
% the response Y_j is regressed on the fully observed variables. The larger
% model also includes R_k and the interactions between R_k and all fully
% observed variables.
q = numel(locmis);
pvalue = ones(q,q);
stat = zeros(q,q);
df1 = zeros(q,q);
df2 = zeros(q,q);

% Matrix of covariates that are observed for every unit. If locfull is
% empty, the reduced model contains only the intercept.
Xfull = Y(:,locfull);
% jj select the partially observed response variable among
% the variables which contain missing values. jj will be put in column jj
% of pvalue matrix
for jj = 1:q
    response = Y(:,locmis(jj));
    % kk loops on the variables which contain missing values which are
    % different from jj
    for kk = 1:q
        if kk == jj
            continue
        end
        % Missingness indicator defining the two comparison groups.
        rkk = R(:,kk);
        % Only units for which the current response is observed can enter
        % this conditional-mean comparison.
        use = ~isnan(response);
        y = response(use);
        % Note that matrix Xfull does not contain missing values
        Xfulluse = Xfull(use,:);
        rkkuse = rkk(use);
        % Reduced model: intercept plus fully observed covariates.
        Xsmall = [ones(numel(y),1) Xfulluse];
        % Full model: reduced terms, group effect and group-by-covariate
        % interactions. Expands rkkuse across the columns of Xfulluse.
        Xlarge = [Xsmall rkkuse Xfulluse.*rkkuse];

        % The question is: after conditioning on all the fully observed
        % variables, knowing whether another partially observed variable is
        % missing provides extra additional information about Y_j?​
        test = nestedLinearF(y,Xsmall,Xlarge,ridge);
        pvalue(kk,jj) = test.pvalue;
        stat(kk,jj) = test.stat;
        df1(kk,jj) = test.df1;
        df2(kk,jj) = test.df2;
    end
end

% There are q(q-1) off-diagonal comparisons. Use a Bonferroni component
% level; max(...,1) avoids division by zero when q=1.
nTests = max(q*(q-1),1);
alphaAdjusted = alpha/nTests;
rejectPairs = pvalue < alphaAdjusted;
rejectPairs(1:q+1:end) = false;
whichReject = any(rejectPairs,2);

if q == 1
    warning('FSDA:mdMAARtest:CCMOneMissingVariable', ...
        ['The CCM procedure requires cross-comparisons between partially ' ...
        'observed variables. With one such variable there is no component test.']);
end

res = struct;
res.pvalue = pvalue;
res.stat = stat;
res.df1 = df1;
res.df2 = df2;
res.alphaAdjusted = alphaAdjusted;
res.rejectPairs = rejectPairs;
res.whichReject = whichReject;
res.reject = any(whichReject);
end

% -------------------------------------------------------------------------
function test = nestedLinearF(y,Xsmall,Xlarge,ridge)
%nestedLinearF compares two nested Gaussian linear models.
%
% Inputs:
%   y       : response vector.
%   Xsmall  : design matrix of the reduced model.
%   Xlarge  : design matrix of the full model; its column space must contain
%             that of Xsmall.
%   ridge   : nonnegative numerical regularization passed to
%             stableLeastSquares.
%
% Output:
%   test    : structure with fields stat, df1, df2 and pvalue.
%
% The statistic compares the reduction in residual sum of squares with the
% residual mean square of the full model.
ns = size(Xsmall,1);
rs = rank(Xsmall);
rl = rank(Xlarge);
dfdiff = rl-rs;
dfres = ns-rl;

if dfdiff <= 0 || dfres <= 0
    test = struct('stat',NaN,'df1',dfdiff,'df2',dfres,'pvalue',NaN);
    return
end

bs = stableLeastSquares(Xsmall,y,ridge);
bl = stableLeastSquares(Xlarge,y,ridge);
rsdSmall = y-Xsmall*bs;
rsdLarge = y-Xlarge*bl;
sseSmall = sum(rsdSmall.^2);
sseLarge = sum(rsdLarge.^2);
num = max(sseSmall-sseLarge,0)/dfdiff;
den = sseLarge/dfres;
if den <= 0
    F = Inf;
    pval = 0;
else
    F = num/den;
    pval = fSurvival(F,dfdiff,dfres);
end

test = struct('stat',F,'df1',dfdiff,'df2',dfres,'pvalue',pval);
end

% -------------------------------------------------------------------------
function b = stableLeastSquares(X,y,ridge)
%stableLeastSquares solves a linear least-squares problem robustly.
%
% Inputs:
%   X     : design matrix.
%   y     : response vector.
%   ridge : nonnegative regularization used when X is rank deficient or
%           numerically ill conditioned.
%
% Output:
%   b     : vector of estimated regression coefficients.
%
% The ordinary backslash solution is used whenever the triangular factor
% from QR decomposition is well conditioned. Otherwise a ridge-stabilized
% normal-equation solution is used.
[~,R] = qr(X,0);
if isempty(R) || size(R,1) ~= size(R,2) || rcond(R) < 1e-12
    scale = max(1,trace(X'*X)/size(X,2));
    ridgeUse = max(ridge,eps(scale));
    b = (X'*X+ridgeUse*scale*eye(size(X,2)))\(X'*y);
else
    b = X\y;
end
end

% -------------------------------------------------------------------------
function [imputations,info] = prepareImputations(Y,imputed,nStockImp,maxiter,tol)
%prepareImputations validates or generates completed data sets.
%
% Inputs:
%   Y          : original n-by-p incomplete data matrix.
%   imputed    : empty, a cell array of completed matrices, one completed
%                matrix, or an n-by-p-by-M numeric array.
%   nStockImp  : number of stochastic imputations to generate when imputed is
%                empty.
%   maxiter : maximum number of iterations used by mdEM.
%   tol     : convergence tolerance used by mdEM.
%
% Outputs:
%   imputations : nStockImp-by-1 cell array of completed n-by-p matrices.
%   info        : structure recording the imputation source, number of
%                 imputations and, when available, EM estimates.
%
% When imputed is empty, mdEM estimates location and covariance and
% mdImputeStochastic generates independent completed data sets.
[n,p] = size(Y);
if isempty(imputed)
    if exist('mdEM','file') ~= 2 || exist('mdImputeStochastic','file') ~= 2
        error('FSDA:mdMAARtest:MissingFSDAFunctions', ...
            ['Functions mdEM and mdImputeStochastic must be on the MATLAB ' ...
            'path when option imputed is empty.']);
    end
    emOut = mdEM(Y,'maxiter',maxiter,'tol',tol);
    imputations = cell(nStockImp,1);
    for m = 1:nStockImp
        imputations{m} = mdImputeStochastic(Y,emOut.loc,emOut.cov);
    end
    info = struct('source','mdEM/mdImputeStochastic', ...
        'nimputations',nStockImp,'loc',emOut.loc,'cov',emOut.cov);
else
    if iscell(imputed)
        imputations = imputed(:);
    elseif isnumeric(imputed) && ndims(imputed) <= 3
        if ismatrix(imputed)
            imputations = {double(imputed)};
        else
            imputations = cell(size(imputed,3),1);
            for m = 1:size(imputed,3)
                imputations{m} = double(imputed(:,:,m));
            end
        end
    else
        error('FSDA:mdMAARtest:WrongImputed', ...
            'Option imputed must be a cell array or an n x p x M numeric array.');
    end

    for m = 1:numel(imputations)
        if istable(imputations{m})
            imputations{m} = table2array(imputations{m});
        end
        if ~isnumeric(imputations{m}) || ...
                ~isequal(size(imputations{m}),[n p])
            error('FSDA:mdMAARtest:WrongImputedSize', ...
                'Every completed data set must have the same size as Y.');
        end
        imputations{m} = double(imputations{m});
        if any(isnan(imputations{m}(:)))
            error('FSDA:mdMAARtest:IncompleteImputation', ...
                'A supplied completed data set still contains NaN values.');
        end
        obs = ~isnan(Y);
        tolerance = 100*eps(max(1,max(abs(Y(obs)))));
        if any(abs(imputations{m}(obs)-Y(obs)) > tolerance)
            warning('FSDA:mdMAARtest:ObservedValuesChanged', ...
                'Supplied imputation %d changes at least one observed value.',m);
        end
    end
    if numel(imputations) < 2
        error('FSDA:mdMAARtest:TooFewImputations', ...
            'At least two completed data sets are required for method dtmm.');
    end
    info = struct('source','user supplied','nimputations',numel(imputations), ...
        'loc',[],'cov',[]);
end
end

% -------------------------------------------------------------------------
function res = dtmmTest(imputations,R,locmis,locfull,alpha,ridge,maxiter,tol)
%dtmmTest directly tests each postulated missingness mechanism.
%
% Inputs:
%   imputations : cell array containing M completed n-by-p data matrices.
%   R           : n-by-q missingness-indicator matrix. Column k is the
%                 binary response for the model of the kth mechanism.
%   locmis      : indices of the q partially observed variables.
%   locfull     : indices of the fully observed variables. It can be empty.
%   alpha       : nominal significance level before Bonferroni correction.
%   ridge       : nonnegative regularization for nearly singular information
%                 matrices.
%   maxiter     : maximum number of IRLS iterations for each logistic fit.
%   tol         : convergence tolerance for the logistic coefficients.
%
% Output:
%   res         : structure containing the pooled D_L statistics, p-values,
%                 degrees of freedom, relative increases in variance and
%                 rejection indicators. The D_L procedure is also commonly
%                 called the D3 procedure.
%
% For missingness indicator R_k, the reduced logistic model contains only
% the fully observed variables. The full logistic model additionally
% contains all the partially observed variables after imputation.
%
% For each missingness mechanism, the response and the fully observed
% variables do not change across imputations. The reduced model is
% therefore identical in all M completed data sets and is fitted only once
% for each value of kk.

M = numel(imputations);

if M < 2
    error('FSDA:dtmmTest:TooFewImputations', ...
        'At least two completed data sets are required for the D_L test.');
end

q = numel(locmis);
if q < 1
    error('FSDA:dtmmTest:NoPartiallyObservedVariables', ...
        'At least one partially observed variable is required.');
end

pvalue = NaN(q,1);
stat = NaN(q,1);
df1 = q*ones(q,1);
df2 = NaN(q,1);
relativeIncrease = NaN(q,1);
dL = NaN(q,1);

% For each missingness mechanism, fullMLE and fullPooled are M-by-1
% vectors. smallMLE and smallPooled are scalars because the reduced model
% is identical across completed data sets.
deviance = repmat(struct('fullMLE',[],'smallMLE',[], ...
    'fullPooled',[],'smallPooled',[]),q,1);

% devfullMLEt=array2table(zeros(M,q),"RowNames",string(1:M));
% devfullPooled=zeros(M,q);
% devsmallMLE=zeros(1,q);
% devsmallPooled=zeros(1,q);


% The reduced model contains only the intercept and the fully observed
% variables. Its response and design matrix are therefore identical
% across completed data sets.
Y1 = imputations{1};
Xsmall = [ones(size(Y1,1),1) Y1(:,locfull)];

% Number of restrictions tested. The full model adds one coefficient
% for each partially observed variable.
% k = size(bFull,1)-size(bSmall,1);
k=q;

for kk = 1:q

    % Binary response for the missingness mechanism of variable locmis(kk).
    response = R(:,kk);


    % Fit the common reduced model once for the current missingness mechanism.
    fitSmall = logisticFit(Xsmall,response,maxiter,tol,ridge);
    % bSmall = fitSmall.beta; (not used)
    devSmallMLE = fitSmall.deviance;

    % The full model includes the imputed partially observed variables and
    % must therefore be fitted separately to every completed data set.
    bFull = zeros(length(locfull)+length(locmis)+1,M);
    devFullMLE = zeros(M,1);

    for m = 1:M
        Ym = imputations{m};
        Xfull = [Xsmall Ym(:,locmis)];

        fitFull = logisticFit(Xfull,response,maxiter,tol,ridge);

        bFull(:,m) = fitFull.beta;
        devFullMLE(m) = fitFull.deviance;
    end

    % The Meng-Rubin D_L procedure requires a common coefficient vector for
    % the full model. It is obtained as the arithmetic mean of the M
    % imputation-specific full-model estimates:
    %
    %   qbarFull = (1/M) * sum_m betaFull_m.
    %
    % No analogous averaging is required for the reduced model because its
    % response and design matrix are identical across completed data sets.
    % Thus, (bSmall) the reduced-model coefficient estimate is already common to all
    % completed data sets.
    qbarFull = mean(bFull,2);

    % The reduced model is identical across imputations. Therefore, its pooled
    % deviance is identical to the reduced-model fitted deviance.
    devSmallPooled = devSmallMLE;

    % Reevaluate the full model in each completed data set at qbarFull,
    % without re-estimating its coefficients.
    %
    % devFullMLE(m) is the full-model deviance evaluated at the coefficient
    % estimate specific to completed data set m.
    %
    % devFullPooled(m) is the full-model deviance for completed data set m
    % evaluated at the common coefficient vector qbarFull.
    devFullPooled = zeros(M,1);

    for m = 1:M
        Ym = imputations{m};
        Xfull = [Xsmall Ym(:,locmis)];

        devFullPooled(m) = ...
            logisticDeviance(Xfull,response,qbarFull);
    end

    % devL is the likelihood-ratio statistic obtained after forcing all
    % completed data sets to use the common full-model coefficient vector:
    %
    %   devL = Dsmall(betaSmallHat)
    %          - mean_m{Dfull_m(qbarFull)},
    %
    % where Dsmall and Dfull denote minus twice the corresponding Bernoulli
    % log likelihoods.
    devL = devSmallPooled - mean(devFullPooled);


    % In the Meng-Rubin formulation, the corresponding average
    % likelihood-ratio statistic based on the imputation-specific
    % full-model estimates is
    %
    %   devM = Dsmall(betaSmallHat)
    %          - mean_m{Dfull_m(betaFullHat_m)}.
    %
    % The common reduced-model deviance cancels from devM-devL:
    %
    %   devM-devL
    %     = mean_m{
    %         Dfull_m(qbarFull)-Dfull_m(betaFullHat_m)
    %       }.
    %
    % Compute this difference directly to avoid subtracting two potentially
    % similar likelihood-ratio statistics.

    % devMminusDevL measures the loss of fit caused by replacing the
    % imputation-specific full-model estimates with the common pooled estimate.
    % For each completed data set, devFullPooled(m)-devFullMLE(m) is the increase
    % in deviance obtained by evaluating the model at qbarFull rather than at
    % its own MLE. Averaging these increases quantifies the disagreement among
    % the imputations and provides a likelihood-based measure of the additional
    % uncertainty due to missing data.
    devMminusDevL = mean(devFullPooled-devFullMLE);

    % Estimate the average relative increase in variance due to missing data.
    % Equivalently, r_L estimates the average odds of missing information:
    %
    %   r_L = lambda/(1-lambda),
    %
    % where lambda denotes the fraction of missing information. This is r_L
    % in equation (3.8) of Meng and Rubin (1992):
    %
    %   r_L = (M+1)/(k*(M-1)) * (devM-devL).    %
    rL = ((M+1)/(k*(M-1))) * devMminusDevL;

    % When fitFull.beta minimizes the same unpenalized Bernoulli deviance
    % evaluated by logisticDeviance, devFullPooled(m) cannot be smaller
    % than devFullMLE(m). Small negative values of devMminusDevL may still
    % occur because of numerical rounding. A substantial negative value
    % can indicate that the fitting and deviance calculations do not use
    % exactly the same objective.
    rL = max(rL,0);

    % Numerator of the pooled D_L statistic. The maximum with zero protects
    % against tiny negative values caused by finite-precision calculations.
    dL(kk) = max(devL,0);

    % Pooled D_L statistic. The numerator dL is divided by its number of
    % restrictions and adjusted for the estimated increase in variance
    % caused by missing information.
    D_L = dL(kk)/(k*(1+rL));

    % Degrees of freedom entering the finite-M reference distribution.
    v = k*(M-1);

    % Denominator degrees of freedom of the F reference distribution,
    % following equation (2.7) of Meng and Rubin (1992). When rL is
    % essentially zero, the denominator degrees of freedom are infinite and
    % the reference distribution reduces to the corresponding chi-square
    % approximation divided by k.
    if rL <= sqrt(eps)
        w = Inf;
    elseif v > 4
        w = 4+(v-4)*(1+(1-2/v)/rL)^2;
    else
        w = v*(1+1/k)*(1+1/rL)^2/2;
    end

    % A small p-value for component kk indicates that, after conditioning on
    % the fully observed variables, at least one partially observed variable
    % provides additional information about whether variable locmis(kk) is
    % missing. This is evidence against the restrictions implied by MAAR or
    % against one of the accompanying assumptions. The test is omnibus and
    % does not by itself identify which partially observed variable is
    % responsible.
    pvalue(kk) = fSurvival(D_L,k,w);
    stat(kk) = D_L;
    df1(kk) = k;
    df2(kk) = w;
    relativeIncrease(kk) = rL;

    % Store the quantities entering the D_L calculation. The reduced-model
    % deviances are scalars, while the full-model deviances contain one
    % value for each completed data set.
    deviance(kk).fullMLE = devFullMLE;
    deviance(kk).smallMLE = devSmallMLE;
    deviance(kk).fullPooled = devFullPooled;
    deviance(kk).smallPooled = devSmallPooled;
end

% Apply a Bonferroni correction across the q missingness mechanisms.
alphaAdjusted = alpha/q;
whichReject = pvalue < alphaAdjusted;

res = struct;
res.pvalue = pvalue;
res.stat = stat;
res.df1 = df1;
res.df2 = df2;
res.relativeIncrease = relativeIncrease;
res.dL = dL;
res.deviance = deviance;
res.alphaAdjusted = alphaAdjusted;
res.whichReject = whichReject;
res.reject = any(whichReject);
end

% The old (not efficient version) which computes the coefficient for each reduced model
% is left below
% function res = dtmmTest(imputations,R,locmis,locfull,alpha,ridge,maxiter,tol)
% %dtmmTest directly tests each postulated missingness mechanism.
% %
% % For missingness indicator R_k, the reduced logistic model contains the
% % fully observed variables and the full model additionally contains all
% % imputed partially observed variables.
% M = numel(imputations);
% q = numel(locmis);
% pvalue = NaN(q,1);
% stat = NaN(q,1);
% df1 = q*ones(q,1);
% df2 = NaN(q,1);
% relativeIncrease = NaN(q,1);
% deviance = repmat(struct('fullMLE',[],'smallMLE',[], ...
%     'fullPooled',[],'smallPooled',[]),q,1);
%
% for kk = 1:q
%     % Binary response for the missingness mechanism of variable locmis(kk).
%     response = R(:,kk);
%     devFullMLE = zeros(M,1);
%     devSmallMLE = zeros(M,1);
%     bSmall = zeros(length(locfull)+1,M);
%     bFull = zeros(length(locfull)+length(locmis)+1,M);
%
%     for m = 1:M
%         Ym = imputations{m};
%         % Reduced model uses only variables observed on all units.
%         Xsmall = [ones(size(Ym,1),1) Ym(:,locfull)];
%         % Full model adds all partially observed variables after imputation.
%         Xfull = [Xsmall Ym(:,locmis)];
%         fitSmall = logisticFit(Xsmall,response,maxiter,tol,ridge); % glmfit(Xsmall(:,2:end),response,"binomial","Link","logit")
%         fitFull = logisticFit(Xfull,response,maxiter,tol,ridge);
%         % if m == 1
%         %     bSmall = zeros(numel(fitSmall.beta),M);
%         %     bFull = zeros(numel(fitFull.beta),M);
%         % end
%         bSmall(:,m) = fitSmall.beta;
%         bFull(:,m) = fitFull.beta;
%         devSmallMLE(m) = fitSmall.deviance;
%         devFullMLE(m) = fitFull.deviance;
%     end
%
%     % The Meng-Rubin D_L likelihood-ratio procedure requires a common
%     % parameter vector for each model. The common vector is the arithmetic
%     % mean of the M imputation-specific estimates:
%     %
%     %   qbarSmall = (1/M) * sum_m betaSmall_m;
%     %   qbarFull  = (1/M) * sum_m betaFull_m.
%     %
%     % In the next loop, the likelihood of every completed dataset is
%     % evaluated at these pooled coefficients, without re-estimation. The
%     % difference between the likelihood-ratio statistic evaluated at the
%     % imputation-specific MLEs and the statistic evaluated at the pooled
%     % coefficients measures disagreement between imputations and is used to
%     % estimate the relative increase in variance due to missing values.
%     %
%     % The reduced and full coefficients must be pooled separately because
%     % the two models have different dimensions and the coefficients of the
%     % common predictors need not be identical after additional predictors
%     % are included in the full logistic model.
%     qbarSmall = mean(bSmall,2);
%     qbarFull = mean(bFull,2);
%
%     % Reevaluate each completed-data model at the pooled coefficients.
%     % Unlike the first loop, no model fitting is performed here.
%     % devSmallPooled(m) and devFullPooled(m) are the deviances of completed
%     % dataset m when all imputations are forced to use the same respective
%     % pooled parameter vectors.
%     devSmallPooled = zeros(M,1);
%     devFullPooled = zeros(M,1);
%     for m = 1:M
%         Ym = imputations{m};
%         Xsmall = [ones(size(Ym,1),1) Ym(:,locfull)];
%         Xfull = [Xsmall Ym(:,locmis)];
%         devSmallPooled(m) = logisticDeviance(Xsmall,response,qbarSmall);
%         devFullPooled(m) = logisticDeviance(Xfull,response,qbarFull);
%     end
%
%
%     % devM is the average likelihood-ratio statistic computed using the
%     % separate maximum-likelihood estimates from each imputation:
%     %
%     %   devM = mean_m{Dsmall_m(betaSmall_m)-Dfull_m(betaFull_m)}.
%     %
%     % devL is the average likelihood-ratio statistic obtained by evaluating
%     % each completed dataset at the common pooled coefficient vectors:
%     %
%     %   devL = mean_m{Dsmall_m(qbarSmall)-Dfull_m(qbarFull)}.
%     %
%     % The difference devM-devL quantifies the loss of fit caused by forcing
%     % all imputations to share common coefficients. After scaling, this
%     % difference estimates the relative increase in variance due to
%     % nonresponse used in the Meng-Rubin D_L statistic.
%     devM = mean(devSmallMLE-devFullMLE);
%     devL = mean(devSmallPooled-devFullPooled);
%     k = size(bFull,1)-size(bSmall,1);
%
%     % Estimate the average odds of the fraction of missing information.
%     % This is r_L in equation (3.8) of Meng and Rubin (1992):
%     %
%     %   r_L = (M+1)/(k*(M-1)) * (dbar_M - dbar_L),
%     %
%     % where devM is the average likelihood-ratio statistic evaluated at the
%     % imputation-specific MLEs and devL is the average likelihood-ratio
%     % statistic evaluated at the pooled coefficient estimates.
%     rL = ((M+1)/(k*(M-1)))*(devM-devL);
%
%
%     % A slightly negative value can occur because of numerical error.
%     rL = max(rL,0);
%     D_L = max(devL,0)/(k*(1+rL));
%     v = k*(M-1);
%     % w = degrees of freedom of the denominator of the F distribution
%     % Equation (2.7) of Meng and Rubin (1992)
%     if rL <= sqrt(eps)
%         w = Inf;
%     elseif v > 4
%         w = 4+(v-4)*(1+(1-2/v)/rL)^2;
%     else
%         w = v*(1+1/k)*(1+1/rL)^2/2;
%     end
%
%     % A small p-value for component kk indicates that, after conditioning on
%     % the fully observed variables, at least one partially observed variable
%     % provides additional information about whether variable locmis(kk) is
%     % missing. This is evidence against the restrictions implied by MAAR or
%     % against one of the accompanying assumptions. The test is omnibus and
%     % does not by itself identify which partially observed variable is
%     % responsible.
%     pvalue(kk) = fSurvival(D_L,k,w);
%     stat(kk) = D_L;
%     df1(kk) = k;
%     df2(kk) = w;
%     relativeIncrease(kk) = rL;
%     deviance(kk).fullMLE = devFullMLE;
%     deviance(kk).smallMLE = devSmallMLE;
%     deviance(kk).fullPooled = devFullPooled;
%     deviance(kk).smallPooled = devSmallPooled;
% end
%
% alphaAdjusted = alpha/q;
% whichReject = pvalue < alphaAdjusted;
% res = struct;
% res.pvalue = pvalue;
% res.stat = stat;
% res.df1 = df1;
% res.df2 = df2;
% res.relativeIncrease = relativeIncrease;
% res.deviance = deviance;
% res.alphaAdjusted = alphaAdjusted;
% res.whichReject = whichReject;
% res.reject = any(whichReject);
% end


% -------------------------------------------------------------------------
function fit = logisticFit(X,y,maxiter,tol,ridge)
%logisticFit computes a binary logistic maximum-likelihood fit.
%
% Inputs:
%   X       : design matrix, including an intercept when required.
%   y       : binary response vector containing zeros and ones.
%   maxiter : maximum number of safeguarded Newton/IRLS iterations.
%   tol     : relative convergence tolerance for the coefficients.
%   ridge   : nonnegative regularization used if the information matrix is
%             nearly singular.
%
% Output:
%   fit     : structure containing beta, loglik, deviance, iter and
%             converged.
%
% Step halving ensures that the log likelihood does not decrease after a
% Newton update.
y = double(y(:));
X = double(X);
k = size(X,2);
beta = zeros(k,1);
ll = logisticLogLik(X,y,beta);
converged = false;

for iter = 1:maxiter
    eta = X*beta;
    mu = logisticCDF(eta);
    w = max(mu.*(1-mu),1e-10);
    gradient = X'*(y-mu);
    information = X'*(w.*X);
    scale = max(1,trace(information)/max(k,1));
    if rcond(information) < 1e-12
        ridgeUse = max(ridge,eps(scale));
        information = information+ridgeUse*scale*eye(k);
    end
    step = information\gradient;

    stepFactor = 1;
    betaNew = beta+step;
    llNew = logisticLogLik(X,y,betaNew);
    while llNew < ll && stepFactor > 2^-20
        stepFactor = stepFactor/2;
        betaNew = beta+stepFactor*step;
        llNew = logisticLogLik(X,y,betaNew);
    end

    if max(abs(betaNew-beta)) <= tol*(1+max(abs(beta)))
        beta = betaNew;
        ll = llNew;
        converged = true;
        break
    end
    beta = betaNew;
    ll = llNew;
end

if ~converged
    warning('FSDA:mdMAARtest:LogisticNoConvergence', ...
        'A logistic model did not converge in %d iterations.',maxiter);
end
fit = struct('beta',beta,'loglik',ll,'deviance',-2*ll, ...
    'iter',iter,'converged',converged);
end

% -------------------------------------------------------------------------
function dev = logisticDeviance(X,y,beta)
%logisticDeviance returns minus twice the Bernoulli log likelihood.
% X is the design matrix, y is the binary response and beta is a fixed
% coefficient vector at which the deviance must be evaluated.
dev = -2*logisticLogLik(X,double(y(:)),beta);
end

% -------------------------------------------------------------------------
function ll = logisticLogLik(X,y,beta)
%logisticLogLik evaluates the Bernoulli logistic log likelihood.
% X is the design matrix, y is a zero-one response and beta is the vector
% of logistic-regression coefficients.
eta = X*beta;
ll = sum(y.*(-softplus(-eta))+(1-y).*(-softplus(eta)));
end

% -------------------------------------------------------------------------
function y = softplus(x)
% softplus evaluates log(1+exp(x)) without numerical overflow.
% see https://en.wikipedia.org/wiki/Softplus
y = max(x,0)+log1p(exp(-abs(x)));
end

% -------------------------------------------------------------------------
function p = logisticCDF(x)
%logisticCDF evaluates 1/(1+exp(-x)) using a stable two-branch formula.
p = zeros(size(x));
pos = x >= 0;
p(pos) = 1./(1+exp(-x(pos)));
ex = exp(x(~pos));
p(~pos) = ex./(1+ex);
end

% -------------------------------------------------------------------------
function res = copulaTest(Y,R,locmis,locfull,alpha,niterMCMC,nburnMCMC,thinMCMC, ...
    pluginthreshold,ridge,savesamples,msg)
%copulaTest performs the semiparametric Gaussian-copula diagnostic.
%
% Inputs:
%
%   Y               : n-by-p incomplete data matrix. Missing entries are
%                     represented by NaN.
%
%   R               : n-by-q missingness-indicator matrix associated with
%                     the variables listed in locmis. In this implementation
%
%                         R(i,k)=1
%
%                     means that Y(i,locmis(k)) is missing. Notice that this
%                     is the reverse of the convention used by Bojinov,
%                     Pillai and Rubin, where R(i,k)=1 means observed. The
%                     reversal may change the sign of an association but
%                     does not affect whether it is equal to zero.
%
%   locmis          : vector containing the column indices of the q
%                     variables of Y having at least one missing value.
%
%   locfull         : vector containing the column indices of the variables
%                     observed for all n units. It can be empty.
%
%   alpha           : nominal significance level used to construct the
%                     simultaneous posterior credible intervals.
%
%   niterMCMC           : total number of MCMC iterations, including burn-in.
%
%   nburnMCMC           : number of initial MCMC iterations discarded as
%                     burn-in.
%
%   thinMCMC            : thinning interval. After burn-in, one posterior draw
%                     is retained every thin iterations.
%
%   pluginthreshold : margins having more than pluginthreshold distinct
%                     observed values use fixed rank-based normal scores for
%                     their observed entries. This approximation is useful
%                     for nearly continuous margins because it avoids
%                     repeatedly sampling a large number of latent values.
%
%   ridge           : nonnegative numerical regularization used when a
%                     covariance matrix is singular or nearly singular.
%
%   savesamples     : if true, the retained posterior draws of the
%                     indicator-outcome conditional covariance matrix are
%                     returned in res.samples. Rows correspond to
%                     missingness indicators and columns to partially
%                     observed outcomes.
%
%   msg             : if true, progress messages are displayed during the
%                     MCMC iterations.
%
% Output:
%
%   res             : structure containing the following fields:
%
%                     pvalue
%                         q-by-q matrix of two-sided posterior tail
%                         probabilities. Row k corresponds to missingness
%                         indicator R_k and column j to partially observed
%                         outcome Y_j. These are posterior sign
%                         probabilities, not classical frequentist
%                         p-values.
%
%                     stat
%                         q-by-q matrix containing posterior medians of the
%                         conditional covariances. Rows correspond to R_k
%                         and columns to Y_j.
%
%                     CI
%                         q-by-q-by-2 array containing the lower and upper
%                         Bonferroni-adjusted credible limits, with rows
%                         corresponding to R_k and columns to Y_j.
%
%                     posteriorMean
%                         q-by-q matrix of posterior means, with rows
%                         corresponding to R_k and columns to Y_j.
%
%                     posteriorMedian
%                         q-by-q matrix of posterior medians, with rows
%                         corresponding to R_k and columns to Y_j.
%
%                     alphaAdjusted
%                         Bonferroni-adjusted component level.
%
%                     rejectPairs
%                         q-by-q logical matrix. Element (k,j) is true when
%                         the credible interval for the latent conditional
%                         covariance between R_k and Y_j excludes zero.
%
%                     whichReject
%                         q-by-1 logical vector. Element k is true when at
%                         least one partially observed outcome is associated
%                         with missingness indicator R_k.
%
%                     reject
%                         true if at least one off-diagonal component is
%                         rejected.
%
%                     nSaved
%                         number of retained posterior draws.
%
%                     pluginMarginal
%                         logical vector identifying margins whose observed
%                         latent values were kept fixed at normal scores.
%
%                     samples
%                         retained posterior draws when savesamples is true.
%                         Rows correspond to R_k and columns to Y_j.
%
% More About:
%
% The method models the joint distribution of the partially observed
% outcomes, their missingness indicators and the fully observed outcomes
% using a semiparametric Gaussian copula.
%
% The observed margins themselves are not assigned parametric
% distributions. Instead, each observed variable W_j is associated with a
% latent Gaussian variable Z_j, and the ordering of the observed values
% imposes ordering restrictions on the corresponding latent values. This is
% the extended rank-likelihood construction of Hoff (2007).
%
% The working matrix is ordered as
%
%        W = [Ymiss, R, Yfull],
%
% where:
%
%        Ymiss = Y(:,locmis),
%        R     = missingness indicators associated with Ymiss,
%        Yfull = Y(:,locfull).
%
% Therefore, the column blocks of W are:
%
%        columns 1:q       = partially observed outcomes;
%        columns q+1:2*q   = their missingness indicators;
%        remaining columns = fully observed outcomes.
%
% Under MAAR and the accompanying assumptions, for j different from k,
%
%        Y_j is conditionally independent of R_k given Yfull.
%
% The copula diagnostic investigates this implication on a latent Gaussian
% scale. At every retained MCMC iteration, the routine computes the
% conditional covariance matrix of
%
%        [Ymiss, R]
%
% given Yfull. The q-by-q block containing the conditional covariances
% between R and Ymiss is stored. Rows correspond to missingness indicators
% and columns to partially observed outcomes, as in the CCM output.
%
% A zero conditional covariance is equivalent to zero conditional
% correlation provided the conditional variances are positive. Consequently,
% testing conditional covariances rather than conditional correlations does
% not change the null hypothesis.
%
% The diagonal association between Y_j and its own missingness indicator
% R_j is not directly estimable from the observed data because Y_j is
% unavailable precisely for the units for which it is missing. The
% diagnostic therefore uses only the q*(q-1) off-diagonal associations.
%
% The original R implementation calls sbgcop::sbgcop.mcmc. The sampler below
% is a self-contained MATLAB implementation inspired by the same extended
% rank-likelihood approach.

% Number of variables having at least one missing value.
q = numel(locmis);

% Construct the working matrix used by the Gaussian-copula sampler.
%
% The ordering of the three blocks is important because the code below
% extracts the outcome-indicator block by its position in the latent
% covariance matrix.
W = [Y(:,locmis) R Y(:,locfull)];

% n is the number of units and d is the total number of margins in the
% copula model:
%
%        d = 2*q + number of fully observed variables.
[n,d] = size(W);

% Number of posterior draws expected after removing burn-in and applying
% thinning. The array used to store the draws is preallocated using this
% value.
nSaved = floor((niterMCMC-nburnMCMC)/thinMCMC);

% Very few retained draws generally produce unstable posterior quantiles.
if nSaved < 20
    warning('FSDA:mdMAARtest:FewCopulaDraws', ...
        'Only %d posterior draws will be retained.',nSaved);
end

% Z is the n-by-d matrix of latent Gaussian variables.
%
% levelIndex(i,j) records the ordered level occupied by observed value
% W(i,j). Observations tied in W receive the same level.
%
% nLevels(j) is the number of distinct observed levels of margin j.
%
% plugin(j) is true when the observed latent values of margin j are kept
% fixed at their initial normal scores.
Z = zeros(n,d);
levelIndex = zeros(n,d);
nLevels = zeros(1,d);
plugin = false(1,d);

% Initialize the latent Gaussian values one margin at a time.
for j = 1:d

    % Identify the observed entries of the current margin. Outcome margins
    % can contain NaNs, whereas the missingness-indicator margins are fully
    % observed binary variables.
    observed = ~isnan(W(:,j));

    % Extract the observed values of margin j.
    x = W(observed,j);

    % Convert the observed values into ordered integer levels.
    %
    % For example, values [10 20 20 40] receive levels [1 2 2 3].
    [~,~,lev] = unique(x,'sorted');
    levelIndex(observed,j) = lev;
    nLevels(j) = max(lev);

    % For a margin with many distinct values, regard it as approximately
    % continuous and keep its observed normal scores fixed during the
    % sampler. Missing latent values are still updated.
    plugin(j) = numel(unique(x)) > pluginthreshold;

    % Compute ranks assigning tied observations their maximum rank.
    rmax = rankMaximum(x);

    % Transform ranks to probabilities strictly inside the interval (0,1).
    % Division by n_j+1 prevents probabilities equal to exactly zero or one.
    u = rmax/(numel(x)+1);

    % Initialize observed latent values using standard normal quantiles.
    Z(observed,j) = normalInv(u);

    % Assign random starting values to latent entries corresponding to
    % missing outcome values. These are only initial values and are updated
    % during the MCMC iterations.
    Z(~observed,j) = randn(sum(~observed),1);
end

% Initialize the latent covariance matrix from the current latent scores.
% The second argument of cov equal to 1 uses normalization by n rather than
% n-1.
S = cov(Z,1);

% Ensure that the starting covariance matrix is symmetric positive
% definite before it is used in conditional Gaussian calculations.
S = stabilizeSPD(S,ridge);

% interestSamples(k,j,s) will contain retained draw s=1, 2, ..., nSaved, of
%
%   Cov(Z_Rk,Z_Yj | Z_Yfull),
%
% where k indexes missingness indicators and j indexes partially observed
% outcomes. This is the same orientation used by the CCM output matrices.
interestSamples = zeros(q,q,nSaved);

% Counter containing the number of posterior draws actually retained.
saveIndex = 0;

% Hyperparameters used for the inverse-Wishart covariance update.
%
% S0 is the prior scale matrix and n0 controls the prior degrees of freedom.
% These choices are implementation-specific and are not prescribed in
% detail in Section 3.4 of Bojinov, Pillai and Rubin.
S0 = eye(d);
n0 = d+2;

% Begin the Gibbs sampler.
for iteration = 1:niterMCMC

    % Update the latent margins in random order. Random ordering avoids
    % favouring the original column order systematically at every sweep.
    order = randperm(d);

    for jj = 1:d

        % j is the index of the margin updated at this step.
        j = order(jj);

        % Indices of all latent margins except margin j.
        other = [1:j-1 j+1:d];

        if isempty(other)
            % This branch is relevant only in the univariate case. The
            % conditional distribution is then the marginal Gaussian
            % distribution of Z_j.
            conditionalMean = zeros(n,1);
            conditionalSD = sqrt(max(S(j,j),eps));
        else
            % Extract and stabilize the covariance matrix of the remaining
            % latent margins.
            Soo = stabilizeSPD(S(other,other),ridge);

            % Regression coefficients in the conditional Gaussian
            % distribution:
            %
            %   S(j,-j) * inv(S(-j,-j)).
            regression = S(j,other)/Soo;

            % Conditional mean of Z_j for every observation:
            %
            %   E(Z_j | Z_-j)
            %       = Z_-j * inv(S_-j,-j) * S_-j,j.
            conditionalMean = Z(:,other)*regression';

            % Conditional variance:
            %
            %   Var(Z_j | Z_-j)
            %       = S(j,j)-S(j,-j)*inv(S(-j,-j))*S(-j,j).
            conditionalVariance = S(j,j)-regression*S(other,j);

            % Numerical rounding may produce a very small negative value.
            % Restrict the conditional variance to be at least eps.
            conditionalSD = sqrt(max(conditionalVariance,eps));
        end

        % Observed-data indicator for the current margin.
        observed = ~isnan(W(:,j));

        if ~plugin(j)
            % For non-plugin margins, resample the observed latent values
            % from their conditional normal distributions subject to the
            % rank restrictions implied by the observed data.

            for lev = 1:nLevels(j)

                % Rows whose observed value belongs to ordered level lev.
                idx = observed & levelIndex(:,j)==lev;

                if ~any(idx)
                    continue
                end

                % The lower truncation point is the largest latent value in
                % the immediately preceding observed level. The first level
                % has no finite lower bound.
                if lev == 1
                    lower = -Inf;
                else
                    lower = max(Z(observed & ...
                        levelIndex(:,j)==lev-1,j));
                end

                % The upper truncation point is the smallest latent value in
                % the immediately following observed level. The last level
                % has no finite upper bound.
                if lev == nLevels(j)
                    upper = Inf;
                else
                    upper = min(Z(observed & ...
                        levelIndex(:,j)==lev+1,j));
                end

                % Sample latent values from their conditional Gaussian
                % distributions while preserving the ordering of the
                % observed margin.
                Z(idx,j) = truncatedNormal(conditionalMean(idx), ...
                    conditionalSD,lower,upper);
            end
        end

        % Latent values corresponding to missing observed outcomes have no
        % rank restrictions. Sample them directly from their conditional
        % Gaussian distribution.
        %
        % Missingness-indicator columns R contain no NaNs, so this update is
        % relevant only to incomplete outcome margins.
        idxMissing = ~observed;

        if any(idxMissing)
            Z(idxMissing,j) = conditionalMean(idxMissing)+ ...
                conditionalSD*randn(sum(idxMissing),1);
        end
    end

    % Update the latent covariance matrix.
    %
    % Conditional on Z, the covariance matrix has an inverse-Wishart full
    % conditional distribution. The code samples the precision matrix from
    % the corresponding Wishart distribution and then inverts it.

    % Posterior scale contribution from the prior and latent cross-products.
    Psi = n0*S0+Z'*Z;

    % Scale matrix of the Wishart distribution for the precision matrix.
    invPsi = stableInverse(Psi,ridge);

    % Draw the latent precision matrix.
    precisionDraw = wishartRandom(invPsi,n0+n,ridge);

    % Convert the precision draw into the new covariance draw.
    S = stableInverse(precisionDraw,ridge);

    % Retain the current covariance draw only after burn-in and at the
    % requested thinning interval.
    if iteration > nburnMCMC && mod(iteration-nburnMCMC,thinMCMC)==0

        saveIndex = saveIndex+1;

        % Convert the latent covariance matrix S into the corresponding
        % correlation matrix C:
        %
        %   C(a,b)=S(a,b)/sqrt(S(a,a)*S(b,b)).
        sdS = sqrt(max(diag(S),eps));
        C = S./(sdS*sdS');

        % Remove minor numerical asymmetry.
        C = (C+C')/2;

        % The first 2*q variables are:
        %
        %   1:q       partially observed outcomes;
        %   q+1:2*q   their missingness indicators.
        %
        % The remaining variables, when present, are fully observed
        % outcomes on which the diagnostic conditions.
        nv = 2*q;

        if nv < d
            % C22 is the latent correlation block associated with the fully
            % observed variables.
            C22 = stabilizeSPD(C(nv+1:end,nv+1:end),ridge);

            % Compute the covariance of [Ymiss,R] conditional on Yfull
            %
            %   Vcond = C11-C12*inv(C22)*C21.
            Vcond = C(1:nv,1:nv)- ...
                C(1:nv,nv+1:end)/C22*C(nv+1:end,1:nv);
            % It removes from the
            % covariance among [Ymiss,R] the component explained through
            % their joint dependence on the fully observed variables.
            %
            % Equivalently, Vcond is the covariance matrix of the residuals
            % obtained after linearly regressing each latent variable in
            % [Z_Ymiss,Z_R] on Z_Yfull.

        else
            % If no fully observed variables are present, no conditioning
            % is possible. In this case the unconditional covariance of
            % [Ymiss,R] is used.
            Vcond = C(1:nv,1:nv);
        end

        % Extract the q-by-q indicator-outcome block:
        %
        %   rows    q+1:2*q   correspond to R;
        %   columns 1:q       correspond to Ymiss.
        %
        % Element (k,j) is the current posterior draw of the latent
        % conditional covariance between indicator R_k and outcome Y_j.
        interestSamples(:,:,saveIndex) = ...
            Vcond(q+1:2*q,1:q);
    end

    % Display progress approximately every 10 percent of the MCMC run.
    if msg && niterMCMC >= 10 && ...
            mod(iteration,max(1,floor(niterMCMC/10)))==0
        fprintf('mdMAARtest copula sampler: %d%% completed.\n', ...
            round(100*iteration/niterMCMC));
    end
end

% Remove any unused preallocated pages. The resulting array has dimensions
%
%             q-by-q-by-nSaved.
%
% Rows identify missingness indicators, columns identify partially
% observed outcomes and pages identify retained posterior draws.
interestSamples = interestSamples(:,:,1:saveIndex);
nSaved = saveIndex;

% Only off-diagonal outcome-indicator associations are used. There are
%
%             q*(q-1)
%
% such comparisons when q>1.
nTests = max(q^2-q,1);

% Apply a Bonferroni correction to obtain simultaneous component credible
% intervals.
alphaAdjusted = alpha/nTests;

% Equal-tail posterior probabilities for the lower and upper limits.
lowerProb = alphaAdjusted/2;
upperProb = 1-lowerProb;

% Preallocate matrices containing posterior summaries.
lower = zeros(q,q);
medianValue = zeros(q,q);
upper = zeros(q,q);
posteriorMean = zeros(q,q);

% The quantity stored in pvalue is a two-sided posterior tail probability:
%
%   2*min{Pr(delta<=0 | data),Pr(delta>=0 | data)}.
%
% It should not be interpreted as a classical frequentist p-value.
pvalue = ones(q,q);

% Compute posterior summaries separately for every indicator-outcome pair.
% As in CCM, rows correspond to missingness indicators R_k and columns to
% partially observed outcomes Y_j.
for k = 1:q
    for j = 1:q

        % Posterior draws for the conditional covariance between missingness
        % indicator k and outcome j.
        draws = squeeze(interestSamples(k,j,:));

        % Bonferroni-adjusted equal-tail credible interval and posterior
        % median.
        qq = quantile(draws,[lowerProb 0.5 upperProb]);

        lower(k,j) = qq(1);
        medianValue(k,j) = qq(2);
        upper(k,j) = qq(3);

        % Posterior mean of the conditional covariance.
        posteriorMean(k,j) = mean(draws);

        % Two-sided posterior sign probability. Small values indicate that
        % nearly all posterior mass lies on one side of zero.
        pvalue(k,j) = min(1, ...
            2*min(mean(draws<=0),mean(draws>=0)));
    end
end

% Flag a pair when its simultaneous posterior credible interval lies
% entirely above or entirely below zero.
rejectPairs = (lower>0 | upper<0);

if q > 1
    % Diagonal comparisons concern Y_j and its own missingness indicator
    % R_j. These associations are not treated as estimable component
    % diagnostics and are therefore excluded from the final decision.
    rejectPairs(1:q+1:end) = false;

    % Assign a neutral value to the diagonal posterior tail probabilities.
    pvalue(1:q+1:end) = 1;
end

% Row k of rejectPairs corresponds to missingness indicator R_k.
% Therefore, whichReject(k) is true when R_k is conditionally associated
% with at least one other partially observed outcome.
whichReject = any(rejectPairs,2);

% Store the posterior summaries and diagnostic decisions.
res = struct;

% Two-sided posterior tail probabilities.
res.pvalue = pvalue;

% For compatibility with the common output structure, stat contains the
% posterior median conditional covariances. Rows correspond to missingness
% indicators R_k and columns to partially observed outcomes Y_j.
res.stat = medianValue;

% Lower and upper credible limits. CI(:,:,1) contains lower limits and
% CI(:,:,2) contains upper limits, with rows corresponding to R_k and
% columns to Y_j.
res.CI = cat(3,lower,upper);

% Posterior location summaries with rows corresponding to R_k and columns
% to Y_j.
res.posteriorMean = posteriorMean;
res.posteriorMedian = medianValue;

% Component significance level after Bonferroni correction.
res.alphaAdjusted = alphaAdjusted;

% Pairwise and indicator-level rejection summaries.
res.rejectPairs = rejectPairs;
res.whichReject = whichReject;
res.reject = any(whichReject);

% Number of posterior draws retained after burn-in and thinning.
res.nSaved = nSaved;

% Logical vector identifying margins treated using fixed observed normal
% scores.
res.pluginMarginal = plugin;

% Store the complete posterior sample only when explicitly requested,
% because this q-by-q-by-nSaved array can be large.
if savesamples
    res.samples = interestSamples;
else
    res.samples = [];
end
end

% -------------------------------------------------------------------------
function ranks = rankMaximum(x)
%rankMaximum computes ranks with ties assigned their maximum rank.
%
% Input:
%   x     : numeric vector.
% Output:
%   ranks : vector having the same size as x. All observations in a tied
%           group receive the largest rank occupied by that group.
x = x(:);
[sorted,order] = sort(x);
ranksSorted = zeros(size(x));
first = 1;
while first <= numel(x)
    last = first;
    while last < numel(x) && sorted(last+1)==sorted(first)
        last = last+1;
    end
    ranksSorted(first:last) = last;
    first = last+1;
end
ranks = zeros(size(x));
ranks(order) = ranksSorted;
end

% -------------------------------------------------------------------------
function x = truncatedNormal(mu,sigma,lower,upper)
%truncatedNormal draws Gaussian variables subject to common truncation.
%
% Inputs:
%   mu    : vector of conditional means.
%   sigma : positive scalar conditional standard deviation.
%   lower : scalar lower truncation limit, possibly -Inf.
%   upper : scalar upper truncation limit, possibly Inf.
%
% Output:
%   x     : vector with size(mu), drawn from N(mu,sigma^2) and restricted
%           to [lower,upper].
%
% Inverse-transform sampling is used first. Extremely narrow or numerically
% degenerate intervals are handled by rejection sampling and a final safe
% fallback.
a = (lower-mu)/sigma;
b = (upper-mu)/sigma;
pLower = normalCDF(a);
pUpper = normalCDF(b);
width = pUpper-pLower;
u = pLower+rand(size(mu)).*width;
u = min(max(u,realmin),1-eps);
x = mu+sigma*normalInv(u);

bad = ~isfinite(x) | width <= 10*eps(max(pUpper,pLower));
if any(bad)
    badIndex = find(bad);
    for ii = 1:numel(badIndex)
        i = badIndex(ii);
        accepted = false;
        for attempt = 1:10000
            proposal = mu(i)+sigma*randn;
            if proposal >= lower && proposal <= upper
                x(i) = proposal;
                accepted = true;
                break
            end
        end
        if ~accepted
            lo = lower;
            hi = upper;
            if isinf(lo), lo = mu(i)-8*sigma; end
            if isinf(hi), hi = mu(i)+8*sigma; end
            x(i) = min(max(mu(i),lo+sqrt(eps)),hi-sqrt(eps));
        end
    end
end
end

% -------------------------------------------------------------------------
function p = normalCDF(x)
%normalCDF evaluates the standard normal cumulative distribution function.
p = 0.5*erfc(-x/sqrt(2));
end

% -------------------------------------------------------------------------
function x = normalInv(p)
%normalInv evaluates the standard normal quantile function. Probabilities
% are bounded away from zero and one to avoid infinite numerical values.
p = min(max(p,realmin),1-eps);
x = -sqrt(2)*erfcinv(2*p);
end

% -------------------------------------------------------------------------
function W = wishartRandom(V,nu,ridge)
%wishartRandom generates a random Wishart matrix by Bartlett decomposition.
%
% Inputs:
%   V     : positive-definite scale matrix.
%   nu    : Wishart degrees of freedom, not smaller than size(V,1).
%   ridge : regularization used to stabilize V before Cholesky factorization.
%
% Output:
%   W     : random matrix distributed as Wishart(V,nu), up to numerical
%           precision.
p = size(V,1);
if nu < p
    error('FSDA:mdMAARtest:WishartDF', ...
        'Wishart degrees of freedom must be at least the matrix dimension.');
end
V = stabilizeSPD(V,ridge);
L = chol(V,'lower');
A = zeros(p,p);
for i = 1:p
    A(i,i) = sqrt(2*randg((nu-i+1)/2));
    if i > 1
        A(i,1:i-1) = randn(1,i-1);
    end
end
W = L*(A*A')*L';
W = (W+W')/2;
end

% -------------------------------------------------------------------------
function A = stabilizeSPD(A,ridge)
%stabilizeSPD symmetrizes and, when necessary, regularizes a matrix.
%
% Inputs:
%   A     : square matrix expected to be symmetric positive definite.
%   ridge : nonnegative eigenvalue floor relative to the scale of A.
% Output:
%   A     : symmetric positive-definite version of the input matrix.
A = (A+A')/2;
if isempty(A)
    return
end
base = max(1,trace(abs(A))/size(A,1));
[~,flag] = chol(A);
if flag ~= 0 || rcond(A) < 1e-12
    [V,D] = eig(A,'vector');
    floorValue = max(ridge*base,eps(base));
    D = max(real(D),floorValue);
    A = V*diag(D)*V';
    A = (A+A')/2;
end
end

% -------------------------------------------------------------------------
function Ainv = stableInverse(A,ridge)
%stableInverse computes a symmetric inverse after SPD stabilization.
% A is the matrix to invert and ridge is passed to stabilizeSPD.
A = stabilizeSPD(A,ridge);
Ainv = A\eye(size(A));
Ainv = (Ainv+Ainv')/2;
end

% -------------------------------------------------------------------------
function p = fSurvival(x,df1,df2)
%fSurvival computes the upper-tail probability of an F distribution.
%
% Inputs:
%   x   : nonnegative F statistic.
%   df1 : numerator degrees of freedom.
%   df2 : denominator degrees of freedom, possibly Inf.
% Output:
%   p   : upper-tail probability.
%
% When df2 is infinite, df1*x has a chi-square distribution with df1
% degrees of freedom.
if isnan(x) || x < 0
    p = NaN;
elseif isinf(x)
    p = 0;
elseif isinf(df2)
    p = gammainc(df1*x/2,df1/2,'upper');
else
    z = df2/(df2+df1*x);
    p = betainc(z,df2/2,df1/2);
end
end

% -------------------------------------------------------------------------
function interpretation = interpretMAAR(out)
%interpretMAAR converts the numerical result into a concise verbal summary.
% out is the complete output structure and interpretation is a character
% vector suitable for display or inclusion in reports.
if out.reject
    implicated = out.missingVariables(out.whichReject);
    if isempty(implicated)
        interpretation = ['At least one component test rejects the restrictions ' ...
            'implied by MAAR and the accompanying assumptions.'];
    else
        interpretation = sprintf(['The diagnostic rejects the restrictions ' ...
            'implied by MAAR. Missingness indicators associated with variables %s ' ...
            'are implicated.'],mat2str(implicated(:)'));
    end
else
    interpretation = ['The diagnostic does not reject the restrictions implied ' ...
        'by MAAR. This result does not prove that the mechanism is MAAR.'];
end
end

% -------------------------------------------------------------------------
function printMAAR(out)
%printMAAR prints the main fields of the output structure in compact form.
fprintf('\nMAAR diagnostic: %s\n',upper(out.method));
fprintf('Observations: %d; variables: %d; partially observed variables: %d\n', ...
    out.n,out.p,out.nvarmiss);
fprintf('Nominal alpha: %.4g; component alpha: %.4g\n', ...
    out.alpha,out.alphaAdjusted);
fprintf('Overall rejection: %s\n',char(string(out.reject)));
if any(out.whichReject)
    fprintf('Implicated variable indices: %s\n', ...
        mat2str(out.missingVariables(out.whichReject)'));
else
    fprintf('No missingness indicator was implicated.\n');
end
fprintf('%s\n\n',out.interpretation);
end

% -------------------------------------------------------------------------
function plotMAAR(out)
%plotMAAR produces the method-specific graphical summary.
%
% Input:
%   out : output structure returned by mdMAARtest. The selected method
%         determines whether a p-value heatmap, a bubble scatter plot or a
%         posterior conditional-covariance heatmap is drawn.
q = out.nvarmiss;
labels = out.variableNames(out.missingVariables);
figure('Name',['mdMAARtest: ' upper(out.method)],'Color','w');

switch out.method

    case 'ccm'

        % Bonferroni-adjusted significance threshold.
        alphaBonf = out.alphaAdjusted;

        % Maximum p-value represented in the colour scale.
        maxPshown = 0.5;
        % Matrix of raw p-values from the pairwise conditional-mean
        % comparisons.
        values = out.pvalue;
        values(values > maxPshown) = maxPshown;

        % Comparisons on the diagonal are not performed.
        values(1:q+1:end) = NaN;

        % Create the heatmap.
        h = heatmap(labels,labels,values);

        h.XLabel = 'Partially observed response variable';
        h.YLabel = 'Missingness indicator';
        h.Title = sprintf(['CCM pairwise p-values; red cells are significant ' ...
            'after Bonferroni correction (p < %.3g)'],alphaBonf);
        h.ColorLimits = [0 maxPshown];
        h.CellLabelFormat = '%.3g';

        % Appearance of cells for which the test is not performed.
        h.MissingDataLabel = 'Not tested';
        h.MissingDataColor = [1 1 1];

        % Number of colours used in the heatmap.
        nColors = 1024;

        % Number of colour levels corresponding to significant p-values.
        nSignificantColors = ceil(alphaBonf/maxPshown*nColors);
        nSignificantColors = max(1,min(nColors-1,nSignificantColors));

        % Exact red used for all Bonferroni-significant cells.
        redSignificant = [0.88 0.08 0.08];
        significantMap = repmat(redSignificant,nSignificantColors,1);

        % Continuous scale for nonsignificant cells:
        % reddish-orange -> yellow -> light green -> dark green.
        anchorPositions = [0 0.30 0.65 1];

        anchorColors = [
            0.95 0.35 0.15   % reddish-orange: just above the threshold
            1.00 0.75 0.20   % orange-yellow
            0.65 0.85 0.35   % light green
            0.10 0.55 0.25   % green: large p-values
            ];

        nNonSignificantColors = nColors-nSignificantColors;
        queryPositions = linspace(0,1,nNonSignificantColors);

        nonSignificantMap = interp1(anchorPositions,anchorColors, ...
            queryPositions,'linear');

        % Low p-values are red; high p-values are green.
        h.Colormap = [significantMap; nonSignificantMap];

    case 'dtmm'
        % Display the two components entering the Meng-Rubin D_L statistic on
        % standardized, interpretable scales. The horizontal coordinate is
        % the likelihood-ratio numerator per restriction, d_L/q. The vertical
        % coordinate is lambda_L=r_L/(1+r_L), the fraction of missing
        % information corresponding to the relative increase r_L.

        % Horizontal coordinate: likelihood-ratio numerator per restriction.
        % xRaw contains the actual values; xPlot is used only for display.
        xRaw = out.DLtable.dL/q;
        rL = out.DLtable.rL;
        y = rL./(1+rL);

        % Keep a common zero origin and ensure that the theoretical reference
        % threshold is visible. Under r_L=0, d_L is asymptotically chi-square with
        % q degrees of freedom. Therefore the Bonferroni-adjusted reference
        % threshold for d_L/q is chi2_{q,1-alphaAdjusted}/q.
        xCritical = 2*gammaincinv(1-out.alphaAdjusted,q/2)/q;

        % Extremely large numerator values can compress the informative part
        % of the plot. Values exceeding six times the reference threshold are
        % therefore displayed at the graphical upper bound. Their exact values
        % remain available in out.DLtable.
        xCap = 6*xCritical;
        isCapped = xRaw > xCap;
        xPlot = min(xRaw,xCap);

        % Leave some space to the right of the most extreme displayed point.
        if any(isCapped)
            xMax = 1.10*xCap;
        else
            xMax = 1.12*max([xRaw(:); xCritical]);
        end

        % Encode the p-value through both bubble size and colour. The
        % transformation -log10(p) makes smaller p-values more prominent.
        % Truncate the graphical strength at 8 so that extremely small
        % p-values do not dominate the display.
        pStrength = -log10(max(out.DLtable.pvalue,realmin));
        pStrengthPlot = min(pStrength,8);
        bubbleSize = 50+45*pStrengthPlot;


        % Main bubble scatter.
        scatter(xPlot,y,bubbleSize,pStrengthPlot,'filled', ...
            'MarkerEdgeColor',[0.35 0.35 0.35])

        hold on

        % Give Bonferroni-significant components a black outer ring.
        significant = out.DLtable.pvalue < out.alphaAdjusted;
        if any(significant)
            scatter(xPlot(significant),y(significant), ...
                bubbleSize(significant)+25,'o', ...
                'MarkerEdgeColor','k','LineWidth',1.5)
        end

        % Colour map: blue for large p-values, through yellow/orange, to red
        % for very small p-values.
        nColors = 256;
        anchorPositions = [0 0.55 1];
        anchorColors = [
            0.20 0.55 0.95
            1.00 0.80 0.15
            0.85 0.05 0.05
            ];
        cmap = interp1(anchorPositions,anchorColors, ...
            linspace(0,1,nColors),'linear');
        colormap(gca,cmap)

        fs=16;

        xlim([0 xMax])
        ylim([0 1])
        xline(xCritical,'--','Bonferroni-adjusted threshold (r_L = 0)', ...
            'LabelVerticalAlignment','bottom', ...
            'LabelHorizontalAlignment','left','FontSize',fs);
        yline(0.5,'--','\lambda_L = 0.5  (r_L = 1)', ...
            'LabelVerticalAlignment','bottom', ...
            'LabelHorizontalAlignment','left','FontSize',fs);

        % Add labels to the four interpretative regions outside the main plotting
        % area. The reference lines are graphical guides only; the formal decision
        % is based on the finite-M Meng-Rubin p-value.
        ax = gca;

        text(ax,-0.02,-0.03, ...
            {'Low numerator','Low missing information'}, ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','top', ...
            'Clipping','off');

        text(ax,0.98,-0.03, ...
            {'High numerator','Low missing information'}, ...
            'Units','normalized', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','top', ...
            'Clipping','off');

        text(ax,-0.02,1.01, ...
            {'Low numerator','High missing information'}, ...
            'Units','normalized', ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom', ...
            'Clipping','off');

        text(ax,0.98,1.01, ...
            {'High numerator','High missing information'}, ...
            'Units','normalized', ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','bottom', ...
            'Clipping','off');

        % Label every bubble with the variable whose missingness indicator is
        % being tested. Use a small offset to avoid drawing text on the marker.
        dx = 0.012*xMax;
        dy = 0.015;


        text(xPlot+dx,y+dy,labels, ...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','left');

        % Mark points whose horizontal coordinate has been truncated.
        if any(isCapped)
            plot(xPlot(isCapped),y(isCapped),'>', ...
                'MarkerSize',8, ...
                'MarkerEdgeColor','k', ...
                'LineWidth',1.2, ...
                'MarkerFaceColor','none');
        end

        % Add a p-value colourbar. The colour variable is -log10(p), but the
        % displayed tick labels are p-values for immediate interpretation.

        % Use a fixed colour scale so that colours have the same meaning
        % across different data sets. The upper limit corresponds to p=1e-8.
        cMax = 8;
        clim([0 cMax])

        cb = colorbar;

        % Standard p-values used as reference labels on the colourbar.
        standardP = [1 0.1 0.05 0.01 0.001 1e-4 1e-6 1e-8];

        % Add the Bonferroni-adjusted significance level. To avoid
        % overlapping labels, remove standard ticks that are too close to
        % the Bonferroni threshold on the -log10(p) scale.
        alphaBonf = out.alphaAdjusted;
        alphaStrength = -log10(alphaBonf);

        standardStrength = -log10(standardP);

        % Minimum separation between colourbar tick labels on the
        % -log10(p) scale.
        minTickDistance = 0.35;

        keep = abs(standardStrength-alphaStrength) >= minTickDistance;

        candidateP = [standardP(keep) alphaBonf];
        candidateP = sort(candidateP,'descend');

        candidateStrength = -log10(candidateP);

        cb.Ticks = candidateStrength;
        cb.Label.String = 'p-value';

        tickLabels = compose('%.3g',candidateP);

        % Values at or below 1e-8 receive the same maximum colour.
        idxMin = candidateP == 1e-8;
        tickLabels(idxMin) = {'<=1e-8'};

        cb.TickLabels = tickLabels;

        grid on
        box on
        xlabel('Likelihood-ratio improvement per restriction, d_L/q','FontSize',fs)
        ylabel('Fraction of missing information, \lambda_L = r_L/(1+r_L)','FontSize',fs)
        title({'Direct missingness-mechanism diagnostic', ...
            sprintf(['Bubble size and colour increase as p decreases; ' ...
            'black ring: p < \\alpha_{\\rm Bonf} = %.4g'], ...
            out.alphaAdjusted)})
        hold off

    case 'cop'

        % The copula output matrices use the same orientation as CCM:
        %
        %   row k    = missingness indicator R_k;
        %   column j = partially observed outcome Y_j.
        %
        % Thus, out.pvalue(k,j) summarizes the posterior evidence about
        % the latent conditional association between R_k and Y_j after
        % conditioning on the fully observed variables.
        %
        % A red cell in displayed row k and column j indicates posterior
        % evidence of a nonzero latent conditional association between Y_j
        % and R_k, after conditioning on the fully observed variables.

        % Bonferroni-adjusted component threshold used for the simultaneous
        % posterior credible intervals.
        alphaBonf = out.alphaAdjusted;

        % Maximum posterior tail probability represented in the colour
        % scale. Larger values receive the same green colour.
        maxPshown = 0.5;

        % out.pvalue already has missingness indicators along the rows and
        % partially observed outcomes along the columns, as in the CCM
        % heatmap.
        values = out.pvalue;

        % Truncate large values only for graphical purposes.
        values(values > maxPshown) = maxPshown;

        % Diagonal associations between an outcome and its own missingness
        % indicator are not used by the diagnostic.
        values(1:q+1:end) = NaN;

        % In the displayed heatmap, columns correspond to partially observed
        % outcomes and rows correspond to missingness indicators.
        h = heatmap(labels,labels,values);

        h.XLabel = 'Partially observed variable Y_j';
        h.YLabel = 'Missingness indicator R_k';

        h.Title = sprintf([ ...
            'Copula posterior tail probabilities for latent conditional ' ...
            'associations; red indicates Bonferroni-adjusted rejection ' ...
            '(value < %.3g)'],alphaBonf);

        h.ColorLimits = [0 maxPshown];
        h.CellLabelFormat = '%.3g';

        % Appearance of diagonal cells, for which no component diagnostic
        % is reported.
        h.MissingDataLabel = 'Not tested';
        h.MissingDataColor = [1 1 1];

        % Number of colours used in the heatmap.
        nColors = 1024;

        % Number of colour levels associated with values below the
        % Bonferroni-adjusted threshold.
        nSignificantColors = ceil(alphaBonf/maxPshown*nColors);
        nSignificantColors = ...
            max(1,min(nColors-1,nSignificantColors));

        % Exact red used for significant components.
        redSignificant = [0.88 0.08 0.08];
        significantMap = repmat(redSignificant, ...
            nSignificantColors,1);

        % Continuous scale for nonsignificant components:
        % reddish-orange -> yellow -> light green -> dark green.
        anchorPositions = [0 0.30 0.65 1];

        anchorColors = [
            0.95 0.35 0.15   % reddish-orange: just above threshold
            1.00 0.75 0.20   % orange-yellow
            0.65 0.85 0.35   % light green
            0.10 0.55 0.25   % dark green: large tail probability
            ];

        nNonSignificantColors = nColors-nSignificantColors;
        queryPositions = linspace(0,1,nNonSignificantColors);

        nonSignificantMap = interp1(anchorPositions,anchorColors, ...
            queryPositions,'linear');

        % Small posterior tail probabilities are red; large values are
        % green.
        h.Colormap = [significantMap; nonSignificantMap];
end
end

%FScategory:MULT-MissingData
