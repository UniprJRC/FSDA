% FSDA
%
% File names, description, category and date last modified
%
%   Name                       - Description                                                                                                                     - Category            - Date last modified
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   add2spm                    - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                   - VIS-Mult            - 2016 May 09
%   addt                       - Produces the t test for an additional expl. variable                                                                            - REG-Regression      - 2016 May 23
%   basicPower                 - Computes the basic power transformation                                                                                         - UTISTAT             - 2016 May 17
%   bc                         - Returns the Binomial coefficient                                                                                                - UTICOMB             - 2016 May 16
%   boxplotb                   - Computes a bivariate boxplot                                                                                                    - VIS-Mult            - 2016 May 23
%   brushFAN                   - Displays a GUI which enables brushing in the fanplot                                                                            - GUI                 - 2016 May 17
%   brushRES                   - Displays a GUI which enables brushing in resfwdplot                                                                             - GUI                 - 2016 May 17
%   brushROB                   - Displays a GUI which enables brushing in resindexplot                                                                           - GUI                 - 2016 May 17
%   cabc                       - Closes all open figures except the one in foreground (the current)                                                              - UTIGEN              - 2016 May 16
%   cascade                    - Is a third party function used in FSDA demos and examples                                                                       - UTIGEN              - 2016 May 16
%   cdsplot                    - Produces the candlestick plot for robust model selection in linear regression                                                   - VIS-Reg             - 2016 May 25
%   clickableMultiLegend       - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                             - UTIGEN              - 2016 May 22
%   combsFS                    - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                      - UTICOMB             - 2016 May 09
%   covplot                    - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                         - VIS-Mult            - 2016 May 23
%   ellipse                    - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT             - 2016 May 17
%   fanplot                    - Plots the fan plot for transformation in linear regression                                                                      - VIS-Reg             - 2016 May 23
%   findDir                    - Finds recursively all directories in root                                                                                       - UTIGEN              - 2016 May 11
%   findFile                   - Finds recursively all files in root                                                                                             - UTIGEN              - 2016 May 09
%   FSM                        - Gives an automatic outlier detection procedure in mult. analysis                                                                - MULT-Multivariate   - 2016 May 25
%   FSMbbm                     - Gives the units belonging to subset at step(s) msel of the forward search                                                       - MULT-Multivariate   - 2016 May 09
%   FSMbonfbound               - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT             - 2016 May 25
%   FSMbsb                     - Gives the units belonging to subset at step(s) msel of the forward search                                                       - MULT-Multivariate   - 2016 May 20
%   FSMeda                     - Performs forward search in multivariate analysis with exploratory data analysis purposes                                        - MULT-Multivariate   - 2016 May 09
%   FSMenvmmd                  - Computes the theoretical envelopes of Minimum MD outside subset during the search                                               - MULT-Multivariate   - 2016 May 26
%   FSMfan                     - Computes confirmatory lrt of a suggested transformation                                                                         - MULT-Transformations- 2016 May 23
%   FSMinvmmd                  - Converts values of minimum Mahalanobis distance into confidence levels                                                          - MULT-Multivariate   - 2016 May 25
%   FSMmmd                     - Monitors minMD                                                                                                                  - MULT-Multivariate   - 2016 May 16
%   FSMmmdeasy                 - Is exactly equal to minMD but much less efficient                                                                               - MULT-Multivariate   - 2016 May 16
%   FSMmmdrs                   - Performs random start monitoring of minimum Mahalanobis distance                                                                - CLUS-RobClaMULT     - 2016 May 25
%   FSMtra                     - Computes MLE of transformation parameters                                                                                       - MULT-Transformations- 2016 May 25
%   FSR                        - Gives an automatic outlier detection procedure in linear regression                                                             - REG-Regression      - 2016 May 25
%   FSRaddt                    - Produces t deletion tests for each expl. variable                                                                               - REG-ModelSelection  - 2016 May 09
%   FSRB                       - Gives an automatic outlier detection procedure in Bayesian linear regression                                                    - REG-Bayes           - 2016 May 25
%   FSRBbsb                    - Returns the units belonging to the subset in each step of the Bayesian forward search                                           - REG-Regression      - 2016 May 25
%   FSRBeda                    - Enables to monitor several quantities in each step of the Bayesian search                                                       - REG-Bayes           - 2016 May 25
%   FSRBmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search             - REG-Bayes           - 2016 May 25
%   FSRbonfbound               - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT             - 2016 May 09
%   FSRBr                      - Bayesian forward search in linear regression reweighted                                                                         - REG-Bayes           - 2016 May 25
%   FSRbsb                     - Returns the units belonging to the subset in each step of the forward search                                                    - REG-Regression      - 2016 May 25
%   FSRcp                      - Monitors Cp and AIC for all models of interest of size smallp                                                                   - REG-ModelSelection  - 2016 May 09
%   FSReda                     - Enables to monitor several quantities in each step of the forward search                                                        - REG-Regression      - 2016 May 25
%   FSRenvmdr                  - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                - REG-Regression      - 2016 May 17
%   FSRfan                     - Monitors the values of the score test statistic for each lambda                                                                 - REG-Transformations - 2016 May 17
%   FSRH                       - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                             - REG-Hetero          - 2016 May 26
%   FSRHbsb                    - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                    - REG-Hetero          - 2016 May 25
%   FSRHeda                    - Enables to monitor several quantities in each step of the forward search                                                        - REG-Hetero          - 2016 May 23
%   FSRHmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search      - REG-Hetero          - 2016 May 26
%   FSRinvmdr                  - Converts values of minimum deletion residual into confidence levels                                                             - REG-Regression      - 2016 May 25
%   FSRmdr                     - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                      - REG-Regression      - 2016 May 25
%   FSRmdrrs                   - Performs random start monitoring of minimum deletion residual                                                                   - CLUS-RobClaREG      - 2016 May 25
%   FSRms                      - Performs robust model selection using flexible trimming in linear regression                                                    - REG-ModelSelection  - 2016 May 25
%   FSRr                       - Forward search in linear regression reweighted                                                                                  - REG-Regression      - 2016 May 09
%   HAbdp                      - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2016 May 09
%   HAc                        - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT             - 2016 May 09
%   HAeff                      - Finds the tuning constant guarrantees a requested asymptotic efficiency                                                         - UTISTAT             - 2016 May 25
%   HApsi                      - Computes psi function  using Hampel proposal                                                                                    - UTISTAT             - 2016 May 25
%   HApsider                   - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT             - 2016 May 25
%   HApsix                     - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT             - 2016 May 25
%   HArho                      - Computes rho function  using Hampel proposal                                                                                    - UTISTAT             - 2016 May 25
%   HAwei                      - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT             - 2016 May 17
%   histFS                     - Plots a histogram with the elements in each bin grouped according to a vector of labels                                         - VIS-Reg             - 2016 May 25
%   HUeff                      - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT             - 2016 May 17
%   HUpsi                      - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT             - 2016 May 11
%   HUpsider                   - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT             - 2016 May 17
%   HUpsix                     - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT             - 2016 May 09
%   HUrho                      - Computes (rho) function for Huber                                                                                               - UTISTAT             - 2016 May 17
%   HUwei                      - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT             - 2016 May 17
%   HYPbdp                     - Finds constant c which is associated to the requested breakdown                                                                 - UTISTAT             - 2016 May 25
%   HYPc                       - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT             - 2016 May 25
%   HYPck                      - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT             - 2016 May 25
%   HYPeff                     - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT             - 2016 May 25
%   HYPk                       - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT             - 2016 May 26
%   HYPpsi                     - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT             - 2016 May 25
%   HYPpsider                  - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT             - 2016 May 25
%   HYPpsix                    - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT             - 2016 May 25
%   HYPrho                     - Computes rho function  using hyperboloc tangent estimator                                                                       - UTISTAT             - 2016 May 25
%   HYPwei                     - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT             - 2016 May 25
%   inversegamcdf              - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT             - 2016 May 09
%   inversegaminv              - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT             - 2016 May 09
%   inversegampdf              - Inverse-gamma probability density function                                                                                      - UTISTAT             - 2016 May 09
%   isfunction                 - Checks if a function exists                                                                                                     - UTIGEN              - 2016 May 09
%   levfwdplot                 - Plots the trajectories of leverage along the search                                                                             - VIS-Reg             - 2016 May 25
%   lexunrank                  - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                    - UTICOMB             - 2016 May 16
%   lga                        - Performs linear grouping analysis                                                                                               - CLUS-RobClaREG      - 2016 May 25
%   logmvnpdfFS                - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT             - 2016 May 26
%   LXS                        - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                            - REG-Regression      - 2016 May 09
%   mahalFS                    - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT             - 2016 May 17
%   makecontentsfileFS         - Extends Matlab function makecontentsfile                                                                                        - UTIHELP             - 2016 May 26
%   malfwdplot                 - Plots the trajectories of scaled Mahalanobis distances along the search                                                         - VIS-Mult            - 2016 May 22
%   malindexplot               - Plots the Mahalanobis distances versus a selected variable                                                                      - VIS-Mult            - 2016 May 25
%   mcd                        - Computes Minimum Covariance Determinant                                                                                         - MULT-Multivariate   - 2016 May 25
%   mdrplot                    - Plots the trajectory of minimum deletion residual (mdr)                                                                         - VIS-Reg             - 2016 May 23
%   MixSim                     - Generates k clusters in v dimensions with given overlap                                                                         - CLUS-MixSim         - 2016 May 25
%   MixSimreg                  - Generates k regression hyperplanes in p dimensions with given overlap                                                           - CLUS-MixSim         - 2016 May 26
%   mmdplot                    - Plots the trajectory of minimum mhalanobis distance (mmd)                                                                       - VIS-Mult            - 2016 May 26
%   MMmult                     - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                          - MULT-Multivariate   - 2016 May 25
%   MMmultcore                 - Computes multivariate MM estimators for a selected fixed scale                                                                  - MULT-Multivariate   - 2016 May 25
%   MMreg                      - Computes MM estimator of regression coefficients                                                                                - REG-Regression      - 2016 May 25
%   MMregcore                  - Computes MM regression estimators for a selected fixed scale                                                                    - REG-Regression      - 2016 May 25
%   Mscale                     - Finds the M estimator of the scale                                                                                              - UTISTAT             - 2016 May 26
%   mve                        - Computes Minimum volume ellipsoid                                                                                               - MULT-Multivariate   - 2016 May 20
%   nchoosekFS                 - Returns the Binomial coefficient or matrix containing all combinations                                                          - UTICOMB             - 2016 May 16
%   ncx2mixtcdf                - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT             - 2016 May 09
%   normBoxCox                 - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT             - 2016 May 17
%   normYJ                     - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT             - 2016 May 17
%   OPTbdp                     - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2016 May 09
%   OPTc                       - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT             - 2016 May 25
%   OPTeff                     - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT             - 2016 May 25
%   OPTpsi                     - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT             - 2016 May 25
%   OPTpsider                  - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT             - 2016 May 25
%   OPTpsix                    - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT             - 2016 May 25
%   OPTrho                     - Computes rho function for optimal weight function                                                                               - UTISTAT             - 2016 May 25
%   OPTwei                     - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT             - 2016 May 25
%   overlap                    - Computes the exact overlap given the parameters of the mixture                                                                  - CLUS-MixSim         - 2016 May 25
%   PoolClose                  - Closes the pool of MATLAB instances opened with PoolPrepare to execute code in parallel                                         - UTIGEN              - 2016 May 09
%   PoolPrepare                - Prepares a pool of MATLAB instances for executing code in parallel                                                              - UTIGEN              - 2016 May 09
%   position                   - Controls the position of the open figures                                                                                       - UTIGEN              - 2016 May 16
%   Powertra                   - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT             - 2016 May 17
%   publishFS                  - Enables to create automatic HELP FILES from structured .m function files                                                        - UTIHELP             - 2016 May 26
%   publishFunctionAlpha       - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 May 26
%   publishFunctionCate        - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 May 23
%   Qn                         - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT             - 2016 May 25
%   quickselectFS              - Finds the k-th order statistic                                                                                                  - UTIGEN              - 2016 May 25
%   RandIndexFS                - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT             - 2016 May 25
%   randsampleFS               - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                         - UTICOMB             - 2016 May 16
%   regressB                   - Computes Bayesian estimates of regression parameters                                                                            - REG-Bayes           - 2016 May 25
%   regressH                   - Fits a multiple linear regression model with heteroskedasticity                                                                 - REG-Hetero          - 2016 May 25
%   regressHart                - Fits a multiple linear regression model using art heteroskedasticity                                                            - REG-Hetero          - 2016 May 26
%   regressHhar                - RegressH fits a multiple linear regression model with Harvey heteroskedasticity                                                 - REG-Hetero          - 2016 May 17
%   resfwdplot                 - Plots the trajectories of the scaled (squared) residuals monitored                                                              - VIS-Reg             - 2016 May 23
%   resindexplot               - Plots the residuals from a regression analysis versus index number or any other variable                                        - VIS-Reg             - 2016 May 25
%   rlga                       - Performs robust linear grouping analysis                                                                                        - CLUS-RobClaREG      - 2016 May 26
%   RobCov                     - Computes covariance matrix of robust regression coefficients                                                                    - REG-Regression      - 2016 May 25
%   RobRegrSize                - Provides proper threshold for robust estimators to obtain an empirical size equal to 1 per cent nominal size                    - REG-Regression      - 2016 May 25
%   Score                      - Computes the score test for transformation                                                                                      - REG-Transformations - 2016 May 09
%   SDest                      - Computes Stahel-Donoho robust estimator of dispersion/location                                                                  - MULT-Multivariate   - 2016 May 25
%   shuffling                  - Does a random permutation of the elements of input vector                                                                       - UTICOMB             - 2016 May 16
%   Smult                      - Computes S estimators in multivariate analysis                                                                                  - MULT-Multivariate   - 2016 May 25
%   Sn                         - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT             - 2016 May 25
%   spmplot                    - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal                                     - VIS-Mult            - 2016 May 25
%   Sreg                       - Computes S estimators in linear regression                                                                                      - REG-Regression      - 2016 May 25
%   subsets                    - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                - UTICOMB             - 2016 May 09
%   tabulateFS                 - Create frequency table of unique values of x, excluding possible 0 counts                                                       - UTISTAT             - 2016 May 09
%   Taureg                     - Computes Tau estimators in linear regression                                                                                    - REG-Regression      - 2016 May 25
%   TBbdp                      - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT             - 2016 May 25
%   TBc                        - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT             - 2016 May 25
%   TBeff                      - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT             - 2016 May 25
%   TBpsi                      - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT             - 2016 May 25
%   TBpsider                   - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT             - 2016 May 25
%   TBpsix                     - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT             - 2016 May 25
%   TBrho                      - Computes (rho) function for Tukey biweight                                                                                      - UTISTAT             - 2016 May 25
%   TBwei                      - Computes weight function psi(u)/u for Tukey biweight                                                                            - UTISTAT             - 2016 May 25
%   tclust                     - Computes trimmed clustering with restricitons on the eigenvalues                                                                - CLUS-RobClaMULT     - 2016 May 18
%   tclustIC                   - Computes tclust for different number of groups k and restriction factors c                                                      - CLUS-RobClaMULT     - 2016 May 26
%   tclustICplot               - Plots information criterion as a function of c and k                                                                            - VIS-Clu             - 2016 May 25
%   tclustICsol                - Extracts a set of best relevant solutions                                                                                       - CLUS-RobClaMULT     - 2016 May 25
%   tclustreg                  - Performs robust linear grouping analysis                                                                                        - CLUS-RobClaREG      - 2016 May 19
%   triu2vec                   - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                       - UTIGEN              - 2016 May 16
%   unibiv                     - Has the purpose of detecting univariate and bivariate outliers                                                                  - MULT-Multivariate   - 2016 May 25
%   upperfracpos               - Positions two figures on the upper part of the screen                                                                           - UTIGEN              - 2016 May 25
%   winsor                     - Returns a winsorized copy of input                                                                                              - UTISTAT             - 2016 May 26
%   yXplot                     - Produces an interactive scatterplot of y against each variable of X in the input dataset                                        - VIS-Reg             - 2016 May 25
%   zscoreFS                   - Computes (robust) standardized z scores                                                                                         - UTIGEN              - 2016 May 25
