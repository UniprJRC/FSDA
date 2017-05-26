% FSDA
%
% File names, description, category and date last modified
%
%   Name                       - Description                                                                                                                     - Category            - Date last modified
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   add2spm                    - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                   - VIS-Mult            - 2017 Mar 07
%   addt                       - Produces the t test for an additional explanatory variable                                                                      - REG-Regression      - 2017 May 26
%   basicPower                 - Computes the basic power transformation                                                                                         - UTISTAT             - 2017 Mar 02
%   bc                         - Returns the Binomial coefficient                                                                                                - UTICOMB             - 2016 Oct 01
%   boxplotb                   - Computes a bivariate boxplot                                                                                                    - VIS-Mult            - 2016 Aug 01
%   brushFAN                   - Displays a GUI which enables brushing in the fanplot                                                                            - GUI                 - 2017 May 17
%   brushRES                   - Displays a GUI which enables brushing in resfwdplot                                                                             - GUI                 - 2017 May 17
%   brushROB                   - Displays a GUI which enables brushing in resindexplot                                                                           - GUI                 - 2017 May 17
%   bwe                        - Estimates the bandwidth smoothing parameter for kernel density estimation                                                       - UTISTAT             - 2017 May 12
%   cabc                       - Closes all open figures except the one in foreground (the current)                                                              - UTIGEN              - 2016 Oct 01
%   carbikeplot                - Produces the carbike plot to find best relevant clustering solutions                                                            - CLUS-RobClaMULT     - 2017 May 12
%   cascade                    - Is a third party function used in FSDA demos and examples                                                                       - UTIGEN              - 2016 Aug 01
%   cdsplot                    - Produces the candlestick plot for robust model selection in linear regression                                                   - VIS-Reg             - 2016 Aug 01
%   clickableMultiLegend       - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                             - UTIGEN              - 2016 Sep 26
%   combsFS                    - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                      - UTICOMB             - 2017 May 26
%   covplot                    - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                         - VIS-Mult            - 2016 Aug 01
%   ellipse                    - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT             - 2016 Aug 01
%   fanplot                    - Plots the fan plot for transformation in linear regression                                                                      - VIS-Reg             - 2017 May 15
%   findDir                    - Finds recursively all directories in root                                                                                       - UTIGEN              - 2017 May 15
%   findFile                   - Finds recursively all files in root                                                                                             - UTIGEN              - 2016 Aug 01
%   FowlkesMallowsIndex        - Computes the Fowlkes and Mallows index                                                                                          - UTISTAT             - 2017 May 15
%   FSM                        - Gives an automatic outlier detection procedure in mult. analysis                                                                - MULT-Multivariate   - 2017 Apr 20
%   FSMbonfbound               - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT             - 2017 May 12
%   FSMbsb                     - Gives the units belonging to subset at step(s) msel of the forward search                                                       - MULT-Multivariate   - 2016 Aug 01
%   FSMeda                     - Performs forward search in multivariate analysis with exploratory data analysis purposes                                        - MULT-Multivariate   - 2017 Apr 20
%   FSMenvmmd                  - Computes the theoretical envelopes of Minimum MD outside subset during the search                                               - MULT-Multivariate   - 2017 May 15
%   FSMfan                     - Computes confirmatory lrt of a suggested transformation                                                                         - MULT-Transformations- 2017 May 12
%   FSMinvmmd                  - Converts values of minimum Mahalanobis distance into confidence levels                                                          - MULT-Multivariate   - 2016 Aug 01
%   FSMmmd                     - Monitors minMD                                                                                                                  - MULT-Multivariate   - 2016 Dec 21
%   FSMmmdeasy                 - Is exactly equal to FSMmmd but it is much less efficient                                                                        - MULT-Multivariate   - 2017 May 12
%   FSMmmdrs                   - Performs random start monitoring of minimum Mahalanobis distance                                                                - CLUS-RobClaMULT     - 2016 Aug 01
%   FSMtra                     - Computes MLE of transformation parameters                                                                                       - MULT-Transformations- 2017 May 12
%   FSR                        - Gives an automatic outlier detection procedure in linear regression                                                             - REG-Regression      - 2017 May 26
%   FSRaddt                    - Produces t deletion tests for each explanatory variable                                                                         - REG-ModelSelection  - 2017 May 26
%   FSRB                       - Gives an automatic outlier detection procedure in Bayesian linear regression                                                    - REG-Bayes           - 2017 May 26
%   FSRBbsb                    - Returns the units belonging to the subset in each step of the Bayesian forward search                                           - REG-Regression      - 2017 May 26
%   FSRBeda                    - Enables to monitor several quantities in each step of the Bayesian search                                                       - REG-Bayes           - 2017 May 26
%   FSRBmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search             - REG-Bayes           - 2017 May 26
%   FSRbonfbound               - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT             - 2017 May 15
%   FSRBr                      - Bayesian forward search in linear regression reweighted                                                                         - REG-Bayes           - 2017 May 26
%   FSRbsb                     - Returns the units belonging to the subset in each step of the forward search                                                    - REG-Regression      - 2017 May 26
%   FSRcp                      - Monitors Cp and AIC for all models of interest of size smallp                                                                   - REG-ModelSelection  - 2017 May 26
%   FSReda                     - Enables to monitor several quantities in each step of the forward search                                                        - REG-Regression      - 2017 May 26
%   FSRenvmdr                  - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                - REG-Regression      - 2017 Feb 17
%   FSRfan                     - Monitors the values of the score test statistic for each lambda                                                                 - REG-Transformations - 2017 May 26
%   FSRH                       - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                             - REG-Hetero          - 2017 May 26
%   FSRHbsb                    - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                    - REG-Hetero          - 2017 May 26
%   FSRHeda                    - Enables to monitor several quantities in each step of the forward search                                                        - REG-Hetero          - 2017 May 26
%   FSRHmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search      - REG-Hetero          - 2017 May 26
%   FSRinvmdr                  - Converts values of minimum deletion residual into confidence levels                                                             - REG-Regression      - 2016 Aug 01
%   FSRmdr                     - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                      - REG-Regression      - 2017 May 26
%   FSRmdrrs                   - Performs random start monitoring of minimum deletion residual                                                                   - CLUS-RobClaREG      - 2017 May 26
%   FSRms                      - Performs robust model selection using flexible trimming in linear regression                                                    - REG-ModelSelection  - 2017 May 12
%   FSRr                       - Forward search in linear regression reweighted                                                                                  - REG-Regression      - 2017 Apr 20
%   HAbdp                      - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2016 Aug 01
%   HAc                        - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT             - 2016 Aug 01
%   HAeff                      - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT             - 2017 May 17
%   HApsi                      - Computes psi function  using Hampel proposal                                                                                    - UTISTAT             - 2017 May 17
%   HApsider                   - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT             - 2017 May 17
%   HApsix                     - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT             - 2017 May 17
%   HArho                      - Computes rho function  using Hampel proposal                                                                                    - UTISTAT             - 2017 May 17
%   HAwei                      - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT             - 2017 May 17
%   histFS                     - Plots a histogram with the elements in each bin grouped according to a vector of labels                                         - VIS-Reg             - 2016 Aug 01
%   htmlwriteFS                - Enables to create automatic HELP FILES from a specific MATLAB structure created with function mreadFS.m                         - UTIHELP             - 2017 May 12
%   HUeff                      - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT             - 2017 May 17
%   HUpsi                      - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT             - 2017 May 17
%   HUpsider                   - Computes derivative of psi function (second derivative of rho function) for Huber                                               - UTISTAT             - 2017 May 17
%   HUpsix                     - Computes psi function (derivative of rho function) times x for Huber                                                            - UTISTAT             - 2017 May 17
%   HUrho                      - Computes rho function for Huber                                                                                                 - UTISTAT             - 2016 Aug 01
%   HUwei                      - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT             - 2016 Aug 01
%   HYPbdp                     - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT             - 2017 May 17
%   HYPc                       - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT             - 2016 Aug 01
%   HYPck                      - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT             - 2016 Oct 02
%   HYPeff                     - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT             - 2017 May 17
%   HYPk                       - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT             - 2016 Aug 01
%   HYPpsi                     - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT             - 2017 May 17
%   HYPpsider                  - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT             - 2017 May 17
%   HYPpsix                    - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT             - 2016 Aug 01
%   HYPrho                     - Computes rho function  using hyperbolic tangent estimator                                                                       - UTISTAT             - 2017 May 17
%   HYPwei                     - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT             - 2017 May 17
%   inversegamcdf              - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT             - 2016 Aug 01
%   inversegaminv              - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT             - 2016 Aug 01
%   inversegampdf              - Computes inverse-gamma probability density function                                                                             - UTISTAT             - 2016 Aug 01
%   isfunction                 - Checks if a function exists                                                                                                     - UTIGEN              - 2016 Aug 01
%   kdebiv                     - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                  - UTISTAT             - 2017 May 15
%   levfwdplot                 - Plots the trajectories of leverage along the search                                                                             - VIS-Reg             - 2017 Feb 14
%   lexunrank                  - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                    - UTICOMB             - 2016 Aug 01
%   lga                        - Performs linear grouping analysis                                                                                               - CLUS-RobClaREG      - 2017 May 22
%   logmvnpdfFS                - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT             - 2016 Aug 01
%   LXS                        - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                            - REG-Regression      - 2017 May 26
%   mahalFS                    - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT             - 2017 May 02
%   makecontentsfileFS         - Extends Matlab function makecontentsfile                                                                                        - UTIHELP             - 2017 May 26
%   malfwdplot                 - Plots the trajectories of scaled Mahalanobis distances along the search                                                         - VIS-Mult            - 2017 May 12
%   malindexplot               - Plots the Mahalanobis distances versus a selected variable                                                                      - VIS-Mult            - 2017 May 12
%   mcd                        - Computes Minimum Covariance Determinant                                                                                         - MULT-Multivariate   - 2017 May 12
%   mdrplot                    - Plots the trajectory of minimum deletion residual (mdr)                                                                         - VIS-Reg             - 2017 Feb 14
%   MixSim                     - Generates k clusters in v dimensions with given overlap                                                                         - CLUS-MixSim         - 2017 May 12
%   MixSimreg                  - Generates k regression hyperplanes in p dimensions with given overlap                                                           - CLUS-MixSim         - 2017 May 12
%   mmdplot                    - Plots the trajectory of minimum Mahalanobis distance (mmd)                                                                      - VIS-Mult            - 2017 May 12
%   MMmult                     - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                          - MULT-Multivariate   - 2017 May 12
%   MMmultcore                 - Computes multivariate MM estimators for a selected fixed scale                                                                  - MULT-Multivariate   - 2016 Aug 01
%   MMmulteda                  - Computes MM estimators in multivariate analysis for a series of values of bdp                                                   - MULT-Multivariate   - 2017 May 04
%   MMreg                      - Computes MM estimator of regression coefficients                                                                                - REG-Regression      - 2017 May 26
%   MMregcore                  - Computes MM regression estimators for a selected fixed scale                                                                    - REG-Regression      - 2017 May 26
%   MMregeda                   - Computes MM estimator in linear regression for a series of values of efficiency                                                 - REG-Regression      - 2017 May 26
%   mreadFS                    - Enables to create a structure with InputArgs/OptArgs/OutArgs ... from .m function files                                         - UTIHELP             - 2017 May 26
%   Mscale                     - Finds the M estimator of the scale                                                                                              - UTISTAT             - 2016 Aug 01
%   mve                        - Computes Minimum volume ellipsoid                                                                                               - MULT-Multivariate   - 2017 May 12
%   mveeda                     - Mve monitors Minimum volume ellipsoid for a series of values of bdp                                                             - MULT-Multivariate   - 2017 May 12
%   nchoosekFS                 - Returns the Binomial coefficient or matrix containing all combinations                                                          - UTICOMB             - 2016 Aug 01
%   ncx2mixtcdf                - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT             - 2016 Aug 01
%   normBoxCox                 - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT             - 2017 Feb 20
%   normYJ                     - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT             - 2017 Mar 02
%   OPTbdp                     - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2017 May 17
%   OPTc                       - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT             - 2016 Aug 01
%   OPTeff                     - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT             - 2017 May 17
%   OPTpsi                     - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT             - 2017 May 17
%   OPTpsider                  - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT             - 2017 May 17
%   OPTpsix                    - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT             - 2017 May 17
%   OPTrho                     - Computes rho function for optimal weight function                                                                               - UTISTAT             - 2017 May 17
%   OPTwei                     - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT             - 2017 May 17
%   overlap                    - Computes the exact overlap given the parameters of the mixture                                                                  - CLUS-MixSim         - 2017 May 12
%   PoolClose                  - Closes the pool of MATLAB instances opened with PoolPrepare to execute code in parallel                                         - UTIGEN              - 2017 May 22
%   PoolPrepare                - Prepares a pool of MATLAB instances for executing code in parallel                                                              - UTIGEN              - 2017 May 22
%   position                   - Controls the position of the open figures                                                                                       - UTIGEN              - 2017 Feb 17
%   Powertra                   - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT             - 2017 Mar 02
%   publishFS                  - Enables to create automatic HELP FILES from structured .m function files                                                        - UTIHELP             - 2017 May 26
%   publishFunctionAlpha       - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 Aug 01
%   publishFunctionCate        - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 Aug 01
%   Qn                         - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT             - 2016 Aug 01
%   quickselectFS              - Finds the k-th order statistic                                                                                                  - UTIGEN              - 2016 Aug 01
%   RandIndexFS                - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT             - 2017 May 02
%   randsampleFS               - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                         - UTICOMB             - 2017 Apr 24
%   regressB                   - Computes Bayesian estimates of regression parameters                                                                            - REG-Bayes           - 2017 May 26
%   regressH                   - Fits a multiple linear regression model with heteroskedasticity                                                                 - REG-Hetero          - 2017 May 26
%   regressHart                - Fits a multiple linear regression model using ART heteroskedasticity                                                            - REG-Hetero          - 2017 May 26
%   regressHhar                - RegressH fits a multiple linear regression model with Harvey heteroskedasticity                                                 - REG-Hetero          - 2017 May 26
%   removeExtraSpacesLF        - Removes extra spaces and selected carriage returns from input string                                                            - UTIGEN              - 2016 Oct 03
%   resfwdplot                 - Plots the trajectories of the scaled (squared) residuals monitored                                                              - VIS-Reg             - 2017 Apr 20
%   resindexplot               - Plots the residuals from a regression analysis versus index number or any other variable                                        - VIS-Reg             - 2017 Feb 14
%   restreigen                 - Computes eigenvalues restriction (without Dykstra algorithm)                                                                    - CLUS-RobClaMULT     - 2017 May 12
%   rlga                       - Performs robust linear grouping analysis                                                                                        - CLUS-RobClaREG      - 2017 May 11
%   RobCov                     - Computes covariance matrix of robust regression coefficients                                                                    - REG-Regression      - 2017 May 26
%   RobRegrSize                - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size                    - REG-Regression      - 2017 Feb 17
%   rthin                      - Applies independent random thinning to a point pattern                                                                          - UTISTAT             - 2017 May 12
%   Score                      - Computes the score test for transformation                                                                                      - REG-Transformations - 2017 May 26
%   ScoreYJ                    - Score computes the score test for Yeo and Johnson transformation                                                                - REG-Transformations - 2017 May 26
%   SDest                      - Computes Stahel-Donoho robust estimator of dispersion-location                                                                  - MULT-Multivariate   - 2016 Aug 01
%   shuffling                  - Does a random permutation of the elements of input vector                                                                       - UTICOMB             - 2016 Dec 19
%   simdataset                 - Simulates and-or contaminates a dataset given the parameters of a finite mixture model with Gaussian components                 - CLUS-MixSim         - 2017 May 12
%   simdatasetreg              - Simulates a regression dataset given the parameters of a mixture regression model                                               - CLUS-MixSim         - 2017 May 02
%   Smult                      - Computes S estimators in multivariate analysis                                                                                  - MULT-Multivariate   - 2017 Apr 20
%   Smulteda                   - Computes S estimators in multivariate analysis for a series of values of bdp                                                    - MULT-Multivariate   - 2017 May 12
%   Sn                         - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT             - 2016 Aug 01
%   spmplot                    - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal                                     - VIS-Mult            - 2017 May 12
%   Sreg                       - Computes S estimators in linear regression                                                                                      - REG-Regression      - 2017 May 26
%   Sregeda                    - Computes S estimators in linear regression for a series of values of bdp                                                        - REG-Regression      - 2017 May 26
%   subsets                    - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                - UTICOMB             - 2017 May 08
%   suplabel                   - Places text as a title, xlabel, or ylabel on a group of subplots                                                                - UTIGEN              - 2017 May 04
%   tabulateFS                 - Create frequency table of unique values of x, excluding possible 0 counts                                                       - UTISTAT             - 2016 Aug 01
%   Taureg                     - Computes Tau estimators in linear regression                                                                                    - REG-Regression      - 2017 May 26
%   TBbdp                      - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT             - 2017 May 17
%   TBc                        - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT             - 2016 Aug 01
%   TBeff                      - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT             - 2017 May 17
%   TBpsi                      - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT             - 2017 May 17
%   TBpsider                   - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT             - 2017 May 10
%   TBpsix                     - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT             - 2017 May 17
%   TBrho                      - Computes rho function for Tukey's biweight                                                                                      - UTISTAT             - 2017 May 17
%   TBwei                      - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT             - 2017 May 17
%   tclust                     - Computes trimmed clustering with restricitons on the eigenvalues                                                                - CLUS-RobClaMULT     - 2017 May 12
%   tclustIC                   - Computes tclust for different number of groups k and restriction factors c                                                      - CLUS-RobClaMULT     - 2017 May 12
%   tclustICplot               - Plots information criterion as a function of c and k                                                                            - VIS-Clu             - 2017 May 12
%   tclustICsol                - Extracts a set of best relevant solutions                                                                                       - CLUS-RobClaMULT     - 2017 May 08
%   tclustreg                  - Performs robust linear grouping analysis                                                                                        - CLUS-RobClaREG      - 2017 May 26
%   tkmeans                    - Computes trimmed k-means                                                                                                        - CLUS-RobClaMULT     - 2017 May 12
%   triu2vec                   - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                       - UTIGEN              - 2016 Aug 01
%   unibiv                     - Has the purpose of detecting univariate and bivariate outliers                                                                  - MULT-Multivariate   - 2017 May 12
%   UnitsSameCluster           - Enables to control the labels of the clusters which contain predefined units                                                    - UTISTAT             - 2017 May 22
%   upperfracpos               - Positions two figures on the upper part of the screen                                                                           - UTIGEN              - 2016 Aug 01
%   winsor                     - Returns a winsorized copy of input                                                                                              - UTISTAT             - 2016 Aug 01
%   WNChygepdf                 - Returns Wallenius' non-central hypergeometric probability density values                                                        - UTISTAT             - 2017 May 12
%   wraptextFS                 - Formats long strings into wrapped text of specified width                                                                       - UTIGEN              - 2017 Mar 27
%   wthin                      - Thin a uni/bi-dimensional dataset                                                                                               - UTISTAT             - 2017 May 15
%   xmlcreateFS                - Create an XML file passing through publishFS                                                                                    - UTIHELP             - 2017 May 03
%   yXplot                     - Produces an interactive scatterplot of y against each variable of X in the input dataset                                        - VIS-Reg             - 2017 May 12
%   zscoreFS                   - Computes (robust) standardized z scores                                                                                         - UTIGEN              - 2017 May 04
