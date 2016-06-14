% FSDA
%
% File names, description, category and date last modified
%
%   Name                       - Description                                                                                                                     - Category            - Date last modified
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   add2spm                    - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                   - VIS-Mult            - 2016 Jun 01
%   addt                       - Produces the t test for an additional explanatory variable                                                                      - REG-Regression      - 2016 Jun 10
%   basicPower                 - Computes the basic power transformation                                                                                         - UTISTAT             - 2016 Jun 01
%   bc                         - Returns the Binomial coefficient                                                                                                - UTICOMB             - 2016 Jun 01
%   boxplotb                   - Computes a bivariate boxplot                                                                                                    - VIS-Mult            - 2016 Jun 01
%   brushFAN                   - Displays a GUI which enables brushing in the fanplot                                                                            - GUI                 - 2016 Jun 12
%   brushRES                   - Displays a GUI which enables brushing in resfwdplot                                                                             - GUI                 - 2016 Jun 12
%   brushROB                   - Displays a GUI which enables brushing in resindexplot                                                                           - GUI                 - 2016 Jun 12
%   cabc                       - Closes all open figures except the one in foreground (the current)                                                              - UTIGEN              - 2016 Jun 01
%   cascade                    - Is a third party function used in FSDA demos and examples                                                                       - UTIGEN              - 2016 Jun 01
%   cdsplot                    - Produces the candlestick plot for robust model selection in linear regression                                                   - VIS-Reg             - 2016 Jun 01
%   clickableMultiLegend       - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                             - UTIGEN              - 2016 Jun 01
%   combsFS                    - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                      - UTICOMB             - 2016 Jun 01
%   covplot                    - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                         - VIS-Mult            - 2016 Jun 01
%   ellipse                    - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT             - 2016 Jun 13
%   fanplot                    - Plots the fan plot for transformation in linear regression                                                                      - VIS-Reg             - 2016 Jun 01
%   findDir                    - Finds recursively all directories in root                                                                                       - UTIGEN              - 2016 Jun 01
%   findFile                   - Finds recursively all files in root                                                                                             - UTIGEN              - 2016 Jun 01
%   FSM                        - Gives an automatic outlier detection procedure in mult. analysis                                                                - MULT-Multivariate   - 2016 Jun 01
%   FSMbbm                     - Gives the units belonging to subset at step(s) msel of the forward search                                                       - MULT-Multivariate   - 2016 Jun 01
%   FSMbonfbound               - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT             - 2016 Jun 01
%   FSMbsb                     - Gives the units belonging to subset at step(s) msel of the forward search                                                       - MULT-Multivariate   - 2016 Jun 01
%   FSMeda                     - Performs forward search in multivariate analysis with exploratory data analysis purposes                                        - MULT-Multivariate   - 2016 Jun 01
%   FSMenvmmd                  - Computes the theoretical envelopes of Minimum MD outside subset during the search                                               - MULT-Multivariate   - 2016 Jun 01
%   FSMfan                     - Computes confirmatory lrt of a suggested transformation                                                                         - MULT-Transformations- 2016 Jun 01
%   FSMinvmmd                  - Converts values of minimum Mahalanobis distance into confidence levels                                                          - MULT-Multivariate   - 2016 Jun 01
%   FSMmmd                     - Monitors minMD                                                                                                                  - MULT-Multivariate   - 2016 Jun 01
%   FSMmmdeasy                 - Is exactly equal to minMD but much less efficient                                                                               - MULT-Multivariate   - 2016 Jun 01
%   FSMmmdrs                   - Performs random start monitoring of minimum Mahalanobis distance                                                                - CLUS-RobClaMULT     - 2016 Jun 01
%   FSMtra                     - Computes MLE of transformation parameters                                                                                       - MULT-Transformations- 2016 Jun 01
%   FSR                        - Gives an automatic outlier detection procedure in linear regression                                                             - REG-Regression      - 2016 Jun 01
%   FSRaddt                    - Produces t deletion tests for each explanatory variable                                                                         - REG-ModelSelection  - 2016 Jun 10
%   FSRB                       - Gives an automatic outlier detection procedure in Bayesian linear regression                                                    - REG-Bayes           - 2016 Jun 01
%   FSRBbsb                    - Returns the units belonging to the subset in each step of the Bayesian forward search                                           - REG-Regression      - 2016 Jun 01
%   FSRBeda                    - Enables to monitor several quantities in each step of the Bayesian search                                                       - REG-Bayes           - 2016 Jun 13
%   FSRBmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search             - REG-Bayes           - 2016 Jun 01
%   FSRbonfbound               - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT             - 2016 Jun 13
%   FSRBr                      - Bayesian forward search in linear regression reweighted                                                                         - REG-Bayes           - 2016 Jun 01
%   FSRbsb                     - Returns the units belonging to the subset in each step of the forward search                                                    - REG-Regression      - 2016 Jun 01
%   FSRcp                      - Monitors Cp and AIC for all models of interest of size smallp                                                                   - REG-ModelSelection  - 2016 Jun 01
%   FSReda                     - Enables to monitor several quantities in each step of the forward search                                                        - REG-Regression      - 2016 Jun 01
%   FSRenvmdr                  - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                - REG-Regression      - 2016 Jun 01
%   FSRfan                     - Monitors the values of the score test statistic for each lambda                                                                 - REG-Transformations - 2016 Jun 01
%   FSRH                       - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                             - REG-Hetero          - 2016 Jun 01
%   FSRHbsb                    - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                    - REG-Hetero          - 2016 Jun 01
%   FSRHeda                    - Enables to monitor several quantities in each step of the forward search                                                        - REG-Hetero          - 2016 Jun 01
%   FSRHmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search      - REG-Hetero          - 2016 Jun 01
%   FSRinvmdr                  - Converts values of minimum deletion residual into confidence levels                                                             - REG-Regression      - 2016 Jun 01
%   FSRmdr                     - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                      - REG-Regression      - 2016 Jun 01
%   FSRmdrrs                   - Performs random start monitoring of minimum deletion residual                                                                   - CLUS-RobClaREG      - 2016 Jun 01
%   FSRms                      - Performs robust model selection using flexible trimming in linear regression                                                    - REG-ModelSelection  - 2016 Jun 01
%   FSRr                       - Forward search in linear regression reweighted                                                                                  - REG-Regression      - 2016 Jun 01
%   HAbdp                      - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2016 Jun 01
%   HAc                        - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT             - 2016 Jun 01
%   HAeff                      - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT             - 2016 Jun 10
%   HApsi                      - Computes psi function  using Hampel proposal                                                                                    - UTISTAT             - 2016 Jun 01
%   HApsider                   - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT             - 2016 Jun 01
%   HApsix                     - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT             - 2016 Jun 01
%   HArho                      - Computes rho function  using Hampel proposal                                                                                    - UTISTAT             - 2016 Jun 01
%   HAwei                      - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT             - 2016 Jun 01
%   histFS                     - Plots a histogram with the elements in each bin grouped according to a vector of labels                                         - VIS-Reg             - 2016 Jun 01
%   HUeff                      - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT             - 2016 Jun 01
%   HUpsi                      - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT             - 2016 Jun 01
%   HUpsider                   - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT             - 2016 Jun 01
%   HUpsix                     - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT             - 2016 Jun 01
%   HUrho                      - Computes rho function for Huber                                                                                                 - UTISTAT             - 2016 Jun 10
%   HUwei                      - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT             - 2016 Jun 01
%   HYPbdp                     - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT             - 2016 Jun 10
%   HYPc                       - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT             - 2016 Jun 01
%   HYPck                      - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT             - 2016 Jun 01
%   HYPeff                     - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT             - 2016 Jun 01
%   HYPk                       - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT             - 2016 Jun 01
%   HYPpsi                     - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT             - 2016 Jun 01
%   HYPpsider                  - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT             - 2016 Jun 01
%   HYPpsix                    - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT             - 2016 Jun 01
%   HYPrho                     - Computes rho function  using hyperboloc tangent estimator                                                                       - UTISTAT             - 2016 Jun 01
%   HYPwei                     - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT             - 2016 Jun 01
%   inversegamcdf              - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT             - 2016 Jun 10
%   inversegaminv              - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT             - 2016 Jun 01
%   inversegampdf              - Computes inverse-gamma probability density function                                                                             - UTISTAT             - 2016 Jun 10
%   isfunction                 - Checks if a function exists                                                                                                     - UTIGEN              - 2016 Jun 01
%   levfwdplot                 - Plots the trajectories of leverage along the search                                                                             - VIS-Reg             - 2016 Jun 01
%   lexunrank                  - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                    - UTICOMB             - 2016 Jun 01
%   lga                        - Performs linear grouping analysis                                                                                               - CLUS-RobClaREG      - 2016 Jun 13
%   logmvnpdfFS                - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT             - 2016 Jun 01
%   LXS                        - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                            - REG-Regression      - 2016 Jun 01
%   mahalFS                    - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT             - 2016 Jun 01
%   makecontentsfileFS         - Extends Matlab function makecontentsfile                                                                                        - UTIHELP             - 2016 Jun 01
%   malfwdplot                 - Plots the trajectories of scaled Mahalanobis distances along the search                                                         - VIS-Mult            - 2016 Jun 10
%   malindexplot               - Plots the Mahalanobis distances versus a selected variable                                                                      - VIS-Mult            - 2016 Jun 01
%   mcd                        - Computes Minimum Covariance Determinant                                                                                         - MULT-Multivariate   - 2016 Jun 01
%   mdrplot                    - Plots the trajectory of minimum deletion residual (mdr)                                                                         - VIS-Reg             - 2016 Jun 01
%   MixSim                     - Generates k clusters in v dimensions with given overlap                                                                         - CLUS-MixSim         - 2016 Jun 01
%   MixSimreg                  - Generates k regression hyperplanes in p dimensions with given overlap                                                           - CLUS-MixSim         - 2016 Jun 01
%   mmdplot                    - Plots the trajectory of minimum Mahalanobis distance (mmd)                                                                      - VIS-Mult            - 2016 Jun 13
%   MMmult                     - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                          - MULT-Multivariate   - 2016 Jun 01
%   MMmultcore                 - Computes multivariate MM estimators for a selected fixed scale                                                                  - MULT-Multivariate   - 2016 Jun 01
%   MMreg                      - Computes MM estimator of regression coefficients                                                                                - REG-Regression      - 2016 Jun 01
%   MMregcore                  - Computes MM regression estimators for a selected fixed scale                                                                    - REG-Regression      - 2016 Jun 01
%   Mscale                     - Finds the M estimator of the scale                                                                                              - UTISTAT             - 2016 Jun 13
%   mve                        - Computes Minimum volume ellipsoid                                                                                               - MULT-Multivariate   - 2016 Jun 01
%   nchoosekFS                 - Returns the Binomial coefficient or matrix containing all combinations                                                          - UTICOMB             - 2016 Jun 01
%   ncx2mixtcdf                - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT             - 2016 Jun 01
%   normBoxCox                 - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT             - 2016 Jun 01
%   normYJ                     - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT             - 2016 Jun 01
%   OPTbdp                     - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2016 Jun 01
%   OPTc                       - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT             - 2016 Jun 01
%   OPTeff                     - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT             - 2016 Jun 01
%   OPTpsi                     - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT             - 2016 Jun 01
%   OPTpsider                  - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT             - 2016 Jun 01
%   OPTpsix                    - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT             - 2016 Jun 01
%   OPTrho                     - Computes rho function for optimal weight function                                                                               - UTISTAT             - 2016 Jun 01
%   OPTwei                     - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT             - 2016 Jun 01
%   overlap                    - Computes the exact overlap given the parameters of the mixture                                                                  - CLUS-MixSim         - 2016 Jun 01
%   PoolClose                  - Closes the pool of MATLAB instances opened with PoolPrepare to execute code in parallel                                         - UTIGEN              - 2016 Jun 01
%   PoolPrepare                - Prepares a pool of MATLAB instances for executing code in parallel                                                              - UTIGEN              - 2016 Jun 01
%   position                   - Controls the position of the open figures                                                                                       - UTIGEN              - 2016 Jun 01
%   Powertra                   - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT             - 2016 Jun 01
%   publishFS                  - Enables to create automatic HELP FILES from structured .m function files                                                        - UTIHELP             - 2016 Jun 13
%   publishFunctionAlpha       - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 Jun 01
%   publishFunctionCate        - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 Jun 01
%   Qn                         - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT             - 2016 Jun 01
%   quickselectFS              - Finds the k-th order statistic                                                                                                  - UTIGEN              - 2016 Jun 01
%   RandIndexFS                - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT             - 2016 Jun 01
%   randsampleFS               - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                         - UTICOMB             - 2016 Jun 01
%   regressB                   - Computes Bayesian estimates of regression parameters                                                                            - REG-Bayes           - 2016 Jun 01
%   regressH                   - Fits a multiple linear regression model with heteroskedasticity                                                                 - REG-Hetero          - 2016 Jun 01
%   regressHart                - Fits a multiple linear regression model using ART heteroskedasticity                                                            - REG-Hetero          - 2016 Jun 10
%   regressHhar                - RegressH fits a multiple linear regression model with Harvey heteroskedasticity                                                 - REG-Hetero          - 2016 Jun 01
%   resfwdplot                 - Plots the trajectories of the scaled (squared) residuals monitored                                                              - VIS-Reg             - 2016 Jun 10
%   resindexplot               - Plots the residuals from a regression analysis versus index number or any other variable                                        - VIS-Reg             - 2016 Jun 01
%   restreigen                 - Computes eigenvalues restriction (without Dykstra algorithm)                                                                    - CLUS-RobClaMULT     - 2016 Jun 01
%   rlga                       - Performs robust linear grouping analysis                                                                                        - CLUS-RobClaREG      - 2016 Jun 01
%   RobCov                     - Computes covariance matrix of robust regression coefficients                                                                    - REG-Regression      - 2016 Jun 01
%   RobRegrSize                - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size                    - REG-Regression      - 2016 Jun 10
%   Score                      - Computes the score test for transformation                                                                                      - REG-Transformations - 2016 Jun 01
%   SDest                      - Computes Stahel-Donoho robust estimator of dispersion-location                                                                  - MULT-Multivariate   - 2016 Jun 08
%   shuffling                  - Does a random permutation of the elements of input vector                                                                       - UTICOMB             - 2016 Jun 01
%   simdataset                 - Simulates and-or contaminates a dataset given the parameters of a finite mixture model with Gaussian components                 - CLUS-MixSim         - 2016 Jun 10
%   simdatasetreg              - Simulates a regression dataset given the parameters of a mixture regression model                                               - CLUS-MixSim         - 2016 Jun 10
%   Smult                      - Computes S estimators in multivariate analysis                                                                                  - MULT-Multivariate   - 2016 Jun 01
%   Sn                         - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT             - 2016 Jun 01
%   spmplot                    - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal                                     - VIS-Mult            - 2016 Jun 01
%   Sreg                       - Computes S estimators in linear regression                                                                                      - REG-Regression      - 2016 Jun 01
%   subsets                    - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                - UTICOMB             - 2016 Jun 01
%   tabulateFS                 - Create frequency table of unique values of x, excluding possible 0 counts                                                       - UTISTAT             - 2016 Jun 01
%   Taureg                     - Computes Tau estimators in linear regression                                                                                    - REG-Regression      - 2016 Jun 01
%   TBbdp                      - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT             - 2016 Jun 01
%   TBc                        - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT             - 2016 Jun 01
%   TBeff                      - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT             - 2016 Jun 10
%   TBpsi                      - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT             - 2016 Jun 01
%   TBpsider                   - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT             - 2016 Jun 01
%   TBpsix                     - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT             - 2016 Jun 01
%   TBrho                      - Computes rho function for Tukey's biweight                                                                                      - UTISTAT             - 2016 Jun 10
%   TBwei                      - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT             - 2016 Jun 10
%   tclust                     - Computes trimmed clustering with restricitons on the eigenvalues                                                                - CLUS-RobClaMULT     - 2016 Jun 06
%   tclustIC                   - Computes tclust for different number of groups k and restriction factors c                                                      - CLUS-RobClaMULT     - 2016 Jun 01
%   tclustICplot               - Plots information criterion as a function of c and k                                                                            - VIS-Clu             - 2016 Jun 10
%   tclustICsol                - Extracts a set of best relevant solutions                                                                                       - CLUS-RobClaMULT     - 2016 Jun 01
%   tclustreg                  - Performs robust clustering in regression                                                                                        - CLUS-RobClaREG      - 2016 Jun 10
%   tkmeans                    - Computes trimmed k-means                                                                                                        - CLUS-RobClaMULT     - 2016 Jun 01
%   triu2vec                   - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                       - UTIGEN              - 2016 Jun 01
%   unibiv                     - Has the purpose of detecting univariate and bivariate outliers                                                                  - MULT-Multivariate   - 2016 Jun 01
%   upperfracpos               - Positions two figures on the upper part of the screen                                                                           - UTIGEN              - 2016 Jun 12
%   winsor                     - Returns a winsorized copy of input                                                                                              - UTISTAT             - 2016 Jun 01
%   yXplot                     - Produces an interactive scatterplot of y against each variable of X in the input dataset                                        - VIS-Reg             - 2016 Jun 12
%   zscoreFS                   - Computes (robust) standardized z scores                                                                                         - UTIGEN              - 2016 Jun 01
