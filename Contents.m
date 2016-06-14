% FSDA
%
% File names, description, category and date last modified
%
%   Name                       - Description                                                                                                                     - Category            - Date last modified
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   add2spm                    - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                   - VIS-Mult            - 2016 Jun 14
%   addt                       - Produces the t test for an additional explanatory variable                                                                      - REG-Regression      - 2016 Jun 14
%   basicPower                 - Computes the basic power transformation                                                                                         - UTISTAT             - 2016 Jun 14
%   bc                         - Returns the Binomial coefficient                                                                                                - UTICOMB             - 2016 Jun 14
%   boxplotb                   - Computes a bivariate boxplot                                                                                                    - VIS-Mult            - 2016 Jun 14
%   brushFAN                   - Displays a GUI which enables brushing in the fanplot                                                                            - GUI                 - 2016 Jun 14
%   brushRES                   - Displays a GUI which enables brushing in resfwdplot                                                                             - GUI                 - 2016 Jun 14
%   brushROB                   - Displays a GUI which enables brushing in resindexplot                                                                           - GUI                 - 2016 Jun 14
%   cabc                       - Closes all open figures except the one in foreground (the current)                                                              - UTIGEN              - 2016 Jun 14
%   cascade                    - Is a third party function used in FSDA demos and examples                                                                       - UTIGEN              - 2016 Jun 14
%   cdsplot                    - Produces the candlestick plot for robust model selection in linear regression                                                   - VIS-Reg             - 2016 Jun 14
%   clickableMultiLegend       - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                             - UTIGEN              - 2016 Jun 14
%   combsFS                    - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                      - UTICOMB             - 2016 Jun 14
%   covplot                    - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                         - VIS-Mult            - 2016 Jun 14
%   ellipse                    - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT             - 2016 Jun 14
%   fanplot                    - Plots the fan plot for transformation in linear regression                                                                      - VIS-Reg             - 2016 Jun 14
%   findDir                    - Finds recursively all directories in root                                                                                       - UTIGEN              - 2016 Jun 14
%   findFile                   - Finds recursively all files in root                                                                                             - UTIGEN              - 2016 Jun 14
%   FSM                        - Gives an automatic outlier detection procedure in mult. analysis                                                                - MULT-Multivariate   - 2016 Jun 14
%   FSMbbm                     - Gives the units belonging to subset at step(s) msel of the forward search                                                       - MULT-Multivariate   - 2016 Jun 14
%   FSMbonfbound               - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT             - 2016 Jun 14
%   FSMbsb                     - Gives the units belonging to subset at step(s) msel of the forward search                                                       - MULT-Multivariate   - 2016 Jun 14
%   FSMeda                     - Performs forward search in multivariate analysis with exploratory data analysis purposes                                        - MULT-Multivariate   - 2016 Jun 14
%   FSMenvmmd                  - Computes the theoretical envelopes of Minimum MD outside subset during the search                                               - MULT-Multivariate   - 2016 Jun 14
%   FSMfan                     - Computes confirmatory lrt of a suggested transformation                                                                         - MULT-Transformations- 2016 Jun 14
%   FSMinvmmd                  - Converts values of minimum Mahalanobis distance into confidence levels                                                          - MULT-Multivariate   - 2016 Jun 14
%   FSMmmd                     - Monitors minMD                                                                                                                  - MULT-Multivariate   - 2016 Jun 14
%   FSMmmdeasy                 - Is exactly equal to minMD but much less efficient                                                                               - MULT-Multivariate   - 2016 Jun 14
%   FSMmmdrs                   - Performs random start monitoring of minimum Mahalanobis distance                                                                - CLUS-RobClaMULT     - 2016 Jun 14
%   FSMtra                     - Computes MLE of transformation parameters                                                                                       - MULT-Transformations- 2016 Jun 14
%   FSR                        - Gives an automatic outlier detection procedure in linear regression                                                             - REG-Regression      - 2016 Jun 14
%   FSRaddt                    - Produces t deletion tests for each explanatory variable                                                                         - REG-ModelSelection  - 2016 Jun 14
%   FSRB                       - Gives an automatic outlier detection procedure in Bayesian linear regression                                                    - REG-Bayes           - 2016 Jun 14
%   FSRBbsb                    - Returns the units belonging to the subset in each step of the Bayesian forward search                                           - REG-Regression      - 2016 Jun 14
%   FSRBeda                    - Enables to monitor several quantities in each step of the Bayesian search                                                       - REG-Bayes           - 2016 Jun 14
%   FSRBmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search             - REG-Bayes           - 2016 Jun 14
%   FSRbonfbound               - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT             - 2016 Jun 14
%   FSRBr                      - Bayesian forward search in linear regression reweighted                                                                         - REG-Bayes           - 2016 Jun 14
%   FSRbsb                     - Returns the units belonging to the subset in each step of the forward search                                                    - REG-Regression      - 2016 Jun 14
%   FSRcp                      - Monitors Cp and AIC for all models of interest of size smallp                                                                   - REG-ModelSelection  - 2016 Jun 14
%   FSReda                     - Enables to monitor several quantities in each step of the forward search                                                        - REG-Regression      - 2016 Jun 14
%   FSRenvmdr                  - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                - REG-Regression      - 2016 Jun 14
%   FSRfan                     - Monitors the values of the score test statistic for each lambda                                                                 - REG-Transformations - 2016 Jun 14
%   FSRH                       - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                             - REG-Hetero          - 2016 Jun 14
%   FSRHbsb                    - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                    - REG-Hetero          - 2016 Jun 14
%   FSRHeda                    - Enables to monitor several quantities in each step of the forward search                                                        - REG-Hetero          - 2016 Jun 14
%   FSRHmdr                    - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search      - REG-Hetero          - 2016 Jun 14
%   FSRinvmdr                  - Converts values of minimum deletion residual into confidence levels                                                             - REG-Regression      - 2016 Jun 14
%   FSRmdr                     - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                      - REG-Regression      - 2016 Jun 14
%   FSRmdrrs                   - Performs random start monitoring of minimum deletion residual                                                                   - CLUS-RobClaREG      - 2016 Jun 14
%   FSRms                      - Performs robust model selection using flexible trimming in linear regression                                                    - REG-ModelSelection  - 2016 Jun 14
%   FSRr                       - Forward search in linear regression reweighted                                                                                  - REG-Regression      - 2016 Jun 14
%   HAbdp                      - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2016 Jun 14
%   HAc                        - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT             - 2016 Jun 14
%   HAeff                      - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT             - 2016 Jun 14
%   HApsi                      - Computes psi function  using Hampel proposal                                                                                    - UTISTAT             - 2016 Jun 14
%   HApsider                   - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT             - 2016 Jun 14
%   HApsix                     - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT             - 2016 Jun 14
%   HArho                      - Computes rho function  using Hampel proposal                                                                                    - UTISTAT             - 2016 Jun 14
%   HAwei                      - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT             - 2016 Jun 14
%   histFS                     - Plots a histogram with the elements in each bin grouped according to a vector of labels                                         - VIS-Reg             - 2016 Jun 14
%   HUeff                      - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT             - 2016 Jun 14
%   HUpsi                      - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT             - 2016 Jun 14
%   HUpsider                   - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT             - 2016 Jun 14
%   HUpsix                     - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT             - 2016 Jun 14
%   HUrho                      - Computes rho function for Huber                                                                                                 - UTISTAT             - 2016 Jun 14
%   HUwei                      - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT             - 2016 Jun 14
%   HYPbdp                     - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT             - 2016 Jun 14
%   HYPc                       - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT             - 2016 Jun 14
%   HYPck                      - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT             - 2016 Jun 14
%   HYPeff                     - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT             - 2016 Jun 14
%   HYPk                       - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT             - 2016 Jun 14
%   HYPpsi                     - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT             - 2016 Jun 14
%   HYPpsider                  - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT             - 2016 Jun 14
%   HYPpsix                    - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT             - 2016 Jun 14
%   HYPrho                     - Computes rho function  using hyperboloc tangent estimator                                                                       - UTISTAT             - 2016 Jun 14
%   HYPwei                     - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT             - 2016 Jun 14
%   inversegamcdf              - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT             - 2016 Jun 14
%   inversegaminv              - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT             - 2016 Jun 14
%   inversegampdf              - Computes inverse-gamma probability density function                                                                             - UTISTAT             - 2016 Jun 14
%   isfunction                 - Checks if a function exists                                                                                                     - UTIGEN              - 2016 Jun 14
%   levfwdplot                 - Plots the trajectories of leverage along the search                                                                             - VIS-Reg             - 2016 Jun 14
%   lexunrank                  - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                    - UTICOMB             - 2016 Jun 14
%   lga                        - Performs linear grouping analysis                                                                                               - CLUS-RobClaREG      - 2016 Jun 14
%   logmvnpdfFS                - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT             - 2016 Jun 14
%   LXS                        - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                            - REG-Regression      - 2016 Jun 14
%   mahalFS                    - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT             - 2016 Jun 14
%   makecontentsfileFS         - Extends Matlab function makecontentsfile                                                                                        - UTIHELP             - 2016 Jun 14
%   malfwdplot                 - Plots the trajectories of scaled Mahalanobis distances along the search                                                         - VIS-Mult            - 2016 Jun 14
%   malindexplot               - Plots the Mahalanobis distances versus a selected variable                                                                      - VIS-Mult            - 2016 Jun 14
%   mcd                        - Computes Minimum Covariance Determinant                                                                                         - MULT-Multivariate   - 2016 Jun 14
%   mdrplot                    - Plots the trajectory of minimum deletion residual (mdr)                                                                         - VIS-Reg             - 2016 Jun 14
%   MixSim                     - Generates k clusters in v dimensions with given overlap                                                                         - CLUS-MixSim         - 2016 Jun 14
%   MixSimreg                  - Generates k regression hyperplanes in p dimensions with given overlap                                                           - CLUS-MixSim         - 2016 Jun 14
%   mmdplot                    - Plots the trajectory of minimum Mahalanobis distance (mmd)                                                                      - VIS-Mult            - 2016 Jun 14
%   MMmult                     - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                          - MULT-Multivariate   - 2016 Jun 14
%   MMmultcore                 - Computes multivariate MM estimators for a selected fixed scale                                                                  - MULT-Multivariate   - 2016 Jun 14
%   MMreg                      - Computes MM estimator of regression coefficients                                                                                - REG-Regression      - 2016 Jun 14
%   MMregcore                  - Computes MM regression estimators for a selected fixed scale                                                                    - REG-Regression      - 2016 Jun 14
%   Mscale                     - Finds the M estimator of the scale                                                                                              - UTISTAT             - 2016 Jun 14
%   mve                        - Computes Minimum volume ellipsoid                                                                                               - MULT-Multivariate   - 2016 Jun 14
%   nchoosekFS                 - Returns the Binomial coefficient or matrix containing all combinations                                                          - UTICOMB             - 2016 Jun 14
%   ncx2mixtcdf                - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT             - 2016 Jun 14
%   normBoxCox                 - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT             - 2016 Jun 14
%   normYJ                     - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT             - 2016 Jun 14
%   OPTbdp                     - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT             - 2016 Jun 14
%   OPTc                       - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT             - 2016 Jun 14
%   OPTeff                     - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT             - 2016 Jun 14
%   OPTpsi                     - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT             - 2016 Jun 14
%   OPTpsider                  - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT             - 2016 Jun 14
%   OPTpsix                    - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT             - 2016 Jun 14
%   OPTrho                     - Computes rho function for optimal weight function                                                                               - UTISTAT             - 2016 Jun 14
%   OPTwei                     - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT             - 2016 Jun 14
%   overlap                    - Computes the exact overlap given the parameters of the mixture                                                                  - CLUS-MixSim         - 2016 Jun 14
%   PoolClose                  - Closes the pool of MATLAB instances opened with PoolPrepare to execute code in parallel                                         - UTIGEN              - 2016 Jun 14
%   PoolPrepare                - Prepares a pool of MATLAB instances for executing code in parallel                                                              - UTIGEN              - 2016 Jun 14
%   position                   - Controls the position of the open figures                                                                                       - UTIGEN              - 2016 Jun 14
%   Powertra                   - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT             - 2016 Jun 14
%   publishFS                  - Enables to create automatic HELP FILES from structured .m function files                                                        - UTIHELP             - 2016 Jun 14
%   publishFunctionAlpha       - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 Jun 14
%   publishFunctionCate        - Enables to create web page which contains the alphabetical list of functions                                                    - UTIHELP             - 2016 Jun 14
%   Qn                         - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT             - 2016 Jun 14
%   quickselectFS              - Finds the k-th order statistic                                                                                                  - UTIGEN              - 2016 Jun 14
%   RandIndexFS                - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT             - 2016 Jun 14
%   randsampleFS               - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                         - UTICOMB             - 2016 Jun 14
%   regressB                   - Computes Bayesian estimates of regression parameters                                                                            - REG-Bayes           - 2016 Jun 14
%   regressH                   - Fits a multiple linear regression model with heteroskedasticity                                                                 - REG-Hetero          - 2016 Jun 14
%   regressHart                - Fits a multiple linear regression model using ART heteroskedasticity                                                            - REG-Hetero          - 2016 Jun 14
%   regressHhar                - RegressH fits a multiple linear regression model with Harvey heteroskedasticity                                                 - REG-Hetero          - 2016 Jun 14
%   resfwdplot                 - Plots the trajectories of the scaled (squared) residuals monitored                                                              - VIS-Reg             - 2016 Jun 14
%   resindexplot               - Plots the residuals from a regression analysis versus index number or any other variable                                        - VIS-Reg             - 2016 Jun 14
%   restreigen                 - Computes eigenvalues restriction (without Dykstra algorithm)                                                                    - CLUS-RobClaMULT     - 2016 Jun 14
%   rlga                       - Performs robust linear grouping analysis                                                                                        - CLUS-RobClaREG      - 2016 Jun 14
%   RobCov                     - Computes covariance matrix of robust regression coefficients                                                                    - REG-Regression      - 2016 Jun 14
%   RobRegrSize                - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size                    - REG-Regression      - 2016 Jun 14
%   Score                      - Computes the score test for transformation                                                                                      - REG-Transformations - 2016 Jun 14
%   SDest                      - Computes Stahel-Donoho robust estimator of dispersion-location                                                                  - MULT-Multivariate   - 2016 Jun 14
%   shuffling                  - Does a random permutation of the elements of input vector                                                                       - UTICOMB             - 2016 Jun 14
%   simdataset                 - Simulates and-or contaminates a dataset given the parameters of a finite mixture model with Gaussian components                 - CLUS-MixSim         - 2016 Jun 14
%   simdatasetreg              - Simulates a regression dataset given the parameters of a mixture regression model                                               - CLUS-MixSim         - 2016 Jun 14
%   Smult                      - Computes S estimators in multivariate analysis                                                                                  - MULT-Multivariate   - 2016 Jun 14
%   Sn                         - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT             - 2016 Jun 14
%   spmplot                    - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal                                     - VIS-Mult            - 2016 Jun 14
%   Sreg                       - Computes S estimators in linear regression                                                                                      - REG-Regression      - 2016 Jun 14
%   subsets                    - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                - UTICOMB             - 2016 Jun 14
%   tabulateFS                 - Create frequency table of unique values of x, excluding possible 0 counts                                                       - UTISTAT             - 2016 Jun 14
%   Taureg                     - Computes Tau estimators in linear regression                                                                                    - REG-Regression      - 2016 Jun 14
%   TBbdp                      - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT             - 2016 Jun 14
%   TBc                        - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT             - 2016 Jun 14
%   TBeff                      - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT             - 2016 Jun 14
%   TBpsi                      - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT             - 2016 Jun 14
%   TBpsider                   - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT             - 2016 Jun 14
%   TBpsix                     - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT             - 2016 Jun 14
%   TBrho                      - Computes rho function for Tukey's biweight                                                                                      - UTISTAT             - 2016 Jun 14
%   TBwei                      - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT             - 2016 Jun 14
%   tclust                     - Computes trimmed clustering with restricitons on the eigenvalues                                                                - CLUS-RobClaMULT     - 2016 Jun 14
%   tclustIC                   - Computes tclust for different number of groups k and restriction factors c                                                      - CLUS-RobClaMULT     - 2016 Jun 14
%   tclustICplot               - Plots information criterion as a function of c and k                                                                            - VIS-Clu             - 2016 Jun 14
%   tclustICsol                - Extracts a set of best relevant solutions                                                                                       - CLUS-RobClaMULT     - 2016 Jun 14
%   tclustreg                  - Performs robust clustering in regression                                                                                        - CLUS-RobClaREG      - 2016 Jun 14
%   tkmeans                    - Computes trimmed k-means                                                                                                        - CLUS-RobClaMULT     - 2016 Jun 14
%   triu2vec                   - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                       - UTIGEN              - 2016 Jun 14
%   unibiv                     - Has the purpose of detecting univariate and bivariate outliers                                                                  - MULT-Multivariate   - 2016 Jun 14
%   upperfracpos               - Positions two figures on the upper part of the screen                                                                           - UTIGEN              - 2016 Jun 14
%   winsor                     - Returns a winsorized copy of input                                                                                              - UTISTAT             - 2016 Jun 14
%   yXplot                     - Produces an interactive scatterplot of y against each variable of X in the input dataset                                        - VIS-Reg             - 2016 Jun 14
%   zscoreFS                   - Computes (robust) standardized z scores                                                                                         - UTIGEN              - 2016 Jun 14
