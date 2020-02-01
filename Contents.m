% FSDA
%
% File names, description, category and date last modified
%
%   Name                         - Description                                                                                                                        - Category            - Date last modified
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   ace                          - Computes alternative conditional expectation                                                                                       - REG-Transformations - 2019 May 20
%   aceplot                      - Produces the aceplot to visualize the results of ace                                                                               - REG-Regression      - 2019 May 20
%   add2spm                      - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                      - VIS-Mult            - 2019 Oct 17
%   add2yX                       - Adds objects (personalized clickable multilegends and text labels) to the yXplot                                                   - VIS-Reg             - 2019 Oct 17
%   addt                         - Produces the t test for an additional explanatory variable                                                                         - REG-Regression      - 2019 Oct 17
%   avas                         - Computes additivity and variance stabilization for regression                                                                      - REG-Transformations - 2019 May 20
%   barnardtest                  - Barnard's unconditional test                                                                                                       - MULT-Categorical    - 2019 May 20
%   basicPower                   - Computes the basic power transformation                                                                                            - UTISTAT             - 2019 Oct 17
%   bc                           - Returns the Binomial coefficient                                                                                                   - UTICOMB             - 2019 Oct 17
%   boxplotb                     - Computes a bivariate boxplot                                                                                                       - VIS-Mult            - 2019 Oct 17
%   boxtest                      - Performs Box test of equality of covariance matrices                                                                               - MULT-Multivariate   - 2019 Oct 17
%   brushFAN                     - Displays a GUI which enables brushing in the fanplot                                                                               - GUI                 - 2019 Oct 17
%   brushRES                     - Displays a GUI which enables brushing in resfwdplot                                                                                - GUI                 - 2019 Oct 17
%   brushROB                     - Displays a GUI which enables brushing in resindexplot                                                                              - GUI                 - 2019 Oct 17
%   bwe                          - Estimates the bandwidth smoothing parameter for kernel density estimation                                                          - UTISTAT             - 2019 Oct 17
%   cabc                         - Closes all open figures except the one in foreground (the current)                                                                 - UTIGEN              - 2019 Oct 17
%   carbikeplot                  - Produces the carbike plot to find best relevant clustering solutions                                                               - VIS-Clu             - 2019 Oct 17
%   cascade                      - Is a third party function used in FSDA demos and examples                                                                          - UTIGEN              - 2019 Oct 17
%   cdsplot                      - Produces the candlestick plot for robust model selection in linear regression                                                      - VIS-Reg             - 2019 Oct 17
%   clickableMultiLegend         - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                                - UTIGEN              - 2019 Oct 17
%   ClusterRelabel               - Enables to control the labels of the clusters which contain predefined units                                                       - UTISTAT             - 2019 Oct 17
%   combsFS                      - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                         - UTICOMB             - 2019 Oct 17
%   CorAna                       - Performs correspondence analysis                                                                                                   - MULT-Categorical    - 2019 Oct 17
%   CorAnaplot                   - Draws the Correspondence Analysis (CA) graphs with confidence ellipses                                                             - VIS-Mult            - 2019 May 20
%   corrNominal                  - Measures strength of association between two unordered (nominal) categorical variables                                             - MULT-Categorical    - 2019 Oct 17
%   corrOrdinal                  - Measures strength of association between two ordered categorical variables                                                         - MULT-Categorical    - 2019 Oct 17
%   covplot                      - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                            - VIS-Mult            - 2019 Oct 17
%   CressieRead                  - Computes the power divergence family                                                                                               - MULT-Categorical    - 2019 Oct 17
%   crosstab2datamatrix          - Recreates the original data matrix X from contingency table N                                                                      - MULT-Categorical    - 2019 Oct 17
%   ctsub                        - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                         - UTISTAT             - 2019 May 20
%   dempk                        - Performs a merging of components found by tkmeans                                                                                  - CLUS-RobClaMULT     - 2019 Oct 17
%   ellipse                      - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                         - UTISTAT             - 2019 Oct 17
%   exactcdf                     - Finds exact p-values                                                                                                               - UTISTAT             - 2019 May 20
%   existFS                      - Check if file exists and puts answer in a cached persistent variable                                                               - UTIGEN              - 2019 Oct 21
%   fanplot                      - Plots the fan plot for transformation in linear regression                                                                         - VIS-Reg             - 2019 Oct 17
%   findDir                      - Finds recursively all directories in root                                                                                          - UTIGEN              - 2019 Oct 17
%   findFile                     - Finds recursively all files in root                                                                                                - UTIGEN              - 2019 Oct 17
%   forecastTS                   - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                                  - REG-Regression      - 2019 Aug 31
%   FowlkesMallowsIndex          - Computes the Fowlkes and Mallows index                                                                                             - UTISTAT             - 2019 Oct 17
%   FSM                          - Gives an automatic outlier detection procedure in multivariate analysis                                                            - MULT-Multivariate   - 2019 Oct 17
%   FSMbonfbound                 - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                         - UTISTAT             - 2019 Oct 17
%   FSMbsb                       - Gives the units belonging to subset at step(s) msel of the forward search                                                          - MULT-Multivariate   - 2019 Oct 17
%   FSMeda                       - Performs forward search in multivariate analysis with exploratory data analysis purposes                                           - MULT-Multivariate   - 2019 Oct 17
%   FSMedaeasy                   - Is exactly equal to FSMeda but it is much less efficient                                                                           - MULT-Multivariate   - 2019 May 20
%   FSMenvmmd                    - Computes the theoretical envelopes of Minimum MD outside subset during the search                                                  - MULT-Multivariate   - 2019 Oct 17
%   FSMfan                       - Computes confirmatory lrt of a suggested transformation                                                                            - MULT-Transformations- 2019 Oct 17
%   FSMinvmmd                    - Converts values of minimum Mahalanobis distance into confidence levels                                                             - MULT-Multivariate   - 2019 Oct 17
%   FSMmmd                       - Monitors minMD                                                                                                                     - MULT-Multivariate   - 2019 Oct 17
%   FSMmmdeasy                   - Is exactly equal to FSMmmd but it is much less efficient                                                                           - MULT-Multivariate   - 2019 Oct 17
%   FSMmmdrs                     - Performs random start monitoring of minimum Mahalanobis distance                                                                   - CLUS-RobClaMULT     - 2019 Oct 17
%   FSMtra                       - Computes MLE of transformation parameters                                                                                          - MULT-Transformations- 2019 Oct 17
%   FSR                          - Gives an automatic outlier detection procedure in linear regression                                                                - REG-Regression      - 2019 Oct 17
%   FSRaddt                      - Produces t deletion tests for each explanatory variable                                                                            - REG-ModelSelection  - 2019 Oct 17
%   FSRB                         - Gives an automatic outlier detection procedure in Bayesian linear regression                                                       - REG-Bayes           - 2019 Oct 28
%   FSRBbsb                      - Returns the units belonging to the subset in each step of the Bayesian forward search                                              - REG-Regression      - 2019 Oct 17
%   FSRBeda                      - Enables to monitor several quantities in each step of the Bayesian search                                                          - REG-Bayes           - 2019 Oct 17
%   FSRBmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search                - REG-Bayes           - 2019 Oct 17
%   FSRbonfbound                 - Computes Bonferroni bounds for each step of the search (in linear regression)                                                      - UTISTAT             - 2019 Oct 17
%   FSRBr                        - Bayesian forward search in linear regression reweighted                                                                            - REG-Bayes           - 2019 Oct 17
%   FSRbsb                       - Returns the units belonging to the subset in each step of the forward search                                                       - REG-Regression      - 2019 Oct 17
%   FSRcp                        - Monitors Cp and AIC for all models of interest of size smallp                                                                      - REG-ModelSelection  - 2019 Oct 17
%   FSReda                       - Enables to monitor several quantities in each step of the forward search                                                           - REG-Regression      - 2019 Oct 17
%   FSRenvmdr                    - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                   - REG-Regression      - 2019 Oct 17
%   FSRfan                       - Monitors the values of the score test statistic for each lambda                                                                    - REG-Transformations - 2019 Oct 17
%   FSRH                         - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                                - REG-Hetero          - 2019 Oct 17
%   FSRHbsb                      - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                       - REG-Hetero          - 2019 Oct 17
%   FSRHeda                      - Enables to monitor several quantities in each step of the forward search                                                           - REG-Hetero          - 2019 Oct 17
%   FSRHmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search         - REG-Hetero          - 2019 Oct 17
%   FSRinvmdr                    - Converts values of minimum deletion residual into confidence levels                                                                - REG-Regression      - 2019 Oct 17
%   FSRmdr                       - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                         - REG-Regression      - 2019 Oct 17
%   FSRmdrrs                     - Performs random start monitoring of minimum deletion residual                                                                      - CLUS-RobClaREG      - 2019 Oct 17
%   FSRms                        - Performs robust model selection using flexible trimming in linear regression                                                       - REG-ModelSelection  - 2019 Oct 17
%   FSRr                         - Forward search in linear regression reweighted                                                                                     - REG-Regression      - 2019 Oct 17
%   FSRts                        - Is an automatic adaptive procedure to detect outliers in time series                                                               - REG-Regression      - 2019 May 20
%   FSRtsbsb                     - Returns the units belonging to the subset in each step of the forward search                                                       - REG-Regression      - 2019 May 20
%   FSRtsmdr                     - Computes minimum deletion residual for time series models in each step of the search                                               - REG-Regression      - 2019 May 20
%   genr8                        - Returns a vector of pseudorandom number sequence                                                                                   - UTISTAT             - 2019 May 04
%   genSigmaGPCM                 - Generates covariance matrix for the 14 Gaussian Parsimonious Clustering Models                                                     - CLUS-RobClaMULT     - 2019 May 20
%   GowerIndex                   - Computes matrix of similarity indexes using Gower metric                                                                           - CLUS-RobClaMULT     - 2019 May 20
%   GYfilt                       - Computes the Gervini-Yohai univariate outlier identifier                                                                           - UTISTAT             - 2019 Oct 17
%   HAbdp                        - Finds the constant c associated to the supplied breakdown point                                                                    - UTISTAT             - 2019 Oct 17
%   HAc                          - Computes breakdown point and efficiency associated with constant c                                                                 - UTISTAT             - 2019 Oct 17
%   HAeff                        - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                       - UTISTAT             - 2019 Oct 17
%   HApsi                        - Computes psi function  using Hampel proposal                                                                                       - UTISTAT             - 2019 Oct 17
%   HApsider                     - Computes derivative of psi function  using Hampel proposal                                                                         - UTISTAT             - 2019 Oct 17
%   HApsix                       - Computes psi function  using Hampel proposal times x                                                                               - UTISTAT             - 2019 Oct 17
%   HArho                        - Computes rho function  using Hampel proposal                                                                                       - UTISTAT             - 2019 Oct 17
%   HAwei                        - Computes weight function psi(u)/u using Hampel proposal                                                                            - UTISTAT             - 2019 Oct 17
%   histFS                       - Plots a histogram with the elements in each bin grouped according to a vector of labels                                            - VIS-Reg             - 2019 Oct 17
%   htmlwriteFS                  - Enables to create automatic HELP FILES from a specific MATLAB structure created with function mreadFS.m                            - UTIHELP             - 2019 Oct 17
%   HUeff                        - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                  - UTISTAT             - 2019 Oct 17
%   HUpsi                        - Computes psi function (derivative of rho function) for Huber                                                                       - UTISTAT             - 2019 Oct 17
%   HUpsider                     - Computes derivative of psi function (second derivative of rho function) for Huber                                                  - UTISTAT             - 2019 Oct 17
%   HUpsix                       - Computes psi function (derivative of rho function) times x for Huber                                                               - UTISTAT             - 2019 Oct 17
%   HUrho                        - Computes rho function for Huber                                                                                                    - UTISTAT             - 2019 Oct 17
%   HUwei                        - Computes weight function psi(u)/u for Huber                                                                                        - UTISTAT             - 2019 Oct 17
%   HYPbdp                       - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                      - UTISTAT             - 2019 Oct 17
%   HYPc                         - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)    - UTISTAT             - 2019 Oct 17
%   HYPck                        - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                            - UTISTAT             - 2019 Oct 17
%   HYPeff                       - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                          - UTISTAT             - 2019 Oct 17
%   HYPk                         - Computes breakdown point and efficiency for hyp. tan. estimator                                                                    - UTISTAT             - 2019 Oct 17
%   HYPpsi                       - Computes psi function for hyperbolic tangent estimator                                                                             - UTISTAT             - 2019 Oct 17
%   HYPpsider                    - Computes derivative of psi function for hyperbolic tangent estimator                                                               - UTISTAT             - 2019 Oct 17
%   HYPpsix                      - Computes psi function for hyperbolic tangent estimator times x                                                                     - UTISTAT             - 2019 Oct 17
%   HYPrho                       - Computes rho function  using hyperbolic tangent estimator                                                                          - UTISTAT             - 2019 Oct 17
%   HYPwei                       - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                 - UTISTAT             - 2019 Oct 17
%   inversegamcdf                - Computes inverse-gamma cumulative distribution function                                                                            - UTISTAT             - 2019 Oct 17
%   inversegaminv                - Inversegampdf Inverse-gamma cumulative distribution function                                                                       - UTISTAT             - 2019 Oct 17
%   inversegampdf                - Computes inverse-gamma probability density function                                                                                - UTISTAT             - 2019 Oct 17
%   isfunction                   - Checks if a function exists                                                                                                        - UTIGEN              - 2019 Oct 17
%   kdebiv                       - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                     - UTISTAT             - 2019 Oct 17
%   levfwdplot                   - Plots the trajectories of leverage along the search                                                                                - VIS-Reg             - 2019 Oct 17
%   lexunrank                    - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                       - UTICOMB             - 2019 Oct 17
%   lga                          - Performs linear grouping analysis                                                                                                  - CLUS-RobClaREG      - 2019 Oct 17
%   logmvnpdfFS                  - Produces log of Multivariate normal probability density function (pdf)                                                             - UTISTAT             - 2019 Oct 17
%   LTSts                        - Extends LTS estimator to time series                                                                                               - REG-Regression      - 2019 Oct 17
%   LTStsVarSel                  - Does variable selection in the robust time series model LTSts                                                                      - REG-Regression      - 2019 Sep 01
%   LXS                          - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                               - REG-Regression      - 2019 Oct 17
%   mahalFS                      - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                         - UTISTAT             - 2019 Oct 17
%   makecontentsfileFS           - Extends Matlab function makecontentsfile                                                                                           - UTIHELP             - 2019 Oct 17
%   malfwdplot                   - Plots the trajectories of scaled Mahalanobis distances along the search                                                            - VIS-Mult            - 2019 Oct 17
%   malindexplot                 - Plots the Mahalanobis distances versus a selected variable                                                                         - VIS-Mult            - 2019 Oct 17
%   mcd                          - Computes Minimum Covariance Determinant                                                                                            - MULT-Multivariate   - 2019 Oct 17
%   mdrplot                      - Plots the trajectory of minimum deletion residual (mdr)                                                                            - VIS-Reg             - 2019 Oct 17
%   MixSim                       - Generates k clusters in v dimensions with given overlap                                                                            - CLUS-MixSim         - 2019 Oct 17
%   MixSimreg                    - Generates k regression hyperplanes in p dimensions with given overlap                                                              - CLUS-MixSim         - 2019 Oct 17
%   mmdplot                      - Plots the trajectory of minimum Mahalanobis distance (mmd)                                                                         - VIS-Mult            - 2019 Oct 17
%   mmdrsplot                    - Plots the trajectories of minimum Mahalanobis distances from different starting points                                             - VIS-Mult            - 2019 Oct 17
%   MMmult                       - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                             - MULT-Multivariate   - 2019 Oct 17
%   MMmultcore                   - Computes multivariate MM estimators for a selected fixed scale                                                                     - MULT-Multivariate   - 2019 Oct 17
%   MMmulteda                    - Computes MM estimators in multivariate analysis for a series of values of eff                                                      - MULT-Multivariate   - 2019 Oct 17
%   MMreg                        - Computes MM estimator of regression coefficients                                                                                   - REG-Regression      - 2019 Oct 17
%   MMregcore                    - Computes MM regression estimators for a selected fixed scale                                                                       - REG-Regression      - 2019 Oct 17
%   MMregeda                     - Computes MM estimator in linear regression for a series of values of efficiency                                                    - REG-Regression      - 2019 Oct 17
%   mreadFS                      - Enables to create a structure with InputArgs/OptArgs/OutArgs ... from .m function files                                            - UTIHELP             - 2019 Oct 17
%   Mscale                       - Finds the M estimator of the scale                                                                                                 - UTISTAT             - 2019 Oct 17
%   mtR                          - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                           - UTISTAT             - 2019 Oct 17
%   mve                          - Computes Minimum volume ellipsoid                                                                                                  - MULT-Multivariate   - 2019 Oct 17
%   mveeda                       - Monitors Minimum volume ellipsoid for a series of values of bdp                                                                    - MULT-Multivariate   - 2019 Oct 17
%   nchoosekFS                   - Returns the Binomial coefficient or matrix containing all combinations                                                             - UTICOMB             - 2019 Oct 17
%   ncpci                        - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                   - UTISTAT             - 2019 Oct 17
%   ncx2mixtcdf                  - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                        - UTISTAT             - 2019 Oct 17
%   normBoxCox                   - Computes (normalized) Box-Cox transformation                                                                                       - UTISTAT             - 2019 Oct 17
%   normYJ                       - Computes (normalized) Yeo-Johnson transformation                                                                                   - UTISTAT             - 2019 Oct 17
%   openMatlabFileFromHTML       - Enables to put in HTML an hypertextual link to a specific MATLAB file                                                              - UTIGEN              - 2019 Oct 17
%   OPTbdp                       - Finds the constant c associated to the supplied breakdown point                                                                    - UTISTAT             - 2019 Oct 17
%   OPTc                         - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                        - UTISTAT             - 2019 Oct 17
%   OPTeff                       - Finds the constant c which is associated to the requested efficiency                                                               - UTISTAT             - 2019 Oct 17
%   OPTpsi                       - Computes psi function (derivative of rho function) for optimal weight function                                                     - UTISTAT             - 2019 Oct 17
%   OPTpsider                    - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                - UTISTAT             - 2019 Oct 17
%   OPTpsix                      - Computes psi function (derivative of rho function) times x                                                                         - UTISTAT             - 2019 Oct 17
%   OPTrho                       - Computes rho function for optimal weight function                                                                                  - UTISTAT             - 2019 Oct 17
%   OPTwei                       - Computes weight function psi(u)/u for optimal weight function                                                                      - UTISTAT             - 2019 Oct 17
%   overlap                      - Computes the exact overlap given the parameters of the mixture                                                                     - CLUS-MixSim         - 2019 Oct 17
%   overlapmap                   - Produce an interactive overlap map                                                                                                 - CLUS-RobClaMULT     - 2019 Oct 17
%   PoolClose                    - Closes the pool of MATLAB instances opened with PoolPrepare to execute code in parallel                                            - UTIGEN              - 2019 Oct 17
%   PoolPrepare                  - Prepares a pool of MATLAB instances for executing code in parallel                                                                 - UTIGEN              - 2019 Oct 17
%   position                     - Controls the position of the open figures                                                                                          - UTIGEN              - 2019 Oct 17
%   Powertra                     - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                            - UTISTAT             - 2019 Oct 17
%   publishBibliography          - Enables to create web page which contains the references inside the input .m files                                                 - UTIHELP             - 2019 Aug 31
%   publishFS                    - Enables to create automatic HELP FILES from structured .m function files                                                           - UTIHELP             - 2019 Nov 03
%   publishFunctionAlpha         - Enables to create web page which contains the alphabetical list of functions                                                       - UTIHELP             - 2019 Oct 17
%   publishFunctionCate          - Enables to create web page which contains the categorical list of functions                                                        - UTIHELP             - 2019 Oct 17
%   Qn                           - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                     - UTISTAT             - 2019 Oct 17
%   qqplotFS                     - Qqplot of studentized residuals with envelopes                                                                                     - VIS-Reg             - 2019 May 20
%   quickselectFS                - Finds the k-th order statistic                                                                                                     - UTIGEN              - 2019 Oct 17
%   RandIndexFS                  - Calculates Rand type Indices to compare two partitions                                                                             - UTISTAT             - 2019 Oct 17
%   randsampleFS                 - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                            - UTICOMB             - 2019 Oct 17
%   rcontFS                      - Generates a random two-way table with given marginal totals                                                                        - MULT-Categorical    - 2019 Oct 17
%   regressB                     - Computes Bayesian estimates of regression parameters                                                                               - REG-Bayes           - 2019 Oct 17
%   regressH                     - Fits a multiple linear regression model with heteroskedasticity                                                                    - REG-Hetero          - 2019 Oct 17
%   regressHart                  - Fits a multiple linear regression model using ART heteroskedasticity                                                               - REG-Hetero          - 2019 Oct 17
%   regressHhar                  - Fits a multiple linear regression model with Harvey heteroskedasticity                                                             - REG-Hetero          - 2019 Oct 17
%   regressts                    - Computes estimates of regression parameters for a time series models                                                               - REG-Regression      - 2019 Oct 03
%   removeExtraSpacesLF          - Removes extra spaces and selected carriage returns from input string                                                               - UTIGEN              - 2019 Oct 17
%   repDupValWithMean            - Replaces values of y which have non unique elements in vector x with local means                                                   - UTIGEN              - 2019 Oct 16
%   resfwdplot                   - Plots the trajectories of the monitored scaled (squared) residuals                                                                 - VIS-Reg             - 2019 Oct 17
%   resindexplot                 - Plots the residuals from a regression analysis versus index number or any other variable                                           - VIS-Reg             - 2019 Oct 17
%   restrdeter                   - Computes determinant restriction                                                                                                   - CLUS-RobClaMULT     - 2019 May 20
%   restrdeterGPCM               - Applies determinat restrictions for the 14 GPCM                                                                                    - CLUS-RobClaMULT     - 2019 Oct 09
%   restreigen                   - Computes eigenvalues restriction (without Dykstra algorithm)                                                                       - CLUS-RobClaMULT     - 2019 Oct 17
%   restreigeneasy               - Restreigen computes eigenvalues restriction (without Dykstra algorithm)                                                            - CLUS-RobClaMULT     - 2019 May 20
%   restrshapeGPCM               - Produces the restricted shape matrix for the 14 GPCM                                                                               - CLUS-RobClaMULT     - 2019 Oct 09
%   RKbdp                        - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                    - UTISTAT             - 2019 May 20
%   RKeff                        - Finds the constants c and M which are associated to the requested efficiency and ARP                                               - UTISTAT             - 2019 May 20
%   RKpsi                        - Computes psi function for Rocke (translated Tukey's) biweight                                                                      - UTISTAT             - 2019 May 20
%   RKpsider                     - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                    - UTISTAT             - 2019 May 20
%   RKpsix                       - Computes psi function times x for Rocke (translated Tukey's) biweight                                                              - UTISTAT             - 2019 May 20
%   RKrho                        - Computes rho function for Rocke (translated Tukey's) biweight                                                                      - UTISTAT             - 2019 May 20
%   RKwei                        - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                          - UTISTAT             - 2019 May 20
%   rlga                         - Performs robust linear grouping analysis                                                                                           - CLUS-RobClaREG      - 2019 Oct 17
%   rlsmo                        - Computes a running-lines smoother with global cross-validation                                                                     - REG-Transformations - 2019 May 20
%   RobCov                       - Computes covariance matrix of robust regression coefficients                                                                       - REG-Regression      - 2019 Oct 17
%   RobRegrSize                  - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size                       - REG-Regression      - 2019 Oct 17
%   rthin                        - Applies independent random thinning to a point pattern                                                                             - UTISTAT             - 2019 Oct 17
%   Score                        - Computes the score test for transformation                                                                                         - REG-Transformations - 2019 Oct 17
%   ScoreYJ                      - Computes the score test for Yeo and Johnson transformation                                                                         - REG-Transformations - 2019 Oct 17
%   ScoreYJpn                    - Computes the score test for YJ transformation separately for pos and neg observations                                              - REG-Transformations - 2019 May 20
%   SDest                        - Computes Stahel-Donoho robust estimator of dispersion-location                                                                     - MULT-Multivariate   - 2019 Oct 17
%   shuffling                    - Does a random permutation of the elements of input vector                                                                          - UTICOMB             - 2019 Oct 17
%   simdataset                   - Simulates and-or contaminates a dataset given the parameters of a finite mixture model with Gaussian components                    - CLUS-MixSim         - 2019 Oct 17
%   simdatasetreg                - Simulates a regression dataset given the parameters of a mixture regression model                                                  - CLUS-MixSim         - 2019 Oct 17
%   simulateTS                   - Simulate a time series with trend, time varying seasonal, level shift and irregular component                                      - REG-Regression      - 2019 Aug 31
%   smothr                       - Produces smoothed values with constraints                                                                                          - REG-Transformations - 2019 May 20
%   Smult                        - Computes S estimators in multivariate analysis                                                                                     - MULT-Multivariate   - 2019 Oct 17
%   Smulteda                     - Computes S estimators in multivariate analysis for a series of values of bdp                                                       - MULT-Multivariate   - 2019 Oct 17
%   Sn                           - Robust estimator of scale (robust version of Gini's average difference)                                                            - UTISTAT             - 2019 Oct 17
%   SparseTableTest              - Computes independence test for large and sparse contingency tables                                                                 - MULT-Categorical    - 2019 Oct 17
%   spmplot                      - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal and possibly robust bivariate contours - VIS-Mult            - 2019 Oct 17
%   Sreg                         - Computes S estimators in linear regression                                                                                         - REG-Regression      - 2019 Oct 17
%   Sregeda                      - Computes S estimators in linear regression for a series of values of bdp                                                           - REG-Regression      - 2019 Oct 17
%   subsets                      - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                   - UTICOMB             - 2019 Oct 17
%   suplabel                     - Places text as a title, xlabel, or ylabel on a group of subplots                                                                   - UTIGEN              - 2019 Oct 17
%   supsmu                       - Smooths scatterplots using Friedman's supersmoother algorithm                                                                      - REG-Transformations - 2019 May 20
%   tabulateFS                   - Creates frequency table of unique values of x, excluding possible 0 counts                                                         - UTISTAT             - 2019 Oct 17
%   Taureg                       - Computes Tau estimators in linear regression                                                                                       - REG-Regression      - 2019 Oct 17
%   TBbdp                        - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                               - UTISTAT             - 2019 Oct 17
%   TBc                          - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                            - UTISTAT             - 2019 Oct 17
%   TBeff                        - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                          - UTISTAT             - 2019 Oct 17
%   tBothSides                   - Allows users to transform both sides of a (nonlinear) regression model                                                             - REG-Regression      - 2019 Oct 09
%   TBpsi                        - Computes psi function (derivative of rho function) for Tukey's biweight                                                            - UTISTAT             - 2019 Oct 17
%   TBpsider                     - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                       - UTISTAT             - 2019 Oct 17
%   TBpsix                       - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                    - UTISTAT             - 2019 Oct 17
%   TBrho                        - Computes rho function for Tukey's biweight                                                                                         - UTISTAT             - 2019 Oct 17
%   TBwei                        - Computes weight function psi(u)/u for Tukey's biweight                                                                             - UTISTAT             - 2019 Oct 17
%   tclust                       - Computes trimmed clustering with scatter restrictions                                                                              - CLUS-RobClaMULT     - 2019 Oct 17
%   tclusteda                    - Computes tclust for a series of values of the trimming factor                                                                      - CLUS-RobClaMULT     - 2019 May 20
%   tclustIC                     - Computes tclust for different number of groups k and restriction factors c                                                         - CLUS-RobClaMULT     - 2019 Oct 17
%   tclustICplot                 - Plots information criterion as a function of c and k                                                                               - VIS-Clu             - 2019 Oct 17
%   tclustICsol                  - Extracts a set of best relevant solutions                                                                                          - CLUS-RobClaMULT     - 2019 Oct 17
%   tclustreg                    - Performs robust linear grouping analysis                                                                                           - CLUS-RobClaREG      - 2019 Oct 17
%   tclustregIC                  - Computes tclustreg for different number of groups k and restriction factors c                                                      - CLUS-RobClaMULT     - 2019 May 20
%   tkmeans                      - Computes trimmed k-means                                                                                                           - CLUS-RobClaMULT     - 2019 Oct 17
%   triu2vec                     - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                          - UTIGEN              - 2019 Oct 17
%   unibiv                       - Has the purpose of detecting univariate and bivariate outliers                                                                     - MULT-Multivariate   - 2019 Oct 17
%   upperfracpos                 - Positions two figures on the upper part of the screen                                                                              - UTIGEN              - 2019 Oct 17
%   verLessThanFS                - Compares version of MATLAB to specified version number                                                                             - UTIGEN              - 2019 May 20
%   vervaatrnd                   - Simulates random variates from the Vervaat perpetuity distribution                                                                 - UTISTAT             - 2019 May 20
%   vervaatsim                   - Returns a Vervaat perpetuity                                                                                                       - UTISTAT             - 2019 May 20
%   vervaatxdf                   - Returns the pdf and cdf of a Vervaat perpetuity                                                                                    - UTISTAT             - 2019 May 20
%   wedgeplot                    - Generates the double wedge plot of a time series                                                                                   - VIS-Reg             - 2019 Oct 17
%   winsor                       - Returns a winsorized copy of input                                                                                                 - UTISTAT             - 2019 Oct 17
%   WNChygepdf                   - Returns Wallenius' non-central hypergeometric probability density values                                                           - UTISTAT             - 2019 Oct 17
%   wraptextFS                   - Formats long strings into wrapped text of specified width                                                                          - UTIGEN              - 2019 Oct 17
%   wthin                        - Thins a uni/bi-dimensional dataset                                                                                                 - UTISTAT             - 2019 Oct 17
%   xmlcreateFS                  - Creates an XML file passing through publishFS                                                                                      - UTIHELP             - 2019 Oct 17
%   yXplot                       - Produces an interactive scatterplot of y against each variable of X in the input dataset                                           - VIS-Reg             - 2019 Oct 17
%   zscoreFS                     - Computes (robust) standardized z scores                                                                                            - UTIGEN              - 2019 Oct 17
