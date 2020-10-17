%% FSDA
%
% FSDA code for any function is open and extensible.
%
% FSDA  is developed by the Robust Statistics Academy of the University of Parma
% (http://rosa.unipr.it) jointly with the
% Joint Research Centre of European Commission
% (https://ec.europa.eu/jrc/en/about/jrc-site/ispra)
%
% The source code is also available on github
% (https://uniprjrc.github.io/FSDA/)
%
% The html documentation of each function in Mathworks style can be found
% in the supplementary software section of MATLAB help system.
% A copy of the documentation can also be found at the web address
% (http://rosa.unipr.it/FSDA/guide.html)
%
% File names, description, category and date last modified
%
%   Name                         - Description                                                                                                                         - Category            - Date last modified
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   ace                          - Computes alternative conditional expectation                                                                                        - REG-Transformations - 2020 Jan 29
%   aceplot                      - Produces the aceplot to visualize the results of ace                                                                                - VIS-Reg             - 2020 May 26
%   add2spm                      - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                       - VIS-Mult            - 2020 Jan 29
%   add2yX                       - Adds objects (personalized clickable multilegends and text labels) to the yXplot                                                    - VIS-Reg             - 2020 Apr 10
%   addt                         - Produces the t test for an additional explanatory variable                                                                          - REG-Regression      - 2020 Feb 20
%   avas                         - Computes additivity and variance stabilization for regression                                                                       - REG-Transformations - 2020 Jan 29
%   barnardtest                  - Barnard's unconditional test                                                                                                        - MULT-Categorical    - 2020 Jan 29
%   basicPower                   - Computes the basic power transformation                                                                                             - UTISTAT             - 2020 Feb 05
%   bc                           - Returns the Binomial coefficient                                                                                                    - UTICOMB             - 2020 Jan 29
%   boxcoxR                      - Finds MLE of lambda in linear regression (and confidence interval) using Box Cox, YJ or extended YJ  transformation                 - REG-Transformations - 2020 Jun 07
%   boxplotb                     - Computes a bivariate boxplot                                                                                                        - VIS-Mult            - 2020 Jan 29
%   boxtest                      - Performs Box test of equality of covariance matrices                                                                                - MULT-Multivariate   - 2020 Jan 29
%   brushFAN                     - Displays a GUI which enables brushing in the fanplot                                                                                - GUI                 - 2020 Jan 29
%   brushRES                     - Displays a GUI which enables brushing in resfwdplot                                                                                 - GUI                 - 2020 Jan 29
%   brushROB                     - Displays a GUI which enables brushing in resindexplot                                                                               - GUI                 - 2020 Jan 29
%   bwe                          - Estimates the bandwidth smoothing parameter for kernel density estimation                                                           - UTISTAT             - 2020 Feb 05
%   cabc                         - Closes all open figures except the one in foreground (the current)                                                                  - UTIGEN              - 2020 Jan 29
%   carbikeplot                  - Produces the carbike plot to find best relevant clustering solutions                                                                - VIS-Clu             - 2020 Oct 12
%   carbikeplotGPCM              - Carbikeplot produces the carbike plot to find best relevant clustering solutions                                                    - VIS-Clu             - 2020 Oct 12
%   cascade                      - Is a third party function used in FSDA demos and examples                                                                           - UTIGEN              - 2020 Jan 29
%   cdsplot                      - Produces the candlestick plot for robust model selection in linear regression                                                       - VIS-Reg             - 2020 Jan 29
%   clickableMultiLegend         - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                                 - UTIGEN              - 2020 Aug 27
%   ClusterRelabel               - Enables to control the labels of the clusters which contain predefined units                                                        - UTISTAT             - 2020 Apr 21
%   combsFS                      - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                          - UTICOMB             - 2020 Jan 29
%   CorAna                       - Performs correspondence analysis                                                                                                    - MULT-Categorical    - 2020 Jan 29
%   CorAnaplot                   - Draws the Correspondence Analysis (CA) graphs with confidence ellipses                                                              - VIS-Mult            - 2020 Feb 18
%   corrNominal                  - Measures strength of association between two unordered (nominal) categorical variables                                              - MULT-Categorical    - 2020 Oct 10
%   corrOrdinal                  - Measures strength of association between two ordered categorical variables                                                          - MULT-Categorical    - 2020 Jan 29
%   covplot                      - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                             - VIS-Mult            - 2020 Jan 29
%   CressieRead                  - Computes the power divergence family                                                                                                - MULT-Categorical    - 2020 Jan 29
%   crosstab2datamatrix          - Recreates the original data matrix X from contingency table N                                                                       - MULT-Categorical    - 2020 Feb 05
%   ctsub                        - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                          - UTISTAT             - 2020 Feb 18
%   dempk                        - Performs a merging of components found by tkmeans                                                                                   - CLUS-RobClaMULT     - 2020 Feb 20
%   ellipse                      - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                          - UTISTAT             - 2020 Feb 05
%   estregimeTAR                 - Estimate a regression model with OLS in one of the regimes of a TAR model                                                           - REG-Regression      - 2020 Aug 07
%   exactcdf                     - Finds exact p-values                                                                                                                - UTISTAT             - 2020 Jan 29
%   existFS                      - Check if file exists and puts answer in a cached persistent variable                                                                - UTIGEN              - 2020 Feb 05
%   fanBIC                       - Uses the output of FSRfan to choose the best value of the transformation parameter in linear regression                             - VIS-Reg             - 2020 Jun 22
%   fanBICpn                     - Uses the output of FSRfan called with input option family 'YJpn' to choose la_P and la_N                                            - VIS-Reg             - 2020 Jun 05
%   fanplot                      - Plots the fan plot for transformation in linear regression                                                                          - VIS-Reg             - 2020 Jun 04
%   findDir                      - Finds recursively all directories in root                                                                                           - UTIGEN              - 2020 Apr 21
%   findFile                     - Finds recursively all files in root                                                                                                 - UTIGEN              - 2020 Apr 21
%   forecastTS                   - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                                   - REG-Regression      - 2020 Jan 29
%   FowlkesMallowsIndex          - Computes the Fowlkes and Mallows index                                                                                              - UTISTAT             - 2020 Jan 29
%   FSM                          - Gives an automatic outlier detection procedure in multivariate analysis                                                             - MULT-Multivariate   - 2020 Jan 29
%   FSMbonfbound                 - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                          - UTISTAT             - 2020 Feb 05
%   FSMbsb                       - Gives the units belonging to subset at step(s) msel of the forward search                                                           - MULT-Multivariate   - 2020 Jan 29
%   FSMeda                       - Performs forward search in multivariate analysis with exploratory data analysis purposes                                            - MULT-Multivariate   - 2020 Jan 29
%   FSMedaeasy                   - Is exactly equal to FSMeda but it is much less efficient                                                                            - MULT-Multivariate   - 2020 Jan 29
%   FSMenvmmd                    - Computes the theoretical envelopes of Minimum MD outside subset during the search                                                   - MULT-Multivariate   - 2020 Jan 29
%   FSMfan                       - Computes confirmatory lrt of a suggested transformation                                                                             - MULT-Transformations- 2020 Jan 29
%   FSMinvmmd                    - Converts values of minimum Mahalanobis distance into confidence levels                                                              - MULT-Multivariate   - 2020 Jan 29
%   FSMmmd                       - Monitors minMD                                                                                                                      - MULT-Multivariate   - 2020 May 23
%   FSMmmdeasy                   - Is exactly equal to FSMmmd but it is much less efficient                                                                            - MULT-Multivariate   - 2020 Jan 29
%   FSMmmdrs                     - Performs random start monitoring of minimum Mahalanobis distance                                                                    - CLUS-RobClaMULT     - 2020 Apr 21
%   FSMtra                       - Computes MLE of transformation parameters                                                                                           - MULT-Transformations- 2020 Jan 29
%   FSR                          - Computes forward search estimator in linear regression                                                                              - REG-Regression      - 2020 Oct 13
%   FSRaddt                      - Produces t deletion tests for each explanatory variable                                                                             - REG-ModelSelection  - 2020 Apr 02
%   FSRB                         - Gives an automatic outlier detection procedure in Bayesian linear regression                                                        - REG-Bayes           - 2020 Jul 08
%   FSRBbsb                      - Returns the units belonging to the subset in each step of the Bayesian forward search                                               - REG-Regression      - 2020 Apr 02
%   FSRBeda                      - Enables to monitor several quantities in each step of the Bayesian search                                                           - REG-Bayes           - 2020 Apr 21
%   FSRBmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search                 - REG-Bayes           - 2020 Apr 02
%   FSRbonfbound                 - Computes Bonferroni bounds for each step of the search (in linear regression)                                                       - UTISTAT             - 2020 Feb 05
%   FSRBr                        - Bayesian forward search in linear regression reweighted                                                                             - REG-Bayes           - 2020 Feb 05
%   FSRbsb                       - Returns the units belonging to the subset in each step of the forward search                                                        - REG-Regression      - 2020 Apr 27
%   FSRcp                        - Monitors Cp and AIC for all models of interest of size smallp                                                                       - REG-ModelSelection  - 2020 Apr 02
%   FSReda                       - Enables to monitor several quantities in each step of the forward search                                                            - REG-Regression      - 2020 May 26
%   FSRenvmdr                    - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                    - REG-Regression      - 2020 Jan 29
%   FSRfan                       - Monitors the values of the score test statistic for each lambda                                                                     - REG-Transformations - 2020 Jun 11
%   FSRH                         - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                                 - REG-Hetero          - 2020 Feb 18
%   FSRHbsb                      - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                        - REG-Hetero          - 2020 Apr 02
%   FSRHeda                      - Enables to monitor several quantities in each step of the forward search                                                            - REG-Hetero          - 2020 Apr 02
%   FSRHmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search          - REG-Hetero          - 2020 Apr 02
%   FSRinvmdr                    - Converts values of minimum deletion residual into confidence levels                                                                 - REG-Regression      - 2020 Jan 29
%   FSRmdr                       - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                          - REG-Regression      - 2020 Jun 06
%   FSRmdrrs                     - Performs random start monitoring of minimum deletion residual                                                                       - CLUS-RobClaREG      - 2020 Jun 06
%   FSRms                        - Performs robust model selection using flexible trimming in linear regression                                                        - REG-ModelSelection  - 2020 Apr 21
%   FSRr                         - Forward search in linear regression reweighted                                                                                      - REG-Regression      - 2020 Feb 20
%   FSRts                        - Is an automatic adaptive procedure to detect outliers in time series                                                                - REG-Regression      - 2020 Feb 20
%   FSRtsbsb                     - Returns the units belonging to the subset in each step of the forward search                                                        - REG-Regression      - 2020 Apr 02
%   FSRtsmdr                     - Computes minimum deletion residual for time series models in each step of the search                                                - REG-Regression      - 2020 Apr 21
%   funnelchart                  - Displays a funnel chart                                                                                                             - VIS-Mult            - 2020 Oct 10
%   genSigmaGPCM                 - Generates covariance matrix for the 14 Gaussian Parsimonious Clustering Models                                                      - CLUS-RobClaMULT     - 2020 Jul 14
%   GowerIndex                   - Computes matrix of similarity indexes using Gower metric                                                                            - CLUS-RobClaMULT     - 2020 Jan 29
%   GYfilt                       - Computes the Gervini-Yohai univariate outlier identifier                                                                            - UTISTAT             - 2020 Jan 29
%   HAbdp                        - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT             - 2020 Jan 29
%   HAc                          - Computes breakdown point and efficiency associated with constant c                                                                  - UTISTAT             - 2020 Feb 05
%   HAeff                        - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                        - UTISTAT             - 2020 Jan 29
%   HApsi                        - Computes psi function  using Hampel proposal                                                                                        - UTISTAT             - 2020 Feb 05
%   HApsider                     - Computes derivative of psi function  using Hampel proposal                                                                          - UTISTAT             - 2020 Feb 05
%   HApsix                       - Computes psi function  using Hampel proposal times x                                                                                - UTISTAT             - 2020 Feb 05
%   HArho                        - Computes rho function  using Hampel proposal                                                                                        - UTISTAT             - 2020 Feb 05
%   HAwei                        - Computes weight function psi(u)/u using Hampel proposal                                                                             - UTISTAT             - 2020 Jan 29
%   histFS                       - Plots a histogram with the elements in each bin grouped according to a vector of labels                                             - VIS-Reg             - 2020 Feb 18
%   htmlwriteFS                  - Is an obsolete function which will be removed in future releases. Use publishFS.m instead                                           - UTIHELP             - 2020 Apr 21
%   HUeff                        - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                   - UTISTAT             - 2020 Jan 29
%   HUpsi                        - Computes psi function (derivative of rho function) for Huber                                                                        - UTISTAT             - 2020 Jan 29
%   HUpsider                     - Computes derivative of psi function (second derivative of rho function) for Huber                                                   - UTISTAT             - 2020 Feb 05
%   HUpsix                       - Computes psi function (derivative of rho function) times x for Huber                                                                - UTISTAT             - 2020 Feb 05
%   HUrho                        - Computes rho function for Huber                                                                                                     - UTISTAT             - 2020 Feb 13
%   HUwei                        - Computes weight function psi(u)/u for Huber                                                                                         - UTISTAT             - 2020 Feb 05
%   HYPbdp                       - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                       - UTISTAT             - 2020 Feb 20
%   HYPc                         - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)     - UTISTAT             - 2020 Feb 05
%   HYPck                        - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                             - UTISTAT             - 2020 Feb 05
%   HYPeff                       - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                           - UTISTAT             - 2020 Feb 05
%   HYPk                         - Computes breakdown point and efficiency for hyp. tan. estimator                                                                     - UTISTAT             - 2020 Feb 05
%   HYPpsi                       - Computes psi function for hyperbolic tangent estimator                                                                              - UTISTAT             - 2020 Jan 29
%   HYPpsider                    - Computes derivative of psi function for hyperbolic tangent estimator                                                                - UTISTAT             - 2020 Jan 29
%   HYPpsix                      - Computes psi function for hyperbolic tangent estimator times x                                                                      - UTISTAT             - 2020 Jan 29
%   HYPrho                       - Computes rho function  using hyperbolic tangent estimator                                                                           - UTISTAT             - 2020 Jan 29
%   HYPwei                       - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                  - UTISTAT             - 2020 Jan 29
%   inversegamcdf                - Computes inverse-gamma cumulative distribution function                                                                             - UTISTAT             - 2020 Jan 29
%   inversegaminv                - Inversegampdf Inverse-gamma cumulative distribution function                                                                        - UTISTAT             - 2020 Jan 29
%   inversegampdf                - Computes inverse-gamma probability density function                                                                                 - UTISTAT             - 2020 Jan 29
%   isfunction                   - Checks if a function exists                                                                                                         - UTIGEN              - 2020 Apr 21
%   kdebiv                       - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                      - UTISTAT             - 2020 Feb 18
%   levfwdplot                   - Plots the trajectories of leverage along the search                                                                                 - VIS-Reg             - 2020 Apr 21
%   lexunrank                    - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                        - UTICOMB             - 2020 May 20
%   lga                          - Performs linear grouping analysis                                                                                                   - CLUS-RobClaREG      - 2020 Apr 28
%   logmvnpdfFS                  - Produces log of Multivariate normal probability density function (pdf)                                                              - UTISTAT             - 2020 Jul 03
%   LTSts                        - Extends LTS estimator to time series                                                                                                - REG-Regression      - 2020 May 26
%   LTStsVarSel                  - Does variable selection in the robust time series model LTSts                                                                       - REG-Regression      - 2020 Apr 21
%   LXS                          - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                                - REG-Regression      - 2020 Apr 02
%   mahalFS                      - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT             - 2020 Feb 05
%   makecontentsfileFS           - Extends Matlab function makecontentsfile                                                                                            - UTIHELP             - 2020 Oct 17
%   malfwdplot                   - Plots the trajectories of scaled Mahalanobis distances along the search                                                             - VIS-Mult            - 2020 Sep 03
%   malindexplot                 - Plots the Mahalanobis distances versus a selected variable                                                                          - VIS-Mult            - 2020 Jan 29
%   mcd                          - Computes Minimum Covariance Determinant                                                                                             - MULT-Multivariate   - 2020 Feb 20
%   mdpd                         - Computes Minimum Distance Power Divergence statistics                                                                               - UTISTAT             - 2020 Apr 21
%   mdpdR                        - Allows to apply Minimum Density Power Divergence criterion to parametric regression problems                                        - REG-Regression      - 2020 Apr 03
%   mdpdReda                     - Allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems                                     - REG-Regression      - 2020 Feb 24
%   mdrplot                      - Plots the trajectory of minimum deletion residual (mdr)                                                                             - VIS-Reg             - 2020 Feb 18
%   mdrrsplot                    - Plots the trajectory of minimum deletion residual from random starts                                                                - CLUS-RobClaREG      - 2020 Jun 15
%   MixSim                       - Generates k clusters in v dimensions with given overlap                                                                             - CLUS-MixSim         - 2020 Sep 28
%   MixSimreg                    - Generates k regression hyperplanes in p dimensions with given overlap                                                               - CLUS-MixSim         - 2020 Apr 21
%   mmdplot                      - Plots the trajectory of minimum Mahalanobis distance (mmd)                                                                          - VIS-Mult            - 2020 Apr 24
%   mmdrsplot                    - Plots the trajectories of minimum Mahalanobis distances from different starting points                                              - VIS-Mult            - 2020 Oct 12
%   MMmult                       - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                              - MULT-Multivariate   - 2020 Feb 20
%   MMmultcore                   - Computes multivariate MM estimators for a selected fixed scale                                                                      - MULT-Multivariate   - 2020 Apr 21
%   MMmulteda                    - Computes MM estimators in multivariate analysis for a series of values of eff                                                       - MULT-Multivariate   - 2020 Apr 25
%   MMreg                        - Computes MM estimator of regression coefficients                                                                                    - REG-Regression      - 2020 Apr 25
%   MMregcore                    - Computes MM regression estimators for a selected fixed scale                                                                        - REG-Regression      - 2020 Apr 25
%   MMregeda                     - Computes MM estimator in linear regression for a series of values of efficiency                                                     - REG-Regression      - 2020 Feb 20
%   mreadFS                      - Enables to create a structure with InputArgs/OptArgs/OutArgs ... from .m function files (OBSOLETE FUNCTION REPLACED BY publishFS.m) - UTIHELP             - 2020 Apr 21
%   Mscale                       - Finds the M estimator of the scale                                                                                                  - UTISTAT             - 2020 Feb 20
%   mtR                          - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                            - UTISTAT             - 2020 Feb 18
%   mve                          - Computes Minimum volume ellipsoid                                                                                                   - MULT-Multivariate   - 2020 Feb 20
%   mveeda                       - Monitors Minimum volume ellipsoid for a series of values of bdp                                                                     - MULT-Multivariate   - 2020 Feb 20
%   nchoosekFS                   - Returns the Binomial coefficient or matrix containing all combinations                                                              - UTICOMB             - 2020 Jan 29
%   ncpci                        - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                    - UTISTAT             - 2020 Feb 05
%   ncx2mixtcdf                  - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                         - UTISTAT             - 2020 Sep 11
%   normBoxCox                   - Computes (normalized) Box-Cox transformation                                                                                        - UTISTAT             - 2020 Feb 05
%   normYJ                       - Computes (normalized) Yeo-Johnson transformation                                                                                    - UTISTAT             - 2020 May 12
%   normYJpn                     - NormYJ computes (normalized) extended Yeo-Johnson transformation                                                                    - UTISTAT             - 2020 May 12
%   openMatlabFileFromHTML       - Enables to put in HTML an hypertextual link to a specific MATLAB file                                                               - UTIGEN              - 2020 Apr 21
%   OPTbdp                       - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT             - 2020 Jan 29
%   OPTc                         - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                         - UTISTAT             - 2020 Feb 20
%   OPTeff                       - Finds the constant c which is associated to the requested efficiency                                                                - UTISTAT             - 2020 Apr 21
%   OPTpsi                       - Computes psi function (derivative of rho function) for optimal weight function                                                      - UTISTAT             - 2020 Jan 29
%   OPTpsider                    - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                 - UTISTAT             - 2020 Jan 29
%   OPTpsix                      - Computes psi function (derivative of rho function) times x                                                                          - UTISTAT             - 2020 Jan 29
%   OPTrho                       - Computes rho function for optimal weight function                                                                                   - UTISTAT             - 2020 Jan 29
%   OPTwei                       - Computes weight function psi(u)/u for optimal weight function                                                                       - UTISTAT             - 2020 Oct 11
%   overlap                      - Computes the exact overlap given the parameters of the mixture                                                                      - CLUS-MixSim         - 2020 Sep 27
%   overlapmap                   - Produces an interactive overlap map                                                                                                 - CLUS-RobClaMULT     - 2020 Feb 20
%   PDbdp                        - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                          - UTISTAT             - 2020 Aug 25
%   PDc                          - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator                - UTISTAT             - 2020 Feb 07
%   PDeff                        - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                     - UTISTAT             - 2020 Feb 11
%   PDpsi                        - Computes psi function (derivative of rho function) for minimum density power divergence estimator                                   - UTISTAT             - 2020 Feb 07
%   PDpsider                     - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                      - UTISTAT             - 2020 Apr 03
%   PDpsix                       - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                           - UTISTAT             - 2020 Feb 07
%   PDrho                        - Computes rho function for minimum density power divergence estimator                                                                - UTISTAT             - 2020 Feb 07
%   PDwei                        - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                               - UTISTAT             - 2020 Aug 25
%   position                     - Controls the position of the open figures                                                                                           - UTIGEN              - 2020 Oct 11
%   Powertra                     - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                             - UTISTAT             - 2020 Feb 05
%   publishBibliography          - Enables to create web page which contains the references inside the input .m files                                                  - UTIHELP             - 2020 Jun 05
%   publishFS                    - Enables to create automatic HELP FILES from structured .m function files                                                            - UTIHELP             - 2020 Jul 08
%   publishFunctionAlpha         - Enables to create web page which contains the alphabetical list of functions                                                        - UTIHELP             - 2020 Apr 21
%   publishFunctionCate          - Enables to create web page which contains the categorical list of functions                                                         - UTIHELP             - 2020 Jun 22
%   Qn                           - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                      - UTISTAT             - 2020 Jan 29
%   qqplotFS                     - Qqplot of studentized residuals with envelopes                                                                                      - VIS-Reg             - 2020 Feb 05
%   quickselectFS                - Finds the k-th order statistic                                                                                                      - UTIGEN              - 2020 Jan 29
%   RandIndexFS                  - Calculates Rand type Indices to compare two partitions                                                                              - UTISTAT             - 2020 Aug 07
%   randsampleFS                 - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                             - UTICOMB             - 2020 May 23
%   rcontFS                      - Generates a random two-way table with given marginal totals                                                                         - MULT-Categorical    - 2020 Feb 29
%   regressB                     - Computes Bayesian estimates of regression parameters                                                                                - REG-Bayes           - 2020 Feb 05
%   regressH                     - Fits a multiple linear regression model with heteroskedasticity                                                                     - REG-Hetero          - 2020 Feb 05
%   regressHart                  - Fits a multiple linear regression model using ART heteroskedasticity                                                                - REG-Hetero          - 2020 Feb 05
%   regressHhar                  - Fits a multiple linear regression model with Harvey heteroskedasticity                                                              - REG-Hetero          - 2020 May 07
%   regressts                    - Computes estimates of regression parameters for a time series models                                                                - REG-Regression      - 2020 Feb 29
%   removeExtraSpacesLF          - Removes extra spaces and selected carriage returns from input string                                                                - UTIGEN              - 2020 Jan 29
%   repDupValWithMean            - Replaces values of y which have non unique elements in vector x with local means                                                    - UTIGEN              - 2020 Jan 29
%   resfwdplot                   - Plots the trajectories of the monitored scaled (squared) residuals                                                                  - VIS-Reg             - 2020 Oct 12
%   resindexplot                 - Plots the residuals from a regression analysis versus index number or any other variable                                            - VIS-Reg             - 2020 Feb 18
%   restrdeter                   - Computes determinant restriction                                                                                                    - CLUS-RobClaMULT     - 2020 Feb 18
%   restrdeterGPCM               - Applies determinat restrictions for the 14 GPCM                                                                                     - CLUS-RobClaMULT     - 2020 Jul 15
%   restreigen                   - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2020 Jul 12
%   restreigeneasy               - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2020 Jul 12
%   restreigenmemopt             - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2020 Jul 12
%   restrshapeExact              - Computes constrained Gamma (shape) matrix with exact constraints                                                                    - CLUS-RobClaMULT     - 2020 Jul 15
%   restrshapeGPCM               - Produces the restricted shape matrix for the 14 GPCM                                                                                - CLUS-RobClaMULT     - 2020 Jul 15
%   restrSigmaGPCM               - Computes constrained covariance matrices for the 14 GPCM specifications                                                             - CLUS-RobClaMULT     - 2020 Aug 31
%   RKbdp                        - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                     - UTISTAT             - 2020 Jan 29
%   RKeff                        - Finds the constants c and M which are associated to the requested efficiency and ARP                                                - UTISTAT             - 2020 Apr 21
%   RKpsi                        - Computes psi function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT             - 2020 Jan 29
%   RKpsider                     - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                     - UTISTAT             - 2020 Jan 29
%   RKpsix                       - Computes psi function times x for Rocke (translated Tukey's) biweight                                                               - UTISTAT             - 2020 Jan 29
%   RKrho                        - Computes rho function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT             - 2020 Jan 29
%   RKwei                        - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                           - UTISTAT             - 2020 Jan 29
%   rlga                         - Performs robust linear grouping analysis                                                                                            - CLUS-RobClaREG      - 2020 Apr 21
%   rlsmo                        - Computes a running-lines smoother with global cross-validation                                                                      - REG-Transformations - 2020 Feb 05
%   RobCov                       - Computes covariance matrix of robust regression coefficients                                                                        - REG-Regression      - 2020 Sep 24
%   RobRegrSize                  - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size                        - REG-Regression      - 2020 Feb 05
%   rthin                        - Applies independent random thinning to a point pattern                                                                              - UTISTAT             - 2020 Jan 29
%   scatterboxplot               - Creates scatter diagram with marginal boxplots                                                                                      - VIS-Mult            - 2020 Oct 17
%   Score                        - Computes the score test for transformation                                                                                          - REG-Transformations - 2020 May 20
%   ScoreYJ                      - Computes the score test for Yeo and Johnson transformation                                                                          - REG-Transformations - 2020 May 20
%   ScoreYJall                   - Computes all the 4 score tests for YJ transformation                                                                                - REG-Transformations - 2020 May 23
%   ScoreYJmle                   - Computes the likelihood ratio test fof H_0=lambdaP=lambdaP0 and lambdaN=lambdaN0                                                    - REG-Transformations - 2020 May 20
%   ScoreYJpn                    - Computes the score test for YJ transformation for pos and neg observations                                                          - REG-Transformations - 2020 May 07
%   SDest                        - Computes Stahel-Donoho robust estimator of dispersion-location                                                                      - MULT-Multivariate   - 2020 Apr 21
%   SETARX                       - Implements Threshold autoregressive models with two regimes                                                                         - REG-Regression      - 2020 Aug 07
%   shuffling                    - Does a random permutation of the elements of input vector                                                                           - UTICOMB             - 2020 Jan 29
%   simdataset                   - Simulates and-or contaminates a dataset given the parameters of a finite mixture model with Gaussian components                     - CLUS-MixSim         - 2020 Feb 18
%   simdatasetreg                - Simulates a regression dataset given the parameters of a mixture regression model                                                   - CLUS-MixSim         - 2020 Feb 18
%   simulateLM                   - Simulates linear regression data with prespecified values of statistical indexes                                                    - REG-Regression      - 2020 May 07
%   simulateTS                   - Simulates a time series with trend, time varying seasonal, level shift and irregular component                                      - REG-Regression      - 2020 May 07
%   smothr                       - Produces smoothed values with constraints                                                                                           - REG-Transformations - 2020 Jan 29
%   Smult                        - Computes S estimators in multivariate analysis                                                                                      - MULT-Multivariate   - 2020 Feb 20
%   Smulteda                     - Computes S estimators in multivariate analysis for a series of values of bdp                                                        - MULT-Multivariate   - 2020 Jan 29
%   Sn                           - Robust estimator of scale (robust version of Gini's average difference)                                                             - UTISTAT             - 2020 Jan 29
%   SparseTableTest              - Computes independence test for large and sparse contingency tables                                                                  - MULT-Categorical    - 2020 Jan 29
%   spmplot                      - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal and possibly robust bivariate contours  - VIS-Mult            - 2020 Oct 12
%   Sreg                         - Computes S estimators in linear regression                                                                                          - REG-Regression      - 2020 Feb 20
%   Sregeda                      - Computes S estimators in linear regression for a series of values of bdp                                                            - REG-Regression      - 2020 Oct 10
%   subsets                      - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                    - UTICOMB             - 2020 May 22
%   suplabel                     - Places text as a title, xlabel, or ylabel on a group of subplots                                                                    - UTIGEN              - 2020 Apr 21
%   supsmu                       - Smooths scatterplots using Friedman's supersmoother algorithm                                                                       - REG-Transformations - 2020 Apr 25
%   tabulateFS                   - Creates frequency table of unique values of x, excluding possible 0 counts                                                          - UTISTAT             - 2020 Jan 29
%   Taureg                       - Computes Tau estimators in linear regression                                                                                        - REG-Regression      - 2020 Feb 20
%   TBbdp                        - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                - UTISTAT             - 2020 Oct 11
%   TBc                          - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                             - UTISTAT             - 2020 Jan 29
%   TBeff                        - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                           - UTISTAT             - 2020 Feb 05
%   tBothSides                   - Allows users to transform both sides of a (nonlinear) regression model                                                              - REG-Regression      - 2020 Feb 18
%   TBpsi                        - Computes psi function (derivative of rho function) for Tukey's biweight                                                             - UTISTAT             - 2020 Aug 30
%   TBpsider                     - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                        - UTISTAT             - 2020 Feb 05
%   TBpsix                       - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                     - UTISTAT             - 2020 Jan 29
%   TBrho                        - Computes rho function for Tukey's biweight                                                                                          - UTISTAT             - 2020 Jan 29
%   TBwei                        - Computes weight function psi(u)/u for Tukey's biweight                                                                              - UTISTAT             - 2020 Jan 29
%   tclust                       - Computes trimmed clustering with scatter restrictions                                                                               - CLUS-RobClaMULT     - 2020 Sep 25
%   tclusteda                    - Computes tclust for a series of values of the trimming factor                                                                       - CLUS-RobClaMULT     - 2020 Jul 16
%   tclustIC                     - Computes tclust for different number of groups k and restriction factors c                                                          - CLUS-RobClaMULT     - 2020 Jun 23
%   tclustICgpcm                 - Computes tclust for different number of groups k and restr. factors $c_{det}$ and $c_{shw}$                                         - CLUS-RobClaMULT     - 2020 Oct 12
%   tclustICplot                 - Plots information criterion as a function of c and k                                                                                - VIS-Clu             - 2020 Oct 10
%   tclustICplotGPCM             - Plots information criterion as a function of c and k                                                                                - VIS-Clu             - 2020 Oct 11
%   tclustICsol                  - Extracts a set of best relevant solutions                                                                                           - CLUS-RobClaMULT     - 2020 Sep 01
%   tclustICsolGPCM              - Extracts a set of best relevant solutions from 3D array computed using function tclustICgpcm                                        - CLUS-RobClaMULT     - 2020 Sep 21
%   tclustreg                    - Performs robust linear grouping analysis                                                                                            - CLUS-RobClaREG      - 2020 Aug 05
%   tclustregeda                 - Performs robust linear grouping analysis for a series of values of the trimming factor                                              - CLUS-RobClaREG      - 2020 Aug 27
%   tclustregIC                  - Computes tclustreg for different number of groups k and restriction factors c                                                       - CLUS-RobClaREG      - 2020 Sep 01
%   tkmeans                      - Computes trimmed k-means                                                                                                            - CLUS-RobClaMULT     - 2020 Apr 21
%   triu2vec                     - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                           - UTIGEN              - 2020 Apr 21
%   twdcdf                       - Computes the cumulative distribution function of the Tweedie distribution                                                           - UTISTAT             - 2020 May 06
%   twdpdf                       - Computes the probability density function of the Tweedie distribution                                                               - UTISTAT             - 2020 Apr 22
%   twdrnd                       - Generates random variates from the Tweedie distribution                                                                             - UTISTAT             - 2020 May 06
%   unibiv                       - Has the purpose of detecting univariate and bivariate outliers                                                                      - MULT-Multivariate   - 2020 Feb 18
%   upperfracpos                 - Positions two figures on the upper part of the screen                                                                               - UTIGEN              - 2020 Jan 29
%   verLessThanFS                - Compares version of MATLAB to specified version number                                                                              - UTIGEN              - 2020 May 26
%   vervaatrnd                   - Simulates random variates from the Vervaat perpetuity distribution                                                                  - UTISTAT             - 2020 Jan 29
%   vervaatsim                   - Returns a Vervaat perpetuity                                                                                                        - UTISTAT             - 2020 Jan 29
%   vervaatxdf                   - Returns the pdf and cdf of a Vervaat perpetuity                                                                                     - UTISTAT             - 2020 Feb 05
%   VIOM                         - Computes weights estimates under Variance-Inflation Model                                                                           - REG-Regression      - 2020 Apr 25
%   waterfallchart               - Creates a waterfall chart                                                                                                           - VIS-Mult            - 2020 Sep 28
%   wedgeplot                    - Generates the double wedge plot of a time series                                                                                    - VIS-Reg             - 2020 Jun 09
%   winsor                       - Returns a winsorized copy of input                                                                                                  - UTISTAT             - 2020 Apr 21
%   WNChygepdf                   - Returns Wallenius' non-central hypergeometric probability density values                                                            - UTISTAT             - 2020 Jan 29
%   wraptextFS                   - Formats long strings into wrapped text of specified width                                                                           - UTIGEN              - 2020 Feb 18
%   wthin                        - Thins a uni/bi-dimensional dataset                                                                                                  - UTISTAT             - 2020 Jul 04
%   xmlcreateFS                  - Creates an XML file passing through publishFS                                                                                       - UTIHELP             - 2020 Apr 21
%   yXplot                       - Produces an interactive scatterplot of y against each variable of X in the input dataset                                            - VIS-Reg             - 2020 Jun 22
%   zscoreFS                     - Computes (robust) standardized z scores                                                                                             - UTIGEN              - 2020 Jan 29
