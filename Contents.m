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
%   ace                          - Computes alternative conditional expectation                                                                                        - REG-Transformations - 2023 Apr 19
%   aceplot                      - Produces the aceplot to visualize the results of ace                                                                                - VIS-Reg             - 2023 May 10
%   add2boxplot                  - Adds labels to the boxplot figure                                                                                                   - VIS-Mult            - 2022 Dec 03
%   add2spm                      - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                       - VIS-Mult            - 2023 Apr 19
%   add2yX                       - Adds objects (personalized clickable multilegends and text labels) to the yXplot                                                    - VIS-Reg             - 2023 Apr 19
%   addt                         - Produces the t test for an additional explanatory variable                                                                          - REG-Regression      - 2022 Nov 10
%   ASbdp                        - Finds the constant c associated to the supplied breakdown point for Andrew's sine function                                          - UTISTAT             - 2023 May 29
%   ASc                          - Computes breakdown point and efficiency associated with constant c for Andrew's rho function                                        - UTISTAT             - 2023 May 28
%   ASeff                        - Finds the constant c which is associated to the requested efficiency for Andrew's sine function                                     - UTISTAT             - 2023 May 28
%   ASpsi                        - Computes psi function (derivative of rho function) for Andrew's sine function                                                       - UTISTAT             - 2023 May 28
%   ASpsider                     - Computes derivative of psi function (second derivative of rho function) for Andrew's sine function                                  - UTISTAT             - 2023 May 29
%   ASpsix                       - Computes psi function (derivative of rho function) times x for Andrew's sine function                                               - UTISTAT             - 2023 May 28
%   ASrho                        - Computes rho function for Andrew's sine function                                                                                    - UTISTAT             - 2023 May 28
%   ASwei                        - Computes weight function psi(u)/u for Andrew's sine function                                                                        - UTISTAT             - 2023 May 28
%   avas                         - Computes additivity and variance stabilization for regression                                                                       - REG-Transformations - 2023 Apr 19
%   avasms                       - Computes avas using a series of alternative options                                                                                 - REG-Transformations - 2023 May 17
%   avasmsplot                   - Produces the augmented star plot and enables interactivity                                                                          - VIS-Reg             - 2023 Apr 19
%   balloonplot                  - Creates a balloon plot of a contingency table                                                                                       - VIS-Mult            - 2023 Apr 19
%   barnardtest                  - Barnard's unconditional test                                                                                                        - MULT-Categorical    - 2023 Apr 19
%   barVariableWidth             - Produces a bar plot with different widths and colors for each bar                                                                   - VIS-Reg             - 2023 Apr 19
%   basicPower                   - Computes the basic power transformation                                                                                             - UTISTAT             - 2023 Apr 19
%   bc                           - Returns the Binomial coefficient                                                                                                    - UTICOMB             - 2022 Nov 10
%   biplotFS                     - Calls biplotAPP.mlapp to show a dynamic biplot                                                                                      - VIS-Mult            - 2023 Mar 09
%   boxcoxR                      - Finds MLE of lambda in linear regression (and confidence interval) using Box Cox, YJ or extended YJ  transformation                 - REG-Transformations - 2023 Apr 19
%   boxplotb                     - Computes a bivariate boxplot                                                                                                        - VIS-Mult            - 2023 Apr 19
%   boxtest                      - Performs Box test of equality of covariance matrices                                                                                - MULT-Multivariate   - 2023 Apr 19
%   brushFAN                     - Displays a GUI which enables brushing in the fanplot                                                                                - GUI                 - 2022 Nov 10
%   brushRES                     - Displays a GUI which enables brushing in resfwdplot                                                                                 - GUI                 - 2022 Nov 10
%   brushROB                     - Displays a GUI which enables brushing in resindexplot                                                                               - GUI                 - 2022 Nov 10
%   bwe                          - Estimates the bandwidth smoothing parameter for kernel density estimation                                                           - UTISTAT             - 2023 Apr 06
%   cabc                         - Closes all open figures except the one in foreground (the current)                                                                  - UTIGEN              - 2022 Nov 10
%   carbikeplot                  - Produces the carbike plot to find best relevant clustering solutions                                                                - VIS-Clu             - 2023 Apr 30
%   carbikeplotGPCM              - Produces the carbike plot to find best relevant clustering solutions                                                                - VIS-Clu             - 2023 Apr 30
%   cascade                      - Is a third party function used in FSDA demos and examples                                                                           - UTIGEN              - 2022 Nov 10
%   cdsplot                      - Produces the candlestick plot for robust model selection in linear regression                                                       - VIS-Reg             - 2023 Apr 19
%   CEVmodel                     - Computes price and instantaneous variance processes from the CEV model                                                              - UTISTAT             - 2022 Nov 10
%   clickableMultiLegend         - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                                 - UTIGEN              - 2022 Nov 10
%   ClusterRelabel               - Enables to control the labels of the clusters which contain predefined units                                                        - UTISTAT             - 2022 Nov 10
%   combsFS                      - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                          - UTICOMB             - 2022 Nov 10
%   CorAna                       - Performs correspondence analysis                                                                                                    - MULT-Categorical    - 2023 Apr 19
%   CorAnaplot                   - Draws the Correspondence Analysis (CA) graphs with confidence ellipses                                                              - VIS-Mult            - 2023 Apr 19
%   corrcdf                      - Computes correlation coefficient probability distribution function                                                                  - ProbDist            - 2022 Dec 02
%   corrNominal                  - Measures strength of association between two unordered (nominal) categorical variables                                              - MULT-Categorical    - 2023 Apr 19
%   corrOrdinal                  - Measures strength of association between two ordered categorical variables                                                          - MULT-Categorical    - 2023 Apr 19
%   corrpdf                      - Computes correlation coefficient probability density function                                                                       - ProbDist            - 2022 Dec 02
%   covplot                      - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                             - VIS-Mult            - 2023 Apr 19
%   CressieRead                  - Computes the power divergence family                                                                                                - MULT-Categorical    - 2023 Apr 19
%   crosstab2datamatrix          - Recreates the original data matrix X from contingency table N                                                                       - MULT-Categorical    - 2023 Apr 19
%   ctlcurves                    - Computes Classification Trimmed Likelihood Curves                                                                                   - CLUS-RobClaMULT     - 2023 Apr 19
%   ctlcurvesplot                - Plots the output of routine ctlcurves                                                                                               - VIS-Clu             - 2023 Apr 19
%   ctsub                        - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                          - UTISTAT             - 2023 Apr 06
%   dempk                        - Performs a merging of components found by tkmeans                                                                                   - CLUS-RobClaMULT     - 2023 Apr 19
%   distribspec                  - Plots a probability density function between specification limits                                                                   - UTISTAT             - 2023 May 03
%   ellipse                      - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                          - UTISTAT             - 2022 Nov 10
%   estregimeTAR                 - Estimate a regression model with OLS in one of the regimes of a TAR model                                                           - REG-Regression      - 2022 Nov 10
%   exactcdf                     - Finds exact p-values                                                                                                                - UTISTAT             - 2022 Nov 10
%   existFS                      - Check if file exists and puts answer in a cached persistent variable                                                                - UTIGEN              - 2022 Nov 10
%   fanBIC                       - Uses the output of FSRfan to choose the best value of the transformation parameter in linear regression                             - VIS-Reg             - 2023 Apr 19
%   fanBICpn                     - Uses the output of FSRfan called with input option family 'YJpn' to choose la_P and la_N                                            - VIS-Reg             - 2023 Apr 19
%   fanplot                      - Plots the fan plot for transformation in linear regression or deletion t stat                                                       - VIS-Reg             - 2023 Apr 28
%   FE_int_vol                   - Computes the integrated variance from a diffusion process via the Fourier estimator using Dirichlet kernel                          - UTISTAT             - 2023 Apr 19
%   FE_int_vol_Fejer             - Computes the integrated variance from a diffusion process via the Fourier estimator using Fejer kernel                              - UTISTAT             - 2023 Apr 19
%   FE_spot_vol                  - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel                          - UTISTAT             - 2023 Apr 19
%   FE_spot_vol_FFT              - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel, using the FFT algorithm - UTISTAT             - 2023 Apr 19
%   findDir                      - Finds recursively all directories in root                                                                                           - UTIGEN              - 2023 Apr 19
%   findFile                     - Finds recursively all files in root                                                                                                 - UTIGEN              - 2023 Apr 19
%   forecastH                    - Produce forecasts with confidence bands for regression model with heteroskedasticity                                                - REG-Hetero          - 2023 May 14
%   forecastTS                   - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                                   - REG-Regression      - 2023 Apr 19
%   FowlkesMallowsIndex          - Computes the Fowlkes and Mallows index                                                                                              - UTISTAT             - 2022 Nov 10
%   FSCorAnaeda                  - Performs forward search in correspondence analysis with exploratory data analysis purposes                                          - MULT-Multivariate   - 2023 Apr 19
%   FSCorAnaenvmmd               - Computes the empirical envelopes of Minimum MD outside subset during the search                                                     - MULT-Multivariate   - 2023 Apr 19
%   FSM                          - Computes forward search estimator in multivariate analysis                                                                          - MULT-Multivariate   - 2023 Apr 19
%   FSMbonfbound                 - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                          - UTISTAT             - 2022 Nov 10
%   FSMbsb                       - Gives the units belonging to subset at step(s) msel of the forward search                                                           - MULT-Multivariate   - 2023 Apr 19
%   FSMeda                       - Performs forward search in multivariate analysis with exploratory data analysis purposes                                            - MULT-Multivariate   - 2023 Apr 19
%   FSMedaeasy                   - Is exactly equal to FSMeda but it is much less efficient                                                                            - MULT-Multivariate   - 2023 Apr 19
%   FSMenvmmd                    - Computes the theoretical envelopes of Minimum MD outside subset during the search                                                   - MULT-Multivariate   - 2023 Apr 19
%   FSMfan                       - Computes confirmatory lrt of a suggested transformation                                                                             - MULT-Transformations- 2023 Apr 19
%   FSMinvmmd                    - Converts values of minimum Mahalanobis distance into confidence levels                                                              - MULT-Multivariate   - 2023 Apr 19
%   FSMmmd                       - Monitors minMD                                                                                                                      - MULT-Multivariate   - 2023 Apr 19
%   FSMmmdeasy                   - Is exactly equal to FSMmmd but it is much less efficient                                                                            - MULT-Multivariate   - 2023 Apr 19
%   FSMmmdrs                     - Performs random start monitoring of minimum Mahalanobis distance                                                                    - CLUS-RobClaMULT     - 2023 Apr 19
%   FSMtra                       - Computes MLE of transformation parameters                                                                                           - MULT-Transformations- 2023 Apr 19
%   FSR                          - Computes forward search estimator in linear regression                                                                              - REG-Regression      - 2023 Apr 19
%   FSRaddt                      - Produces t deletion tests for each explanatory variable                                                                             - REG-ModelSelection  - 2023 Apr 19
%   FSRB                         - Gives an automatic outlier detection procedure in Bayesian linear regression                                                        - REG-Bayes           - 2023 Apr 19
%   FSRBbsb                      - Returns the units belonging to the subset in each step of the Bayesian forward search                                               - REG-Regression      - 2023 Apr 19
%   FSRBeda                      - Enables to monitor several quantities in each step of the Bayesian search                                                           - REG-Bayes           - 2023 Apr 19
%   FSRBmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search                 - REG-Bayes           - 2023 Apr 19
%   FSRbonfbound                 - Computes Bonferroni bounds for each step of the search (in linear regression)                                                       - UTISTAT             - 2023 Apr 06
%   FSRBr                        - Bayesian forward search in linear regression reweighted                                                                             - REG-Bayes           - 2023 Apr 07
%   FSRbsb                       - Returns the units belonging to the subset in each step of the forward search                                                        - REG-Regression      - 2023 Apr 19
%   FSRcp                        - Monitors Cp and AIC for all models of interest of size smallp                                                                       - REG-ModelSelection  - 2023 Apr 19
%   FSReda                       - Enables to monitor several quantities in each step of the forward search                                                            - REG-Regression      - 2023 Apr 19
%   FSRenvmdr                    - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                    - REG-Regression      - 2023 Apr 19
%   FSRfan                       - Monitors the values of the score test statistic for each lambda                                                                     - REG-Transformations - 2023 Apr 19
%   FSRH                         - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                                 - REG-Hetero          - 2023 May 13
%   FSRHbsb                      - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                        - REG-Hetero          - 2023 May 13
%   FSRHeda                      - Enables to monitor several quantities in each step of the forward search                                                            - REG-Hetero          - 2023 May 14
%   FSRHmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search          - REG-Hetero          - 2023 May 13
%   FSRinvmdr                    - Converts values of minimum deletion residual into confidence levels                                                                 - REG-Regression      - 2023 Apr 19
%   FSRmdr                       - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                          - REG-Regression      - 2023 Apr 19
%   FSRmdrrs                     - Performs random start monitoring of minimum deletion residual                                                                       - CLUS-RobClaREG      - 2023 Apr 06
%   FSRms                        - Performs robust model selection using flexible trimming in linear regression                                                        - REG-ModelSelection  - 2023 Apr 19
%   FSRr                         - Forward search in linear regression reweighted                                                                                      - REG-Regression      - 2022 Nov 10
%   FSRts                        - Is an automatic adaptive procedure to detect outliers in time series                                                                - REG-Regression      - 2023 Apr 19
%   FSRtsbsb                     - Returns the units belonging to the subset in each step of the forward search                                                        - REG-Regression      - 2023 Apr 19
%   FSRtsmdr                     - Computes minimum deletion residual for time series models in each step of the search                                                - REG-Regression      - 2023 Apr 19
%   funnelchart                  - Displays a funnel chart                                                                                                             - VIS-Mult            - 2023 Apr 19
%   genSigmaGPCM                 - Generates covariance matrix for the 14 Gaussian Parsimonious Clustering Models                                                      - CLUS-RobClaMULT     - 2022 Nov 10
%   GowerIndex                   - Computes matrix of similarity indexes using Gower metric                                                                            - CLUS-RobClaMULT     - 2023 Apr 19
%   GUIautocorr                  - Shows the necessary calculations to obtain the autocorrelation in a GUI                                                             - GUI                 - 2023 Apr 19
%   GUIconcentration             - Shows the necessary calculations to obtain the GINI concentration index in a GUI                                                    - GUI                 - 2023 Apr 19
%   GUIcorr                      - Shows the necessary calculations to obtain the correlation in a GUI                                                                 - GUI                 - 2023 Apr 06
%   GUIcov                       - Shows the necessary calculations to obtain the covariance in a GUI                                                                  - GUI                 - 2023 Apr 06
%   GUImad                       - Shows the necessary calculations to obtain MAD, S_M or S_Me in a GUI                                                                - GUI                 - 2023 Apr 19
%   GUIpowermean                 - Shows the necessary calculations to obtain the power (generalized) mean in a GUI                                                    - GUI                 - 2023 Apr 06
%   GUIquantile                  - Shows the necessary calculations to obtain $x_z$ quantile                                                                           - GUI                 - 2023 Apr 19
%   GUIregress                   - Shows the necessary calculations to obtain simple linear regression statistics in a GUI                                             - GUI                 - 2023 May 01
%   GUIskewness                  - Shows the necessary calculations to obtain the variance in a GUI                                                                    - GUI                 - 2023 Apr 06
%   GUIstd                       - Shows the necessary calculations to obtain the standard deviation in a GUI                                                          - GUI                 - 2023 Apr 06
%   GUItrimmean                  - Shows the necessary calculations to obtain the trimmed mean in a GUI                                                                - GUI                 - 2023 Apr 06
%   GUIvar                       - Shows the necessary calculations to obtain the variance in a GUI                                                                    - GUI                 - 2023 Apr 06
%   GYfilt                       - Computes the Gervini-Yohai univariate outlier identifier                                                                            - UTISTAT             - 2023 Apr 19
%   HAbdp                        - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT             - 2022 Nov 10
%   HAc                          - Computes breakdown point and efficiency associated with constant c                                                                  - UTISTAT             - 2023 Apr 19
%   HAeff                        - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                        - UTISTAT             - 2022 Nov 10
%   HApsi                        - Computes psi function  using Hampel proposal                                                                                        - UTISTAT             - 2022 Nov 10
%   HApsider                     - Computes derivative of psi function  using Hampel proposal                                                                          - UTISTAT             - 2022 Nov 10
%   HApsix                       - Computes psi function  using Hampel proposal times x                                                                                - UTISTAT             - 2022 Nov 10
%   HArho                        - Computes rho function  using Hampel proposal                                                                                        - UTISTAT             - 2022 Nov 10
%   HAwei                        - Computes weight function psi(u)/u using Hampel proposal                                                                             - UTISTAT             - 2022 Nov 10
%   histFS                       - Plots a histogram with the elements in each bin grouped according to a vector of labels                                             - VIS-Reg             - 2022 Nov 10
%   htmlwriteFS                  - Is an obsolete function which will be removed in future releases. Use publishFS.m instead                                           - UTIHELP             - 2023 Apr 19
%   HUc                          - Computes breakdown point and efficiency associated with constant c for Huber link                                                   - UTISTAT             - 2023 May 28
%   HUeff                        - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                   - UTISTAT             - 2023 May 28
%   HUpsi                        - Computes psi function (derivative of rho function) for Huber                                                                        - UTISTAT             - 2022 Nov 10
%   HUpsider                     - Computes derivative of psi function (second derivative of rho function) for Huber                                                   - UTISTAT             - 2022 Nov 10
%   HUpsix                       - Computes psi function (derivative of rho function) times x for Huber                                                                - UTISTAT             - 2022 Nov 10
%   HUrho                        - Computes rho function for Huber                                                                                                     - UTISTAT             - 2022 Nov 10
%   HUwei                        - Computes weight function psi(u)/u for Huber                                                                                         - UTISTAT             - 2023 Apr 06
%   HYPbdp                       - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                       - UTISTAT             - 2023 May 28
%   HYPc                         - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)     - UTISTAT             - 2023 Apr 19
%   HYPck                        - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                             - UTISTAT             - 2022 Nov 10
%   HYPeff                       - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                           - UTISTAT             - 2022 Nov 10
%   HYPk                         - Computes breakdown point and efficiency for hyp. tan. estimator                                                                     - UTISTAT             - 2023 Apr 19
%   HYPpsi                       - Computes psi function for hyperbolic tangent estimator                                                                              - UTISTAT             - 2022 Nov 10
%   HYPpsider                    - Computes derivative of psi function for hyperbolic tangent estimator                                                                - UTISTAT             - 2022 Nov 10
%   HYPpsix                      - Computes psi function for hyperbolic tangent estimator times x                                                                      - UTISTAT             - 2022 Nov 10
%   HYPrho                       - Computes rho function  using hyperbolic tangent estimator                                                                           - UTISTAT             - 2022 Nov 10
%   HYPwei                       - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                  - UTISTAT             - 2022 Nov 10
%   inversegamcdf                - Computes inverse-gamma cumulative distribution function                                                                             - ProbDist            - 2022 Dec 02
%   inversegaminv                - Computes the inverse of the inverse-gamma cumulative distribution function                                                          - ProbDist            - 2022 Dec 02
%   inversegampdf                - Computes inverse-gamma probability density function                                                                                 - ProbDist            - 2022 Dec 02
%   isfunction                   - Checks if a function exists                                                                                                         - UTIGEN              - 2022 Nov 10
%   kdebiv                       - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                      - UTISTAT             - 2022 Nov 10
%   levfwdplot                   - Plots the trajectories of leverage along the search                                                                                 - VIS-Reg             - 2023 Apr 19
%   lexunrank                    - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                        - UTICOMB             - 2022 Nov 10
%   lga                          - Performs linear grouping analysis                                                                                                   - CLUS-RobClaREG      - 2023 Apr 19
%   logmvnpdfFS                  - Produces log of Multivariate normal probability density function (pdf)                                                              - ProbDist            - 2022 Dec 02
%   LTSts                        - Extends LTS estimator to time series                                                                                                - REG-Regression      - 2023 Apr 19
%   LTStsLSmult                  - Extends LTSts to the detection of multiple Level Shifts in time series                                                              - REG-Regression      - 2023 Apr 19
%   LTStsVarSel                  - Does variable selection in the robust time series model LTSts                                                                       - REG-Regression      - 2023 Apr 19
%   LXS                          - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                                - REG-Regression      - 2023 Apr 19
%   mahalCorAna                  - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT             - 2022 Dec 03
%   mahalFS                      - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT             - 2022 Nov 10
%   makecontentsfileFS           - Extends Matlab function makecontentsfile                                                                                            - UTIHELP             - 2023 Apr 19
%   malfwdplot                   - Plots the trajectories of scaled Mahalanobis distances along the search                                                             - VIS-Mult            - 2023 Apr 19
%   malindexplot                 - Plots the Mahalanobis distances versus a selected variable                                                                          - VIS-Mult            - 2023 May 01
%   mcd                          - Computes Minimum Covariance Determinant                                                                                             - MULT-Multivariate   - 2023 May 10
%   mcdCorAna                    - Computes Minimum Covariance Determinant in correspondence analysis                                                                  - MULT-Multivariate   - 2023 Apr 19
%   mcdeda                       - Monitors Minimum Covariance Determinant for a series of values of bdp                                                               - MULT-Multivariate   - 2023 Apr 19
%   mdpattern                    - Finds and plots missing data patterns                                                                                               - UTISTAT             - 2023 Apr 19
%   mdpd                         - Computes Minimum Distance Power Divergence statistics                                                                               - UTISTAT             - 2023 Apr 19
%   mdpdR                        - Allows to apply Minimum Density Power Divergence criterion to parametric regression problems                                        - REG-Regression      - 2023 Apr 19
%   mdpdReda                     - Allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems                                     - REG-Regression      - 2023 Apr 19
%   mdrplot                      - Plots the trajectory of minimum deletion residual (mdr)                                                                             - VIS-Reg             - 2023 Apr 19
%   mdrrsplot                    - Plots the trajectory of minimum deletion residual from random starts                                                                - CLUS-RobClaREG      - 2023 Apr 19
%   medcouple                    - Computes the medcouple, a robust skewness estimator                                                                                 - UTISTAT             - 2022 Nov 10
%   MixSim                       - Generates k clusters in v dimensions with given overlap                                                                             - CLUS-MixSim         - 2023 Apr 19
%   MixSimreg                    - Generates k regression hyperplanes in p dimensions with given overlap                                                               - CLUS-MixSim         - 2023 Apr 19
%   Mlocation                    - Finds the M estimator of location in a univariate sample                                                                            - UTISTAT             - 2022 Nov 10
%   Mlocsca                      - Computes M estimator of location and scale in univariate samples                                                                    - UTISTAT             - 2022 Nov 10
%   mmdplot                      - Plots the trajectory of minimum Mahalanobis distance (mmd)                                                                          - VIS-Mult            - 2023 Apr 19
%   mmdrsplot                    - Plots the trajectories of minimum Mahalanobis distances from different starting points                                              - VIS-Mult            - 2023 Apr 19
%   MMmult                       - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                              - MULT-Multivariate   - 2023 Apr 19
%   MMmultcore                   - Computes multivariate MM estimators for a selected fixed scale                                                                      - MULT-Multivariate   - 2023 Apr 19
%   MMmulteda                    - Computes MM estimators in multivariate analysis for a series of values of eff                                                       - MULT-Multivariate   - 2023 Apr 19
%   MMreg                        - Computes MM estimator of regression coefficients                                                                                    - REG-Regression      - 2023 May 28
%   MMregcore                    - Computes MM regression estimators for a selected fixed scale                                                                        - REG-Regression      - 2023 May 29
%   MMregeda                     - Computes MM estimator in linear regression for a series of values of efficiency                                                     - REG-Regression      - 2023 May 29
%   moonplot                     - Draws the Correspondence Analysis (CA) moonplot                                                                                     - VIS-Mult            - 2023 Apr 19
%   mreadFS                      - Enables to create a structure with InputArgs/OptArgs/OutArgs ... from .m function files (OBSOLETE FUNCTION REPLACED BY publishFS.m) - UTIHELP             - 2023 Apr 19
%   Mscale                       - Finds the M estimator of the scale                                                                                                  - UTISTAT             - 2023 May 29
%   msdcutoff                    - Mahalanobis Squared Distance cutoff                                                                                                 - UTISTAT             - 2023 May 10
%   mtR                          - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                            - ProbDist            - 2022 Dec 02
%   mve                          - Computes Minimum volume ellipsoid                                                                                                   - MULT-Multivariate   - 2023 Apr 19
%   mveeda                       - Monitors Minimum volume ellipsoid for a series of values of bdp                                                                     - MULT-Multivariate   - 2023 Apr 19
%   nchoosekFS                   - Returns the Binomial coefficient or matrix containing all combinations                                                              - UTICOMB             - 2022 Nov 10
%   ncpci                        - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                    - UTISTAT             - 2023 Apr 19
%   ncx2mixtcdf                  - Computes cumulative distribution function of a linear combination of non-central chi-square +N(0,s)                                 - ProbDist            - 2023 Apr 19
%   normBoxCox                   - Computes (normalized) Box-Cox transformation                                                                                        - UTISTAT             - 2023 Apr 19
%   normYJ                       - Computes (normalized) Yeo-Johnson transformation                                                                                    - UTISTAT             - 2023 Apr 19
%   normYJpn                     - Computes (normalized) extended Yeo-Johnson transformation                                                                           - UTISTAT             - 2023 Apr 19
%   openMatlabFileFromHTML       - Enables to put in HTML an hypertextual link to a specific MATLAB file                                                               - UTIGEN              - 2022 Nov 10
%   OPTbdp                       - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT             - 2022 Nov 10
%   OPTc                         - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                         - UTISTAT             - 2022 Nov 10
%   OPTeff                       - Finds the constant c which is associated to the requested efficiency                                                                - UTISTAT             - 2023 May 28
%   OptimalCuttingFrequency      - Computes the optimal cutting frequency for the Fourier estimator of integrated variance                                             - UTISTAT             - 2023 Apr 03
%   OPTpsi                       - Computes psi function (derivative of rho function) for optimal weight function                                                      - UTISTAT             - 2022 Nov 10
%   OPTpsider                    - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                 - UTISTAT             - 2022 Nov 10
%   OPTpsix                      - Computes psi function (derivative of rho function) times x                                                                          - UTISTAT             - 2022 Nov 10
%   OPTrho                       - Computes rho function for optimal weight function                                                                                   - UTISTAT             - 2022 Nov 10
%   OPTwei                       - Computes weight function psi(u)/u for optimal weight function                                                                       - UTISTAT             - 2022 Nov 10
%   overlap                      - Computes the exact overlap given the parameters of the mixture                                                                      - CLUS-MixSim         - 2022 Nov 10
%   overlapmap                   - Produces an interactive overlap map                                                                                                 - CLUS-RobClaMULT     - 2023 Apr 19
%   pcaFS                        - Performs Principal Component Analysis (PCA) on raw data                                                                             - MULT-Multivariate   - 2023 Apr 19
%   PDbdp                        - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                          - UTISTAT             - 2023 Apr 07
%   PDc                          - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator                - UTISTAT             - 2023 May 28
%   PDeff                        - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                     - UTISTAT             - 2023 May 28
%   PDpsi                        - Computes psi function (derivative of rho function) for minimum density power divergence estimator                                   - UTISTAT             - 2023 Apr 07
%   PDpsider                     - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                      - UTISTAT             - 2023 Apr 07
%   PDpsix                       - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                           - UTISTAT             - 2023 Apr 07
%   PDrho                        - Computes rho function for minimum density power divergence estimator                                                                - UTISTAT             - 2023 Apr 07
%   PDwei                        - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                               - UTISTAT             - 2023 Apr 07
%   position                     - Controls the position of the open figures                                                                                           - UTIGEN              - 2023 Jan 09
%   Powertra                     - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                             - UTISTAT             - 2023 Apr 19
%   publishBibliography          - Enables to create web page which contains the references inside the input .m files                                                  - UTIHELP             - 2023 Apr 19
%   publishFS                    - Enables to create automatic HELP FILES from structured .m function files                                                            - UTIHELP             - 2023 Apr 19
%   publishFunctionAlpha         - Enables to create web page which contains the alphabetical list of functions                                                        - UTIHELP             - 2023 Apr 19
%   publishFunctionCate          - Enables to create web page which contains the categorical list of functions                                                         - UTIHELP             - 2023 Apr 19
%   Qn                           - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                      - UTISTAT             - 2023 Apr 06
%   qqplotFS                     - Creates qqplot of studentized residuals with envelopes                                                                              - VIS-Reg             - 2023 Apr 19
%   quickselectFS                - Finds the k-th order statistic                                                                                                      - UTIGEN              - 2023 May 10
%   quickselectFS_demo           - Illustrates the functioning of quickselectFS                                                                                        - GUI                 - 2023 May 10
%   quickselectFSw               - Finds the 100*p-th weighted order statistic for $0<p<1$                                                                             - UTIGEN              - 2023 May 10
%   quickselectFSw_demo          - Illustrates the functioning of quickselectFSw                                                                                       - GUI                 - 2023 May 10
%   RandIndexFS                  - Calculates Rand type Indices to compare two partitions                                                                              - UTISTAT             - 2022 Nov 10
%   randsampleFS                 - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                             - UTICOMB             - 2023 Apr 06
%   rcontFS                      - Generates a random two-way table with given marginal totals                                                                         - MULT-Categorical    - 2023 Apr 19
%   regressB                     - Computes Bayesian estimates of regression parameters                                                                                - REG-Bayes           - 2023 Apr 19
%   regressH                     - Fits a multiple linear regression model with heteroskedasticity                                                                     - REG-Hetero          - 2023 May 13
%   regressHart                  - Fits a multiple linear regression model using ART heteroskedasticity                                                                - REG-Hetero          - 2023 May 14
%   regressHhar                  - Fits a multiple linear regression model with Harvey heteroskedasticity                                                              - REG-Hetero          - 2023 May 14
%   regressts                    - Computes estimates of regression parameters for a time series models                                                                - REG-Regression      - 2023 Apr 19
%   removeExtraSpacesLF          - Removes extra spaces and selected carriage returns from input string                                                                - UTIGEN              - 2022 Nov 10
%   repDupValWithMean            - Replaces values of y which have non unique elements in vector x with local means                                                    - UTIGEN              - 2023 Apr 19
%   resfwdplot                   - Plots the trajectories of the monitored scaled (squared) residuals                                                                  - VIS-Reg             - 2023 Apr 19
%   resindexplot                 - Plots the residuals from a regression analysis versus index number or any other variable                                            - VIS-Reg             - 2023 Apr 19
%   restrdeter                   - Computes determinant restriction                                                                                                    - CLUS-RobClaMULT     - 2022 Nov 10
%   restrdeterGPCM               - Applies determinat restrictions for the 14 GPCM                                                                                     - CLUS-RobClaMULT     - 2023 Apr 06
%   restreigen                   - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2022 Nov 10
%   restreigeneasy               - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2022 Nov 10
%   restreigenmemopt             - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2022 Nov 10
%   restrshapeExact              - Computes constrained Gamma (shape) matrix with exact constraints                                                                    - CLUS-RobClaMULT     - 2022 Nov 10
%   restrshapeGPCM               - Produces the restricted shape matrix for the 14 GPCM                                                                                - CLUS-RobClaMULT     - 2023 Apr 06
%   restrSigmaGPCM               - Computes constrained covariance matrices for the 14 GPCM specifications                                                             - CLUS-RobClaMULT     - 2023 Apr 06
%   RhoPsiWei                    - Finds rho, psi, psi', w functions given bdp, or eff or tuning constant c                                                            - UTISTAT             - 2023 Apr 19
%   RKbdp                        - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                     - UTISTAT             - 2023 Apr 06
%   RKeff                        - Finds the constants c and M which are associated to the requested efficiency and ARP                                                - UTISTAT             - 2022 Nov 10
%   RKpsi                        - Computes psi function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT             - 2023 Apr 06
%   RKpsider                     - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                     - UTISTAT             - 2023 Apr 06
%   RKpsix                       - Computes psi function times x for Rocke (translated Tukey's) biweight                                                               - UTISTAT             - 2023 Apr 06
%   RKrho                        - Computes rho function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT             - 2023 Apr 06
%   RKwei                        - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                           - UTISTAT             - 2023 Apr 06
%   rlga                         - Performs robust linear grouping analysis                                                                                            - CLUS-RobClaREG      - 2023 Apr 19
%   rlsmo                        - Computes a running-lines smoother with global cross-validation                                                                      - REG-Transformations - 2023 Apr 06
%   RobCov                       - Computes covariance matrix of robust regression coefficients                                                                        - REG-Regression      - 2023 Apr 19
%   RobRegrSize                  - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size                        - REG-Regression      - 2022 Nov 10
%   rthin                        - Applies independent random thinning to a point pattern                                                                              - UTISTAT             - 2022 Nov 10
%   scatterboxplot               - Creates scatter diagram with marginal boxplots                                                                                      - VIS-Mult            - 2023 Apr 19
%   Score                        - Computes the score test for transformation                                                                                          - REG-Transformations - 2023 Apr 19
%   ScoreYJ                      - Computes the score test for Yeo and Johnson transformation                                                                          - REG-Transformations - 2023 Apr 19
%   ScoreYJall                   - Computes all the 4 score tests for YJ transformation                                                                                - REG-Transformations - 2023 Apr 19
%   ScoreYJmle                   - Computes the likelihood ratio test fof H_0=lambdaP=lambdaP0 and lambdaN=lambdaN0                                                    - REG-Transformations - 2023 Apr 19
%   ScoreYJpn                    - Computes the score test for YJ transformation for pos and neg observations                                                          - REG-Transformations - 2023 Apr 19
%   SDest                        - Computes Stahel-Donoho robust estimator of dispersion-location                                                                      - MULT-Multivariate   - 2023 Apr 19
%   SETARX                       - Implements Threshold autoregressive models with two regimes                                                                         - REG-Regression      - 2023 Apr 06
%   setdiffFS                    - Finds the positive integers in a which are not present in the positive integers in b                                                - UTIGEN              - 2022 Nov 10
%   setToolboxStartEnd           - Sets release compatibility in ToolboxPackagingConfiguration.prj file                                                                - UTIHELP             - 2023 Apr 19
%   shuffling                    - Does a random permutation of the elements of input vector                                                                           - UTICOMB             - 2023 Apr 06
%   simdataset                   - Simulates and-or contaminates a dataset given the parameters of a finite mixture model with Gaussian components                     - CLUS-MixSim         - 2023 Apr 19
%   simdatasetreg                - Simulates a regression dataset given the parameters of a mixture regression model                                                   - CLUS-MixSim         - 2023 Apr 19
%   simulateLM                   - Simulates linear regression data with prespecified values of statistical indexes                                                    - REG-Regression      - 2023 Apr 19
%   simulateTS                   - Simulates a time series with trend, time varying seasonal, level shift and irregular component                                      - REG-Regression      - 2023 Apr 19
%   smothr                       - Produces smoothed values with constraints                                                                                           - REG-Transformations - 2023 Apr 06
%   Smult                        - Computes S estimators in multivariate analysis                                                                                      - MULT-Multivariate   - 2023 Apr 19
%   Smulteda                     - Computes S estimators in multivariate analysis for a series of values of bdp                                                        - MULT-Multivariate   - 2023 Apr 19
%   Sn                           - Robust estimator of scale (robust version of Gini's average difference)                                                             - UTISTAT             - 2023 Apr 06
%   SparseTableTest              - Computes independence test for large and sparse contingency tables                                                                  - MULT-Categorical    - 2023 Apr 19
%   spmplot                      - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal and possibly robust bivariate contours  - VIS-Mult            - 2023 Apr 19
%   Sreg                         - Computes S estimators in linear regression                                                                                          - REG-Regression      - 2023 May 29
%   Sregeda                      - Computes S estimators in linear regression for a series of values of bdp                                                            - REG-Regression      - 2023 May 29
%   subsets                      - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                    - UTICOMB             - 2022 Nov 10
%   suplabel                     - Places text as a title, xlabel, or ylabel on a group of subplots                                                                    - UTIGEN              - 2022 Nov 10
%   supsmu                       - Smooths scatterplots using Friedman's supersmoother algorithm                                                                       - REG-Transformations - 2022 Nov 10
%   tabledisp                    - Displays in good format a table or array in the command line                                                                        - UTIGEN              - 2023 Mar 27
%   tabulateFS                   - Creates frequency table of unique values of x, excluding possible 0 counts                                                          - UTISTAT             - 2022 Nov 10
%   Taureg                       - Computes Tau estimators in linear regression                                                                                        - REG-Regression      - 2023 May 29
%   TBbdp                        - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                - UTISTAT             - 2023 May 27
%   TBc                          - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                             - UTISTAT             - 2023 May 28
%   TBeff                        - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                           - UTISTAT             - 2023 May 28
%   tBothSides                   - Allows users to transform both sides of a (nonlinear) regression model                                                              - REG-Regression      - 2023 Apr 19
%   TBpsi                        - Computes psi function (derivative of rho function) for Tukey's biweight                                                             - UTISTAT             - 2023 Apr 06
%   TBpsider                     - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                        - UTISTAT             - 2022 Nov 10
%   TBpsix                       - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                     - UTISTAT             - 2022 Nov 10
%   TBrho                        - Computes rho function for Tukey's biweight                                                                                          - UTISTAT             - 2023 Apr 06
%   TBwei                        - Computes weight function psi(u)/u for Tukey's biweight                                                                              - UTISTAT             - 2023 Apr 06
%   tclust                       - Computes trimmed clustering with scatter restrictions                                                                               - CLUS-RobClaMULT     - 2023 Apr 19
%   tclusteda                    - Computes tclust for a series of values of the trimming factor                                                                       - CLUS-RobClaMULT     - 2023 Apr 19
%   tclustIC                     - Computes tclust for different number of groups k and restriction factors c                                                          - CLUS-RobClaMULT     - 2023 Apr 19
%   tclustICgpcm                 - Computes tclust for different number of groups k and restr. factors $c_{det}$ and $c_{shw}$                                         - CLUS-RobClaMULT     - 2023 Apr 19
%   tclustICplot                 - Plots information criterion as a function of c and k                                                                                - VIS-Clu             - 2023 Apr 19
%   tclustICplotGPCM             - Plots information criterion as a function of  $c_{det}$, $c_{shw}$,  $c_{shb}$ and $k$                                              - VIS-Clu             - 2023 Apr 19
%   tclustICsol                  - Extracts a set of best relevant solutions                                                                                           - CLUS-RobClaMULT     - 2023 Apr 19
%   tclustICsolGPCM              - Extracts a set of best relevant solutions from 3D array computed using function tclustICgpcm                                        - CLUS-RobClaMULT     - 2023 Apr 19
%   tclustreg                    - Performs robust linear grouping analysis                                                                                            - CLUS-RobClaREG      - 2023 Apr 07
%   tclustregeda                 - Performs robust linear grouping analysis for a series of values of the trimming factor                                              - CLUS-RobClaREG      - 2023 Apr 07
%   tclustregIC                  - Computes tclustreg for different number of groups k and restriction factors c                                                       - CLUS-RobClaREG      - 2023 Apr 19
%   tkmeans                      - Computes trimmed k-means                                                                                                            - CLUS-RobClaMULT     - 2023 Apr 19
%   triu2vec                     - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                           - UTIGEN              - 2022 Nov 10
%   twdcdf                       - Computes the cumulative distribution function of the Tweedie distribution                                                           - ProbDist            - 2022 Dec 02
%   twdpdf                       - Computes the probability density function of the Tweedie distribution                                                               - ProbDist            - 2022 Dec 02
%   twdrnd                       - Generates random variates from the Tweedie distribution                                                                             - ProbDist            - 2022 Dec 02
%   txmerge                      - Performs a (hierarchical) merging of the inflated number of components found by tkmeans or TCLUST                                   - CLUS-RobClaMULT     - 2023 May 29
%   unibiv                       - Has the purpose of detecting univariate and bivariate outliers                                                                      - MULT-Multivariate   - 2023 Apr 19
%   upperfracpos                 - Positions two figures on the upper part of the screen                                                                               - UTIGEN              - 2022 Nov 10
%   verLessThanFS                - Compares version of MATLAB to specified version number                                                                              - UTIGEN              - 2022 Nov 10
%   vervaatrnd                   - Simulates random variates from the Vervaat perpetuity distribution                                                                  - ProbDist            - 2022 Dec 02
%   vervaatsim                   - Returns a Vervaat perpetuity                                                                                                        - ProbDist            - 2022 Dec 02
%   vervaatxdf                   - Returns the pdf and cdf of a Vervaat perpetuity                                                                                     - ProbDist            - 2022 Dec 02
%   VIOM                         - Computes weights estimates under Variance-Inflation Model                                                                           - REG-Regression      - 2023 Apr 19
%   waterfallchart               - Creates a waterfall chart                                                                                                           - VIS-Mult            - 2023 Apr 19
%   wedgeplot                    - Generates the double wedge plot of a time series                                                                                    - VIS-Reg             - 2023 Apr 19
%   winsor                       - Returns a winsorized copy of input                                                                                                  - UTISTAT             - 2022 Nov 10
%   WNChygecdf                   - Returns Wallenius' non-central hypergeometric cumulative distribution function                                                      - ProbDist            - 2023 Apr 06
%   WNChygepdf                   - Returns Wallenius' non-central hypergeometric probability density function                                                          - ProbDist            - 2023 Apr 06
%   wraptextFS                   - Formats long strings into wrapped text of specified width                                                                           - UTIGEN              - 2023 Apr 19
%   wthin                        - Thins a uni/bi-dimensional dataset                                                                                                  - UTISTAT             - 2022 Nov 10
%   xmlcreateFS                  - Creates an XML file passing through publishFS                                                                                       - UTIHELP             - 2023 Apr 19
%   yXplot                       - Produces an interactive scatterplot of y against each variable of X in the input dataset                                            - VIS-Reg             - 2023 May 10
%   zscoreFS                     - Computes (robust) standardized z scores                                                                                             - UTIGEN              - 2022 Nov 10
