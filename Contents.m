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
%   ace                          - Computes alternative conditional expectation                                                                                        - REG-Transformations - 2021 Oct 16
%   aceplot                      - Produces the aceplot to visualize the results of ace                                                                                - VIS-Reg             - 2022 Feb 28
%   add2boxplot                  - Add labels to the boxplot figure                                                                                                    - VIS-Mult            - 2022 Aug 21
%   add2spm                      - Adds objects (personalized clickable multilegends and text labels) to the scatter plot matrix                                       - VIS-Mult            - 2021 Feb 01
%   add2yX                       - Adds objects (personalized clickable multilegends and text labels) to the yXplot                                                    - VIS-Reg             - 2021 Feb 01
%   addt                         - Produces the t test for an additional explanatory variable                                                                          - REG-Regression      - 2021 Aug 10
%   avas                         - Computes additivity and variance stabilization for regression                                                                       - REG-Transformations - 2022 Feb 23
%   avasms                       - Computes avas using a series of alternative options                                                                                 - REG-Transformations - 2022 Mar 14
%   balloonplot                  - Creates a balloon plot of a contingency table                                                                                       - VIS-Mult            - 2022 Aug 23
%   barnardtest                  - Barnard's unconditional test                                                                                                        - MULT-Categorical    - 2021 Feb 01
%   barVariableWidth             - Produces a bar plot with different widths and colors for each bar                                                                   - VIS-Reg             - 2021 Feb 19
%   basicPower                   - Computes the basic power transformation                                                                                             - UTISTAT             - 2021 Feb 01
%   bc                           - Returns the Binomial coefficient                                                                                                    - UTICOMB             - 2021 Feb 01
%   biplotFS                     - Calls biplotAPP.mlapp to show a dynamic biplot                                                                                      - VIS-Mult            - 2021 Feb 01
%   boxcoxR                      - Finds MLE of lambda in linear regression (and confidence interval) using Box Cox, YJ or extended YJ  transformation                 - REG-Transformations - 2021 Feb 22
%   boxplotb                     - Computes a bivariate boxplot                                                                                                        - VIS-Mult            - 2022 Sep 26
%   boxtest                      - Performs Box test of equality of covariance matrices                                                                                - MULT-Multivariate   - 2021 Feb 01
%   brushFAN                     - Displays a GUI which enables brushing in the fanplot                                                                                - GUI                 - 2021 Feb 01
%   brushRES                     - Displays a GUI which enables brushing in resfwdplot                                                                                 - GUI                 - 2021 Feb 01
%   brushROB                     - Displays a GUI which enables brushing in resindexplot                                                                               - GUI                 - 2021 Feb 01
%   bwe                          - Estimates the bandwidth smoothing parameter for kernel density estimation                                                           - UTISTAT             - 2021 Feb 01
%   cabc                         - Closes all open figures except the one in foreground (the current)                                                                  - UTIGEN              - 2021 Feb 01
%   carbikeplot                  - Produces the carbike plot to find best relevant clustering solutions                                                                - VIS-Clu             - 2021 Jun 10
%   carbikeplotGPCM              - Carbikeplot produces the carbike plot to find best relevant clustering solutions                                                    - VIS-Clu             - 2021 Feb 01
%   cascade                      - Is a third party function used in FSDA demos and examples                                                                           - UTIGEN              - 2021 Feb 01
%   cdsplot                      - Produces the candlestick plot for robust model selection in linear regression                                                       - VIS-Reg             - 2021 Feb 01
%   CEVmodel                     - Computes price and instantaneous variance processes from the CEV model                                                              - UTISTAT             - 2021 Mar 01
%   clickableMultiLegend         - Hides/shows symbols inside all gplotmatrix subplots (or similar multi-plots) clicking on the legend                                 - UTIGEN              - 2021 Feb 01
%   ClusterRelabel               - Enables to control the labels of the clusters which contain predefined units                                                        - UTISTAT             - 2021 Feb 01
%   combsFS                      - Is an iterative algorithm equivalent to the MATLAB combs.m                                                                          - UTICOMB             - 2021 Jul 06
%   CorAna                       - Performs correspondence analysis                                                                                                    - MULT-Categorical    - 2021 Dec 04
%   CorAnaplot                   - Draws the Correspondence Analysis (CA) graphs with confidence ellipses                                                              - VIS-Mult            - 2021 Dec 13
%   corrcdf                      - Correlation coefficient probability distribution function                                                                           - UTISTAT             - 2022 Nov 09
%   corrNominal                  - Measures strength of association between two unordered (nominal) categorical variables                                              - MULT-Categorical    - 2021 Jun 24
%   corrOrdinal                  - Measures strength of association between two ordered categorical variables                                                          - MULT-Categorical    - 2021 Feb 01
%   corrpdf                      - Correlation coefficient probability density function                                                                                - UTISTAT             - 2022 Oct 31
%   covplot                      - Plots the trajectories of the elements of the covariance (correlation) matrix monitored                                             - VIS-Mult            - 2022 Apr 24
%   CressieRead                  - Computes the power divergence family                                                                                                - MULT-Categorical    - 2021 Feb 01
%   crosstab2datamatrix          - Recreates the original data matrix X from contingency table N                                                                       - MULT-Categorical    - 2022 Mar 29
%   ctlcurves                    - Computes Classification Trimmed Likelihood Curves                                                                                   - CLUS-RobClaMULT     - 2022 Jun 22
%   ctlcurvesplot                - Plots the output of routine ctlcurves                                                                                               - VIS-Clu             - 2022 Apr 24
%   ctsub                        - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                          - UTISTAT             - 2021 Oct 11
%   dempk                        - Performs a merging of components found by tkmeans                                                                                   - CLUS-RobClaMULT     - 2022 Jun 23
%   distribspec                  - Plots a probability density function between specification limits                                                                   - UTISTAT             - 2022 Oct 18
%   ellipse                      - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                          - UTISTAT             - 2022 Aug 31
%   estregimeTAR                 - Estimate a regression model with OLS in one of the regimes of a TAR model                                                           - REG-Regression      - 2021 Feb 01
%   exactcdf                     - Finds exact p-values                                                                                                                - UTISTAT             - 2021 Feb 01
%   existFS                      - Check if file exists and puts answer in a cached persistent variable                                                                - UTIGEN              - 2021 Feb 01
%   fanBIC                       - Uses the output of FSRfan to choose the best value of the transformation parameter in linear regression                             - VIS-Reg             - 2021 Oct 16
%   fanBICpn                     - Uses the output of FSRfan called with input option family 'YJpn' to choose la_P and la_N                                            - VIS-Reg             - 2021 Oct 16
%   fanplot                      - Plots the fan plot for transformation in linear regression                                                                          - VIS-Reg             - 2022 Apr 24
%   FE_int_vol                   - Computes the integrated variance from a diffusion process via the Fourier estimator using Dirichlet kernel                          - UTISTAT             - 2021 Mar 01
%   FE_int_vol_Fejer             - Computes the integrated variance from a diffusion process via the Fourier estimator using Fejer kernel                              - UTISTAT             - 2021 Mar 01
%   FE_spot_vol                  - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel                          - UTISTAT             - 2021 Mar 24
%   FE_spot_vol_FFT              - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel, using the FFT algorithm - UTISTAT             - 2021 Mar 24
%   findDir                      - Finds recursively all directories in root                                                                                           - UTIGEN              - 2021 Feb 01
%   findFile                     - Finds recursively all files in root                                                                                                 - UTIGEN              - 2021 May 31
%   forecastTS                   - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                                   - REG-Regression      - 2022 Oct 31
%   FowlkesMallowsIndex          - Computes the Fowlkes and Mallows index                                                                                              - UTISTAT             - 2021 Feb 01
%   FSCorAnaeda                  - Performs forward search in correspondence analysis with exploratory data analysis purposes                                          - MULT-Multivariate   - 2021 Aug 30
%   FSCorAnaenvmmd               - Computes the empirical envelopes of Minimum MD outside subset during the search                                                     - MULT-Multivariate   - 2021 Feb 01
%   FSM                          - Computes forward search estimator in multivariate analysis                                                                          - MULT-Multivariate   - 2022 Jul 11
%   FSMbonfbound                 - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                          - UTISTAT             - 2021 Mar 24
%   FSMbsb                       - Gives the units belonging to subset at step(s) msel of the forward search                                                           - MULT-Multivariate   - 2021 Mar 18
%   FSMeda                       - Performs forward search in multivariate analysis with exploratory data analysis purposes                                            - MULT-Multivariate   - 2021 Feb 01
%   FSMedaeasy                   - Is exactly equal to FSMeda but it is much less efficient                                                                            - MULT-Multivariate   - 2021 Feb 01
%   FSMenvmmd                    - Computes the theoretical envelopes of Minimum MD outside subset during the search                                                   - MULT-Multivariate   - 2021 Jul 17
%   FSMfan                       - Computes confirmatory lrt of a suggested transformation                                                                             - MULT-Transformations- 2021 Feb 01
%   FSMinvmmd                    - Converts values of minimum Mahalanobis distance into confidence levels                                                              - MULT-Multivariate   - 2021 Dec 09
%   FSMmmd                       - Monitors minMD                                                                                                                      - MULT-Multivariate   - 2021 Oct 16
%   FSMmmdeasy                   - Is exactly equal to FSMmmd but it is much less efficient                                                                            - MULT-Multivariate   - 2021 Feb 01
%   FSMmmdrs                     - Performs random start monitoring of minimum Mahalanobis distance                                                                    - CLUS-RobClaMULT     - 2021 Nov 30
%   FSMtra                       - Computes MLE of transformation parameters                                                                                           - MULT-Transformations- 2021 Feb 01
%   FSR                          - Computes forward search estimator in linear regression                                                                              - REG-Regression      - 2022 Oct 15
%   FSRaddt                      - Produces t deletion tests for each explanatory variable                                                                             - REG-ModelSelection  - 2022 May 30
%   FSRB                         - Gives an automatic outlier detection procedure in Bayesian linear regression                                                        - REG-Bayes           - 2021 Feb 22
%   FSRBbsb                      - Returns the units belonging to the subset in each step of the Bayesian forward search                                               - REG-Regression      - 2021 Apr 29
%   FSRBeda                      - Enables to monitor several quantities in each step of the Bayesian search                                                           - REG-Bayes           - 2021 Feb 22
%   FSRBmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search                 - REG-Bayes           - 2021 Feb 22
%   FSRbonfbound                 - Computes Bonferroni bounds for each step of the search (in linear regression)                                                       - UTISTAT             - 2021 Mar 15
%   FSRBr                        - Bayesian forward search in linear regression reweighted                                                                             - REG-Bayes           - 2021 Feb 22
%   FSRbsb                       - Returns the units belonging to the subset in each step of the forward search                                                        - REG-Regression      - 2021 Jul 16
%   FSRcp                        - Monitors Cp and AIC for all models of interest of size smallp                                                                       - REG-ModelSelection  - 2021 Feb 22
%   FSReda                       - Enables to monitor several quantities in each step of the forward search                                                            - REG-Regression      - 2022 Mar 23
%   FSRenvmdr                    - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                                    - REG-Regression      - 2021 Jul 19
%   FSRfan                       - Monitors the values of the score test statistic for each lambda                                                                     - REG-Transformations - 2021 Oct 16
%   FSRH                         - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                                 - REG-Hetero          - 2021 Feb 22
%   FSRHbsb                      - Returns the units belonging to the subset in each step of the heteroskedastic forward search                                        - REG-Hetero          - 2021 Feb 22
%   FSRHeda                      - Enables to monitor several quantities in each step of the forward search                                                            - REG-Hetero          - 2021 Feb 22
%   FSRHmdr                      - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search          - REG-Hetero          - 2021 Feb 22
%   FSRinvmdr                    - Converts values of minimum deletion residual into confidence levels                                                                 - REG-Regression      - 2021 Mar 28
%   FSRmdr                       - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                          - REG-Regression      - 2021 Aug 12
%   FSRmdrrs                     - Performs random start monitoring of minimum deletion residual                                                                       - CLUS-RobClaREG      - 2021 Feb 01
%   FSRms                        - Performs robust model selection using flexible trimming in linear regression                                                        - REG-ModelSelection  - 2021 Feb 22
%   FSRr                         - Forward search in linear regression reweighted                                                                                      - REG-Regression      - 2021 Feb 22
%   FSRts                        - Is an automatic adaptive procedure to detect outliers in time series                                                                - REG-Regression      - 2022 May 30
%   FSRtsbsb                     - Returns the units belonging to the subset in each step of the forward search                                                        - REG-Regression      - 2021 Feb 22
%   FSRtsmdr                     - Computes minimum deletion residual for time series models in each step of the search                                                - REG-Regression      - 2021 Feb 22
%   funnelchart                  - Displays a funnel chart                                                                                                             - VIS-Mult            - 2022 Aug 23
%   genSigmaGPCM                 - Generates covariance matrix for the 14 Gaussian Parsimonious Clustering Models                                                      - CLUS-RobClaMULT     - 2021 Jun 10
%   GowerIndex                   - Computes matrix of similarity indexes using Gower metric                                                                            - CLUS-RobClaMULT     - 2021 Feb 01
%   GUIconcentration             - Shows the necessary calculations to obtain the GINI concentration index in a GUI                                                    - GUI                 - 2022 Apr 29
%   GUIcorr                      - Shows the necessary calculations to obtain the correlation in a GUI                                                                 - GUI                 - 2022 Nov 08
%   GUIcov                       - Shows the necessary calculations to obtain the covariance in a GUI                                                                  - GUI                 - 2022 Oct 09
%   GUImad                       - Shows the necessary calculations to obtain MAD, S_M or S_Me in a GUI                                                                - GUI                 - 2022 Jun 19
%   GUIpowermean                 - Shows the necessary calculations to obtain the power (generalized) mean in a GUI                                                    - GUI                 - 2022 Apr 29
%   GUIquantile                  - Shows the necessary calculations to obtain $x_z$ quantile                                                                           - GUI                 - 2022 Aug 19
%   GUIregress                   - Shows the necessary calculations to obtain simple linear regression statistics in a GUI                                             - GUI                 - 2022 Jun 22
%   GUIskewness                  - Shows the necessary calculations to obtain the variance in a GUI                                                                    - GUI                 - 2022 Apr 29
%   GUIstd                       - Shows the necessary calculations to obtain the standard deviation in a GUI                                                          - GUI                 - 2022 Apr 29
%   GUItrimmean                  - Shows the necessary calculations to obtain the trimmed mean in a GUI                                                                - GUI                 - 2022 May 30
%   GUIvar                       - Shows the necessary calculations to obtain the variance in a GUI                                                                    - GUI                 - 2022 Apr 29
%   GYfilt                       - Computes the Gervini-Yohai univariate outlier identifier                                                                            - UTISTAT             - 2021 Jun 27
%   HAbdp                        - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT             - 2022 Aug 23
%   HAc                          - Computes breakdown point and efficiency associated with constant c                                                                  - UTISTAT             - 2021 Feb 01
%   HAeff                        - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                        - UTISTAT             - 2022 Aug 03
%   HApsi                        - Computes psi function  using Hampel proposal                                                                                        - UTISTAT             - 2022 Aug 10
%   HApsider                     - Computes derivative of psi function  using Hampel proposal                                                                          - UTISTAT             - 2021 Feb 01
%   HApsix                       - Computes psi function  using Hampel proposal times x                                                                                - UTISTAT             - 2021 Feb 01
%   HArho                        - Computes rho function  using Hampel proposal                                                                                        - UTISTAT             - 2021 Jul 28
%   HAwei                        - Computes weight function psi(u)/u using Hampel proposal                                                                             - UTISTAT             - 2021 Jul 28
%   histFS                       - Plots a histogram with the elements in each bin grouped according to a vector of labels                                             - VIS-Reg             - 2021 Feb 01
%   htmlwriteFS                  - Is an obsolete function which will be removed in future releases. Use publishFS.m instead                                           - UTIHELP             - 2021 Feb 01
%   HUc                          - Computes breakdown point and efficiency associated with constant c for Huber link                                                   - UTISTAT             - 2022 Jul 29
%   HUeff                        - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                   - UTISTAT             - 2022 Aug 23
%   HUpsi                        - Computes psi function (derivative of rho function) for Huber                                                                        - UTISTAT             - 2021 Feb 01
%   HUpsider                     - Computes derivative of psi function (second derivative of rho function) for Huber                                                   - UTISTAT             - 2022 Jul 08
%   HUpsix                       - Computes psi function (derivative of rho function) times x for Huber                                                                - UTISTAT             - 2021 Feb 01
%   HUrho                        - Computes rho function for Huber                                                                                                     - UTISTAT             - 2022 Jul 04
%   HUwei                        - Computes weight function psi(u)/u for Huber                                                                                         - UTISTAT             - 2021 Feb 01
%   HYPbdp                       - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                       - UTISTAT             - 2022 Aug 03
%   HYPc                         - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)     - UTISTAT             - 2021 Feb 01
%   HYPck                        - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                             - UTISTAT             - 2022 Aug 03
%   HYPeff                       - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                           - UTISTAT             - 2022 Aug 03
%   HYPk                         - Computes breakdown point and efficiency for hyp. tan. estimator                                                                     - UTISTAT             - 2021 Feb 01
%   HYPpsi                       - Computes psi function for hyperbolic tangent estimator                                                                              - UTISTAT             - 2021 Feb 01
%   HYPpsider                    - Computes derivative of psi function for hyperbolic tangent estimator                                                                - UTISTAT             - 2022 Jul 20
%   HYPpsix                      - Computes psi function for hyperbolic tangent estimator times x                                                                      - UTISTAT             - 2021 Feb 01
%   HYPrho                       - Computes rho function  using hyperbolic tangent estimator                                                                           - UTISTAT             - 2021 Jul 23
%   HYPwei                       - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                  - UTISTAT             - 2021 Jul 21
%   inversegamcdf                - Computes inverse-gamma cumulative distribution function                                                                             - UTISTAT             - 2021 Feb 01
%   inversegaminv                - Inversegampdf Inverse-gamma cumulative distribution function                                                                        - UTISTAT             - 2021 Feb 01
%   inversegampdf                - Computes inverse-gamma probability density function                                                                                 - UTISTAT             - 2022 Oct 14
%   isfunction                   - Checks if a function exists                                                                                                         - UTIGEN              - 2021 Feb 01
%   kdebiv                       - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                      - UTISTAT             - 2022 Aug 24
%   levfwdplot                   - Plots the trajectories of leverage along the search                                                                                 - VIS-Reg             - 2022 Apr 24
%   lexunrank                    - Gives the the $k$-combination of $n$ elements of position $N$ in the lexicographic order of all combinations                        - UTICOMB             - 2021 Mar 13
%   lga                          - Performs linear grouping analysis                                                                                                   - CLUS-RobClaREG      - 2021 Feb 01
%   logmvnpdfFS                  - Produces log of Multivariate normal probability density function (pdf)                                                              - UTISTAT             - 2022 Jul 01
%   LTSts                        - Extends LTS estimator to time series                                                                                                - REG-Regression      - 2022 Nov 07
%   LTStsVarSel                  - Does variable selection in the robust time series model LTSts                                                                       - REG-Regression      - 2022 Nov 08
%   LXS                          - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                                - REG-Regression      - 2022 Oct 08
%   mahalCorAna                  - MahalFS computes Mahalanobis distances (in squared units) for each row of matrix Y                                                  - UTISTAT             - 2021 Feb 01
%   mahalFS                      - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT             - 2021 Nov 14
%   makecontentsfileFS           - Extends Matlab function makecontentsfile                                                                                            - UTIHELP             - 2021 Feb 01
%   malfwdplot                   - Plots the trajectories of scaled Mahalanobis distances along the search                                                             - VIS-Mult            - 2022 Apr 24
%   malindexplot                 - Plots the Mahalanobis distances versus a selected variable                                                                          - VIS-Mult            - 2022 May 04
%   mcd                          - Computes Minimum Covariance Determinant                                                                                             - MULT-Multivariate   - 2022 Jun 19
%   mcdCorAna                    - Computes Minimum Covariance Determinant in correspondence analysis                                                                  - MULT-Multivariate   - 2021 Sep 01
%   mcdeda                       - Monitors Minimum Covariance Determinant for a series of values of bdp                                                               - MULT-Multivariate   - 2021 Aug 07
%   mdpattern                    - Finds and plots missing data patterns                                                                                               - UTISTAT             - 2022 Jan 12
%   mdpd                         - Computes Minimum Distance Power Divergence statistics                                                                               - UTISTAT             - 2021 Feb 01
%   mdpdR                        - Allows to apply Minimum Density Power Divergence criterion to parametric regression problems                                        - REG-Regression      - 2021 Feb 22
%   mdpdReda                     - Allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems                                     - REG-Regression      - 2021 Feb 22
%   mdrplot                      - Plots the trajectory of minimum deletion residual (mdr)                                                                             - VIS-Reg             - 2022 Jul 11
%   mdrrsplot                    - Plots the trajectory of minimum deletion residual from random starts                                                                - CLUS-RobClaREG      - 2022 Apr 24
%   medcouple                    - Computes the medcouple, a robust skewness estimator                                                                                 - UTISTAT             - 2022 Aug 03
%   MixSim                       - Generates k clusters in v dimensions with given overlap                                                                             - CLUS-MixSim         - 2021 Feb 01
%   MixSimreg                    - Generates k regression hyperplanes in p dimensions with given overlap                                                               - CLUS-MixSim         - 2021 Feb 01
%   Mlocation                    - Finds the M estimator of location in a univariate sample                                                                            - UTISTAT             - 2022 Jul 03
%   Mlocsca                      - Computes M estimator of location and scale in univariate samples                                                                    - UTISTAT             - 2022 Jul 20
%   mmdplot                      - Plots the trajectory of minimum Mahalanobis distance (mmd)                                                                          - VIS-Mult            - 2021 Feb 01
%   mmdrsplot                    - Plots the trajectories of minimum Mahalanobis distances from different starting points                                              - VIS-Mult            - 2022 Apr 24
%   MMmult                       - Computes MM estimators in multivariate analysis with auxiliary S-scale                                                              - MULT-Multivariate   - 2022 Apr 25
%   MMmultcore                   - Computes multivariate MM estimators for a selected fixed scale                                                                      - MULT-Multivariate   - 2021 Feb 01
%   MMmulteda                    - Computes MM estimators in multivariate analysis for a series of values of eff                                                       - MULT-Multivariate   - 2022 Mar 17
%   MMreg                        - Computes MM estimator of regression coefficients                                                                                    - REG-Regression      - 2022 Apr 25
%   MMregcore                    - Computes MM regression estimators for a selected fixed scale                                                                        - REG-Regression      - 2022 Jul 20
%   MMregeda                     - Computes MM estimator in linear regression for a series of values of efficiency                                                     - REG-Regression      - 2022 Mar 17
%   moonplot                     - Draws the Correspondence Analysis (CA) moonplot                                                                                     - VIS-Mult            - 2021 Sep 02
%   mreadFS                      - Enables to create a structure with InputArgs/OptArgs/OutArgs ... from .m function files (OBSOLETE FUNCTION REPLACED BY publishFS.m) - UTIHELP             - 2021 Feb 01
%   Mscale                       - Finds the M estimator of the scale                                                                                                  - UTISTAT             - 2022 Jul 02
%   mtR                          - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                            - UTISTAT             - 2021 Feb 01
%   mve                          - Computes Minimum volume ellipsoid                                                                                                   - MULT-Multivariate   - 2022 Aug 23
%   mveeda                       - Monitors Minimum volume ellipsoid for a series of values of bdp                                                                     - MULT-Multivariate   - 2021 Aug 06
%   nchoosekFS                   - Returns the Binomial coefficient or matrix containing all combinations                                                              - UTICOMB             - 2021 Feb 01
%   ncpci                        - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                    - UTISTAT             - 2021 Mar 24
%   ncx2mixtcdf                  - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                         - UTISTAT             - 2021 Mar 24
%   normBoxCox                   - Computes (normalized) Box-Cox transformation                                                                                        - UTISTAT             - 2021 Mar 24
%   normYJ                       - Computes (normalized) Yeo-Johnson transformation                                                                                    - UTISTAT             - 2021 Mar 24
%   normYJpn                     - NormYJ computes (normalized) extended Yeo-Johnson transformation                                                                    - UTISTAT             - 2021 Mar 24
%   openMatlabFileFromHTML       - Enables to put in HTML an hypertextual link to a specific MATLAB file                                                               - UTIGEN              - 2021 Feb 01
%   OPTbdp                       - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT             - 2022 Jul 20
%   OPTc                         - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                         - UTISTAT             - 2022 Jul 20
%   OPTeff                       - Finds the constant c which is associated to the requested efficiency                                                                - UTISTAT             - 2022 Aug 03
%   OptimalCuttingFrequency      - Computes the optimal cutting frequency for the Fourier estimator of integrated variance                                             - UTISTAT             - 2021 Mar 01
%   OPTpsi                       - Computes psi function (derivative of rho function) for optimal weight function                                                      - UTISTAT             - 2022 Jul 19
%   OPTpsider                    - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                 - UTISTAT             - 2022 Jul 19
%   OPTpsix                      - Computes psi function (derivative of rho function) times x                                                                          - UTISTAT             - 2022 Jul 19
%   OPTrho                       - Computes rho function for optimal weight function                                                                                   - UTISTAT             - 2022 Jul 19
%   OPTwei                       - Computes weight function psi(u)/u for optimal weight function                                                                       - UTISTAT             - 2022 Jul 19
%   overlap                      - Computes the exact overlap given the parameters of the mixture                                                                      - CLUS-MixSim         - 2021 Feb 01
%   overlapmap                   - Produces an interactive overlap map                                                                                                 - CLUS-RobClaMULT     - 2022 Jun 23
%   pcaFS                        - Performs Principal Component Analysis (PCA) on raw data                                                                             - MULT-Multivariate   - 2021 Oct 16
%   PDbdp                        - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                          - UTISTAT             - 2021 Feb 01
%   PDc                          - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator                - UTISTAT             - 2021 Feb 01
%   PDeff                        - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                     - UTISTAT             - 2021 Feb 01
%   PDpsi                        - Computes psi function (derivative of rho function) for minimum density power divergence estimator                                   - UTISTAT             - 2021 Feb 01
%   PDpsider                     - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                      - UTISTAT             - 2021 Feb 01
%   PDpsix                       - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                           - UTISTAT             - 2022 Jul 28
%   PDrho                        - Computes rho function for minimum density power divergence estimator                                                                - UTISTAT             - 2022 Jul 03
%   PDwei                        - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                               - UTISTAT             - 2021 Aug 12
%   position                     - Controls the position of the open figures                                                                                           - UTIGEN              - 2021 Feb 08
%   Powertra                     - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                             - UTISTAT             - 2022 Apr 29
%   publishBibliography          - Enables to create web page which contains the references inside the input .m files                                                  - UTIHELP             - 2021 Feb 01
%   publishFS                    - Enables to create automatic HELP FILES from structured .m function files                                                            - UTIHELP             - 2022 Jun 22
%   publishFunctionAlpha         - Enables to create web page which contains the alphabetical list of functions                                                        - UTIHELP             - 2021 Feb 01
%   publishFunctionCate          - Enables to create web page which contains the categorical list of functions                                                         - UTIHELP             - 2021 Sep 01
%   Qn                           - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                      - UTISTAT             - 2021 Mar 27
%   qqplotFS                     - Qqplot of studentized residuals with envelopes                                                                                      - VIS-Reg             - 2022 Jul 30
%   quickselectFS                - Finds the k-th order statistic                                                                                                      - UTIGEN              - 2022 Mar 29
%   quickselectFS_demo           - Illustrates the functioning of quickselectFS                                                                                        - GUI                 - 2022 Aug 03
%   quickselectFSw               - Finds the 100*p-th weighted order statistic for $0<p<1$                                                                             - UTIGEN              - 2022 May 30
%   quickselectFSw_demo          - Illustrates the functioning of quickselectFSw                                                                                       - GUI                 - 2022 Aug 03
%   RandIndexFS                  - Calculates Rand type Indices to compare two partitions                                                                              - UTISTAT             - 2021 Feb 01
%   randsampleFS                 - Generates a random sample of k elements from the integers 1 to n (k<=n)                                                             - UTICOMB             - 2021 Feb 01
%   rcontFS                      - Generates a random two-way table with given marginal totals                                                                         - MULT-Categorical    - 2021 Aug 28
%   regressB                     - Computes Bayesian estimates of regression parameters                                                                                - REG-Bayes           - 2021 Feb 22
%   regressH                     - Fits a multiple linear regression model with heteroskedasticity                                                                     - REG-Hetero          - 2021 Feb 22
%   regressHart                  - Fits a multiple linear regression model using ART heteroskedasticity                                                                - REG-Hetero          - 2021 Feb 22
%   regressHhar                  - Fits a multiple linear regression model with Harvey heteroskedasticity                                                              - REG-Hetero          - 2021 Feb 22
%   regressts                    - Computes estimates of regression parameters for a time series models                                                                - REG-Regression      - 2021 Oct 16
%   removeExtraSpacesLF          - Removes extra spaces and selected carriage returns from input string                                                                - UTIGEN              - 2021 Feb 01
%   repDupValWithMean            - Replaces values of y which have non unique elements in vector x with local means                                                    - UTIGEN              - 2021 Apr 21
%   resfwdplot                   - Plots the trajectories of the monitored scaled (squared) residuals                                                                  - VIS-Reg             - 2021 Mar 11
%   resindexplot                 - Plots the residuals from a regression analysis versus index number or any other variable                                            - VIS-Reg             - 2022 Jun 19
%   restrdeter                   - Computes determinant restriction                                                                                                    - CLUS-RobClaMULT     - 2021 Feb 01
%   restrdeterGPCM               - Applies determinat restrictions for the 14 GPCM                                                                                     - CLUS-RobClaMULT     - 2021 Feb 01
%   restreigen                   - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2021 Oct 16
%   restreigeneasy               - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2022 May 05
%   restreigenmemopt             - Computes eigenvalues restriction (without Dykstra algorithm)                                                                        - CLUS-RobClaMULT     - 2021 Oct 16
%   restrshapeExact              - Computes constrained Gamma (shape) matrix with exact constraints                                                                    - CLUS-RobClaMULT     - 2021 Oct 16
%   restrshapeGPCM               - Produces the restricted shape matrix for the 14 GPCM                                                                                - CLUS-RobClaMULT     - 2021 Feb 01
%   restrSigmaGPCM               - Computes constrained covariance matrices for the 14 GPCM specifications                                                             - CLUS-RobClaMULT     - 2021 Oct 16
%   RhoPsiWei                    - Finds rho, psi, psi', w functions given bdp, or eff or tuning constant c                                                            - UTISTAT             - 2022 Aug 19
%   RKbdp                        - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                     - UTISTAT             - 2021 Feb 01
%   RKeff                        - Finds the constants c and M which are associated to the requested efficiency and ARP                                                - UTISTAT             - 2021 Feb 01
%   RKpsi                        - Computes psi function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT             - 2021 Feb 01
%   RKpsider                     - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                     - UTISTAT             - 2021 Feb 01
%   RKpsix                       - Computes psi function times x for Rocke (translated Tukey's) biweight                                                               - UTISTAT             - 2021 Feb 01
%   RKrho                        - Computes rho function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT             - 2021 Feb 01
%   RKwei                        - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                           - UTISTAT             - 2021 Feb 01
%   rlga                         - Performs robust linear grouping analysis                                                                                            - CLUS-RobClaREG      - 2021 Feb 01
%   rlsmo                        - Computes a running-lines smoother with global cross-validation                                                                      - REG-Transformations - 2022 Jun 27
%   RobCov                       - Computes covariance matrix of robust regression coefficients                                                                        - REG-Regression      - 2022 Jul 20
%   RobRegrSize                  - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size                        - REG-Regression      - 2021 Aug 11
%   rthin                        - Applies independent random thinning to a point pattern                                                                              - UTISTAT             - 2021 Feb 01
%   scatterboxplot               - Creates scatter diagram with marginal boxplots                                                                                      - VIS-Mult            - 2021 Dec 04
%   Score                        - Computes the score test for transformation                                                                                          - REG-Transformations - 2021 Feb 28
%   ScoreYJ                      - Computes the score test for Yeo and Johnson transformation                                                                          - REG-Transformations - 2021 Feb 28
%   ScoreYJall                   - Computes all the 4 score tests for YJ transformation                                                                                - REG-Transformations - 2021 Feb 28
%   ScoreYJmle                   - Computes the likelihood ratio test fof H_0=lambdaP=lambdaP0 and lambdaN=lambdaN0                                                    - REG-Transformations - 2021 Feb 22
%   ScoreYJpn                    - Computes the score test for YJ transformation for pos and neg observations                                                          - REG-Transformations - 2021 Mar 01
%   SDest                        - Computes Stahel-Donoho robust estimator of dispersion-location                                                                      - MULT-Multivariate   - 2021 Feb 01
%   SETARX                       - Implements Threshold autoregressive models with two regimes                                                                         - REG-Regression      - 2021 Mar 02
%   setdiffFS                    - Finds the positive integers in a which are not present in the positive integers in b                                                - UTIGEN              - 2021 Aug 24
%   setToolboxStartEnd           - Sets release compatibility in ToolboxPackagingConfiguration.prj file                                                                - UTIHELP             - 2022 Feb 04
%   shuffling                    - Does a random permutation of the elements of input vector                                                                           - UTICOMB             - 2021 Feb 01
%   simdataset                   - Simulates and-or contaminates a dataset given the parameters of a finite mixture model with Gaussian components                     - CLUS-MixSim         - 2021 Feb 01
%   simdatasetreg                - Simulates a regression dataset given the parameters of a mixture regression model                                                   - CLUS-MixSim         - 2021 Feb 01
%   simulateLM                   - Simulates linear regression data with prespecified values of statistical indexes                                                    - REG-Regression      - 2022 Apr 04
%   simulateTS                   - Simulates a time series with trend, time varying seasonal, level shift and irregular component                                      - REG-Regression      - 2022 Oct 23
%   smothr                       - Produces smoothed values with constraints                                                                                           - REG-Transformations - 2021 Apr 13
%   Smult                        - Computes S estimators in multivariate analysis                                                                                      - MULT-Multivariate   - 2021 Feb 01
%   Smulteda                     - Computes S estimators in multivariate analysis for a series of values of bdp                                                        - MULT-Multivariate   - 2021 Feb 01
%   Sn                           - Robust estimator of scale (robust version of Gini's average difference)                                                             - UTISTAT             - 2021 Feb 01
%   SparseTableTest              - Computes independence test for large and sparse contingency tables                                                                  - MULT-Categorical    - 2021 Feb 01
%   spmplot                      - Produces an interactive scatterplot matrix with boxplots or histograms on the main diagonal and possibly robust bivariate contours  - VIS-Mult            - 2022 Aug 31
%   Sreg                         - Computes S estimators in linear regression                                                                                          - REG-Regression      - 2022 Oct 31
%   Sregeda                      - Computes S estimators in linear regression for a series of values of bdp                                                            - REG-Regression      - 2022 Aug 03
%   subsets                      - Creates a matrix of indexes where rows are distinct p-subsets extracted from a set of n elements                                    - UTICOMB             - 2021 Aug 11
%   suplabel                     - Places text as a title, xlabel, or ylabel on a group of subplots                                                                    - UTIGEN              - 2021 Feb 01
%   supsmu                       - Smooths scatterplots using Friedman's supersmoother algorithm                                                                       - REG-Transformations - 2021 Mar 01
%   tabulateFS                   - Creates frequency table of unique values of x, excluding possible 0 counts                                                          - UTISTAT             - 2021 Feb 01
%   Taureg                       - Computes Tau estimators in linear regression                                                                                        - REG-Regression      - 2022 Aug 03
%   TBbdp                        - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                - UTISTAT             - 2021 Jul 20
%   TBc                          - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                             - UTISTAT             - 2021 Feb 01
%   TBeff                        - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                           - UTISTAT             - 2021 Jul 23
%   tBothSides                   - Allows users to transform both sides of a (nonlinear) regression model                                                              - REG-Regression      - 2021 Feb 22
%   TBpsi                        - Computes psi function (derivative of rho function) for Tukey's biweight                                                             - UTISTAT             - 2021 Feb 01
%   TBpsider                     - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                        - UTISTAT             - 2021 Feb 01
%   TBpsix                       - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                     - UTISTAT             - 2021 Feb 01
%   TBrho                        - Computes rho function for Tukey's biweight                                                                                          - UTISTAT             - 2022 Jul 28
%   TBwei                        - Computes weight function psi(u)/u for Tukey's biweight                                                                              - UTISTAT             - 2021 Dec 04
%   tclust                       - Computes trimmed clustering with scatter restrictions                                                                               - CLUS-RobClaMULT     - 2022 Jul 03
%   tclusteda                    - Computes tclust for a series of values of the trimming factor                                                                       - CLUS-RobClaMULT     - 2021 Oct 16
%   tclustIC                     - Computes tclust for different number of groups k and restriction factors c                                                          - CLUS-RobClaMULT     - 2022 Sep 22
%   tclustICgpcm                 - Computes tclust for different number of groups k and restr. factors $c_{det}$ and $c_{shw}$                                         - CLUS-RobClaMULT     - 2022 Mar 10
%   tclustICplot                 - Plots information criterion as a function of c and k                                                                                - VIS-Clu             - 2022 Apr 24
%   tclustICplotGPCM             - Plots information criterion as a function of  $c_{det}$, $c_{shw}$,  $c_{shb}$ and $k$                                              - VIS-Clu             - 2022 Apr 24
%   tclustICsol                  - Extracts a set of best relevant solutions                                                                                           - CLUS-RobClaMULT     - 2021 Feb 01
%   tclustICsolGPCM              - Extracts a set of best relevant solutions from 3D array computed using function tclustICgpcm                                        - CLUS-RobClaMULT     - 2021 Feb 01
%   tclustreg                    - Performs robust linear grouping analysis                                                                                            - CLUS-RobClaREG      - 2022 Aug 23
%   tclustregeda                 - Performs robust linear grouping analysis for a series of values of the trimming factor                                              - CLUS-RobClaREG      - 2022 May 09
%   tclustregIC                  - Computes tclustreg for different number of groups k and restriction factors c                                                       - CLUS-RobClaREG      - 2022 Jun 23
%   tkmeans                      - Computes trimmed k-means                                                                                                            - CLUS-RobClaMULT     - 2022 Jun 23
%   triu2vec                     - Extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix                           - UTIGEN              - 2021 Feb 01
%   twdcdf                       - Computes the cumulative distribution function of the Tweedie distribution                                                           - UTISTAT             - 2021 Dec 04
%   twdpdf                       - Twopdf computes the probability density function of the Tweedie distribution                                                        - UTISTAT             - 2021 Mar 01
%   twdrnd                       - Generates random variates from the Tweedie distribution                                                                             - UTISTAT             - 2021 Feb 01
%   unibiv                       - Has the purpose of detecting univariate and bivariate outliers                                                                      - MULT-Multivariate   - 2021 Jul 26
%   upperfracpos                 - Positions two figures on the upper part of the screen                                                                               - UTIGEN              - 2021 Feb 01
%   verLessThanFS                - Compares version of MATLAB to specified version number                                                                              - UTIGEN              - 2022 Mar 17
%   vervaatrnd                   - Simulates random variates from the Vervaat perpetuity distribution                                                                  - UTISTAT             - 2021 Feb 01
%   vervaatsim                   - Returns a Vervaat perpetuity                                                                                                        - UTISTAT             - 2021 Feb 01
%   vervaatxdf                   - Returns the pdf and cdf of a Vervaat perpetuity                                                                                     - UTISTAT             - 2021 Feb 01
%   VIOM                         - Computes weights estimates under Variance-Inflation Model                                                                           - REG-Regression      - 2021 Feb 22
%   waterfallchart               - Creates a waterfall chart                                                                                                           - VIS-Mult            - 2022 Oct 18
%   wedgeplot                    - Generates the double wedge plot of a time series                                                                                    - VIS-Reg             - 2022 May 09
%   winsor                       - Returns a winsorized copy of input                                                                                                  - UTISTAT             - 2021 Feb 01
%   WNChygepdf                   - Returns Wallenius' non-central hypergeometric probability density values                                                            - UTISTAT             - 2021 Feb 01
%   wraptextFS                   - Formats long strings into wrapped text of specified width                                                                           - UTIGEN              - 2021 Feb 01
%   wthin                        - Thins a uni/bi-dimensional dataset                                                                                                  - UTISTAT             - 2021 Mar 24
%   xmlcreateFS                  - Creates an XML file passing through publishFS                                                                                       - UTIHELP             - 2021 Feb 01
%   yXplot                       - Produces an interactive scatterplot of y against each variable of X in the input dataset                                            - VIS-Reg             - 2022 Apr 24
%   zscoreFS                     - Computes (robust) standardized z scores                                                                                             - UTIGEN              - 2021 May 11
