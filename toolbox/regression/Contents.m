% REGRESSION
%
% File names, description, category and date last modified
%
%   Name                 - Description                                                                                                                - Category           - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   ace                  - Computes alternative conditional expectation                                                                               - REG-Transformations- 2025 May 30
%   addt                 - Produces the t test for an additional explanatory variable                                                                 - REG-Regression     - 2025 May 30
%   avas                 - Computes additivity and variance stabilization for regression                                                              - REG-Transformations- 2025 May 30
%   avasms               - Computes avas using a series of alternative options                                                                        - REG-Transformations- 2025 May 30
%   boxcoxR              - Finds MLE of lambda in linear regression (and confidence interval) using Box Cox, YJ or extended YJ  transformation        - REG-Transformations- 2025 May 30
%   estregimeTAR         - Estimate a regression model with OLS in one of the regimes of a TAR model                                                  - REG-Regression     - 2025 May 30
%   fanBIC               - Uses the output of FSRfan to choose the best value of the transformation parameter in linear regression                    - VIS-Reg            - 2025 May 30
%   fanBICpn             - Uses the output of FSRfan called with input option family 'YJpn' to choose la_P and la_N                                   - VIS-Reg            - 2025 May 30
%   forecastH            - Produce forecasts with confidence bands for regression model with heteroskedasticity                                       - REG-Hetero         - 2025 May 30
%   forecastTS           - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                          - REG-Regression     - 2025 Aug 09
%   FSR                  - Computes forward search estimator in linear regression                                                                     - REG-Regression     - 2025 Oct 14
%   FSRaddt              - Produces t deletion tests for each explanatory variable                                                                    - REG-ModelSelection - 2025 May 30
%   FSRB                 - Gives an automatic outlier detection procedure in Bayesian linear regression                                               - REG-Bayes          - 2025 May 30
%   FSRBbsb              - Returns the units belonging to the subset in each step of the Bayesian forward search                                      - REG-Regression     - 2025 May 30
%   FSRBeda              - Enables to monitor several quantities in each step of the Bayesian search                                                  - REG-Bayes          - 2025 May 30
%   FSRBmdr              - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search        - REG-Bayes          - 2025 May 30
%   FSRBr                - Bayesian forward search in linear regression reweighted                                                                    - REG-Bayes          - 2025 May 30
%   FSRbsb               - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2025 May 30
%   FSRcp                - Monitors Cp and AIC for all models of interest of size smallp                                                              - REG-ModelSelection - 2025 Sep 14
%   FSReda               - Enables to monitor several quantities in each step of the forward search                                                   - REG-Regression     - 2025 May 30
%   FSRedaCens           - Enables to monitor several quantities in each step of the forward search                                                   - REG-Regression     - 2025 May 30
%   FSRenvmdr            - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                           - REG-Regression     - 2025 May 30
%   FSRfan               - Monitors the values of the score test statistic for each lambda                                                            - REG-Transformations- 2025 May 30
%   FSRfanCens           - Monitors the values of the signed sqrt LR test for each lambda in the tobit model                                          - REG-Transformations- 2025 May 30
%   FSRH                 - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                        - REG-Hetero         - 2025 May 30
%   FSRHbsb              - Returns the units belonging to the subset in each step of the heteroskedastic forward search                               - REG-Hetero         - 2025 May 30
%   FSRHeda              - Enables to monitor several quantities in each step of the forward search                                                   - REG-Hetero         - 2025 May 30
%   FSRHmdr              - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search - REG-Hetero         - 2025 May 30
%   FSRinvmdr            - Converts values of mdr into confidence levels and mdr in normal coordinates                                                - REG-Regression     - 2025 May 30
%   FSRmdr               - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                 - REG-Regression     - 2025 May 30
%   FSRms                - Performs robust model selection using flexible trimming in linear regression                                               - REG-ModelSelection - 2025 May 30
%   FSRr                 - Forward search in linear regression reweighted                                                                             - REG-Regression     - 2025 May 30
%   FSRts                - Is an automatic adaptive procedure to detect outliers in time series                                                       - REG-Regression     - 2025 May 30
%   FSRtsbsb             - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2025 May 30
%   FSRtsmdr             - Computes minimum deletion residual for time series models in each step of the search                                       - REG-Regression     - 2025 May 30
%   GYfilt               - Computes the Gervini-Yohai univariate outlier identifier                                                                   - UTISTAT            - 2025 May 30
%   LTSts                - Extends LTS estimator to time series                                                                                       - REG-Regression     - 2025 May 30
%   LTStsLSmult          - Extends LTSts to the detection of multiple Level Shifts in time series                                                     - REG-Regression     - 2025 May 30
%   LTStsVarSel          - Does variable selection in the robust time series model LTSts                                                              - REG-Regression     - 2025 May 30
%   LXS                  - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                       - REG-Regression     - 2025 May 30
%   mdpdR                - Allows to apply Minimum Density Power Divergence criterion to parametric regression problems                               - REG-Regression     - 2025 May 30
%   mdpdReda             - Allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems                            - REG-Regression     - 2025 May 30
%   MMreg                - Computes MM estimator of regression coefficients                                                                           - REG-Regression     - 2025 May 30
%   MMregcore            - Computes MM regression estimators for a selected fixed scale                                                               - REG-Regression     - 2025 May 30
%   MMregeda             - Computes MM estimator in linear regression for a series of values of efficiency                                            - REG-Regression     - 2025 May 30
%   regressB             - Computes Bayesian estimates of regression parameters                                                                       - REG-Bayes          - 2025 May 30
%   regressCens          - Computes estimates of regression parameters under the censored (Tobit) model                                               - REG-Regression     - 2025 May 30
%   regressCensTra       - Computes signed sqrt LR test for lambda in the censored (Tobit) model                                                      - REG-Regression     - 2025 May 30
%   regressH             - Fits a multiple linear regression model with heteroskedasticity                                                            - REG-Hetero         - 2025 May 30
%   regressHart          - Fits a multiple linear regression model using ART heteroskedasticity                                                       - REG-Hetero         - 2025 May 30
%   regressHhar          - Fits a multiple linear regression model with Harvey heteroskedasticity                                                     - REG-Hetero         - 2025 May 30
%   regressts            - Computes estimates of regression parameters for a time series models                                                       - REG-Regression     - 2025 May 30
%   rlsmo                - Computes a running-lines smoother with global cross-validation                                                             - REG-Transformations- 2025 May 30
%   RobCov               - Computes covariance matrix of robust regression coefficients                                                               - REG-Regression     - 2025 May 30
%   RobRegrSize          - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size               - REG-Regression     - 2025 May 30
%   Score                - Computes the score test for Box-Cox transformation                                                                         - REG-Transformations- 2025 May 30
%   ScoreT               - Computes Tukey's one df test for non-additivity on Box-Cox transformed data                                                - REG-Transformations- 2025 May 30
%   ScoreYJ              - Computes the score test for Yeo and Johnson transformation                                                                 - REG-Transformations- 2025 May 30
%   ScoreYJall           - Computes all the 4 score tests for YJ transformation                                                                       - REG-Transformations- 2025 May 30
%   ScoreYJmle           - Computes the likelihood ratio test for H_0=lambdaP=lambdaP0 and lambdaN=lambdaN0                                           - REG-Transformations- 2025 May 30
%   ScoreYJpn            - Computes the score test for YJ transformation for pos and neg observations                                                 - REG-Transformations- 2025 May 30
%   SETARX               - Implements Threshold autoregressive models with two regimes                                                                - REG-Regression     - 2025 May 30
%   simulateLM           - Simulates linear regression data with pre-specified values of statistical indexes                                          - REG-Regression     - 2025 May 30
%   simulateTS           - Simulates a time series with trend, time varying seasonal, level shift and irregular component                             - REG-Regression     - 2025 May 30
%   smothr               - Produces smoothed values with constraints                                                                                  - REG-Transformations- 2025 May 30
%   Sreg                 - Computes S estimators in linear regression                                                                                 - REG-Regression     - 2025 Dec 02
%   Sregeda              - Computes S estimators in linear regression for a series of values of bdp                                                   - REG-Regression     - 2025 May 30
%   supsmu               - Smooths scatterplots using Friedman's supersmoother algorithm                                                              - REG-Transformations- 2025 May 30
%   Taureg               - Computes Tau estimators in linear regression                                                                               - REG-Regression     - 2025 May 30
%   tBothSides           - Allows users to transform both sides of a (nonlinear) regression model                                                     - REG-Regression     - 2025 May 30
%   univariatems         - Performs preliminary univariate robust model selection in linear regression                                                - REG-Regression     - 2025 May 30
%   VIOM                 - Computes weights estimates under Variance-Inflation Model                                                                  - REG-Regression     - 2025 May 30
