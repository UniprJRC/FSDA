% REGRESSION
%
% File names, description, category and date last modified
%
%   Name                 - Description                                                                                                                - Category           - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   ace                  - Computes alternative conditional expectation                                                                               - REG-Transformations- 2022 Nov 10
%   addt                 - Produces the t test for an additional explanatory variable                                                                 - REG-Regression     - 2022 Nov 10
%   avas                 - Computes additivity and variance stabilization for regression                                                              - REG-Transformations- 2022 Nov 10
%   avasms               - Computes avas using a series of alternative options                                                                        - REG-Transformations- 2022 Nov 10
%   boxcoxR              - Finds MLE of lambda in linear regression (and confidence interval) using Box Cox, YJ or extended YJ  transformation        - REG-Transformations- 2022 Nov 10
%   estregimeTAR         - Estimate a regression model with OLS in one of the regimes of a TAR model                                                  - REG-Regression     - 2022 Nov 10
%   fanBIC               - Uses the output of FSRfan to choose the best value of the transformation parameter in linear regression                    - VIS-Reg            - 2022 Dec 02
%   fanBICpn             - Uses the output of FSRfan called with input option family 'YJpn' to choose la_P and la_N                                   - VIS-Reg            - 2022 Nov 10
%   forecastTS           - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                          - REG-Regression     - 2022 Nov 14
%   FSR                  - Computes forward search estimator in linear regression                                                                     - REG-Regression     - 2022 Nov 10
%   FSRaddt              - Produces t deletion tests for each explanatory variable                                                                    - REG-ModelSelection - 2022 Nov 10
%   FSRB                 - Gives an automatic outlier detection procedure in Bayesian linear regression                                               - REG-Bayes          - 2022 Nov 10
%   FSRBbsb              - Returns the units belonging to the subset in each step of the Bayesian forward search                                      - REG-Regression     - 2022 Nov 10
%   FSRBeda              - Enables to monitor several quantities in each step of the Bayesian search                                                  - REG-Bayes          - 2022 Nov 10
%   FSRBmdr              - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search        - REG-Bayes          - 2022 Nov 10
%   FSRBr                - Bayesian forward search in linear regression reweighted                                                                    - REG-Bayes          - 2022 Nov 10
%   FSRbsb               - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2022 Nov 10
%   FSRcp                - Monitors Cp and AIC for all models of interest of size smallp                                                              - REG-ModelSelection - 2022 Nov 10
%   FSReda               - Enables to monitor several quantities in each step of the forward search                                                   - REG-Regression     - 2022 Nov 10
%   FSRenvmdr            - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                           - REG-Regression     - 2022 Nov 10
%   FSRfan               - Monitors the values of the score test statistic for each lambda                                                            - REG-Transformations- 2022 Nov 10
%   FSRH                 - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                        - REG-Hetero         - 2022 Nov 10
%   FSRHbsb              - Returns the units belonging to the subset in each step of the heteroskedastic forward search                               - REG-Hetero         - 2022 Nov 10
%   FSRHeda              - Enables to monitor several quantities in each step of the forward search                                                   - REG-Hetero         - 2022 Nov 10
%   FSRHmdr              - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search - REG-Hetero         - 2022 Nov 10
%   FSRinvmdr            - Converts values of minimum deletion residual into confidence levels                                                        - REG-Regression     - 2022 Nov 10
%   FSRmdr               - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                 - REG-Regression     - 2022 Nov 10
%   FSRms                - Performs robust model selection using flexible trimming in linear regression                                               - REG-ModelSelection - 2022 Nov 10
%   FSRr                 - Forward search in linear regression reweighted                                                                             - REG-Regression     - 2022 Nov 10
%   FSRts                - Is an automatic adaptive procedure to detect outliers in time series                                                       - REG-Regression     - 2022 Nov 10
%   FSRtsbsb             - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2022 Nov 10
%   FSRtsmdr             - Computes minimum deletion residual for time series models in each step of the search                                       - REG-Regression     - 2022 Nov 10
%   GYfilt               - Computes the Gervini-Yohai univariate outlier identifier                                                                   - UTISTAT            - 2022 Dec 02
%   LTSts                - Extends LTS estimator to time series                                                                                       - REG-Regression     - 2022 Dec 02
%   LTStsLSmult          - Extends LTSts to the detection of multiple Level Shifts in time series                                                     - REG-Regression     - 2022 Nov 28
%   LTStsVarSel          - Does variable selection in the robust time series model LTSts                                                              - REG-Regression     - 2022 Nov 21
%   LXS                  - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                       - REG-Regression     - 2022 Nov 21
%   mdpdR                - Allows to apply Minimum Density Power Divergence criterion to parametric regression problems                               - REG-Regression     - 2022 Nov 10
%   mdpdReda             - Allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems                            - REG-Regression     - 2022 Nov 10
%   MMreg                - Computes MM estimator of regression coefficients                                                                           - REG-Regression     - 2022 Nov 10
%   MMregcore            - Computes MM regression estimators for a selected fixed scale                                                               - REG-Regression     - 2022 Nov 10
%   MMregeda             - Computes MM estimator in linear regression for a series of values of efficiency                                            - REG-Regression     - 2022 Nov 10
%   regressB             - Computes Bayesian estimates of regression parameters                                                                       - REG-Bayes          - 2022 Nov 10
%   regressH             - Fits a multiple linear regression model with heteroskedasticity                                                            - REG-Hetero         - 2022 Nov 10
%   regressHart          - Fits a multiple linear regression model using ART heteroskedasticity                                                       - REG-Hetero         - 2022 Nov 10
%   regressHhar          - Fits a multiple linear regression model with Harvey heteroskedasticity                                                     - REG-Hetero         - 2022 Nov 10
%   regressts            - Computes estimates of regression parameters for a time series models                                                       - REG-Regression     - 2022 Nov 10
%   rlsmo                - Computes a running-lines smoother with global cross-validation                                                             - REG-Transformations- 2022 Nov 10
%   RobCov               - Computes covariance matrix of robust regression coefficients                                                               - REG-Regression     - 2022 Nov 10
%   RobRegrSize          - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size               - REG-Regression     - 2022 Nov 10
%   Score                - Computes the score test for transformation                                                                                 - REG-Transformations- 2022 Nov 10
%   ScoreYJ              - Computes the score test for Yeo and Johnson transformation                                                                 - REG-Transformations- 2022 Nov 10
%   ScoreYJall           - Computes all the 4 score tests for YJ transformation                                                                       - REG-Transformations- 2022 Nov 10
%   ScoreYJmle           - Computes the likelihood ratio test fof H_0=lambdaP=lambdaP0 and lambdaN=lambdaN0                                           - REG-Transformations- 2022 Nov 10
%   ScoreYJpn            - Computes the score test for YJ transformation for pos and neg observations                                                 - REG-Transformations- 2022 Nov 10
%   SETARX               - Implements Threshold autoregressive models with two regimes                                                                - REG-Regression     - 2022 Nov 10
%   simulateLM           - Simulates linear regression data with prespecified values of statistical indexes                                           - REG-Regression     - 2022 Nov 10
%   simulateTS           - Simulates a time series with trend, time varying seasonal, level shift and irregular component                             - REG-Regression     - 2022 Nov 10
%   smothr               - Produces smoothed values with constraints                                                                                  - REG-Transformations- 2022 Nov 10
%   Sreg                 - Computes S estimators in linear regression                                                                                 - REG-Regression     - 2022 Nov 10
%   Sregeda              - Computes S estimators in linear regression for a series of values of bdp                                                   - REG-Regression     - 2022 Nov 10
%   supsmu               - Smooths scatterplots using Friedman's supersmoother algorithm                                                              - REG-Transformations- 2022 Nov 10
%   Taureg               - Computes Tau estimators in linear regression                                                                               - REG-Regression     - 2022 Nov 10
%   tBothSides           - Allows users to transform both sides of a (nonlinear) regression model                                                     - REG-Regression     - 2022 Nov 10
%   VIOM                 - Computes weights estimates under Variance-Inflation Model                                                                  - REG-Regression     - 2022 Nov 10
