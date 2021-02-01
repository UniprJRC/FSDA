% REGRESSION
%
% File names, description, category and date last modified
%
%   Name             - Description                                                                                                                - Category           - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   ace              - Computes alternative conditional expectation                                                                               - REG-Transformations- 2021 Feb 01
%   addt             - Produces the t test for an additional explanatory variable                                                                 - REG-Regression     - 2021 Feb 01
%   avas             - Computes additivity and variance stabilization for regression                                                              - REG-Transformations- 2021 Feb 01
%   boxcoxR          - Finds MLE of lambda in linear regression (and confidence interval) using Box Cox, YJ or extended YJ  transformation        - REG-Transformations- 2021 Feb 01
%   estregimeTAR     - Estimate a regression model with OLS in one of the regimes of a TAR model                                                  - REG-Regression     - 2021 Feb 01
%   fanBIC           - Uses the output of FSRfan to choose the best value of the transformation parameter in linear regression                    - VIS-Reg            - 2021 Feb 01
%   fanBICpn         - Uses the output of FSRfan called with input option family 'YJpn' to choose la_P and la_N                                   - VIS-Reg            - 2021 Feb 01
%   forecastTS       - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                          - REG-Regression     - 2021 Feb 01
%   FSR              - Computes forward search estimator in linear regression                                                                     - REG-Regression     - 2021 Feb 01
%   FSRaddt          - Produces t deletion tests for each explanatory variable                                                                    - REG-ModelSelection - 2021 Feb 01
%   FSRB             - Gives an automatic outlier detection procedure in Bayesian linear regression                                               - REG-Bayes          - 2021 Feb 01
%   FSRBbsb          - Returns the units belonging to the subset in each step of the Bayesian forward search                                      - REG-Regression     - 2021 Feb 01
%   FSRBeda          - Enables to monitor several quantities in each step of the Bayesian search                                                  - REG-Bayes          - 2021 Feb 01
%   FSRBmdr          - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search        - REG-Bayes          - 2021 Feb 01
%   FSRBr            - Bayesian forward search in linear regression reweighted                                                                    - REG-Bayes          - 2021 Feb 01
%   FSRbsb           - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2021 Feb 01
%   FSRcp            - Monitors Cp and AIC for all models of interest of size smallp                                                              - REG-ModelSelection - 2021 Feb 01
%   FSReda           - Enables to monitor several quantities in each step of the forward search                                                   - REG-Regression     - 2021 Feb 01
%   FSRenvmdr        - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                           - REG-Regression     - 2021 Feb 01
%   FSRfan           - Monitors the values of the score test statistic for each lambda                                                            - REG-Transformations- 2021 Feb 01
%   FSRH             - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                        - REG-Hetero         - 2021 Feb 01
%   FSRHbsb          - Returns the units belonging to the subset in each step of the heteroskedastic forward search                               - REG-Hetero         - 2021 Feb 01
%   FSRHeda          - Enables to monitor several quantities in each step of the forward search                                                   - REG-Hetero         - 2021 Feb 01
%   FSRHmdr          - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search - REG-Hetero         - 2021 Feb 01
%   FSRinvmdr        - Converts values of minimum deletion residual into confidence levels                                                        - REG-Regression     - 2021 Feb 01
%   FSRmdr           - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                 - REG-Regression     - 2021 Feb 01
%   FSRms            - Performs robust model selection using flexible trimming in linear regression                                               - REG-ModelSelection - 2021 Feb 01
%   FSRr             - Forward search in linear regression reweighted                                                                             - REG-Regression     - 2021 Feb 01
%   FSRts            - Is an automatic adaptive procedure to detect outliers in time series                                                       - REG-Regression     - 2021 Feb 01
%   FSRtsbsb         - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2021 Feb 01
%   FSRtsmdr         - Computes minimum deletion residual for time series models in each step of the search                                       - REG-Regression     - 2021 Feb 01
%   GYfilt           - Computes the Gervini-Yohai univariate outlier identifier                                                                   - UTISTAT            - 2021 Feb 01
%   LTSts            - Extends LTS estimator to time series                                                                                       - REG-Regression     - 2021 Feb 01
%   LTStsVarSel      - Does variable selection in the robust time series model LTSts                                                              - REG-Regression     - 2021 Feb 01
%   LXS              - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                       - REG-Regression     - 2021 Feb 01
%   mdpdR            - Allows to apply Minimum Density Power Divergence criterion to parametric regression problems                               - REG-Regression     - 2021 Feb 01
%   mdpdReda         - Allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems                            - REG-Regression     - 2021 Feb 01
%   MMreg            - Computes MM estimator of regression coefficients                                                                           - REG-Regression     - 2021 Feb 01
%   MMregcore        - Computes MM regression estimators for a selected fixed scale                                                               - REG-Regression     - 2021 Feb 01
%   MMregeda         - Computes MM estimator in linear regression for a series of values of efficiency                                            - REG-Regression     - 2021 Feb 01
%   regressB         - Computes Bayesian estimates of regression parameters                                                                       - REG-Bayes          - 2021 Feb 01
%   regressH         - Fits a multiple linear regression model with heteroskedasticity                                                            - REG-Hetero         - 2021 Feb 01
%   regressHart      - Fits a multiple linear regression model using ART heteroskedasticity                                                       - REG-Hetero         - 2021 Feb 01
%   regressHhar      - Fits a multiple linear regression model with Harvey heteroskedasticity                                                     - REG-Hetero         - 2021 Feb 01
%   regressts        - Computes estimates of regression parameters for a time series models                                                       - REG-Regression     - 2021 Feb 01
%   rlsmo            - Computes a running-lines smoother with global cross-validation                                                             - REG-Transformations- 2021 Feb 01
%   RobCov           - Computes covariance matrix of robust regression coefficients                                                               - REG-Regression     - 2021 Feb 01
%   RobRegrSize      - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size               - REG-Regression     - 2021 Feb 01
%   Score            - Computes the score test for transformation                                                                                 - REG-Transformations- 2021 Feb 01
%   ScoreYJ          - Computes the score test for Yeo and Johnson transformation                                                                 - REG-Transformations- 2021 Feb 01
%   ScoreYJall       - Computes all the 4 score tests for YJ transformation                                                                       - REG-Transformations- 2021 Feb 01
%   ScoreYJmle       - Computes the likelihood ratio test fof H_0=lambdaP=lambdaP0 and lambdaN=lambdaN0                                           - REG-Transformations- 2021 Feb 01
%   ScoreYJpn        - Computes the score test for YJ transformation for pos and neg observations                                                 - REG-Transformations- 2021 Feb 01
%   SETARX           - Implements Threshold autoregressive models with two regimes                                                                - REG-Regression     - 2021 Feb 01
%   simulateLM       - Simulates linear regression data with prespecified values of statistical indexes                                           - REG-Regression     - 2021 Feb 01
%   simulateTS       - Simulates a time series with trend, time varying seasonal, level shift and irregular component                             - REG-Regression     - 2021 Feb 01
%   smothr           - Produces smoothed values with constraints                                                                                  - REG-Transformations- 2021 Feb 01
%   Sreg             - Computes S estimators in linear regression                                                                                 - REG-Regression     - 2021 Feb 01
%   Sregeda          - Computes S estimators in linear regression for a series of values of bdp                                                   - REG-Regression     - 2021 Feb 01
%   supsmu           - Smooths scatterplots using Friedman's supersmoother algorithm                                                              - REG-Transformations- 2021 Feb 01
%   Taureg           - Computes Tau estimators in linear regression                                                                               - REG-Regression     - 2021 Feb 01
%   tBothSides       - Allows users to transform both sides of a (nonlinear) regression model                                                     - REG-Regression     - 2021 Feb 01
%   VIOM             - Computes weights estimates under Variance-Inflation Model                                                                  - REG-Regression     - 2021 Feb 01
