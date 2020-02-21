% REGRESSION
%
% File names, description, category and date last modified
%
%   Name             - Description                                                                                                                - Category           - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   ace              - Computes alternative conditional expectation                                                                               - REG-Transformations- 2020 Jan 29
%   addt             - Produces the t test for an additional explanatory variable                                                                 - REG-Regression     - 2020 Feb 20
%   avas             - Computes additivity and variance stabilization for regression                                                              - REG-Transformations- 2020 Jan 29
%   forecastTS       - Forecast for a time series with trend, time varying seasonal, level shift and irregular component                          - REG-Regression     - 2020 Jan 29
%   FSR              - Gives an automatic outlier detection procedure in linear regression                                                        - REG-Regression     - 2020 Jan 29
%   FSRaddt          - Produces t deletion tests for each explanatory variable                                                                    - REG-ModelSelection - 2020 Feb 05
%   FSRB             - Gives an automatic outlier detection procedure in Bayesian linear regression                                               - REG-Bayes          - 2020 Feb 05
%   FSRBbsb          - Returns the units belonging to the subset in each step of the Bayesian forward search                                      - REG-Regression     - 2020 Feb 05
%   FSRBeda          - Enables to monitor several quantities in each step of the Bayesian search                                                  - REG-Bayes          - 2020 Feb 20
%   FSRBmdr          - Computes minimum deletion residual and other basic linear regression quantities in each step of the Bayesian search        - REG-Bayes          - 2020 Feb 05
%   FSRBr            - Bayesian forward search in linear regression reweighted                                                                    - REG-Bayes          - 2020 Feb 05
%   FSRbsb           - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2020 Feb 05
%   FSRcp            - Monitors Cp and AIC for all models of interest of size smallp                                                              - REG-ModelSelection - 2020 Feb 05
%   FSReda           - Enables to monitor several quantities in each step of the forward search                                                   - REG-Regression     - 2020 Feb 05
%   FSRenvmdr        - Computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search                           - REG-Regression     - 2020 Jan 29
%   FSRfan           - Monitors the values of the score test statistic for each lambda                                                            - REG-Transformations- 2020 Feb 18
%   FSRH             - Gives an automatic outlier detection procedure in heteroskedastic linear regression                                        - REG-Hetero         - 2020 Feb 18
%   FSRHbsb          - Returns the units belonging to the subset in each step of the heteroskedastic forward search                               - REG-Hetero         - 2020 Feb 05
%   FSRHeda          - Enables to monitor several quantities in each step of the forward search                                                   - REG-Hetero         - 2020 Feb 18
%   FSRHmdr          - Computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search - REG-Hetero         - 2020 Feb 05
%   FSRinvmdr        - Converts values of minimum deletion residual into confidence levels                                                        - REG-Regression     - 2020 Jan 29
%   FSRmdr           - Computes minimum deletion residual and other basic linear regression quantities in each step of the search                 - REG-Regression     - 2020 Feb 05
%   FSRms            - Performs robust model selection using flexible trimming in linear regression                                               - REG-ModelSelection - 2020 Feb 18
%   FSRr             - Forward search in linear regression reweighted                                                                             - REG-Regression     - 2020 Feb 20
%   FSRts            - Is an automatic adaptive procedure to detect outliers in time series                                                       - REG-Regression     - 2020 Feb 20
%   FSRtsbsb         - Returns the units belonging to the subset in each step of the forward search                                               - REG-Regression     - 2020 Jan 29
%   FSRtsmdr         - Computes minimum deletion residual for time series models in each step of the search                                       - REG-Regression     - 2020 Feb 18
%   GYfilt           - Computes the Gervini-Yohai univariate outlier identifier                                                                   - UTISTAT            - 2020 Jan 29
%   LTSts            - Extends LTS estimator to time series                                                                                       - REG-Regression     - 2020 Feb 21
%   LTStsVarSel      - Does variable selection in the robust time series model LTSts                                                              - REG-Regression     - 2020 Feb 18
%   LXS              - Computes the Least Median of Squares (LMS) or Least Trimmed Squares (LTS) estimators                                       - REG-Regression     - 2020 Feb 20
%   mdpdR            - Allows to apply Minimum Density Power Divergence criterion to parametric regression problems                               - REG-Regression     - 2020 Feb 20
%   mdpdReda         - Allows to monitor  Minimum Density Power Divergence criterion to parametric regression problems                            - REG-Regression     - 2020 Feb 20
%   MMreg            - Computes MM estimator of regression coefficients                                                                           - REG-Regression     - 2020 Feb 13
%   MMregcore        - Computes MM regression estimators for a selected fixed scale                                                               - REG-Regression     - 2020 Feb 13
%   MMregeda         - Computes MM estimator in linear regression for a series of values of efficiency                                            - REG-Regression     - 2020 Feb 20
%   regressB         - Computes Bayesian estimates of regression parameters                                                                       - REG-Bayes          - 2020 Feb 05
%   regressH         - Fits a multiple linear regression model with heteroskedasticity                                                            - REG-Hetero         - 2020 Feb 05
%   regressHart      - Fits a multiple linear regression model using ART heteroskedasticity                                                       - REG-Hetero         - 2020 Feb 05
%   regressHhar      - Fits a multiple linear regression model with Harvey heteroskedasticity                                                     - REG-Hetero         - 2020 Feb 20
%   regressts        - Computes estimates of regression parameters for a time series models                                                       - REG-Regression     - 2020 Feb 18
%   rlsmo            - Computes a running-lines smoother with global cross-validation                                                             - REG-Transformations- 2020 Feb 05
%   RobCov           - Computes covariance matrix of robust regression coefficients                                                               - REG-Regression     - 2020 Feb 20
%   RobRegrSize      - Provides proper threshold for robust estimators to obtain an empirical size close to 1 per cent nominal size               - REG-Regression     - 2020 Feb 05
%   Score            - Computes the score test for transformation                                                                                 - REG-Transformations- 2020 Feb 20
%   ScoreYJ          - Computes the score test for Yeo and Johnson transformation                                                                 - REG-Transformations- 2020 Feb 20
%   ScoreYJpn        - Computes the score test for YJ transformation separately for pos and neg observations                                      - REG-Transformations- 2020 Feb 05
%   simulateLM       - Simulate linear regression data                                                                                            - REG-Regression     - 2020 Feb 18
%   simulateTS       - Simulate a time series with trend, time varying seasonal, level shift and irregular component                              - REG-Regression     - 2020 Feb 18
%   smothr           - Produces smoothed values with constraints                                                                                  - REG-Transformations- 2020 Jan 29
%   Sreg             - Computes S estimators in linear regression                                                                                 - REG-Regression     - 2020 Feb 20
%   Sregeda          - Computes S estimators in linear regression for a series of values of bdp                                                   - REG-Regression     - 2020 Feb 20
%   supsmu           - Smooths scatterplots using Friedman's supersmoother algorithm                                                              - REG-Transformations- 2020 Feb 05
%   Taureg           - Computes Tau estimators in linear regression                                                                               - REG-Regression     - 2020 Feb 20
%   tBothSides       - Allows users to transform both sides of a (nonlinear) regression model                                                     - REG-Regression     - 2020 Feb 18
%   VIOM             - Computes weights estimates under Variance-Inflation Model                                                                  - REG-Regression     - 2020 Feb 20
