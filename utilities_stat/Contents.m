% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                - Description                                                                                                                     - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   basicPower          - Computes the basic power transformation                                                                                         - UTISTAT         - 2017 Nov 17
%   bwe                 - Estimates the bandwidth smoothing parameter for kernel density estimation                                                       - UTISTAT         - 2017 Nov 17
%   crosstab2datamatrix - Recreates the original data matrix X from contingency table N                                                                   - MULT-Categorical- 2017 Nov 17
%   ellipse             - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT         - 2017 Nov 17
%   FowlkesMallowsIndex - Computes the Fowlkes and Mallows index                                                                                          - UTISTAT         - 2017 Nov 17
%   FSMbonfbound        - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT         - 2017 Nov 17
%   FSRbonfbound        - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT         - 2017 Nov 17
%   HAbdp               - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2017 Nov 17
%   HAc                 - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT         - 2017 Nov 17
%   HAeff               - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT         - 2017 Nov 17
%   HApsi               - Computes psi function  using Hampel proposal                                                                                    - UTISTAT         - 2017 Nov 17
%   HApsider            - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT         - 2017 Nov 17
%   HApsix              - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT         - 2017 Nov 17
%   HArho               - Computes rho function  using Hampel proposal                                                                                    - UTISTAT         - 2017 Nov 17
%   HAwei               - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT         - 2017 Nov 17
%   HUeff               - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT         - 2017 Nov 17
%   HUpsi               - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT         - 2017 Nov 17
%   HUpsider            - Computes derivative of psi function (second derivative of rho function) for Huber                                               - UTISTAT         - 2017 Nov 17
%   HUpsix              - Computes psi function (derivative of rho function) times x for Huber                                                            - UTISTAT         - 2017 Nov 17
%   HUrho               - Computes rho function for Huber                                                                                                 - UTISTAT         - 2017 Nov 17
%   HUwei               - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT         - 2017 Nov 17
%   HYPbdp              - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT         - 2017 Nov 17
%   HYPc                - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT         - 2017 Nov 17
%   HYPck               - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT         - 2017 Nov 17
%   HYPeff              - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT         - 2017 Nov 17
%   HYPk                - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT         - 2017 Nov 17
%   HYPpsi              - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT         - 2017 Nov 17
%   HYPpsider           - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT         - 2017 Nov 17
%   HYPpsix             - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT         - 2017 Nov 17
%   HYPrho              - Computes rho function  using hyperbolic tangent estimator                                                                       - UTISTAT         - 2017 Nov 17
%   HYPwei              - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT         - 2017 Nov 17
%   inversegamcdf       - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT         - 2017 Nov 17
%   inversegaminv       - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT         - 2017 Nov 17
%   inversegampdf       - Computes inverse-gamma probability density function                                                                             - UTISTAT         - 2017 Nov 17
%   kdebiv              - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                  - UTISTAT         - 2017 Nov 17
%   logmvnpdfFS         - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT         - 2017 Nov 17
%   mahalFS             - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT         - 2017 Nov 17
%   Mscale              - Finds the M estimator of the scale                                                                                              - UTISTAT         - 2017 Nov 17
%   ncpci               - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                - UTISTAT         - 2017 Nov 19
%   ncx2mixtcdf         - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT         - 2017 Nov 17
%   normBoxCox          - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT         - 2017 Nov 17
%   normYJ              - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT         - 2017 Nov 17
%   OPTbdp              - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2017 Nov 17
%   OPTc                - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT         - 2017 Nov 17
%   OPTeff              - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT         - 2017 Nov 17
%   OPTpsi              - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT         - 2017 Nov 17
%   OPTpsider           - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT         - 2017 Nov 17
%   OPTpsix             - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT         - 2017 Nov 17
%   OPTrho              - Computes rho function for optimal weight function                                                                               - UTISTAT         - 2017 Nov 17
%   OPTwei              - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT         - 2017 Nov 17
%   Powertra            - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT         - 2017 Nov 17
%   Qn                  - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT         - 2017 Nov 17
%   RandIndexFS         - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT         - 2017 Nov 17
%   rthin               - Applies independent random thinning to a point pattern                                                                          - UTISTAT         - 2017 Nov 17
%   Sn                  - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT         - 2017 Nov 17
%   tabulateFS          - Create frequency table of unique values of x, excluding possible 0 counts                                                       - UTISTAT         - 2017 Nov 17
%   TBbdp               - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT         - 2017 Nov 17
%   TBc                 - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT         - 2017 Nov 17
%   TBeff               - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT         - 2017 Nov 17
%   TBpsi               - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT         - 2017 Nov 17
%   TBpsider            - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT         - 2017 Nov 17
%   TBpsix              - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT         - 2017 Nov 17
%   TBrho               - Computes rho function for Tukey's biweight                                                                                      - UTISTAT         - 2017 Nov 17
%   TBwei               - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT         - 2017 Nov 17
%   winsor              - Returns a winsorized copy of input                                                                                              - UTISTAT         - 2017 Nov 17
%   WNChygepdf          - Returns Wallenius' non-central hypergeometric probability density values                                                        - UTISTAT         - 2017 Nov 17
%   wthin               - Thin a uni/bi-dimensional dataset                                                                                               - UTISTAT         - 2017 Nov 17
