% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                      - Description                                                                                                                     - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   basicPower                - Computes the basic power transformation                                                                                         - UTISTAT         - 2018 Jun 08
%   bwe                       - Estimates the bandwidth smoothing parameter for kernel density estimation                                                       - UTISTAT         - 2018 Sep 15
%   crosstab2datamatrix       - Recreates the original data matrix X from contingency table N                                                                   - MULT-Categorical- 2018 May 31
%   ellipse                   - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT         - 2018 Sep 15
%   FowlkesMallowsIndex       - Computes the Fowlkes and Mallows index                                                                                          - UTISTAT         - 2018 Jun 08
%   FSMbonfbound              - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT         - 2018 Jun 08
%   FSRbonfbound              - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT         - 2018 Jun 08
%   HAbdp                     - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2018 Jun 08
%   HAc                       - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT         - 2018 Jun 08
%   HAeff                     - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT         - 2018 Sep 20
%   HApsi                     - Computes psi function  using Hampel proposal                                                                                    - UTISTAT         - 2018 Jun 08
%   HApsider                  - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT         - 2018 Jun 08
%   HApsix                    - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT         - 2018 Jun 08
%   HArho                     - Computes rho function  using Hampel proposal                                                                                    - UTISTAT         - 2018 Jun 08
%   HAwei                     - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT         - 2018 Jun 08
%   HUeff                     - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT         - 2018 Jun 08
%   HUpsi                     - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT         - 2018 Jun 08
%   HUpsider                  - Computes derivative of psi function (second derivative of rho function) for Huber                                               - UTISTAT         - 2018 Jul 13
%   HUpsix                    - Computes psi function (derivative of rho function) times x for Huber                                                            - UTISTAT         - 2018 Jun 08
%   HUrho                     - Computes rho function for Huber                                                                                                 - UTISTAT         - 2018 Jun 08
%   HUwei                     - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT         - 2018 Jun 08
%   HYPbdp                    - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT         - 2018 Sep 15
%   HYPc                      - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT         - 2018 Jun 08
%   HYPck                     - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT         - 2018 Jun 08
%   HYPeff                    - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT         - 2018 Sep 15
%   HYPk                      - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT         - 2018 Jun 08
%   HYPpsi                    - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT         - 2018 Jun 08
%   HYPpsider                 - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT         - 2018 Jun 08
%   HYPpsix                   - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT         - 2018 Jun 08
%   HYPrho                    - Computes rho function  using hyperbolic tangent estimator                                                                       - UTISTAT         - 2018 Jun 08
%   HYPwei                    - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT         - 2018 Jun 08
%   inversegamcdf             - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT         - 2018 May 31
%   inversegaminv             - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT         - 2018 May 31
%   inversegampdf             - Computes inverse-gamma probability density function                                                                             - UTISTAT         - 2018 Sep 15
%   kdebiv                    - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                  - UTISTAT         - 2018 Jun 08
%   logmvnpdfFS               - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT         - 2018 Jun 19
%   mahalFS                   - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT         - 2018 May 31
%   Mscale                    - Finds the M estimator of the scale                                                                                              - UTISTAT         - 2018 Sep 15
%   ncpci                     - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                - UTISTAT         - 2018 Sep 15
%   ncx2mixtcdf               - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT         - 2018 Jun 08
%   normBoxCox                - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT         - 2018 Jun 08
%   normYJ                    - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT         - 2018 Sep 15
%   OPTbdp                    - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2018 Jun 08
%   OPTc                      - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT         - 2018 Jun 08
%   OPTeff                    - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT         - 2018 Jun 08
%   OPTpsi                    - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT         - 2018 Jun 08
%   OPTpsider                 - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT         - 2018 Jun 08
%   OPTpsix                   - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT         - 2018 Jun 08
%   OPTrho                    - Computes rho function for optimal weight function                                                                               - UTISTAT         - 2018 Sep 15
%   OPTwei                    - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT         - 2018 Sep 15
%   Powertra                  - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT         - 2018 Sep 15
%   Qn                        - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT         - 2018 Sep 15
%   qqplotFS                  - Qqplot of studentized residuals with envelopes                                                                                  - VIS-Reg         - 2018 Sep 20
%   RandIndexFS               - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT         - 2018 Jun 08
%   RKbdp                     - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                 - UTISTAT         - 2018 Jul 08
%   RKeff                     - ROeff finds the constants c and M which are associated to the requested efficiency and ARP                                      - UTISTAT         - 2018 Jul 09
%   RKpsi                     - Computes psi function for Rocke (translated Tukey's) biweight                                                                   - UTISTAT         - 2018 Jul 06
%   RKpsider                  - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                 - UTISTAT         - 2018 Jul 06
%   RKpsix                    - Computes psi function times x for Rocke (translated Tukey's) biweight                                                           - UTISTAT         - 2018 Jul 06
%   RKrho                     - Computes rho function for Rocke (translated Tukey's) biweight                                                                   - UTISTAT         - 2018 Sep 17
%   RKwei                     - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                       - UTISTAT         - 2018 Sep 03
%   rthin                     - Applies independent random thinning to a point pattern                                                                          - UTISTAT         - 2018 Sep 20
%   Sn                        - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT         - 2018 Sep 15
%   tabulateFS                - Creates frequency table of unique values of x, excluding possible 0 counts                                                      - UTISTAT         - 2018 May 31
%   TBbdp                     - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT         - 2018 Jul 02
%   TBc                       - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT         - 2018 Jun 08
%   TBeff                     - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT         - 2018 Jun 08
%   TBpsi                     - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT         - 2018 Jun 08
%   TBpsider                  - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT         - 2018 Jun 08
%   TBpsix                    - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT         - 2018 Jun 08
%   TBrho                     - Computes rho function for Tukey's biweight                                                                                      - UTISTAT         - 2018 Jun 08
%   TBwei                     - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT         - 2018 Sep 20
%   winsor                    - Returns a winsorized copy of input                                                                                              - UTISTAT         - 2018 Jun 08
%   WNChygepdf                - Returns Wallenius' non-central hypergeometric probability density values                                                        - UTISTAT         - 2018 Jun 08
%   wthin                     - Thins a uni/bi-dimensional dataset                                                                                              - UTISTAT         - 2018 Jun 08
