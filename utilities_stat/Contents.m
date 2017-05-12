% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                - Description                                                                                                                     - Category- Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   basicPower          - Computes the basic power transformation                                                                                         - UTISTAT- 2017 Mar 02
%   bwe                 - Estimates the bandwidth smoothing parameter for kernel density estimation                                                       - UTISTAT- 2017 May 12
%   ellipse             - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT- 2016 Aug 01
%   FowlkesMallowsIndex - Computes the Fowlkes and Mallows index                                                                                          - UTISTAT- 2017 Mar 13
%   FSMbonfbound        - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT- 2016 Aug 01
%   FSRbonfbound        - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT- 2016 Aug 01
%   HAbdp               - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT- 2016 Aug 01
%   HAc                 - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT- 2016 Aug 01
%   HAeff               - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT- 2016 Aug 01
%   HApsi               - Computes psi function  using Hampel proposal                                                                                    - UTISTAT- 2016 Aug 01
%   HApsider            - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT- 2016 Aug 01
%   HApsix              - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT- 2016 Aug 01
%   HArho               - Computes rho function  using Hampel proposal                                                                                    - UTISTAT- 2016 Aug 01
%   HAwei               - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT- 2016 Aug 01
%   HUeff               - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT- 2017 May 03
%   HUpsi               - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT- 2016 Aug 01
%   HUpsider            - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT- 2016 Aug 01
%   HUpsix              - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT- 2016 Aug 01
%   HUrho               - Computes rho function for Huber                                                                                                 - UTISTAT- 2016 Aug 01
%   HUwei               - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT- 2016 Aug 01
%   HYPbdp              - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT- 2016 Aug 01
%   HYPc                - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT- 2016 Aug 01
%   HYPck               - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT- 2016 Oct 02
%   HYPeff              - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT- 2016 Aug 01
%   HYPk                - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT- 2016 Aug 01
%   HYPpsi              - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT- 2016 Aug 01
%   HYPpsider           - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT- 2017 Feb 09
%   HYPpsix             - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT- 2016 Aug 01
%   HYPrho              - Computes rho function  using hyperboloc tangent estimator                                                                       - UTISTAT- 2016 Aug 01
%   HYPwei              - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT- 2017 Feb 17
%   inversegamcdf       - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT- 2016 Aug 01
%   inversegaminv       - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT- 2016 Aug 01
%   inversegampdf       - Computes inverse-gamma probability density function                                                                             - UTISTAT- 2016 Aug 01
%   kdebiv              - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                  - UTISTAT- 2017 May 12
%   logmvnpdfFS         - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT- 2016 Aug 01
%   mahalFS             - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT- 2017 May 02
%   Mscale              - Finds the M estimator of the scale                                                                                              - UTISTAT- 2016 Aug 01
%   ncx2mixtcdf         - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT- 2016 Aug 01
%   normBoxCox          - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT- 2017 Feb 20
%   normYJ              - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT- 2017 Mar 02
%   OPTbdp              - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT- 2016 Aug 01
%   OPTc                - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT- 2016 Aug 01
%   OPTeff              - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT- 2016 Aug 01
%   OPTpsi              - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT- 2016 Aug 01
%   OPTpsider           - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT- 2016 Aug 01
%   OPTpsix             - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT- 2016 Aug 01
%   OPTrho              - Computes rho function for optimal weight function                                                                               - UTISTAT- 2016 Aug 01
%   OPTwei              - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT- 2016 Aug 01
%   Powertra            - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT- 2017 Mar 02
%   Qn                  - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT- 2016 Aug 01
%   RandIndexFS         - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT- 2017 May 02
%   rthin               - Applies independent random thinning to a point pattern                                                                          - UTISTAT- 2017 May 12
%   Sn                  - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT- 2016 Aug 01
%   tabulateFS          - Create frequency table of unique values of x, excluding possible 0 counts                                                       - UTISTAT- 2016 Aug 01
%   TBbdp               - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT- 2016 Aug 01
%   TBc                 - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT- 2016 Aug 01
%   TBeff               - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT- 2016 Aug 01
%   TBpsi               - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT- 2016 Aug 01
%   TBpsider            - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT- 2017 May 10
%   TBpsix              - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT- 2016 Aug 01
%   TBrho               - Computes rho function for Tukey's biweight                                                                                      - UTISTAT- 2016 Aug 01
%   TBwei               - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT- 2016 Aug 01
%   winsor              - Returns a winsorized copy of input                                                                                              - UTISTAT- 2016 Aug 01
%   WNChygepdf          - Returns Wallenius' non-central hypergeometric probability density values                                                        - UTISTAT- 2017 May 04
%   wthin               - Thin a uni/bi-dimensional dataset                                                                                               - UTISTAT- 2017 May 12
