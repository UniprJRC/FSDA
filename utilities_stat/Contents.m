% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                - Description                                                                                                                     - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   basicPower          - Computes the basic power transformation                                                                                         - UTISTAT         - 2020 May 06
%   bwe                 - Estimates the bandwidth smoothing parameter for kernel density estimation                                                       - UTISTAT         - 2020 May 06
%   crosstab2datamatrix - Recreates the original data matrix X from contingency table N                                                                   - MULT-Categorical- 2020 May 06
%   ctsub               - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                      - UTISTAT         - 2020 May 06
%   ellipse             - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT         - 2020 May 06
%   exactcdf            - Finds exact p-values                                                                                                            - UTISTAT         - 2020 May 06
%   FowlkesMallowsIndex - Computes the Fowlkes and Mallows index                                                                                          - UTISTAT         - 2020 May 06
%   FSMbonfbound        - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT         - 2020 May 06
%   FSRbonfbound        - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT         - 2020 May 06
%   HAbdp               - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2020 May 06
%   HAc                 - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT         - 2020 May 06
%   HAeff               - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT         - 2020 May 06
%   HApsi               - Computes psi function  using Hampel proposal                                                                                    - UTISTAT         - 2020 May 06
%   HApsider            - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT         - 2020 May 06
%   HApsix              - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT         - 2020 May 06
%   HArho               - Computes rho function  using Hampel proposal                                                                                    - UTISTAT         - 2020 May 06
%   HAwei               - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT         - 2020 May 06
%   HUeff               - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT         - 2020 May 06
%   HUpsi               - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT         - 2020 May 06
%   HUpsider            - Computes derivative of psi function (second derivative of rho function) for Huber                                               - UTISTAT         - 2020 May 06
%   HUpsix              - Computes psi function (derivative of rho function) times x for Huber                                                            - UTISTAT         - 2020 May 06
%   HUrho               - Computes rho function for Huber                                                                                                 - UTISTAT         - 2020 May 06
%   HUwei               - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT         - 2020 May 06
%   HYPbdp              - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT         - 2020 May 06
%   HYPc                - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT         - 2020 May 06
%   HYPck               - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT         - 2020 May 06
%   HYPeff              - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT         - 2020 May 06
%   HYPk                - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT         - 2020 May 06
%   HYPpsi              - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT         - 2020 May 06
%   HYPpsider           - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT         - 2020 May 06
%   HYPpsix             - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT         - 2020 May 06
%   HYPrho              - Computes rho function  using hyperbolic tangent estimator                                                                       - UTISTAT         - 2020 May 06
%   HYPwei              - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT         - 2020 May 06
%   inversegamcdf       - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT         - 2020 May 06
%   inversegaminv       - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT         - 2020 May 06
%   inversegampdf       - Computes inverse-gamma probability density function                                                                             - UTISTAT         - 2020 May 06
%   kdebiv              - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                  - UTISTAT         - 2020 May 06
%   logmvnpdfFS         - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT         - 2020 May 27
%   mahalFS             - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT         - 2020 May 06
%   mdpd                - Computes Minimum Distance Power Divergence statistics                                                                           - UTISTAT         - 2020 May 06
%   Mscale              - Finds the M estimator of the scale                                                                                              - UTISTAT         - 2020 May 06
%   mtR                 - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                        - UTISTAT         - 2020 May 06
%   ncpci               - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                - UTISTAT         - 2020 May 06
%   ncx2mixtcdf         - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT         - 2020 May 06
%   normBoxCox          - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT         - 2020 May 06
%   normYJ              - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT         - 2020 May 17
%   normYJpn            - NormYJ computes (normalized) extended Yeo-Johnson transformation                                                                - UTISTAT         - 2020 May 17
%   OPTbdp              - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2020 May 06
%   OPTc                - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT         - 2020 May 06
%   OPTeff              - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT         - 2020 May 06
%   OPTpsi              - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT         - 2020 May 06
%   OPTpsider           - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT         - 2020 May 06
%   OPTpsix             - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT         - 2020 May 06
%   OPTrho              - Computes rho function for optimal weight function                                                                               - UTISTAT         - 2020 May 06
%   OPTwei              - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT         - 2020 May 06
%   PDbdp               - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                      - UTISTAT         - 2020 May 06
%   PDc                 - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator            - UTISTAT         - 2020 May 06
%   PDeff               - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                 - UTISTAT         - 2020 May 06
%   PDpsi               - Computes psi function (derivative of rho function) for minimum density power divergence estimator                               - UTISTAT         - 2020 May 06
%   PDpsider            - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                  - UTISTAT         - 2020 May 06
%   PDpsix              - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                       - UTISTAT         - 2020 May 06
%   PDrho               - Computes rho function for minimum density power divergence estimator                                                            - UTISTAT         - 2020 May 06
%   PDwei               - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                           - UTISTAT         - 2020 May 06
%   Powertra            - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT         - 2020 May 06
%   Qn                  - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT         - 2020 May 06
%   qqplotFS            - Qqplot of studentized residuals with envelopes                                                                                  - VIS-Reg         - 2020 May 06
%   RandIndexFS         - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT         - 2020 May 06
%   RKbdp               - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                 - UTISTAT         - 2020 May 06
%   RKeff               - Finds the constants c and M which are associated to the requested efficiency and ARP                                            - UTISTAT         - 2020 May 06
%   RKpsi               - Computes psi function for Rocke (translated Tukey's) biweight                                                                   - UTISTAT         - 2020 May 06
%   RKpsider            - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                 - UTISTAT         - 2020 May 06
%   RKpsix              - Computes psi function times x for Rocke (translated Tukey's) biweight                                                           - UTISTAT         - 2020 May 06
%   RKrho               - Computes rho function for Rocke (translated Tukey's) biweight                                                                   - UTISTAT         - 2020 May 06
%   RKwei               - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                       - UTISTAT         - 2020 May 06
%   rthin               - Applies independent random thinning to a point pattern                                                                          - UTISTAT         - 2020 May 06
%   Sn                  - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT         - 2020 May 06
%   tabulateFS          - Creates frequency table of unique values of x, excluding possible 0 counts                                                      - UTISTAT         - 2020 May 06
%   TBbdp               - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT         - 2020 May 06
%   TBc                 - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT         - 2020 May 06
%   TBeff               - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT         - 2020 May 06
%   TBpsi               - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT         - 2020 May 06
%   TBpsider            - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT         - 2020 May 06
%   TBpsix              - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT         - 2020 May 06
%   TBrho               - Computes rho function for Tukey's biweight                                                                                      - UTISTAT         - 2020 May 06
%   TBwei               - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT         - 2020 May 06
%   twdcdf              - Computes the cumulative distribution function of the Tweedie distribution                                                       - UTISTAT         - 2020 May 06
%   twdpdf              - Computes the probability density function of the Tweedie distribution                                                           - UTISTAT         - 2020 May 06
%   twdrnd              - Generates random variates from the Tweedie distribution                                                                         - UTISTAT         - 2020 May 06
%   vervaatrnd          - Simulates random variates from the Vervaat perpetuity distribution                                                              - UTISTAT         - 2020 May 06
%   vervaatsim          - Returns a Vervaat perpetuity                                                                                                    - UTISTAT         - 2020 May 06
%   vervaatxdf          - Returns the pdf and cdf of a Vervaat perpetuity                                                                                 - UTISTAT         - 2020 May 06
%   winsor              - Returns a winsorized copy of input                                                                                              - UTISTAT         - 2020 May 06
%   WNChygepdf          - Returns Wallenius' non-central hypergeometric probability density values                                                        - UTISTAT         - 2020 May 06
%   wthin               - Thins a uni/bi-dimensional dataset                                                                                              - UTISTAT         - 2020 May 06
