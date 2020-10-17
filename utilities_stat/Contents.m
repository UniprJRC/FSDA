% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                - Description                                                                                                                     - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   basicPower          - Computes the basic power transformation                                                                                         - UTISTAT         - 2020 Feb 05
%   bwe                 - Estimates the bandwidth smoothing parameter for kernel density estimation                                                       - UTISTAT         - 2020 Feb 05
%   crosstab2datamatrix - Recreates the original data matrix X from contingency table N                                                                   - MULT-Categorical- 2020 Feb 05
%   ctsub               - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                      - UTISTAT         - 2020 Feb 18
%   ellipse             - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                      - UTISTAT         - 2020 Feb 05
%   exactcdf            - Finds exact p-values                                                                                                            - UTISTAT         - 2020 Jan 29
%   FowlkesMallowsIndex - Computes the Fowlkes and Mallows index                                                                                          - UTISTAT         - 2020 Jan 29
%   FSMbonfbound        - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                      - UTISTAT         - 2020 Feb 05
%   FSRbonfbound        - Computes Bonferroni bounds for each step of the search (in linear regression)                                                   - UTISTAT         - 2020 Feb 05
%   HAbdp               - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2020 Jan 29
%   HAc                 - Computes breakdown point and efficiency associated with constant c                                                              - UTISTAT         - 2020 Feb 05
%   HAeff               - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                    - UTISTAT         - 2020 Jan 29
%   HApsi               - Computes psi function  using Hampel proposal                                                                                    - UTISTAT         - 2020 Feb 05
%   HApsider            - Computes derivative of psi function  using Hampel proposal                                                                      - UTISTAT         - 2020 Feb 05
%   HApsix              - Computes psi function  using Hampel proposal times x                                                                            - UTISTAT         - 2020 Feb 05
%   HArho               - Computes rho function  using Hampel proposal                                                                                    - UTISTAT         - 2020 Feb 05
%   HAwei               - Computes weight function psi(u)/u using Hampel proposal                                                                         - UTISTAT         - 2020 Jan 29
%   HUeff               - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                               - UTISTAT         - 2020 Jan 29
%   HUpsi               - Computes psi function (derivative of rho function) for Huber                                                                    - UTISTAT         - 2020 Jan 29
%   HUpsider            - Computes derivative of psi function (second derivative of rho function) for Huber                                               - UTISTAT         - 2020 Feb 05
%   HUpsix              - Computes psi function (derivative of rho function) times x for Huber                                                            - UTISTAT         - 2020 Feb 05
%   HUrho               - Computes rho function for Huber                                                                                                 - UTISTAT         - 2020 Feb 13
%   HUwei               - Computes weight function psi(u)/u for Huber                                                                                     - UTISTAT         - 2020 Feb 05
%   HYPbdp              - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                   - UTISTAT         - 2020 Feb 20
%   HYPc                - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC) - UTISTAT         - 2020 Feb 05
%   HYPck               - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                         - UTISTAT         - 2020 Feb 05
%   HYPeff              - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                       - UTISTAT         - 2020 Feb 05
%   HYPk                - Computes breakdown point and efficiency for hyp. tan. estimator                                                                 - UTISTAT         - 2020 Feb 05
%   HYPpsi              - Computes psi function for hyperbolic tangent estimator                                                                          - UTISTAT         - 2020 Jan 29
%   HYPpsider           - Computes derivative of psi function for hyperbolic tangent estimator                                                            - UTISTAT         - 2020 Jan 29
%   HYPpsix             - Computes psi function for hyperbolic tangent estimator times x                                                                  - UTISTAT         - 2020 Jan 29
%   HYPrho              - Computes rho function  using hyperbolic tangent estimator                                                                       - UTISTAT         - 2020 Jan 29
%   HYPwei              - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                              - UTISTAT         - 2020 Jan 29
%   inversegamcdf       - Computes inverse-gamma cumulative distribution function                                                                         - UTISTAT         - 2020 Jan 29
%   inversegaminv       - Inversegampdf Inverse-gamma cumulative distribution function                                                                    - UTISTAT         - 2020 Jan 29
%   inversegampdf       - Computes inverse-gamma probability density function                                                                             - UTISTAT         - 2020 Jan 29
%   kdebiv              - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                  - UTISTAT         - 2020 Feb 18
%   logmvnpdfFS         - Produces log of Multivariate normal probability density function (pdf)                                                          - UTISTAT         - 2020 Jul 03
%   mahalFS             - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                      - UTISTAT         - 2020 Feb 05
%   mdpd                - Computes Minimum Distance Power Divergence statistics                                                                           - UTISTAT         - 2020 Apr 21
%   Mscale              - Finds the M estimator of the scale                                                                                              - UTISTAT         - 2020 Feb 20
%   mtR                 - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                        - UTISTAT         - 2020 Feb 18
%   ncpci               - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                - UTISTAT         - 2020 Feb 05
%   ncx2mixtcdf         - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                     - UTISTAT         - 2020 Sep 11
%   normBoxCox          - Computes (normalized) Box-Cox transformation                                                                                    - UTISTAT         - 2020 Feb 05
%   normYJ              - Computes (normalized) Yeo-Johnson transformation                                                                                - UTISTAT         - 2020 May 12
%   normYJpn            - NormYJ computes (normalized) extended Yeo-Johnson transformation                                                                - UTISTAT         - 2020 May 12
%   OPTbdp              - Finds the constant c associated to the supplied breakdown point                                                                 - UTISTAT         - 2020 Jan 29
%   OPTc                - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                     - UTISTAT         - 2020 Feb 20
%   OPTeff              - Finds the constant c which is associated to the requested efficiency                                                            - UTISTAT         - 2020 Apr 21
%   OPTpsi              - Computes psi function (derivative of rho function) for optimal weight function                                                  - UTISTAT         - 2020 Jan 29
%   OPTpsider           - Computes derivative of psi function (second derivative of rho function) for optimal weight function                             - UTISTAT         - 2020 Jan 29
%   OPTpsix             - Computes psi function (derivative of rho function) times x                                                                      - UTISTAT         - 2020 Jan 29
%   OPTrho              - Computes rho function for optimal weight function                                                                               - UTISTAT         - 2020 Jan 29
%   OPTwei              - Computes weight function psi(u)/u for optimal weight function                                                                   - UTISTAT         - 2020 Oct 11
%   PDbdp               - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                      - UTISTAT         - 2020 Aug 25
%   PDc                 - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator            - UTISTAT         - 2020 Feb 07
%   PDeff               - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                 - UTISTAT         - 2020 Feb 11
%   PDpsi               - Computes psi function (derivative of rho function) for minimum density power divergence estimator                               - UTISTAT         - 2020 Feb 07
%   PDpsider            - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                  - UTISTAT         - 2020 Apr 03
%   PDpsix              - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                       - UTISTAT         - 2020 Feb 07
%   PDrho               - Computes rho function for minimum density power divergence estimator                                                            - UTISTAT         - 2020 Feb 07
%   PDwei               - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                           - UTISTAT         - 2020 Aug 25
%   Powertra            - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                         - UTISTAT         - 2020 Feb 05
%   Qn                  - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                  - UTISTAT         - 2020 Jan 29
%   qqplotFS            - Qqplot of studentized residuals with envelopes                                                                                  - VIS-Reg         - 2020 Feb 05
%   RandIndexFS         - Calculates Rand type Indices to compare two partitions                                                                          - UTISTAT         - 2020 Aug 07
%   RKbdp               - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                 - UTISTAT         - 2020 Jan 29
%   RKeff               - Finds the constants c and M which are associated to the requested efficiency and ARP                                            - UTISTAT         - 2020 Apr 21
%   RKpsi               - Computes psi function for Rocke (translated Tukey's) biweight                                                                   - UTISTAT         - 2020 Jan 29
%   RKpsider            - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                 - UTISTAT         - 2020 Jan 29
%   RKpsix              - Computes psi function times x for Rocke (translated Tukey's) biweight                                                           - UTISTAT         - 2020 Jan 29
%   RKrho               - Computes rho function for Rocke (translated Tukey's) biweight                                                                   - UTISTAT         - 2020 Jan 29
%   RKwei               - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                       - UTISTAT         - 2020 Jan 29
%   rthin               - Applies independent random thinning to a point pattern                                                                          - UTISTAT         - 2020 Jan 29
%   Sn                  - Robust estimator of scale (robust version of Gini's average difference)                                                         - UTISTAT         - 2020 Jan 29
%   tabulateFS          - Creates frequency table of unique values of x, excluding possible 0 counts                                                      - UTISTAT         - 2020 Jan 29
%   TBbdp               - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                            - UTISTAT         - 2020 Oct 11
%   TBc                 - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                         - UTISTAT         - 2020 Jan 29
%   TBeff               - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                       - UTISTAT         - 2020 Feb 05
%   TBpsi               - Computes psi function (derivative of rho function) for Tukey's biweight                                                         - UTISTAT         - 2020 Aug 30
%   TBpsider            - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                    - UTISTAT         - 2020 Feb 05
%   TBpsix              - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                 - UTISTAT         - 2020 Jan 29
%   TBrho               - Computes rho function for Tukey's biweight                                                                                      - UTISTAT         - 2020 Jan 29
%   TBwei               - Computes weight function psi(u)/u for Tukey's biweight                                                                          - UTISTAT         - 2020 Jan 29
%   twdcdf              - Computes the cumulative distribution function of the Tweedie distribution                                                       - UTISTAT         - 2020 May 06
%   twdpdf              - Computes the probability density function of the Tweedie distribution                                                           - UTISTAT         - 2020 Apr 22
%   twdrnd              - Generates random variates from the Tweedie distribution                                                                         - UTISTAT         - 2020 May 06
%   vervaatrnd          - Simulates random variates from the Vervaat perpetuity distribution                                                              - UTISTAT         - 2020 Jan 29
%   vervaatsim          - Returns a Vervaat perpetuity                                                                                                    - UTISTAT         - 2020 Jan 29
%   vervaatxdf          - Returns the pdf and cdf of a Vervaat perpetuity                                                                                 - UTISTAT         - 2020 Feb 05
%   winsor              - Returns a winsorized copy of input                                                                                              - UTISTAT         - 2020 Apr 21
%   WNChygepdf          - Returns Wallenius' non-central hypergeometric probability density values                                                        - UTISTAT         - 2020 Jan 29
%   wthin               - Thins a uni/bi-dimensional dataset                                                                                              - UTISTAT         - 2020 Jul 04
