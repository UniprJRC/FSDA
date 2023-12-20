% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                    - Description                                                                                                                         - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   add2boxplot             - Adds labels to the boxplot figure                                                                                                   - VIS-Mult        - 2023 Dec 13
%   ASbdp                   - Finds the constant c associated to the supplied breakdown point for Andrew's sine function                                          - UTISTAT         - 2023 Dec 13
%   ASc                     - Computes breakdown point and efficiency associated with constant c for Andrew's rho function                                        - UTISTAT         - 2023 Dec 13
%   ASeff                   - Finds the constant c which is associated to the requested efficiency for Andrew's sine function                                     - UTISTAT         - 2023 Dec 13
%   ASpsi                   - Computes psi function (derivative of rho function) for Andrew's sine function                                                       - UTISTAT         - 2023 Dec 13
%   ASpsider                - Computes derivative of psi function (second derivative of rho function) for Andrew's sine function                                  - UTISTAT         - 2023 Dec 13
%   ASpsix                  - Computes psi function (derivative of rho function) times x for Andrew's sine function                                               - UTISTAT         - 2023 Dec 13
%   ASrho                   - Computes rho function for Andrew's sine function                                                                                    - UTISTAT         - 2023 Dec 13
%   ASwei                   - Computes weight function psi(u)/u for Andrew's sine function                                                                        - UTISTAT         - 2023 Dec 13
%   basicPower              - Computes the basic power transformation                                                                                             - UTISTAT         - 2023 Dec 13
%   bwe                     - Estimates the bandwidth smoothing parameter for kernel density estimation                                                           - UTISTAT         - 2023 Dec 13
%   CEVmodel                - Computes price and instantaneous variance processes from the CEV model                                                              - UTISTAT         - 2023 Dec 13
%   corrcdf                 - Computes correlation coefficient probability distribution function                                                                  - ProbDist        - 2023 Dec 13
%   corrpdf                 - Computes correlation coefficient probability density function                                                                       - ProbDist        - 2023 Dec 13
%   crosstab2datamatrix     - Recreates the original data matrix X from contingency table N                                                                       - MULT-Categorical- 2023 Dec 13
%   ctsub                   - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                          - UTISTAT         - 2023 Dec 13
%   ellipse                 - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                          - UTISTAT         - 2023 Dec 13
%   exactcdf                - Finds exact p-values                                                                                                                - UTISTAT         - 2023 Dec 13
%   FE_int_vol              - Computes the integrated variance from a diffusion process via the Fourier estimator using Dirichlet kernel                          - UTISTAT         - 2023 Dec 13
%   FE_int_vol_Fejer        - Computes the integrated variance from a diffusion process via the Fourier estimator using Fejer kernel                              - UTISTAT         - 2023 Dec 13
%   FE_spot_vol             - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel                          - UTISTAT         - 2023 Dec 13
%   FE_spot_vol_FFT         - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel, using the FFT algorithm - UTISTAT         - 2023 Dec 13
%   FowlkesMallowsIndex     - Computes the Fowlkes and Mallows index                                                                                              - UTISTAT         - 2023 Dec 13
%   FSMbonfbound            - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                          - UTISTAT         - 2023 Dec 13
%   FSRbonfbound            - Computes Bonferroni bounds for each step of the search (in linear regression)                                                       - UTISTAT         - 2023 Dec 13
%   HAbdp                   - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT         - 2023 Dec 13
%   HAc                     - Computes breakdown point and efficiency associated with constant c                                                                  - UTISTAT         - 2023 Dec 13
%   HAeff                   - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                        - UTISTAT         - 2023 Dec 13
%   HApsi                   - Computes psi function  using Hampel proposal                                                                                        - UTISTAT         - 2023 Dec 13
%   HApsider                - Computes derivative of psi function  using Hampel proposal                                                                          - UTISTAT         - 2023 Dec 13
%   HApsix                  - Computes psi function  using Hampel proposal times x                                                                                - UTISTAT         - 2023 Dec 13
%   HArho                   - Computes rho function  using Hampel proposal                                                                                        - UTISTAT         - 2023 Dec 13
%   HAwei                   - Computes weight function psi(u)/u using Hampel proposal                                                                             - UTISTAT         - 2023 Dec 13
%   HUc                     - Computes breakdown point and efficiency associated with constant c for Huber link                                                   - UTISTAT         - 2023 Dec 13
%   HUeff                   - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                   - UTISTAT         - 2023 Dec 13
%   HUpsi                   - Computes psi function (derivative of rho function) for Huber                                                                        - UTISTAT         - 2023 Dec 13
%   HUpsider                - Computes derivative of psi function (second derivative of rho function) for Huber                                                   - UTISTAT         - 2023 Dec 13
%   HUpsix                  - Computes psi function (derivative of rho function) times x for Huber                                                                - UTISTAT         - 2023 Dec 13
%   HUrho                   - Computes rho function for Huber                                                                                                     - UTISTAT         - 2023 Dec 13
%   HUwei                   - Computes weight function psi(u)/u for Huber                                                                                         - UTISTAT         - 2023 Dec 13
%   HYPbdp                  - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                       - UTISTAT         - 2023 Dec 13
%   HYPc                    - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)     - UTISTAT         - 2023 Dec 13
%   HYPck                   - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                             - UTISTAT         - 2023 Dec 13
%   HYPeff                  - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                           - UTISTAT         - 2023 Dec 13
%   HYPk                    - Computes breakdown point and efficiency for hyp. tan. estimator                                                                     - UTISTAT         - 2023 Dec 13
%   HYPpsi                  - Computes psi function for hyperbolic tangent estimator                                                                              - UTISTAT         - 2023 Dec 13
%   HYPpsider               - Computes derivative of psi function for hyperbolic tangent estimator                                                                - UTISTAT         - 2023 Dec 13
%   HYPpsix                 - Computes psi function for hyperbolic tangent estimator times x                                                                      - UTISTAT         - 2023 Dec 13
%   HYPrho                  - Computes rho function  using hyperbolic tangent estimator                                                                           - UTISTAT         - 2023 Dec 13
%   HYPwei                  - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                  - UTISTAT         - 2023 Dec 13
%   inversegamcdf           - Computes inverse-gamma cumulative distribution function                                                                             - ProbDist        - 2023 Dec 13
%   inversegaminv           - Computes the inverse of the inverse-gamma cumulative distribution function                                                          - ProbDist        - 2023 Dec 13
%   inversegampdf           - Computes inverse-gamma probability density function                                                                                 - ProbDist        - 2023 Dec 13
%   kdebiv                  - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                      - UTISTAT         - 2023 Dec 13
%   logmvnpdfFS             - Produces log of Multivariate normal probability density function (pdf)                                                              - ProbDist        - 2023 Dec 13
%   mahalCorAna             - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT         - 2023 Dec 13
%   mahalFS                 - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT         - 2023 Dec 13
%   mdpattern               - Finds and plots missing data patterns                                                                                               - UTISTAT         - 2023 Dec 13
%   mdpd                    - Computes Minimum Distance Power Divergence statistics                                                                               - UTISTAT         - 2023 Dec 13
%   medcouple               - Computes the medcouple, a robust skewness estimator                                                                                 - UTISTAT         - 2023 Dec 13
%   Mlocation               - Finds the M estimator of location in a univariate sample                                                                            - UTISTAT         - 2023 Dec 13
%   Mlocsca                 - Computes M estimator of location and scale in univariate samples                                                                    - UTISTAT         - 2023 Dec 13
%   Mscale                  - Finds the M estimator of the scale                                                                                                  - UTISTAT         - 2023 Dec 13
%   msdcutoff               - Mahalanobis Squared Distance cutoff                                                                                                 - UTISTAT         - 2023 Dec 13
%   mtR                     - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                            - ProbDist        - 2023 Dec 13
%   ncpci                   - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                    - UTISTAT         - 2023 Dec 13
%   ncx2mixtcdf             - Computes cumulative distribution function of a linear combination of non-central chi-square +N(0,s)                                 - ProbDist        - 2023 Dec 13
%   normBoxCox              - Computes (normalized) Box-Cox transformation                                                                                        - UTISTAT         - 2023 Dec 13
%   normYJ                  - Computes (normalized) Yeo-Johnson transformation                                                                                    - UTISTAT         - 2023 Dec 13
%   normYJpn                - Computes (normalized) extended Yeo-Johnson transformation                                                                           - UTISTAT         - 2023 Dec 13
%   OPTbdp                  - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT         - 2023 Dec 13
%   OPTc                    - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                         - UTISTAT         - 2023 Dec 13
%   OPTeff                  - Finds the constant c which is associated to the requested efficiency                                                                - UTISTAT         - 2023 Dec 13
%   OptimalCuttingFrequency - Computes the optimal cutting frequency for the Fourier estimator of integrated variance                                             - UTISTAT         - 2023 Dec 13
%   OPTpsi                  - Computes psi function (derivative of rho function) for optimal weight function                                                      - UTISTAT         - 2023 Dec 13
%   OPTpsider               - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                 - UTISTAT         - 2023 Dec 13
%   OPTpsix                 - Computes psi function (derivative of rho function) times x                                                                          - UTISTAT         - 2023 Dec 13
%   OPTrho                  - Computes rho function for optimal weight function                                                                                   - UTISTAT         - 2023 Dec 13
%   OPTwei                  - Computes weight function psi(u)/u for optimal weight function                                                                       - UTISTAT         - 2023 Dec 13
%   PDbdp                   - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                          - UTISTAT         - 2023 Dec 13
%   PDc                     - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator                - UTISTAT         - 2023 Dec 13
%   PDeff                   - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                     - UTISTAT         - 2023 Dec 13
%   PDpsi                   - Computes psi function (derivative of rho function) for minimum density power divergence estimator                                   - UTISTAT         - 2023 Dec 13
%   PDpsider                - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                      - UTISTAT         - 2023 Dec 13
%   PDpsix                  - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                           - UTISTAT         - 2023 Dec 13
%   PDrho                   - Computes rho function for minimum density power divergence estimator                                                                - UTISTAT         - 2023 Dec 13
%   PDwei                   - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                               - UTISTAT         - 2023 Dec 13
%   Powertra                - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                             - UTISTAT         - 2023 Dec 13
%   Qn                      - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                      - UTISTAT         - 2023 Dec 13
%   qqplotFS                - Creates qqplot of studentized residuals with envelopes                                                                              - VIS-Reg         - 2023 Dec 13
%   RandIndexFS             - Calculates Rand type Indices to compare two partitions                                                                              - UTISTAT         - 2023 Dec 13
%   RhoPsiWei               - Finds rho, psi, psi', w functions given bdp, or eff or tuning constant c                                                            - UTISTAT         - 2023 Dec 13
%   RKbdp                   - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                     - UTISTAT         - 2023 Dec 13
%   RKeff                   - Finds the constants c and M which are associated to the requested efficiency and ARP                                                - UTISTAT         - 2023 Dec 13
%   RKpsi                   - Computes psi function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT         - 2023 Dec 13
%   RKpsider                - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                     - UTISTAT         - 2023 Dec 13
%   RKpsix                  - Computes psi function times x for Rocke (translated Tukey's) biweight                                                               - UTISTAT         - 2023 Dec 13
%   RKrho                   - Computes rho function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT         - 2023 Dec 13
%   RKwei                   - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                           - UTISTAT         - 2023 Dec 13
%   rthin                   - Applies independent random thinning to a point pattern                                                                              - UTISTAT         - 2023 Dec 13
%   Sn                      - Robust estimator of scale (robust version of Gini's average difference)                                                             - UTISTAT         - 2023 Dec 13
%   tabulateFS              - Creates frequency table of unique values of x, excluding possible 0 counts                                                          - UTISTAT         - 2023 Dec 13
%   TBbdp                   - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                - UTISTAT         - 2023 Dec 13
%   TBc                     - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                             - UTISTAT         - 2023 Dec 13
%   TBeff                   - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                           - UTISTAT         - 2023 Dec 13
%   TBpsi                   - Computes psi function (derivative of rho function) for Tukey's biweight                                                             - UTISTAT         - 2023 Dec 13
%   TBpsider                - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                        - UTISTAT         - 2023 Dec 13
%   TBpsix                  - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                     - UTISTAT         - 2023 Dec 13
%   TBrho                   - Computes rho function for Tukey's biweight                                                                                          - UTISTAT         - 2023 Dec 13
%   TBwei                   - Computes weight function psi(u)/u for Tukey's biweight                                                                              - UTISTAT         - 2023 Dec 13
%   twdcdf                  - Computes the cumulative distribution function of the Tweedie distribution                                                           - ProbDist        - 2023 Dec 13
%   twdpdf                  - Computes the probability density function of the Tweedie distribution                                                               - ProbDist        - 2023 Dec 13
%   twdrnd                  - Generates random variates from the Tweedie distribution                                                                             - ProbDist        - 2023 Dec 13
%   vervaatrnd              - Simulates random variates from the Vervaat perpetuity distribution                                                                  - ProbDist        - 2023 Dec 13
%   vervaatsim              - Returns a Vervaat perpetuity                                                                                                        - ProbDist        - 2023 Dec 13
%   vervaatxdf              - Returns the pdf and cdf of a Vervaat perpetuity                                                                                     - ProbDist        - 2023 Dec 13
%   winsor                  - Returns a winsorized copy of input                                                                                                  - UTISTAT         - 2023 Dec 13
%   WNChygecdf              - Returns Wallenius' non-central hypergeometric cumulative distribution function                                                      - ProbDist        - 2023 Dec 13
%   WNChygepdf              - Returns Wallenius' non-central hypergeometric probability density function                                                          - ProbDist        - 2023 Dec 13
%   wthin                   - Thins a uni/bi-dimensional dataset                                                                                                  - UTISTAT         - 2023 Dec 13
