% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                    - Description                                                                                                                         - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   add2boxplot             - Adds labels to the boxplot figure                                                                                                   - VIS-Mult        - 2022 Dec 03
%   ASbdp                   - Finds the constant c associated to the supplied breakdown point for Andrew's sine function                                          - UTISTAT         - 2023 May 29
%   ASc                     - Computes breakdown point and efficiency associated with constant c for Andrew's rho function                                        - UTISTAT         - 2023 May 28
%   ASeff                   - Finds the constant c which is associated to the requested efficiency for Andrew's sine function                                     - UTISTAT         - 2023 May 28
%   ASpsi                   - Computes psi function (derivative of rho function) for Andrew's sine function                                                       - UTISTAT         - 2023 May 28
%   ASpsider                - Computes derivative of psi function (second derivative of rho function) for Andrew's sine function                                  - UTISTAT         - 2023 May 29
%   ASpsix                  - Computes psi function (derivative of rho function) times x for Andrew's sine function                                               - UTISTAT         - 2023 May 28
%   ASrho                   - Computes rho function for Andrew's sine function                                                                                    - UTISTAT         - 2023 May 28
%   ASwei                   - Computes weight function psi(u)/u for Andrew's sine function                                                                        - UTISTAT         - 2023 May 28
%   basicPower              - Computes the basic power transformation                                                                                             - UTISTAT         - 2023 Apr 19
%   bwe                     - Estimates the bandwidth smoothing parameter for kernel density estimation                                                           - UTISTAT         - 2023 Apr 06
%   CEVmodel                - Computes price and instantaneous variance processes from the CEV model                                                              - UTISTAT         - 2022 Nov 10
%   corrcdf                 - Computes correlation coefficient probability distribution function                                                                  - ProbDist        - 2022 Dec 02
%   corrpdf                 - Computes correlation coefficient probability density function                                                                       - ProbDist        - 2022 Dec 02
%   crosstab2datamatrix     - Recreates the original data matrix X from contingency table N                                                                       - MULT-Categorical- 2023 Apr 19
%   ctsub                   - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                          - UTISTAT         - 2023 Apr 06
%   ellipse                 - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                          - UTISTAT         - 2022 Nov 10
%   exactcdf                - Finds exact p-values                                                                                                                - UTISTAT         - 2022 Nov 10
%   FE_int_vol              - Computes the integrated variance from a diffusion process via the Fourier estimator using Dirichlet kernel                          - UTISTAT         - 2023 Apr 19
%   FE_int_vol_Fejer        - Computes the integrated variance from a diffusion process via the Fourier estimator using Fejer kernel                              - UTISTAT         - 2023 Apr 19
%   FE_spot_vol             - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel                          - UTISTAT         - 2023 Apr 19
%   FE_spot_vol_FFT         - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel, using the FFT algorithm - UTISTAT         - 2023 Apr 19
%   FowlkesMallowsIndex     - Computes the Fowlkes and Mallows index                                                                                              - UTISTAT         - 2022 Nov 10
%   FSMbonfbound            - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                          - UTISTAT         - 2022 Nov 10
%   FSRbonfbound            - Computes Bonferroni bounds for each step of the search (in linear regression)                                                       - UTISTAT         - 2023 Apr 06
%   HAbdp                   - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT         - 2022 Nov 10
%   HAc                     - Computes breakdown point and efficiency associated with constant c                                                                  - UTISTAT         - 2023 Apr 19
%   HAeff                   - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                        - UTISTAT         - 2022 Nov 10
%   HApsi                   - Computes psi function  using Hampel proposal                                                                                        - UTISTAT         - 2022 Nov 10
%   HApsider                - Computes derivative of psi function  using Hampel proposal                                                                          - UTISTAT         - 2022 Nov 10
%   HApsix                  - Computes psi function  using Hampel proposal times x                                                                                - UTISTAT         - 2022 Nov 10
%   HArho                   - Computes rho function  using Hampel proposal                                                                                        - UTISTAT         - 2022 Nov 10
%   HAwei                   - Computes weight function psi(u)/u using Hampel proposal                                                                             - UTISTAT         - 2022 Nov 10
%   HUc                     - Computes breakdown point and efficiency associated with constant c for Huber link                                                   - UTISTAT         - 2023 May 28
%   HUeff                   - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                   - UTISTAT         - 2023 May 28
%   HUpsi                   - Computes psi function (derivative of rho function) for Huber                                                                        - UTISTAT         - 2022 Nov 10
%   HUpsider                - Computes derivative of psi function (second derivative of rho function) for Huber                                                   - UTISTAT         - 2022 Nov 10
%   HUpsix                  - Computes psi function (derivative of rho function) times x for Huber                                                                - UTISTAT         - 2022 Nov 10
%   HUrho                   - Computes rho function for Huber                                                                                                     - UTISTAT         - 2022 Nov 10
%   HUwei                   - Computes weight function psi(u)/u for Huber                                                                                         - UTISTAT         - 2023 Apr 06
%   HYPbdp                  - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                       - UTISTAT         - 2023 May 28
%   HYPc                    - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)     - UTISTAT         - 2023 Apr 19
%   HYPck                   - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                             - UTISTAT         - 2022 Nov 10
%   HYPeff                  - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                           - UTISTAT         - 2022 Nov 10
%   HYPk                    - Computes breakdown point and efficiency for hyp. tan. estimator                                                                     - UTISTAT         - 2023 Apr 19
%   HYPpsi                  - Computes psi function for hyperbolic tangent estimator                                                                              - UTISTAT         - 2022 Nov 10
%   HYPpsider               - Computes derivative of psi function for hyperbolic tangent estimator                                                                - UTISTAT         - 2022 Nov 10
%   HYPpsix                 - Computes psi function for hyperbolic tangent estimator times x                                                                      - UTISTAT         - 2022 Nov 10
%   HYPrho                  - Computes rho function  using hyperbolic tangent estimator                                                                           - UTISTAT         - 2022 Nov 10
%   HYPwei                  - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                  - UTISTAT         - 2022 Nov 10
%   inversegamcdf           - Computes inverse-gamma cumulative distribution function                                                                             - ProbDist        - 2022 Dec 02
%   inversegaminv           - Computes the inverse of the inverse-gamma cumulative distribution function                                                          - ProbDist        - 2022 Dec 02
%   inversegampdf           - Computes inverse-gamma probability density function                                                                                 - ProbDist        - 2022 Dec 02
%   kdebiv                  - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                      - UTISTAT         - 2022 Nov 10
%   logmvnpdfFS             - Produces log of Multivariate normal probability density function (pdf)                                                              - ProbDist        - 2022 Dec 02
%   mahalCorAna             - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT         - 2022 Dec 03
%   mahalFS                 - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT         - 2022 Nov 10
%   mdpattern               - Finds and plots missing data patterns                                                                                               - UTISTAT         - 2023 Apr 19
%   mdpd                    - Computes Minimum Distance Power Divergence statistics                                                                               - UTISTAT         - 2023 Apr 19
%   medcouple               - Computes the medcouple, a robust skewness estimator                                                                                 - UTISTAT         - 2022 Nov 10
%   Mlocation               - Finds the M estimator of location in a univariate sample                                                                            - UTISTAT         - 2022 Nov 10
%   Mlocsca                 - Computes M estimator of location and scale in univariate samples                                                                    - UTISTAT         - 2022 Nov 10
%   Mscale                  - Finds the M estimator of the scale                                                                                                  - UTISTAT         - 2023 May 29
%   msdcutoff               - Mahalanobis Squared Distance cutoff                                                                                                 - UTISTAT         - 2023 May 10
%   mtR                     - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                            - ProbDist        - 2022 Dec 02
%   ncpci                   - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                    - UTISTAT         - 2023 Apr 19
%   ncx2mixtcdf             - Computes cumulative distribution function of a linear combination of non-central chi-square +N(0,s)                                 - ProbDist        - 2023 Apr 19
%   normBoxCox              - Computes (normalized) Box-Cox transformation                                                                                        - UTISTAT         - 2023 Apr 19
%   normYJ                  - Computes (normalized) Yeo-Johnson transformation                                                                                    - UTISTAT         - 2023 Apr 19
%   normYJpn                - Computes (normalized) extended Yeo-Johnson transformation                                                                           - UTISTAT         - 2023 Apr 19
%   OPTbdp                  - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT         - 2022 Nov 10
%   OPTc                    - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                         - UTISTAT         - 2022 Nov 10
%   OPTeff                  - Finds the constant c which is associated to the requested efficiency                                                                - UTISTAT         - 2023 May 28
%   OptimalCuttingFrequency - Computes the optimal cutting frequency for the Fourier estimator of integrated variance                                             - UTISTAT         - 2023 Apr 03
%   OPTpsi                  - Computes psi function (derivative of rho function) for optimal weight function                                                      - UTISTAT         - 2022 Nov 10
%   OPTpsider               - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                 - UTISTAT         - 2022 Nov 10
%   OPTpsix                 - Computes psi function (derivative of rho function) times x                                                                          - UTISTAT         - 2022 Nov 10
%   OPTrho                  - Computes rho function for optimal weight function                                                                                   - UTISTAT         - 2022 Nov 10
%   OPTwei                  - Computes weight function psi(u)/u for optimal weight function                                                                       - UTISTAT         - 2022 Nov 10
%   PDbdp                   - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                          - UTISTAT         - 2023 Apr 07
%   PDc                     - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator                - UTISTAT         - 2023 May 28
%   PDeff                   - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                     - UTISTAT         - 2023 May 28
%   PDpsi                   - Computes psi function (derivative of rho function) for minimum density power divergence estimator                                   - UTISTAT         - 2023 Apr 07
%   PDpsider                - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                      - UTISTAT         - 2023 Apr 07
%   PDpsix                  - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                           - UTISTAT         - 2023 Apr 07
%   PDrho                   - Computes rho function for minimum density power divergence estimator                                                                - UTISTAT         - 2023 Apr 07
%   PDwei                   - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                               - UTISTAT         - 2023 Apr 07
%   Powertra                - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                             - UTISTAT         - 2023 Apr 19
%   Qn                      - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                      - UTISTAT         - 2023 Apr 06
%   qqplotFS                - Creates qqplot of studentized residuals with envelopes                                                                              - VIS-Reg         - 2023 Apr 19
%   RandIndexFS             - Calculates Rand type Indices to compare two partitions                                                                              - UTISTAT         - 2022 Nov 10
%   RhoPsiWei               - Finds rho, psi, psi', w functions given bdp, or eff or tuning constant c                                                            - UTISTAT         - 2023 Apr 19
%   RKbdp                   - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                     - UTISTAT         - 2023 Apr 06
%   RKeff                   - Finds the constants c and M which are associated to the requested efficiency and ARP                                                - UTISTAT         - 2022 Nov 10
%   RKpsi                   - Computes psi function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT         - 2023 Apr 06
%   RKpsider                - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                     - UTISTAT         - 2023 Apr 06
%   RKpsix                  - Computes psi function times x for Rocke (translated Tukey's) biweight                                                               - UTISTAT         - 2023 Apr 06
%   RKrho                   - Computes rho function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT         - 2023 Apr 06
%   RKwei                   - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                           - UTISTAT         - 2023 Apr 06
%   rthin                   - Applies independent random thinning to a point pattern                                                                              - UTISTAT         - 2022 Nov 10
%   Sn                      - Robust estimator of scale (robust version of Gini's average difference)                                                             - UTISTAT         - 2023 Apr 06
%   tabulateFS              - Creates frequency table of unique values of x, excluding possible 0 counts                                                          - UTISTAT         - 2022 Nov 10
%   TBbdp                   - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                - UTISTAT         - 2023 May 27
%   TBc                     - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                             - UTISTAT         - 2023 May 28
%   TBeff                   - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                           - UTISTAT         - 2023 May 28
%   TBpsi                   - Computes psi function (derivative of rho function) for Tukey's biweight                                                             - UTISTAT         - 2023 Apr 06
%   TBpsider                - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                        - UTISTAT         - 2022 Nov 10
%   TBpsix                  - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                     - UTISTAT         - 2022 Nov 10
%   TBrho                   - Computes rho function for Tukey's biweight                                                                                          - UTISTAT         - 2023 Apr 06
%   TBwei                   - Computes weight function psi(u)/u for Tukey's biweight                                                                              - UTISTAT         - 2023 Apr 06
%   twdcdf                  - Computes the cumulative distribution function of the Tweedie distribution                                                           - ProbDist        - 2022 Dec 02
%   twdpdf                  - Computes the probability density function of the Tweedie distribution                                                               - ProbDist        - 2022 Dec 02
%   twdrnd                  - Generates random variates from the Tweedie distribution                                                                             - ProbDist        - 2022 Dec 02
%   vervaatrnd              - Simulates random variates from the Vervaat perpetuity distribution                                                                  - ProbDist        - 2022 Dec 02
%   vervaatsim              - Returns a Vervaat perpetuity                                                                                                        - ProbDist        - 2022 Dec 02
%   vervaatxdf              - Returns the pdf and cdf of a Vervaat perpetuity                                                                                     - ProbDist        - 2022 Dec 02
%   winsor                  - Returns a winsorized copy of input                                                                                                  - UTISTAT         - 2022 Nov 10
%   WNChygecdf              - Returns Wallenius' non-central hypergeometric cumulative distribution function                                                      - ProbDist        - 2023 Apr 06
%   WNChygepdf              - Returns Wallenius' non-central hypergeometric probability density function                                                          - ProbDist        - 2023 Apr 06
%   wthin                   - Thins a uni/bi-dimensional dataset                                                                                                  - UTISTAT         - 2022 Nov 10
