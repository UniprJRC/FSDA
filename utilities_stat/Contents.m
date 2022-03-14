% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                    - Description                                                                                                                         - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   basicPower              - Computes the basic power transformation                                                                                             - UTISTAT         - 2021 Feb 01
%   bwe                     - Estimates the bandwidth smoothing parameter for kernel density estimation                                                           - UTISTAT         - 2021 Feb 01
%   CEVmodel                - Computes price and instantaneous variance processes from the CEV model                                                              - UTISTAT         - 2021 Mar 01
%   crosstab2datamatrix     - Recreates the original data matrix X from contingency table N                                                                       - MULT-Categorical- 2021 Feb 01
%   ctsub                   - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                          - UTISTAT         - 2021 Oct 11
%   ellipse                 - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                          - UTISTAT         - 2021 Oct 31
%   exactcdf                - Finds exact p-values                                                                                                                - UTISTAT         - 2021 Feb 01
%   FE_int_vol              - Computes the integrated variance from a diffusion process via the Fourier estimator using Dirichlet kernel                          - UTISTAT         - 2021 Mar 01
%   FE_int_vol_Fejer        - Computes the integrated variance from a diffusion process via the Fourier estimator using Fejer kernel                              - UTISTAT         - 2021 Mar 01
%   FE_spot_vol             - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel                          - UTISTAT         - 2021 Mar 24
%   FE_spot_vol_FFT         - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel, using the FFT algorithm - UTISTAT         - 2021 Mar 24
%   FowlkesMallowsIndex     - Computes the Fowlkes and Mallows index                                                                                              - UTISTAT         - 2021 Feb 01
%   FSMbonfbound            - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                          - UTISTAT         - 2021 Mar 24
%   FSRbonfbound            - Computes Bonferroni bounds for each step of the search (in linear regression)                                                       - UTISTAT         - 2021 Mar 15
%   HAbdp                   - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT         - 2021 Jul 21
%   HAc                     - Computes breakdown point and efficiency associated with constant c                                                                  - UTISTAT         - 2021 Feb 01
%   HAeff                   - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                        - UTISTAT         - 2021 Jul 23
%   HApsi                   - Computes psi function  using Hampel proposal                                                                                        - UTISTAT         - 2021 Feb 01
%   HApsider                - Computes derivative of psi function  using Hampel proposal                                                                          - UTISTAT         - 2021 Feb 01
%   HApsix                  - Computes psi function  using Hampel proposal times x                                                                                - UTISTAT         - 2021 Feb 01
%   HArho                   - Computes rho function  using Hampel proposal                                                                                        - UTISTAT         - 2021 Jul 28
%   HAwei                   - Computes weight function psi(u)/u using Hampel proposal                                                                             - UTISTAT         - 2021 Jul 28
%   HUeff                   - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                   - UTISTAT         - 2021 Feb 01
%   HUpsi                   - Computes psi function (derivative of rho function) for Huber                                                                        - UTISTAT         - 2021 Feb 01
%   HUpsider                - Computes derivative of psi function (second derivative of rho function) for Huber                                                   - UTISTAT         - 2021 Feb 01
%   HUpsix                  - Computes psi function (derivative of rho function) times x for Huber                                                                - UTISTAT         - 2021 Feb 01
%   HUrho                   - Computes rho function for Huber                                                                                                     - UTISTAT         - 2021 Feb 01
%   HUwei                   - Computes weight function psi(u)/u for Huber                                                                                         - UTISTAT         - 2021 Feb 01
%   HYPbdp                  - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                       - UTISTAT         - 2021 Jul 21
%   HYPc                    - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)     - UTISTAT         - 2021 Feb 01
%   HYPck                   - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                             - UTISTAT         - 2021 Jul 21
%   HYPeff                  - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                           - UTISTAT         - 2021 Feb 01
%   HYPk                    - Computes breakdown point and efficiency for hyp. tan. estimator                                                                     - UTISTAT         - 2021 Feb 01
%   HYPpsi                  - Computes psi function for hyperbolic tangent estimator                                                                              - UTISTAT         - 2021 Feb 01
%   HYPpsider               - Computes derivative of psi function for hyperbolic tangent estimator                                                                - UTISTAT         - 2021 Feb 01
%   HYPpsix                 - Computes psi function for hyperbolic tangent estimator times x                                                                      - UTISTAT         - 2021 Feb 01
%   HYPrho                  - Computes rho function  using hyperbolic tangent estimator                                                                           - UTISTAT         - 2021 Jul 23
%   HYPwei                  - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                  - UTISTAT         - 2021 Jul 21
%   inversegamcdf           - Computes inverse-gamma cumulative distribution function                                                                             - UTISTAT         - 2021 Feb 01
%   inversegaminv           - Inversegampdf Inverse-gamma cumulative distribution function                                                                        - UTISTAT         - 2021 Feb 01
%   inversegampdf           - Computes inverse-gamma probability density function                                                                                 - UTISTAT         - 2021 Feb 01
%   kdebiv                  - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                      - UTISTAT         - 2021 Feb 01
%   logmvnpdfFS             - Produces log of Multivariate normal probability density function (pdf)                                                              - UTISTAT         - 2021 Dec 04
%   mahalCorAna             - MahalFS computes Mahalanobis distances (in squared units) for each row of matrix Y                                                  - UTISTAT         - 2021 Feb 01
%   mahalFS                 - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                          - UTISTAT         - 2021 Nov 14
%   mdpattern               - Finds and plots missing data patterns                                                                                               - UTISTAT         - 2022 Jan 12
%   mdpd                    - Computes Minimum Distance Power Divergence statistics                                                                               - UTISTAT         - 2021 Feb 01
%   Mscale                  - Finds the M estimator of the scale                                                                                                  - UTISTAT         - 2021 Aug 13
%   mtR                     - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                            - UTISTAT         - 2021 Feb 01
%   ncpci                   - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                    - UTISTAT         - 2021 Mar 24
%   ncx2mixtcdf             - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                         - UTISTAT         - 2021 Mar 24
%   normBoxCox              - Computes (normalized) Box-Cox transformation                                                                                        - UTISTAT         - 2021 Mar 24
%   normYJ                  - Computes (normalized) Yeo-Johnson transformation                                                                                    - UTISTAT         - 2021 Mar 24
%   normYJpn                - NormYJ computes (normalized) extended Yeo-Johnson transformation                                                                    - UTISTAT         - 2021 Mar 24
%   OPTbdp                  - Finds the constant c associated to the supplied breakdown point                                                                     - UTISTAT         - 2021 Jul 20
%   OPTc                    - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                         - UTISTAT         - 2021 Mar 02
%   OPTeff                  - Finds the constant c which is associated to the requested efficiency                                                                - UTISTAT         - 2021 Jul 23
%   OptimalCuttingFrequency - Computes the optimal cutting frequency for the Fourier estimator of integrated variance                                             - UTISTAT         - 2021 Mar 01
%   OPTpsi                  - Computes psi function (derivative of rho function) for optimal weight function                                                      - UTISTAT         - 2021 Feb 01
%   OPTpsider               - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                 - UTISTAT         - 2021 Feb 01
%   OPTpsix                 - Computes psi function (derivative of rho function) times x                                                                          - UTISTAT         - 2021 Feb 01
%   OPTrho                  - Computes rho function for optimal weight function                                                                                   - UTISTAT         - 2021 Jul 28
%   OPTwei                  - Computes weight function psi(u)/u for optimal weight function                                                                       - UTISTAT         - 2021 Jul 28
%   PDbdp                   - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                          - UTISTAT         - 2021 Feb 01
%   PDc                     - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator                - UTISTAT         - 2021 Feb 01
%   PDeff                   - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                     - UTISTAT         - 2021 Feb 01
%   PDpsi                   - Computes psi function (derivative of rho function) for minimum density power divergence estimator                                   - UTISTAT         - 2021 Feb 01
%   PDpsider                - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                      - UTISTAT         - 2021 Feb 01
%   PDpsix                  - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                           - UTISTAT         - 2021 Feb 01
%   PDrho                   - Computes rho function for minimum density power divergence estimator                                                                - UTISTAT         - 2021 Aug 12
%   PDwei                   - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                               - UTISTAT         - 2021 Aug 12
%   Powertra                - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                             - UTISTAT         - 2021 Feb 01
%   Qn                      - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                      - UTISTAT         - 2021 Mar 27
%   qqplotFS                - Qqplot of studentized residuals with envelopes                                                                                      - VIS-Reg         - 2021 Oct 01
%   RandIndexFS             - Calculates Rand type Indices to compare two partitions                                                                              - UTISTAT         - 2021 Feb 01
%   RKbdp                   - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                     - UTISTAT         - 2021 Feb 01
%   RKeff                   - Finds the constants c and M which are associated to the requested efficiency and ARP                                                - UTISTAT         - 2021 Feb 01
%   RKpsi                   - Computes psi function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT         - 2021 Feb 01
%   RKpsider                - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                     - UTISTAT         - 2021 Feb 01
%   RKpsix                  - Computes psi function times x for Rocke (translated Tukey's) biweight                                                               - UTISTAT         - 2021 Feb 01
%   RKrho                   - Computes rho function for Rocke (translated Tukey's) biweight                                                                       - UTISTAT         - 2021 Feb 01
%   RKwei                   - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                           - UTISTAT         - 2021 Feb 01
%   rthin                   - Applies independent random thinning to a point pattern                                                                              - UTISTAT         - 2021 Feb 01
%   Sn                      - Robust estimator of scale (robust version of Gini's average difference)                                                             - UTISTAT         - 2021 Feb 01
%   tabulateFS              - Creates frequency table of unique values of x, excluding possible 0 counts                                                          - UTISTAT         - 2021 Feb 01
%   TBbdp                   - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                - UTISTAT         - 2021 Jul 20
%   TBc                     - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                             - UTISTAT         - 2021 Feb 01
%   TBeff                   - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                           - UTISTAT         - 2021 Jul 23
%   TBpsi                   - Computes psi function (derivative of rho function) for Tukey's biweight                                                             - UTISTAT         - 2021 Feb 01
%   TBpsider                - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                        - UTISTAT         - 2021 Feb 01
%   TBpsix                  - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                     - UTISTAT         - 2021 Feb 01
%   TBrho                   - Computes rho function for Tukey's biweight                                                                                          - UTISTAT         - 2021 Jul 28
%   TBwei                   - Computes weight function psi(u)/u for Tukey's biweight                                                                              - UTISTAT         - 2021 Dec 04
%   twdcdf                  - Computes the cumulative distribution function of the Tweedie distribution                                                           - UTISTAT         - 2021 Dec 04
%   twdpdf                  - Twopdf computes the probability density function of the Tweedie distribution                                                        - UTISTAT         - 2021 Mar 01
%   twdrnd                  - Generates random variates from the Tweedie distribution                                                                             - UTISTAT         - 2021 Feb 01
%   vervaatrnd              - Simulates random variates from the Vervaat perpetuity distribution                                                                  - UTISTAT         - 2021 Feb 01
%   vervaatsim              - Returns a Vervaat perpetuity                                                                                                        - UTISTAT         - 2021 Feb 01
%   vervaatxdf              - Returns the pdf and cdf of a Vervaat perpetuity                                                                                     - UTISTAT         - 2021 Feb 01
%   winsor                  - Returns a winsorized copy of input                                                                                                  - UTISTAT         - 2021 Feb 01
%   WNChygepdf              - Returns Wallenius' non-central hypergeometric probability density values                                                            - UTISTAT         - 2021 Feb 01
%   wthin                   - Thins a uni/bi-dimensional dataset                                                                                                  - UTISTAT         - 2021 Mar 24
