% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name                    - Description                                                                                                                                                     - Category        - Date last modified
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   add2boxplot             - Adds labels to the boxplot figure                                                                                                                               - VIS-Mult        - 2024 Aug 18
%   ASbdp                   - Finds the constant c associated to the supplied breakdown point for Andrew's sine function                                                                      - UTISTAT         - 2024 Feb 05
%   ASc                     - Computes breakdown point and efficiency associated with constant c for Andrew's rho function                                                                    - UTISTAT         - 2024 Feb 05
%   ASeff                   - Finds the constant c which is associated to the requested efficiency for Andrew's sine function                                                                 - UTISTAT         - 2024 Feb 05
%   ASpsi                   - Computes psi function (derivative of rho function) for Andrew's sine function                                                                                   - UTISTAT         - 2024 Feb 05
%   ASpsider                - Computes derivative of psi function (second derivative of rho function) for Andrew's sine function                                                              - UTISTAT         - 2024 Feb 05
%   ASpsix                  - Computes psi function (derivative of rho function) times x for Andrew's sine function                                                                           - UTISTAT         - 2024 Feb 05
%   ASrho                   - Computes rho function for Andrew's sine function                                                                                                                - UTISTAT         - 2024 Feb 05
%   ASwei                   - Computes weight function psi(u)/u for Andrew's sine function                                                                                                    - UTISTAT         - 2024 Feb 05
%   basicPower              - Computes the basic power transformation                                                                                                                         - UTISTAT         - 2024 Feb 05
%   bwe                     - Estimates the bandwidth smoothing parameter for kernel density estimation                                                                                       - UTISTAT         - 2024 Apr 03
%   CEVmodel                - Computes price and instantaneous variance processes from the CEV model                                                                                          - UTISTAT         - 2024 Feb 05
%   corrcdf                 - Computes correlation coefficient probability distribution function                                                                                              - ProbDist        - 2024 Feb 05
%   corrpdf                 - Computes correlation coefficient probability density function                                                                                                   - ProbDist        - 2024 Feb 05
%   crosstab2datamatrix     - Recreates the original data matrix X from contingency table N                                                                                                   - MULT-Categorical- 2024 Feb 05
%   ctsub                   - Computes numerical integration from x(1) to z(i) of y=f(x)                                                                                                      - UTISTAT         - 2024 Feb 05
%   ellipse                 - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                                                      - UTISTAT         - 2024 Feb 05
%   exactcdf                - Finds exact p-values                                                                                                                                            - UTISTAT         - 2024 Feb 05
%   FE_int_vol              - Computes the integrated variance from a diffusion process via the Fourier estimator using Dirichlet kernel                                                      - FMvol           - 2024 Feb 05
%   FE_int_vol_Fejer        - Computes the integrated variance from a diffusion process via the Fourier estimator using Fejer kernel                                                          - FMvol           - 2024 Feb 05
%   FE_spot_vol             - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel                                                      - FMvol           - 2025 Jan 30
%   FE_spot_vol_FFT         - Estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel, using the FFT algorithm                             - FMvol           - 2024 Feb 05
%   FM_cov_matrix           - Computes the integrated variance-covariance matrix of any number of processes from discrete observations, via the Fourier-Malliavin estimator with Fejer kernel - FMvol           - 2025 Mar 25
%   FM_int_cov              - Computes the integrated covariance of a bivariate diffusion process via the Fourier-Malliavin estimator                                                         - FMvol           - 2024 Feb 05
%   FM_int_lev              - Computes the integrated leverage of a diffusion process via the Fourier-Malliavin estimator                                                                     - FMvol           - 2024 Feb 05
%   FM_int_quart            - Computes the integrated quarticity of a diffusion process via the Fourier-Malliavin estimator                                                                   - FMvol           - 2024 Feb 05
%   FM_int_vol              - Computes the integrated variance of a diffusion process via the Fourier-Malliavin estimator                                                                     - FMvol           - 2024 Feb 05
%   FM_int_volvol           - Computes the integrated volatility of volatility of a diffusion process via the Fourier-Malliavin estimator                                                     - FMvol           - 2024 Feb 05
%   FM_spot_cov             - Computes the spot covariance of a bivariate diffusion process via the Fourier-Malliavin estimator                                                               - FMvol           - 2024 Feb 05
%   FM_spot_lev             - Computes the spot leverage of a diffusion process via the Fourier-Malliavin estimator                                                                           - FMvol           - 2024 Feb 05
%   FM_spot_quart           - Computes the spot quarticity of a diffusion process via the Fourier-Malliavin estimator                                                                         - FMvol           - 2024 Feb 05
%   FM_spot_vol             - Computes the spot volatility of a diffusion process via the Fourier-Malliavin estimator                                                                         - FMvol           - 2024 Feb 05
%   FM_spot_volvol          - Computes the spot volatiity of volatility of a diffusion process via the Fourier-Malliavin estimator                                                            - FMvol           - 2024 Feb 05
%   FNChygecdf              - Returns Fisher non-central hypergeometric cumulative distribution function                                                                                      - ProbDist        - 2024 Apr 04
%   FNChygeinv              - Computes the inverse of the Fisher non central hypergeometric cumulative distribution function (cdf)                                                            - ProbDist        - 2024 Apr 19
%   FNChygepdf              - Returns Fisher non-central hypergeometric probability density function                                                                                          - ProbDist        - 2024 Apr 16
%   FNChygernd              - Random arrays from the Fisher non central hypergeometric distribution                                                                                           - ProbDist        - 2024 Apr 19
%   FowlkesMallowsIndex     - Computes the Fowlkes and Mallows index                                                                                                                          - UTISTAT         - 2024 Feb 05
%   FSMbonfbound            - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                                                      - UTISTAT         - 2024 Feb 05
%   FSRbonfbound            - Computes Bonferroni bounds for each step of the search (in linear regression)                                                                                   - UTISTAT         - 2024 Feb 05
%   genr8next               - Returns a vector of pseudorandom number sequence                                                                                                                - UTISTAT         - 2024 Apr 03
%   grpstatsFS              - Calls grpstats and reshapes the output in a much better way                                                                                                     - UTISTAT         - 2025 Mar 25
%   HAbdp                   - Finds the constant c associated to the supplied breakdown point                                                                                                 - UTISTAT         - 2024 Feb 05
%   HAc                     - Computes breakdown point and efficiency associated with constant c                                                                                              - UTISTAT         - 2024 Feb 05
%   HAeff                   - Finds the tuning constant that guarrantees a requested asymptotic efficiency                                                                                    - UTISTAT         - 2024 Feb 05
%   HApsi                   - Computes psi function  using Hampel proposal                                                                                                                    - UTISTAT         - 2024 Feb 05
%   HApsider                - Computes derivative of psi function  using Hampel proposal                                                                                                      - UTISTAT         - 2024 Feb 05
%   HApsix                  - Computes psi function  using Hampel proposal times x                                                                                                            - UTISTAT         - 2024 Feb 05
%   HArho                   - Computes rho function  using Hampel proposal                                                                                                                    - UTISTAT         - 2024 Feb 05
%   HAwei                   - Computes weight function psi(u)/u using Hampel proposal                                                                                                         - UTISTAT         - 2024 Feb 05
%   Heston1D                - Simulates observations and instantaneous variances from the Heston model                                                                                        - FMvol           - 2024 Feb 05
%   Heston2D                - Simulates observations and instantaneous variances from the bivariate Heston model                                                                              - FMvol           - 2024 Feb 05
%   HUc                     - Computes breakdown point and efficiency associated with constant c for Huber link                                                                               - UTISTAT         - 2024 Feb 05
%   HUeff                   - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                                               - UTISTAT         - 2024 Feb 05
%   HUpsi                   - Computes psi function (derivative of rho function) for Huber                                                                                                    - UTISTAT         - 2024 Feb 05
%   HUpsider                - Computes derivative of psi function (second derivative of rho function) for Huber                                                                               - UTISTAT         - 2024 Feb 05
%   HUpsix                  - Computes psi function (derivative of rho function) times x for Huber                                                                                            - UTISTAT         - 2024 Feb 05
%   HUrho                   - Computes rho function for Huber                                                                                                                                 - UTISTAT         - 2024 Feb 05
%   HUwei                   - Computes weight function psi(u)/u for Huber                                                                                                                     - UTISTAT         - 2024 Feb 05
%   HYPbdp                  - Finds constant c which is associated to the requested breakdown point for hyp. tan. estimator                                                                   - UTISTAT         - 2024 Feb 05
%   HYPc                    - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)                                 - UTISTAT         - 2024 Feb 05
%   HYPck                   - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                                                         - UTISTAT         - 2024 Feb 05
%   HYPeff                  - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                                                       - UTISTAT         - 2024 Feb 05
%   HYPk                    - Computes breakdown point and efficiency for hyp. tan. estimator                                                                                                 - UTISTAT         - 2024 Feb 05
%   HYPpsi                  - Computes psi function for hyperbolic tangent estimator                                                                                                          - UTISTAT         - 2024 Feb 05
%   HYPpsider               - Computes derivative of psi function for hyperbolic tangent estimator                                                                                            - UTISTAT         - 2024 Feb 05
%   HYPpsix                 - Computes psi function for hyperbolic tangent estimator times x                                                                                                  - UTISTAT         - 2024 Feb 05
%   HYPrho                  - Computes rho function  using hyperbolic tangent estimator                                                                                                       - UTISTAT         - 2024 Feb 05
%   HYPwei                  - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                                              - UTISTAT         - 2024 Feb 05
%   inversegamcdf           - Computes inverse-gamma cumulative distribution function                                                                                                         - ProbDist        - 2024 Feb 05
%   inversegaminv           - Computes the inverse of the inverse-gamma cumulative distribution function                                                                                      - ProbDist        - 2024 Feb 05
%   inversegampdf           - Computes inverse-gamma probability density function                                                                                                             - ProbDist        - 2024 Feb 05
%   kdebiv                  - Computes (and optionally plots) a kernel smoothing estimate for bivariate data                                                                                  - UTISTAT         - 2024 Feb 05
%   logfactorial            - Returns the logarithm of the factorial                                                                                                                          - UTISTAT         - 2024 Apr 19
%   logmvnpdfFS             - Produces log of Multivariate normal probability density function (pdf)                                                                                          - ProbDist        - 2024 Feb 05
%   mahalCorAna             - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                                                      - UTISTAT         - 2024 Feb 05
%   mahalFS                 - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                                                      - UTISTAT         - 2024 Feb 05
%   mdpattern               - Finds and plots missing data patterns                                                                                                                           - UTISTAT         - 2024 Feb 05
%   mdpd                    - Computes Minimum Distance Power Divergence statistics                                                                                                           - UTISTAT         - 2024 Feb 05
%   mFNChygepdf             - Returns Fisher multivariate non-central hypergeometric probability density function                                                                             - ProbDist        - 2024 Apr 19
%   mFNChygernd             - Returns Fisher multivariate non-central hypergeometric random variate generation                                                                                - ProbDist        - 2025 Jan 12
%   Mlocation               - Finds the M estimator of location in a univariate sample                                                                                                        - UTISTAT         - 2024 Feb 05
%   Mlocsca                 - Computes M estimator of location and scale in univariate samples                                                                                                - UTISTAT         - 2024 Feb 05
%   Mscale                  - Finds the M estimator of the scale                                                                                                                              - UTISTAT         - 2024 Feb 05
%   msdcutoff               - Mahalanobis Squared Distance cutoff                                                                                                                             - UTISTAT         - 2024 Mar 04
%   mtR                     - Generates the same random numbers produced by R software with Mersenne Twister mt19937ar                                                                        - ProbDist        - 2024 Feb 05
%   mWNChygepdf             - Returns Wallenius' multivariate non-central hypergeometric probability density function                                                                         - ProbDist        - 2024 Apr 16
%   mWNChygernd             - Returns random arrays from the Wallenius non central hypergeometric distribution                                                                                - ProbDist        - 2024 Apr 16
%   ncpci                   - Non centrality parameter confidence interval (taken from effect_of_size_toolbox)                                                                                - UTISTAT         - 2024 Feb 05
%   ncx2mixtcdf             - Computes cumulative distribution function of a linear combination of non-central chi-square +N(0,s)                                                             - ProbDist        - 2024 Feb 05
%   normBoxCox              - Computes (normalized) Box-Cox transformation                                                                                                                    - UTISTAT         - 2024 Jun 06
%   normYJ                  - Computes (normalized) Yeo-Johnson transformation                                                                                                                - UTISTAT         - 2024 Jun 03
%   normYJpn                - Computes (normalized) extended Yeo-Johnson transformation                                                                                                       - UTISTAT         - 2024 Feb 05
%   OPTbdp                  - Finds the constant c associated to the supplied breakdown point                                                                                                 - UTISTAT         - 2024 Feb 05
%   OPTc                    - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                                                     - UTISTAT         - 2024 Feb 05
%   OPTeff                  - Finds the constant c which is associated to the requested efficiency                                                                                            - UTISTAT         - 2024 Feb 05
%   OptimalCuttingFrequency - Computes the optimal cutting frequency for the Fourier estimator of integrated variance                                                                         - FMvol           - 2025 Mar 25
%   OPTpsi                  - Computes psi function (derivative of rho function) for optimal weight function                                                                                  - UTISTAT         - 2024 Feb 05
%   OPTpsider               - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                                             - UTISTAT         - 2024 Feb 05
%   OPTpsix                 - Computes psi function (derivative of rho function) times x                                                                                                      - UTISTAT         - 2024 Feb 05
%   OPTrho                  - Computes rho function for optimal weight function                                                                                                               - UTISTAT         - 2024 Feb 05
%   OPTwei                  - Computes weight function psi(u)/u for optimal weight function                                                                                                   - UTISTAT         - 2024 Feb 05
%   OSILA                   - Finds the k-th order statistic                                                                                                                                  - UTISTAT         - 2024 Apr 18
%   PDbdp                   - Finds the constant alpha associated to the supplied breakdown point for minimum power divergence estimator                                                      - UTISTAT         - 2024 Feb 05
%   PDc                     - Computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator                                            - UTISTAT         - 2024 Feb 05
%   PDeff                   - Finds the constant alpha which is associated to the requested efficiency for minimum power divergence estimator                                                 - UTISTAT         - 2024 Feb 05
%   PDpsi                   - Computes psi function (derivative of rho function) for minimum density power divergence estimator                                                               - UTISTAT         - 2024 Feb 05
%   PDpsider                - Computes derivative of psi function (second derivative of rho function) for minimum power divergence estimator                                                  - UTISTAT         - 2024 Feb 05
%   PDpsix                  - Computes psi function (derivative of rho function) times x for minimum density power divergence estimator                                                       - UTISTAT         - 2024 Feb 05
%   PDrho                   - Computes rho function for minimum density power divergence estimator                                                                                            - UTISTAT         - 2024 Feb 05
%   PDwei                   - Computes weight function psi(u)/u for  for minimum density power divergence estimator                                                                           - UTISTAT         - 2024 Feb 05
%   pivotCoord              - Transforms into pivot coordinates as a special case of isometric logratio coordinates                                                                           - UTISTAT         - 2024 Aug 04
%   pivotCoordInv           - Transforms back the output of pivotCoord                                                                                                                        - UTISTAT         - 2024 Apr 03
%   Powertra                - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                                                         - UTISTAT         - 2024 Feb 05
%   Qn                      - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                                                  - UTISTAT         - 2024 Feb 05
%   qqplotFS                - Creates qqplot of studentized residuals with envelopes                                                                                                          - VIS-Reg         - 2024 Feb 05
%   RandIndexFS             - Calculates Rand type Indices to compare two partitions                                                                                                          - UTISTAT         - 2024 Nov 04
%   RhoPsiWei               - Finds rho, psi, psi', w functions given bdp, or eff or tuning constant c                                                                                        - UTISTAT         - 2024 Feb 05
%   RKbdp                   - Finds the constants c associated to the supplied breakdown point and asymptotic rejection point                                                                 - UTISTAT         - 2024 Feb 05
%   RKeff                   - Finds the constants c and M which are associated to the requested efficiency and ARP                                                                            - UTISTAT         - 2024 Feb 05
%   RKpsi                   - Computes psi function for Rocke (translated Tukey's) biweight                                                                                                   - UTISTAT         - 2024 Feb 05
%   RKpsider                - Computes derivative of psi function (second derivative of rho function) for Rocke (translated Tukey's) biweight                                                 - UTISTAT         - 2024 Feb 05
%   RKpsix                  - Computes psi function times x for Rocke (translated Tukey's) biweight                                                                                           - UTISTAT         - 2024 Feb 05
%   RKrho                   - Computes rho function for Rocke (translated Tukey's) biweight                                                                                                   - UTISTAT         - 2024 Feb 05
%   RKwei                   - Computes weight function psi(u)/u for Rocke (translated Tukey's) biweight                                                                                       - UTISTAT         - 2024 Feb 05
%   rthin                   - Applies independent random thinning to a point pattern                                                                                                          - UTISTAT         - 2024 Feb 05
%   Sn                      - Robust estimator of scale (robust version of Gini's average difference)                                                                                         - UTISTAT         - 2024 Feb 05
%   tabulateFS              - Creates frequency table of unique values of x, excluding possible 0 counts                                                                                      - UTISTAT         - 2024 Feb 05
%   TBbdp                   - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                                            - UTISTAT         - 2024 Feb 05
%   TBc                     - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                                                         - UTISTAT         - 2024 Feb 05
%   TBeff                   - Finds the constant c which is associated to the requested efficiency for Tukey's biweight                                                                       - UTISTAT         - 2024 Feb 05
%   TBpsi                   - Computes psi function (derivative of rho function) for Tukey's biweight                                                                                         - UTISTAT         - 2024 Feb 05
%   TBpsider                - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                                                    - UTISTAT         - 2024 Feb 05
%   TBpsix                  - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                                                 - UTISTAT         - 2024 Feb 05
%   TBrho                   - Computes rho function for Tukey's biweight                                                                                                                      - UTISTAT         - 2024 Feb 05
%   TBwei                   - Computes weight function psi(u)/u for Tukey's biweight                                                                                                          - UTISTAT         - 2024 Feb 05
%   tobitcdf                - Returns cumulative distribution function from the tobit model                                                                                                   - ProbDist        - 2024 May 07
%   tobitinv                - Computes the inverse of the tobit cumulative distribution function                                                                                              - ProbDist        - 2024 May 07
%   tobitpdf                - Returns probability density function from the tobit model                                                                                                       - ProbDist        - 2024 May 07
%   tobitrnd                - Random arrays from the tobit distribution                                                                                                                       - ProbDist        - 2024 May 13
%   twdcdf                  - Computes the cumulative distribution function of the Tweedie distribution                                                                                       - ProbDist        - 2024 Feb 05
%   twdpdf                  - Computes the probability density function of the Tweedie distribution                                                                                           - ProbDist        - 2024 Feb 05
%   twdrnd                  - Generates random variates from the Tweedie distribution                                                                                                         - ProbDist        - 2024 Feb 05
%   vervaatrnd              - Simulates random variates from the Vervaat perpetuity distribution                                                                                              - ProbDist        - 2024 Feb 05
%   vervaatsim              - Returns a Vervaat perpetuity                                                                                                                                    - ProbDist        - 2024 Feb 05
%   vervaatxdf              - Returns the pdf and cdf of a Vervaat perpetuity                                                                                                                 - ProbDist        - 2024 Feb 05
%   winsor                  - Returns a winsorized copy of input                                                                                                                              - UTISTAT         - 2024 Feb 05
%   WNChygecdf              - Returns Wallenius' non-central hypergeometric cumulative distribution function                                                                                  - ProbDist        - 2024 Apr 04
%   WNChygeinv              - Computes the inverse of the Wallenius non central hypergeometric cumulative distribution function (cdf)                                                         - ProbDist        - 2024 Apr 19
%   WNChygepdf              - Returns Wallenius' non-central hypergeometric probability density function                                                                                      - ProbDist        - 2024 Apr 16
%   WNChygernd              - Random arrays from the Wallenius non central hypergeometric distribution                                                                                        - ProbDist        - 2024 Apr 19
%   wthin                   - Thins a uni/bi-dimensional dataset                                                                                                                              - UTISTAT         - 2024 Apr 03
