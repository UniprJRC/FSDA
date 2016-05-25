% UTILITIES_STAT
%
% File names, description, category and date last modified
%
%   Name          - Description                                                                                                                                                         - Category- Date last modified
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   basicPower    - Computes the basic power transformation                                                                                                                             - :UTISTAT- 2016 May 17
%   eigs_sigma    - Sets a non zero value for the optional parameter sigma of function eigs                                                                                             -         - 2016 May 17
%   ellipse       - Generates an ellipse given mu (location vector) and Sigma (scatter matrix)                                                                                          - :UTISTAT- 2016 May 25
%   FSMbonfbound  - Computes Bonferroni bounds for each step of the  search (in mult analysis)                                                                                          - :UTISTAT- 2016 May 09
%   FSRbonfbound  - Computes Bonferroni bounds for each step of the search (in linear regression)                                                                                       - :UTISTAT- 2016 May 09
%   HAbdp         - Finds the constant c associated to the supplied breakdown point                                                                                                     - :UTISTAT- 2016 May 09
%   HAc           - Computes breakdown point and efficiency associated with constant c                                                                                                  - :UTISTAT- 2016 May 25
%   HAeff         - Finds the tuning constant guarrantees a requested asymptotic efficiency                                                                                             - :UTISTAT- 2016 May 25
%   HApsi         - Computes psi function  using Hampel proposal                                                                                                                        - :UTISTAT- 2016 May 25
%   HApsider      - Computes derivative of psi function  using Hampel proposal                                                                                                          - :UTISTAT- 2016 May 25
%   HApsix        - Computes psi function  using Hampel proposal times x                                                                                                                - :UTISTAT- 2016 May 25
%   HArho         - Computes rho function  using Hampel proposal                                                                                                                        - :UTISTAT- 2016 May 17
%   HAwei         - Computes weight function psi(u)/u using Hampel proposal                                                                                                             - :UTISTAT- 2016 May 17
%   HUeff         - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                                                   - :UTISTAT- 2016 May 11
%   HUpsi         - Computes psi function (derivative of rho function) for Huber                                                                                                        - :UTISTAT- 2016 May 17
%   HUpsider      - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                                                        - :UTISTAT- 2016 May 09
%   HUpsix        - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                                                     - :UTISTAT- 2016 May 17
%   HUrho         - Computes (rho) function for Huber                                                                                                                                   - :UTISTAT- 2016 May 17
%   HUwei         - Computes weight function psi(u)/u for Huber                                                                                                                         - :UTISTAT- 2016 May 25
%   HYPbdp        - Finds constant c which is associated to the requested breakdown                                                                                                     - :UTISTAT- 2016 May 25
%   HYPc          - Computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)                                     - :UTISTAT- 2016 May 25
%   HYPck         - Computes values of the scalars A, B, d for hyperbolic tangent estimator                                                                                             - :UTISTAT- 2016 May 25
%   HYPeff        - Finds constant c which is associated to the requested efficiency for hyperbolic estimator                                                                           - :UTISTAT- 2016 May 25
%   HYPk          - Computes breakdown point and efficiency associated with constant k=sup CVC for hyperbolic tangent estimator (for a given value of c)                                - :UTISTAT- 2016 May 25
%   HYPpsi        - Computes psi function for hyperbolic tangent estimator                                                                                                              - :UTISTAT- 2016 May 25
%   HYPpsider     - Computes derivative of psi function for hyperbolic tangent estimator                                                                                                - :UTISTAT- 2016 May 25
%   HYPpsix       - Computes psi function for hyperbolic tangent estimator times x                                                                                                      - :UTISTAT- 2016 May 25
%   HYPrho        - Computes rho function  using hyperboloc tangent estimator                                                                                                           - :UTISTAT- 2016 May 25
%   HYPwei        - Computes weight function psi(u)/u for hyperbolic tangent estimator                                                                                                  - :UTISTAT- 2016 May 09
%   inversegamcdf - Inversegampdf Inverse-gamma cumulative distribution function                                                                                                        - :UTISTAT- 2016 May 09
%   inversegaminv - Inversegampdf Inverse-gamma cumulative distribution function                                                                                                        - :UTISTAT- 2016 May 09
%   inversegampdf - Inverse-gamma probability density function                                                                                                                          - :UTISTAT- 2016 May 25
%   logmvnpdfFS   - Produces log of Multivariate normal probability density function (pdf)                                                                                              - :UTISTAT- 2016 May 17
%   mahalFS       - Computes Mahalanobis distances (in squared units) for each row of matrix Y                                                                                          - :UTISTAT- 2016 May 25
%   minscale      - Finds the M estimator of the scale for TB                                                                                                                           -         - 2016 May 09
%   Mscale        - Finds the M estimator of the scale                                                                                                                                  - :UTISTAT- 2016 May 17
%   Mscale1       - Finds the M estimator of the scale                                                                                                                                  -         - 2016 May 17
%   ncx2mixtcdf   - Cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ sigma * N(0,1))                                                         - :UTISTAT- 2016 May 09
%   normBoxCox    - Computes (normalized) Box-Cox transformation                                                                                                                        - :UTISTAT- 2016 May 25
%   normYJ        - Computes (normalized) Yeo-Johnson transformation                                                                                                                    - :UTISTAT- 2016 May 25
%   OPTbdp        - Finds the constant c associated to the supplied breakdown point                                                                                                     - :UTISTAT- 2016 May 25
%   OPTc          - Computes breakdown point and efficiency associated with constant c for Optimal rho function                                                                         - :UTISTAT- 2016 May 25
%   OPTeff        - Finds the constant c which is associated to the requested efficiency                                                                                                - :UTISTAT- 2016 May 25
%   OPTpsi        - Computes psi function (derivative of rho function) for optimal weight function                                                                                      - :UTISTAT- 2016 May 25
%   OPTpsider     - Computes derivative of psi function (second derivative of rho function) for optimal weight function                                                                 - :UTISTAT- 2016 May 25
%   OPTpsix       - Computes psi function (derivative of rho function) times x                                                                                                          - :UTISTAT- 2016 May 17
%   OPTrho        - Computes rho function for optimal weight function                                                                                                                   - :UTISTAT- 2016 May 25
%   OPTwei        - Computes weight function psi(u)/u for optimal weight function                                                                                                       - :UTISTAT- 2016 May 25
%   Powertra      - Computes power transformation (Box-Cox or  Yeo-Johnson)                                                                                                             - :UTISTAT- 2016 May 25
%   Qn            - Robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)                                                                                      - :UTISTAT- 2016 May 09
%   RandIndexFS   - Calculates Rand type Indices to compare two partitions                                                                                                              - :UTISTAT- 2016 May 25
%   Sn            - Robust estimator of scale (robust version of Gini's average difference)                                                                                             - :UTISTAT- 2016 May 25
%   tabulateFS    - Create frequency table of unique values of x, excluding possible 0 counts                                                                                           - :UTISTAT- 2016 May 25
%   TBbdp         - Finds the constant c associated to the supplied breakdown point for Tukey's biweight                                                                                - :UTISTAT- 2016 May 25
%   TBc           - Computes breakdown point and efficiency associated with constant c for Tukey's biweight                                                                             - :UTISTAT- 2016 May 25
%   TBeff         - Finds the constant c which is associated to the requested efficiency for Tukey biweight estimator                                                                   - :UTISTAT- 2016 May 25
%   TBpsi         - Computes psi function (derivative of rho function) for Tukey's biweight                                                                                             - :UTISTAT- 2016 May 25
%   TBpsider      - Computes derivative of psi function (second derivative of rho function) for Tukey's biweight                                                                        - :UTISTAT- 2016 May 25
%   TBpsix        - Computes psi function (derivative of rho function) times x for Tukey's biweight                                                                                     - :UTISTAT- 2016 May 11
