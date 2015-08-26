% UTILITIES_STAT
%
% Files :
% basicPower                   - basicPower computes the basic power transformation
% eigs_sigma                   - eigs_sigma sets a non zero value for the optional parameter sigma of function eigs
% ellipse                      - ellipse generates an ellipse given mu (location vector) and Sigma (scatter matrix)
% FSMbonfbound                 - FSMbonfbound computes Bonferroni bounds for each step of the  search (in mult analysis)
% FSRbonfbound                 - FSRbonfbound computes Bonferroni bounds for each step of the search (in linear regression)
% HAbdp                        - HAbdp finds the constant c associated to the supplied breakdown point
% HAc                          - HAc computes breakdown point and efficiency associated with constant c 
% HAeff                        - HAeff finds the tuning constant guarrantees a requested asymptotic efficiency
% HApsi                        - HApsi computes psi function  using Hampel proposal
% HApsider                     - HApsider computes derivative of psi function  using Hampel proposal
% HApsix                       - HApsix computes psi function  using Hampel proposal times x
% HArho                        - HArho computes rho function  using Hampel proposal
% HAwei                        - HAwei computes weight function psi(u)/u using Hampel proposal
% HYPbdp                       - HYPbdp finds constant c which is associated to the requested breakdown
% HYPc                         - HYPc computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)
% HYPck                        - HYPck computes values of the scalars A, B, d for hyperbolic tangent estimator
% HYPeff                       - HYPeff finds constant c which is associated to the requested efficiency for hyperbolic estimator
% HYPk                         - HYPk computes breakdown point and efficiency associated with constant k=sup CVC for hyperbolic tangent estimator (for a given value of c)
% HYPpsi                       - HYPpsi computes psi function for hyperbolic tangent estimator
% HYPpsider                    - HYPpsider computes derivative of psi function for hyperbolic tangent estimator
% HYPpsix                      - HYPpsix computes psi function for hyperbolic tangent estimator times x
% HYPrho                       - HYPrho computes rho function  using hyperboloc tangent estimator
% HYPwei                       - HYPwei computes weight function psi(u)/u for hyperbolic tangent estimator
% inversegamcdf                - inversegampdf Inverse-gamma cumulative distribution function.
% inversegaminv                - inversegampdf Inverse-gamma cumulative distribution function.
% inversegampdf                - inversegampdf Inverse-gamma probability density function.
% logmvnpdfFS                  - logmvnpdfFS log of Multivariate normal probability density function (pdf)
% mahalFS                      - mahalFS computes Mahalanobis distances (in squared units) for each row of matrix Y 
% minscale                     - minscale finds the M estimator of the scale for TB
% Mscale                       - Mscale finds the M estimator of the scale
% Mscale1                      - Mscale1 finds the M estimator of the scale
% ncx2mixtcdf                  - ncx2mixtcdf cumulative distribution function (cdf) of a linear combination of non-central chi-square (+ \sigma * N(0,1))
% normBoxCox                   - normBoxCox computes (normalized) Box-Cox transformation
% normYJ                       - normYJ computes (normalized) Yeo-Johnson transformation
% OPTbdp                       - OPTbdp finds the constant c associated to the supplied breakdown point
% OPTc                         - OPTc computes breakdown point and efficiency associated with constant c for Optimal rho function
% OPTeff                       - OPTeff finds the constant c which is associated to the requested efficiency
% OPTpsi                       - OPTpsi computes psi function (derivative of rho function) for optimal weight function
% OPTpsider                    - OPTpsider computes derivative of psi function (second derivative of rho function) for optimal weight function
% OPTpsix                      - OPTpsix computes psi function (derivative of rho function) times x
% OPTrho                       - OPTrho computes rho function for optimal weight function
% OPTwei                       - OPTwei computes weight function psi(u)/u for optimal weight function
% Powertra                     - Powertra computes power transformation (Box-Cox or  Yeo-Johnson)
% rescale                      - rescale rescales numeric array to have specified minimum (a)  and maximum (b)
% TBbdp                        - TBbdp finds the constant c associated to the supplied breakdown point for Tukey's biweight
% TBc                          - TBc computes breakdown point and efficiency associated with constant c for Tukey's biweight
% TBeff                        - Tbeff finds the constant c which is associated to the requested efficiency for Tukey biweight estimator
% TBpsi                        - TBpsi computes psi function (derivative of rho function) for Tukey's biweight  
% TBpsider                     - TBpsider computes derivative of psi function (second derivative of rho function) for Tukey's biweight  
% TBpsix                       - TBpsix computes psi function (derivative of rho function) times x for Tukey's biweight  
% TBrho                        - TBrho computes (rho) function for Tukey biweight
% TBwei                        - TBwei computes weight function psi(u)/u for Tukey biweight  
% triu2vec                     - triu2vec extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix
% zscoreFS                     - zscoresFS computes (robust) standardized z scores
