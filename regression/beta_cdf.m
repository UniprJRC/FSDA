function cdf = beta_cdf(x, a, b)
% PURPOSE: cdf of the beta distribution
%--------------------------------------------------------------
% USAGE: cdf = beta_cdf(x,a,b)
% where:   x = prob[beta(a,b) <= x], x = vector
%          a = beta distribution parameter, a = scalar 
%          b = beta distribution parameter  b = scalar 
% NOTE: mean [beta(a,b)], variance = ab/((a+b)*(a+b)*(a+b+1))
%--------------------------------------------------------------
% RETURNS: cdf at each element of x of the beta distribution
%--------------------------------------------------------------
% SEE ALSO: beta_d, beta_pdf, beta_inv, beta_rnd
%--------------------------------------------------------------
  

if any(any((a<=0)|(b<=0)))
   error('Parameter a or b is nonpositive');
end

cdf = betainc(x,a,b);

