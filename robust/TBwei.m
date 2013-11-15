function w = TBwei(u,c)
%TBwei computes weight function psi(u)/u for Tukey biweight  
%
%<a href="matlab: docsearch('tbwei')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    u:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator 
%               (beta in regression or mu in the location case ...) 
%
% Function TBwei transforms vector u as follows 
% TBwei(u)=
% (c^2/6) psi(u)/u = (c^2/6) [ 1-(u/c)]^2 if |u/k|<=1 0 otherwise
% See p. 30 of Maronna et al. (2006)
%
%
% Remark: Tukey's biweight  psi-function is almost linear around u = 0 in accordance with
% Winsor's principle that all distributions are normal in the middle.
% This means that  \psi (u)/u is approximately constant over the linear region of \psi,
% so the points in that region tend to get equal weight.
%
%
% References:
%
% ``Robust Statistics, Theory and Methods'' by Maronna, Martin and Yohai;
% Wiley 2006.
%
%
% Copyright 2008-2013.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti
%
%
%<a href="matlab: docsearch('tbwei')">Link to the help page for this function</a>
% Last modified 02-May-2013
%
% Examples:

%{

x=-6:0.01:6;
weiTB=TBwei(x,2);
plot(x,weiTB)
xlabel('x','Interpreter','Latex')
ylabel('$W (x) =\psi(x)/x$','Interpreter','Latex')

%}

%% beginning of code

w = (1 - (u/c).^2).^2;

% The following instruction is unnecessary
% however it is the proper expression for the weights
% if we start with the normalized \rho (\infty)=1
% w = w .* (c^2/6);

w( abs(u/c) > 1 )= 0;
end
