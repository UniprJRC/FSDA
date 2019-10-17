function sc = minscale(u, c, kc, initialsc, tol, maxiter)
%minscale finds the M estimator of the scale for TB
%
%  REMARK: THIS FUNCTION HAS BEEN REPLACED BY Mscale because minscale just
%  refers to Tukey's biweight At present this routine is just called by
%  Smult and MMmult or MMmultcore
%
% u = residuals or Mahalanobis distances 
% (note that u is kept fixed in each iteration)
% Remark: the scale must satisfy the following equation
% $(1/n) \sum_{i=1}^n \rho((u_i/c)/s) = kc$
% This routine computes the minimum value of s which satisfies the above
% equation%
%
%
%
% Required input arguments:
%
%    u:       : n x 1 vector which contains the residuals 
%    c        : scalar, tuning constant of the equation for Tukey biweight
%   kc        : kc= E(\rho) scalar, tuning constant linked to Tukey's biweight
%               E(\rho) = bdp* \rho(c)
%               For Tukey's biweight \rho(c) = (c^2/6) therefore 
%               E(\rho) = bdp* c^2/6
%
%  Optional input arguments:
%
%    initialsc: scalar. The initial estimate of the scale.
%               If not defined, scaled MAD of vector |u| is used.
%     maxiter : scalar. The maximum number of iterations to find the scale.
%               If not defined, maxiter is fixed to 200.
%     tol     : scalar. The tolerance for controlling convergence.
%               If not defined, tol is fixed to 1e-7.
%
%  Output:
%
%  sc : scalar, minimum value of the scale
%
% Remark: this routine is called by Smult.m
%
% References
%
% Huber, P.J. and Ronchetti, E.M. (2009), "Robust Statistics, 2nd Edition",
% Wiley. [equation 7.119,  p. 176].
%
% See also Mscale1, minscale
%
%<a href="matlab: docsearchFS('minscale')">Link to the help function</a>
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%$LastChangedDate::                      $: Date of the last commit
%


%% Beginning of code
if nargin<4
    initialsc = median(abs(u))/.6745;
end
if nargin<5
    tol = 1e-7;
end

if nargin<6
    maxiter = 200;
end

sc = initialsc;
loop = 0;
err = 1;
while  (( loop < maxiter ) && (err > tol))
    % scale step: see equation 7.119 of Huber and Ronchetti, p. 176
    % scalenew = scaleold *(1/n)*\sum  \rho(u_i/scaleold) / kc
    scnew = sc*sqrt( mean(TBrho(u/sc,c)) / kc);
    
    % Note that when there is convergence 
    % sqrt( mean(TBrho(u/sc,c)) / kc) tends to 1 (from below)
    % disp([loop sc sqrt( mean(TBrho(u/sc,c)) / kc)])
    
    err = abs(scnew/sc - 1);
    sc = scnew;
    % disp(sc)
    loop = loop+1;
end
% disp(loop)
% sc=sc;
end
