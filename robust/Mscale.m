function sc = Mscale(u, psifunc, initialsc, tol, maxiter)
%Mscale finds the M estimator of the scale
%
% u = residuals or Mahalanobis distances
% (note that u is kept fixed in each iteration)
% Remark: the M estimator of scale must satisfy the following equation
% (1/n) \sum_{i=1}^n \rho((u_i/c)/s) = kc
% This routine computes the value of s which satisfies the above
% equation
%
% Required input arguments:
%
%    u:       : n x 1 vector which contains the residuals
%     psifunc : a structure specifying the class of rho function to use, the
%               consistency factor, and the value associated with the Exp
%               of rho in correspondence of the consistency factor
%               psifunc must contain the following fields
%               c1(1) = consistency factor associated to required
%                     breakdown point of nominal efficiency
%               c1(2:end) = other paramters associated with the rho (psi)
%                     function. For example if psifunc.class='hampel'
%                     c1(2:4) must contain parameters (a, b and c) of
%                     Hampel rho function
%               kc1= Expectation for rho associated with c1
%               class = string identyfing the rho (psi) function to use.
%                    Admissible values for class are 'bisquare', 'optimal',
%                    'hyperboloc' and 'hampel
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
% Remark: this routine is called by Taureg, Sreg.m and Smult.m
%
% References
% Huber P. and Ronchetti E. (2009), Robust Statistics, Wiley (equation 7.119,  p.
% 176).
%
% Copyright 2008-2014.
% Written by FSDA team
%


%% Beginning of code
c=psifunc.c1;
kc=psifunc.kc1;

XXrho=strcat(psifunc.class,'rho');
hrho=str2func(XXrho);

% M-estimator of scale using the requested rho function.

if nargin<5
    maxiter = 200;
end

if nargin<4
    tol = 1e-7;
end

if nargin<3
    sc=median(abs(u))/.6745;
else
    sc=initialsc;
end

loop = 0;
err = 1;
while  (( loop < maxiter ) && (err > tol))
    % scale step: see equation 7.119 of Huber and Ronchetti, p. 176
    % scalenew = scaleold *(1/n)*\sum  \rho(u_i/scaleold) / kc
    % scnew = sc*sqrt( mean(TBrho(u/sc,c)) / kc);
    scnew = sc*sqrt( mean( feval(hrho,u/sc,c)) /kc);

    
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