function sc = Mscale1(u, psifunc, initialsc, tol, maxiter)
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

XXpsix=strcat(psifunc.class,'psix');
hpsix=str2func(XXpsix);

% M-estimator of scale using the requested rho function.

if nargin<5
    maxiter = 200;
end

if (initialsc==0)
    s=median(abs(u))/.6745;
else
    s=initialsc;
end


% rhoold = mean(OPTrho(u/s,c)) - kc ;
rhoold = mean(feval(hrho,u/s,c)) - kc ;

iter = 0;
while (abs(rhoold) > tol) && (iter < maxiter)
    
    % delta = rhoold / mean( OPTpsix(u/s,c) ) / s;
    delta = rhoold / mean(feval(hpsix,u/s,c)) / s;
    
    
    isqu = 1; ok = 0;
    while  (isqu < maxiter && ok~=1)
        
        % rhonew = mean( OPTrho(u/(s+delta),c)) - kc;
        rhonew = mean( feval(hrho,u/(s+delta),c)) - kc;
        
        if abs(rhonew) < abs(rhoold)
            s = s + delta; ok = 1;
        else
            delta = delta/2; isqu = isqu + 1 ;
        end
    end
%     if isqu==30
%         maxiter = iter; % we tell it to stop, but we keep the iter for info
%     end
    rhoold = rhonew;
    iter = iter + 1;
end
sc = abs(s);
 disp(iter)
end