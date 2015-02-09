function sc = Mscale1(u, psifunc, initialsc, tol, maxiter)
%Mscale1 finds the M estimator of the scale
%
% REMARK: Mscale1 has the same structure of Mscale but a different (less efficient)
% algorithm is used to find the solution
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
%               kc1= Expectation of rho associated with c1 to get a consistent
%                    estimator at the model distribution kc1=E(rho)
%               class = string identyfing the rho (psi) function to use.
%                    Admissible values for class are 'bisquare' (TB),
%                    'optimal', (OPT) 'hyperbolic' (HYP) and 'hampel' (HA)
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
%  sc : scalar, M-estimate of the scale
%
% Remark: this routine in the regression context is never called
% because the more efficient function Mscale is used
%
% References
% Huber P. and Ronchetti E. (2009), Robust Statistics, Wiley 
%
% See also Mscale, minscale
%
% Copyright 2008-2015.
% Written by FSDA team
%

% Last modified 06-Feb-2015

% Examples

%{
    % M estimate of the scale using Tukey biweight rho function with a
    % value of c associated to a breakdown point of 0.5
    psifunc=struct;
    psifunc.class='TB';
    bdp=0.5;
    c=TBbdp(bdp,1);
    % kc = E(rho) = sup(rho)*bdp
    kc=c^2/6*bdp;
    psifunc.c1=c;
    psifunc.kc1=kc;
    n=10000;
    shift=5;
    u=2*randn(n,1);
    u(1:10)=u(1:10)+shift;
    Mscale1(u,psifunc);
%}

%{
    % M estimate of the scale using Hampel biweight rho function with a
    % value of c associated to a breakdown point of 0.5
    psifunc=struct;
    psifunc.class='HA'
    abc=[1.5 3.5 8];
    bdp=0.5;
    c=HAbdp(bdp,1,abc);
    % kc = E(rho) = sup(rho)*bdp
    kc=HArho(c*abc(3),[c, abc])*bdp;
    psifunc.c1=[c abc];
    psifunc.kc1=kc;
    n=10000;
    shift=5;
    u=3*randn(n,1);
    u(1:10)=u(1:10)+shift;
    s=Mscale1(u,psifunc);
%}

%% Beginning of code
c=psifunc.c1;
kc=psifunc.kc1;

XXrho=strcat(psifunc.class,'rho');
hrho=str2func(XXrho);

XXpsix=strcat(psifunc.class,'psix');
hpsix=str2func(XXpsix);

% M-estimator of scale using the requested rho function.

if nargin<5
    maxiter = 1200;
end

if nargin<4
    tol = 1e-7;
end


if nargin<3
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
    % Note that (mean(feval(hrho,u/s,c)) - kc) / mean(feval(hpsix,u/s,c)) 
    % = rhoold / mean(feval(hpsix,u/s,c)) = IF
    
    isqu = 1; ok = 0;
    while  (isqu < maxiter && ok~=1)
        
        % rhonew = mean( OPTrho(u/(s+delta),c)) - kc;
        rhonew = mean( feval(hrho,u/(s+delta),c)) - kc;
        
        if abs(rhonew) < abs(rhoold)
            s = s + delta; ok = 1;
        else
            delta = delta/2; isqu = isqu + 1 ;
        end
        % delta should become smalle and smaller
        % disp(delta)
    end
%     if isqu==30
%         maxiter = iter; % we tell it to stop, but we keep the iter for info
%     end
    rhoold = rhonew;
    disp(rhonew)
    iter = iter + 1;
end
sc = abs(s);
 disp(['Number of iterations=' num2str(iter)])
end