function outIRWLS = Mlocsca(y,psifunc,initialmu,initialscale,tol,maxiter)
%Mlocsca computes M estimator of location and scale in univariate samples
%
%<a href="matlab: docsearchFS('Mlocsca')">Link to the help function</a>
%
% This routines compute simultaneous M estimates of location and scale in
% univariate samples.
%
%  Required input arguments:
%
%    y:         A vector with n elements which contains the response variable.
%               It can be both a row or column vector.
%
%     psifunc : a structure specifying the class of rho function to use, the
%               consistency factor, and the value associated with the
%               Expectation of rho in correspondence of the consistency
%               factor. struct.
%               psifunc is a struct which contains the following fields
%               psifunc.class = string identyfing the rho (psi) function to use.
%                    Admissible values for class are 'bisquare', 'optimal'
%                    'hyperbolic' and 'hampel'
%               psifunc.c1 = consistency factor associated to required
%                    breakdown point
%                   More precisely, psifunc.c1(1) contains consistency
%                   factor associated to required breakdown point or
%                   nominal efficiency psifunc.c1(2:end) may contain other
%                   parameters associated with the rho (psi) function.
%                   For example, if psifunc.class is 'hampel', c1(2:4) must
%                   contain parameters (a, b and c) of Hampel rho function.
%                   If length(psifunc.c1)=1 psifunc.class is HA, then
%                   default parameters for Hampel rho function are used.
%               Remark: if class is 'hyperbolic' it is also possible to
%                   specify parameters k (sup CVC), A, B and d
%               psifunc.kc1= Expectation of rho associated with c1
%               Data Types - struct
%
%
%  Optional input arguments:
%
%    initialmu : starting value of the location estimate. Scalar.
%               The initial estimate of location to use in the first
%               iteration.
%               If not defined, initialmu is set equal to median(y).
%               Example - 0.34
%               Data Types - double
% initialscale : starting value of the dispersion estimate. Scalar.
%               The initial estimate of the scale to use in the first
%               iteration.
%               If not defined, initialscale is set equal to 1.4826*mad(y,1).
%               Example - 1.32
%               Data Types - double
%
%   tol     : scalar. The tolerance for controlling convergence.
%               If not defined, tol is fixed to 1e-7.
%               Example - 1e-10
%               Data Types - double
%
%     maxiter : scalar. Maximum number of iterations to find the location estimate.
%               If not defined, maxiter is fixed to 200.
%               Example - 100
%               Data Types - double
%
%
%  Output:
%
%   outIRWLS : a structure containing the following fields:
%      outIRWLS.location  = Location estimate. Estimate of location after
%           refsteps refining steps
%     outIRWLS.scale  = scale estimate. Estimate of scale after refsteps
%           refining step
%     outIRWLS.weights = n-by-1 vector. Weights assigned to each
%           observation
%
% More About:
%
%
% In the IRWLS (iterative reweighted least square) procedure the value of
% location and the value of the scale are updated in each step. In order to
% find simultaneous estimates of location and dispersion we need to solve
% the system of two equations
% \[
% \sum_{i=1}^n \psi_\text{loc}\left( \frac{y_i - \hat{\mu}}{\hat{\sigma}}\right)  =  0
% \]
% \[
% \sum_{i=1}^n \rho_\text{scale}\left( \frac{y_i - \hat{\mu}}{\hat{\sigma}}\right)  =  K.
% \]

% In the two equations above we distinguish between $\rho_\text{scale}$ used for scale estimation
% and $\psi_\text{loc}$ and its derivatives used for location. The two need
% not be different and are, indeed, often the same. In this routine we
% assume that they are the same and are specified in input parameter psifunc.
% $K$ corresponds to input parameter psifunc.kc1.
%
% Given starting values $\hat{\mu}_0$ and $\hat{\sigma}_0$ the pair of reweighting equations moves forward from stage $k$ using the calculations
% [1] Find the location weights
% $ w_{i,k} = w\{(y_i - \hat{\mu}_k)/\hat{\sigma}_k\}$,
% Note that in order to find the weights we need psifunc.c1, the tuning
% constant associated to a nominal value of breakdown point or efficiency.
% [2] Calculate the new location estimate as
% \[
% \hat{\mu}_{k+1} =
% \sum_{i=1}^n w_{i,k}y_i / \sum_{i=1}^n w_{i,k}.
% \]
% The new (squared) scale estimate is
% [3]
% \[
% \hat{\sigma}^2_{k+1} = \hat{\sigma}^2_k\{{1}/(n K)\}\sum_{i=1}^n
% \rho_\text{scale}\{(y_i - \hat{\mu}_{k+1})/\hat{\sigma}_k\}.
%\]
% Note that in order to compute $\rho_\text{scale}$ we need psifunc.c1, the
% tuning constant associated to a nominal value of breakdown point or
% efficiency.  $K$ corresponds to input parameter psifunc.kc1.
% [4] Return to Step 1 until the change in the estimates defined as
% \[
%  |\mu_{k+1}-\mu_{k}|/|\mu_{k}| +  |\sigma_{k+1}-\sigma_{k}|/\sigma_{k}
% \]
% is less than a
% prespecified tolerance (optional input argument tol) or $k$
% is equal to optional input parameter maxiter (maximum number of
% iterations).
%
% This alternating algorithm converges to a point with zero derivatives,
% which may be a minimum, a maximum or a saddle point.
%
%
% See also: Mlocation, Mscale, Sreg
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Mlocsca')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Example of use of Mlocsca with two input arguments.
    % Example of M estimate of location and scale using TB function.
    psifuncTB=struct;
    psifuncTB.class='TB';
    bdp=0.5;
    c=TBbdp(bdp,1);
    % kc = E(rho) = sup(rho)*bdp
    kc=c^2/6*bdp;
    psifuncTB.c1=c;
    psifuncTB.kc1=kc;
    % trueloc and truescale are the true values of the two parameters
    trueloc=5;
    truescale=2;
    
    % fraction of units to contaminate
    fraccontamination=0.3;
    % shift in contamination
    shift=100;
    % n is sample size
    n=2000;
    rng(1000)
    
    u=truescale*randn(n,1)+trueloc;
    n1=round(n*fraccontamination);
    % shift contamination
    u(1:n1)=u(1:n1)+shift;
    % Specify initial estimate of location.
    initialmu=10;
    out=Mlocsca(u,psifuncTB,initialmu);
    disp('Robust estimate of location')
    disp(out.location)
    disp('Robust estimate of scale')
    disp(out.scale)
%}

%{
    % Example of use of Mlocsca with three input arguments.
    % M estimate of location ans scale using Hampel rho function with a
    % value of c associated to a breakdown point of 0.5
    psifunc=struct;
    psifunc.class='HA';
    abc=[1.5 3.5 8];
    bdp=0.5;
    c=HAbdp(bdp,1,abc);
    psifunc.c1=[c abc];
    % kc = E(rho) = sup(rho)*bdp
    kc=HArho(c*abc(3),[c, abc])*bdp;
    psifunc.c1=[c, abc];
    psifunc.kc1=kc;
     % trueloc and truescale are the true values of the two parameters
    trueloc=5;
    truescale=2;
    
    % n is sample size
    n=20000;
    u=truescale*randn(n,1)+trueloc;
    % Just 2 outliers with point mass contamination
    u(1:2)=10000000;
    out=Mlocsca(u,psifunc);
    disp('Robust estimate of location')
    disp(out.location)
    disp('Robust estimate of scale')
    disp(out.scale)
%}

%{
    % Example of use of Mlocsca with four input arguments.
    % M estimate of location and scale using optimal rho function with a
    % value of c associated to a breakdown point of 0.3
    psifuncOPT=struct;
    psifuncOPT.class='OPT';
    bdp=0.3;
    cOPT=OPTbdp(bdp,1);
    rhoOPTsup=OPTrho(200000,1);
    % rhoOPTsup=1;
    psifuncOPT.c1=cOPT;
    % kc = E(rho) = sup(rho)*bdp
    psifuncOPT.kc1=rhoOPTsup*bdp;
    % trueloc and truescale are the true values of the two parameters
    trueloc=5;
    truescale=2;
    % n is sample size
    n=200;
    u=truescale*randn(n,1)+trueloc;
    % Just 10 outliers with point mass contamination
    u(1:10)=10000000;
    % 10 and 8 are our initial estimates of location and scale
    out=Mlocsca(u,psifuncOPT,10,8);
    disp('Robust estimate of location')
    disp(out.location)
    disp('Robust estimate of scale')
    disp(out.scale)
%}

%{
    %% Example of use of Mlocsca with five input arguments.
    % M estimate of location ans scale using power divergence rho function with a
    % value of c associated to a breakdown point of 0.2
    psifuncPD=struct;
    psifuncPD.class='PD';
    bdp=0.2;
    c1=PDbdp(bdp);
    psifuncPD.c1=c1;
    psifuncPD.kc1=bdp;
    % trueloc and truescale are the true values of the two parameters
    trueloc=5;
    truescale=2;
    
    % n is sample size
    n=200;
    u=truescale*randn(n,1)+trueloc;
    % Just 10 outliers with point mass contamination
    u(1:20)=10000000;
    % 10 and 8 are our initial estimates of location and scale
    % Set the tolerance
    tol=1e-20;
    outIRWLS=Mlocsca(u,psifuncPD,10,8,tol);
    disp('Robust estimate of location')
    disp(outIRWLS.location)
    disp('Robust estimate of scale')
    disp(outIRWLS.scale)
%}

%{
    % Example of use of Mlocsca with six input arguments.
    % M estimate of location and scale using hyperbolic rho function with a
    % value of c associated to a breakdown point of 0.11
    psifuncHYP=struct;
    psifuncHYP.class='HYP';
    bdp=0.11;
    [cHYP,A,B,d]=HYPbdp(bdp,1);
    k=6;
    rhoHYPsup=HYPrho(200000,[cHYP,k,A,B,d]);
    % rhoHAsup=1;
    psifuncHYP.c1=[cHYP,k,A,B,d];
    % kc = E(rho) = sup(rho)*bdp
    psifuncHYP.kc1=rhoHYPsup*bdp;
    % trueloc and truescale are the true values of the two parameters
    trueloc=5;
    truescale=2;
    
    % n is sample size
    n=200;
    u=truescale*randn(n,1)+trueloc;
    % Just 10 outliers with point mass contamination
    u(1:20)=10000000;
    % 10 and 8 are our initial estimates of location and scale
    % Set the tolerance
    tol=1e-20;
    % Set maximum number of iterations
    maxiter=500;
    out=Mlocsca(u,psifuncHYP,10,8,tol,maxiter);
    disp('Robust estimate of location')
    disp(out.location)
    disp('Robust estimate of scale')
    disp(out.scale)
%}

%% Beginning of code
c=psifunc.c1;
kc=psifunc.kc1;
mediany=median(y);
if nargin < 3
    initialmu = mediany;
end

% Residuals for the initialbeta
res = y - initialmu;

% The scaled MAD of residuals is the initial scale estimate default value
if nargin < 4
    initialscale = median(abs(y-mediany))/.6745;
end

if nargin < 5
    tol = 1e-07;
end

if nargin < 6
    maxiter = 200;
end

mu = initialmu;
scale = initialscale;
newmu =initialmu;


XXrho=strcat(psifunc.class,'rho');
hrho=str2func(XXrho);

XXwei=strcat(psifunc.class,'wei');
hwei=str2func(XXwei);

iter = 0;
betadiff = 9999;

while ( (betadiff > tol) && (iter < maxiter) )
    iter = iter + 1;


    % Compute n x 1 vector of weights (using TB)
    weights = feval(hwei,res/scale,c);
    % weights = TBwei(res/scale,c);
    % Compute new estimate of location
    newmu=sum(y.*weights)/sum(weights);

    % Solve for the scale
    meanrho=mean(feval(hrho,res/scale,c));
    oldscale=scale;
    scale = scale * sqrt(meanrho / kc );

    % betadiff is linked to the tolerance (specified in scalar reftol)
    betadiff = abs(mu - newmu) / abs(mu)+ abs(scale - oldscale) / oldscale;

    % update residuals
    res = y - newmu;
    mu = newmu;

end

% store final estimate of location
outIRWLS.location = newmu;
% store final estimate of scale
outIRWLS.scale = scale;
% store final estimate of the weights for each observation
outIRWLS.weights=weights;

end
%FScategory:UTISTAT