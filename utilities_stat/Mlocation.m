function newloc = Mlocation(y, psifunc, scaleest, initialmu, tol, maxiter)
%Mlocation finds the M estimator of location in a univariate sample
%
%
%<a href="matlab: docsearchFS('Mlocation')">Link to the help function</a>
%
% Required input arguments:
%
%    y:       : Response variable. Vector.
%               n-by-1 vector of observations
%               Data Types - single | double
%     psifunc : rho (psi) function. Structure.
%               A structure specifying the class of rho (psi) function to use and the
%               consistency factor.
%               psifunc must contain the following fields:
%               psifunc.class = string identyfing the rho (psi) function to use.
%                    Admissible values for class are 'bisquare' (TB),
%                    'optimal', (OPT) 'hyperbolic' (HYP), 'hampel' (HA)
%                    'power divergence' (PD)
%               psifunc.c1 = consistency factor (and other parameters)
%                   associated to required breakdown point or nominal
%                   efficiency.
%                   More precisely, psifunc.c1(1) contains consistency
%                   factor associated to required breakdown point or
%                   nominal efficiency psifunc.c1(2:end) may contain other
%                   parameters associated with the rho (psi) function.
%                   For example, if psifunc.class='hampel', c1(2:4) may
%                   contain parameters (a, b and c) of Hampel rho function.
%               Example - psifunc.class='TB';psifunc.c1=1.5476;
%               Data Types - struct
%
%  Optional input arguments:
%
%    scaleest : value of sigma. Scalar. The estimate of the scale to use in the iterative loop.
%               If not defined, scaled MAD of vector y is used.
%               Example - 0.34
%               Data Types - double
%
%    initialmu : starting value of the location estimate. Scalar.
%               The initial estimate of location to use in the first
%               iteration.
%               If not defined, initialmu is set equal to median(y).
%               Example - trimmean(y,30)
%               Data Types - double
%
%     tol     : scalar. The tolerance for controlling convergence.
%               If not defined, tol is fixed to 1e-7.
%               Example - 1e-10
%               Data Types - double
%
%     maxiter : scalar. Maximum number of iterations to find the location estimate.
%               If not defined, maxiter is fixed to 200.
%               Example - 100
%               Data Types - double
%
%  Output:
%
%  newloc : M-estimate of location. Scalar.
%           Robust M estimate of location.
%
% More About:
%
% Remark: the M estimator of location must satisfy the following equation
% \[
% \hat{\mu} = \underset{\mu}{\text{argmin}}\sum_{i=1}^n \rho \left( \frac{y_i-\mu}{\sigma} \right).
% \]
% If $\rho$ is differentiable, with
% $\psi(u) = \rho'(\psi) = \text{d}\psi(u)/\text{d}u, $
% $\hat{\mu}$ is the solution of the equation (estimating equation)
% \[
% \sum_{i=1}^n \psi \left( \frac{y_i -  \hat{\mu}}{\sigma} \right) = 0.
% \]
%
% This estimating equation shows the importance of $\psi(y)$. The equation
% can be rewritten to provide an algorithm for estimation of $\mu$ as
% \[
% \sum_{i=1}^n w_i(y_i -  \hat{\mu}) = 0 \quad \text{where} \quad w_i = \psi(y_i -  \hat{\mu})/(y_i -  \hat{\mu}).
% \]
%
% This routine computes the value of $\mu$ which satisfies the above
% equation. $\sigma$ corresponds to optional input parameter scaleest.
% Note that the value of $\sigma$ is kept fixed in each iteration.
% The iterative procedure starts with some $\hat{\mu}_0$, for example the sample median.
% In this routine $\hat{\mu}_0$ corresponds to optional input parameter
% initialmu.
% Given the estimate  $\hat{\mu}_k$ at stage $k$, compute
% \[
% w_{i,k}  =  w(y_i -  \hat{\mu}_k), \quad i = 1,\ldots,n
% \]
% \[
% \hat{\mu}_{k+1} =  \sum_{i=1}^n w_{i,k}y_i/\sum_{i=1}^n w_{i,k}.
% \]
%
%
% See also: Mscale, Mlocsca, Sreg
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Mlocation')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Example of M estimate of location.
    % M estimate of location using Tukey biweight rho function with a
    % value of c associated to a breakdown point of 0.5.
    psifunc=struct;
    psifunc.class='TB';
    bdp=0.5;
    c=TBbdp(bdp,1);
    psifunc.c1=c;
    n=100;
    shift=100;
    u=randn(n,1);
    u(1:10)=u(1:10)+shift;
    % Function Mlocation is just called with two input arguments and
    % therefore the fixed estimate of the scale which is used is the scaled
    % MAD
    loc=Mlocation(u,psifunc)
    disp('Non robust estimate of location')
    disp(mean(u))
    disp('Robust estimate of location')
    disp(loc)
%}

%{
    % Example of use of option scaleest.
    % M estimate of the location using Hampel rho function with a
    % value of c associated to a breakdown point of 0.5
    psifunc=struct;
    psifunc.class='HA'
    abc=[1.5 3.5 8];
    bdp=0.5;
    c=HAbdp(bdp,1,abc);
    psifunc.c1=[c abc];
    n=10000;
    shift=50;
    truescale=6;
    u=truescale*randn(n,1);
    u(1:10)=u(1:10)+shift;
    % Mlocation is called with three arguments. 
    % The third argument is the fixed value of the scale to use in the loop.
    locHA=Mlocation(u,psifunc,truescale+1)
%}

%{
    % Example of use of option initialmu.
    % Estimate of location using Hampel rho function.
    % M estimate of the location using Hampel rho function with a
    % value of c associated to a breakdown point of 0.5
    psifunc=struct;
    psifunc.class='HA'
    abc=[1.5 3.5 8];
    bdp=0.5;
    c=HAbdp(bdp,1,abc);
    psifunc.c1=[c abc];
    n=10000;
    shift=50;
    truescale=3;
    u=truescale*randn(n,1);
    u(1:10)=u(1:10)+shift;
    % Mlocation is called with four arguments.
    % The fourth argument is the initial estimate of mu to use.
    % Use as initial estimate of mu a trimmed mean.
    initialmu=trimmean(u,60);
    locHA=Mlocation(u,psifunc,truescale,initialmu)
%}

%{
    % Example of use of options tol and maxiter.
    % M estimate of the scale using Tukey biweight rho function with a
    % value of c associated to a breakdown point of 0.5.
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
    % Mlocation is called with 6 arguments.
    % The fifth and sixth arguments are respectively tol and maxiter
    loc=Mlocation(u,psifunc,3,median(u),1e-7,20)
%}

%{
    % Compare location estimate using two different link functions.
    psifunc=struct;
    psifunc.class='HA'
    abc=[1.5 3.5 8];
    bdp=0.5;
    c=HAbdp(bdp,1,abc);
    psifunc.c1=[c, abc];
    n=1000;
    shift=100;
    u=3*randn(n,1);
    u(100:200)=u(100:200)+shift;
    sHA=Mlocation(u,psifunc)
    psifunc=struct;
    psifunc.class='TB';
    c=TBbdp(bdp,1);
    psifunc.c1=c;
    sTB=Mlocation(u,psifunc);
    sMLE=mean(u);
    cate=categorical({'Robust location (Hampel)' 'Robust location (TB)' 'Non robust location'});
    bar(cate,[sHA sTB sMLE])
%}

%% Beginning of code
c=psifunc.c1;

% M-estimator of location using the requested rho function.
% Use the appropriate weight function
if coder.target('MATLAB')
    XXrho=strcat(psifunc.class,'wei');
    hrho=str2func(XXrho);
end

medy=median(y);
if nargin<3
    sc=median(abs(y-medy))/.6745;
else
    sc=scaleest;
end

if nargin<4
    oldloc=medy;
else
    oldloc=initialmu;
end


if nargin<5
    tol = 1e-7;
end

if nargin<6
    maxiter = 200;
end


loop = 0;
err=Inf;
while  (( loop < maxiter ) && (err > tol))


    wei=feval(hrho,(y-oldloc)/sc,c);

    % Find new estimate of location
    newloc= sum(wei.*y)/sum(wei);

    % Display
    % disp(['Iteration n.' num2str(i)   ' Diff= ' num2str(abs(newloc-oldloc))])

    % If the difference between old and new estimate of location is smaller
    % than tol get out of the loop
    err=abs(newloc-oldloc)/sc;
    % disp(err)
    oldloc=newloc;

    loop = loop+1;
end

end
%FScategory:UTISTAT