function Nopt = OptimalCuttingFrequency(x,t)
%OptimalCuttingFrequency computes the optimal cutting frequency for the Fourier estimator of integrated variance
%
%
%<a href="matlab: docsearchFS('OptimalCuttingFrequency')">Link to the help function</a>
%
% OptimalCuttingFrequency computes the optimal cutting frequency for
% running  the Fourier estimator of the integrated variance
% on noisy timeseries data.
%
% Required input arguments:
%
%   x   :    Observation values. Vector. Row or column vector containing
%              the observed values.
%   t    :     Observation times. Vector.  Row or column vector with the same
%              length of x containing the observation times
%
%
% Optional input arguments:
%
%
% Output:
%
%   Nopt  :  Optimal cutting frequency. Scalar. Integer representing the
%           optimal cutting frequency.
%
%
%
% More About:
%
% We assume our timeseries data are noisy observations $\tilde x$ from a diffusion
% process following the Ito stochastic differential equation
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$
% where $W$ is a Brownian motion on a filtered probability space. Let
% $\sigma$ and $b$ be random processes, adapted to the Brownian filtration.
% The integrated variance of the process over the time interval $[0,T]$ is defined as
% $$\int_0^T \sigma^2(t) dt.$$

% For any positive integer $n$, let ${\cal S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times.
% The observations are affected by i.i.d. noise terms $\eta(t_i)$  with
% mean zero and finite variance
% $$\tilde x(t_i)=x(t_i)+\eta(t_i).$$
% See the Reference for further  mathematical details.
% Moreover, let $\delta_i(\tilde x):= \tilde x(t_{i+1})-\tilde x(t_i)$ be
% the increments of $\tilde x$.
%
% The optimal cutting frequency $N$ for computing  the Fourier estimator of
% the integrated variance from noisy timeseries data is obtained by
% minimization of the estimated MSE.
%
% The Fourier estimator of the integrated variance over $[0,T]$, is
% then defined as
% $$\widehat\sigma^{2}_{n,N}:= {T^2 \over {2N+1}}\sum_{|s|\leq N} c_s(d\tilde x_n)
% c_{-s}(d\tilde x_n),$$
% where for any integer $k$, $|k|\leq N$, the discretized Fourier
% coefficients of the increments are
% $$c_k(d\tilde x_{n}):= {1\over {T}} \sum_{i=0}^{n-1} e^{-{\rm i} {{2\pi}\over {T}}
% kt_i}\delta_i(\tilde x).$$
%
% See also: FE_int_vol.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin
% Volatility Estimation. Theory and Practice, "Springer Briefs in
% Quantitative Finance", Springer.
%
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('OptimalCuttingFrequency')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Computation of the optimal cutting frequency for estimating the integrated variance from a
    % vector x of noisy observations of a univariate diffusion process.
    % Generate data.
    n=1000;
    dt=1/n;
    t=0:dt:1;
    x=randn(n,1)*sqrt(dt);
    % generate the diffusion process
    x=[0;cumsum(x)];
    % Add noise:
     noise_to_signal =0.5; % noise-to-signal ratio
     sigma_eps = noise_to_signal*std(diff(x));
     noise=sigma_eps*randn(size(x));
     % add noise, which is i.i.d. N(0,sigma_eps^2)
     x=x+noise;
     Nopt = OptimalCuttingFrequency(x,t); % optimal cutting frequency
     ivar=FE_int_vol(x,t,'N',Nopt);
    disp(['The optimal cutting frequency is: ' num2str(Nopt)])
    disp(['The value of the integrated variance is: ' num2str(ivar)])
%}



%% Beginning of code

% Make sure that x and t are both column vector.
x=x(:);
t=t(:);

n=length(x);
if n ~= length(t)
    error('FSDA:OptimalCuttingFrequency:WrongInputOpt','Input arguments x and t must have the same length.');
end

% Sparse sampling and preliminary computation of integrated variance and
% quarticity. Sparse sampling removes the spurious first order
% autocorrelation of the increments due to noise contamintion:
[acf,~,bounds] = autocorr(diff(x));
step=1;
while (acf(2) <  bounds(2))
    step=step+1; xsparse=x(1:step:end);
    [acf,~,bounds] = autocorr(diff(xsparse));
end
RV=sum(diff(xsparse).^2); % integrated realized variance
Q=IntegratedQuarticity(xsparse,length(xsparse)-1,(2*pi)); % integrated quarticity over [0,2 pi]

% Compute moments of noise
[alpha,beta,gamma,E2,Eeta4]=noise_moments(x,n-1,RV);
% BIAS and MSE estimates
[~,MSE] = estimates(length(x)-1,floor((n-1)/2),Q,alpha,beta,gamma,E2,Eeta4);
% optimal cutting frequency Nopt
[~,Nopt]=min(MSE);
Nopt=min(Nopt,floor((n-1)/2));

end


function [BIAS,MSE] = estimates(n,N,Q,alpha,beta,gamma,E2,Eeta4)
%estimates computes BIAS and MSE
%
% Input Arguments:
%
%   n          :    total number of increments
%   N          :   maximum Fourier frequency
%   Q          :   integrated quarticity
%   alpha     :   estimations of alpha, see Reference
%   beta       :   estimations of beta, see Reference
%   gamma  :   estimations of gamma, see Reference
%   E2         :   estimations of E2, see Reference
%   Eeta4    :    estimations of Eeta4, see Reference
%
%
% Output arguments:
%
%   BIAS    row vector of bias estimates as a function
%               of the cutting frequencies
%   MSE     row vector of MSE estimates as a function
%               of the cutting frequencies

h=2*pi/n;
MSE=zeros(1,N); BIAS=zeros(1,N);
alphaFE=zeros(1,N); betaFE=zeros(1,N);
gammaFE=zeros(1,N);

for k=1:N
    trunc=min(n/2,k);
    rD=rDirichletKernel(trunc,h);
    BIAS(k)=n*E2*(1-rD);
    
    alphaFE(k)=alpha*(1+rD^2-2*rD);
    
    betaFE(k)=beta*(1+rD^2-2*rD);
    
    gammaFE(k)=gamma+4*(Eeta4+E2^2)*(2*rD-rD^2)+4*pi*Q/(2*trunc+1);
    
    MSE(k) = 2*Q*h+betaFE(k)*n+alphaFE(k)*n^2+gammaFE(k);
end

end

function d = rDirichletKernel(N,t)
% Rescaled Dirichlet kernel:
d=1;
for s=1:N
    d=d+2*cos(s*t);
end
d=d/(2*N+1);
end

function Q=IntegratedQuarticity(P,M,h)
% Integrated quarticity
Q=sum((diff(P)).^4)*M/(3*h);
end

function [alpha,beta,gamma,E2,Eeta4] = noise_moments(P,n,V)
% Computes the sample moments of the noise

% Input variables:
%   P       column vector of the observed log-prices
%   n       total number of increments
%   V       integrated variance

r=diff(P); % log-returns

E2=sum(r.^2)/n-V/n;
E4=sum(r.^4)/n-(6*E2*V)/n;
Eeta2=E2/2;
Eeta4=E4/2-3*E2^2/4;

alpha=E2^2;
beta=4*Eeta4;
gamma=8*Eeta2*V+alpha/2-2*Eeta4;

end