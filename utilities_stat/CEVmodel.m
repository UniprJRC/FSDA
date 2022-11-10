function [S,A]=CEVmodel(t,x)
% CEVmodel computes price and instantaneous variance processes from the CEV model 
%
%
%<a href="matlab: docsearchFS('CEVmodel')">Link to the help function</a>
%
% CEVmodel computes price and instantaneous variance for the Constant
% Elasticity of Variance model [S. Beckers, The Journal of Finance, Vol.
% 35, No. 3, 1980] via Euler method
%
% Required input arguments:
%
%   t   :    Discrete time grid. Vector.  Row or column vector. 
%
%   x   :    Initial price value. Scalar. 
%
%
% Optional input arguments:
%
%
% Output:
%
%   S   :    Spot prices. Vector. Column vector with the same length of t. 
%
%   A   :    Spot variance values. Vector. Column vector with the same length of t.
%
% More About:
%
% The Constant Elasticity of Variance model [S. Beckers, The Journal of
% Finance, Vol. 35, No. 3, 1980] is given by the following stochastic
% differential equation
% $$\left\{\begin{array}{l}
%dS_t= \sigma \, S_t^{\delta} \, dW_t \\
%S_0=x, 
%\end{array}\right.$$
% where $\sigma$  and $\delta$ are positive constants and $W$ is a
% Brownian motion on a filtered probability space. We assume $\sigma=0.3$
% and $\delta=1.5$. The instanteneous variance is given by
% $$A_t=\sigma^2S_t^{2(\delta-1)}.$$
%
% See also: FE_spot_vol.m, FE_spot_vol_FFT.m
%
% References:
%
% S. Beckers, "The Constant Elasticity of Variance Model and Its
% Implications For Option Pricing", The Journal of Finance, Vol. 35, No. 3,
% 1980.
%
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('CEVmodel')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of call of CEVmodel providing only price values.
    % Generates spot prices for the CEV model at times t.
    n=1000; dt=1/n; 
    t=0:dt:1; % discrete time grid
    x=100; % initial price value
    S=CEVmodel(t,x); % spot prices
    plot(t,S)
    xlabel('Time')
    ylabel('Spot price')
    title('CEV model')
%}

%{
    %% Example of call of CEVmodel providing both price and variance values.
    % Generates price and instantaneous variance values for the CEV model at times t.
    n=1000; dt=1/n; 
    t=0:dt:1; % discrete time grid
    x=100; % initial price value
    [S,A]=CEVmodel(t,x); % spot prices and variance
    subplot(2,1,1)
    plot(t,S)
    xlabel('Time')
    ylabel('Spot price')
    title('CEV model')
    subplot(2,1,2)
    plot(t,A)
    xlabel('Time')
    ylabel('Spot variance')
    title('CEV model')
%}

%% Beginning of code

n=length(t); S=zeros(n,1); S(1)=x; sigma=0.3; delta=1.5;

for i=1:n-1
    dt=t(i+1)-t(i); W=randn;
    S(i+1)=S(i)+sigma*S(i)^delta*sqrt(dt)*W; % spot price
end
A=sigma^2*S.^(2*(delta-1)); % spot variance

end
%FScategory:UTISTAT