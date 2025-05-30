function [x,V,t] = Heston1D(T,n,parameters,rho,x0,V0) 
% Heston1D simulates observations and instantaneous variances from the Heston model 
%
%
%<a href="matlab: docsearchFS('Heston1D')">Link to the help function</a>
%
% Heston1D simulates, using the Euler–Maruyama method, observations and instantaneous variances from the
% model by [S. Heston, The Review of Financial Studies, Vol. 6, No. 2, 1993].
%
% Required input arguments:
%
%   T   :    Time horizon.  Scalar. 
%
%   n   :    Number of simulated observation. Scalar. 
%
%   parameters : Model parameters. Vector. Vector of dimension 4.
%
%   rho :  Leverage parameter. Scalar.  
%
%   x0  :  Initial observation value. Scalar.
%
%   V0  :  Initial variance value. Scalar.
%
%
% Optional input arguments:
%
%
% Output:
%
%   x   :    Process observations. Vector. Vector of dimension n+1.
%
%   V   :    Spot variance values. Vector. Vector of dimension n+1.
%
%   t   :    Observation times. Vector. Vector of dimension n+1.
%
% More About:
%
% The Heston model [S. Heston, The Review of Financial Studies, Vol. 6, No. 2, 1993] is given by the following stochastic
% differential equations:
% $$\left\{\begin{array}{l}
% dx_t=\left(  \mu-\frac{1}{2} \sigma^2_t \right) \,  dt +  \sigma_t \, dW_t \\
% d\sigma^2_t=\theta \, \left(  \alpha-\sigma^2_t \right) \,  dt +  \gamma \, \sigma_t \, dZ_t   
% \end{array}\right. ,$$
% where $\mu$ is a real-valued constant, $\theta, \, \alpha$  and $\gamma$ are
% positive constants, and $W$ and $Z$ are
% Brownian motions with correlation $\rho$.  
%
% See also: Heston2D.m
%

%
% References:
%
% Heston, S. (1993), A closed-form solution for options with stochastic volatility 
% with applications to bond and currency options, The Review of Financial Studies, Vol. 6, No. 2.
%
%
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('Heston1D')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of call of Heston1D for obtaining process observations only.
    % Generates observations from the Heston model.
    T=1;  
    n=23400;   
    parameters=[0,0.4,2,1];
    rho=-0.5;
    x0=log(100); 
    V0=0.4;
    x = Heston1D(T,n,parameters,rho,x0,V0);  
    figure
    plot(x)
    ylabel('Observations')
    title('Heston model')
    
%}

%{
    %% Example of call of Heston1D for obtaning process observations and volatility values.
    % Generates observations and volatilities from the Heston model.
    T=1; 
    n=23400; 
    parameters=[0,0.4,2,1];
    rho=-0.5;
    x0=log(100); 
    V0=0.4;
    [x,V] = Heston1D(T,n,parameters,rho,x0,V0);  
    figure
    subplot(2,1,1)
    plot(x)
    ylabel('Observations')
    title('Heston model')
    subplot(2,1,2)
    plot(V)
    ylabel('Spot Variances')
    title('Heston model')
%}

%{
    %% Example of call of Heston1D for obtaning process observations, volatility values and sampling times.
    % Generates observations, volatilities and sampling times from the Heston model.
    T=1; 
    n=23400; 
    parameters=[0,0.4,2,1];
    rho=-0.5;
    x0=log(100); 
    V0=0.4;
    [x,V,t] = Heston1D(T,n,parameters,rho,x0,V0);  
    figure
    subplot(2,1,1)
    plot(t,x)
    xlabel('Time')
    ylabel('Observations')
    title('Heston model')
    subplot(2,1,2)
    plot(t,V)
    xlabel('Time')
    ylabel('Spot Variances')
    title('Heston model')
%}

%% Beginning of code

dt=T/n; 

t=0:dt:T;

m=parameters(1); a=parameters(2); k=parameters(3); g=parameters(4);  

if 2*a*k-g^2 < 0 
    error('FSDA:Heston1D:WrongInputOpt','The process does not satisfy the Feller condition.');
end

mu=zeros(2,1);   
S=[1 rho; rho  1]; 

dW=sqrt(dt)*mvnrnd(mu,S,n+1); 
 
x=zeros(1,n+1);
V=zeros(1,n+1);
 
if V0 <= 0
    error('FSDA:Heston1D:WrongInputOpt','The initial variance must be positive.');
end

x(1)=x0;  
V(1)=V0;

for nn = 2 :  n+1 
    V(nn) = V(nn-1)   +   k*(a-V(nn-1))*dt     +   g*sqrt(V(nn-1))*dW(nn-1,1);
    x(nn) = x(nn-1)    +   (m-0.5*V(nn-1))*dt   +   sqrt(V(nn-1))*dW(nn-1,2);
end

end

%FScategory:FMvol