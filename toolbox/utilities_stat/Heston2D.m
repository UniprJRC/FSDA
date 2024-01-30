function [x,V,t] = Heston2D(T,n,parameters,Rho,x0,V0) 
% Heston2D simulates observations and instantaneous variances from the bivariate Heston model 
%
%
%<a href="matlab: docsearchFS('Heston2D')">Link to the help function</a>
%
% Heston2D simulates, using the Euler–Maruyama method, observations and instantaneous variances from 
% a 2-dimensional version of the model by [S. Heston, The Review of Financial Studies, Vol. 6, No. 2, 1993].
%
% Required input arguments:
%
%   T   :    Time horizon.  Scalar. 
%
%   n   :    Number of simulated observation. Scalar.  
%
%   parameters : Model parameters. Matrix. Matrix with 4 rows and 2 columns.
%
%   Rho : Correlation vector. Vector. Vector of dimension 6.
%
%   x0  :  Initial process observation values. Vector. Vector of dimension 2.
%
%   V0  :  Initial variance values. Vector. Vector of dimension 2.
%
%
% Optional input arguments:
%
%
% Output:
%
%   x   :    Process observations. Matrix. Matrix with n+1 rows and 2 columns. 
%
%   V   :    Spot variance values. Matrix. Matrix with n+1 rows and 2 columns. 
%
%   t   :    Observation times. Vector. Vector of dimension n+1.
%
% More About:
%
% A bivariate version of the Heston model [S. Heston, The Review of Financial Studies, Vol. 6, No. 2, 1993] is given by the following stochastic
% differential equations:
% $$\left\{\begin{array}{l}
% dx^i(t)=\left(  \mu_{i}-\frac{1}{2} v^i(t) \right) \,  dt +  \sqrt{v^i(t)} \, dW^i_t \\
% dv^i(t)=\theta_{i} \, \left(  \alpha_{i}- v^i(t) \right) \,  dt +  \gamma_{i} \,  \sqrt{v^i(t)}\, dZ^i_t   
% \end{array}\right. ,\quad i=1,2,$$
% where, for each i, $v^i:=(\sigma^i)^2$, $\mu_i$ is a real-valued constant, $\theta_i, \, \alpha_i$  and $\gamma_i$ are
% positive constants, and $W_i$ and $Z_i$ are Brownian motions that are
% correlated with each other and with $W_j$ and $Z_j$, $j \neq i$.
%
% See also: Heston1D.m
%
% References:
%
% Heston, S. (1993), A closed-form solution for options with stochastic volatility 
% with applications to bond and currency options, The Review of Financial Studies, Vol. 6, No. 2.
%
%
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('Heston2D')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of call of Heston2D for obtaining process observations only.
    % Generates observations from the bivariate Heston model.
    T=1;  
    n=23400;   
    parameters=[0,0;0.4,0.4;2,2;1,1];
    Rho=[0.5;-0.5;0;0;-0.5;0.5];
    x0=[log(100),log(100)]; 
    V0=[0.4,0.4];
    x = Heston2D(T,n,parameters,Rho,x0,V0);  
    figure
    subplot(2,1,1)
    plot(x(:,1))
    ylabel('Observed values of process 1')
    title('Heston model')
    subplot(2,1,2)
    plot(x(:,2))
    ylabel('Observed values of process 2')
    title('Heston model')
%}

%{
    %% Example of call of Heston2D for obtaining process observations and volatility values.
    % Generates observations and volatilities from the bivariate Heston model.
    T=1; 
    n=23400; 
    parameters=[0,0;0.4,0.4;2,2;1,1];
    Rho=[0.5,-0.5,0,0,-0.5,0.5];
    x0=[log(100),log(100)]; V0=[0.4,0.4];
    [x,V] = Heston2D(T,n,parameters,Rho,x0,V0); 
    figure
    subplot(2,1,1)
    plot(x(:,1))
    ylabel('Observed values of process 1')
    title('Heston model')
    subplot(2,1,2)
    plot(x(:,2))
    ylabel('Observed values of process 2')
    title('Heston model')
    figure
    subplot(2,1,1)
    plot(V(:,1))
    ylabel('Spot variance of process 1')
    title('Heston model')
    subplot(2,1,2)
    plot(V(:,2))
    ylabel('Spot variance of process 2')
    title('Heston model')
%}

%{
    %% Example of call of Heston2D for obtaining process observations, volatility values and sampling times.
    % Generates observations, volatilities and sampling times from the bivariate Heston model.
    T=1; 
    n=23400; 
    parameters=[0,0;0.4,0.4;2,2;1,1];
    Rho=[0.5,-0.5,0,0,-0.5,0.5];
    x0=[log(100),log(100)]; V0=[0.4,0.4];
    [x,V,t] = Heston2D(T,n,parameters,Rho,x0,V0); 
    figure
    subplot(2,1,1)
    plot(t,x(:,1))
    xlabel('Time')
    ylabel('Observed values of process 1')
    title('Heston model')
    subplot(2,1,2)
    plot(t,x(:,2))
    xlabel('Time')
    ylabel('Observed values of process 2')
    title('Heston model')
    figure
    subplot(2,1,1)
    plot(t,V(:,1))
    xlabel('Time')
    ylabel('Spot variance of process 1')
    title('Heston model')
    subplot(2,1,2)
    plot(t,V(:,2))
    xlabel('Time')
    ylabel('Spot variance of process 2')
    title('Heston model')
%}

%% Beginning of code

dt=T/n;  

t=0:dt:T;

m=parameters(1,1); a=parameters(2,1); k=parameters(3,1); g=parameters(4,1); % parameters of process 1
m2=parameters(1,2); a2=parameters(2,2); k2=parameters(3,2) ; g2=parameters(4,2); % parameters of process 2  

if 2*a*k-g^2 < 0 
    error('FSDA:Heston2D:WrongInputOpt','Process 1 does not satisfy the Feller condition.');
end

if 2*a2*k2-g2^2 < 0 
    error('FSDA:Heston2D:WrongInputOpt','Process 2 does not satisfy the Feller condition.');
end

mu=zeros(4,1);  
S=[1     Rho(1)  Rho(2)  Rho(3); ...
   Rho(1)  1     Rho(4)  Rho(5); ...
   Rho(2)  Rho(4)  1     Rho(6); ...
   Rho(3)  Rho(5)  Rho(6)  1];   

dW=sqrt(dt)*mvnrnd(mu,S,n+1); 
 
x=zeros(n+1,2);
V=zeros(n+1,2);

if V0(1) <= 0
    error('FSDA:Heston2D:WrongInputOpt','The initial variance of process 1 must be positive.');
end

if V0(2) <= 0
    error('FSDA:Heston2D:WrongInputOpt','The initial variance of process 2 must be positive.');
end

x(1,1)=x0(1); V(1,1)=V0(1);  
x(1,2)=x0(2); V(1,2)=V0(2);  
 
for nn = 2 :  n+1 
    V(nn,1) = V(nn-1,1)   +   k*(a-V(nn-1,1))*dt      +   g*sqrt(V(nn-1,1))*dW(nn-1,1);
    
    V(nn,2) = V(nn-1,2)   +   k2*(a2-V(nn-1,2))*dt    +   g2*sqrt(V(nn-1,2))*dW(nn-1,2);
    
    x(nn,1) = x(nn-1,1)   +   (m-0.5*V(nn-1,1))*dt    +   sqrt(V(nn-1,1))*dW(nn-1,3);
    
    x(nn,2) = x(nn-1,2)   +   (m2-0.5*V(nn-1,2))*dt   +   sqrt(V(nn-1,2))*dW(nn-1,4);   
end

end
 
%FScategory:FMvol