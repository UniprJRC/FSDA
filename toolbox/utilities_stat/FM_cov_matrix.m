function var_cov_matrix = FM_cov_matrix(N,T,data)
% FM_cov_matrix computes the integrated variance-covariance matrix of any number of processes from discrete observations, via the Fourier-Malliavin estimator with Fejer kernel
%
%
%<a href="matlab: docsearchFS('FM_cov_matrix')">Link to the help function</a>
%  
% FM_cov_matrix computes the integrated variance-covariance matrix of any number of processes from discrete observations over the time horizon [0,T], 
% using the Fourier-Malliavin estimator with Fejer kernel
%
% Required input arguments:
%
%   N     :    Cutting frequency. Scalar.   
%
%   T     :    Estimation horizon. Scalar. 
%
%   data     :   Observed process values and observation times. Cell array.
%   Cell array containing vectors of process values and observation times of different sizes.
%
%
%
% Optional input arguments:
%
% Output:
%
% var_cov_matrix  : Integrated variances-covariance matrix. Matrix. Returns a
%                   data matrix d x d of double.
%
% More About:
%
% We assume that vectors xi contain discrete observations from a d-dimensional diffusion
% process $(x_1,..,x_d)$ following the Ito stochastic differential equation 
% $$dx_i(t)= \sigma_i(t) \ dW_i(t) + b_i(t) \ dt, \quad i=1,...,d,$$ 
% where, for all i, $W_i$ is a Brownian motion defined on the filtered probability space 
% $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, while $\sigma_i$ and $b_i$ are random processes,
% adapted to  $\mathcal{F}_t$.
% See the Reference for further  mathematical details.
% The integrated covariance of the process $(x_i,x_j)$ on $[0,T]$ is defined as
% $$IC_{ij,[0,T]}:=\rho_{ij} \, \int_0^T \sigma_i(t) \sigma_j(t) \, dt$$
% where $\rho_{ij}$ denotes the correlation between $W_i$ and $W_j$.
% 
% For any positive integer $n_i$, let $\mathcal{S}^i_{n_i}:=\{ 0=t^i_{0}\leq \cdots
% \leq t^i_{n_i}=T  \}$ be the observation times for the $i$-th asset. Moreover, let $\delta_l(x^i):=
% x^i(t^i_{l+1})-x^i(t^i_l)$ be the increments of $x^i$. 
% The Fourier estimator of the integrated covariance with fejer kernel on $[0,T]$ is
% given by 
% $$T c_0(c_{n_i,n_j,N})={T^2\over {N+1}} \sum_{|k|\leq N} \left( 1- \frac{|k|}{N+1}\right)c_{k}(dx^i_{n_1})c_{-k}(dx^j_{n_2}),$$
%
% $$c_k(dx^i_{n_i})= {1\over {T}} \sum_{l=0}^{n_i-1} e^{-{\rm i}\frac{2\pi}{T}kt^i_l}\delta_{l}(x_i).$$
% 
% See also: FE_int_vol.m, FE_int_vol_Fejer.m, FM_int_vol, FM_int_cov.m, FM_int_quart.m, FM_int_volvol.m, FM_int_lev.m, Heston1D.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin Volatility Estimation. Theory and Practice, "Springer Briefs in Quantitative Finance", Springer.
% 
% Sanfelici, S., Toscano, G. (2024), The Fourier-Malliavin Volatility (FMVol) MATLAB toolbox, available on ArXiv.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FM_cov_matrix')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:  

%{ 
    %% Example of call of FM_cov_matrix.
    % The following example estimates the daily integrated variance of a univariate 
    % Heston model, from a discrete sample. The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1; parameters=[0,0.4,2,1]; rho=-0.5;
x0=log(100); V0=0.4;
n=10000; [x1,V1,t1] = Heston1D(T,n,parameters,rho,x0,V0); 

% Integrated variance estimation 
data={x1,t1}; % Create the cell array of values-times data; x1 is a row or column vector containing the observed process values.
%   t1 is a row or column vector with the same length of x1 containing the observation times over [0,T].
var_cov_matrix = FM_cov_matrix(floor(n/2),T,data) 
%}

%{ 
    %% Example of call of FM_cov_matrix with two processes.
    % The following example estimates the integrated variance-covariance of two processes from discrete samples of different sizes. 

T=1; 
n1=1000; dt=1/n1; t1=0:dt:T; x1=randn(n1,1)*sqrt(dt); x1=[0;cumsum(x1)];
n2=2000; dt=1/n2; t2=0:dt:T; x2=randn(n2,1)*sqrt(dt); x2=[0;cumsum(x2)];

% Integrated covariance estimation 
N=min(floor(n1/2),floor(n2/2)); % cutting frequency
data={x1,t1,x2,t2}; % Create the cell array of values-times data; xi is a row or column vector containing the observed process values.
%   ti is a row or column vector with the same length of xi containing the corresponding observation times over [0,T].
var_cov_matrix = FM_cov_matrix(N,T,data) 
%}

%{ 
    %% Example of call of FM_cov_matrix with three processes.
    % The following example estimates the integrated variance-covariance of three processes from discrete samples of different sizes. 

T=1; 
n1=10000; dt=1/n1; t1=0:dt:T; x1=randn(n1,1)*sqrt(dt); x1=[0;cumsum(x1)];
n2=5000; dt=1/n2; t2=0:dt:T; x2=randn(n2,1)*sqrt(dt); x2=[0;cumsum(x2)];
n3=9000; dt=1/n3; t3=0:dt:T; x3=randn(n3,1)*sqrt(dt); x3=[0;cumsum(x3)];

data={x1,t1,x2,t2,x3,t3}; % Create the cell array of values-times data; xi is a row or column vector containing the observed process values.
%   ti is a row or column vector with the same length of xi containing the corresponding observation times over [0,T].

N=2500; % cutting freq. 
var_cov_matrix = FM_cov_matrix(N,T,data)

%}



%% Beginning of code

 n_vectors=length(data);
    if round(2*round(n_vectors/2)) ~= n_vectors
        error('FSDA:FM_cov_matrix:WrongInputOpt','The number of input vectors in input data must be even.');
    end

 n_columns=n_vectors/2; n=length(data{:,2});
 for i=1:n_columns
     if length(data{:,1+2*(i-1)}) ~= length(data{:,2*i})
         error('FSDA:FM_cov_matrix:WrongInputOpt','The vectors of observation values and times must have the same length.');
     end
     t_i=data{:,2*i}; n=max(n,length(t_i));
     if (abs(t_i(end)-T) > 1e-6)
         error('FSDA:FM_cov_matrix:WrongInputOpt','The time horizon must be the same for all processes.');
     end

    R=zeros(n-1,n_columns); t=zeros(n-1,n_columns); 
    for i=1:n_columns
        R(1:length(data{:,1+2*(i-1)})-1,i)=diff(data{:,1+2*(i-1)});
        t_i=data{:,2*i}; t(1:length(data{:,2*i})-1,i)=t_i(1:end-1);
    end

    const=2*pi/T;
    c_p=zeros(n_columns,2*N+1);

    for k=1:(2*N+1)
        s=k-N-1;
        c_p(:,k)=diag(exp(-1i*s*const*t)'*R*sqrt(1-abs(s)/(N+1)));
    end

    var_cov_matrix=real((c_p*c_p')/(N+1)); % already includes conjugation
end