function C_int = FM_int_cov(x1,x2,t1,t2,T,varargin) 
%FM_int_cov computes the integrated covariance of a bivariate diffusion process via the Fourier-Malliavin estimator  
%
%<a href="matlab: docsearchFS('FM_int_cov')">Link to the help function</a>
%  
% Required input arguments:
%
%   x1     :   Observation values of process 1. Vector. Row or column vector containing
%              the observed values.
%   x2     :   Observation values of process 2. Vector. Row or column vector containing
%              the observed values.
%   t1     :   Observation times of process 1. Vector. Row or column vector containing
%              the observation times.
%   t2     :   Observation times of process 2. Vector. Row or column vector containing
%              the observation times.
%   T      :   Estimation horizon. Scalar.  
%
% Optional input arguments:
%
%   N      :   Cutting frequency. Scalar. If N is not specified, it is set
%              equal to floor((min(length(x1),length(x2))-1)/2).
%                 Example - 'N',400
%                 Data Types - single | double
%
% Output:
%
%  C_int  :    Integrated covariance of processes 1 and 2 on [0,T]. Scalar. Value of the integrated
%              covariance of processes 1 and 2.
%
%
% More About:
%
% We assume that vectors x1 and x2 contain discrete observations from a bivariate diffusion
% process $(x_1,x_2)$ following the Ito stochastic differential equation 
% $$dx_i(t)= \sigma_i(t) \ dW_i(t) + b_i(t) \ dt, \quad i=1,2,$$ 
% where $W_1$ and $W_2$ are two Brownian motions defined on the filtered probability space 
% $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, with correlation $\rho$, 
% while $\sigma_1, \sigma_2, b_1$ and $b_2$ are random processes, adapted to  $\mathcal{F}_t$.
% See the References for further  mathematical details.
% The integrated covariance of the process $(x_1,x_2)$ on $[0,T]$ is defined as
% $$IC_{[0,T]}:=\rho \, \int_0^T \sigma_1(t) \sigma_2(t) \, dt$$
% 
% Let $i=1,2$. For any positive integer $n_i$, let $\mathcal{S}^i_{n_i}:=\{ 0=t^i_{0}\leq \cdots
% \leq t^i_{n_i}=T  \}$ be the observation times for the $i$-th asset. Moreover, let $\delta_l(x^i):=
% x^i(t^i_{l+1})-x^i(t^i_l)$ be the increments of $x^i$. 
% The Fourier estimator of the integrated covariance on $[0,T]$ is
% given by 
% $$T c_0(c_{n_1,n_2,N})={T\over {2N+1}} \sum_{|k|\leq N} c_{k}(dx^1_{n_1})c_{-k}(dx^2_{n_2}),$$
%
% $$c_k(dx^i_{n_i})= {1\over {T}} \sum_{l=0}^{n_i-1} e^{-{\rm i}\frac{2\pi}{T}kt^i_l}\delta_{l}(x_i).$$
% 
% See also: FE_int_vol.m, FE_int_vol_Fejer.m, FM_int_vol, FM_int_quart.m, FM_int_volvol.m, FM_int_lev.m, Heston2D.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin Volatility Estimation. Theory and Practice, "Springer Briefs in Quantitative Finance", Springer. 
% 
% Sanfelici, S., Toscano, G. (2024), The Fourier-Malliavin Volatility (FMVol) MATLAB toolbox, available on ArXiv.
%
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FM_int_cov')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:  

%{ 
    %% Example of call of FM_int_cov with default value of N.
    % The following example estimates the daily integrated covariance of a bivariate 
    % Heston model from a discrete sample. The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1; 
n=23400; 
parameters=[0,0;0.4,0.4;2,2;1,1];
Rho=[0.5,-0.5,0,0,-0.5,0.5];
x0=[log(100); log(100)]; 
V0=[0.4; 0.4];
[x,V,t]=Heston2D(T,n,parameters,Rho,x0,V0); 


% Integrated covariance estimation 
t1=t; 
t2=t;
C_int=FM_int_cov(x(:,1),x(:,2),t1,t2,T);
 
%}

%{ 
    %% Example of call of FM_int_cov with custom choice of N.
    % The following example estimates the integrated covariance of a bivariate 
    % Heston model from a discrete sample. The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1; % horizon of the trajectory (in days)
n=23400; % number of observations simulated in one trajectory 
parameters=[0 , 0; 0.4 , 0.4; 2 , 2; 1 , 1 ];
Rho=[0.5 ; -0.5 ; 0 ; 0 ; -0.5 ; 0.5];
x0=[log(100); log(100)]; V0=[0.4; 0.4];
[x,V,t] = Heston2D(T,n,parameters,Rho,x0,V0); 


% Integrated covariance estimation 
t1=t; 
t2=t;
C_int=FM_int_cov(x(:,1),x(:,2),t1,t2,T,'N',5000);


%}

%% Beginning of code

% Make sure that x1, x2, t1 and t2 are column vectors.

x1=x1(:);
x2=x2(:);
t1=t1(:);
t2=t2(:);

if length(x1) ~= length(t1)
    error('FSDA:FM_int_cov:WrongInputOpt','Input arguments x1 and t1 must have the same length.');
end

if length(x2) ~= length(t2)
    error('FSDA:FM_int_cov:WrongInputOpt','Input arguments x2 and t2 must have the same length.');
end

const=2*pi/T;
   
r1=diff(x1);  
r2=diff(x2);  
        
n1=length(r1); n2=length(r2);
N=floor((min(n1,n2)-1)/2);

if nargin>2
    options=struct('N',N);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FM_int_cov:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FM_int_cov:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end  
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    N=options.N;
end

if N  >= min(n1,n2)
    error('FSDA:FM_int_cov:WrongInputOpt','N must be strictly smaller than min(n1,n2).');
end
 
K=1:1:N; 
tt1=-1i*const*t1(1:end-1)'; tt2=-1i*const*t2(1:end-1)';
        
c_r1a=zeros(N,1); 
c_r2a=zeros(N,1); 
        
for j = 1:N  
    c_r1a(j)= exp(K(j)*tt1)*r1;  
    c_r2a(j)= exp(K(j)*tt2)*r2;              
end 
                
c_r1=1/T* [ flip(conj(c_r1a))  ; sum(r1)  ; c_r1a];  
c_r2=1/T* [ flip(conj(c_r2a))  ; sum(r2)  ; c_r2a];  

c_r1b=zeros(2*N+1,1);
c_r2b=zeros(2*N+1,1);

for j = 1 : 2*N+1
c_r1b(j)= (1 - abs(j-N-1 )/(N+1)) *c_r1 (j) ;
c_r2b(j)= (1 - abs(j-N-1 )/(N+1))*c_r2 (j)  ;      
end

C_int  = T^2/(N+1) * real(c_r2b.' * flip(c_r1)) ;  


end

%FScategory:FMvol