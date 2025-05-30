function [C_spot, tau_out] = FM_spot_cov(x1,x2,t1,t2,T,varargin) 
%FM_spot_cov computes the spot covariance of a bivariate diffusion process via the Fourier-Malliavin estimator  
%
%<a href="matlab: docsearchFS('FM_spot_cov')">Link to the help function</a>
%  
% Required input arguments:
%
%   x1     :   Observation values of process 1. Vector. Row or column vector containing
%              the observed values.
%   x2     :   Observation values of process 2. Vector. Row or column vector containing
%              the observed values.
%   t1     :   Observation times of process 1. Vector.  Row or column vector with the same
%              length of x1 containing the observation times.
%   t2     :   Observation times of process 2. Vector.  Row or column vector with the same
%              length of x2 containing the observation times.
%   T      :   Estimation horizon. Scalar.  
%
% Optional input arguments:
%
%   N      :   Cutting frequency. Scalar. If N is not specified, it is set
%              equal to floor(min(length(x1),length(x2))-1)/2).
%                 Example - 'N',400
%                 Data Types - single | double
%   M      :   Cutting frequency. Scalar. If M is not specified, it is set
%              equal to floor(floor(min(length(x1),length(x2))-1)/2)^0.5).
%                 Example - 'M',20
%                 Data Types - single | double
%   tau    :   Estimation times. Vector. If tau is not specified, it is set
%                equal to 0:T/(2*M):T.
%                 Example - 'tau', 0:T/100:T 
%                 Data Types - single | double
%
% Output:
%
%  C_spot  :    Spot covariance of processes 1 and 2. Vector. Values of the
%               spot covariance of processes 1 and 2.
%
% tau_out :    Estimation times. Vector. Coincides with the input vector tau unless the length 
%              of the latter is larger than 2M+1.
%
%
% More About:
%
% We assume that vectors x1 and x2 contain discrete observations from a bivariate diffusion
% process $(x_1,x_2)$ following the Ito stochastic differential equation 
% $$dx^i(t)= \sigma^i(t) \ dW^i(t) + b^i(t) \ dt, \quad i=1,2,$$ 
% where $W^1$ and $W^2$ are two Brownian motions defined on the filtered probability space 
% $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, with correlation $\rho$, while  
% $\sigma^1, \sigma^2, b^1$ and $b^2$ are random processes, adapted to  $\mathcal{F}_t$.
% See the References for further  mathematical details.
% The spot covariance $c$ between the processes $x^1$ and $x^2$ at time $t
% \in [0,T]$ is defined as
% $$c(t):= \frac{d \langle x^1, x^2 \rangle_t}{dt} =  \rho \sigma^1(t) \sigma^2(t).$$
% 
% Let $i=1,2$. For any positive integer $n_i$, let $\mathcal{S}^i_{n_i}:=\{ 0=t^i_{0}\leq \cdots
% \leq t^i_{n_i}=T  \}$ be the observation times for the $i$-th asset. Moreover, let $\delta_l(x^i):=
% x^i(t^i_{l+1})-x^i(t^i_l)$ be the increments of $x^i$. 
% The Fourier estimator of the spot covariance at time $t \in [0,T]$ is
% given by 
% $$\widehat c_{n_1,n_2,N,M}(\tau)= \sum_{|k|\leq M} \left(1-{|k|\over
% M}\right)c_k(c_{n_1,n_2,N}) \, e^{{\rm i}\frac{2\pi}{T}k\tau},$$
% where:
% $$c_k(c_{n_1,n_2,N})={T\over {2N+1}} \sum_{|s|\leq N} c_{s}(dx^1_{n_1})c_{k-s}(dx^2_{n_2}),$$
%
% $$c_k(dx^i_{n_i})= {1\over {T}} \sum_{l=0}^{n_i-1} e^{-{\rm i}\frac{2\pi}{T}kt^i_l}\delta_{l}(x_i).$$
% 
% See also: FE_spot_vol.m, FE_spot_vol_FFT.m, FM_spot_vol.m, FM_spot_quart.m, FM_spot_volvol.m, FM_spot_lev.m, Heston2D.m
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
%<a href="matlab: docsearchFS('FM_spot_cov')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:  

%{ 
    %% Example of call of FM_spot_cov with default values of N, M and tau.
    % The following example estimates the path of the spot covariance  of a bivariate 
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


% Spot covariance estimation
t1=t; 
t2=t;
[C_spot, tau_out]=FM_spot_cov(x(:,1),x(:,2),t1,t2,T);
M=(length(C_spot)-1)/2;   

figure
Cov=Rho(6)*sqrt(V(:,1)).*sqrt(V(:,2));
plot(tau_out,Cov(1:round(n/(2*M)):end));
hold on
plot(tau_out,C_spot); 
xlabel('tau'); 
title('Spot co-volatility estimates Vs Actual values')
legend('Actual values','Estimated values')

 
%}

%{ 
    %% Example of call of FM_spot_cov with custom choices of N, M and tau.
    % The following example estimates the path of the spot covariance of a bivariate 
    % Heston model from a discrete sample. The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0 , 0; 0.4 , 0.4; 2 , 2; 1 , 1 ];
Rho=[0.5 ; -0.5 ; 0 ; 0 ; -0.5 ; 0.5];
x0=[log(100); log(100)]; V0=[0.4; 0.4];
[x,V,t] = Heston2D(T,n,parameters,Rho,x0,V0); 



% Estimation of the spot covariance 
t1=t; 
t2=t;
tau=0:T/50:T;
[C_spot,tau_out]=FM_spot_cov(x(:,1),x(:,2),t1,t2,T,'N',5000,'M',100,'tau',tau);

figure
Cov=Rho(6)*sqrt(V(:,1)).*sqrt(V(:,2));
plot(tau_out,Cov(1:round(n/50):end));
hold on
plot(tau_out,C_spot); 
xlabel('tau'); 
title('Spot co-volatility estimates Vs Actual values')
legend('Actual values','Estimated values')

 

%}

%{ 
    %% Example of call of FM_spot_cov when tau has length larger than 2M &plus; 1.
    % The following example estimates the path of the spot covariance of a bivariate 
    % Heston model from a discrete sample. The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0 , 0; 0.4 , 0.4; 2 , 2; 1 , 1 ];
Rho=[0.5 ; -0.5 ; 0 ; 0 ; -0.5 ; 0.5];
x0=[log(100); log(100)]; V0=[0.4; 0.4];
[x,V,t] = Heston2D(T,n,parameters,Rho,x0,V0); 



% Estimation of the spot covariance 
t1=t; 
t2=t;
tau=0:T/1000:T;
[C_spot,tau_out]=FM_spot_cov(x(:,1),x(:,2),t1,t2,T,'N',5000,'M',100,'tau',tau);

figure
Cov=Rho(6)*sqrt(V(:,1)).*sqrt(V(:,2));
M=100;
plot(tau_out,Cov(1:round(n/(2*M)):end));
hold on
plot(tau_out,C_spot); 
xlabel('tau'); 
title('Spot co-volatility estimates Vs Actual values')
legend('Actual values','Estimated values')

 

%} 

%% Beginning of code

% Make sure that x1, x2, t1 and t2 are column vectors.

x1=x1(:);
x2=x2(:);
t1=t1(:);
t2=t2(:);

if length(x1) ~= length(t1)
    error('FSDA:FM_spot_cov:WrongInputOpt','Input arguments x1 and t1 must have the same length.');
end

if length(x2) ~= length(t2)
    error('FSDA:FM_spot_cov:WrongInputOpt','Input arguments x2 and t2 must have the same length.');
end

const=2*pi/T;
   
r1=diff(x1);  
r2=diff(x2);  
        
n1=length(r1); n2=length(r2);
N=floor(min(n1,n2)/2);
M=floor(N^0.5);   
tau=0:T/(2*M):T;

if nargin>2
    options=struct('N',N,'M',M,'tau',tau);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FM_spot_cov:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FM_spot_cov:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    N=options.N;
    M=options.M;
    tau=options.tau;
end

if N  >= min(n1,n2)
    error('FSDA:FM_spot_cov:WrongInputOpt','N must be strictly smaller than min(n1,n2).');
end

if M  >= N
    error('FSDA:FM_spot_cov:WrongInputOpt','M must be strictly smaller than N.');
end


if length(tau) > 2*M+1
    disp('WARNING: estimation will be performed on the equally-spaced grid with mesh size equal to T/(2*M), provided as an output variable.')
    tau_out=0:T/(2*M):T;
else 
    tau_out=tau;
end

k=1:1:N+M; 
tt1=-1i*const*t1(1:end-1)'; tt2=-1i*const*t2(1:end-1)';
        
c_r1a=zeros(length(k),1); 
c_r2a=zeros(length(k),1); 
        
for j = 1:length(k)  
    c_r1a(j)= exp(k(j)*tt1)*r1;  
    c_r2a(j)= exp(k(j)*tt2)*r2;              
end 
                
c_r1=1/T* [ flip(conj(c_r1a))  ; sum(r1)  ; c_r1a];  
c_r2=1/T* [ flip(conj(c_r2a))  ; sum(r2)  ; c_r2a]; 
        
         
c_c=zeros(2*M+1,1);   
c_r_aux1=zeros(2*M+1,2*N+1);  
c_c_spot=zeros(2*M+1,length(tau_out)); 
       

center=N+M+1;

for j = 1 : 2*M+1        
    c_r_aux1(j,1:2*N+1)=   c_r1  (  center -N +( j-M-1 ) : center +N +( j-M-1 ) ) ;
    c_c(j)= T/(2*N+1) * c_r2  (  center -N   : center +N   ).' * flip(c_r_aux1(j, : )).';
    for ii=1 : length(tau_out)
          c_c_spot(j,ii)=c_c(j)*(1-abs(j-M-1)/(M+1))*exp(const*tau_out(ii)*(j-M-1)*1i);  
    end         
end
               
C_spot=real(sum(c_c_spot));  

end

%FScategory:FMvol