function [V_spot,tau_out] = FM_spot_vol(x,t,T,varargin) 
%FM_spot_vol computes the spot volatility of a diffusion process via the Fourier-Malliavin estimator  
%
%<a href="matlab: docsearchFS('FM_spot_vol')">Link to the help function</a>
%  
% Required input arguments:
%
%   x     :    Observed process values. Vector. Row or column vector containing
%              the observed values.
%   t     :    Observation times. Vector.  Row or column vector with the same
%              length of x containing the observation times.
%   T      :   Estimation horizon. Scalar.  
%
% Optional input arguments:
%
%   N      :   Cutting frequency. Scalar. If N is not specified, it is set
%              equal to floor((length(x)-1)/2).
%                 Example - 'N',400
%                 Data Types - single | double
%   M      :   Cutting frequency. Scalar. If M is not specified, it is set
%              equal to floor(floor((length(x)-1)/2)^0.5).
%                 Example - 'M',20
%                 Data Types - single | double
%   tau    :   Estimation times. Vector. If tau is not specified, it is set
%                equal to 0:T/(2*M):T.
%                 Example - 'tau', 0:T/100:T 
%                 Data Types - single | double
%
% Output:
%
% V_spot :     Spot variance estimates. Vector. Estimated values of the spot
%              variance of the process.
% tau_out :    Estimation times. Vector. Coincides with the input vector tau unless the length 
%              of the latter is larger than 2M+1.
% More About:
%
% We assume that the vector x contains discrete observations from a diffusion
% process $x$ following the Ito stochastic differential equation 
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$ 
% where $W$ is a Brownian motion defined on the filtered probability space $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, 
% while $\sigma$ and $b$ are random processes, adapted to  $\mathcal{F}_t$.
% See the Reference for further  mathematical details.
% The spot volatility of the process $x$ at time $t \in [0,T]$ is defined as $V(t):=\sigma^2(t)$.
% 
% For any positive integer $n$, let $\mathcal{S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times. Moreover, let $\delta_l(x):=
% x(t_{l+1})-x(t_l)$ be the increments of $x$. 
% The Fourier estimator of the spot volatility at time $t \in [0,T]$ is
% given by 
% $$\widehat V_{n,N,M}(\tau)= \sum_{|k|\leq M} \left(1-{|k|\over
% M}\right)c_k(\sigma_{n,N}) \, e^{{\rm i}\frac{2\pi}{T}k\tau},$$
% where:
% $$c_k(\sigma_{n,N})={T\over {2N+1}} \sum_{|s|\leq N} c_{s}(dx_{n})c_{k-s}(dx_{n}),$$
%
% $$c_k(dx_{n})= {1\over {T}} \sum_{l=0}^{n-1} e^{-{\rm i}\frac{2\pi}{T}kt_l}\delta_{l}(x).$$
% 
% See also: FE_spot_vol.m, FE_spot_vol_FFT.m,  FM_spot_quart.m, FM_spot_volvol.m, FM_spot_lev.m, Heston1D.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin Volatility Estimation. Theory and Practice, "Springer Briefs in Quantitative Finance", Springer. 
% 
% Sanfelici, S., Toscano, G. (2024), The Fourier-Malliavin Volatility (FMVol) MATLAB toolbox, available on ArXiv.
%
%
%
%
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FM_spot_vol')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:  


%{ 
%% Example of call of FM_spot_vol with default values of M,N and tau.
% The following example estimates the path of the spot volatility of a
% Heston model from a discrete sample. The Heston model assumes that the spot variance follows 
% a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0,0.4,2,1];
Rho=-0.5;
x0=log(100); 
V0=0.4;
[x,V,t]=Heston1D(T,n,parameters,Rho,x0,V0); 


% Spot volatility estimation 
[V_spot, tau_out]=FM_spot_vol(x,t,T);
M=(length(V_spot)-1)/2;   

figure
plot(tau_out,V(1:round(n/(2*M)):end));
hold on
plot(tau_out,V_spot); 
xlabel('tau'); 
title('Spot volatility estimates Vs Actual values')
legend('Actual values','Estimated values')
 
%}

%{ 
%% Example of call of FM_spot_vol with custom choices of M,N and tau.
% The following example estimates the path of the spot volatility of a
% Heston model from a discrete sample. The Heston model assumes that the spot variance follows 
% a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1; 
n=23400;  
parameters=[0,0.4,2,1];
Rho=-0.5;
x0=log(100); 
V0=0.4;
[x,V,t]=Heston1D(T,n,parameters,Rho,x0,V0); 


% Spot volatility estimation
tau=0:T/50:T;
[V_spot, tau_out]= FM_spot_vol(x,t,T,'N',5000,'M',100,'tau',tau);

figure
plot(tau_out,V(1:round(n/50):end));
hold on
plot(tau_out,V_spot); 
xlabel('tau'); 
title('Spot volatility estimates Vs Actual values')
legend('Actual values','Estimated values')
%}

%{ 
%% Example of call of FM_spot_vol when tau has length larger than 2M &plus; 1.
% The following example estimates the path of the spot volatility of a
% Heston model from a discrete sample. The Heston model assumes that the spot variance follows 
% a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1; 
n=23400;  
parameters=[0,0.4,2,1];
Rho=-0.5;
x0=log(100); 
V0=0.4;
[x,V,t]=Heston1D(T,n,parameters,Rho,x0,V0); 


% Spot volatility estimation
tau=0:T/1000:T;
[V_spot, tau_out]= FM_spot_vol(x,t,T,'N',5000,'M',100,'tau',tau);

figure
M=100;
plot(tau_out,V(1:round(n/(2*M)):end));
hold on
plot(tau_out,V_spot); 
xlabel('tau'); 
title('Spot volatility estimates Vs Actual values')
legend('Actual values','Estimated values')
%}


%% Beginning of code

% Make sure that x and t are column vectors.

x=x(:);
t=t(:);


if length(x) ~= length(t)
    error('FSDA:FM_spot_vol:WrongInputOpt','Input arguments x and t must have the same length.');
end



const=2*pi/T;
   
r=diff(x);  

        
n=length(r);
N=floor(n/2);
M=floor(N^0.5);   
tau=0:T/(2*M):T;


if nargin>2
    options=struct('N',N,'M',M,'tau',tau);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FM_spot_vol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FM_spot_vol:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
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

if N  >= n
    error('FSDA:FM_spot_vol:WrongInputOpt','N must be strictly smaller than n.');
end

if M  >= N
    error('FSDA:FM_spot_vol:WrongInputOpt','M must be strictly smaller than N.');
end

if length(tau) > 2*M+1
    disp('WARNING: estimation will be performed on the equally-spaced grid with mesh size equal to T/(2*M), provided as an output variable.')
    tau_out=0:T/(2*M):T;
else 
    tau_out=tau;
end


k=1:1:N+M; 
tt=-1i*const*t(1:end-1)';  
        
c_ra=zeros(length(k),1); 

        
for j = 1:length(k)  
    c_ra(j)= exp(k(j)*tt)*r;           
end 
                
c_r=1/T* [ flip(conj(c_ra))  ; sum(r)  ; c_ra]; % Fourier coefficients of dx

        
        
c_v=zeros(2*M+1,1);
c_r_aux=zeros(2*M+1,2*N+1); 
c_v_spot=zeros(2*M+1,length(tau_out)); 
       

center=N+M+1;

for j = 1 : 2*M+1        
    c_r_aux(j,1:2*N+1)=   c_r  (  center -N +( j-M-1 ) : center +N +( j-M-1 ) ) ;
    c_v(j)= T/(2*N+1) * c_r  (  center -N   : center +N   ).' * flip(c_r_aux(j, : )).';
    for ii=1 : length(tau_out)
          c_v_spot(j,ii)=c_v(j)*(1-abs(j-M-1)/(M+1))*exp(const*tau_out(ii)*(j-M-1)*1i);  % spot variance of process x
    end         
end
               
V_spot = real(sum(c_v_spot));
end

%FScategory:FMvol