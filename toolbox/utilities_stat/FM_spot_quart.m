function [Q_spot, tau_out] = FM_spot_quart(x,t,T,varargin) 
%FM_spot_quart computes the spot quarticity of a diffusion process via the Fourier-Malliavin estimator  
%
%<a href="matlab: docsearchFS('FM_spot_quart')">Link to the help function</a>
%  
% Required input arguments:
%
%   x     :   Observed process values. Vector. Row or column vector containing
%             the observed values.
%   t     :   Observation times. Vector.  Row or column vector with the same
%             length of x containing the observation times.
%   T     :   Estimation horizon. Scalar.  
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
%   L      :   Cutting frequency. Scalar. If L is not specified, it is set
%              equal to floor(floor(floor((length(x)-1)/2)^0.5)^0.5).
%                 Example - 'L',5
%                 Data Types - single | double
%   tau    :   Estimation times. Vector. If tau is not specified, it is set
%                equal to 0:T/(2*L):T.
%                 Example - 'tau', 0:T/100:T 
%                 Data Types - single | double
%
% Output:
%
%  Q_spot  :    Spot quarticity estimates. Vector. Estimated values of the spot
%               quarticity of the process.
%  tau_out :    Estimation times. Vector. Coincides with the input vector tau unless the length 
%                of the latter is larger than 2L+1.
% More About:
%
% We assume that the vectors x contains discrete observations from a diffusion
% process $x$ following the Ito stochastic differential equation 
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$  
% where $W$ is a Brownian motions defined on the filtered probability space 
% $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, while  
% $\sigma$ and $b$ are random processes, adapted to  $\mathcal{F}_t$.
% See the References for further  mathematical details.
% The spot quarticity at time $t \in [0,T]$ is defined as $Q(t):=\sigma^4(t).$
% 
% 
% For any positive integer $n$, let $\mathcal{S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times. Moreover, let $\delta_l(x):=
% x(t_{l+1})-x(t_l)$ be the increments of $x$. 
% The Fourier estimator of the spot quarticity at time $t \in [0,T]$ is
% given by
% $$\widehat Q_{n,N,M,L}(\tau)= \sum_{|k|\leq L} \left(1-{|k|\over
% L}\right)c_k(Q_{n,N,M}) \, e^{{\rm i}\frac{2\pi}{T}k\tau},$$
% where:
% $$c_k(Q_{n,N,M})= {T \over {2M+1}} \sum_{|s|\leq M} c_s( \sigma_{{n},N})c_{k-s}(\sigma_{{n},N}),$$
%
% $$c_k(\sigma_{n,N})={T\over {2N+1}} \sum_{|s|\leq N} c_{s}(dx_{n})c_{k-s}(dx_{n}),$$
%
% $$c_k(dx_{n})= {1\over {T}} \sum_{l=0}^{n-1} e^{-{\rm i}\frac{2\pi}{T}kt_l}\delta_{l}(x).$$
% 
% See also: FM_spot_vol.m, FM_spot_volvol.m, FM_spot_lev.m, Heston1D.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin Volatility Estimation. Theory and Practice, "Springer Briefs in Quantitative Finance", Springer. 
% 
% Sanfelici, S., Toscano, G. (2024), The Fourier-Malliavin Volatility (FMVol) MATLAB toolbox, available on ArXiv.
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FM_spot_quart')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:  

%{
    %% Example of call of FM_spot_quart with default values of N,M,L and tau.
    % The following example estimates the path of the spot quarticity
    % of a random process following the Heston model from a discrete sample. 
    % The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0,0.4,5,1];
rho=-0.5;
x0=log(100); 
V0=0.4;
[x,V,t] = Heston1D(T,n,parameters,rho,x0,V0);


% Spot qaurticity estimation 
[Q_spot, tau_out] = FM_spot_quart(x,t,T);
L=(length(Q_spot)-1)/2;
 
figure
Q=V.^2;
plot(tau_out,Q(1:round(n/(2*L)):end));
hold on
plot(tau_out,Q_spot);  
title('Spot quarticity estimates Vs Actual values')
legend('Actual values','Estimated values')
 
%}
%{
    %% Example of call of FM_spot_quart with custom choices of N,M,L and tau.
    % The following example estimates the path of the spot quarticity
    % of a random process following the Heston model from a discrete sample. 
    % The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0,0.4,5,1];
rho=-0.5;
x0=log(100); 
V0=0.4;
[x,V,t]=Heston1D(T,n,parameters,rho,x0,V0);

% Spot quarticity estimation  
tau=0:T/10:T;
[Q_spot,tau_out]= FM_spot_quart(x,t,T,'N',11000,'M',120,'L',15,'tau',tau);

figure
Q=V.^2;
plot(tau_out,Q(1:round(n/10):end));
hold on
plot(tau_out,Q_spot); 
xlabel('tau'); 
title('Spot quarticity estimates Vs Actual values')
legend('Actual values','Estimated values')

 
%}

%{
    %% Example of call of FM_spot_quart when tau has length larger than 2L &plus; 1.
    % The following example estimates the path of the spot quarticity
    % of a random process following the Heston model from a discrete sample. 
    % The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0,0.4,5,1];
rho=-0.5;
x0=log(100); 
V0=0.4;
[x,V,t]=Heston1D(T,n,parameters,rho,x0,V0);

% Spot quarticity estimation  
tau=0:T/500:T;
[Q_spot,tau_out]= FM_spot_quart(x,t,T,'N',11000,'M',120,'L',10,'tau',tau);

figure
Q=V.^2;
L=10;
plot(tau_out,Q(1:round(n/(2*L)):end));
hold on
plot(tau_out,Q_spot); 
xlabel('tau'); 
title('Spot quarticity estimates Vs Actual values')
legend('Actual values','Estimated values')
%}

%% Beginning of code

% Make sure that x and t are column vectors.

x=x(:);
t=t(:);

if length(x) ~= length(t)
    error('FSDA:FM_spot_quart:WrongInputOpt','Input arguments x and t must have the same length.');
end


const=2*pi/T;

r=diff(x);  
  
        
n=length(r);  
N=floor(n/2);
M=floor(N^0.5); 
L=floor(M^0.5);

tau=0:T/(2*L):T;

if nargin>2
    options=struct('N',N,'M',M,'L',L,'tau',tau);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FM_spot_quart:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FM_spot_quart:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    N=options.N;
    M=options.M;
    L=options.L;
    tau=options.tau;
end

if N  >= n
    error('FSDA:FM_spot_quart:WrongInputOpt','N must be strictly smaller than min(n1,n2).');
end

if M  >= N
    error('FSDA:FM_spot_quart:WrongInputOpt','M must be strictly smaller than N.');
end

if L  >= M
    error('FSDA:FM_spot_quart:WrongInputOpt','L must be strictly smaller than M.');
end        

if length(tau) > 2*L+1
    disp('WARNING: estimation will be performed on the equally-spaced grid with mesh size equal to T/(2*L), provided as an output variable.')
    tau_out=0:T/(2*L):T;
else 
    tau_out=tau;
end

k=1:1:N+M+L; 
tt=-1i*const*t(1:end-1)';  
 

c_ra=zeros(length(k),1); 
         
for j = 1:length(k)  
    c_ra(j)= exp(k(j)*tt)*r;  
 end 
                
c_r1=1/T* [ flip(conj(c_ra))  ; sum(r)  ; c_ra];  
         
         
c_v=zeros(2*M+1,1);  
c_r_aux1=zeros(2*M+1,2*N+1);       
 
center=N+M+L+1;
for j = 1 : 2*M+1        
    c_r_aux1(j,1:2*N+1)=   c_r1  (  center -N +( j-M-1 ) : center +N +( j-M-1 ) ) ;
    c_v(j)= T/(2*N+1) * c_r1  (  center -N   : center +N   ).' * flip(c_r_aux1(j, : )).';
end

MM=M+L;
 
c_v2=zeros(2*MM+1,1); 
c_r_aux2=zeros(2*MM+1,2*N+1); 

 

for j = 1 : 2*MM+1    
c_r_aux2(j,1:2*N+1)= c_r1(center-N+(j-MM-1) : center+N+(j-MM-1));    
c_v2(j)= T/(2*N+1)*c_r1(center-N : center+N).'*flip(c_r_aux2(j, : )).';
end

c_q=zeros(2*L+1,1); c_q_spot=zeros(2*L+1,length(tau_out)); c_v3=zeros(2*L+1, 2*M +1);
 
center2=MM+1;

for j = 1 : 2*L+1
c_v3(j,1:2*M +1)= c_v2 (  center2 -M  +(j-L-1 ) : center2 +M  +(j-L-1 ) );
c_q(j)=  T  *  c_v.'*flip(c_v3(j,:)).';

  for ii=1 : length(tau_out)
  c_q_spot(j,ii)=c_q(j)*(1-abs(j-L-1)/L)*exp(const*tau_out(ii)*(j-L-1)*1i);
  end
end
 
Q_spot=real(sum(c_q_spot));  
 
end
 
%FScategory:FMvol
