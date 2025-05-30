function VV_int = FM_int_volvol(x,t,T,varargin) 
%FM_int_volvol computes the integrated volatility of volatility of a diffusion process via the Fourier-Malliavin estimator  
%
%<a href="matlab: docsearchFS('FM_int_volvol')">Link to the help function</a>
%  
% Required input arguments:
%
%   x     :    Observed process values. Vector. Row or column vector containing
%              the observed values.
%   t     :    Observation times. Vector.  Row or column vector with the same
%              length of x containing the observation times.
%   T     :    Estimation horizon. Scalar. 
%
% Optional input arguments:
%
%   N      :   Cutting frequency. Scalar. If N is not specified, it is set
%              equal to floor((length(x)-1)/2).
%                 Example - 'N',400
%                 Data Types - single | double
%   M      :   Cutting frequency. Scalar. If M is not specified, it is set
%              equal to floor(floor((length(x)-1)/2)^0.4).
%                 Example - 'M',20
%                 Data Types - single | double
%                 Data Types - single | double
%
% Output:
%
%  VV_int  :    Integrated volatility of volatility. Scalar. Value of the integrated volatility of volatility.
%
%
% More About:
%
% We assume that the vectors x contains discrete observations from a diffusion
% process $x$ following the Ito stochastic differential equation 
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$ 
% $$d\sigma^2(t)= \gamma(t) \ dZ(t) + a(t) \ dt,$$ 
% where $W$ and $Z$ are two Brownian motions defined on the filtered probability space 
% $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, with correlation $\rho$, while 
% $\sigma, \gamma, b$ and $a$ are random processes, adapted to  $\mathcal{F}_t$.
% See the References for further  mathematical details.
% The integrated volatility of volatility on $[0,T]$ is defined as 
% $$ IVV_{[0,T]}:= \langle  \sigma^2 ,\sigma^2 \rangle_T =\int_0^T\gamma^2(s)ds.$$
% 
% 
%For any positive integer $n$, let $\mathcal{S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times. Moreover, let $\delta_l(x):=
% x(t_{l+1})-x(t_l)$ be the increments of $x$. 
% The Fourier estimator of the integrated volatility of volatility on $[0,T]$ is
% given by
% $$T c_0(V_{{n},N,M})={T^2\over {2M+1}} \sum_{|k|\leq M} c_{k}(d\sigma_{n,N})c_{-k}(d\sigma_{n,N}),$$
% 
% where:
% $$c_k(d\sigma_{n,N})= i k \frac{2\pi}{T} c_k(\sigma_{n,N}), \quad   c_k(\sigma_{{n},N})={T\over {2N+1}} \sum_{|s|\leq N} c_{s}(dx_{n})c_{k-s}(dx_{n}),$$
% $$c_k(dx_{n})= {1\over {T}} \sum_{l=0}^{n-1} e^{-{\rm i}\frac{2\pi}{T}kt_l}\delta_{l}(x).$$
% 
% See also: FM_int_vol.m, FM_int_cov.m, FM_int_quart.m, FM_int_lev.m, Heston1D.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin Volatility Estimation. Theory and Practice, "Springer Briefs in Quantitative Finance", Springer. 
% 
% Sanfelici, S., Toscano, G. (2024), The Fourier-Malliavin Volatility (FMVol) MATLAB toolbox, available on ArXiv.
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FM_int_volvol')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Example:  

%{
    %% Example of call of FM_int_volvol with default values of N and M.
    % The following example estimates the integrated volatility of volatility
    % of a random process following the Heston model from a discrete sample. 
    % The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0,0.8,10,3.25];
rho=-0.3;
x0=log(100); 
V0=0.8;
[x,V,t] = Heston1D(T,n,parameters,rho,x0,V0);
 

% Integrated volatility of volatility estimation
VV_int = FM_int_volvol(x,t,T);
 
%}

%{
    %% Example of call of FM_int_volvol with custom choices of N and M.
    % The following example estimates the integrated volatility of volatility
    % of a random process following the Heston model from a discrete sample. 
    % The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1; 
n=23400; % number of observations simulated in one trajectory 
parameters=[0,0.8,10,3.25];
rho=-0.3;
x0=log(100); 
V0=0.8;
[x,V,t]=Heston1D(T,n,parameters,rho,x0,V0);
 

% Integrated volatility of volatility estimation 
VV_int = FM_int_volvol(x,t,T,'N',11200, 'M', 50);



 
%}


%% Beginning of code

% Make sure that x and t are column vectors.

x=x(:);
t=t(:);

if length(x) ~= length(t)
    error('FSDA:FM_int_volvol:WrongInputOpt','Input arguments x and t must have the same length.');
end


const=2*pi/T;

r=diff(x);  
  
        
n=length(r);  
N=floor(n/2);
M=floor(N^0.4); 

if nargin>2
    options=struct('N',N,'M',M);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FM_int_volvol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FM_int_volvol:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    N=options.N;
    M=options.M;
end

if N  >= n
    error('FSDA:FM_int_volvol:WrongInputOpt','N must be strictly smaller than min(n1,n2).');
end

if M  >= N
    error('FSDA:FM_int_volvol:WrongInputOpt','M must be strictly smaller than N.');
end

 

k=1:1:N+M; 
tt=-1i*const*t(1:end-1)';  
 

c_ra=zeros(length(k),1); 
         
for j = 1:length(k)  
    c_ra(j)= exp(k(j)*tt)*r;  
 end 
                
c_r1=1/T* [ flip(conj(c_ra))  ; sum(r)  ; c_ra]; % Fourier coefficients of dx1
         
         
c_v=zeros(2*M+1,1);  
c_r_aux1=zeros(2*M+1,2*N+1);  
c_dv=zeros(2*M+1,1);       
c_dv2=zeros(2*M+1,1);

center=N+M+1;

for j = 1 : 2*M+1        
    c_r_aux1(j,1:2*N+1)=   c_r1  (  center -N +( j-M-1 ) : center +N +( j-M-1 ) ) ;
    c_v(j)= T/(2*N+1) * c_r1  (  center -N   : center +N   ).' * flip(c_r_aux1(j, : )).';
    c_dv(j)=const*1i*(j-M-1) *(1-abs(j-M-1)/(M+1))*c_v(j); 
    c_dv2(j)=const*1i*(j-M-1) *c_v(j);
end
 
VV_int= T^2/(M+1)*   real(  c_dv.'*flip(c_dv2)) ;

end

%FScategory:FMvol