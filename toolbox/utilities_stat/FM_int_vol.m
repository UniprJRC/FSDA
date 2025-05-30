function V_int = FM_int_vol(x,t,T,varargin) 
%FM_int_vol computes the integrated variance of a diffusion process via the Fourier-Malliavin estimator  
%
%<a href="matlab: docsearchFS('FM_int_vol')">Link to the help function</a>
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
%   N     :   Cutting frequency. Scalar. If N is not specified, it is set
%              equal to floor(length(x)-/2).
%                 Example - 'N',400
%                 Data Types - single | double
%
% Output:
%
%  V_int :     Integrated variance. Scalar. Values of the integrated
%              variance.
%
%
% More About:
%
% We assume that the vector x contains discrete observations from a diffusion
% process following the Ito stochastic differential equation 
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$ 
% where $W$ is a Brownian motions defined on the filtered probability space 
% $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, while $\sigma$ and $b$ are random processes, adapted to $\mathcal{F}_t$.
% See the References for further  mathematical details.
% The integrated variance of the process $x$ on $[0,T]$ is defined as
% $$IV_{[0,T]}:=\int_0^T \sigma^2(t)\, dt.$$
% 
% For any positive integer $n$, let $\mathcal{S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times. Moreover, let $\delta_l(x):=
% x(t_{l+1})-x(t_l)$ be the increments of $x$. 
% The Fourier estimator of the integrated volatility on $[0,T]$ is
% given by
% 
% $$T c_0(\sigma_{{n},N})={T^2\over {2N+1}} \sum_{|k|\leq N} c_{k}(dx_{n})c_{-k}(dx_{n}),$$
% 
% where
% $$c_k(dx_{n})= {1\over {T}} \sum_{l=0}^{n-1} e^{-{\rm i}\frac{2\pi}{T}kt_l}\delta_{l}(x).$$
% 
% See also: FE_int_vol.m, FE_int_vol_Fejer.m,  OptimalCuttingFrequency.m, FM_int_cov.m, FM_int_quart.m, FM_int_volvol.m, FM_int_lev.m, Heston1D.m
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
%<a href="matlab: docsearchFS('FM_int_vol')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:  

%{ 
%% Example of call of FM_int_vol with default value of N.
% The following example estimates the integrated volatility of a
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


% Integrated volatility estimation
V_int=FM_int_vol(x,t,T);
   
%}

%{ 
%% Example of call of FM_int_vol with custom choice of N.
% The following example estimates the integrated volatility of a
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


% Integrated volatility estimation
V_int= FM_int_vol(x,t,T,'N',5000);


%}


 


%% Beginning of code

% Make sure that x1, x2, t1 and t2 are column vectors.

x=x(:);
t=t(:);

if length(x) ~= length(t)
    error('FSDA:FM_int_vol:WrongInputOpt','Input arguments x and t must have the same length.');
end


const=2*pi/T;
   
r=diff(x);  
  
        
n=length(r); 
N=floor(n/2);


if nargin>2
    options=struct('N',N);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FM_int_vol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FM_int_vol:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    N=options.N;
end

if N  >= n
    error('FSDA:FM_int_vol:WrongInputOpt','N must be strictly smaller than n.');
end
 
K=1:1:N; 
tt=-1i*const*t(1:end-1)';
        
c_ra=zeros(N,1); 

        
for j = 1:N  
    c_ra(j)= exp(K(j)*tt)*r;  
                  
end 
                
c_r=1/T* [ flip(conj(c_ra))  ; sum(r)  ; c_ra]; % Fourier coefficients of dx
        
V_int  = T^2/(2*N+1) * real(c_r.' * flip(c_r)) ;  
             
end

%FScategory:FMvol