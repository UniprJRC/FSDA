function L_int = FM_int_lev(x,t,T,varargin) 
%FM_int_lev computes the integrated leverage of a diffusion process via the Fourier-Malliavin estimator  
%
%<a href="matlab: docsearchFS('FM_int_lev')">Link to the help function</a>
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
%              equal to floor(floor((length(x)-1)/2)^0.5).
%                 Example - 'M',20
%                 Data Types - single | double
%                 Data Types - single | double
%
% Output:
%
%  L_int  :    Integrated leverage. Scalar. Value of the integrated
%              leverage.
%
%
% More About:
%
% We assume that the vectors x contains discrete observations from a diffusion
% process $x$ following the Ito stochastic differential equation 
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$ 
% $$d\sigma^2(t)= \gamma(t) \ dZ(t) + a(t) \ dt,$$ 
% where  $W$ and $Z$ are two Brownian motions defined on the filtered probability space 
% $(\Omega, (\mathcal{F}_t)_{t \in [0,T]}, P)$, with correlation $\rho$, while 
% $\sigma, \gamma, b$ and $a$ are random processes, adapted to  $\mathcal{F}_t$.
% See the References for further  mathematical details.
% The integrated leverage on $[0,T]$ is defined as 
% $$IL_{[0,T]}:=  \langle  x,\sigma^2 \rangle_T =\rho\int_0^T\sigma(s)\gamma(s)ds.$$
% 
% 
%For any positive integer $n$, let $\mathcal{S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times. Moreover, let $\delta_l(x):=
% x(t_{l+1})-x(t_l)$ be the increments of $x$. 
% The Fourier estimator of the integrated leverage on $[0,T]$ is
% given by 
% $$T c_0(B_{{n},N})={T^2\over {2M+1}} \sum_{|k|\leq M} c_{k}(d\sigma_{n,N})c_{-k}(dx_{n}),$$
% 
% where:
% $$c_k(d\sigma_{n,N})= i k \frac{2\pi}{T} c_k(\sigma_{n,N}), \quad   c_k(\sigma_{{n},N})={T\over {2N+1}} \sum_{|s|\leq N} c_{s}(dx_{n})c_{k-s}(dx_{n}),$$
% $$c_k(dx_{n})= {1\over {T}} \sum_{l=0}^{n-1} e^{-{\rm i}\frac{2\pi}{T}kt_l}\delta_{l}(x).$$
% 
% See also: FM_int_vol.m, FM_int_cov.m, FM_int_quart.m, FM_int_volvol.m, Heston1D.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin Volatility Estimation. Theory and Practice, "Springer Briefs in Quantitative Finance", Springer. 
% 
% Sanfelici, S., Toscano, G. (2024), The Fourier-Malliavin Volatility (FMVol) MATLAB toolbox, available on ArXiv.
%
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FM_int_lev')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Example:  

%{
    %% Example of call of FM_int_lev with default values of N and M.
    % The following example estimates the integrated leverage
    % of a random process following the Heston model from a discrete sample. 
    % The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0; 0.5; 7; 2];
rho=-0.9;
x0=log(100); 
V0=0.5;
[x,V,t] = Heston1D(T,n,parameters,rho,x0,V0);
 

% Integrated leverage estimation
L_int = FM_int_lev(x,t,T);
 
%}

%{
    %% Example of call of FM_int_lev with custom choices of N and M.
    % The following example estimates the integrated leverage
    % of a random process following the Heston model from a discrete sample. 
    % The Heston model assumes that the spot variance follows 
    % a Cox-Ingersoll-Ross model.

% Heston model simulation
T=1;  
n=23400;  
parameters=[0; 0.5; 7; 2];
rho=-0.9;
x0=log(100); 
V0=0.5;
[x,V,t]=Heston1D(T,n,parameters,rho,x0,V0);
 

% Integrated leverage estimation
L_int= FM_int_lev(x,t,T,'N',10000,'M',70);
 
%}

%% Beginning of code

% Make sure that x and t are column vectors.

x=x(:);
t=t(:);

if length(x) ~= length(t)
    error('FSDA:FM_int_lev:WrongInputOpt','Input arguments x and t must have the same length.');
end


const=2*pi/T;

r=diff(x);  
  
        
n=length(r);  
N=floor(n/2);
M=floor(N^0.5); 

if nargin>2
    options=struct('N',N,'M',M);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FM_int_lev:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FM_int_lev:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
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
    error('FSDA:FM_int_lev:WrongInputOpt','N must be strictly smaller than min(n1,n2).');
end

if M  >= N
    error('FSDA:FM_int_lev:WrongInputOpt','M must be strictly smaller than N.');
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
 
center=N+M+1;

for j = 1 : 2*M+1        
    c_r_aux1(j,1:2*N+1)=   c_r1  (  center -N +( j-M-1 ) : center +N +( j-M-1 ) ) ;
    c_v(j)= T/(2*N+1) * c_r1  (  center -N   : center +N   ).' * flip(c_r_aux1(j, : )).';
    c_dv(j)=const*1i*(j-M-1 ) *(1-abs(j-M-1)/(M+1))*c_v(j);       
end



 
L_int= T^2/(M+1)*   real(  c_dv.'*flip(c_r1(center-M: center+M)));

end

%FScategory:FMvol