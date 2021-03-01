function [spotvar,tau] = FE_spot_vol_FFT(x,t, varargin)
% FE_spot_vol estimates the instantaneous variance of a diffusion process by the Fourier estimator with Dirichlet kernel, using the FFT algorithm
%
%
%<a href="matlab: docsearchFS('FE_spot_vol_FFT')">Link to the help function</a>
%
%
% Required input arguments:
%
%   x   :    Observation values. Vector. Row or column vector containing
%            the observed values.
%   t   :    Observation times. Vector.  Row or column vector with the same
%            length of x containing the observation times
%
% Optional input arguments:
%
%   N    :   cutting frequency. Scalar. If N is not specified, it is set
%            equal to (length(x)-1)/2.
%                 Example - 'N',500
%                 Data Types - single | double
%
%   M    :   maximum Fourier frequency for spot variance reconstruction.
%            Scalar. If M is not specified, it is set equal to
%            sqrt(length(x))*log(length(x)).
%                 Example - 'M',218
%                 Data Types - single | double
%
% Output:
%
%   spotvar : Spot variance. Vector. Row vector with length 2M+1 of spot variance values at the time grid tau.
%
%   tau :    Estimation times. Vector.  Row vector with length 2M+1
%            where the spot variance is estimated.
%
%
% More About:
%
% We assume our timeseries data are discrete observations from a diffusion
% process $x$ following the It\^o stochastic differential equation 
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$ 
% where $W$ is a Brownian motion on a filtered probability space. Let
% $\sigma$ and $b$ be random processes, adapted to the Brownian filtration.
% See the Reference for further  mathematical details.
% The squared diffusion coefficient $\sigma^2(t)$ is called
% instantaneous variance of the process $x$.
%
% For any positive integer $n$, let ${\cal S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times. Moreover, let $\delta_i(x):=
% x(t_{i+1})-x(t_i)$ be the increments of $x$.
% 
% The Fourier estimator of the instantaneous variance for $t \in (0,T)$,
% is defined as 
% $$\widehat \sigma^2_{n,N,M}(t):= \sum_{|k|\leq M} \left(1- {|k|\over
% M}\right) e^{{\rm i} {{2\pi}\over {T}} tk} c_k(\sigma^2_{n,N}),$$
% where
% $$c_k(\sigma^2_{n,N}):= {T \over {2N+1}}\sum_{|s|\leq N} c_s(dx_n)
% c_{k-s}(dx_n),$$
% where for any integer $s$, $|s|\leq 2N$, the discretized Fourier
% coefficients of the increments are
% $$c_s(dx_{n}):= {1\over {T}} \sum_{i=0}^{n-1} e^{-{\rm i} {{2\pi}\over {T}}
% st_i}\delta_i(x).$$
%
%
% See also: FE_spot_vol.m, FE_int_vol.m, FE_int_vol_Fejer.m, CEVmodel.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin
% Volatility Estimation. Theory and Practice, "Springer Briefs in
% Quantitative Finance", Springer.
%
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FE_spot_vol_FFT')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of call of FE_spot_vol with just two input arguments.
    % Generates price and instantaneous variance from the Constant Elasticity of Variance model [S. Beckers, The Journal of Finance, Vol. 35, No. 3, 1980] and estimates the instantaneous variance via the Fourier method.
    n=1000;  T=1; t=0:T/n:T; S0=100;
    [S,sigma]=CEVmodel(t,S0); % data generation
    x=log(S); % log-price
    spotvar=FE_spot_vol_FFT(x,t);
    plot(spotvar(2:end-1))
    ylabel('Spot variance')
%}

%{
    %% FE_spot_vol called with optional input argument N and M.
    % Generates price and instantaneous variance from the Constant Elasticity of Variance model [S. Beckers, The Journal of Finance, Vol. 35, No. 3, 1980] and estimates the instantaneous variance via the Fourier method.
    n=21600;  T=1; t=0:T/n:T; S0=100;
    [S,sigma]=CEVmodel(t,S0); % data generation
    x=log(S); % log-price
    N=floor(n/2); M=floor(n^0.5); % cutting frequencies
    spotvar=FE_spot_vol_FFT(x,t,'N',N,'M',M);
    plot(spotvar(2:end-1))
    ylabel('Spot variance')
%}

%{
    %% Example of call of FE_spot_vol providing both estimation values and times.
    % Generates price and instantaneous variance from the Constant Elasticity of Variance model [S. Beckers, The Journal of Finance, Vol. 35, No. 3, 1980] and estimates the instantaneous variance via the Fourier method.
    n=21600;  T=1; t=0:T/n:T; S0=100;
    [S,sigma]=CEVmodel(t,S0); % data generation
    x=log(S); % log-price
    N=floor(n/2); M=150; % cutting frequencies
    [spotvar,tau]=FE_spot_vol_FFT(x,t,'N',N,'M',M);
    plot(tau,spotvar)
    xlabel('time')
    hold on; plot(tau,sigma(1:n/(2*M):end));
    ylabel('Spot variance')
    title('Estimated value (blue), True value (red)')
%}

%% Beginning of code

% Make sure that x and t are both column vector.
x=x(:);
t=t(:);
N=floor((length(x)-1)/2); M=floor(sqrt(length(x))*log(length(x)));

if length(x) ~= length(t)
    error('FSDA:FE_spot_vol:WrongInputOpt','Input arguments x and t must have the same length.');
end

if nargin>2
    options=struct('N',N,'M',M);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FE_spot_vol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FE_spot_vol:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    N=options.N; M=options.M;
end

r=diff(x)'; % must be a row vector
T=t(end)-t(1); tau=t(1):T/(2*M):t(end);
fft_v=fft(r); 
idx=M+N+1:-1:2; ff=fft_v(idx);
fft_def=[conj(ff) fft_v(1:M+N+1)];
fft_def=fft_def./T; % Fourier coeff. of log-returns

idxk=-N:1:N; nshift=M+N+1;
for kk=-M:M
    idxx=idxk+nshift+kk;
    Capp=fft_def(idxx);
    coeff(M+kk+1)=Capp*fft_def(nshift-N:nshift+N)';
end
C=coeff.*(T/(2*N+1)); % Fourier coeff. of variance

k=(-M:1:M);
f=C.*(1-abs(k)/M);
Fsum=(2*M+1)*ifft(f);
Fsum=exp(-1i*2*pi*M*(k+M)/(2*M+1)).*Fsum;
spotvar=real(Fsum);

end
%FScategory:UTISTAT