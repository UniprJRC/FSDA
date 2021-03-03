function ivar = FE_int_vol(x,t, varargin)
%FE_int_vol computes the integrated variance from a diffusion process via the Fourier estimator using Dirichlet kernel
%
%
%<a href="matlab: docsearchFS('FE_int_vol')">Link to the help function</a>
%
% FE_int_vol computes the integrated variance of univariate timeseries data
% from a diffusion process by the Fourier estimator with Dirichlet kernel
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
%           equal to (length(x)-1)/2.
%                 Example - 'N',500
%                 Data Types - single | double
%
% Output:
%
%   ivar    integrated variance. Scalar. Value of the integrated variance.
%
%
% More About:
%
% We assume our timeseries data are discrete observations from a diffusion
% process $x$ following the Ito stochastic differential equation 
% $$dx(t)= \sigma(t) \ dW(t) + b(t) \ dt,$$ 
% where $W$ is a Brownian motion on a filtered probability space. Let
% $\sigma$ and $b$ be random processes, adapted to the Brownian filtration.
% See the Reference for further  mathematical details.
% The integrated variance of the process over the time interval $[0,T]$ is defined as
% $$\int_0^T \sigma^2(t) dt.$$
%
% For any positive integer $n$, let ${\cal S}_{n}:=\{ 0=t_{0}\leq \cdots
% \leq t_{n}=T  \}$ be the observation times. Moreover, let $\delta_i(x):=
% x(t_{i+1})-x(t_i)$ be the increments of $x$.
% 
% The Fourier estimator of the integrated variance over $[0,T]$, is
% defined as 
% $$\widehat\sigma^{2}_{n,N}:= {T^2 \over {2N+1}}\sum_{|s|\leq N} c_s(dx_n)
% c_{-s}(dx_n),$$
% where for any integer $k$, $|k|\leq N$, the discretized Fourier
% coefficients of the increments are
% $$c_k(dx_{n}):= {1\over {T}} \sum_{i=0}^{n-1} e^{-{\rm i} {{2\pi}\over {T}}
% kt_i}\delta_i(x).$$
% The cutting frequency $N$ is a scalar integer. If not specified, $N$ is
% set equal to $n/2$ (Nyquist frequency).
%
% See also: FE_int_vol_Fejer.m, OptimalCuttingFrequency.m, FE_spot_vol.m
%
% References:
%
% Mancino, M.E., Recchioni, M.C., Sanfelici, S. (2017), Fourier-Malliavin
% Volatility Estimation. Theory and Practice, "Springer Briefs in
% Quantitative Finance", Springer.
%
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('FE_int_vol')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of call of FE_int_vol with just two input arguments.
    % The following example calculates the integrated variance from a
    % vector x of discrete observations of a univariate diffusion process
    n=1000;
    dt=1/n; t=0:dt:1;
    x=randn(n,1)*sqrt(dt); x=[0;cumsum(x)];
    ivar=FE_int_vol(x,t);
    disp(['The value of the integrated variance is: ' num2str(ivar)])
%}

%{
    %% FE_int_vol called with optional input argument N. 
    % Analysis of the change of the integrated variance estimates of a
    % univariate diffusion process as a function of the different values of
    % the cutting frequency.
    n=1000;
    dt=1/n; t=0:dt:1;
    x=randn(n,1)*sqrt(dt); x=[0;cumsum(x)];
    cuttingfreq=(50:500)';
    l=length(cuttingfreq);
    Ivar=zeros(l,1);
    for i=1:l
        ivar=FE_int_vol(x,t,'N',cuttingfreq(i));
        Ivar(i)=ivar;
    end
    plot(cuttingfreq,Ivar)
    xlabel('Cutting frequency')
    ylabel('Integrated variance')
%}



%% Beginning of code

% Make sure that x and t are both column vector.
x=x(:);
t=t(:);
N=floor((length(x)-1)/2);

if length(x) ~= length(t)
    error('FSDA:FE_int_vol:WrongInputOpt','Input arguments x and t must have the same length.');
end

if nargin>2
    options=struct('N',N);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FE_int_vol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:FE_int_vol:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    N=options.N;
end

r=diff(x); T=t(end)-t(1); const=2*pi/T;
c_p=zeros(N,1);
c_0=sum(r);
for k=1:N
    c_p(k)=sum(exp(1i*const*k*t(1:end-1)).*r);
end
ivar=(c_0.*conj(c_0)+2*sum(c_p.*conj(c_p)))/(2*N+1);

end
%FScategory:UTISTAT
