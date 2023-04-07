function [res] = mdpd(y, alpha, varargin)
%mdpd computes Minimum Distance Power Divergence statistics
%
%<a href="matlab: docsearchFS('mdpd')">Link to the help function</a>
%
%
% The Power Divergence for a density function $f$ and observations
% $y_1 , \ldots , y_n$ is defined as
%
% \[
% \Delta(f,\alpha) = \int_R f^{1+\alpha}(y) dy - (1+1/\alpha) \sum_{i=1}^n f^\alpha(y_i)/n
% \]
%
% for $\alpha > 0$
%
% \[
% \Delta(f,0) = -\sum_{i=1}^n \log f(y_i)/n
% \]
%
% for $\alpha = 0$.
%
%
%  Required input arguments:
%
%         y:    Response variable. Vector. A vector with n elements that
%               contains the response variable.  It can be either a row or
%               a column vector.
%                 Data Types - double
%
%    alpha :    Numeric for the power divergence parameter. Non negative
%               scalar. It can be shown that as the tuning parameter
%               $\alpha$ increases the robustness of the Minimum Density
%               Power Divergence estimator increases while its efficiency
%               decreases (Basu et al., 1998).
%                 Data Types - double
%
% Optional input arguments:
%
% densfunc  :   handle to the function computing the theoretical density.
%               Function handle. Function handle which defines the function
%               to be integrated from lower to upper. The default density
%               function is the standard normal distribution.
%                 Example - 'densfunc', @tpdf
%                 Data Types - handle
%
% lower :       Lower bound of the domain of the density function.
%               Scalar. The default value of lower is 1.
%                 Example - 'lower',0
%                 Data Types - double
%
% upper :       Upper bound of the domain of the density function.
%               Scalar. The default value of upper is Inf.
%                 Example - 'upper',10
%                 Data Types - double
%
%  theta :      The parameters of the distribution given as a vector.
%               Numeric vector. The default values of theta is [0 1] given
%               that the default density function is the standard normal
%               distribution.
%                 Example - 'theta', [10 100]
%                 Data Types - double
%
%  RelTol :     Relative error tolerance. Non negative scalar.
%               Relative error tolerance, specified as the comma-separated
%               pair consisting of 'RelTol' and a nonnegative real number.
%               mpdm uses the relative error tolerance to limit an estimate
%               of the relative error, |q - Q| / min(|q|,|Q|), where q is the computed
%               value of the integral and Q is the (unknown) exact value.
%               The default value of RelTol is 1e-6.
%               Example - 'RelTol', 1e-12
%                 Data Types - double
%
%  AbsTol :     Absolute error tolernace. Non negative scalar.
%               Absolute error tolerance, specified as the comma-separated
%               pair consisting of 'AbsTol' and a nonnegative real number.
%               mpdm uses the absolute error tolerance to limit an estimate
%               of the absolute error, |q - Q|, where q is the computed
%               value of the integral and Q is the (unknown) exact value.
%               This can be useful when q or Q becomes close to zero and
%               the relative tolerance risks to go to infinity.
%               The default value of AbsTol is 1e-12.
%               Example - 'AbsTol', 1e-9
%                 Data Types - double
%
%
%  Output:
%
%   res :       power divergence against the density function densfunc.
%               Scalar. Value of power divergence wrt the density specified
%               in densfunc.
%                  Data Types - double.
%
%
% See also: normpdf, mdpdR, mdpdReda, PDrho 
%
%  References:
%
%  Basu, A., Harris, I.R., Hjort, N.L. and Jones, M.C., (1998), Robust
%  and efficient estimation by minimizing a density power divergence,
%  Biometrika, 85, pp. 549-559.
%
%  Riani, M. Atkinson, A.C., Corbellini A. and Perrotta A. (2020), Robust
%  Regression with Density Power Divergence: Theory, Comparisons and Data
%  Analysis, Entropy, Vol. 22, 399. 
%  https://www.mdpi.com/1099-4300/22/4/399 
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdpd')">Link to the help function</a>
%
%$LastChangedDate:: 2019-05-14 16:04:25 #$: Date of the last commit


% Examples:

%{
    %% mdpd with all default arguments.
    rng('default')
    y =randn(10,1);
    out=mdpd(y,0.2);
%}


%{
    %% Student T with 5 degrees of freedom.
    rng('default')
    y =randn(10,1);
    mdpd(y,0.1,'densfunc',@tpdf,'theta',5)
%}


%{
    % Change lower and upper integration limits.
    rng('default')
    % Generate 10 numbers from Uniform.
    y =mtR(10,0);
    out=mdpd(y,0.1,'densfunc',@tpdf,'theta',5, 'lower',-Inf,'upper',Inf);
    expectedRes = -8.870156744130417;
    assert(isequal(round(out,7), round(expectedRes,7)), 'Error: MATLAB did not output the expected result!')
%}

%{
    % Change lower and upper integration limits and relatve tolerance.
    % Change lower and upper integration limits.
    rng('default')
    % Generate 10 numbers from Uniform.
    y =mtR(10,0);
    out=mdpd(y,0.1,'densfunc',@tpdf,'theta',5, 'lower',-Inf,'upper',Inf,'RelTol',1e-15);
    expectedRes =  -8.870156744130275;
    assert(isequal(round(out,7), round(expectedRes,7)), 'Error: MATLAB did not output the expected result!')

%}

%% Beginning of code
if ~isscalar(alpha)
    error('FSDA:mdpd:WrongInputOpt','alpha should be a non negative scalar')
else
    if alpha<0
        error('FSDA:mdpd:WrongInputOpt','alpha should be a non negative scalar')
    end
end

if nargin<2
    error('FSDA:mdpd:missingInputs','y or alpha missing')
end

% Absolute tolerance of the integral
AbsTol = 1e-12;
% Relative tolerance of the integral
RelTol = 1e-6;

% Lower and upper limit
lower = 1;
upper = Inf;

% Default density is normal with parameters 0, 1 (i.e. standardized normal)
densfunc=@normpdf;
% Default parameter of the (normal) distribution given as a vector.
theta=[0 1];

if nargin > 2
    options=struct('densfunc',densfunc,'theta',theta, ...
        'lower',lower,'upper',upper','RelTol',RelTol,'AbsTol',AbsTol);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdpd:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    RelTol = options.RelTol;
    AbsTol = options.AbsTol;
    lower  = options.lower;
    upper  = options.upper;
    densfunc=options.densfunc;
    theta  = options.theta;
end

n      = length(y);
ltheta = length(theta);

if (alpha == 0)
    
    switch ltheta
        case 1
            res = -sum(log(densfunc(y, theta)))/n;
        case 2
            res = -sum(log(densfunc(y, theta(1), theta(2))))/n;
        case 3
            res = -sum(log(densfunc(y, theta(1), theta(2), theta(3))))/n;
        otherwise
            res = -sum(log(densfunc(y, theta(1), theta(2), theta(3), theta(4))))/n;
    end
    
else
    
    switch ltheta
        case 1
            fun = @(x) densfunc(x, theta).^(alpha+1);
            res = -(1 + 1/alpha) * sum(densfunc(y, theta(1)).^(alpha))/n;
        case 2
            fun = @(x) densfunc(x, theta(1), theta(2)).^(alpha+1);
            res = -(1 + 1/alpha) * sum(densfunc(y, theta(1), theta(2)).^(alpha))/n;
        case 3
            fun = @(x) densfunc(x, theta(1), theta(2), theta(3)).^(alpha+1);
            res = -(1 + 1/alpha) * sum(densfunc(y, theta(1), theta(2), theta(3)).^(alpha))/n;
        otherwise
            fun = @(x) densfunc(x, theta(1), theta(2), theta(3), theta(4)).^(alpha+1);
            res = -(1 + 1/alpha) * sum(densfunc(y, theta(1), theta(2), theta(3), theta(4)).^(alpha))/n;
    end
    
    I = integral(fun, lower,upper,'RelTol',RelTol,'AbsTol',AbsTol);
    res = res + I;
end

end
%FScategory:UTISTAT