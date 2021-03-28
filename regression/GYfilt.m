function weights = GYfilt(x,varargin)
%GYfilt computes the Gervini-Yohai univariate outlier identifier
%
%<a href="matlab: docsearchFS('GYfilt')">Link to the help function</a>
%
%  Required input arguments:
%
%    x:         Input vector. Vector. A vector with n elements that
%               contains the univariate data.
%
%
%  Optional input arguments:
%
%      alpha :  coverage probability. Scalar.
%               Scalar in the interval [0.5 1). The default coverage
%               probability is 0.95.
%                 Example - 'alpha',0.99
%                 Data Types - double
%
%   centering:  centering the data. Boolean.
%               If centering is true input data are preliminarly centered.
%               The defalt value of centering is true.
%                 Example - 'centering',false
%                 Data Types - logical
%
%   iterating:  iterative procedure. Boolean.
%               If Boolean is true then an iterative adaptive procedure is
%               applied.  The defalt value of iterating is true.
%                 Example - 'iterating',false
%                 Data Types - logical
%
%      niter :  maximum number of iterations in the iterative adaptive
%               procedure. Positive integer. This option is used just if previous iterating
%               is true. The default value of niter is 10.
%                 Example - 'niter',20
%                 Data Types - double
%
%  Output:
%
%    weights:   Boolean vector of weights. Logical.
%               A boolean vector with n elements that contains false in
%               correspondence of the units declared as outliers.
%
%
% See also: LTSts
%
% References:
%
% Gervini, D. and Yohai, V.J. (2002), A class of robust and fully efficient
% regression estimators, "Annals of Statistics", Vol. 30, pp. 583-616.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('GYfilt')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % GYfilt with all the default options.
    weights=GYfilt(randn(100,1));
%}

%{
    % GYfilt with option alpha.
    alpha=0.999;
    weights=GYfilt(randn(100,1),'alpha',alpha);
%}

%% Beginning of code

if ~isvector(x)
    error('FSDA:GYfilt:WrongInputOpt','The data should be a vector')
end

alphadef=0.95; % default coverage probability
centering=true;
iterating=true;
niterdef=10;

if coder.target('MATLAB')
    options=struct('alpha',alphadef,'centering',centering,'iterating',iterating,'niter',niterdef);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:GYfilt:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

alpha=options.alpha;
centering=options.centering;
iterating=options.iterating;
niter=options.niter;

if centering==true
    mu=median(x);
else
    mu=0;
end
madx=1.4826*mad(x,1);

xs= (x-mu)/madx;
xs2 = xs.^2;
n=length(x);

if iterating
    xs2na = gyiterate(xs2,alpha,niter);
else
    xs2na = gyfiltaux(xs2,alpha);
end

weights =true(n,1);
weights( isnan(xs2na)) =0;


end

% Inner fucntion gyfiltaux
function vna=gyfiltaux(v, alpha)

[v,vorder] = sort(v);
n=length(v);
seq=1:length(v);
i0 = seq(v < chi2inv(alpha, 1));
n0 = 0;
if ~isempty(i0)
    i0 = i0(end);
    t1=chi2cdf(v(i0:n),1)-(((i0:n) - 1)')./n;
    boo=t1>0;
    dnt=zeros(length(t1),1);
    dnt(boo)=t1(boo);
    
    dn = max(dnt);
    n0 = round(dn * n);
end
[~,order]=sort(vorder);
v =v(order);
vna=v;
if (n0 > 0)
    vna(vorder((n - n0 + 1):n)) =NaN;
end
end

% Inner function gyiterate
function vout=gyiterate(v, alpha, niter)
if nargin<3
    niter=10;
end
converge =0;
iter =0;
n=length(v);
id =1:n;
% vold = v;
while converge == 0 &&  iter < niter
    iter = iter + 1;
    v = gyfiltaux(v, alpha);
    id = id(~isnan(v));
    if (~any(isnan(v)))
        converge = 1;
    end
    v =v(~isnan(v));
end
vout =nan(n,1);
vout(id) = v;
% disp([' iter = ' num2str(iter)])
fprintf(' iter = %.0f\n ',iter)
end

%FScategory:UTISTAT