function tpdf = tobitpdf(x,   mu,sigma, left, right)
%tobitpdf returns probability density function from the tobit model
%
%<a href="matlab: docsearchFS('tobitpdf')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    x:         Value at which the pdf must be evaluated.
%               Scalar, vector or matrix 3D array of the same size of mu,
%               sigma, Lower and Upper
%               A scalar input functions as a constant matrix of the
%               same size as the other input.
%               See "More About:" for details about the tobit
%               distribution.
%               Data Types - single | double
%
%  Optional input arguments:
%
%    mu :       location parameter of the tobit distribution.
%               Scalar, vector or matrix 3D array of the same size of x and sigma, Lower, Upper.
%               A scalar input functions as a constant matrix of the same
%               size as the other input. Default value of mu is 0.
%               See "More About:" for details about the tobit
%               distribution.
%               Example - 'mu',1
%               Data Types - single | double
%
%    sigma :    scale parameter of the tobit distribution.
%               Scalar, vector or matrix 3D array of the same size of x and sigma, Lower, Upper.
%               A scalar input functions as a constant matrix of the same
%               size as the other input. Default value of sigma is 1
%               See "More About:" for details about the tobit
%               distribution.
%               Example - 'sigma',10
%               Data Types - single | double
%
%    left :     lower limit for the censored random variable. Scalar.
%               If set to -Inf, the random variable is assumed to be not
%               left-censored; default value of left is zero (classical
%               tobit model).
%               Example - 'left',1
%               Data Types - double
%
%    right :    right limit for the censored random variable. Scalar.
%               If set to Inf, the random variable is assumed to be not
%               right-censored; default value of left is Inf (classical
%               tobit model).
%               Example - 'right',800
%               Data Types - double
%
%  Output:
%
%           tpdf : tobit pdf values. 
%                  The size of tpdf is the common size of the input
%                  arguments. A scalar input functions as a constant matrix
%                  of the same size as the other inputs.
%
%
% See also: tobitcdf, tobitrnd, tobitinv
%
% References:
%
% Greene, W.H. (2008), "Econometric Analysis, Sixth Edition", Prentice Hall, pp. 871-875.
%
% Tobin, J. (1958), Estimation of Relationships for Limited Dependent
% Variables, "Econometrica", 26, pp. 24-36.
%
%
% Copyright 2008-2025.
%
%<a href="matlab: docsearchFS('tobitpdf')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%

%{
    % Example when x is a scalar.
    x=0;
    % The default values for mu and sigma are 0, 1
    % The default values for left and right are 0 Inf
    y=tobitpdf(x); 
    ynorm=normcdf(x);
    assert(y==ynorm,'Tobit density is not correct')
    % Pr(Tobit=0, mu=0, sigma=1, left=0, right=Inf) = Pr (Z<0)
%}


%{
    % In this example mu, and sigma are specified.
    x=0; mu=2; sigma=1.5;
    % The default values for left and right are 0 Inf
    y=tobitpdf(x,mu,sigma); 
    ynorm=normcdf(x,mu,sigma);
    assert(y==ynorm,'Tobit cdf is not correct')
    % Pr(Tobit=0, mu=2, sigma=1.5, left=0, right=Inf) = Pr (Z<0)
%}

%{
    % In this example mu, sigma, left and right are specified.
    x=0; mu=2; sigma=1.5; left=1; right=3;
    y=tobitpdf(x,mu,sigma,left,right); 
    ynorm=normcdf(x,mu,sigma);
    % Pr(Tobit=0, mu=2, sigma=1.5, left=0, right=Inf) = Pr (Z<0)
%}

%{
    % Example where x is not a scalar.
    x=0:10; mu=2; sigma=1.5; left=1; right=3;
    y=tobitcdf(x,mu,sigma,left,right); 
    ynorm=normpdf(x,mu,sigma);
    
%}


%{
    %% Example of tobit density.
    close all
    x=(-3:0.0001:3)';
    left=0;
    right=2;
    mu=0.5;
    sigma=1;
    x(find(x<left,1,'last'))=NaN;
    x(find(x>left,1,'first'))=NaN;
    x(find(x<right,1,'last'))=NaN;
    x(find(x>right,1,'first'))=NaN;
    y=tobitpdf(x,mu,sigma,left,right);
    plot(x,y,'LineWidth',2)
    hold('on')
    stem(left,y(x==left),'Color','b')
    stem(right,y(x==right),'Color','b')
    title(['Tobit density when \mu=' num2str(mu) ', \sigma=' num2str(sigma) ', ' ...
        'left=' num2str(left)  ', right='  num2str(right)])
    text(left,tobitpdf(left,mu,sigma,left,right)-0.01, ...
        'Pr(Tobit(\mu,\sigma^2,left,right)=left)=\Phi(left,\mu, \sigma^2) ','HorizontalAlignment','right')
    text(right,tobitpdf(right,mu,sigma,left,right)+0.01, ...
        ' Pr(Tobit(\mu,\sigma^2,left,right)=right)=1-\Phi(right, \mu, \sigma^2)','HorizontalAlignment','left')
    ylim([0 1])
%}


%% Beginning of code

if nargin<5
    right =  Inf;
end

if nargin<4
    left =  0;
end

if nargin<3
    sigma =  1;
end

if nargin<2
    mu =  0;
end

if nargin < 1
    error(message('FSDA:tobitpdf:TooFewInputs'));
end

[errorcode, x, mu, sigma, left, right] = distchck(5,x,mu,sigma,left,right);

if errorcode > 0
    error(message('FSDA:tobitpdf:InputSizeMismatch'));
end

% Initialize Y to zero.
outType = internal.stats.dominantType(x,mu,sigma,left,right);
tpdf = zeros(size(x),"like",outType);

for i=1:length(x)
    if x(i)>left(i) && x(i)<right(i)
        tpdf(i) = normpdf(x(i),mu(i),sigma(i));
    elseif x(i)==left(i)
        tpdf(i) = normcdf(x(i),mu(i),sigma(i));
    elseif x(i)==right(i)
        tpdf(i) = normcdf(x(i),mu(i),sigma(i),'upper');
    else
    end  
end

end

%FScategory:ProbDist
