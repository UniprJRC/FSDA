function y = inversegamcdf(x,a,b,nocheck)
%inversegamcdf computes inverse-gamma cumulative distribution function.
%
%<a href="matlab: docsearchFS('inversegamcdf')">Link to the help function</a>
%
%  Required input arguments:
%
%    x:         Value at which the cdf must be evaluated.
%               Scalar, vector or matrix 3D array of the same size of x and b.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%               Data Types - single | double
%    a :        shape parameter of the inverse-gamma distribution.
%               Scalar, vector or matrix 3D array of the same size of x and b.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%               Data Types - single | double
%    b :        scale parameter b of the inverse-gamma distribution.
%               Scalar, vector or matrix 3D array of the same size of x and a.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               Unlike the Gamma distribution, which contains a somewhat
%               similar exponential term, $b$ is a scale parameter as the
%               distribution function satisfies:
%               \[
%                   f_{IG}(x,a,b)=\frac{f(x/b,a,1)}{b}
%               \]
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%               Data Types - single | double
%
%  Optional input arguments:
%
%      nocheck: Check input arguments. Scalar. If nocheck is equal to 1 no
%               check is performed and input and the cdf is evaluated
%               directly through MATLAB buit in function gammainc
%               else we use MATLAB function gamcdf
%               Example - 'nocheck',1
%               Data Types - double
%
%  Output:
%
%    y:         CDF value. Scalar, vector or matrix or 3D array of the same size
%               of input arguments x, a and b. $y=\int_0^x f_{IG}(t | a,b) dt$ is the
%               value of the cdf of the inverse gamma distribution
%               evaluated at x
%
% More About:
%
%
% The cdf of the inverse gamma distribution defined over the support
% $x>0$ with shape parameter $a$ and scale parameter $b$ is
%  \[
%  F_{IG}(x, a, b)  =\int_0^x t^{-a -1} \exp (-b/t)
%  \frac{b^a}{\Gamma(a)} dt
%  \]
%
%
% inversegamcdf computes the inverse-gamma cdf at each of the values in x using the
% corresponding shape parameters in a and scale parameters in b. x, a, and
% b can be vectors, matrices, or multidimensional arrays that all have the
% same size. A scalar input is expanded to a constant array with the same
% dimensions as the other inputs. The parameters in a and b must all be
% positive, and the values in x must lie on the interval $[0,  \infty)$.
%
% Note that $F_{IG}(x,a,b)=\frac{\Gamma(a,b/x)}{\Gamma(\alpha)}$ therefore  
% Therefore, the CDF for an inverse Gamma distribution can be computed
% using the incomplete gamma function (also called regularized gamma
% function, i.e. MATLAB function gammainc) of course  keeping into account
% that we need the upper tail.
%
% The chief use of the inverse gamma distribution is in Bayesian
% statistics, where the distribution arises as the marginal posterior
% distribution for the unknown variance of a normal distribution if an
% uninformative prior is used; and as an analytically tractable conjugate
% prior if an informative prior is required.
% Relation with the Gamma distribution.
% If $X \sim Gamma(a,b)$ then $\frac{1}{X} \sim$ inverse-gamma distribution
% with paramters $a$ and $1/b$.
%
% See also: gampdf,inversegampdf
%
% References:
%
% https://en.wikipedia.org/wiki/Inverse-gamma_distribution
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('inversegamcdf')">Link to the help function</a>
% Last modified 11-06-2016
%

% Examples:
%
%{
    %% Plot the cdf for 4 different combinations of parameter values.
    x=(0:0.001:3)';
    a=[1,2,3,3];
    b=[1,1,1,0.5];
    Y=zeros(length(x),length(a));
    for i=1:length(x)
        Y(i,:)=inversegamcdf(x(i),a,b);
    end

    for j=1:4
        subplot(2,2,j)
        plot(x,Y(:,j))
        title(['CDF with a=' num2str(a(j)) ' b=' num2str(b(j))])
        xlabel('x')
    end
%}

%{
    %% Compare the results using option nocheck=1.
    x=(0:0.001:3)';
    a=[1,2,3,50,100,10000];
    b=[1,10,100,0.05,10,800];
    Y=zeros(length(x),length(a));
    Ychk=Y;
    for i=1:length(x)
        Y(i,:)=inversegamcdf(x(i),a,b);
        Ychk(i,:)=inversegamcdf(x(i),a,b,1);
    end
    disp('Maximum absolute difference is:');
    disp(max(max(abs(Y-Ychk))));
%}


%% Beginning of code

if nargin<4
    nocheck=0;
end

if nocheck==1
    y = gammainc(b./x,a,'upper');
else
    y = gamcdf(1./x,a,1./b,'upper');
end

% inversegamcdf(10,2,3)
% pinvgamma(10, 2, scale = 3, log = FALSE)
end
%FScategory:UTISTAT