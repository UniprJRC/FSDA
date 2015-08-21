function y = inversegampdf(x,a,b,nocheck)
%inversegampdf Inverse-gamma probability density function.
%
%<a href="matlab: docsearchFS('inversegampdf')">Link to the help function</a>
%
%  Required input arguments:
%
%    x:         Value at which the pdf must be evaluated.
%               Scalar, vector or matrix 3D array of the same size of x and b.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%    a :        shape parameter of the inverse-gamma distribution. 
%               Scalar, vector or matrix 3D array of the same size of x and b.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%    b :        scale parameter b of the inverse-gamma distribution. 
%               Scalar, vector or matrix 3D array of the same size of x and a.
%               A scalar input functions as a constant matrix of the same
%               size as the other input.
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%
%  Optional input arguments:
%
%      nocheck: Check input arguments. Scalar. If nocheck is equal to 1 no
%               check is performed and input and the density is evaluated
%               directly through the expression
%               y = (b.^a).*(x.^(-a-1)).*exp(-b./x)./gamma(a); else 
%               we use MATLAB function gampdf
%               Example - 'nocheck',1 
%               Data Types - double
%
%
%  Output:
%
%    y:         pdf value. Scalar, vector or matrix or 3D array of the same size
%               of input arguments x, a and b. $y=f_{IG}(x | a,b)$ is the
%               value of the pdf of the inverse gamma distribution
%               evaluated at x
%
% More About:
%
%
% The density of the inverse gamma distribution defined over the support
% $x>0$ with shape parameter $a$ and scale parameter $b$ is
%  \[
%  f_{IG}(x, a, b)  \propto x^{-a -1} \exp (-b/x)
%  \frac{b^a}{\Gamma(a)}
%  \]
%
% inversegampdf computes the gamma pdf at each of the values in x using the corresponding
% shape parameters in a and scale parameters in b. x, a, and b can be
% vectors, matrices, or multidimensional arrays that all have the same
% size. A scalar input is expanded to a constant array with the same
% dimensions as the other inputs. The parameters in a and b must all be
% positive, and the values in x must lie on the interval $[0,  \infty)$.
% The chief use of the inverse gamma distribution is in Bayesian
% statistics, where the distribution arises as the marginal posterior
% distribution for the unknown variance of a normal distribution if an
% uninformative prior is used; and as an analytically tractable conjugate
% prior if an informative prior is required.
% Relation with the Gamma distribution.
% If $X \sim Gamma(a,b)$ then $\frac{1}{X} \sim$ inverse-gamma distribution
% with paramters $a$ and $1/b$
% 
% See also: gampdf
%
% References:
% 
% https://en.wikipedia.org/wiki/Inverse-gamma_distribution
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('inversegampdf')">Link to the help function</a>
% Last modified 06-Feb-2015
%

% Examples:
%
%{
    %% Plot the pdf for 4 different combinations of parameter values.
    x=(0:0.001:3)';
    a=[1,2,3,3];
    b=[1,1,1,0.5];
    Y=zeros(length(x),length(a));
    for i=1:length(x)
        Y(i,:)=inversegampdf(x(i),a,b);
    end

    for j=1:4
        subplot(2,2,j)
        plot(x,Y(:,j))
        title(['PDF with a=' num2str(a(j)) ' b=' num2str(b(j))])
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
        Y(i,:)=inversegampdf(x(i),a,b);
        Ychk(i,:)=inversegampdf(x(i),a,b,1);
    end
    disp('Maximum absolute difference is:')
    disp(max(max(abs(Y-Ychk))))
%}


%% Beginning of code

if nargin<4
    nocheck=0;
end

if nocheck==1
    % This expression is much faster but (less accurate in the extreme cases) and without
    % checks on the input. Our experience (10 milions of simulations) is
    % that the difference is never greater than 1e-9
    % y = (b.^a).*(x.^(-a-1)).*exp(-b./x)./gamma(a); 
    y=exp(a.*log(b)-(a+1).*log(x)-b./x - gammaln(a));
else
% Standard of finding the pdf of the inverse Gamma distribution using the gampdf function
%  and the Jacobian (1/x^2)
    y = gampdf(1./x,a,1./b)./(x.^2);
end
 
% Return NaN for out of range parameters or missing values.
y(a < 0) = NaN;
y(b <= 0) = NaN;

% remark
% dinvgamma(10, 2, scale = 3, log = FALSE)
% inversegampdf(10,2,3)
end