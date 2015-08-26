function x = inversegaminv(p,a,b,nocheck)
%inversegampdf Inverse-gamma cumulative distribution function.
%
%<a href="matlab: docsearchFS('inversegaminv')">Link to the help function</a>
%
%  Required input arguments:
%
%    p:         Probability at which the inverse of the cdf must be
%               evaluated $0 \leq p \leq 1$.
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
%               Unlike the Gamma distribution, which contains a somewhat
%               similar exponential term, $b$ is a scale parameter as the
%               distribution function satisfies:
%               \[
%                   f_{IG}(x,a,b)=\frac{f(x/b,a,1)}{b}
%               \]
%               See the section called "More About:" for more details about
%               the inverse gamma distribution.
%
%  Optional input arguments:
%
%      nocheck: Check input arguments. Scalar. If nocheck is equal to 1 no
%               check is performed and input and the inverse cdf is evaluated
%               directly through MATLAB buit in function gammaincinv
%               else we use MATLAB function gaminv
%               Example - 'nocheck',1
%               Data Types - double
%
%
%  Output:
%
%    x:         inverse CDF value. Scalar, vector or matrix or 3D array of the same size
%               of input arguments p, a and b. $p=\int_0^x f_{IG}(t | a,b) dt$ is the
%               inverse of the inverse-gamma cdf with shape parameters in a
%               and scale parameters in b for the corresponding
%               probabilities in p
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
% inversegaminv computes the inverse of the inverse-gamma cdf with shape
% parameters in a and scale parameters in b for the corresponding
% probabilities in p. p, a, and b can be vectors, matrices, or
% multidimensional arrays that all have the same size. A scalar input is
% expanded to a constant array with the same dimensions as the other
% inputs. The parameters in a and b must all be positive, and the values in
% x must lie on the interval $[0,  \infty)$.
%
% Note that $F_{IG}(x,a,b)=\frac{\Gamma(a,b/x)}{\Gamma(\alpha)}$ therefore
% Therefore, the CDF for an inverse Gamma distribution can be computed
% using the incomplete gamma function (also called regularized gamma
% function, i.e. MATLAB function gammainc) of course  keeping into account
% that we need the upper tail.
%
%
% The chief use of the inverse gamma distribution is in Bayesian
% statistics, where the distribution arises as the marginal posterior
% distribution for the unknown variance of a normal distribution if an
% uninformative prior is used; and as an analytically tractable conjugate
% prior if an informative prior is required.
% Relation with the Gamma distribution.
% If $X \sim Gamma(a,b)$ then $\frac{1}{X} \sim$ inverse-gamma distribution
% with paramters $a$ and $1/b$
%
% See also: gampdf, inversegampdf, inversegamcdf
%
% References:
%
% https://en.wikipedia.org/wiki/Inverse-gamma_distribution
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('inversegaminv')">Link to the help function</a>
% Last modified 06-Feb-2015
%

% Examples:
%

%{
    %% Compare the results using option nocheck=1.
    x=(0:0.0001:0.9999)';
    a=[1,2,3,50,100,10000];
    b=[1,10,100,0.05,10,800];
    Y=zeros(length(x),length(a));
    Ychk=Y;
    for i=1:length(x)
        Y(i,:)=inversegaminv(x(i),a,b);
        Ychk(i,:)=inversegaminv(x(i),a,b,1);
    end
    disp('Maximum absolute difference is:')
    disp(max(max(abs(Y-Ychk))))
%}

%{
    %% Check accuracy of results monitoring $|p-F_{IG} (F_{IG}^{-1}(p))|$.
    a=[1,2,3,50,100,10000];
    b=[1,10,100,0.05,10,800];

    x=(0:0.0001:0.9999)';
    Y=zeros(length(x),length(a));
    Ychk=Y;

    for i=1:length(x)
        Y(i,:)=x(i)-inversegamcdf(inversegaminv(x(i),a,b),a,b);
        Ychk(i,:)=x(i)-inversegamcdf(inversegaminv(x(i),a,b,1),a,b,1);
    end
    disp('Maximum deviation from 0 passing through routine gaminv:')
    disp(max(max(abs(Y))))
    disp('Maximum deviation from 0 using fast routine:')
    disp(max(max(abs(Ychk))))
%}

%{
    % Compare results with R (library actuar).
    % The example below assumes that the Connection with R has already been setup
    % For more information on how to connect R and Matlab see file
    % Connect_Matlab_with_R_HELP
    % in folder
    % disp(which('Connect_Matlab_with_R_HELP'))
    chkMatlab_With_R_connection=exist('openR','file');
    if chkMatlab_With_R_connection==0
        disp('Connection with R has not been setup yet')
        examp=which('Connect_Matlab_with_R_HELP.m');
        examp1=strrep(examp,'\','\\');
        stri=['See file <a href="matlab: opentoline(' examp1 ',27)">Connect_Matlab_with_R_HELP.m</a>  for more information'];
        disp(stri)
    else
        openR
        evalR('library(actuar)');
        evalR('x=seq(0,0.9999,0.0001)')
        evalR('a=2')
        evalR('b=3')
        yfromR=evalR('qinvgamma(x,a,scale=b)');
        x=(0:0.0001:0.9999);
        yfromMatlab=inversegaminv(x,2,3);
        disp(max(abs(yfromR-yfromMatlab)))
        closeR
    end
%}

%% Beginning of code

if nargin<4
    nocheck=0;
end

if nocheck==1
    % This code is much faster but there are no checks
    x = b./gammaincinv(p,a,'upper');
else
    x = 1./gaminv(1 - p,a,1./b);
end

end
