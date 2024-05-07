function r = tobitrnd(mu,sigma, left, right,   mm,nn,oo)
%tobitrnd random arrays from the tobit distribution
%
%   returns an array of random numbers of size mm-by-nn-by-oo chosen from the
%   tobit distribution with parameters mu, sigma, left, right
%
%<a href="matlab: docsearchFS('tobitrnd')">Link to the help page for this function</a>
%
%  Required input arguments:
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
%               Example - 'sigma',1
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
% Optional input arguments:
%
%         mm    : Length of first dimension. Scalar. Number of rows of the
%                 array which contains the random numbers.
%               Example - 3
%               Data Types - double
%         nn    : Length of second dimension. Scalar. Number of columns of
%                 the array which contains the random numbers.
%               Example - 2
%               Data Types - double
%         oo    : Length of third dimension. Scalar. Number of 3D slides of
%                 the array which contains the random numbers.
%               Example - 5
%               Data Types - double
%       accuracy : accuracy of the calculations. Scalar. The default value
%                  of accuracy is 1e-10.
%                  Data Types - single|double
%                  Example - 1e-06
%
% Output:
%
%         r    : Random numbers. Array of random numbers from the
%                tobit distribution with paramters, mu, sigma, left and right.
%                The size of rr is determined by the optional input
%                parameters mm, nn, oo.
%
%
% See also: tobitpdf, tobitcdf, tobitinv
%
%
% References:
%
% Greene, W.H. (2008), "Econometric Analysis, Sixth Edition", Prentice Hall, pp. 871-875.
%
% Tobin, J. (1958), Estimation of Relationships for Limited Dependent
% Variables, "Econometrica", 26, pp. 24-36.
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('tobitrnd')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Generate a random number from  the tobit  distribution(mu,sigma,left,right).
    mu=1; % location parameter
    sigma=2;  % dispersion parameter
    left=-0.5; % Lower limit
    right=5; % Upper limit
    % Generate a random number from the distribution above
    x=tobitrnd(mu,sigma,left,right);
    disp(x)
%}




%{
    % Generate a 2D array of size mmxnn of random number from  the Tobit distribution.
    % total number of balls
    mu=1; % location parameter
    sigma=2;  % dispersion parameter
    left=-0.5; % Lower limit
    right=5; % Upper limit
    % Generate an array of size 2x3x5 of random numbers from the distribution above
    mm=2;   nn=3;
    X=tobitrnd(mu,sigma,left,right, mm,nn);
    disp(X)
%}

%{
    % Generate a 3D array of size mmxnnxoo of random number from  the Tobit distribution.
    % total number of balls
    mu=1; % location parameter
    sigma=2;  % dispersion parameter
    left=-0.5; % Lower limit
    right=5; % Upper limit
    % Generate an array of size 2x3x5 of random numbers from the distribution above
    mm=2;   nn=3; oo=5;
    X=tobitrnd(mu,sigma,left,right, mm,nn,oo);
    disp(X)
%}

%{
    %% Compare relative frequencies with probabilities.
    % Generate data from a Tobit distribution with parameters 
    % mu sigma left right 
    mu=3; sigma=1;  left=2; right=4.2;
    % Generate n observations
    n=10000;
    x=tobitrnd(mu,sigma,left,right,n);
    % Classes are close to the left and open to the right that is [ ), except
    % for the last one which is closed to the left and to the right,
    %  that is [ ]
    edges=[left linspace(left+1e-10, right-1e-10,10) right];
    % First class goes from left (included) to left+1e-10
    % ...
    % Last class goes from right-1e-10 to right (right is included in the last class)
    % The length of edges is 12 so there are 11 classes
    % Set up string for the title of the figures
    para=['with parameters \mu,\sigma,left,right= (' num2str(mu) ',' num2str(sigma) ',' num2str(left) ',' num2str(right) ')']; 
    
    h= histogram(x,edges,'Normalization','probability');
    title('Histogram of relative frequencies from Tobit distribution',para)
    xlabel('Classes')
    
    
    freq=h.Values;
    % Create a new figure to compare theoretical and empirical frequencies
    figure
    % Compute the probabilities for each class in vector yy
    y=tobitcdf(edges,mu,sigma,left,right);
    yy=(y(2:end)-y(1:end-1)); 
    % First element of yy is the probability of tobit distribution is equal to left
    yy(1)=tobitcdf(left,mu,sigma,left,right); 
    % set the labels for the classes
    % From MATLAB 2023b the  x coordinates of bar can be a string array
    lab=string(edges(1:end-1))+"-"+string(edges(2:end));
    % The instruction below is just for compatibility with older versions of
    % MATLAB
    lab=categorical(lab,lab);
    hold('on')
    bar(lab,[freq' yy'])
    legend(["Emprical relative frequencies" "Theoretical probabilities"])
    title(['Data from Tobit distribution', para])
    xlabel('Classes')
%}



%% Beginning of code
if nargin < 4
    error(message('FSDA:tobitrnd:TooFewInputs'));
end

if nargin<7
    oo=1;
end
if nargin <6
    nn=1;
end
if nargin <5
    mm=1;
end

r=zeros(mm,nn,oo);
for i=1:numel(r)
    x=normrnd(mu,sigma);
    if x<=left
        r(i)=left;
    elseif x>=right
        r(i)=right;
    else
        r(i)=x;
    end
end


%FScategory:ProbDist

