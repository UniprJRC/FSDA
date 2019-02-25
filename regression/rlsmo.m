function [smo,span]=rlsmo(x,y,w,span)
% rlsmo computes a running-lines smoother with global cross-validation.
% 
%<a href="matlab: docsearchFS('rlsmo')">Link to the help page for this function</a>
%
%
%   This function is called in each step of the avas function but it
%   can be called directly when it is necessary to smooth a set of values
%   using local regressions. Note that the x values must be non decreasing.
%
%  Required input arguments:
%
%    x   :      Predictor variable sorted. Vector.  Ordered abscissa values.
%               Note that the x values are assumed non decreasing.
%    y  :       Response variable. Vector. Response variable which has to
%               be smoothed, specified as
%               a vector of length n, where n is the number of
%               observations.
%
%
%  Optional input arguments:
%
%       w  : weights for the observations. Vector. Row or column vector of
%           length n containing the weights associated to each
%           observations. If w is not specified we assum $w=1$ for $i=1,
%           2, \ldots, n$.
%           Example - 1:n
%           Data Types - double
%  span    : length of the local regressions. Scalar. Scalar in the
%           interval [0, 1] which specifies the length of the local
%           regressions. If span is 0 (default value) the fractions of
%           observations which are considered for computing the local
%           regerssions are roughly $cvspan=n*[0.3,0.4,0.5,0.6,0.7,1.0]$.
%           The element of $cvspan$ which is associated with the smallest
%           cross validation residual sum of squares is chosen. The
%           smoothing procedure is called using the best value of cvspan
%           and the smoothed values are found without cross validation.
%           If span is not 0 but is a value in
%           the interval (0, 1], the local regression have length n*span and
%           the smoothed values are found without cross validation.
%           Example - 0.4
%           Data Types - double
%
%  Output:
%
%         ysmo:  smoothed values. Vector. A vector with the same dimension
%               of y containing smoothed values, that is the y values on
%               the fitted curve. The smoothed values come from linear
%               local linear regressions whose length is specified by input
%               parameter span.
%         span: length of the local regressions. Scalar. Scalar in the
%               interval [0, 1] which specifies the length of the local
%               regressions which has been used. For example if span=0.3
%               approximately 30 per cent of consecutive obsrvations are
%               used in order to compute the local regressions.
%
%
% More About:
%
% This function makes use of subroutine smth.
% The sintax of $smth$ is $[smo] = smth(x,y,w,span,cross)$. $x$, $y$ and
% $w$ are 3 vectors of length $n$ containing respectively the $x$
% coordinates, the $y$ coordinates and the weights. Input paramter $span$ is
% a scalar in the interval (0 1] which defines the length of the elements
% in the local regressions.
% More precisely, if $span$ is in (0 1), the length of elements in the
% local regressions is $m*2+1$, where $m$ is defined as the $\max([(n
% \times span)/2],1)$ to ensure that minimum length of the local
% regression is 3. Symbol $[ \cdot ]$ denotes the integer part.
% 
% Parameter $cross$ is a Boolean scalar. If it is set to true it specifies
% that, to compute the local regression centered on unit $i$, unit $i$ must
% be deleted. Therefore for example, 
% [1] if $m$ is 3 and $cross$ is true, the
% smoothed value for observation $i$ uses a local regression with $x$
% coordinates $(x(i-1), x(i+1))$, $y$ coordinates $(y(i-1), y(i+1))$ and
% $w$ coordinates  $(w(i-1), w(i+1))$, $i=2, \ldots, n-1$. The smoothed
% values for observation 1 is $y(2)$ and the smoothed value for observation
% $n$ is $y(n-1)$. 
% [2] If $m$ is 3 and $cross$ is false, the smoothed value for
% observations $i$ is based on a local regression with $x$ coordinates
% $(x(i-1), x(i), x(i+1))$, $y$ coordinates $(y(i-1), y(i), y(i+1))$ and
% $w$ coordinates  $(w(i-1), w(1), w(i+1))$, $i=2, \ldots, n-1$. The
% smoothed values for observation 1 uses a local regression based on
% $(x(1), x(2))$, $(y(1), y(2))$, and  $(w(1), w(2))$ while the smoothed
% value for observation $n$ uses a local regression based on $(x(n-1),
% x(n))$, $(y(n-1), y(n))$, and  $(w(n-1), w(n))$.
% 
% [3] If $m=5$ and $cross$ is true, the smoothed value for observations $i$
% uses a local regression based on observations $(i-2), (i-1), (i+1),
% (i+2)$, for $i=3, \ldots, n-2$.  The smoothed values for observation 1
% uses observations 2 and 3, the smoothed value for observations 2 uses
% observations 1, 3 and 4 ... 
% [4] If $m$ is 5 and $cross$ is false, the
% smoothed value for observations $i$ uses a local regression based on
% observations $(i-2), (i-1), i, (i+1), (i+2)$, for $i=3, \ldots, n-2$.
% The smoothed values for observation 1 uses observations 1, 2 and 3, the
% smoothed value for observations 2 uses observations 1, 2, 3 and 4 ...
%
%
% See also: avas.m, smothr.m, ace.m, supsmu.m
%
% References:
% 
% Tibshirani R. (1987), Estimating optimal transformations for regression,
% "Journal of the American Statistical Association", Vol. 83, 394-405.
%  Hastie, T., and Tibshirani, R. (1986), Generalized Additive Models (with
%  discussion), "Statistical Science", Vol 1, pp. 297-318
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('rlsmo')">Link to the help function</a>
%
%$LastChangedDate:: 2018-10-10 19:53:34 #$: Date of the last commit


% Examples:

%{
    % rlsmo with all the default arguments.
    n=200;
    x=sort(randn(n,1));
    y=2*x.^2+-3*x+2*randn(n,1);
    [ysmo,span]=rlsmo(x,y);
    plot(x,[y ysmo])
    title(['span chosen by cross validation= ' num2str(span)])
%}

%{
    % rlsmo with weights.
    n=200;
    x=sort(randn(n,1));
    y=2*x.^2+-3*x+2*randn(n,1);
    w=1:n; w=w(:);
    [ysmo,span]=rlsmo(x,y,w);
    plot(x,[y ysmo])
    title(['span chosen by cross validation= ' num2str(span)])
%}


%{
    % rlsmo with span value supplied as input.
    n=200;
    x=sort(randn(n,1))*10;
    y=3*x.^3-2*x.^2+-4*x+10000*randn(n,1);
    [ysmo,span]=rlsmo(x,y,[],0.5);
    plot(x,[y ysmo])
    title(['Fixed value of span = ' num2str(span)])
%}


%% Beginning of code

n=length(x);
if nargin<3 || isempty(w)
    w=ones(n,1);
end

if nargin<4
    span=0;
end

if span == 0
    cross = 1;
else
    cross = 0;
end

penal=0.01;
cvmin=1e15;

if cross==1
    cvspan=[0.3,0.4,0.5,0.6,0.7,1.0];

    % cvrss = vector which will contain the cross validation residual sum
    % of squares for each value of vector cvspan
    cvrss=zeros(6,1);
    
    for  k=1:6
        % [smo,rss]=smth(x,y,w,span,dof,n,cross);
        % disp(['k=' num2str(k)])
        % Apply the smoothing procedure for a given value of cvpsan
        % Note that smth is called with input parameter cross set to 1
        % (that is cross validation is used to find smoothed values)
        [smo,cvrss(k)]=smth(x,y,w,cvspan(k),1);
        
        if cvrss(k)<=cvmin
            cvmin=cvrss(k);
            idmin=k;
        end
    end
    % Store initial best value of span
   span=cvspan(idmin);
   
   % Penalize cvmin with a penalty
        cvmin= (1.+penal)*cvmin;
        
        % If it is found that a cv residual sum of squares is smaller than
        % penalized cvmin then use that value of cvspan
        for k=6:-1:1
            if cvrss(k) <= cvmin
                span=cvspan(k);
                break
            end
        end
end

% Call smth function with cross validation parameter set to 0
[smo,~,meay]=smth(x,y,w,span,0);
smo=smo+meay;
end


function [smo,rss,meay] = smth(x,y,w,span,cross)
% smoothing function for aperiodic data, uses weights.

% The output of smth is 
% smo = the vector which contains the smoothed values
% coming from local regressions. Note that this vector is given in terms of
% deviations from the mean 
% rss = the vector of residual sum of squares or
% cross validation residual sum of squares (if smth is call with input
% parameter cross set to true). meay = the mean of the y values.

% Get the dimensions of the input.
n = length(y);

wy = w.*y;
wxy = wy.*x;
data = [cumprod([w,x,x],2),wy,wxy];

% FS: The alternative code is
% data =[w,x, x.^2, wy, wxy]
sw=sum(w);
meay=sum(y.*w)/sw;


if span <1
    ispan=n.*span;
    m=fix(ispan./2);
    
    if m<1
        m=1;
    end
    k = 2*m + 1;
    
    % Compute sum(w), sum(w*x), sum(w*x^2), sum(w*y), sum(w*x*y) over k points.
    sums = zeros(n,5);
    % ----------------------------------------------
    % Slower, more accurate code:
    % temp = filter(ones(it,1),1,data);
    % sums(m+1:n-m,:) = temp(k:end,:);
    % ----------------------------------------------
    % Faster, slightly less accurate code:
    cs = [0 0 0 0 0;cumsum(data)];
    sums(m+1:n-m,:) = cs(k+1:end,:) - cs(1:end-k,:);
    % ----------------------------------------------
    % Repeat row m+1 in rows 1:m
    % and repeat row n-m in the last m rows
    sums(1:m,:) = sums((m+1)*ones(1,m),:);
    sums(n-m+1:end,:) = sums((n-m)*ones(1,m),:);
    
    remfromtop=flipud(cumsum(flipud(data(m+2:2*m+1,:)),1));
    sums(1:m,:)=sums(1:m,:)-remfromtop;
    
    remfrombottom=(cumsum((data(n-2*m:n-m-1,:)),1));
    sums(n-m+1:end,:)=sums(n-m+1:end,:)-remfrombottom;
    
else
    sums=repmat(sum(data,1),n,1);
end

% To remove from each row of sums observation i (that is to do cross
% validation)
if cross==1
    sums=sums-data;
end

denom = sums(:,1).*sums(:,3) - sums(:,2).^2;
a = (sums(:,4).*sums(:,3) - sums(:,2).*sums(:,5))./denom;
b = (sums(:,1).*sums(:,5) - sums(:,2).*sums(:,4))./denom;
smo = a + b.*x;
if span<1
    if m==1
        smo(1)=y(2);
        smo(end)=y(end-1);
    end
    
    % Check the presence of NaN inside smo
    % NaN are due to constant values of x over the span
    NaNsmo=isnan(smo);
    
    % The smoothed values are simply equal to the weighted average of y over
    % the span, if x is constant over the span
    if sum(NaNsmo)>0
        smo(NaNsmo)=sums(NaNsmo,4)./sums(NaNsmo,1);
    end
    % Return smoothed values in terms of deviation from the overall mean of y
    smo=smo-meay;
else
    smo=smo-sums(:,4)./sums(:,1);
end

 rss=sum((y-meay-smo).^2.*w)/sw;

% rss=sum((y-smo).^2.*w)/sw;
% dd=1;

end



