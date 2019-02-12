function smo = supsmu(x,y,varargin)
%supsmu: Smoothing of scatterplots using Friedman's supersmoother algorithm.
%
% 
%<a href="matlab: docsearchFS('supsmu')">Link to the help page for this function</a>
%
%
%   This function implements the supersmoother. This is basically the 
%   function supsmu written in MALTAB by Douglas M. Schwarz.
%    Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
%    Real_email = regexprep(Email,{'=','*'},{'@','.'})
%    See the section "More About" of this file for the details of 
%    of the modifications which have been made.
%
%
% Required input arguments:
%
%    x   :      x values for smoothing. Vector.  
%               Predictor variable containing the abscissa values specified
%               as a vector of length n, where n is the number of
%               observations.
%    y  :       y values for smoothing. Vector. 
%               Response variable which has to be smoothed, specified as a
%               vector of length n, where n is the number of observations.
%
%  Optional input arguments:
%
%       Weights  : weights for the observations. Vector. Row or column vector of
%                positive numbers of length n containing the weights associated to each
%                observations. If w is not specified we assum $w=1$ for $i=1,
%                2, \ldots, n$.
%                Example - 'Weights',1:n
%                Data Types - double
%       Span   :  fraction of the observations in the span. Scalar.
%               This option sets the width of a fixed-width smoothing operation
%               relative to the number of data points, 0 < Span < 1.
%               Setting this to be non-zero disables the supersmoother
%               algorithm.  Default is 0 (use supersmoother).
%                Example - 'Span',0.2
%                Data Types - double
%     Period  : Sets the period of periodic data.  Default is Inf
%               (infinity) which implies that the data is not periodic.
%               Can also be set to zero for the same effect.
%                Example - 'Period',1
%                Data Types - double
%   Alpha       Sets a small-span penalty to produce a greater smoothing
%               effect.  0 < Alpha < 10, where 0 does nothing and 10
%               produces the maximum effect.  Default = 0. This paramter
%               controls the smoothness of the fitted curve. Values up to
%               10 indicate increasing smoothness.
%                Example - 'Alpha',5
%                Data Types - double
%   Unsorted    Sorted or unsorted data. Boolean. 
%               If the data points are not already sorted in order of the x
%               values then setting this to true will sort them.
%               Default is false.
%                Example - 'Unsorted',5
%                Data Types - true
%
%  Output:
%
%         smo:  smoothed values. Vector. A vector with the same dimension
%               of y containing smoothed values, that is the y values on
%               the fitted curve.
%
%
% More About:
%
%
% The supersmoother algorithm computes three separate smooth curves from
% the input data with symmetric spans of 0.05*n, 0.2*n and 0.5*n, where n
% is the number of data points.  The best of the three smooth curves is
% chosen for each predicted point using leave-one-out cross validation. The
% best spans are then smoothed by a fixed-span smoother (span = 0.2*n) and
% the prediction is computed by linearly interpolating between the three
% smooth curves.  This final smooth curve is then smmothed again with a
% fixed-span smoother (span = 0.05*n).
%
% According to comments by Friedman, "For small samples (n < 40) or if
% there are substantial serial correlations between observations close in
% x-value, then a prespecified fixed span smoother (span > 0) should be
% used.  Reasonable span values are 0.2 to 0.4."
%
%    This function is basically equal to the 
%   function supsmu written in MALTAB by Douglas M. Schwarz.
%    Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
%    Real_email = regexprep(Email,{'=','*'},{'@','.'})
%
%    The following modifications with respect to the original function have
%    been made:
%    [1] In case of constant values of x over a span program was producng
%    NA. Modifications have been done in the subroutinesto cope with this
%    case. 
%    [2] The formula for the cross validation residuals in presence of
%    constant values has been introduced.
%    [3] All the Mlint suggestions have been incorporated. 
%    [4] The help has been put inside the FSDA style.
%    [5] A series of examples which explore the different options has been
%    added.
%    
%
% 
% 
% See also: ace.m
%
% References:
%
% Friedman, J. H. (1984). A Variable Span Smoother. Tech. Rep. No. 5,
% Laboratory for Computational Statistics, Dept. of Statistics, Stanford
% Univ., California.
%
%

%{
    % Example of use of supsmu with all the default options.
    x = linspace(0,1,201);
    y = sin(2.5*x) + 0.05*randn(1,201);
    smo = supsmu(x,y);
    plot(x,y,'o',x,smo)
%}



%% Beginning of code

% Input checks.
narginchk(2,Inf)

% x and y must be vectors with same number of points (at least 5).
sy = size(y);
n = length(y);
if ~isvector(x) || ~isvector(y) || length(x) ~= n || n < 5
    error('FSDA:supsmu:WrongInputOpt','x and y must be equal-length vectors of at least 5 points.')
end

% Define properties and set default values.
prop.weights = [];
prop.span = 0;
prop.period = Inf;
prop.alpha = 0;
prop.unsorted = false;

% Process inputs and set prop fields.
properties = fieldnames(prop);
arg_index = 1;
while arg_index <= length(varargin)
    arg = varargin{arg_index};
    if ischar(arg)
        prop_index = find(strncmpi(arg,properties,length(arg)));
        if length(prop_index) == 1
            prop.(properties{prop_index}) = varargin{arg_index + 1};
        else
            error('FSDA:supsmo:WrongInputOpt','Property ''%s'' does not exist or is ambiguous.',arg)
        end
        arg_index = arg_index + 2;
    elseif isstruct(arg)
        arg_fn = fieldnames(arg);
        for i = 1:length(arg_fn)
            prop_index = find(strncmpi(arg_fn{i},properties,...
                length(arg_fn{i})));
            if length(prop_index) == 1
                prop.(properties{prop_index}) = arg.(arg_fn{i});
            else
                error('FSDA:supsmo:WrongInputOpt','Property ''%s'' does not exist or is ambiguous.',...
                    arg_fn{i})
            end
        end
        arg_index = arg_index + 1;
    else
        error('FSDA:supsmo:WrongInputOpt',['Properties must be specified by property/value pairs',...
            ' or structures.'])
    end
end

% Validate Weights property.
if isempty(prop.weights)
elseif length(prop.weights) == 1
    prop.weights = [];
elseif isvector(prop.weights) && length(prop.weights) == n
    prop.weights = prop.weights(:);
else
    error('FSDA:supsmo:WrongInputOpt','Weights property must be a vector of the same length as X and Y.')
end

% Validate Span property.
if ~isscalar(prop.span) || prop.span < 0 || prop.span >= 1
    error('FSDA:supsmu:WrongInputOpt','Span property must be a scalar, 0 <= span < 1.')
end

% Validate Periodic property.
if ~isscalar(prop.period) || prop.period < 0
    error('FSDA:supsmu:WrongInputOpt','Periodic property must be a scalar >= 0.')
end
if isinf(prop.period)
    prop.period = 0;
end

% Validate Alpha property.
if isscalar(prop.alpha)
    prop.alpha = min(max(prop.alpha,0),10);
else
    error('FSDA:supsmo:WrongInputOpt','Alpha property must be a scalar.')
end

% Validate Unsorted property.
if ~isscalar(prop.unsorted)
    error('FSDA:supsmo:WrongInputOpt','Unsorted property must be a scalar.')
end


% Select one of four smooth functions.  Each smooth function has been
% speed optimized for the specific conditions.
smooth_fcn_selector = 2*isempty(prop.weights) + (prop.period == 0);
switch smooth_fcn_selector
    case 0 % use weights, periodic
        smooth = @smooth_wt_per;
    case 1 % use weights, aperiodic
        smooth = @smooth_wt_aper;
    case 2 % no weights, periodic
        smooth = @smooth_per;
    case 3 % no weights, aperiodic
        smooth = @smooth_aper;
end

% Make x and y into column vectors and sort if necessary.
x = x(:);
y = y(:);
if prop.unsorted
    [x,order] = sort(x);
    y = y(order);
    if ~isempty(prop.weights)
        prop.weights = prop.weights(order);
    end
end

% If prop.span > 0 then we have a fixed span smooth.
if prop.span > 0
    smo = smooth(x,y,prop.weights,prop.span,prop.period);
    smo = reshape(smo,sy);
    if prop.unsorted
        smo(order) = smo;
    end
    return
end

spans = [0.05;0.2;0.5];
nspans = length(spans);


% i=fix(n/4);
% j=3*i;
%     scale=x(j)-x(i);
% while scale==0
%     if scale==0.0
%         if j<n
%             j=j+1;
%         end
%         if i>1
%             i=i-1;
%         end
%         scale=x(j)-x(i);
%     end
% end
% % Tolerance for the deviance of x over the span 
% % If the deviance of x over the span is smaller than vsmlsq then x is
% % considered constant and the smoothed valued is simply the weighted
% % arithmetic average over the span
% eps=0.001;
% vsmlsq=(eps*scale)^2;

% Compute three smooth curves.
smo_n = zeros(n,nspans);
acvr_smo = zeros(n,nspans);
for i = 1:nspans % row 627 of avas.f
    % The output abs_cv_res corresponds to acvr
    % abs_cv_res = absolute value of cross validated residuals
    [smo_n(:,i),abs_cv_res] = smooth(x,y,prop.weights,spans(i),prop.period);
    acvr_smo(:,i) = smooth(x,abs_cv_res,prop.weights,spans(2),prop.period);
end

% Select which smooth curve has smallest error using cross validation.
[resmin,index] = min(acvr_smo,[],2);
span_cv = spans(index);

% Apply alpha.
if prop.alpha ~= 0
    small = 1e-7;
    tf = resmin < acvr_smo(:,3) & resmin > 0;
    span_cv(tf) = span_cv(tf) + (spans(3) - span_cv(tf)).* ...
        max(small,resmin(tf)./acvr_smo(tf,3)).^(10 - prop.alpha);
end
% if (alpha.gt.0.0.and.alpha.le.10.0.and.resmin.lt.sc(j,6))
%      1     sc(j,7) = sc(j,7) +
%      2        (spans(3)-sc(j,7))*max(sml,resmin/sc(j,6))**
%      3        (10.0-alpha)
%  90   continue


% Smooth span_cv and clip at spans(1) and spans(end).
smo_span = smooth(x,span_cv,prop.weights,spans(2),prop.period);
smo_span = max(min(smo_span,spans(end)),spans(1));

% Interpolate each point.
% The block of code below does the same thing as this, but much faster:
% smo_raw = zeros(n,1);
% for i = 1:n
% 	smo_raw(i) = interp1(spans,smo_n(i,:),smo_span(i));
% end
try
    bin = sum(bsxfun(@ge,smo_span,spans(1:end-1)'),2);
catch % if bsxfun does not exist.
    bin = sum(repmat(smo_span,1,nspans-1) >= repmat(spans(1:end-1)',n,1),2);
end
dspans = diff(spans);
t = (smo_span - spans(bin))./dspans(bin);
index = (1:n)' + n*bin;
smo_raw = (1-t).*smo_n(index - n) + t.*smo_n(index);

% Apply final smooth using spans(1).
%  call smooth (n,x,sc(1,4),w,spans(1),-jper,vsmlsq,smo,h)
smo = smooth(x,smo_raw,prop.weights,spans(1),prop.period);
smo = reshape(smo,sy);
if prop.unsorted
    smo(order) = smo;
end



%-------------------------------------------------------------------------
% Subfunctions
%-------------------------------------------------------------------------

function [smo,acvr] = smooth_wt_aper(x,y,w,span,~)
% smoothing function for aperiodic data, uses weights.

% if nargin<5
%     vsmlsq=0;
% end

    
% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2); % m is  ibw in the Fortran program
k = 2*m + 1;

wy = w.*y;
wxy = wy.*x;
data = [cumprod([w,x,x],2),wy,wxy];

% FS: The alternative code is
% data =[w,x, x.^2, wy, wxy]

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
sums(1:m,:) = sums((m+1)*ones(1,m),:);
sums(n-m+1:end,:) = sums((n-m)*ones(1,m),:);

denom = sums(:,1).*sums(:,3) - sums(:,2).^2;
a = (sums(:,4).*sums(:,3) - sums(:,2).*sums(:,5))./denom;
b = (sums(:,1).*sums(:,5) - sums(:,2).*sums(:,4))./denom;
smo = a + b.*x;

% Check the presence of NaN inside smo
% NaN are due to constant values of x over the span
NaNsmo=isnan(smo);

% The smoothed values are simply equal to the weighted average of y over
% the span, if x is constant over the span
if sum(NaNsmo)>0
    smo(NaNsmo)=sums(NaNsmo,4)./sums(NaNsmo,1);
end


if nargout > 1
    sums_cv = sums - data;
    denom_cv = sums_cv(:,1).*sums_cv(:,3) - sums_cv(:,2).^2;
    a_cv = (sums_cv(:,4).*sums_cv(:,3) - sums_cv(:,2).*sums_cv(:,5))./denom_cv;
    b_cv = (sums_cv(:,1).*sums_cv(:,5) - sums_cv(:,2).*sums_cv(:,4))./denom_cv;
    smo_cv = a_cv + b_cv.*x;
    % Check the presence of NaN inside smo_cv
    % NaN are due to constant values of x over the span
    NaNsmo_cv=isnan(smo_cv);
    
    acvr = abs(smo_cv - y);
    
    if sum(NaNsmo_cv)>0
        acvr(NaNsmo_cv)=abs((y(NaNsmo_cv)-smo(NaNsmo_cv)))./(1-data(NaNsmo_cv,1)./sums(NaNsmo_cv,1));
    end
    
end

% If there are consecutive values of x which are equal than take the
% weighted arithmetic average of the corresponding smoothed values
consval=x(2:end)-x(1:end-1)==0;
if max(consval)>0
    seq=2:n-1;
    indexes=seq(consval(1:n-2)> consval (2:n-1));
    
    for i=1:length(indexes)
        if i==1
            seli=1:indexes(1);
        else
            seli=indexes(i)-sum(consval(indexes(i-1):indexes(i)-1)):indexes(i);
        end
        smo(seli)=sum(smo(seli).*w(seli))/sum(w(seli));
    end
end

%-------------------------------------------------------------------------

function [smo,acvr] = smooth_wt_per(x,y,w,span,period)
% smoothing function for periodic data, uses weights.

% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2);
k = 2*m + 1;

x = [x;x(1:k-1)+period];
y = [y;y(1:k-1)];
w = [w;w(1:k-1)];

wy = w.*y;
wxy = wy.*x;
data = [cumprod([w,x,x],2),wy,wxy];

% Compute sum(w), sum(w*x), sum(w*x^2), sum(w*y), sum(w*x*y) over k points.
% ----------------------------------------------
% Slower, more accurate code:
% temp = filter(ones(k,1),1,data);
% sums = temp(k:end,:);
% ----------------------------------------------
% Faster, slightly less accurate code:
cs = [0 0 0 0 0;cumsum(data)];
sums = cs(k+1:n+k,:) - cs(1:n,:);
% ----------------------------------------------

denom = sums(:,1).*sums(:,3) - sums(:,2).^2;
a = (sums(:,4).*sums(:,3) - sums(:,2).*sums(:,5))./denom;
b = (sums(:,1).*sums(:,5) - sums(:,2).*sums(:,4))./denom;
% smo = circshift(a + b.*x(m+1:n+m),m); % slow
smo([m+1:n,1:m],:) = a + b.*x(m+1:n+m); % fast

if nargout > 1
    sums_cv = sums - data(m+1:n+m,:);
    denom_cv = sums_cv(:,1).*sums_cv(:,3) - sums_cv(:,2).^2;
    a_cv = (sums_cv(:,4).*sums_cv(:,3) - sums_cv(:,2).*sums_cv(:,5))./denom_cv;
    b_cv = (sums_cv(:,1).*sums_cv(:,5) - sums_cv(:,2).*sums_cv(:,4))./denom_cv;
    % 	smo_cv = circshift(a_cv + b_cv.*x(m+1:n+m),m); % slow
    smo_cv([m+1:n,1:m],:) = a_cv + b_cv.*x(m+1:n+m); % fast
    acvr = abs(smo_cv - y(1:n));
end

%-------------------------------------------------------------------------

function [smo,acvr] = smooth_aper(x,y,~,span,~)
% smooth_aper(x,y,w,span,period)
% smoothing function for aperiodic data, does not weights 
% (third input argument) or period (fifth input argument).

% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2);
k = 2*m + 1;

xy = y.*x;
data = [cumprod([x,x],2),y,xy];

% Compute sum(x), sum(x^2), sum(y), sum(x*y) over k points.
sums = zeros(n,4);
% ----------------------------------------------
% Slower, more accurate code:
% temp = filter(ones(k,1),1,data);
% sums(m+1:n-m,:) = temp(k:end,:);
% ----------------------------------------------
% Faster, slightly less accurate code:
cs = [0 0 0 0;cumsum(data)];
sums(m+1:n-m,:) = cs(k+1:end,:) - cs(1:end-k,:);
% ----------------------------------------------
sums(1:m,:) = sums((m+1)*ones(1,m),:);
sums(n-m+1:end,:) = sums((n-m)*ones(1,m),:);

denom = k.*sums(:,2) - sums(:,1).^2;
a = (sums(:,3).*sums(:,2) - sums(:,1).*sums(:,4))./denom;
b = (k.*sums(:,4) - sums(:,1).*sums(:,3))./denom;
smo = a + b.*x;

% Check the presence of NaN inside smo
% NaN are due to constant values of x over the span
NaNsmo=isnan(smo);

% The smoothed values are simply equal to the weighted average of y over
% the span, if x is constant over the span
if sum(NaNsmo)>0
    smo(NaNsmo)=sums(NaNsmo,4)./sums(NaNsmo,1);
end

if nargout > 1
    sums_cv = sums - data;
    denom_cv = (k-1).*sums_cv(:,2) - sums_cv(:,1).^2;
    a_cv = (sums_cv(:,3).*sums_cv(:,2) - sums_cv(:,1).*sums_cv(:,4))./denom_cv;
    b_cv = ((k-1).*sums_cv(:,4) - sums_cv(:,1).*sums_cv(:,3))./denom_cv;
    smo_cv = a_cv + b_cv.*x;
    
        % Check the presence of NaN inside smo_cv
    % NaN are due to constant values of x over the span
    NaNsmo_cv=isnan(smo_cv);
    
    acvr = abs(smo_cv - y);
    
    if sum(NaNsmo_cv)>0
        acvr(NaNsmo_cv)=abs((y(NaNsmo_cv)-smo(NaNsmo_cv)))./(1-data(NaNsmo_cv,1)./sums(NaNsmo_cv,1));
    end

end

%-------------------------------------------------------------------------

function [smo,acvr] = smooth_per(x,y,~,span,period)
% smooth_per(x,y,w,span,period)
% smoothing function for periodic data, does not use weights (third
% argument).

% Get the dimensions of the input.
n = length(y);

m = max(round(0.5*span*n),2);
k = 2*m + 1;

x = [x;x(1:k-1)+period];
y = [y;y(1:k-1)];

xy = y.*x;
data = [cumprod([x,x],2),y,xy];

% Compute sum(x), sum(x^2), sum(y), sum(x*y) over k points.
% ----------------------------------------------
% Slower, more accurate code:
% temp = filter(ones(k,1),1,data);
% sums = temp(k:end,:);
% ----------------------------------------------
% Faster, slightly less accurate code:
cs = [0 0 0 0;cumsum(data)];
sums = cs(k+1:n+k,:) - cs(1:n,:);
% ----------------------------------------------

denom = k.*sums(:,2) - sums(:,1).^2;
a = (sums(:,3).*sums(:,2) - sums(:,1).*sums(:,4))./denom;
b = (k.*sums(:,4) - sums(:,1).*sums(:,3))./denom;
% smo = circshift(a + b.*x(m+1:n+m),m); % slow
smo([m+1:n,1:m],:) = a + b.*x(m+1:n+m); % fast

if nargout > 1
    sums_cv = sums - data(m+1:n+m,:);
    denom_cv = (k-1).*sums_cv(:,2) - sums_cv(:,1).^2;
    a_cv = (sums_cv(:,3).*sums_cv(:,2) - sums_cv(:,1).*sums_cv(:,4))./denom_cv;
    b_cv = ((k-1).*sums_cv(:,4) - sums_cv(:,1).*sums_cv(:,3))./denom_cv;
    % 	smo_cv = circshift(a_cv + b_cv.*x(m+1:n+m),m); % slow
    smo_cv([m+1:n,1:m],:) = a_cv + b_cv.*x(m+1:n+m); % fast
    acvr = abs(smo_cv - y(1:n));
end
