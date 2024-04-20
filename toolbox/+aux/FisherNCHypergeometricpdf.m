function Wpdf = FisherNCHypergeometricpdf(x,n,m,N,omega,accuracy)
%FisherNCHypergeometricpdf returns Fisher non-central hypergeometric probability density function

if (n < 0 || n > N || m < 0 || m > N || omega < 0)
    error('FSDA:FisherNCHypergeometricpdf:WrgInp',"Parameter out of range in CWalleniusNCHypergeometric");
end

xmin = m + n - N;
% calculate xmin
if (xmin < 0)
    xmin = 0;
end
xmax = n;
%  calculate xmax
if (xmax > m)
    xmax = m;
end


if nargin <6
    accuracy=1e-10;
end
accur = accuracy * 0.1;

%% Beginning of code

if x<xmin || x>xmax
    if x<xmin
        disp(['x=' num2str(x) ' is smaller than xmin=' num2str(xmin) ' Density is 0'])
    end
    if x>xmax
        disp(['x=' num2str(x) ' is greater than xmax=' num2str(xmax) ' Density is 0'])
    end

    Wpdf=0;
    return
end

if n==0
    Wpdf=1;
    return
end

if xmin==xmax
    Wpdf=1;
    return
end

if omega==1 % Call the central hypergeometric
    Wpdf=hygepdf(x,N,m,n);
    return
elseif omega==0
    if n > N-m
        errorr("Not enough items with nonzero weight in CFishersNCHypergeometric::probability");
    else
        Wpdf=0;
        return
    end
end



x1 = CFishersNCHypergeometricmean(n,m,N,omega);           %  // start at mean
x1=floor(x1);

if (x1 < xmin)
    x1 = xmin;
end

x2 = x1 + 1;
scale = 0.;
scale = CFishersNCHypergeometric_lng(x1,n,m,N,omega,scale);     %/ calculate scale to avoid overflow
rsum = 1.;                       % // = exp(lng(x1)) with this scale
for xi=x1-1:-1:xmin
    y = exp(CFishersNCHypergeometric_lng(xi,n,m,N,omega,scale));   %  // sum from x1 and down
    rsum = rsum +y;
    if (y < accur)
        break;        % // until value becomes negligible
    end
end
for xi=x2:xmax
    % // sum from x2 and up
    y = exp(CFishersNCHypergeometric_lng(xi,n,m,N,omega,scale));
    rsum = rsum+y;
    if (y < accur)
        break      %  // until value becomes negligible
    end
end
rsum = 1. / rsum;             %   // save reciprocal sum

Wpdf=exp(CFishersNCHypergeometric_lng(x,n,m,N,omega,scale)) * rsum;         % // function value
end


function mea=CFishersNCHypergeometricmean(n,m,N,odds)
% // Find approximate mean

if odds == 1.                  %  // simple hypergeometric
    mea=m*n/N;
    return
end


% // calculate Cornfield mean
a = (m+n)*odds + (N-m-n);
b = a*a - 4.*odds*(odds-1.)*m*n;
if b>0
    b=sqrt(b);
else
    b=0;
end
mea = (a-b)/(2.*(odds-1.));
end

function vari=CFishersNCHypergeometricvariance(n,m,N,odds)
% // find approximate variance (poor approximation)
my= CFishersNCHypergeometricmean(n,m,N,odds);
%// find approximate variance from Fisher's noncentral hypergeometric approximation
r1 = my * (m-my);
r2 = (n-my)*(my+N-n-m);
if (r1 <= 0. || r2 <= 0.)
    vari=0;
    return
end

vari = N*r1*r2/((N-1)*(m*r2+(N-m)*r1));
if vari < 0
    vari=0;
end

end


function  logprop=CFishersNCHypergeometric_lng(x,n,m,N,odds,scale)
%// natural log of proportional function
% // returns lambda = log(m!*x!/(m-x)!*m2!*x2!/(m2-x2)!*odds^x)
x2 = n - x;
m2 = N - m;
mFac = logfactorial(m) + logfactorial(m2);

xFac = logfactorial(x) + logfactorial(x2) + logfactorial(m-x) + logfactorial(m2-x2);

logprop=mFac - xFac + x * log(odds) - scale;
end

%FScategory:ProbDist