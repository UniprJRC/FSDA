function Wcdf = FisherNCHypergeometriccdf(x,n,m,N,omega,accuracy)
%WNChygepdf returns Fisher non-central hypergeometric probability density function

if (n < 0 || n > N || m < 0 || m > N || omega < 0)
    error("Parameter out of range in CWalleniusNCHypergeometric");
end


if nargin <6
    accuracy=1e-10;
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

if x>=xmax
    Wcdf=1;
elseif x<xmin
    Wcdf=0;
else
    Wcdf=0;
    for i=xmin:1:x
        Wcdf=Wcdf+FisherNCHypergeometricpdf(x,n,m,N,omega,accuracy);
    end
end

%FScategory:ProbDist