function r = CMultiFisherNCHypergeometricrnd(m,w,n, accuracy)

% /*
%    This function generates a vector of random variates with the
%    multivariate Fisher's noncentral hypergeometric distribution.
%
%    This distribution is defined as the conditional distribution of 'colors'
%    independent binomial variates
%    x[i] = binomial(source[i], p[i])
%    on the condition that the sum of all x[i] is n.
%    p[i] = r * weights[i] / (1 + r * weights[i]),
%    r is an arbitrary scale factor.
%
%    Parameters:
%    destination:    An output array to receive the number of balls of each
%    color. Must have space for at least 'colors' elements.
%    source:         An input array containing the number of balls of each
%    color in the urn. Must have 'colors' elements.
%    All elements must be non-negative.
%    weights:        The odds of each color. Must have 'colors' elements.
%    All elements must be non-negative.
%    n:              The number of balls drawn from the urn.
%    Can't exceed the total number of balls with nonzero weight
%    in the urn.
%    colors:         The number of possible colors.
%
%    Method: The conditional method is used for generating a sample with the
%    approximate distribution. This sample is used as a starting point for
%    a Gibbs sampler. The accuracy depends on the number of scans with the
%    Gibbs sampler.
%
%    The function will reduce the number of colors, if possible, by eliminating
%    colors with zero weight or zero number and pooling together colors with the
%    same weight. A symmetry transformation is used if more than half the balls
%    are taken. The problem thus reduced is handled in the arrays osource,
%    oweights and osample of dimension colors2.
%    */

if nargin<4
    accuracy=1e-08;
end


m=m(:);
w=w(:);


% // variables
% int order1[MAXCOLORS];              // sort order, index into source and destination
% int order2[MAXCOLORS];              // corresponding index into arrays when equal weights pooled together
% int order3[MAXCOLORS];              // secondary index for sorting by variance
% int32_t osource[MAXCOLORS];         // contents of source, sorted by weight with equal weights pooled together
% int32_t urn[MAXCOLORS];             // balls from osource not taken yet
% int32_t osample[MAXCOLORS];         // balls sampled
% double oweights[MAXCOLORS];         // sorted list of weights
% double wcum[MAXCOLORS];             // list of accumulated probabilities
% double var[MAXCOLORS];              // sorted list of variance
% double w = 0.;                      // weight of balls of one color
% double w1, w2;                      // odds within group; mean weight in group
% double wsum;                        // total weight of all balls of several or all colors
% double p;                           // probability
% double f0, f1;                      // multivariate probability function
% double g0, g1;                      // conditional probability function
% double r1, r2;                      // temporaries in calculation of variance
% int32_t nn;                         // number of balls left to sample
% int32_t m;                          // number of balls of one color
% int32_t msum;                       // total number of balls of several or all colors
% int32_t N;                          // total number of balls with nonzero weight
% int32_t x0, x = 0;                  // sample of one color
% int32_t n1, n2, ng;                 // size of weight group sample or partial sample
% int32_t m1, m2;                     // size of weight group
% int i, j, k;                        // loop counters
% int c, c1, c2;                      // color index
% int colors2;                        // reduced number of colors
% int a, b;                           // color index delimiting weight group
% int nHastings;                      // number of scans in Metropolis-Hastings sampling


colors=length(m);
MAXCOLORS=32;

source=m;

% // check validity of parameters
if (n < 0 || colors < 0 || colors > MAXCOLORS)
    error('FSDA:CMultiFisherNCHypergeometricrnd:WrgInpt',"Parameter out of range in function MultiFisherNCHyp");
end

destination=zeros(colors,1);

if (n == 0)
    return
end

% // check validity of array parameters
if (min(m) < 0 || min(w) < 0)
    error('FSDA:CMultiFisherNCHypergeometricrnd:WrgInpt',"Parameter negative in function MultiFisherNCHyp");
end


N = sum(m);

%   // sort colors by weight, heaviest first
order1=1:colors;
order3=order1;
weights=w;

 % // sort by weight, heaviest first
for i=1:colors
    c = order1(i);
    k = i;
    w = weights(c);
    if (source(c)==0)
        w = 0;       %  // zero number treated as zero weight
    end

    for j=i+1:colors
        c2 = order1(j);
        if weights(c2) > w && source(c2)
            w = weights(c2);  k = j;
        end
    end
    order1(i) = order1(k);
    order1(k) = c;
end

urn=m(order1);
osource=urn;



% // skip any colors with zero weight or zero number.
% // this solves all problems with zero weights
% while (colors && (weights(c=~order1(colors-1))==0 || source(c)==0))
%     colors=colors-1;  destination(c) = 0;
% end
colors=colors-sum(weights==0);

% // check if there are more than n balls with nonzero weight
if (n >= N)
    if (n > N)
        error('FSDA:CMultiFisherNCHypergeometricrnd:WrgInpt',"Taking more items than there are in function MultiFisherNCHyp");
    end

    for i = 1:colors
        c = order1(i);
        destination(c) = source(c);
    end
    return
end


if n > (N /2) 
    % // improve accuracy by symmetry transformation
    j=colors;
    for i=1:(j-1)
        % // reverse order list
        c = order1(i);

        order1(i) = order1(j);
        order1(j) = c;
        j=j-1;

        n = N - n;
        invert = true;
    end
else
    invert=false;
end


%// copy source and weights into ordered lists
%// and pool together colors with same weight
% colors2= number of colors left
osample=zeros(colors,1);
order2=osample;
oweights=osample;

c2=0;
for i=1:(colors)
    c = order1(i);
    if (i==1 || weights(c) ~= w)
        c2=c2+1;
        x = source(c);
        if invert == true
            oweights(c2) = 1/weights(c);
        else
             oweights(c2)=weights(c);
        end

        w = weights(c);
    else
        x = x+ source(c);              %  // join colors with same weight
    end
    urn(c2) =  x;
    osource(c2) = x;
    order2(i) = c2;
    %  osample(c2) = 0;
end
colors2 = c2;
osample=zeros(colors2,1);
oweights=oweights(1:colors2,1);
osource=osource(1:colors2);
urn=urn(1:colors2);

% // check number of colors left
if (colors2 < 3)
    % // simple cases
    if (colors2 == 1)
        osample(1) = n;
    end

    if (colors2 == 2)
        x = aux.FisherNCHypergeometricrnd(n, osource(1), N, oweights(1)/oweights(2));
        osample(1) = x;
        osample(2) = n - x;
    end
else
    %  more than 2 colors
    osample=zeros(colors2,1);

    % // divide weights into two groups, heavy and light
    a = 0;  b = colors2-1;
    w = sqrt(oweights(1) * oweights(colors2));
    while (b > a + 1)
        c = floor((a + b) / 2);
        if (oweights(c+1) > w)
            a = c;
        else
            b = c;
        end
    end
    %// heavy group goes from 0 to b-1, light group goes from b to colors2-1

    % // calculate mean weight for heavy color group
    m1=0; wsum=0;
    for i=0:(b-1)
        m1 = m1+ urn(i+1);
        wsum = wsum+ oweights(i+1) * urn(i+1);
    end
    w1 = wsum / m1;

    % // calculate mean weight for light color group
    m2=0; wsum=0;
    for i=(b+1):colors2
        m2 = m2 + urn(i);
        wsum = wsum  + oweights(i) * urn(i);
    end
    w2 = wsum / m2;

    % // split partial sample n into heavy (n1) and light (n2)
    n1 = aux.FisherNCHypergeometricrnd(n, m1, m1+m2, w1/w2);
    % n1=FNChygernd(m1+m2,m1, n, w1/w2);
    n2 = n - n1;

    % // set parameters for first group (heavy)
    a = 0;  n0 = n1;


    % // loop twice, for the two groops
    for k=1:2
        % // split group into single colors by calling univariate distribution b-a-1 times
        for i = a:b-2
            m = urn(i+1);  w = oweights(i+1);

            msum=0; wsum=0;
            % // calculate mean weight of remaining colors
            for j=i+1:b-1
                m1 = urn(j+1);
                w1 = oweights(j+1);
                msum = msum +m1;
                wsum = wsum +m1 * w1;
            end

            % // split out color i
            if    w==w1
                x=aux.FisherNCHypergeometricrnd(n0, m, msum + m, w * msum / wsum);
            else
                if wsum==0
                    x=n0;
                else
                    odds= w * msum / wsum;
                    x = aux.FisherNCHypergeometricrnd(n0, m, msum + m, odds);
                end
            end

            osample(i+1) = osample(i+1) + x;
            n0 = n0-x;
        end

        if isempty(i)
            osample(a+1)=n0;
        else
            osample(i+2)=n0;
        end

        % // set parameters for second group (light)
        a = b;  b = colors2;  n0 = n2;
    end



    %// finished with conditional method.
    %// osample contains starting point for Metropolis-Hastings sampling
    % mu=CMultiFishersNCHypergeometricmean(osource,oweights,n)

    % // calculate variance
    vari=CMultiFisherNCHypergeometricvariance(osource,oweights,n); %

    % if length(vari)<length(source)
    %     vari=[vari; zeros(length(source)-length(vari),1)];
    % end


    % // sort again, this time by variance
    for i=1:colors2
        c = order3(i);  k = i;
        w = vari(c);
        for j=i+2:colors2
            c2 = order3(j);
            if (vari(c2) > w)
                w = vari(c2);  k = j;
            end
        end
        order3(i) = order3(k);  order3(k) = c;
    end

    % // number of scans (not been fine-tuned)

    ngibbs = 4;

    if (accuracy < 1E-6)
        ngibbs = 6;
    end

    if (colors2 > 5)
        ngibbs=ngibbs+1;
    end


    for k = 1:ngibbs
        for i = 1:colors2
            c1 = order3(i);
            j = i + 1;
            if (j == colors2+1)
                j = 1;
            end

            c2 = order3(j);
            n1 = osample(c1) + osample(c2);
            x = aux.FisherNCHypergeometricrnd(n1, osource(c1), osource(c1)+osource(c2), oweights(c1)/oweights(c2));
            osample(c1) = x;
            osample(c2) = n1 - x;
        end
    end

end



if invert == true
    % // reverse symmetry transformation on result
    for i=1:colors2
        osample(i) = osource(i) - osample(i);
    end
end


%// finished sampling
% // un-sort sample into destination and untangle re-orderings
for i=1:colors
    c1 = order1(i);  c2 = order2(i);
    if source(c1) == osource(c2)
        destination(c1) = osample(c2);
    else
        %// split colors with same weight that have been treated as one
        % x = aux.FisherNCHypergeometricrnd(osample(c2), source(c1), osource(c2));

        % M. K. N
    x = hygernd(osource(c2), source(c1), osample(c2));

        destination(c1) = x;
        osample(c2) = osample(c2) - x;
        osource(c2) = osource(c2) - source(c1);
    end
end

% (n,m,N,odds, accuracy)
r=destination;

end


function mu=CMultiWalleniusNCHypergeometricmean(m,omega,n)

% omega= ordered weights omega(1) is largest weight
% m corresponding populations based on ordered weights

colors=length(m);
LN2=log(2);
N=sum(m);
mu=zeros(colors,1);

if (n == 0)
    %  // needs special case
    return
end

omr=0;
% // calculate mean weight
for i=1:colors
    omr = omr+ omega(i) * m(i);
end

omr = N / omr;
% // scale weights to make mean = 1
omeg= omega * omr;


%  // Newton Raphson iteration
iter = 0;  t = -1.;                  %// first guess
H=Inf;
while (abs(H - n) > 1E-3)

    t1 = t;
    H =0;
    HD = 0;
    % // calculate H and HD
    for i = 1:colors
        if (omeg(i)~= 0.)
            [To1,To] = pow2_1(t * (1./LN2) * omeg(i));
            H = H + m(i) * To1;
            HD = HD - m(i) * omeg(i) * To;
        end
    end
    t = t -(H-n) / HD;
    if (t >= 0)
        t = 0.5 * t1;
    end
    iter=iter+1;
    if (iter > 20)
        error('FSDA:CMultiFisherNCHypergeometricrnd:WrgIter',"Search for mean failed in function CMultiWalleniusNCHypergeometric::mean");
    end
end

% // finished iteration. Get all mu[i]
for i=1:colors
    if (omeg(i) ~= 0.)
        To1 = pow2_1(t * (1./LN2) * omeg(i));
        mu(i) = m(i) * To1;
    else
        mu(i) = 0.;
    end
end
end

function [rt, r2]=pow2_1(r)
r2=2^r;
rt=1-r2;
end

function mu=CMultiFishersNCHypergeometricmean(m,odds,n)
% // calculates approximate mean of multivariate Fisher's noncentral
% // hypergeometric distribution. Result is returned in mu[0..colors-1].
% // The calculation is reasonably fast.
% // Note: The version in BiasedUrn package deals with unused colors

N=sum(m); % TOCHECK
colors=length(m);
mu=zeros(colors,1);

if (colors < 3)
    % // simple cases
    if (colors == 1)
        mu(1) = n;
    end

    if (colors == 2)
        mu(1) = CFishersNCHypergeometric(n,m(1),m(1)+m(2),odds(1)/odds(2)).mean;

        mu(2) = n - mu(0);
        return
    end
end

if (n == N)
    % // Taking all balls
    mu=m;
    return
end

% // initial guess for r
W=0;
for i=1:colors
    W = W+ m(i) * odds(i);
    r = n * N / ((N-n)*W);
end

% // iteration loop to find r
r1=Inf; iter=0;
while abs(r-r1) > 1E-5
    r1 = r;
    q=0;
    for i=1:colors
        q = q+ m(i) * r * odds(i) / (r * odds(i) + 1.);
    end
    r = r* n * (N-q) / (q * (N-n));
    iter=iter+1;
    if (iter > 100)
        error('FSDA:CMultiFisherNCHypergeometricrnd:WrgInpt',"convergence problem in function CMultiFishersNCHypergeometric::mean");
    end
end

% // store result
for i=1:colors
    mu(i) = m(i) * r * odds(i) / (r * odds(i) + 1.);
end
end


function vari=CMultiFisherNCHypergeometricvariance(m,odds,n)
% // calculates approximate variance of multivariate Fisher's noncentral
% // hypergeometric distribution (accuracy is not too good).
% // Result is returned in variance[0..colors-1].
% // The calculation is reasonably fast.
% // Note: The version in BiasedUrn package deals with unused colors
colors=length(m);
vari=zeros(colors,1);
N=sum(m);
mu=CMultiFishersNCHypergeometricmean(m,odds,n);
for i=1:colors
    r1 = mu(i) * (m(i)-mu(i));
    r2 = (n-mu(i))*(mu(i)+N-n-m(i));
    if (r1 <= 0. || r2 <= 0.)
        vari(i) = 0.;

    else
        vari(i) = N*r1*r2/((N-1)*(m(i)*r2+(N-m(i))*r1));
    end
end
end