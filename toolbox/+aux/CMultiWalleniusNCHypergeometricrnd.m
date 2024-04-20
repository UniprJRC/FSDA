function r = CMultiWalleniusNCHypergeometricrnd(m,w,n, accuracy)

% /*
% This function generates a vector of random variables with the
% multivariate Wallenius noncentral hypergeometric distribution.
%
% The multivariate Wallenius noncentral hypergeometric distribution is
% the distribution you get when drawing colored balls from an urn
% with any number of colors, without replacement, and with bias.
%
% The weights are defined so that the probability of taking a particular
% ball is proportional to its weight.
%
% Parameters:
% destination:    An output array to receive the number of balls of each
% color. Must have space for at least 'colors' elements.
% source:         An input array containing the number of balls of each
% color in the urn. Must have 'colors' elements.
% All elements must be non-negative.
% weights:        The odds of each color. Must have 'colors' elements.
% All elements must be non-negative.
% n:              The number of balls to draw from the urn.
% Cannot exceed the total number of balls with nonzero weight
% in source.
% colors:         The number of possible colors.
%
% MAXCOLORS  (defined in stocc.h): You may adjust MAXCOLORS to the maximum
% number of colors you need.
%
% The function will reduce the number of colors, if possible, by eliminating
% colors with zero weight or zero number and pooling together colors with the
% same weight. The problem thus reduced is handled in the arrays osource,
% urn, oweights and osample of size colors2.
%
% The sampling proceeds by either of two methods: simulating urn experiment,
% or conditional method followed by Metropolis-Hastings sampling.
%
% Simulating the urn experiment is simply taking one ball at a time, requiring
% n uniform random variates. The problem is reduced whenever a color has been
% exhausted.
%
% The conditional method divides the colors into groups where the number of
% balls in each group is determined by sampling from the marginal distribution
% which is approximated by the univariate Wallenius distribution. Each group
% is then subdivided by sampling one color at a time until all colors have
% been sampled.
%
% The sample from the conditional method does not have the exact distribution,
% but it is used as a starting point for the Metropolis-Hastings sampling,
% which proceeds as follows: colors c1 and c2 are re-sampled using the
% univariate Wallenius distribution, keeping the samples of all other colors
% constant. The new sample is accepted or the old sample retained, according
% to the Metropolis formula which corrects for the slight error introduced
% by not using the true conditional distribution. c1 and c2 are rotated in
% an order determined by the variance of each color. This rotation (scan) is
% repeated nHastings times.
% */

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
    error('FSDA:CMultiWalleniusNCHypergeometricrnd:WrgInp',"Parameter out of range in function MultiWalleniusNCHyp");
end

destination=zeros(colors,1);

if (n == 0)
    return
end

% // check validity of array parameters
if (min(m) < 0 || min(w) < 0)
    error('FSDA:CMultiWalleniusNCHypergeometricrnd:WrgInp',"Parameter negative in function MultiWalleniusNCHyp");
end


N = sum(m);

%   // sort colors by weight, heaviest first
order1=1:colors;
order3=order1;
weights=w;

% for i=1:colors
%     c = order1(i);
%     k = i;
%     w = weights(c);
%     if (source(c)==0)
%         w = 0;       %  // zero number treated as zero weight
%     end
%
%     for j=i+1:colors
%         c2 = order1(j);
%         if weights(c2) > w && source(c2)
%             w = weights(c2);  k = j;
%         end
%     end
%     order1(i) = order1(k);
%     order1(k) = c;
% end
[~,order1]=sort(weights,'descend');

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
        error('FSDA:CMultiWalleniusNCHypergeometricrnd:WrgInp',"Taking more items than there are in function MultiWalleniusNCHyp");
    end

    for i = 1:colors
        c = order1(i);
        destination(c) = source(c);
    end
    return
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
        oweights(c2) = weights(c);
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
        x = aux.WalleniusNCHypergeometricrnd(n, osource(1), N, oweights(1)/oweights(2));
        osample(1) = x;
        osample(2) = n - x;
    end
else
    %  more than 2 colors
    nn = n;
    osample=zeros(colors2,1);

    % // decide which method to use
    if (nn < 5000 * colors2)

        %// Simulate urn experiment

        %// Make list of accumulated probabilities of each color
        wcum=cumsum(urn.*oweights);


        % // take one item nn times
        j = colors2;
        while nn>0

            % // get random color according to probability distribution wcum
            randomnumber=unifrnd(0,1);
            p =  randomnumber * wcum(colors2);
            % // get color from search in probability distribution wcum
            for i=1:j
                if (p < wcum(i))
                    notbreak=0;
                    break
                else
                    notbreak=1;
                end
            end

            i=i+notbreak;
            % // sample one ball of color i
            osample(i)= osample(i) +1;
            urn(i)=urn(i)-1;  nn=nn-1;

            % // check if this color has been exhausted
            if (urn(i) == 0)
                if (i ~= (j+1))
                    %// put exhausted color at the end of lists so that colors2 can be reduced
                    m = osource(i); osource(i) = osource(j); osource(j) = m;
                    m = urn(i); urn(i) = urn(j); urn(j) = m;
                    m = osample(i); osample(i) = osample(j); osample(j) = m;
                    w = oweights(i); oweights(i) = oweights(j); oweights(j) = w;
                    % // update order2 list (no longer sorted by weight)
                    for k=1:colors
                        if order2(k) == i
                            order2(k) = j;
                        else
                            if (order2(k) == j)
                                order2(k) = i;
                            end
                        end
                    end
                end
                colors2=colors2-1;  j = colors2-1;
                % // decrement number of colors left in urn

                if (colors2 == 2 && nn > 50)
                    % // two colors left. use univariate distribution for the rest
                    x = WalleniusNCHyp(nn, urn(1), urn(1)+urn(2), oweights(1)/oweights(2));
                    osample(1) = osample(1) +x;
                    osample(2) = osample(2)+ nn - x;
                    break
                end

                if (colors2 == 1)
                    %// only one color left. The rest is deterministic
                    osample(1) = osample(1) +nn;
                    break
                end

                % // make sure wcum is re-calculated from beginning
                i = 1;
            end


            % // update list of accumulated probabilities
            if i > 1
                wsum = wcum(i-1);
            else
                wsum=0;
            end

            for k=i:colors2
                wsum = wsum +urn(k) * oweights(k);
                wcum(k) = wsum;
            end
        end

    else
        %// use conditional method to make starting point for
        %// Metropolis-Hastings sampling

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
        n1 = aux.WalleniusNCHypergeometricrnd(n, m1, m1+m2, w1/w2);
        % n1=WNChygernd(m1+m2,m1, n, w1/w2);
        n2 = n - n1;

        % // set parameters for first group (heavy)
        a = 0;  ng = n1;


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

                % // sample color i in group
                if    wsum ~=0
                    x=aux.WalleniusNCHypergeometricrnd(ng, m, msum + m, w * msum / wsum);


                    % if k==1
                    %     x=8062; % TODO
                    % else
                    %     x=13778;
                    % end


                else
                    x=ng;
                end

                osample(i+1) = x;
                ng = ng-x;
            end

            if isempty(i)
                osample(1)=ng;
            else
                osample(i+2)=ng;
            end

            % // set parameters for second group (light)
            a = b;  b = colors2;  ng = n2;
        end



        %// finished with conditional method.
        %// osample contains starting point for Metropolis-Hastings sampling

        % // calculate mean
        vari=CMultiWalleniusNCHypergeometricmean(osource,oweights,n); %

        if length(vari)<length(source)
            vari=[vari; zeros(length(source)-length(vari),1)];
        end
        % // calculate approximate variance from mean
        for i=1:colors2

            r1 = vari(i) * (osource(i)-vari(i));
            r2 = (n-vari(i))*(vari(i)+N-n-osource(i));
            if (r1 <= 0. || r2 <= 0.)
                vari(i) = 0;
            else
                vari(i) = N*r1*r2/((N-1)*(osource(i)*r2+(N-osource(i))*r1));
            end

        end

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

        % // number of scans (this value of nHastings has not been fine-tuned)
        nHastings = 4;
        if (accuracy < 1E-6)
            nHastings = 6;
        end

        if (colors2 > 5)
            nHastings=nHastings+1;
        end

        % // Metropolis-Hastings sampler
        f0 = -1.;
        for k =1:nHastings
            for i = 1:colors2
                j = i+1;
                if (j >= colors2)
                    j = 1;
                end
                c1 = order3(i);  c2 = order3(j);
                w = oweights(c1) / oweights(c2);
                n1 = osample(c1) + osample(c2);
                x0 = osample(c1);
                x = aux.WalleniusNCHypergeometricrnd(n1, osource(c1), osource(c1)+osource(c2), w);
                if (x == x0)
                    continue % // accepted
                end

                if (f0 < 0)
                    f0 = aux.CMultiWalleniusNCHypergeometricpdf(osample, osource,oweights, accuracy);
                end


                g0 = aux.WalleniusNCHypergeometricpdf(x0, n1, osource(c1), osource(c1)+osource(c2), w, accuracy);
                g1 = aux.WalleniusNCHypergeometricpdf(x,  n1, osource(c1), osource(c1)+osource(c2), w, accuracy);
                osample(c1) = x;
                osample(c2) = n1 - x;
                f1 = aux.CMultiWalleniusNCHypergeometricpdf(osample, osource,oweights, accuracy);
                g0 = f1 * g0;  g1 = f0 * g1;
                if (g0 >= g1 || g0 > g1 * unifrnd(0,1))
                    % / new state accepted
                    f0 = -1.;
                else
                    % // rejected. restore old sample
                    osample(c1) = x0;
                    osample(c2) = n1 - x0;
                end
            end
        end
    end
end


%// finished sampling by either method
% // un-sort sample into destination and untangle re-orderings
for i=1:colors
    c1 = order1(i);  c2 = order2(i);
    if source(c1) == osource(c2)
        destination(c1) = osample(c2);
    else
        %// split colors with same weight that have been treated as one
        x = aux.WalleniusNCHypergeometricrnd(osample(c2), source(c1), osource(c2));
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
        error('FSDA:CMultiWalleniusNCHypergeometricrnd:WrgIter',"Search for mean failed in function CMultiWalleniusNCHypergeometric::mean");
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
