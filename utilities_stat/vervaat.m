function y = vervaat(betav,steps,d)
%vervaat returns a Vervaat perpetuity.
%
%<a href="matlab: docsearchFS('vervaat')">Link to the help function</a>
%
% This function allows to simulate exactly from a Vervaat perpetuity
% distribution with parameter $\beta \in (0,+\infty)$. The simulation
% method, by K. Cloud, M. Huber (2018), is accurate and very efficient: it
% is $\mathcal{O}(\beta \log(\beta))$, in the sense that it uses $T$
% uniform random variates, where $\mathbf{E}[T] \le \mathcal{O}(\beta
% \log(\beta))$. Therefore the algorithm is good even for $\beta \gg 1$.
%
% The FSDA code has been translated from the R version of the authors.
% Below we give briefly some definitions, background and motivations. In
% the references, we also include several key works on the distribution.
%
%
% Required input arguments:
% 
%
% betav : Distribution parameter value. Positive integer. The parameter of
%         the Vervaat family. Default is betav = 1.
% steps : Markov chain step. Positive integer. The inital number of steps
%         to run to move the chain forward. Default is steps = 1.
%
% Optional input arguments:
%
% d     : Chain value at time 0. Positive integer or -1. The value of the 
%         dominating chain at time 0. Default is $D_0 \leftarray x_0−1+G$,
%         where $G \sim \mbox{Geo}(1/2)$ and 
%         $x_0 = \frac{1 + (2/3)^{1/betav}}[1 - (2/3)^(1/betav)}$;
%             Example - d=2
%
% Output:
%
% y     : The Vervaat perpetuity value. Scalar. The estimated distribution 
%         value. 
%
%
% More About:
%
% A perpetuity is a random variable of the form: 
% \[ 
%  Y = W_1 + W_1 \cdot W_2 + W_1 \cdot W_2 \cdot W_3 + \ldots 
% \]
% where the $W_i$ are an independent, identically distributed (iid)
% sequence of random variables. If each $W_i$ has the same distribution,
% say $W_i \sim W$, then $Y \sim W(1 + Y)$ for $Y$ and $W$ independent.
%
% We are interested in this distribution because the running time of the
% Quickselect algorithm of Hoare, for finding the order statistics in a
% numerical array, approaches asymptotically the Dickman distribution,
% which is a perpetuity with $W \sim Unif([0,1])$. Unfortunately there is
% no known closed form for the Dickman.
%
% The Dickman distribution can be also seen as a special case of Vervaat
% perpetuity, which is such that $W_i \sim U^{1/\beta}$ for some $\beta
% \in (0,\infty)$ for $U \sim Unif([0,1])$. In other words, the Dickman
% distributon is a Vervaat perpetutiy with $\beta = 1$ (although some
% authors define it as $1 + Y$ for $Y$ being a Vervaat perpetutiy with
% $\beta = 1$).
%
%
% See also: sort.m
%
% References:
%
% K. Cloud, Huber M. (2018), Fast Perfect Simulation of Vervaat
% Perpetuities, "Journal of Complexity", Vol. 42, pp. 19-30.
%
% Devroye L. (2001), Simulating  perpetuities, "Methodology And Computing
% In Applied Probability", Vol. 3, pp. 97-115.
%
% Fill J.A.  and  Huber M.L. (2010),  Perfect  simulation  of  Vervaat
% perpetuities, "Electronic Journal of Probability", Vol. 15, pp. 96-109.
%
% Devroye L. and Fawzi O. (2010), Simulating the Dickman distribution,
% "Statistics and Probability Letters", Vol. 80, pp. 242-247.
%
% Blanchet J.H.  and Sigman K. (2011). On exact sampling of stochastic
% perpetuities, "Journal of Applied Probability", Vol. 48, pp. 165-182.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('vervaat')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-02-01 15:26:55 #$: Date of the last commit
%
% Examples:
%

%{  
    % Density of Vervaat perpetuity with default parameters.
    % First example.
    clear all;
    close all;
    betav = 1;
    y=vervaat(betav);
%}

%{
    %% Density of Vervaat perpetuity with beta = 1.
    betav = 1;
    steps = 10;
    d = 3;
    % This is just to check coherence of results with the R implementation.
    % See below.
    rn = mtR(1);    
    y=vervaat(betav,steps,d);
    % To verify coherence of results with the R implementation, execute
    % in R the following four lines, which allow to replicate the same
    % random numbers; then execute the vervaat after sourcing the following 
    % R code:
    %
    %    RNGkind("Mersenne-Twister") # set "Mersenne-Twister" "Inversion"
    %    set.seed(0)
    %    state = .Random.seed
    %    runif(1)
    %
    %    vervaat(beta = 1,steps = 10,d = 3)
    
%}

%{  
    % Densities of Vervaat perpetuity.
    % Densities are estimated for $\beta = 1$ and $\beta = 10$ using 10000 samples.
    % With this algorithm, only 23 steps per sample are required.
    betav10 = 10;
    betav01 = 1;
    steps = 10;
    d = -1;
    N = 10000;
    y10 = zeros(N,1);
    y01 = zeros(N,1);
    for i=1:N
        y10(i)=vervaat(betav10,steps,d);
        y01(i)=vervaat(betav01,steps,d);
    end
    %histogram(y); 
    %hold on;
    histfit(y10,N/100,'kernel');
    hold on;
    histfit(y01,N/100,'kernel');
%}



%% 
if nargin<3 
    d = -1;
end

if nargin<2
    steps = 1;
end

x0 = (1 + (2/3)^(1/betav)) / (1 - (2/3)^(1/betav));
if (d == -1)
    d = x0 - 1 + geornd(1/2);  %geornd generates a random number from a  
                               %geometric distribution, with probability p
end

d = [zeros(1,steps),d];
a = rand(1,steps);
for t = steps:-1:1
    d(t) = d(t + 1) + (a(t) > 2/3) - (a(t) <= 2/3) * (d(t + 1) >= x0);
end
m  = 0;
M  = d(1);
u1 = zeros(1, steps);
u2 = rand(1,steps);

for t = 2:(steps + 1)
    up = d(t) > d(t - 1);
    a = 2/3*up; b = 2/3+1/3*up;
    u1(t - 1) = a + (b-a).*rand; % random number in the interval (a,b) 
                                 % with r = a + (b-a).*rand(N,1)
end

for t = 2:(steps + 1)
    m = (1 + m) * u2(t - 1)^(1/betav);
    a = u1(t - 1)^(1/betav);
    s = (a < ((1 + m) / (1 + M)));
    M = s * m + (1 - s) * a * (1 + M);
end

if (m == M)
    y=m;
else
    % call the algorithm recursively
    y = vervaat(betav, 2 * steps, d(1));
    m = 0;
    M = d(1); %#ok<NASGU>
    for t = 2:(steps + 1)
        r = (u1(t - 1) < ((1 + m) / (1 + y))^betav);
        m = (1 + m) * u2(t - 1)^(1/betav);
        y = r * m + (1 - r) * u1(t - 1)^(1 / betav)*(1 + y);
    end
end

end

% 
% 
% # This code implements the Vervaat perpetuity algorithm for R version 3.0.2
% # (2013-09-25), following Section 2 of Cloud, Huber (2018), Fast Perfect
% # Simulation of Vervaat Perpetuities, Journal of Complexity, Volume 42,
% # October 2017, Pages 19-30.
% 
% vervaat <- function(beta = 1,
%                     steps = 1,
%                     d = -1) {
%   x0 <- (1 + (2 / 3) ^ (1 / beta)) / (1 - (2 / 3) ^ (1 / beta))
%   if (d == -1)
%     d <- x0 - 1 + rgeom(1, prob = 1 / 2)
%   d <- c(rep(0, steps), d)
%   a <- runif(steps)
%   for (t in steps:1)
%     d[t] <- d[t + 1] + (a[t] > 2 / 3) - (a[t] <= 2 / 3) * (d[t + 1] >= x0)
%   m <- 0
%   M <- d[1]
%   u1 <- rep(0, steps)
%   u2 <- runif(steps)
%   for (t in 2:(steps + 1)) {
%     up <- d[t] > d[t - 1]
%     u1[t - 1] <- runif(1, min = 2 / 3 * up, max = 2 / 3 + 1 / 3 * up)
%   }
%   for (t in 2:(steps + 1)) {
%     m <- (1 + m) * u2[t - 1] ^ (1 / beta)
%     a <- u1[t - 1] ^ (1 / beta)
%     s <- (a < ((1 + m) / (1 + M)))
%     M <- s * m + (1 - s) * a * (1 + M)
%   }
%   if (m == M)
%     return(m)
%   else {
%     y <- vervaat(beta, 2 * steps, d[1])
%     m <- 0
%     M <- d[1]
%     for (t in 2:(steps + 1)) {
%       r <- (u1[t - 1] < ((1 + m) / (1 + y)) ^ beta)
%       m <- (1 + m) * u2[t - 1] ^ (1 / beta)
%       y <- r * m + (1 - r) * u1[t - 1] ^ (1 / beta) * (1 + y)
%     }
%     return(y)
%   }
% }

%FScategory:UTISTAT
