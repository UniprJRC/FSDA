function y = vervaatsim(betav,steps,d)
%vervaatsim returns a Vervaat perpetuity.
%
%<a href="matlab: docsearchFS('vervaatsim')">Link to the help function</a>
%
% This function allows to simulate exactly from a Vervaat perpetuity
% distribution with parameter $\beta \in (0,+\infty)$. The simulation
% method, by K. Cloud, M. Huber (2018), is accurate and very efficient: it
% is $\mathcal{O}(\beta \log(\beta))$, in the sense that it uses $T$
% uniform random variates, where $\mathbf{E}[T] \le \mathcal{O}(\beta
% \log(\beta))$. Therefore the algorithm is good even for $\beta \gg 1$.
% Remark: this function is recursive and cannot be easily extended to
% the generation of several vervaat numbers.
%
% The FSDA code has been translated from the R version of the authors.
% Below we give briefly some definitions, background and motivations. In
% the references, we also include several key works on the distribution.
%
%
% Required input arguments:
%
% betav : Distribution parameter value. Positive integer. The parameter of
%         the Vervaat family. Dickman function is obtained for betav = 1.
%
% Optional input arguments:
%
% steps : Markov chain step. Positive integer. The inital number of steps
%         to run to move the chain forward. Default is steps = 1.
%         Example - steps=5.
%         Data Types - Scalar.
%
%
% d     : Chain value at time 0. Positive integer or -1. The value of the
%         dominating chain at time 0. Default is $D_0 \leftarrow x_0−1+G$,
%         where $G \sim \mbox{Geo}(1/2)$ and
%         $x_0 = \frac{1 + (2/3)^{1/\beta}}{1 - (2/3)^{1/\beta}}$
%         ($\mbox{Geo}$ represents the geometric distribution).
%         Example - d=10.
%         Data Types - Scalar.
%
% Output:
%
% y     : The Vervaat perpetuity value. Scalar. The estimated distribution
%         value. Data Types - Double.
%
%
% More About:
%
% A perpetuity is a random variable of the form:
% \[ Y = W_1 + W_1 \cdot W_2 + W_1 \cdot W_2 \cdot W_3 + \ldots \]
% where the $W_i$ are an independent, identically distributed (iid)
% sequence of random variables. If each $W_i$ has the same distribution,
% say $W_i \sim W$, then $Y \sim W(1 + Y)$ for $Y$ and $W$ independent.
%
% We are interested in this distribution because the running time of the
% Quickselect algorithm of Hoare, for finding the order statistics in a
% numerical array, approaches asymptotically a particular perpetuity,
% called Dickman distribution, with $W \sim Unif([0,1])$. Unfortunately
% such distribution has no closed form.
%
% The Dickman distribution can be also seen as a special case of Vervaat
% perpetuity, which is such that $W_i \sim U^{1/\beta}$ for some $\beta
% \in (0,\infty)$ for $U \sim Unif([0,1])$. In other words, the Dickman
% distributon is a Vervaat perpetutiy with $\beta = 1$.
%
% A generalization of the perpetuity takes the form \[ Y =
% \sum_{n=0}^{\infty}A_n \prod_{i}^{n} W_i\] with $A_n$ not necessarily
% equal to $1$, which is known as Takacs distribution. 
%
%
% See also: vervaatrnd, vervaatxdf, quickselectFS
%
%
%
% References:
%
% Cloud, K. and Huber, M. (2018), Fast Perfect Simulation of Vervaat
% Perpetuities, "Journal of Complexity", Vol. 42, pp. 19-30.
%
% Devroye, L. (2001), Simulating  perpetuities, "Methodology And Computing
% In Applied Probability", Vol. 3, Num. 1, pp. 97–115.
%
% Fill, J. A. and  Huber, M. (2010), Perfect  simulation  of  Vervaat
% perpetuities, "Electronic Journal of Probability", Vol. 15, pp. 96–109.
%
% Devroye, L. and Fawzi, O. (2010), Simulating the Dickman distribution,
% "Statistics and Probability Letters", Vol. 80, pp. 242–247.
%
% Blanchet, J. H. and Sigman, K. (2011), On exact sampling of stochastic
% perpetuities, "Journal of Applied Probability", Vol. 48A, pp. 165–182.
%
% Takacs, L. (1955), On stochastic processes connected with certain
% physical recording apparatuses. "Acta Mathematica Academiae
% Scientificarum Hungarica", Vol. 6, pp 363-379.
%
% Barabesi, L. and Pratelli, L. (2019), On the properties of a Takacs
% distribution, "Statistics and Probability Letters", Vol. 148, pp. 66-73.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('vervaatsim')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-02-01 15:26:55 #$: Date of the last commit
%
% Examples:
%

%{
    % Density of Vervaat perpetuity with $\beta=1$ and no optional parameters.
    clear all;
    close all;
    betav = 1;
    y=vervaatsim(betav)
%}

%{
    % Density of Vervaat perpetuity with $\beta=1$ and 10 chain steps.
    clear all;
    close all;
    betav = 1;
    steps = 10;
    y=vervaatsim(betav,steps)
%}

%{
    % Density of Vervaat perpetuity with $\beta=1$, 10 chain steps and initial chain value d=2.
    clear all;
    close all;
    betav = 1;
    steps = 10;
    d = 2;
    y=vervaatsim(betav,steps,d)
%}

%{
    %% N=5000 random values extracted from two Vervaat perpetuities.
    % Vervaat parameters ate: $\beta = 1$ and $\beta = 10$.
    % The superimposed normal kernel density is just for illustration.

    clear all; close all;
    betav10 = 10;
    betav01 = 1;
    N = 5000;

    for i=1:N
        y10(i) = vervaatsim(betav10);
        y01(i) = vervaatsim(betav01);
    end

    pdy10 = fitdist(y10(:),'Kernel','Kernel','normal','Support','positive');
    x10   = (1:N)*(5*betav10)/N;
    pdf10 = pdf(pdy10,x10);

    pdy01 = fitdist(y01(:),'Kernel','Kernel','normal','Support','positive');
    x01   = (1:N)*(5*betav01)/N;
    pdf01 = pdf(pdy01,x01);

    figure;
    h1=subplot(2,1,1);
    plot(x10,pdf10,'-','LineWidth',2); ylim([0,1]);
    if ~verLessThan('matlab','1.7.0')
        hold on
        h=histogram(y10);
        h.Normalization='pdf';
        h.BinWidth=0.02;
        h.EdgeColor='none'; 
        hold off
    end

    h2=subplot(2,1,2);
    plot(x01,pdf01,'-','LineWidth',2); ylim([0,1]);
    if ~verLessThan('matlab','1.7.0')
        hold on
        h=histogram(y01);
        h.Normalization='pdf';
        h.BinWidth=0.02;
        h.EdgeColor='none'; 
        hold off
    end

    title(h1,'$$\beta = 10$$' ,'Fontsize',20,'interpreter','latex');
    title(h2,'$$\beta = 1 $$ (Dickman)'  ,'Fontsize',20,'interpreter','latex');

%}

%{ 
    % A rough time test on the method.

    clear all;
    close all;
    betav = 1;
    n = 100;
    nmax    = 500;
    betamax = 10;
    rng(123);
    N     = randi(nmax,n,1);
    betav = randi(betamax,n,1);

    ti = tic;
    for i=1:n;
        for j=1:N(i)
            y1(j) = vervaatsim(betav(i));
        end
    end
    tf = toc(ti);

    disp(['Cloud-Huber: etime = ' num2str(tf)]);

%}

%{
    %% Density of Vervaat perpetuity: comparison with original R function.
    
    clear all;
    close all;

    % set common parameters
    betav = 1;
    steps = 10;
    d = 3;

    % ensure same random numbers in MATLAB and R 
    rn = mtR(1);    

    % Now compute the veervaat in MATLAB
    y=vervaatsim(betav,steps,d)

    % And finally compute the vervaat in R
    disp('To verify coherence of results with the R implementation,');
    disp('execute in R the following four lines,'); 
    disp('which allow to replicate the same random numbers;');
    disp('then execute the vervaat after sourcing the following R code: ') ;
    disp(' ');
    disp('NGkind("Mersenne-Twister") #set "Mersenne-Twister" "Inversion"');
    disp('set.seed(0) ');
    disp('state = .Random.seed ');
    disp('runif(1) ');
    disp('vervaat(beta = 1,steps = 10,d = 3) ');
    disp(' ');

%}



%% Beginning of code

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
    y = vervaatsim(betav, 2 * steps, d(1));
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

% (although some authors define it as $1 + Y$ for $Y$ being a Vervaat
% perpetutiy with $\beta = 1$).

%FScategory:UTISTAT
