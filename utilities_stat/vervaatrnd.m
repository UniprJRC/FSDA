function X = vervaatrnd(betav,n,method)
%vervaatrnd simulates random variates from the Vervaat perpetuity distribution.
%
%<a href="matlab: docsearchFS('vervaatrnd')">Link to the help function</a>
%
% This function generates random numbers from a Vervaat perpetuity
% distribution with parameter $\beta \in (0,+\infty)$. Two alternative
% methods are implemented. The first method calls FSDA function vervaatsim,
% which is based on an accurate and efficient recursive simulation
% algorithm by Cloud and Huber (2018). The second method calls FSDA
% function vervaatxdf, which uses an algorithm by Barabesi and Pratelli
% (2019) to generate the pdf and cdf of the distribution. Below we give
% briefly some definitions, background, motivations and some key
% references.
%
%
% Required input arguments:
%
% betav : Distribution parameter value. Positive integer. The parameter of
%         the Vervaat family. Default is betav = 1.
%
% Optional input arguments:
%
%
% n     : Number of random numbers to extract. Positive integer. The sample
%         size. Default is n = 1.
%         Example - n=100.
%         Data Types - Scalar.
%
% method: Computation method. Integer in {1,2} or structure. The method
%         used to generate the random numbers: method = 1 (default) is for
%         Barabesi and Pratelli (2019); method = 2 is for Cloud and
%         Huber (2018). Method can be also a structure to pass the optional
%         parameters for method 2, which are:
%         - method.steps: Markov chain step. Positive integer. The inital
%                  number of steps to run to move the chain forward.
%                  Default is method.steps = 1.
%         - method.d    : Chain value at time 0. Positive integer or -1.
%                  The value  of the dominating chain at time 0. Default is
%                  $D_0 \leftarrow x_0−1+G$, where $G \sim \mbox{Geo}(1/2)$
%                  and $x_0 = \frac{1 + (2/3)^{1/\beta}}{1 -
%                  (2/3)^{1/\beta}}$ ($\mbox{Geo}$ represents the geometric
%                  distribution).
%         Example - method = 2;
%         Example - method = struct; method.steps=10; method.d=3;
%         Data Types - Struct.
%
% Output:
%
% X     : Random values extracted from the Vervaat perpetuity. Scalar
%         or array. The random value(s). Data Types - Double.
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
% numerical array, approaches asymptotically the Dickman distribution,
% which is a perpetuity with $W \sim Unif([0,1])$. Unfortunately there is
% no known closed form for the Dickman.
%
% The Dickman distribution can be also seen as a special case of Vervaat
% perpetuity, which is such that $W_i \sim U^{1/\beta}$ for some $\beta
% \in (0,\infty)$ for $U \sim Unif([0,1])$. In other words, the Dickman
% distributon is a Vervaat perpetutiy with $\beta = 1$.
%
%
% See also: vervaatsim, vervaatxdf, quickselectFS
%
%
% References:
%
% Cloud, K. and Huber, M. (2018), Fast Perfect Simulation of Vervaat
% Perpetuities, "Journal of Complexity", Vol. 42, pp. 19-30.
%
% Devroye, L. (2001), Simulating  perpetuities, "Methodology And Computing
% In Applied Probability", Vol. 3, Num. 1, pp. 97-115.
%
% Fill, J. A. and  Huber, M. (2010), Perfect  simulation  of  Vervaat
% perpetuities, "Electronic Journal of Probability", Vol. 15, pp. 96-109.
%
% Devroye, L. and Fawzi, O. (2010), Simulating the Dickman distribution,
% "Statistics and Probability Letters", Vol. 80, pp. 242-247.
%
% Blanchet, J. H. and Sigman, K. (2011), On exact sampling of stochastic
% perpetuities, "Journal of Applied Probability", Vol. 48A, pp. 165-182.
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
%<a href="matlab: docsearchFS('vervaatrnd')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-02-01 15:26:55 #$: Date of the last commit
%
% Examples:
%

%{
    % A random value from a Vervaat perpetuity with a default parameter.
    % Argument betav is set to 1 (Dickman) if not specified.
    clear all;
    close all;
    y=vervaatrnd
%}

%{
    % n random values from a Vervaat perpetuity with a specified betav parameter.
    clear all;
    close all;
    n = 5;
    betav = 2;
    y=vervaatrnd(betav,5);
%}

%{ 
    %% n Vervaat perpetuity values extracterd with the two implemented methods.
    % In the example, we set $n = 5000$ and $\beta = 1$.
    % The superimposed normal kernel density is just for illustration: 
    % a more precise density can be simulated with function vervaatxdf.

    clear all;
    close all
    betav = 1;
    N = 5000;

    y1 = vervaatrnd(betav,N,1);
    y2 = vervaatrnd(betav,N,2);

    x    = (1:N)*(5*betav)/N;
    pd1  = fitdist(y1(:),'Kernel','Kernel','normal','Support','positive');
    pdf1 = pdf(pd1,x);

    pd2  = fitdist(y2(:),'Kernel','Kernel','normal','Support','positive');
    pdf2 = pdf(pd2,x);

    figure;
    h1 = subplot(2,1,1);
    p1h = plot(x,pdf1,'.','LineWidth',2); ylim([0,1]);
    if ~verLessThan('matlab','1.7.0')
        hold on
        h=histogram(y1);
        h.Normalization='pdf';
        h.BinWidth=0.02;
        h.EdgeColor='none'; 
        hold off
    end

    h2 = subplot(2,1,2);
    plot(x,pdf2,'.','LineWidth',2); ylim([0,1]);
    if ~verLessThan('matlab','1.7.0')
        hold on
        h=histogram(y2);
        h.Normalization='pdf';
        h.BinWidth=0.02;
        h.EdgeColor='none'; 
        hold off
    end

    title(h1,'Dickman ($$\beta = 1 $$) values, from Barabesi-Pratelli' ,'Fontsize',20,'interpreter','latex');
    title(h2,'Dickman ($$\beta = 1 $$) values, from Cloud-Huber'       ,'Fontsize',20,'interpreter','latex');

%}

%{ 
    % A time test for the two methods used to extract the Vervaat values.
    % Results can vary considerably with rep n and N (in this simulation
    % method 1 computes the pascal matrix at each replicate).

    clear all;
    close all;
    betav = 1;
    rep = 10;
    n = 10000;
    N = randi(n,rep,1);

    t1 = tic;
    for i=1:rep;
        y1 = vervaatrnd(betav,N(i),1);
    end
    t1end = toc(t1);

    t2 = tic;
    for i=1:rep;
        y2 = vervaatrnd(betav,N(i),2);
    end
    t2end = toc(t2);

    disp(['Barabesi-Pratelli: etime = ' num2str(t1end)]);
    disp(['Cloud-Huber: etime = ' num2str(t2end)]);

%}

%{
    %% N=5000 random values extracted from two Vervaat perpetuities.
    % Parameters are: $\beta = 1$ and $\beta = 10$.
    % The superimposed normal kernel density is for illustration: 
    % a more precise density can be simulated with function vervaatxdf.

    betav10 = 10;
    betav01 = 1;
    N = 5000;

    y10 = vervaatrnd(betav10,N);
    y01 = vervaatrnd(betav01,N);

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




if nargin<3 || isempty(method)
    method=1;
end
if nargin==3 && isstruct(method)
    steps   = method.steps;
    d       = method.d;
    method  = 2;
end
if nargin==3 && method==2
    steps   = 1;
    d       = -1;    
end
if nargin<2 || isempty(n)
    n=1;
end
if nargin==0 || isempty(betav) || betav < 1
    betav = 1;
end

switch method

    case 1
        % F is the cdf of the vervaat; x is the set of evaluation points;
        % the pdf is not used (~).
        [~ , F , x] = vervaatxdf(betav,500);
        
        % remove inf and nan
        maskinf = not(isinf(F));
        masknan = not(isnan(F));
        mask = and(maskinf,masknan);
        x = x(mask); F = F(mask);
        
        % remove non-unique elements
        [F, mask] = unique(F);
        x = x(mask);

        % extract n random values from the uniform
        rvu = rand(1, n);

        % inverse interpolation to get the projection of the random values
        X = interp1(F, x, rvu,'linear','extrap');  % 'spline'

    case 2
        X  = zeros(1,n);
        for k = 1:n
            X(k)= vervaatsim(betav,steps,d);
        end
end

end

%FScategory:UTISTAT

