function cdf = twdcdf(x,alpha,theta,delta)
% TWDCDF computes the cumulative distribution function of the Tweedie distribution.
%
%<a href="matlab: docsearchFS('twdcdf')">Link to the help function</a>
%
% This function returns the cdf of a Tweedie distribution evaluated on an
% array of values $x$. The description of the Tweedie distribution and its
% parameters is detailed in fuction twdrnd.
%
% Required input arguments:
%
% x     : Values at which to evaluate cdf. Array of scalar values.
%         To evaluate the pdf at multiple values, specify x using an array.
%
% alpha : Distribution parameter value. Non-zero value smaller than 1. The
%         location parameter of the Tweedie distribution. Default is alpha = 1.
%
% theta : Distribution parameter value. Positive value. The dispersion
%         parameter of the Tweedie distribution. Default is theta = 0
%         (Dirac).
%
% delta : Distribution parameter value. Positive value. The power parameter
%         of the Tweedie distribution. delta is such that the variance is
%         $var(Y) = \theta * \alpha^\delta$. delta is greater than or equal
%         to one, or less than or equal to zero. Default is delta = 1.
%         Interesting special cases are: the normal (delta=0), Poisson
%         (delta=1 with theta=1), gamma (delta=2) and inverse Gaussian
%         (delta=3). Other values of delta lead to cases that cannot be
%         written in closed form, and the computation becomes difficult.
%         When 1 < delta < 2, the distribution is continuous for Y>0 and
%         has a positive mass at Y=0. For delta > 2, the distribution
%         is continuous for Y>0.
%
%
% Optional input arguments:
%
%
% Output:
%
%   cdf : Cumulative distribution function values. Scalar value or array of
%         scalar values. pdf values, evaluated at the values in x,
%         returned as a scalar value or an array of scalar values.
%         Data Types - Double.
%
%
% More About: 
%
% Detailed information can be found in function twdrnd.
%
% See also: twdrnd, twdpdf
%
% References:
%
% Tweedie, M. C. K. (1984), An index which distinguishes between some important
% exponential families, "in Statistics: Applications and New Directions,
% Proceedings of the Indian Statistical Institute Golden Jubilee
% International Conference (J.K. Ghosh and J. Roy, eds.), Indian Statistical
% Institute, Calcutta", pp. 579-604.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('twdcdf')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-02-01 15:26:55 #$: Date of the last commit
%
% Examples:
%
%{
        % Estimating the empirical CDF passing from the integral of the pdf. 
        % Data taken from example by Barabesi et. al (2016).
        clear all
        close all
        n = 500;
        pMushrooms = [-0.0936 , 1.27 *10^(-6)  , 0.6145];
        param1 = pMushrooms;
        tit1 = 'Mushrooms';
        al = param1(1) ; th = param1(2) ; de = param1(3) ;
        x = twdrnd(al,th,de,n);

        pdf = twdpdf(x,al,th,de);
        figure;
        subplot(2,1,1);
        plot(x,pdf,'r.');
        title([tit1 ' - Tweedie PDF'],'interpreter','latex','Fontsize',16);

        cdf = twdcdf(x,al,th,de);
        subplot(2,1,2);
        plot(x,cdf,'b.');
        title([tit1 ' - Tweedie CDF, found by integrating the PDF'],'interpreter','latex','Fontsize',16);
        
%}

%{
    % Estimating the empirical CDF by its definition. 
    % F(x) = (Number of observations in X<=x)/(Total number of observations).
    % x = vector specifying random variables.
    % y = vector specifying points for which CDF and PDF has to be evaluated.  

    % The CDF is evaluated using the data in the vector x above, 
    % at the points in y.

    close all
    clear all

    n = 500;
    param1 = [-0.0936 , 1.27 *10^(-6)  , 0.6145];
    al = param1(1) ; th = param1(2) ; de = param1(3) ;
    x = twdrnd(al,th,de,n);
    y  = linspace(0,max(x));

    ny = length(y);
    nx = length(x);

    for i = 1:ny
        p = 0;              % True  Probability
        q = 0;              % False Probability
        for j = 1:nx
            if x(j)<=y(i)   % Definition of CDF
                p = p + 1;
            else
                q = q + 1;
            end
        end
        F(i) = p/(p + q);   % Calulating Probability
    end
    plot(y,F,'o');

    % Now a faster and more compact version of the previous definition

    x1  = unique(x); % in principle duplicates should be removed
    nx = length(x1);
    y1  = linspace(0,max(x1));
    ny = length(y1);
    sump = sum(x1<=y1 , 1);
    FF = sump/nx;
    hold on;
    plot(y1,FF,'-r')

    legend('CDF using loops' , 'CDF using vectorization');

    title({'CDF computed using definition' , '$F(x) = (\sum (I(X<=x)))/(n)$'},'interpreter','latex','Fontsize',16);

    % both definitions above are along the matlab fiunction cdfplot(x)

%}

%% Beginning of code

f= @(x) twdpdf(x,alpha,theta,delta);
cdf = zeros(size(x));
for j=1:length(x)
   cdf(j) = integral(f,0,x(j),'ArrayValued',true);
end

end
%FScategory:UTISTAT
