function pdf = twdpdf(x,alpha,theta,delta)
%twopdf computes the probability density function of the Tweedie distribution.
%
%<a href="matlab: docsearchFS('twdpdf')">Link to the help function</a>
%
% This function returns the pdf of a Tweedie distribution evaluated on an
% array of values $x$. The description of the Tweedie distribution and its
% parameters is detailed in fuction twdrnd.
%
% Required input arguments:
%
% x     : Values at which to evaluate pdf. Array of scalar values.
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
%   pdf : Probability density function values. Scalar value or array of
%         scalar values. pdf values, evaluated at the values in x,
%         returned as a scalar value or an array of scalar values.
%         Data Types - Double.
%
%
% More About: 
%
% Detailed information can be found in function twdrnd.
%
% See also: twdrnd
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
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('twdpdf')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2019-02-01 15:26:55 #$: Date of the last commit
%
% Examples:
%
%{
        % Example from Barabesi et. al (2016).
        n = 1000;
        pMushrooms = [-0.0936 , 1.27 *10^(-6)  , 0.6145];
        param1 = pMushrooms;
        tit1 = 'Mushrooms';
        al = param1(1) ; th = param1(2) ; de = param1(3) ;
        x = twdrnd(al,th,de,n);
        pdf = twdpdf(x,al,th,de);
        figure;
        plot(x,pdf,'r.');
        title(tit1);
%}

%{
        %% Check consistency with R theedie package.
        % Note that the package adopts the parametrization of Jorgensen
        % (1987), which introduces the Tweedie distribution as a special
        % case of the exponential dispersion model.
        % In the R tweedie package, the parameters are:
        % p     = vector of probabilities;
        % n     = the number of observations;
        % xi    = the value of xi such that the variance is $var(Y)=\phi \mu^{xi}$;
        % power = a synonym for xi
        % mu    = the mean
        % phi	= the dispersion
        
        % # R code
        % power <- 2.5
        % mu    <- 1
        % phi   <- 1
        % y     <- seq(0, 6, length=500)
        % fy    <- dtweedie( y=y, power=power, mu=mu, phi=phi)
        % plot(y, fy, type="l", lwd=2, ylab="Density")
        % # Compare to the saddlepoint density
        % f.saddle <- dtweedie.saddle( y=y, power=power, mu=mu, phi=phi)
        % lines( y, f.saddle, col=2 )
        % legend("topright", col=c(1,2), lwd=c(2,1),
        %     legend=c("Actual","Saddlepoint") )
        
        % parameter values in Jorgensen parametrization
        power = 2.5;
        mu    = 1 ;
        phi   = 1 ;

        % reparametrization
        alpha = (power-2)/(power-1);
        theta = phi*mu;
        delta = (phi/(power-1))^(1/(power-1));

        % pdf
        y = linspace(0,6,500);
        pdf = twdpdf(y,alpha,theta,delta);
        figure;
        subplot(2,1,1);
        plot(y,pdf,'r.');
        title({'Jorgensen parametrization: power = 2.5, mu = 1, phi = 1' , ...
               're-parametrization: alpha = 0.3333, theta = 1, delta = 0.7631'});
        
        Y = twdrnd(alpha,theta,delta,500);
        subplot(2,1,2);
        histogram(Y);

%}

%% Beginning of code

% ensure that x is in column vector
x = x(:);
        
pdf   = zeros(size(x));
x0    = (x==0);
xnot0 = (x>0);
c     = (delta*theta^alpha)/alpha;
if alpha>0
    pdf(x0) = 0;      % When alpha>0, Prob(x=0)=0
else
    pdf(x0) = exp(c); % When alpha<0, this is the Prob(x=0)
end
pdf(xnot0)  = exp(-theta*x(xnot0)+c) .* InfSum(x(xnot0),alpha,delta);


%% Subfunction that computes the series that approximate the pdf

    function s = InfSum(x,al,de)
        
        N   = 100;   % number of addends
        seq = (1:N);
        
        m = length(x);
        
        log_fact_n  = cumsum(log(seq));  % vector 1xN
        exp_x       = -seq*al-1;         % vector 1xN (exponent of x)
        exp_x_log_x = log(x)*exp_x;      % matrix mxN
        m_ones      = ones(m,1);         % vector mx1 of 1
        if al>0
            n_log_de = seq*log(de);                                                   % vector 1xN
            n_log_al = seq*log(al);                                                   % vector 1xN
            M1       = exp(m_ones*(n_log_de - log_fact_n - n_log_al) + exp_x_log_x);  % matrix mxN
            minus1n  = (-1) .^ seq;                                                   % vector 1xN
            M2       = m_ones * (minus1n./gamma(-seq*al));                            % matrix mxN
            M        = M1 .* M2;
        else
            log_gamma_nal = gammaln(-seq*al);                                               % vector 1xN
            n_log_de_al   = seq*log(-de/al);                                                % vector 1xN
            M             = exp(m_ones*(n_log_de_al-log_fact_n-log_gamma_nal)+exp_x_log_x); % matrix mxN
        end
        s = max(sum(M,2),realmin);
    end

end
%FScategory:UTISTAT
