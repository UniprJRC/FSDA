function ctun = HAbdp(bdp,p,abc)
%HAbdp finds the constant c associated to the supplied breakdown point
% The constant is found through a dichotomic search
%
%
%<a href="matlab: docsearchFS('HAbdp')">Link to the help function</a>
%
%  Required input arguments:
%
%      bdp    : breakdown point. Scalar. Scalar defining breakdown point
%               (i.e a number between 0 and 0.5)
%               Data Types - single|double
%        p    : number of response variables. Scalar. e.g. in regression p=1
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
%     abc     : parameters of Hampel estimator. Vector. Vector of length 3
%               which contains the parameters of Hampel estimator. If
%               vector abc is not specified it is set equal to [2, 4, 8]
%               Example - [1.5,3.5,8]
%               Data Types - double
%
% Output:
%
%  ctun : Requested tuning constant. Scalar. Tuning constatnt of Hampel rho
%         function associated to requested breakdown point
%
%
% More About:
%
% Function HApsi transforms vector u as follows.
%  \[
%  HApsi(u)  = \left\{   
%  \begin{array}{cc}
%    u & |u| <= a                                       \\
%    a \times sign(u) & a <= |u| < b                    \\
%    a \frac{c-|u|}{c-b} \times sign(u) & b <= |u| <  c \\
%    0 & |u| >= c 
%  \end{array} \right.
% \]
%
%             where $a$= ctun *param(1).
%                   $b$= ctun *param(2).
%                   $c$= ctun *param(3).
%
%             The default is
%                   $a$= 2*ctun. 
%                   $b$= 4*ctun. 
%                   $c$= 8*ctun. 
%
%	It is necessary to have 0 <= a <= b <= c
%
% See also: TBbdp, HYPbdp, OPTbdp
%
% References:
%
% D. C. Hoaglin, F. Mosteller, J. W. Tukey (1982), Understanding Robust and
% Exploratory Data Analysis Wiley, New York.

% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HAbdp')">Link to the help page for this function</a>
% Last modified 31-05-2016
%
%
%
% Examples:
%
%{
    %% Find constant c for bdp=0.5.
    % The constant c associated to a breakdown point of 50 per cent in regression is
    % 0.198131771596856
    c=HAbdp(0.5,1);
    disp(c);
%}

%{
    %% Find constant c for bdp=0.5 when abc=[1.5 3.5 8].
    % The constant c associated to a breakdown point of 50 per cent in regression is
    c=HAbdp(0.5,1,[1.5 3.5 8]);
    disp(c);
%}



%% Beginning of code


if (nargin >2),
    if ((abc(1) < 0) || (abc(2) < abc(1)) || (abc(3) < abc(2))),
        error('FSDA:HAbdp:WrongAbc',[' illegal choice of parameters in Hampel: ' ...
            num2str(abc) ]')
        
    end
    a0 = abc(1);
    b0 = abc(2);
    c0 = abc(3);
else
    a0 = 2;
    b0 = 4;
    c0 = 8;
    %     a0 = 1.5;
    %     b0 = 3.5;
    %     c0 = 8;
end

% step = width of the dichotomic search (it decreases by half at each
% iteration). Generally it can be smaller. A large value ensures converge
% when bdp is very small and p is very large.
step=200;

% ctun = starting point of the iteration
ctun=0.4;

% Convergence condition is E(\rho) = \rho(\infty) bdp
Erho1=10;
eps=1e-10;
while abs(Erho1-1)>eps
    
    
    a=(a0*ctun);
    b=(b0*ctun);
    c=(c0*ctun);
    
    a2=a.^2/2;
    b2=b.^2/2;
    c2=c.^2/2;
    
    
    phic=a*b-0.5*a^2+0.5*(c-b)*a;
    
    % |u| <a
    % Erhoa=  \int_-a^a u^2/2
    Erhoa=0.5*p*gammainc(a2,(p+2)/2);
    % rhoa = @(u,a,b,c)(0.5*u.^2.*(1/sqrt(2*pi)).*exp(-0.5*u.^2));
    % Erhoack=integral(@(u)rhoa(u,a,b,c),-a,a);
    
    % a< |u| <b
    % Erhoab= 2 * \int_a^b (au-0.5a^2)
    Erhoab=2*a*(normpdf(a)-normpdf(b))-a^2*(normcdf(b)-normcdf(a));
    % rhoab = @(u,a,b,c)((a*u-0.5*a^2).*(1/sqrt(2*pi)).*exp(-0.5*u.^2));
    % Erhoabck=2*integral(@(u)rhoab(u,a,b,c),a,b);
    
    % b< |u| <c
    % Erhobc = \int_b^c \rho(x) \Phi(x)
    Erhobc1=2*(a*b-0.5*a^2+0.5*(c-b)*a*(1 -c^2/((c-b)^2)))*(normcdf(c)-normcdf(b));
    Erhobc2=0.5*a*p*(gammainc(b2,(p+2)/2) -gammainc(c2,(p+2)/2)) /(c-b);
    Erhobc3=2*a*c*(normpdf(b)-normpdf(c))/(c-b);
    Erhobc=Erhobc1+Erhobc2+Erhobc3;
    
    %     psi2 = @(u,a,b,c) (a*b-0.5*a^2+0.5*(c-b)*a*(1 -(c-u).^2/((c-b)^2))) .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
    %     Erhobcck=2*integral(@(u)psi2(u,a,b,c),b,c);
    %
    %     psi2 = @(u,a,b,c) (a*b-0.5*a^2+0.5*(c-b)*a*(1 -(c.^2)/((c-b)^2))) .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
    %     Erhobc1ck=2*integral(@(u)psi2(u,a,b,c),b,c);
    %
    %     psi2 = @(u,a,b,c) (0.5*(c-b)*a*( -u.^2/((c-b)^2))) .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
    %     Erhobc2ck=2*integral(@(u)psi2(u,a,b,c),b,c);
    %
    %     psi2 = @(u,a,b,c) (0.5*(c-b)*a*( 2*c*u/((c-b)^2))) .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
    %     Erhobc3ck=2*integral(@(u)psi2(u,a,b,c),b,c);
    
    
    % |u| >c
    Erhoc=phic*( 1-gammainc(c2,p/2) );
    
    
    Erho= Erhoa+Erhoab+Erhobc+Erhoc;
    
    Erho1=Erho/(phic*bdp);
    
    step=step/2;
    if Erho1>1
        ctun=ctun+step;
    else
        ctun=max(ctun-step,0.1);
    end
    % disp([step ctun Erho1])
end
% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);


end
%FScategory:UTISTAT