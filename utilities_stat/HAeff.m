function ceff = HAeff(eff,v,abc)
%HAeff finds the tuning constant that guarrantees a requested asymptotic efficiency
%
%<a href="matlab: docsearchFS('HAeff')">Link to the help function</a>
%
%  Required input arguments:
%
%    eff:       efficiency. Scalar.  Scalar which contains the required
%               efficiency (of location or scale estimator).
%               Generally eff=0.85, 0.9 or 0.95.
%    v :        number of response variables. Scalar. Number of variables of
%               the  dataset (for regression v=1)
%               UP TO NOW v=1 (JUST REGRESSION) TO DO FOR MULTIVARIATE
%               ANALYSIS
%
%
%  Optional input arguments:
%
%     abc     : parameters of Hampel estimator. Vector. Vector of length 3
%               which contains the parameters of Hampel estimator. If
%               vector abc is not specified it is set equal to [2, 4, 8]
%               Example - [1.5,3.5,8]
%               Data Types - double
%
%
% Output:
%
%  ceff : Requested tuning constant. Scalar. Tuning constatnt of Hampel rho
%         function associated to requested value of efficiency
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
%	It is necessary to have 0 <= a <= b <= c.
%
% Parameter ctun multiplies parameters (a,b,c) of Hampel estimator.
%
%
% See also: TBeff, HYPeff, OPTeff, RKeff, HUeff
%
% References:
%
% Hoaglin, D.C., Mosteller, F., Tukey, J.W. (1982), "Understanding Robust and
% Exploratory Data Analysis", Wiley, New York.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HAeff')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
%
% Examples:
%
%{
    % Find c for fixed efficiency.
    % The constant c associated to a nominal location efficiency of 95% in regression is
    % c = 0.690998716841394
    c=HAeff(0.95,1)
%}

%{
    % Example where three input parameters are supplied.
    % Find constant c associated to a nominal location efficiency of 95 per
    % cent in regression when tun=[1.5,3.5,8].
     tun=[1.5,3,8];
    c=HAeff(0.95,1,tun);
%}
%
%


%% Beginning of code

if (nargin >2)
    if ((abc(1) < 0) || (abc(2) < abc(1)) || (abc(3) < abc(2)))
        error('FSDA:HAeff:WrongAbc','Illegal choice of parameters in Hampel: %f,%f,%f', abc(1),abc(2),abc(3));
        % error('FSDA:HAeff:WrongAbc',[' illegal choice of parameters in Hampel: ' ...
        %    num2str(abc) ]')
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



% LOCATION EFFICIENCY

% ctun = starting point of the iteration
ctun=0.57;
% c = starting point of the iteration
% Remark: the more refined approximation for the starting value of
% sqrt(chi2inv(eff,p))+2; does not seem to be necessary in the case of
% location efficiency

% step = width of the dichotomic search (it decreases by half at each
% iteration).
step=0.5;


% Convergence condition is
%  .......
empeff=1;

eps=1e-14;
% ctun=(0:0.01:1)';
while abs(empeff-eff)> eps
    
    a=(a0*ctun);
    b=(b0*ctun);
    c=(c0*ctun);
    
    a2=a.^2/2;
    b2=b.^2/2;
    c2=c.^2/2;
    
    % bet  = \int  \psi'(x) d \Phi(x)
    % bet = \int_-a^a d \Phi(x) +2* \int_b^c -a/(c-b)
    bet= gammainc(a2,v/2)+(gammainc(b2,v/2)-gammainc(c2,v/2))*a/(c-b);
    
    % alph = \int \psi^2(x) d \Phi(x)
    alph= v*gammainc(a2,(v+2)/2)...                                        % 2* \int_0^a x^2 f(x) dx
        +a.^2 .*(gammainc(b2,v/2)-gammainc(a2,v/2))...                     % 2* a^2 \int_a^b f(x) dx
        +(a./(c-b)).^2 .*(c.^2.*(gammainc(c2,v/2)-gammainc(b2,v/2)) ...    %(a./(c-b)).^2 (2 c^2 \int_b^c f(x) dx
        + v*(gammainc(c2,(v+2)/2)-gammainc(b2,(v+2)/2)) ...                %   + 2*  \int_b^c x^2 f(x) dx
        -2*c*v*sqrt(2/pi)*(gammainc(c2,(v+1)/2)-gammainc(b2,(v+1)/2)));        % +2 *2* \int_b^c |x| f(x)
    
    
    % Remark: if v=1
    % -2*c*v*sqrt(2/pi)*(gammainc(c2,(v+1)/2)-gammainc(b2,(v+1)/2)));
    %     -4*c.*(normpdf(b)-normpdf(c))  );
    
    % empeff = bet^2/alph = 1 / [var (robust estimator of location)]
    empeff=(real(bet)^2)/real(alph);
    
    step=step/2;
    if empeff<eff
        ctun=ctun+step;
    elseif empeff>eff
        ctun=max(ctun-step,0.01);
    else
    end
    
    % disp([empeff eff ctun])
end


ceff=ctun;
% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);

end
%FScategory:UTISTAT
