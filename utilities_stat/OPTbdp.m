function c = OPTbdp(bdp,v)
%OPTbdp finds the constant c associated to the supplied breakdown point
% The constant is found through a dichotomic search
%
%
%<a href="matlab: docsearchFS('OPTbdp')">Link to the help function</a>
%
%  Required input arguments:
%
%      bdp    : breakdown point. Scalar. Scalar defining breakdown point
%               (i.e a number between 0 and 0.5)
%               Data Types - single|double
%        v    : number of response variables. Scalar. e.g. in regression p=1
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
% Output:
%
%  c : Requested tuning constant. Scalar. Tuning constatnt of optimal rho
%         function associated to requested breakdown point
%
%
% More About:
%
%
%    $\rho$ ($\psi$) function which is considered is standardized 
%    using intervals 0---(2/3)c , (2/3)c---c, >c.   
%    $\rho$ function is
%
% \[
% TBrho(u)= \left\{
%    \begin{array}{lr}
%     1.3846 \left(\frac{u}{c}\right)^2                      &                                      |\frac{u}{c}| \leq  \frac{2}{3} \\
%    0.5514-2.6917 \left(\frac{u}{c}\right)^2 +10.7668\left(\frac{u}{c}\right)^4-11.6640\left(\frac{u}{c}\right)^6+4.0375\left(\frac{u}{c}\right)^8   & \qquad \frac{2}{3} \leq  |\frac{u}{c}|\leq  1 \\
%    1                                                    &                      |\frac{u}{c}|>1 \\
% \end{array}
%    \right.
%  \]
%
%   Therefore, to obtain the value of c for the (rho) psi function defined in the
%   interval 0---2c, 2c---3c, >3c it is necessary to divide the output of
%   function OPTbdp by 3.
%
%
% See also: TBbdp, HYPbdp, HAbdp
%
%
% References:
% 
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('OPTbdp')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

%
% Examples:
%
%{
    % Find c given bdp.
    % The constant c associated to a breakdown point of 50% in regression is 
    % c= 1.213897063441983
    c=OPTbdp(0.5,1)
%}
%

%% Beginning of code

% c = starting point of the iteration
c=5;
% step = width of the dichotomic search (it decreases by half at each
% iteration). Generally it can be smaller. A large value ensures converge
% when bdp is very small and p is very large.
step=200;

% Convergence condition is E(\rho) = \rho(c) bdp
% For simplicity below we put the standardized rho function so that
% \rho(c)=1
Erho1=10;
eps=1e-10;
while abs(Erho1-1)>eps
    
    b=c.^2/2;
    a=2*c.^2/9;
    Erho= (1.3846./c.^2)*v*gammainc(a,0.5*(v+2)) ...
          + 0.5514.*1*(gammainc(b,v/2) -gammainc(a,v/2) ) ...
          -(2.6917./c^2)*v*(gammainc(b,(v+2)/2) -gammainc(a,(v+2)/2) )+...
          +(10.7668./c^4)*v*(v+2)* (gammainc(b,(v+4)/2) -gammainc(a,(v+4)/2) )+...
          -(11.664./c^6)*v*(v+2)*(v+4)*(gammainc(b,(v+6)/2) -gammainc(a,(v+6)/2) )+...
          +(4.0375./c^8)*v*(v+2)*(v+4)*(v+6)*(gammainc(b,(v+8)/2) -gammainc(a,(v+8)/2) )+...
          +(1-gammainc(b,v/2));
    Erho1=(Erho/bdp);
    
    step=step/2;
    if Erho1>1
        c=c+step;
    else
        c=max(c-step,0.1);
    end
   %  disp([step c Erho1])
end
% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);


end
%FScategory:UTISTAT