function c = OPTbdp(bdp,p)
%OPTbdp finds the constant c associated to the supplied breakdown point
% The constant is found through a dichotomic search
%
%
%<a href="matlab: docsearchFS('OPTbdp')">Link to the help function</a>
%
%  Required input arguments:
%
%      bdp    : scalar defining breakdown point (i.e a number between 0 and 0.5)
%        p    : number of response variables (e.g. in regression p=1) 
%
% Output:
%
%  c = scalar associated to that particular breakdown point
%
%    REMARK: \rho (\psi) function which is considered is standardized 
%    using intervals 0---(2/3)c , (2/3)c---c, >c   
%    Rho function is
%
%               |   1.3846 |r/c|^2                                                            |t/c|<=2/3
%               |   
%     \rho(r) = |   0.5514-2.6917|r/c|^2+10.7668|r/c|^4-11.6640|r/c|^6+4.0375|r/c|^8     2/3<=|t/c|<=1
%               |
%               |   1                                                                          |t/c|>1                              
%
%   Therefore, to obtain the value of c for the (rho) psi function defined in the
%   interval 0---2c, 2c---3c, >3c it is necessary to divide the output of
%   function OPTbdp by 3.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('OPTbdp')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%

%
% Examples:
%
%{
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
    Erho= (1.3846./c.^2)*p*gammainc(a,0.5*(p+2)) ...
          + 0.5514.*1*(gammainc(b,p/2) -gammainc(a,p/2) ) ...
          -(2.6917./c^2)*p*(gammainc(b,(p+2)/2) -gammainc(a,(p+2)/2) )+...
          +(10.7668./c^4)*p*(p+2)* (gammainc(b,(p+4)/2) -gammainc(a,(p+4)/2) )+...
          -(11.664./c^6)*p*(p+2)*(p+4)*(gammainc(b,(p+6)/2) -gammainc(a,(p+6)/2) )+...
          +(4.0375./c^8)*p*(p+2)*(p+4)*(p+6)*(gammainc(b,(p+8)/2) -gammainc(a,(p+8)/2) )+...
          +(1-gammainc(b,p/2));
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
