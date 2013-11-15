function c = TBbdp(bdp,p)
%TBbdp finds the constant c associated to the supplied breakdown point
% The constant is found through a dichotomic search
%
%
%<a href="matlab: docsearch('tbbdp')">Link to the help function</a>
%
%  Required input arguments:
%
%      bdp    : scalar defining breakdown point (i.e a number between 0 and 0.5)
%        p    : number of response variables (e.g. in regression p=1) 
%
% Output:
%
%  c = scalar of Tukey Biweight associated to that particular breakdown point
%
% Copyright 2008-2013.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti
%
%
%<a href="matlab: docsearch('tbbdp')">Link to the help page for this function</a>
% Last modified 02-May-2013
%
%
%
% Examples:
%
%{
    % The constant c associated to a breakdown point of 50% in regression is 
    % c=1.5476449809
    c=TBbdp(0.5,1)
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
%  where \rho(c) for TBW is c^2/6
Erho1=10;
eps=1e-11;
p1=(p^2+2*p);       % p1 = p(p+2)
p2=(p^3+6*p^2+8*p); % p2 =p(p+2)(p+4)

p2h=0.5*(p+2);
p4h=0.5*(p+4);
p6h=0.5*(p+6);

while abs(Erho1-1)>eps
    
    c2=c.^2/2;
    Erho= (p*gammainc(c2,p2h)/2-p1*gammainc(c2,p4h)./(4*c2)+...
        +p2*gammainc(c2,p6h)./(6*(c.^4))+ ((c.^2)/6).*(1-gammainc(c2,p/2))  );
    Erho1=(Erho./(c.^2))*(6/bdp);
    
    step=step/2;
    if Erho1>1
        c=c+step;
    else
        c=max(c-step,0.1);
    end
    % disp([step c Erho1])
end
% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);
end
