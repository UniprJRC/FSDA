function c = TBbdp(bdp,v)
%TBbdp finds the constant c associated to the supplied breakdown point for Tukey's biweight
% The constant is found through a dichotomic search
%
%
%<a href="matlab: docsearchFS('TBbdp')">Link to the help function</a>
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
%  c : Requested tuning constant. Scalar. Tuning constatnt of Tukey Biweight
%         function associated to requested breakdown point
%
%
% See also: OPTbdp, HYPbdp, HAbdp
%
% References:
% 
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
%
% Copyright 2008-2017.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('TBbdp')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
%
% Examples:
%
%{
    % Find c given bdp.
    % The constant c associated to a breakdown point of 50% in regression is 
    % c=1.547644980928226
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
while abs(Erho1-1)>eps
    
    c2=c.^2/2;
    Erho= (v*gammainc(c2,0.5*(v+2))/2-(v^2+2*v)*gammainc(c2,0.5*(v+4))./(4*c2)+...
        +(v^3+6*v^2+8*v)*gammainc(c2,0.5*(v+6))./(6*(c.^4))+ ((c.^2)/6).*(1-gammainc(c2,v/2))  );
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
%FScategory:UTISTAT