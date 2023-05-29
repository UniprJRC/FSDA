function c = ASbdp(bdp,v)
%ASbdp finds the constant c associated to the supplied breakdown point for Andrew's sine function
% The constant is found through a dichotomic search
%
%
%<a href="matlab: docsearchFS('ASbdp')">Link to the help function</a>
%
%  Required input arguments:
%
%      bdp    : breakdown point. Scalar. Scalar defining breakdown point
%               (i.e a number between 0 and 0.5)
%               Data Types - single|double
%        v    : number of response variables. Scalar. e.g. in regression v=1
%               UP TO NOW v=1 (JUST REGRESSION) TO DO FOR MULTIVARIATE
%               ANALYSIS
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
% Output:
%
%  c : Requested tuning constant. Scalar. Tuning constatnt of Andrew's sine
%         function associated to requested breakdown point
%
%
% See also: OPTbdp, HYPbdp, HAbdp, PDbdp, TBbdp
%
% References:
%
% Andrews, D.F., Bickel, P.J., Hampel, F.R., Huber, P.J., Rogers, W.H., and
% Tukey, J.W. (1972), "Robust Estimates of Location: Survey and Advances",
% Princeton Univ. Press, Princeton, NJ. [p. 203]
%
% Andrews, D. F. (1974). A Robust Method for Multiple Linear Regression,
% "Technometrics", V. 16, pp. 523-531, https://doi.org/10.1080/00401706.1974.10489233
%
% Riani, M., Cerioli, A. and Torti, F. (2014), On consistency factors and
% efficiency of robust S-estimators, "TEST", Vol. 23, pp. 356-387,
% http://dx.doi.org/10.1007/s11749-014-0357-7
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('ASbdp')">Link to the help page for this function</a>
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
    % c=0.4495
    c=ASbdp(0.5,1)
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
%  where \rho(c) for AS is 2*c
Erho1=10;
eps=1e-11;

if nargin<=2 || v==1

    while abs(Erho1-1)>eps

        Erho =integral(@(u)(ASrho(u,c)).*normpdf(u),-c*pi,c*pi);
        c2=2*c;
        Erho1=(real(Erho)+2*c2*normcdf(c*pi,'upper'))/(c2*bdp);

        step=step/2;
        if Erho1>1
            c=c+step;
        else
            c=max(c-step,0.1);
        end
        % disp([step c Erho1])
    end
else
    error('FSDA:ASbdp:Wrongv','Not yet implemented for v>1')
end

end
%FScategory:UTISTAT