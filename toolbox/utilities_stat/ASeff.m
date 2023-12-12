function ceff = ASeff(eff, v)
%ASeff finds the constant c which is associated to the requested efficiency for Andrew's sine function
%
%
%
%<a href="matlab: docsearchFS('ASeff')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    eff:       required efficiency. Scalar.
%               Scalar which contains the required efficiency (of location
%               or scale estimator).
%               Data Types - single|double
%               Generally eff=0.85, 0.9 or 0.95
%    v :        Number of response variables. Scalar. e.g. in regression p=1
%               UP TO NOW v=1 (JUST REGRESSION) TO DO FOR MULTIVARIATE
%               ANALYSIS
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
%
% Output:
%
%  ceff : Requested tuning constant. Scalar. Tuning constatnt of Andrew's rho
%         function associated to requested value of efficiency
%
% See also: OPTeff, HYPeff, HAeff, TBeff, PDeff
%
% References:
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
%<a href="matlab: docsearchFS('ASeff')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

%
% Examples:
%
%{
    % Find c for a given efficiency.
    % The constant c associated to a nominal location efficiency of 95% in regression is
    % c = 1.3387
    c=ASeff(0.95,1)
%}


%% Beginning of code

eps=1e-12;
if nargin<=2 || v==1
    % LOCATION EFFICIENCY
    c= 2;
    % c = starting point of the iteration
    % Remark: the more refined approximation for the starting value of
    % sqrt(chi2inv(eff,p))+2; does not seem to be necessary in the case of
    % location efficiency
    
    % step = width of the dichotomic search (it decreases by half at each
    % iteration).
    step=30;
    
    
    % Convergence condition is 
    %  .......
    empeff=10;
     
    
    % bet  = \int_{-c}^c  \psi'(x) d \Phi(x)
    % alph = \int_{-c}^c  \psi^2(x) d \Phi(x)
    % empeff = bet^2/alph = 1 / [var (robust estimator of location)]
    while abs(empeff-eff)> eps
        
       
        bet= integral(@(u)(ASpsider(u,c)).*normpdf(u),-c*pi,c*pi);
        
        alph= integral(@(u)((ASpsi(u,c)).^2).*normpdf(u),-c*pi,c*pi);
        empeff=(real(bet)^2)/real(alph);
        
        step=step/2;
        if empeff<eff
            c=c+step;
        elseif empeff>eff
            c=max(c-step,0.1);
        else
        end
        
    end
else
    error('FSDA:ASeff:Wrongv','Not yet implemented for v>1')
end

ceff=c;

end
%FScategory:UTISTAT
