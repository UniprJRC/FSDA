function ceff = HUeff(eff,v,shapeeff)
%HUeff finds the constant c which is associated to the requested efficiency for Tukey biweight estimator
%
%
%
%<a href="matlab: docsearchFS('HUeff')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    eff:       required efficienty. Scalar.
%               Scalar which contains the required efficiency (of location
%               or scale estimator).
%               Data Types - single|double
%               Generally eff=0.85, 0.9 or 0.95
%    v :        Number of response variables. Scalar. e.g. in regression v=1
%               Now it is implemented just for v=1
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
%   shapeeff : Location or shape efficiency. Scalar.
%              If 1, the efficiency is referred to shape else (default)
%              is referred to location (not implemented yet)
%               Example - 'shapeeff',1
%               Data Types - double
%
% Output:
%
%  ceff : Requested tuning constant. Scalar. Tuning constatnt of Tukey Biweigh rho
%         function associated to requested value of efficiency
%
% See also: OPTeff, HYPeff, HAeff
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
%<a href="matlab: docsearchFS('HUeff')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

%
% Examples:
%
%{
    % Find c in regression for 95 per cent efficiency.
    % The constant c associated to a nominal location efficiency of 95% in regression is
    % c = 1.344997508513144
    c=HUeff(0.95,1)
%}
%
%
%
%{
    % Analyze constant c as a function of eff.
    % Initialize the matrix which stores the values of c for the two
    % methods
    eff=[0.70:0.0001:0.9999];
    cc=[eff' zeros(length(eff),1)];

    for i=1:length(eff)
        % Use exact formula for finding the value of c associated to a fixed
        % level of shape efficiency
        cc(i,2)=TBeff(eff(i),1);
    end
    figure
    plot(cc(:,1),cc(:,2),'LineStyle','-','LineWidth',2)
    xlabel('Effciency')
    ylabel('Value of c')
%}


%% Beginning of code

eps=1e-12;


    % LOCATION EFFICIENCY
    c= 2;
    % c = starting point of the iteration
    % Remark: the more refined approximation for the starting value of
    % sqrt(chi2inv(eff,p))+2; does not seem to be necessary in the case of
    % location efficiency
    
    % step = width of the dichotomic search (it decreases by half at each
    % iteration).
    step=10;
    
    
    % Convergence condition is 
    %  .......
    empeff=20;
    
    
    % bet  = \int_{-c}^c  \psi'(x) d \Phi(x)
    % alph = \int_{-c}^c  \psi^2(x) d \Phi(x)
    % empeff = bet^2/alph = 1 / [var (robust estimator of location)]
    while abs(empeff-eff)> eps
        
        bet= normcdf(c)-normcdf(-c);
        
        alph= 2*(c.^2*(1-normcdf(c))+normcdf(c)-0.5-c*normpdf(c));
        empeff=(bet^2)/alph;
        % disp(empeff)
        step=step/2;
        if empeff<eff
            c=c+step;
        elseif empeff>eff
            c=max(c-step,1e-10);
        else
        end
        
    end


ceff=c;
end

%FScategory:UTISTAT