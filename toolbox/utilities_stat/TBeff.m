function ceff = TBeff(eff, v, shapeeff, approxsheff)
%Tbeff finds the constant c which is associated to the requested efficiency for Tukey's biweight
%
%
%
%<a href="matlab: docsearchFS('TBeff')">Link to the help function</a>
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
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
%   shapeeff : Location or shape efficiency. Scalar.
%              If 1, the efficiency is referred to shape else (default)
%              is referred to location
%               Example - 'shapeeff',1
%               Data Types - double
%
% approxsheff: approximate or exact calculations. Dummy scalar.
%              If 1, when p > 1 the approximate formula for scale
%              efficiency is used else (default) the exact formula of the
%              variance of the robust estimator of the scale is used
%               Example - 'approxsheff',1
%               Data Types - double
%
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
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('TBeff')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

%
% Examples:
%
%{
    % Find c for a given efficiency.
    % The constant c associated to a nominal location efficiency of 95% in regression is
    % c = 4.685064948559557
    c=TBeff(0.95,1)
%}

%{
    % Find c associated to a fixed shape efficiency 
    % c= 5.490249208447494
    c=TBeff(0.95,1,1)
%}
%
%
%
%{
    % Find shape efficiency.
    % Compare constant c for the range of values of p from 1 to 200 using
    % approximate and exact formula for fixed shape efficiency
    % Number of variables
    p=200;
    pp=1:1:p;
    % Initialize the matrix which stores the values of c for the two
    % methods
    cc=[pp' zeros(p,2)];

    eff=0.75;
    % eff=0.99;
    for i=pp
        % Use exact formula for finding the value of c associated to a fixed
        % level of shape efficiency
        cc(i,2)=TBeff(eff,i,1);
        % Use approximate formula for finding the value of c associated to a fixed
        % level of shape efficiency
        cc(i,3)=TBeff(eff,i,1,1);
    end
    figure
    plot(cc(:,1),cc(:,2),'LineStyle','-','LineWidth',2)
    hold('on')
    plot(cc(:,1),cc(:,3),'LineStyle','-.','LineWidth',2)
    legend('True value of c', 'Approximate value of c','Location','best')
    xlabel('Number of variables')
    ylabel('Value of c')
%}

%

%% Beginning of code

eps=1e-12;
if nargin<=2 || shapeeff  ~=1
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
    p4=(v+4)*(v+2);
    p6=(v+6)*p4;
    p8=(v+8)*p6;
    
    
    % bet  = \int_{-c}^c  \psi'(x) d \Phi(x)
    % alph = \int_{-c}^c  \psi^2(x) d \Phi(x)
    % empeff = bet^2/alph = 1 / [var (robust estimator of location)]
    while abs(empeff-eff)> eps
        
        cs=c.^2/2;
        
        bet= p4*gammainc(cs,(v+4)/2)./((c.^4)) ...
            -2*(v+2)*gammainc(cs,(v+2)/2)./(c.^2)+...
            + gammainc(cs,v/2)  ;
        
        alph= p8*gammainc(cs,(v+10)/2)./(c.^8)-4*p6*gammainc(cs,(v+8)/2)./(c.^6)+...
            6*p4*gammainc(cs,(v+6)/2)./(c.^4)-2*(v+2)*gammainc(cs,(v+4)/2)./cs+...
            gammainc(cs,(v+2)/2);
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
    % SHAPE EFFICIENCY
    if nargin<=3 || approxsheff ~=1
        % approxsheff 0 ==> use exact formulae
        approxsheff=0;
    else
        % approxsheff 1 ==> use Tyler approximation
        approxsheff=1;
    end
    % constant for second Tukey Biweight rho-function for MM, for fixed shape-efficiency
    % c = starting point of the iteration
    c=sqrt(chi2inv(eff,v))+7;
    % step = width of the dichotomic search (it decreases by half at each
    % iteration).
    if eff<0.92 && approxsheff==0
        step=5;
    else
        step=15;
    end
    varrobestsc=10;
    p4=(v+4);
    p6=p4*(v+6);
    p8=p6*(v+8);
    p10=p8*(v+10);
    % alphsc = E[ \psi^2(x) x^2] /{(v(v+2)]^2}
    % betsc =  E [ \psi'(x) x^2+(v+1)  \psi^2(x) x ]/[v(v+2)]
    % res = [var (robust estimator of scale)] = alphsc/(betsc^2)
    while abs(1-eff*varrobestsc)> eps
        cs=c.^2/2;
        alphsc = gammainc(cs,(v+4)/2) ...
            -4*p4*gammainc(cs,(v+6)/2)./(c.^2)+...
            +6*p6*gammainc(cs,(v+8)/2)./(c.^4)+...
            -4*p8*gammainc(cs,(v+10)/2)./(c.^6)...
            +p10*gammainc(cs,(v+12)/2)./(c.^8);
        
        betsc= gammainc(cs,(v+2)/2) ...
            -2*p4*gammainc(cs,(v+4)/2)./(c.^2)+...
            +p6*gammainc(cs,(v+6)/2)./(c.^4);
        
        varrobestsc=alphsc./(betsc.^2);
        
        if v>1 && approxsheff==0
            
            Erho2=v*(v+2)*gammainc(cs,0.5*(v+4))./4-0.5*v*(v+2)*p4*gammainc(cs,0.5*(v+6))./(c.^2)...
                +(5/12)*v*(v+2)*p6*gammainc(cs,0.5*(v+8))./(c^4)-...
                -(1/6)*v*(v+2)*p8*gammainc(cs,0.5*(v+10))./(c^6)-...
                +(1/36)*v*(v+2)*p10*gammainc(cs,0.5*(v+12))./(c^8);
            Erho= v*gammainc(cs,0.5*(v+2))/2-(v*(v+2))*0.5*gammainc(cs,0.5*(v+4))./(c^2)+...
                +v*(v+2)*(v+4)*gammainc(cs,0.5*(v+6))./(6*(c.^4))+ (cs.*(1-gammainc(cs,v/2))  );
            Epsixx= v*gammainc(cs,0.5*(v+2))-2*(v*(v+2))*gammainc(cs,0.5*(v+4))./(c^2)+...
                +v*(v+2)*(v+4)*gammainc(cs,0.5*(v+6))./(c.^4);
            k3=(Erho2-Erho^2)/(Epsixx^2);
            varrobestsc=((v-1)/v)*varrobestsc+2*k3;
        end
        
        % disp(1-eff*varrobestsc)
        
        step=step/2;
        if (1-eff*varrobestsc)<0
            c=c+step;
        elseif (1-eff*varrobestsc)>0
            c=max(c-step,0.1);
        end
        %  disp([c 1-eff*res])
        
    end
end

ceff=c;
% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);

end
%FScategory:UTISTAT
