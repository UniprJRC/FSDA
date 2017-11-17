function [bdp,eff,approxsheff] = TBc(c,v,shapeeff)
%TBc computes breakdown point and efficiency associated with constant c for Tukey's biweight
%
%<a href="matlab: docsearchFS('TBc')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    c :     tuning constant c. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%    v :        number of response variables. Scalar. Number of variables of
%               the  dataset (for regression v=1)
%
%  Optional input arguments:
%
%   shapeeff : location or shape efficiency. Scalar. 
%              If shapeeff=1, the efficiency is referred to the shape else
%              (default) is referred to the location estimator
%               Example - 1 
%               Data Types - double
%
% Output:
%
%     bdp      :  bdp. Scalar. Breakdown point associated to the supplied
%                 value of c
%     eff      :  eff. Scalar. Efficiency associated to the supplied
%                 value of c
%                 Remark: if c is a vector bdp and eff will also be vectors
%                 with the same size of c. For example bdp(3) and eff(3)
%                 are associated to c(3) ....
%  approxsheff :  Approximate value of efficiency. Scalar. 
%                 Approximate value of efficiency in case input
%                 option shapeeff=1 and v>1.
%                 This output is left for comparability with the value
%                 which comes out from R library robustbase.
%
%
% See also: HYPc, HAc, OPTc
%
% References:
%
%
% Frank R. Hampel, Peter J. Rousseeuw and Elvezio Ronchetti (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% Journal of the American Statistical Association , Vol. 76, No. 375,
% pp. 643-648 (HRR)
%
%
% Copyright 2008-2017.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('TBc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %Tbc with just one output argument.
    [bdp]=TBc(2,1)
    disp('Break down point')
    disp(bdp)
%}

%{
    %Tbc with 2 output arguments.
    [bdp,eff]=TBc(2,1)
    disp('Break down point and efficienty')
    disp(bdp)
    disp(eff)
%}

%{
    % Find also approximate value of scale efficienty (for R comparability).
    [bdp,eff,approxeff]=TBc(2,2,1)
%}

%{
    %Breakdown point and efficiency.
    %Analysis of breakdown point and asymptotic efficiency
    %at the normal distribution as a function of c in regression.
    c=1:0.01:6;
    [bdp,eff]=TBc(c,1);
    subplot(2,1,1)
    plot(c,bdp)
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Breakdown point','Interpreter','none')
    subplot(2,1,2)
    plot(c,eff)
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Asymptotic efficiency','Interpreter','none')

%}

%% Beginning of code

c2=c.^2/2;

Erho= (v*gammainc(c2,0.5*(v+2))/2-(v^2+2*v)*gammainc(c2,0.5*(v+4))./(4*c2)+...
    +(v^3+6*v^2+8*v)*gammainc(c2,0.5*(v+6))./(6*(c.^4))+ ((c.^2)/6).*(1-gammainc(c2,v/2))  );
% Convergence condition is E(\rho) = \rho(c) bdp
%  where \rho(c) for TBW is c^2/6
bdp=(Erho./(c.^2))*(6);



if nargin<=2 || shapeeff ~=1
    % LOCATION EFFICIENCY
    p4=(v+4)*(v+2);
    p6=(v+6)*(v+4)*(v+2);
    p8=(v+8)*(v+6)*(v+4)*(v+2);
    
    bet= p4*gammainc(c2,(v+4)/2)./((c.^4)) ...
        -2*(v+2)*gammainc(c2,(v+2)/2)./(c.^2)+...
        + gammainc(c2,v/2)  ;
    
    alph= p8*gammainc(c2,(v+10)/2)./(c.^8)-4*p6*gammainc(c2,(v+8)/2)./(c.^6)+...
        6*p4*gammainc(c2,(v+6)/2)./(c.^4)-2*(v+2)*gammainc(c2,(v+4)/2)./c2+...
        gammainc(c2,(v+2)/2);
    eff=(bet.^2)./alph;
    
    
else
    % SCALE EFFICIENCY
    p4=(v+4);
    p6=p4*(v+6);
    p8=p6*(v+8);
    p10=p8*(v+10);
    
    
    betsc= gammainc(c2,(v+2)/2) ...
        -2*p4*gammainc(c2,(v+4)/2)./(c.^2)+...
        +p6*gammainc(c2,(v+6)/2)./(c.^4);
    
    % In the paper alphsc is A+B
    alphsc = gammainc(c2,(v+4)/2) ...
        -4*p4*gammainc(c2,(v+6)/2)./(c.^2)+...
        +6*p6*gammainc(c2,(v+8)/2)./(c.^4)+...
        -4*p8*gammainc(c2,(v+10)/2)./(c.^6)...
        +p10*gammainc(c2,(v+12)/2)./(c.^8);
    
    
    eff=(betsc.^2)./alphsc;
    
    if v>1
        % approxsheff is the approximate value of efficieny using Tyler
        % approximation.
        approxsheff=eff;
        
        % Now we compute the true value of efficiency
        k1=1/eff;
        
        Erho2=v*(v+2)*gammainc(c2,0.5*(v+4))./4-0.5*v*(v+2)*p4*gammainc(c2,0.5*(v+6))./(c.^2)...
            +(5/12)*v*(v+2)*p6*gammainc(c2,0.5*(v+8))./(c^4)-...
            -(1/6)*v*(v+2)*p8*gammainc(c2,0.5*(v+10))./(c^6)-...
            +(1/36)*v*(v+2)*p10*gammainc(c2,0.5*(v+12))./(c^8);
        Erho= v*gammainc(c2,0.5*(v+2))/2-(v*(v+2))*0.5*gammainc(c2,0.5*(v+4))./(c^2)+...
            +v*(v+2)*(v+4)*gammainc(c2,0.5*(v+6))./(6*(c.^4))+ (c2.*(1-gammainc(c2,v/2))  );
        Epsixx= v*gammainc(c2,0.5*(v+2))-2*(v*(v+2))*gammainc(c2,0.5*(v+4))./(c^2)+...
            +v*(v+2)*(v+4)*gammainc(c2,0.5*(v+6))./(c.^4);
        k3=(Erho2-Erho^2)/(Epsixx^2);
        
        effcor=((v-1)/v)*k1+2*k3;
        eff=1/effcor;
        
    end
end


end
%FScategory:UTISTAT