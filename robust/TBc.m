function [bdp,eff,approxsheff] = TBc(c,p,varargin)
%TBc computes breakdown point and efficiency associated with constant c for Tukey's biweight
%
%<a href="matlab: docsearch('tbc')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    c :        scalar greater than 0 which controls the robustness/efficiency of the estimator
%    p :        number of response variables of the dataset (for regression p=1)
%
%  Optional input arguments:
%
%   shapeeff : If 1, the efficiency is referred to the shape else (default)
%              is referred to the location
%
% Output:
%
%     bdp      :  scalar, breakdown point associated to the supplied
%                 value of c
%     eff      :  scalar, efficiency associated to the supplied
%                 value of c
%  approxsheff :  scalar, approximate value of efficiency in case input
%                 option shapeeff=1.
%                 This output is left for comparability with the value
%                 which comes out from R library robustbase.
%
%
% Copyright 2008-2011.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti
%
% References:
%
%
% Frank R. Hampel, Peter J. Rousseeuw and Elvezio Ronchetti (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% Journal of the American Statistical Association , Vol. 76, No. 375,
% pp. 643-648 (HRR)
%
%<a href="matlab: docsearch('tbc')">Link to the help page for this function</a>
% Last modified 15-Nov-2011
%
% Examples:

%{

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

%% beginning of code

c2=c.^2/2;

Erho= (p*gammainc(c2,0.5*(p+2))/2-(p^2+2*p)*gammainc(c2,0.5*(p+4))./(4*c2)+...
    +(p^3+6*p^2+8*p)*gammainc(c2,0.5*(p+6))./(6*(c.^4))+ ((c.^2)/6).*(1-gammainc(c2,p/2))  );
% Convergence condition is E(\rho) = \rho(c) bdp
%  where \rho(c) for TBW is c^2/6
bdp=(Erho./(c.^2))*(6);



if nargin<=2 || varargin{1} ~=1
    % LOCATION EFFICIENCY
    p4=(p+4)*(p+2);
    p6=(p+6)*(p+4)*(p+2);
    p8=(p+8)*(p+6)*(p+4)*(p+2);
    
    bet= p4*gammainc(c2,(p+4)/2)./((c.^4)) ...
        -2*(p+2)*gammainc(c2,(p+2)/2)./(c.^2)+...
        + gammainc(c2,p/2)  ;
    
    alph= p8*gammainc(c2,(p+10)/2)./(c.^8)-4*p6*gammainc(c2,(p+8)/2)./(c.^6)+...
        6*p4*gammainc(c2,(p+6)/2)./(c.^4)-2*(p+2)*gammainc(c2,(p+4)/2)./c2+...
        gammainc(c2,(p+2)/2);
    eff=(bet^2)/alph;
    
    
else
    % SCALE EFFICIENCY
    p4=(p+4);
    p6=p4*(p+6);
    p8=p6*(p+8);
    p10=p8*(p+10);
    
    
    betsc= gammainc(c2,(p+2)/2) ...
        -2*p4*gammainc(c2,(p+4)/2)./(c.^2)+...
        +p6*gammainc(c2,(p+6)/2)./(c.^4);
    
    % In the paper alphsc is A+B
    alphsc = gammainc(c2,(p+4)/2) ...
        -4*p4*gammainc(c2,(p+6)/2)./(c.^2)+...
        +6*p6*gammainc(c2,(p+8)/2)./(c.^4)+...
        -4*p8*gammainc(c2,(p+10)/2)./(c.^6)...
        +p10*gammainc(c2,(p+12)/2)./(c.^8);
    
    
    eff=(betsc.^2)/alphsc;
    
    if p>1
        % approxsheff is the approximate value of efficieny using Tyler
        % approximation.
        approxsheff=eff;
        
        % Now we compute the true value of efficiency
        k1=1/eff;
        
        Erho2=p*(p+2)*gammainc(c2,0.5*(p+4))./4-0.5*p*(p+2)*p4*gammainc(c2,0.5*(p+6))./(c.^2)...
            +(5/12)*p*(p+2)*p6*gammainc(c2,0.5*(p+8))./(c^4)-...
            -(1/6)*p*(p+2)*p8*gammainc(c2,0.5*(p+10))./(c^6)-...
            +(1/36)*p*(p+2)*p10*gammainc(c2,0.5*(p+12))./(c^8);
        Erho= p*gammainc(c2,0.5*(p+2))/2-(p*(p+2))*0.5*gammainc(c2,0.5*(p+4))./(c^2)+...
            +p*(p+2)*(p+4)*gammainc(c2,0.5*(p+6))./(6*(c.^4))+ (c2.*(1-gammainc(c2,p/2))  );
        Epsixx= p*gammainc(c2,0.5*(p+2))-2*(p*(p+2))*gammainc(c2,0.5*(p+4))./(c^2)+...
            +p*(p+2)*(p+4)*gammainc(c2,0.5*(p+6))./(c.^4);
        k3=(Erho2-Erho^2)/(Epsixx^2);
        
        effcor=((p-1)/p)*k1+2*k3;
        eff=1/effcor;
        
    end
end


end