function [bdp,eff,approxsheff] = OPTc(c,p,varargin)
%OPTc computes breakdown point and efficiency associated with constant c for Optimal rho function
%
%<a href="matlab: docsearchFS('optc')">Link to the help function</a>
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
%  The output consists of a structure 'out' containing the following fields:
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
%    REMARK: \rho (\psi) function which is considered is standardized 
%    using intervals 0---(2/3)c , (2/3)c---c, >c   
%    Rho function is
%
%               |   1.3846 |r/c|^2                                                         |t/c|<=2/3
%               |   
%     \rho(r) = |   0.5514-2.6917|r/c|^2+10.7668|r/c|^4-11.6640|r/c|^6+4.0375|r/c|^8     2/3<=|t/c|<=1
%               |
%               |   1                                                                          |t/c|>1                              
%
%   Therefore, the input c for the (rho) psi function above corresponds to c/3
%   in the rho (psi) function defined in the
%   interval 0---2c, 2c---3c, >3c 
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('OPTc')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%{

    %Analysis of breakdown point and asymptotic efficiency
    %at the normal distribution as a function of c in regression.
    c=1:0.01:4;
    CC=[c' zeros(length(c),1)];
    jk=0;
    for j=c
        jk=jk+1;
         [bdp,eff]=OPTc(j,1);
        CC(jk,2:3)=[bdp,eff];
    end

    
    subplot(2,1,1)
    plot(c',CC(:,2))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Breakdown point','Interpreter','none')
    subplot(2,1,2)
    plot(c',CC(:,3))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Asymptotic efficiency','Interpreter','none')

%}

%% Beginning of code



    b=c.^2/2;
    a=2*c.^2/9;
    bdp= (1.3846./(c.^2)).*p.*gammainc(a,0.5*(p+2)) ...
          + 0.5514.*(gammainc(b,p/2) -gammainc(a,p/2) ) ...
          -(2.6917./c.^2).*p.*(gammainc(b,(p+2)/2) -gammainc(a,(p+2)/2) )+...
          +(10.7668./c.^4).*p.*(p+2).*(gammainc(b,(p+4)/2) -gammainc(a,(p+4)/2) )+...
          -(11.664./c.^6).*p.*(p+2).*(p+4).*(gammainc(b,(p+6)/2) -gammainc(a,(p+6)/2) )+...
          +(4.0375./c.^8).*p.*(p+2).*(p+4).*(p+6).*(gammainc(b,(p+8)/2) -gammainc(a,(p+8)/2) )+...
          +(1-gammainc(b,p/2));

% Convergence condition is E(\rho) = \rho(c) bdp
%  where \rho(c) for standardized optimal rho function is 1


if nargin<=2 || varargin{1} ~=1
    
        % Coefficients of optimal psi standardized using intervals 
    % 0---(2/3)c , (2/3)c---c, >c  
    p0=2*1.3846;
    p1=-2*2.6917;
    p2=4*10.7668;
    p3=-6*11.664;
    p4=8*4.0375;

    % LOCATION EFFICIENCY
        b=c.^2/2;
        a=2*c.^2/9;
        
        Epsisq=(p0^2/(c^4))*p*gammainc(a,(p+2)/2)+(p1^2/(c^4))*p*(gammainc(b,(p+2)/2)-gammainc(a,(p+2)/2))...
            + (p2^2/(c^8))*(p+4)*(p+2)*p*(gammainc(b,(p+6)/2)-gammainc(a,(p+6)/2))... 
            + ((p3^2+2*p2*p4)/(c^12))*(p+8)*(p+6)*(p+4)*(p+2)*p*(gammainc(b,(p+10)/2)-gammainc(a,(p+10)/2))... 
            + (p4^2/(c^16))*(p+12)*(p+10)*(p+8)*(p+6)*(p+4)*(p+2)*p*(gammainc(b,(p+14)/2)-gammainc(a,(p+14)/2))... 
            + (2*p1*p2/(c^6))*(p+2)*p*(gammainc(b,(p+4)/2)-gammainc(a,(p+4)/2))...
            + (2*p1*p3/(c^8))*(p+4)*(p+2)*p*(gammainc(b,(p+6)/2)-gammainc(a,(p+6)/2))...
            + (2*(p1*p4+p2*p3)/(c^10))*(p+6)*(p+4)*(p+2)*p*(gammainc(b,(p+8)/2)-gammainc(a,(p+8)/2))...
            + (2*p3*p4/(c^14))*(p+10)*(p+8)*(p+6)*(p+4)*(p+2)*p*(gammainc(b,(p+12)/2)-gammainc(a,(p+12)/2));
           
        Epsidivx=(p0/(c^2))*gammainc(a,p/2)+(p1/(c^2))*(gammainc(b,p/2)-gammainc(a,p/2))...
            + (p2/(c^4))*p*(gammainc(b,(p+2)/2)-gammainc(a,(p+2)/2))... 
            + (p3/(c^6))*(p+2)*p*(gammainc(b,(p+4)/2)-gammainc(a,(p+4)/2))... 
            + (p4/(c^8))*(p+4)*(p+2)*p*(gammainc(b,(p+6)/2)-gammainc(a,(p+6)/2)); 

        Epsider=(p0/(c^2))*gammainc(a,p/2)+(p1/(c^2))*(gammainc(b,p/2)-gammainc(a,p/2))...
            + (3*p2/(c^4))*p*(gammainc(b,(p+2)/2)-gammainc(a,(p+2)/2))... 
            + (5*p3/(c^6))*(p+2)*p*(gammainc(b,(p+4)/2)-gammainc(a,(p+4)/2))... 
            + (7*p4/(c^8))*(p+4)*(p+2)*p*(gammainc(b,(p+6)/2)-gammainc(a,(p+6)/2)); 
        
        bet=(1-1/p)*Epsidivx+(1/p)*Epsider;   
        
        eff=(bet^2)/(Epsisq/p);    
    
else
%TODO:OPTc:shapeff
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