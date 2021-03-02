function [bdp,eff,approxsheff] = OPTc(c, v, shapeeff)
%OPTc computes breakdown point and efficiency associated with constant c for Optimal rho function
%
%<a href="matlab: docsearchFS('OPTc')">Link to the help function</a>
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
%  approxsheff :  Approximate value of efficiency. Scalar. 
%                 Approximate value of efficiency in case input
%                 option shapeeff=1 and v>1.
%                 This output is left for comparability with the value
%                 which comes out from R library robustbase.
%
%
%
% More About:
%
%
%    $\rho$ ($\psi$) function which is considered is standardized 
%    using intervals 0---(2/3)c , (2/3)c---c, >c.   
%    $\rho$ function is
%
% \[
% OPTrho(u)= \left\{
%    \begin{array}{lr}
%     1.3846 \left(\frac{u}{c}\right)^2                      &                                      |\frac{u}{c}| \leq  \frac{2}{3} \\
%    0.5514-2.6917 \left(\frac{u}{c}\right)^2 +10.7668\left(\frac{u}{c}\right)^4-11.6640\left(\frac{u}{c}\right)^6+4.0375\left(\frac{u}{c}\right)^8   & \qquad \frac{2}{3} \leq  |\frac{u}{c}|\leq  1 \\
%    1                                                    &                      |\frac{u}{c}|>1 \\
% \end{array}
%    \right.
%  \]
%                                                                      |t/c|>1                              
%   Therefore, the input c for the (rho) psi function above corresponds to c/3
%   in the rho (psi) function defined in the
%   interval 0---2c, 2c---3c, >3c 
%
%
% See also: TBc, HYPc, HAc
%
% References:
% 
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.                                                                            
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('OPTc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % bdp associated with a particular c.
    c=5;
    bdp = OPTc(c,1)
%}

%{
    % bdp and eff associated with a particular c.
    c=5;
    [bdp, eff] = OPTc(c,1)
%}

%{

    %% Breakdown vs efficiency.
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

%{
    % An example with 3 output arguments.
    c=5;
    % third input argument is 1, that is shape efficiency
    [bdp, eff, approxeff] = OPTc(c,2,1);
%}

%% Beginning of code



    b=c.^2/2;
    a=2*c.^2/9;
    bdp= (1.3846./(c.^2)).*v.*gammainc(a,0.5*(v+2)) ...
          + 0.5514.*(gammainc(b,v/2) -gammainc(a,v/2) ) ...
          -(2.6917./c.^2).*v.*(gammainc(b,(v+2)/2) -gammainc(a,(v+2)/2) )+...
          +(10.7668./c.^4).*v.*(v+2).*(gammainc(b,(v+4)/2) -gammainc(a,(v+4)/2) )+...
          -(11.664./c.^6).*v.*(v+2).*(v+4).*(gammainc(b,(v+6)/2) -gammainc(a,(v+6)/2) )+...
          +(4.0375./c.^8).*v.*(v+2).*(v+4).*(v+6).*(gammainc(b,(v+8)/2) -gammainc(a,(v+8)/2) )+...
          +(1-gammainc(b,v/2));

% Convergence condition is E(\rho) = \rho(c) bdp
%  where \rho(c) for standardized optimal rho function is 1


if nargin<=2 || shapeeff ~=1
    
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
        
        Epsisq=(p0^2/(c^4))*v*gammainc(a,(v+2)/2)+(p1^2/(c^4))*v*(gammainc(b,(v+2)/2)-gammainc(a,(v+2)/2))...
            + (p2^2/(c^8))*(v+4)*(v+2)*v*(gammainc(b,(v+6)/2)-gammainc(a,(v+6)/2))... 
            + ((p3^2+2*p2*p4)/(c^12))*(v+8)*(v+6)*(v+4)*(v+2)*v*(gammainc(b,(v+10)/2)-gammainc(a,(v+10)/2))... 
            + (p4^2/(c^16))*(v+12)*(v+10)*(v+8)*(v+6)*(v+4)*(v+2)*v*(gammainc(b,(v+14)/2)-gammainc(a,(v+14)/2))... 
            + (2*p1*p2/(c^6))*(v+2)*v*(gammainc(b,(v+4)/2)-gammainc(a,(v+4)/2))...
            + (2*p1*p3/(c^8))*(v+4)*(v+2)*v*(gammainc(b,(v+6)/2)-gammainc(a,(v+6)/2))...
            + (2*(p1*p4+p2*p3)/(c^10))*(v+6)*(v+4)*(v+2)*v*(gammainc(b,(v+8)/2)-gammainc(a,(v+8)/2))...
            + (2*p3*p4/(c^14))*(v+10)*(v+8)*(v+6)*(v+4)*(v+2)*v*(gammainc(b,(v+12)/2)-gammainc(a,(v+12)/2));
           
        Epsidivx=(p0/(c^2))*gammainc(a,v/2)+(p1/(c^2))*(gammainc(b,v/2)-gammainc(a,v/2))...
            + (p2/(c^4))*v*(gammainc(b,(v+2)/2)-gammainc(a,(v+2)/2))... 
            + (p3/(c^6))*(v+2)*v*(gammainc(b,(v+4)/2)-gammainc(a,(v+4)/2))... 
            + (p4/(c^8))*(v+4)*(v+2)*v*(gammainc(b,(v+6)/2)-gammainc(a,(v+6)/2)); 

        Epsider=(p0/(c^2))*gammainc(a,v/2)+(p1/(c^2))*(gammainc(b,v/2)-gammainc(a,v/2))...
            + (3*p2/(c^4))*v*(gammainc(b,(v+2)/2)-gammainc(a,(v+2)/2))... 
            + (5*p3/(c^6))*(v+2)*v*(gammainc(b,(v+4)/2)-gammainc(a,(v+4)/2))... 
            + (7*p4/(c^8))*(v+4)*(v+2)*v*(gammainc(b,(v+6)/2)-gammainc(a,(v+6)/2)); 
        
        bet=(1-1/v)*Epsidivx+(1/v)*Epsider;   
        
        eff=(bet^2)/(Epsisq/v);    
    
else
%TODO:OPTc:shapeff
    c2=c.^2/2;
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
    
    
    eff=(betsc.^2)/alphsc;
    
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