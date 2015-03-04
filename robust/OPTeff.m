function ceff = OPTeff(eff,p,varargin)
%OPTeff finds the constant c which is associated to the requested efficiency
%
%  Required input arguments:
%
%    eff:       scalar which contains the required efficiency (of location
%               of scale estimator)
%               Generally eff=0.85, 0.9 or 0.95
%    p :        scalar, number of response variables
%
%TODO:OPTeff_INPUT_OPTIONS
%
% Output:
%
%  c = scalar of Optimal rho associated to the nominal (location or
%  shape) efficiency
%
% Copyright 2008-2015.
% Written by FSDA team
% Last modified 06-Feb-2015
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
%   Therefore, to obtain the value of c for the (rho) psi function defined in the
%   interval 0---2c, 2c---3c, >3c it is necessary to divide the output of
%   function OPTeff by 3.
%
%<a href="matlab: docsearchFS('OPTeff')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
%
%
% Examples:
%
%{
    % The constant c associated to a nominal location efficiency of 95% in regression is
    % c = 3.180662196584308
    c=OPTeff(0.95,1)
%}
%
%{
    % Find the value of c for efficiency which goes from 0.75 to 0.99 with step 0.01
    ef=0.75:0.01:0.99;
    CC=[ef' zeros(length(ef),1)];
    jk=0;
    for j=ef
        jk=jk+1;
        CC(jk,2)=OPTeff(j,1)
    end

%}
%


%

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
    
    
    empeff=10;
%     p4=(p+4)*(p+2);
%     p6=(p+6)*(p+4)*(p+2);
%     p8=(p+8)*(p+6)*(p+4)*(p+2);
    
 
    % Coefficients of optimal psi standardized using intervals 
    % 0---(2/3)c , (2/3)c---c, >c  
    p0=2*1.3846;
    p1=-2*2.6917;
    p2=4*10.7668;
    p3=-6*11.664;
    p4=8*4.0375;
    
    
    % Epsisq = E( \psi(x)^2) 
    % Epsidivx = E( \psi(x)/x)
    % Epsider = E( \psi'(x))
    % bet=(1-1/p)*Epsidivx+(1/p)*Epsider;   
    % [var (robust estimator of location using optimal rho function)] = (Epsisq/p) / (bet^2)
     
    while abs(empeff-eff)> eps
        
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
        
        empeff=(bet^2)/(Epsisq/p);
        
        step=step*0.5;
         % disp([step c empeff])
        if empeff<eff
            c=c+step;
        elseif empeff>eff
            c=max(c-step,0.1);
        else
        end
        % disp([step c empeff])
    end

ceff=c;
% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);

end