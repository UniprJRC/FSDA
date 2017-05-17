function ceff = OPTeff(eff,v)
%OPTeff finds the constant c which is associated to the requested efficiency
%
%<a href="matlab: docsearchFS('OPTeff')">Link to the help function</a>
%
%  Required input arguments:
%
%    eff:       required efficienty. Scalar.
%               Scalar which contains the required efficiency (of location
%               or scale estimator).
%               Generally eff=0.85, 0.9 or 0.95
%               Data Types - single|double
%    v :        Number of response variables. Scalar. e.g. in regression p=1
%               Data Types - single|double|int32|int64
%               
%
%  Optional input arguments: TODO_OPTeff_INPUT_OPTIONS
%  
%
%
% Output:
%
%  ceff : Requested tuning constant. Scalar. Tuning constatnt of optimal rho
%         function associated to requested value of efficiency
%
% More About:
%
%
%    $\rho$ ($\psi$) function which is considered is standardized 
%    using intervals 0---(2/3)c , (2/3)c---c, >c.   
%    $\rho$ function is
%
% \[
% TBrho(u)= \left\{
%    \begin{array}{lr}
%     1.3846 \left(\frac{u}{c}\right)^2                      &                                      |\frac{u}{c}| \leq  \frac{2}{3} \\
%    0.5514-2.6917 \left(\frac{u}{c}\right)^2 +10.7668\left(\frac{u}{c}\right)^4-11.6640\left(\frac{u}{c}\right)^6+4.0375\left(\frac{u}{c}\right)^8   & \qquad \frac{2}{3} \leq  |\frac{u}{c}|\leq  1 \\
%    1                                                    &                      |\frac{u}{c}|>1 \\
% \end{array}
%    \right.
%  \]
%                                                                      |t/c|>1                              
%
%   Therefore, to obtain the value of c for the (rho) psi function defined in the
%   interval 0---2c, 2c---3c, >3c it is necessary to divide the output of
%   function OPTeff by 3.
%
% See also: TBeff, HYPeff, HAeff
%
% References:
% 
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('OPTeff')">Link to the help page for this function</a>
% Last modified 31-05-2016
%
% Examples:
%
%{
    % Find c given a value of efficiency.
    % The constant c associated to a nominal location efficiency of 95% in regression is
    % c = 3.180662196584308
    c=OPTeff(0.95,1)
%}
%
%{
    % Find the value of c for efficiency which goes to 1.
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
        
        empeff=(bet^2)/(Epsisq/v);
        
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
%FScategory:UTISTAT