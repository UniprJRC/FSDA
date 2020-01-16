function [c, M] = RKeff(eff,v,ARP)
%RKeff finds the constants c and M which are associated to the requested efficiency and ARP
%
%<a href="matlab: docsearchFS('RKeff')">Link to the help function</a>
%
%  Required input arguments:
%
%    eff:       required efficiency. Scalar.
%               Scalar which contains the required efficiency (of location
%               or scale estimator).
%               Generally eff=0.85, 0.9 or 0.95
%               Data Types - single|double
%    v :        Number of response variables. Scalar. e.g. in regression v=1
%               Data Types - single|double|int32|int64
%
%  Optional input arguments:
%
%    ARP      :  asymptotic rejection probability. Scalar.
%                The asymptotic rejection probability of an estimator is
%                defined as the probability in large sample under a
%                reference distribution that a Malanobis distance excees
%                $c_0$, where $c_0=inf \{ u_0 | w(u)=0, \forall u>u_0 \}$.
%                $w(u)$ is the weight function (defined in RKwei.m). In
%                other words, given $c_0=sup(\rho(u))$,if an estimator is normed
%                to the normal distribution ARP is $1-F_{\chi^2_v}(c_0^2)$.
%                The default value of ARP is 0.05.
%                 Example - 0.04
%               Data Types - double
%
% Output:
%
%  c : Requested tuning constant. Scalar. Tuning constatnt of Rocke rho
%       function (translated Tukey Biweight) associated to requested
%       efficiency and asymptotic rejection probability
%  M : Requested tuning constant. Scalar. Tuning constant of Rocke rho
%       function (translated Tukey Biweight) associated to requested
%       efficiency and asymptotic rejection probability
%
%
% More About:
%
%
%
% See also: TBeff, HYPeff, HAeff, HUeff, OPTeff
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('RKeff')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit
%
% Examples:
%
%{
    % Find c given a value of efficiency.
    % The constant c associated to a nominal location efficiency of 95% in regression is
    % c = 3.180662196584308
    [c,M]=RKeff(0.95,5,0.05)
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


if nargin<3
    ARP=0.05;
end


eps=1e-12;

% deps = tolerance to compute the numerical derivative
deps = 1e-4;

a=1;
prod1=exp(gammaln((v+a)/2)-gammaln(v/2))*2^(a/2);
prod2=v;
a=3;
prod3=exp(gammaln((v+a)/2)-gammaln(v/2))*2^(a/2);
% a=4;
% prod4=prod(v:2:(v+a-2));
prod4=v*(v+2);
a=5;
prod5=exp(gammaln((v+a)/2)-gammaln(v/2))*2^(a/2);
% a=6;
% prod6=prod(v:2:(v+a-2));
prod6=v*(v+2)*(v+4);
a=7;
prod7=exp(gammaln((v+a)/2)-gammaln(v/2))*2^(a/2);
% prod8=prod(v:2:(v+a-2));
prod8=v*(v+2)*(v+4)*(v+6);
a=9;
prod9=exp(gammaln((v+a)/2)-gammaln(v/2))*2^(a/2);

% prod10=prod(v:2:(v+a-2));
prod10=v*(v+2)*(v+4)*(v+6)*(v+8);



% cplusM = initial value of c+M
cplusM = sqrt(chi2inv(1-ARP,v));

M = sqrt(v);
% Initial value for c
c = cplusM - M;

iter = 1;
crit = 100;

% maxiter = maximum number of iterations
maxiter=10000;

while (crit > eps && iter<maxiter)
    
    Mold = M;
    Erho = EmpEff(c,v,M);
    % iterate up to when crit=|F| is <=eps
    % crit= |F| = |E(rho)-bdp*sup(rho)| <=eps
    % sup(rho)=(M^2/2+c*(5*c+16*M)/30)
    F = Erho - eff;
    
    % Fderc = numerical partial derivative of E(rho) with respect to c
    % Fderc =(E(rho(c1+deps))-E(rho(c1)))/1e-4
    Fderc =(EmpEff(c+deps,v,M) - Erho)/deps;
    
    % FderM=numerical partial derivative of E(rho) with respect to M
    % FderM =(E(rho(M+deps))-E(rho(c1)))/1e-4
    FderM  = (EmpEff(c,v,M+deps) - Erho)/deps;
    
    
    % Fder = derivative of F with respect to M, minus derivative of F
    % with respect to c
    Fder = FderM - Fderc;
    
    % Update equation for M
    % Note that this is the typical equation x_new=x_old-f(x_old)/f'(x_old)
    M = M - F/Fder;
    
    if M >= cplusM
        M = (Mold + cplusM)/2;
    end
    c= cplusM - M;
    crit = abs(F);
    iter=iter+1;
    disp([iter crit])
    disp('----')
end

if iter==maxiter
    error('FSDA:RKeff:NoConv','No convergence in the computation of c and M')
end

% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);

    function EF = EmpEff(ct,p,M)
        % This routine computes the empirical efficiency
        % given the values of c, p and M
        
        % Epsisq = E( \psi(x)^2)
        % Epsidivx = E( \psi(x)/x)
        % Epsider = E( \psi'(x))
        % bet=(1-1/p)*Epsidivx+(1/p)*Epsider;
        % [var (robust estimator of location using optimal rho function)] = (Epsisq/p) / (bet^2)
        
        % Notation:
        % Mct2=(M+c)^2, M2=M^2, ct2=c^2, ct4=c^4
        
        Mct2= (M+ct)^2;
        M2=M^2;
        c2=ct^2;
        %         c4=ct^4;
        %         M2mc2=3*M2-c2;
        %         c2mM2=(c2-M2)^2;
        
        %         Epsisq=prod2*chi2cdf(M2,p+2)...
        %               +(16*prod4*(chi2cdf(Mct2,p+4)-chi2cdf(M2,p+4))*M2*c2mM2 ...
        %             +16*prod5*(chi2cdf(Mct2,p+5)-chi2cdf(M2,p+5))*M*(c2-M2)*(3*M2-c2) ...
        %             +4*prod6*(chi2cdf(Mct2,p+6)-chi2cdf(M2,p+6))*(M2mc2^2-8*M2*(c2^2-M2^2))...
        %             +4*prod7*(chi2cdf(Mct2,p+7)-chi2cdf(M2,p+7))*M*(2*(c2-M2) -4*M2mc2)...
        %             +prod8*(chi2cdf(Mct2,p+8)-chi2cdf(M2,p+8))*(28*M2 -c2)...
        %             -8*prod9*(chi2cdf(Mct2,p+9)-chi2cdf(M2,p+9))*M...
        %             +prod10*(chi2cdf(Mct2,p+10)-chi2cdf(M2,p+10)))/(c^8);
        
        
        Epsisq=prod2*chi2cdf(M2,p+2)...
            + prod10*Euj(Mct2,M2,p+10)/ct^8 + (-(8*M)/ct^8)*prod9*Euj(Mct2,M2,p+9)...
            + ((28*M^2)/ct^8 - 4/ct^6)*prod8*Euj(Mct2,M2,p+8)...
            + ((24*M)/ct^6 - (56*M^3)/ct^8)*prod7*Euj(Mct2,M2,p+7)...
            + (6/ct^4 - (60*M^2)/ct^6 + (70*M^4)/ct^8)*prod6*Euj(Mct2,M2,p+6)...
            + ((80*M^3)/ct^6 - (24*M)/ct^4 - (56*M^5)/ct^8)*prod5*Euj(Mct2,M2,p+5)...
            + ((36*M^2)/ct^4 - 4/ct^2 - (60*M^4)/ct^6 + (28*M^6)/ct^8)*prod4*Euj(Mct2,M2,p+4)...
            + ((8*M)/ct^2 - (24*M^3)/ct^4 + (24*M^5)/ct^6 - (8*M^7)/ct^8)*prod3*Euj(Mct2,M2,p+3)...
            + ((6*M^4)/ct^4 - (4*M^2)/ct^2 - (4*M^6)/ct^6 + M^8/ct^8 + 1)*prod2*Euj(Mct2,M2,p+2);
  
        %         Epsidivx=chi2cdf(M2,p)+(chi2cdf(Mct2,p)-chi2cdf(M2,p))*(1-2*M2/c2) ...
        %             + prod4*(chi2cdf(Mct2,p+4)-chi2cdf(M2,p+4))*(1/c4) ...
        %             -2*prod3*(chi2cdf(Mct2,p+3)-chi2cdf(M2,p+3))*(M/c4) ...
        %             + prod2*(chi2cdf(Mct2,p+2)-chi2cdf(M2,p+2))*(6*M2/c2-2) ...
        %             -4*prod1*(chi2cdf(Mct2,p+1)-chi2cdf(M2,p+1))*4*(M/c2)*(M2/c2-1);
        
        Epsidivx= chi2cdf(M2,p)+prod4*Euj(Mct2,M2,p+4)/ct^4 ...
            + (-(4*M)/ct^4)*prod3*Euj(Mct2,M2,p+3) ...
            + ((6*M^2)/ct^4 - 2/ct^2)*prod2*Euj(Mct2,M2,p+2) ...
            + ((4*M)/ct^2 - (4*M^3)/ct^4)*prod1*Euj(Mct2,M2,p+1) ...
            + (M^4/ct^4 - (2*M^2)/ct^2 + 1)*( chi2cdf(Mct2,p)-chi2cdf(M2,p) );
        
        %         Epsider=Epsidivx+prod4*( chi2cdf(Mct2,p+4)-chi2cdf(M2,p+4))*(1/c4) ...
        %             -4*prod2*(chi2cdf(Mct2,p+2)-chi2cdf(M2,p+2))*(M2+ c2+2*M/c2)*(1/c2) ...
        %             +4*prod3*(chi2cdf(Mct2,p+3)-chi2cdf(M2,p+3))*(4*M/c4) ...
        %             +4*prod1*(chi2cdf(Mct2,p+1)-chi2cdf(M2,p+1))*(M/c2)*(M2+c2);
        Epsider=chi2cdf(M2,p)+ (5/ct^4)*prod4*Euj(Mct2,M2,p+4)...
            -(16*M)/(ct^4)*prod3* Euj(Mct2,M2,p+3) ...
            +((6*(M^2/ct^2 - 1))/ct^2) *prod2* Euj(Mct2,M2,p+2)  ...
            +(-(8*M*(M^2/ct^2 - 1))/ct^2)*prod1*Euj(Mct2,M2,p+1) ...
            +(M^2/ct^2 - 1)^2*Euj(Mct2,M2,p) ;
  
        bet=(1-1/p)*Epsidivx+(1/p)*Epsider;
        
        EF=(bet^2)/(Epsisq/p);
    end

    function res=Euj(b,a,j)
        res=gammainc(0.5*b,0.5*j)-gammainc(0.5*a,0.5*j);
        % res=chi2cdf(b,j)-chi2cdf(a,j);
    end
end
%FScategory:UTISTAT