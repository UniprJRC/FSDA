function [c,M] = RKbdp(bdp,v,ARP)
%RKbdp finds the constants c associated to the supplied breakdown point and asymptotic rejection point
%
%
%<a href="matlab: docsearchFS('RKbdp')">Link to the help function</a>
%
%  Required input arguments:
%
%      bdp    : breakdown point. Scalar. Scalar defining breakdown point
%               (i.e a number in the interval [0 0.5). Please notice that the
%               maximum achievable breakdown point is (n-p)/(2*n), and
%               therefore the value 0.5 is reached only when the sample
%               size goes to infinity. However, this routine assume a
%               sample of size infinity and allows you to specify a bdp
%               equal to 0.5.
%               Data Types - single|double
%        v    : number of response variables. Scalar. e.g. in regression v=1
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
%       breakdown point and asymptotic rejection probability
%  M : Requested tuning constant. Scalar. Tuning constant of Rocke rho
%       function (translated Tukey Biweight) associated to requested
%       breakdown point and asymptotic rejection probability
%
%
% See also: TBbdp, OPTbdp, HYPbdp, HAbdp
%
% References:
%
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
% Rocke D.M. (1996), Robustness properties of S-estimators of multivariate
% location and shape in high dimension, The Annals of Statistics, Vol. 24,
% pp. 1327-1345.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('RKbdp')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit
%
%
%
% Examples:
%
%{
    % Find c and M given bdp and ARP.
    bdp=0.5;
    % Suppose the dimension is 3
    v=3;
    % The constants c and M associated to a breakdown point of 50 per cent
    % and an ARP of 0.05 when there are three variables are
    % c=1.133183024897769 and M= 1.662300458017338
    [c,M]=RKbdp(bdp,v)
%}

%{
    %% Computation of c and M for a series of values of bdp (v=3).
    v=3;
    bdp=0.01:0.01:0.5;
    cM=zeros(length(bdp),2);
    for i=1:length(bdp)
     [c,M]=RKbdp(bdp(i),v);
     cM(i,:)=[c M];
    end

    subplot(2,1,1)
    plot(bdp,cM(:,1))
    subplot(2,1,2)
    plot(bdp,cM(:,2))
%}

%{
    %% Computation of c and M for a series of values of bdp (v=10).
    % ARP fixed to 0.01.
    v=10;
    bdp=0.01:0.01:0.5;
    cM=zeros(length(bdp),2);
    for i=1:length(bdp)
     [c,M]=RKbdp(bdp(i),v,0.01);
     cM(i,:)=[c M];
    end

    subplot(2,1,1)
    plot(bdp,cM(:,1))
    subplot(2,1,2)
    plot(bdp,cM(:,2))
%}

%% Beginning of code

if nargin<3
    ARP=0.05;
end


iter = 1;
crit = 100;
% eps = tolerance for the main loop
eps = 1e-5;
% deps = tolerance to compute the numerical derivative
deps = 1e-4;

% Constants prod2, ..., prod6 are fixed once and for all and
% can be computed out of the main loop
% a=2;
% prod2=prod(v:2:(v+a-2));
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

% The preliminary loop below is just to check whether ARP is reachable
% When c goes to 0 the estimator becomes the lease Winsorized squares (LWS)
% estimtor
% \rho_{LWS}= d^2/2 if      0 \leq d \leq c
%             c^2/2          d>c

c=2*v;
%     .chiInt(p,2,M)/2 + M^2/2*.chiInt2(p,0,M))
% maxiter = maximum number of iterations
maxiter=10000;

while (crit > eps && iter<maxiter)
    cold=c;
    c2=c^2;
    % Note that erhoLim computes E(d^2) and not E(d^2/2)
    % This is fine because Flim is E(d^2)-bdp*c^2 rather than
    %  E(d^2)/2-bdp*c^2/2
    % Slow implementation calling function chi2cdf
    % ErhoLim=prod2*chi2cdf(c2,v+2)+ c2*(1-chi2cdf(c2,v));
    % Fast implementation calling directly function gammainc
    % Remark: chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);
    ErhoLim=prod2*gammainc(0.5*c2,0.5*(v+2))+ c2*(1-gammainc(0.5*c2,0.5*v));
    
    FLim=ErhoLim-c2*bdp;
    % ErhoLimder = derivative of Erholim with respect to c
    % Slow implementation calling function chi2cdf
    % ErhoLimder=prod2*chi2pdf(c2,v+2)*2*c +2*c*(1-chi2cdf(c2,v))-c^2*chi2pdf(c2,v)*2*c;
    % Fast implementation calling directly function gammainc
    % Remark: chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);
    ErhoLimder=prod2*gampdf(c2,0.5*(v+2),2)*2*c +2*c*(1-gammainc(0.5*c2,0.5*v))-c^2*gampdf(c2,0.5*v,2)*2*c;
    
    Flimder=ErhoLimder-2*c*bdp;
    
    % Update equation for c
    % Note that this is the typical Newton Raphson equation x_new=x_old-f(x_old)/f'(x_old)
    % https://en.wikipedia.org/wiki/Newton%27s_method
    c = c - FLim/Flimder;
    
    if c < 0
        c = cold/2;
    end
    crit = abs(FLim);
    iter=iter+1;
end

if iter==maxiter
    error('FSDA:RKbdp:NoConv','No convergence in the computation of maxARP')
end

% If the estimator is normed to the normal distribution the asymptotic
% rejection probability, that is the probability that (under multivariate
% normality) a Mahalanobis distance exceeds c_1 is $1-F_\chi^2_v(c_1^2)$
maxARP=1-chi2cdf(c^2,v);

if maxARP <= ARP
    % if rejection is not achievable, use c=0 and best rejection
    c=0;
    M=sqrt(chi2inv(1-ARP, v));
else
    
    % cplusM = initial value of c+M
    cplusM = sqrt(chi2inv(1-ARP,v));
    
    M = sqrt(v);
    % Initial value for c
    c = cplusM - M;
    
    iter = 1;
    crit = 100;
    
    while (crit > eps && iter<maxiter)
        
        Mold = M;
        Erho = erho(c,v,M);
        % iterate up to when crit=|F| is <=eps
        % crit= |F| = |E(rho)-bdp*sup(rho)| <=eps
        % sup(rho)=(M^2/2+c*(5*c+16*M)/30)
        F = Erho - bdp*(M^2/2+c*(5*c+16*M)/30);
        
        % Fderc = numerical partial derivative of E(rho) with respect to c
        % Fderc =(E(rho(c1+deps))-E(rho(c1)))/1e-4
        Fderc =(erho(c+deps,v,M) - Erho)/deps;
        
        % FderM=numerical partial derivative of E(rho) with respect to M
        % FderM =(E(rho(M+deps))-E(rho(c1)))/1e-4
        FderM  = (erho(c,v,M+deps) - Erho)/deps;
        
        % (14M+6c)/30 is the derivative of suprho=M^2/2+c*(5*c+16*M)/30 with
        % respect to M minus the derivative of suprho with respect to c
        
        % Fder = derivative of F with respect to M, minus derivative of F
        % with respect to c -(derivative of suprho with respect to M -
        % derivative of suprho with respect to c)
        Fder = FderM - Fderc - bdp*(14*M+6*c)/30;
        
        % Update equation for M
        % Note that this is the typical equation x_new=x_old-f(x_old)/f'(x_old)
        M = M - F/Fder;
        
        if M >= cplusM
            M = (Mold + cplusM)/2;
        end
        c= cplusM - M;
        crit = abs(F);
        iter=iter+1;
        %         disp([iter crit])
        %         disp('----')
    end
    
    if iter==maxiter
        error('FSDA:RKbdp:NoConv','No convergence in the computation of c and M')
    end
    
end

    function Expectationrho = erho(ct,p,M)
        % This routine computes the expectation of Rocke rho function for
        % a given set of values of ct, p and M
        % Notation:
        % Mct2=(M+c)^2, M2=M^2, ct2=c^2, ct4=c^4
        Mct2= (M+ct)^2;
        M2=M^2;
        ct2=ct^2;
        ct4=ct^4;
        
        %       %  Slow implementation using function chi2cdf (which does a series
        %       %  of checks than calls gamcdf, which finally calls gammainc)
        %                 if (ct == 0)
        %                     %     .chiInt(p,2,M)/2 + M^2/2*.chiInt2(p,0,M))
        %                     Expectationrho=prod2*chi2cdf(M2,p+2)/2+ M2/2*(1-chi2cdf(M2,p));
        %                 else
        %                     Expectationrho=prod2*chi2cdf(M2,p+2)/2 ...
        %                         +(M2/2+ct*(5*ct+16*M)/30)*(1-chi2cdf(Mct2,p)) ...
        %                         +(M2/2-M2*(M^4-5*M2*ct2+15*ct4)/(30*ct4))*(chi2cdf(Mct2,p)-chi2cdf(M2,p)) ...
        %                         +(1/2+M^4/(2*ct4)-M2/ct2)*prod2*(chi2cdf(Mct2,p+2)-chi2cdf(M2,p+2)) ...
        %                         +(4*M/(3*ct2)-4*M^3/(3*ct4))*prod3*(chi2cdf(Mct2,p+3)-chi2cdf(M2,p+3)) ...
        %                         +(3*M2/(2*ct4)-1/(2*ct2))*prod4*(chi2cdf(Mct2,p+4)-chi2cdf(M2,p+4))...
        %                         -(4*M/(5*ct4))*prod5*(chi2cdf(Mct2,p+5)-chi2cdf(M2,p+5))...
        %                         +(1/(6*ct4))*prod6*(chi2cdf(Mct2,p+6)-chi2cdf(M2,p+6));
        %                 end
        
        % Fast implementation using directly function gammainc
        if (ct == 0)
            %     .chiInt(p,2,M)/2 + M^2/2*.chiInt2(p,0,M))
            Expectationrho=prod2*chi2cdf(M2,p+2)/2+ M2/2*(1-chi2cdf(M2,p));
        else
            % Expection of rho calling directly gammmainc
            %             Expectationrho=prod2*gammainc(0.5*M2,0.5*(p+2))/2 ...
            %                 +(M2/2+ct*(5*ct+16*M)/30)*(1-gammainc(0.5*Mct2,0.5*p)) ...
            %                 +(M2/2-M2*(M^4-5*M2*ct2+15*ct4)/(30*ct4))*(gammainc(0.5*Mct2,0.5*p)-gammainc(0.5*M2,0.5*p)) ...
            %                 +(1/2+M^4/(2*ct4)-M2/ct2)*prod2*(gammainc(0.5*Mct2,0.5*(p+2))-gammainc(0.5*M2,0.5*(p+2))) ...
            %                 +(4*M/(3*ct2)-4*M^3/(3*ct4))*prod3*(gammainc(0.5*Mct2,0.5*(p+3))-gammainc(0.5*M2,0.5*(p+3))) ...
            %                 +(3*M2/(2*ct4)-1/(2*ct2))*prod4*(gammainc(0.5*Mct2,0.5*(p+4))-gammainc(0.5*M2,0.5*(p+4)))...
            %                 -(4*M/(5*ct4))*prod5*(gammainc(0.5*Mct2,0.5*(p+5))-gammainc(0.5*M2,0.5*(p+5)))...
            %                 +(1/(6*ct4))*prod6*(gammainc(0.5*Mct2,0.5*(p+6))-gammainc(0.5*M2,0.5*(p+6)));
            
            % Expectation of rho calling function Euj
            Expectationrho=prod2*gammainc(0.5*M2,0.5*(p+2))/2 ...
                +(M2/2+ct*(5*ct+16*M)/30)*(1-gammainc(0.5*Mct2,0.5*p)) ...
                +(M2/2-M2*(M^4-5*M2*ct2+15*ct4)/(30*ct4))*(gammainc(0.5*Mct2,0.5*p)-gammainc(0.5*M2,0.5*p)) ...
                +(1/2+M^4/(2*ct4)-M2/ct2)*prod2*Euj(Mct2,M2,p+2) ...
                +(4*M/(3*ct2)-4*M^3/(3*ct4))*prod3*Euj(Mct2,M2,p+3) ...
                +(3*M2/(2*ct4)-1/(2*ct2))*prod4*Euj(Mct2,M2,p+4)...
                -(4*M/(5*ct4))*prod5*Euj(Mct2,M2,p+5)...
                +(1/(6*ct4))*prod6*Euj(Mct2,M2,p+6);
        end
        
    end


    function res=Euj(b,a,j)
        res=gammainc(0.5*b,0.5*j)-gammainc(0.5*a,0.5*j);
        % res=chi2cdf(b,j)-chi2cdf(a,j);
    end

% Remark:
% chi2cdf(x,v) = gamcdf(x,v/2,2) = gammainc(x ./ 2, v/2);
end
%FScategory:UTISTAT