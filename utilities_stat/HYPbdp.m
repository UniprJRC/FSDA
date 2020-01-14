function [c,A,B,d] = HYPbdp(bdp, v,k,traceiter)
%HYPbdp finds constant c which is associated to the requested breakdown point for hyp. tan. estimator
%
%<a href="matlab: docsearchFS('HYPbdp')">Link to the help page for this function</a>
%
%  Required input arguments:
%
%      bdp    : requested breakdown point. Scalar.
%               Scalar defining breakdown point (i.e a number between 0 and
%               0.5)
%               Data Types - single|double
%        v    : number of response variables. Scalar.
%               (e.g. in regression p=1) 
%               Data Types - single|double|int32|int64
%
%
%
%  Optional input arguments:
%
%   k        : supremum of the change of variance curve. Scalar.
%              $\sup CVC(psi,x) x \in R$.
%              Default value is k=4.5.
%                 Example - 'k',5
%                 Data Types - double
%  traceiter : Level of display. Scalar.
%               If traceiter = 1 it is possible to monitor
%              how the value of the objective function E(rho)/\rho(\infty)
%              gets closer to the target (bdp) during the iterations
%               Example - 'traceiter',0
%               Data Types - single|double|int32|int64
% Output:
%
%  c    : parameter c of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%   A   : parameter A of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%   B   : parameter B of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%   d   : parameter d of hyperbolic tangent estimator. Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%
% More About:
%
%  \[
%   HYPpsi(u) =
% \left\{
%   \begin{array}{cc}
%  	 u &        |u| \leq  d \\
%                  \sqrt{A (k - 1)}  \tanh \left( \sqrt{(k - 1) B^2/A} (c -|u|)/2 \right) sign(u) &
% 		         	                 d \leq |u| <  c, \\
%                 0 &                      |u| \geq c.
% \end{array}
%    \right.
%  \]
%  	It is necessary to have $0 < A < B < 2 normcdf(c)-1- 2 c \times normpdf(c) <1$
%
%
% See also TBbdp, HAbdp, OPTbdp
%
%
% References:
%
%
% Hampel, F.R., Rousseeuw, P.J. and  Ronchetti E. (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% "Journal of the American Statistical Association", Vol. 76,
% pp. 643-648. [HRR]
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('HYPbdp')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:

%{
    % Find value of c, A, B, for fixed break down point. 
    % Find required values when k=4.5 and bdp=0.5.
 
    ktuning=4.5;
    [c,A,B,d]=HYPbdp(0.5,1,ktuning);
    % In this case
    % c = 2.010311082005501
    % A = 0.008931591866092
    % B = 0.051928487236632
    % d=  0.132017481327058
%}

%{
    % Efficienty monitoring.
    % Analysis of efficiency and of paramters A, B abd k as function of bdp
    % for a given value of sup CVC=4
    seqi=0.1:0.1:0.5;
    eff=[seqi' zeros(length(seqi),4)];

    iter=0;
    k=4;

    for i=seqi
        [c,A,B,d] = HYPbdp(i,1,k);
        iter=iter+1;
        eff(iter,2:5)=[B^2/A A B d];
    end

    subplot(2,2,1)
    plot(eff(:,1),eff(:,2))
    title('efficiency')


    subplot(2,2,2)
    plot(eff(:,1),eff(:,3))
    title('A')

    subplot(2,2,3)
    plot(eff(:,1),eff(:,4))
    title('B')

    subplot(2,2,4)
    plot(eff(:,1),eff(:,5))
    title('d')

%}

%% Beginning of code


if (nargin >2)
    if (k < 0) 
        error('FSDA:HYPbdp:WrongK',[' Illegal choice of parameters in hyperbolic tangent estimator: ' ...
            num2str(k) ]')
    end
else
    k=4.5;
end

if nargin <4
    traceiter =0;
end


% c = starting point of the iteration
% Note that the procedure is not sensitive at all to the starting point of
% the iteration
c=2.3;

% step = width of the dichotomic search (it decreases by half at each
% iteration). Generally it can be smaller. A large value ensures converge
% when bdp is very small and p is very large.
step=3;

% Convergence condition is E(\rho) = \rho(c) bdp
%  where \rho(c) 
Erho1=10;

if v==1
    eps=1e-7;
    iter=0;
    while abs(Erho1-1)>eps
        % Find value of A, B and d given c and k
        [A,B,d]=HYPck(c,k);
        
        iter=iter+1;
        
        if iter==80
            disp(['Effective tolerance in routine HYPbdp=' num2str(abs(Erho1-1))])
            break
        end
        
        c2=c.^2/2;
        Erhoa= v*gammainc(d.^2/2,0.5*(v+2))/2;
        
        % Erhoa is also equal to
        % tmpu=@(u) (u.^2) .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
        % integral(@(u)tmpu(u),-d,d)/2
        
        % Rho function inside interval d----c
        rhodc = @(u,c,A,B,k) -2*(A/B) * log(cosh(0.5*sqrt((k - 1) * B^2/A)...
            *(c - u))) .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
        
        Erhob= 2*integral(@(u)rhodc(u,c,A,B,k),d,c)...
            +(d^2/2 + 2*(A/B)*log(cosh(0.5*sqrt((k - 1) * B^2/A)*(c -d))))*(gammainc(c2,0.5*v)-gammainc(d.^2/2,0.5*v));
        
        rhoc=d^2/2 +2*(A/B)*log(cosh(0.5*sqrt((k - 1) * B^2/A)*(c -d)));
        Erhoc=rhoc*(1-gammainc(c2,0.5*v));
        
        % Eho = E [ rho]
        Erho= Erhoa+Erhob+Erhoc;
        
        Erho1=Erho/(rhoc*bdp);
        
        if traceiter==1
            disp([iter c Erho/rhoc])
        end
        
        step=step/2;
        if Erho1>1
            c=c+step;
        else
            % 1.7 is the value which usually guarantes a breakdown point greater
            % than 0.5
            c=max(c-step,1.7);
        end
        % disp([step c Erho1 iter])
    end
else
    error('FSDA:HYPeff:Wrongv','Not yet implemented for v>1')
end

end
%FScategory:UTISTAT