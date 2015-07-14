function [Anew,Bnew,d]=HYPck(c,k,A,B,d)
%HYPck computes values of the scalars A, B, d for hyperbolic tangent estimator
%
%
%<a href="matlab: docsearchFS('hypck')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    c :     tuning constant c. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%    k :     sup of change of variance curve (CVC). Scalar. $k= supCVC(psi,x) x \in R$
%
%  Optional input arguments:
%
%         A   : A parameter. Scalar. Starting value for parameter A
%         B   : B parameter. Scalar. Starting value for parameter B
%         d   : d parameter. scalar. Starting value for parameter d
%
%
% Output:
%
%  Anew : Value of parameter A.  Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%  Bnew : Value of parameter B. Scalar. 
%         For more details see the  methodological details inside "More
%         About" below
%  d    : Value of parameter d.  Scalar.
%         For more details see the  methodological details inside "More
%         About" below
%
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
% See also HYPc
%
% References:
%
%
% Hampel,F.R.,  Rousseeuw P.J. and  Ronchetti E.(1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% Journal of the American Statistical Association , Vol. 76, No. 375,
% pp. 643-648 (HRR)
%
%
%<a href="matlab: docsearchFS('hypck')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
% Examples:
%
%{
    % Reconstuct columns 3:6 of Table 2 of HRR
    cc=3:6;
    kk=4:0.5:5;

    ABD=zeros(length(cc)*length(kk),4);
    ij=1;

    for c=cc
        for k=kk
            [A,B,d]=HYPck(c,k);
            eff=B^2/A;
            ABD(ij,:)=[A,B,d,eff];
            ij=ij+1;
        end
    end
%}
%

% Remark: this function is called by HYPrho, HYPpsi, HYPwei

%% Beginning of code

if nargin==2
    % Set initial values for A, B and d
    d=0.1;
    B=(2*normcdf(c)-1-2*c*normpdf(c))/2;
    A=B*(3/4);
end
dd=1;

Anew=0.5;
Bnew=0.5;
iter=0;
while max(abs([A-Anew,B-Bnew,dd-d]))>1e-8
    iter=iter+1;
    %disp(iter)
    if iter==1000
        disp(['Effective tolerance in routine HYPck=' num2str(max(abs([A-Anew,B-Bnew,dd-d])))])
        break
    end
    
    
    % disp([A-Anew B-Bnew dd-d max(abs([A-Anew,B-Bnew,dd-d]))])
    if iter>1
        A=Anew;
        B=Bnew;
    end
    
    % In the loop below we find the value of d
    % keeping into account the constraint
    
    % Check if the solution is feasible
    dchx=(1e-15:0.1:10)';
    dchy=sqrt(A*(k - 1))*tanh(0.5*sqrt((k - 1) * B^2/A)*(c - dchx));
    difdd=dchx-dchy;
    if min(difdd)<0 && max(difdd)>0
        % Find best starting value of d
        d=dchx(find(difdd<0,1,'last'));
        dinit=d;
    else
        
        error('FSDA:HYPck:NoConvergence','Solution for k is impossible: increase the value of k');
        %         B=(2*normcdf(c)-1-2*c*normpdf(c))*rand(1,1);
        %         A=B*rand(1,1);
    end
    
    step=0.8;
    while abs(dd-d)>1e-8
        
        dd=sqrt(A*(k - 1))*tanh(0.5*sqrt((k - 1) * B^2/A)*(c - d));
        % disp(dd)
        step=step/2;
        
        
        
        if dd/d>1
            d=d+step;
        else
            ep=dinit;
            d=max(d-step,ep);
            % d=(0:0.1:6)';
            % dd= sqrt(A*(k - 1))*tanh(0.5*sqrt((k - 1) * B^2/A)*(c - d));
            % plot(d,dd)
        end
        
    end
    
    % disp(d)
    % Find new value of A = \int psi^2
    psi2 = @(u,c,A,B,k) (A * (k - 1)) * (tanh(sqrt((k - 1) * B^2/A)...
        *(c - u)/2)).^2 .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
    
    Anew = 2*integral(@(u)psi2(u,c,A,B,k),d,c)+gammainc(d^2/2,3/2);
    
    % Derivative of psi
    psider = @(u,c,A,B,k) (k - 1)*B*(-0.5)* (cosh(sqrt((k - 1) * B^2/A)...
        *(c - u)/2)).^(-2) .*(1/sqrt(2*pi)).*exp(-0.5*u.^2);
    
    
    % Find new value of B = \int psi'
    Bnew = 2*integral(@(u)psider(u,c,A,B,k),d,c)+gammainc(d^2/2,0.5);
    
    if Anew>Bnew
        % This a pathological situation which may happen when the value of
        % c is small (<2) therefore bpd >0.5 A and B and very small (close
        % to zero
        Bnew=abs(Bnew);
        Anew=Bnew*3/4;
    end
    
    dd=sqrt(Anew*(k - 1))*tanh(0.5*sqrt((k - 1) * Bnew^2/Anew)*(c - d));
    
end


end
