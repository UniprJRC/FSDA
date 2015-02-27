function psiHYPder = HYPpsider(u, cktuning)
%HYPpsider computes derivative of psi function for hyperbolic tangent estimator
%
%<a href="matlab: docsearchFS('hyppsider')">Link to the help function</a>
%
%  Required input arguments:
%
%    u:         n x 1 vector containing residuals or Mahalanobis distances
%               for the n units of the sample
%    cktuning :  vector of length 2 or of length 5 which specifies the value of the tuning
%                constant c (scalar greater than 0 which controls the
%                robustness/efficiency of the estimator)
%                and the prefixed value k (sup of the
%                change-of-variance sensitivity) and the values of
%                parameters A, B and d
%                cktuning(1) = c
%                cktuning(2) = k = supCVC(psi,x) x \in R
%                cktuning(3)=A;
%                cktuning(4)=B;
%                cktuning(5)=d;
%                Remark: if length(cktuning)==2 values of A, B and d will be
%                computed automatically
%
%
% Function HYPpsi transforms vector u as follows
%
% HYPpsider(u)= { 1			                               |u| <= d,
%               {
%		        { 0.5*B*(1-k) * ( 1/cosh(sqrt((k - 1) * B^2/A)*(c -|u|)/2)  )^2
%		        { 	                 d <= |u| <  c,
%               {
%		        { 0,			                         |u| >= c.
%
%	It is necessary to have 0 < A < B < 2 *normcdf(c)-1- 2*c*normpdf(c) <1
%
%
%
%
% References:
%
%
% Frank R. Hampel, Peter J. Rousseeuw and Elvezio Ronchetti (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% Journal of the American Statistical Association , Vol. 76, No. 375,
% pp. 643-648 (HRR)
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('hyppsider')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{

    % Plot of derivative of hyperbolic psi function
    x=-9:0.1:9;
    ctuning=6;
    ktuning=4.5;
    psiHYPder=HYPpsider(x,[ctuning,ktuning]);
    plot(x,psiHYPder)
    xlabel('x','Interpreter','Latex')
    ylabel(' Hyperbolic $\psi''(x) $','Interpreter','Latex')

%}

%% Comparison among four derivatives of psi function
% TB, Optimal, Hampel, Hyperbolic

%{

    bdp=0.5;
    c=TBbdp(bdp,1);
    subplot(2,2,1)
    % 1st panel is Tukey biweight
    x=-4:0.01:4;
    psiTBder=TBpsider(x,c);
    plot(x,psiTBder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')
    title('Tukey biweight')


    subplot(2,2,2)
    % 2nd panel is optimal
    c=OPTbdp(bdp,1);
    c=c/3;
    % Remark: in this case it is necessary to multiply by 3.25*c^2 because the
    % optimal has been standardized in such a way that sup(rho(x))=1
    psiOPTder=(3.25*c^2)*OPTpsider(x,c);
    plot(x,psiOPTder)
    xlim([-4 4])
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x)$','Interpreter','Latex')
    title('Optimal')

    subplot(2,2,3)
    % 3rd panel is Hampel
    % Obtain bottom panel of Figure 11.10 p. 375 of
    % Hoaglin et al. (1987)
    c=HAbdp(bdp,1);
    psiHA=HApsider(x,c);
    plot(x,psiHA)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x) $','Interpreter','Latex')
    title('Hampel')

    subplot(2,2,4)
    % 4th panel is hyperbolic
    % c=HYPbdp(0.5,1);
    %ktuning=4.5;
    ktuning=4.5;
    % Precalculated values
    c = 2.010311082005501;
    A = 0.008931591866092;
    B = 0.051928487236632;
    d=  0.132017481327058;
    % Alternatively the values can be found using
    %[c,A,B,d]=HYPbdp(0.5,1,ktuning);

    psiHYPder=HYPpsider(x,[c,ktuning,A,B,d]);

    plot(x,psiHYPder)
    xlabel('x','Interpreter','Latex')
    ylabel('$\psi''(x) $','Interpreter','Latex')
    title('Hyperbolic')
%}

%% Beginning of code

c = cktuning(1);
k = cktuning(2);
if length(cktuning)>2

        A=cktuning(3);
        B=cktuning(4);
        d=cktuning(5);

    if ((A < 0) || (B < A) || (B>1)),
        error('FSDA:HYPpsider:WrongAorB',[' Illegal choice of parameters in hyperbolic tangent estimator: ' ...
            num2str(param) ]')
    else   
    end
    
else
    % Find parameters A, B and d using routine HYPck
    [A,B,d]=HYPck(c,k);
    
    % For example if c=4 and k=5
    %     A = 0.857044;
    %     B = 0.911135;
    %     d =1.803134;
    % see Table 2 of HRR
end


psiHYPder = zeros(size(u));
absu=abs(u);

%  1,		   |u| <=d
psiHYPder(absu<=d) = 1;


%                d <= |u| < c,
%  (k - 1)*(-1/2) * (cosh(sqrt((k - 1) * B^2/A)*(c -|u|)/2))^-2 
psiHYPder(absu > d & absu <=c) = 0.5*B* (1-k) *( cosh(sqrt((k - 1) * B^2/A)...
    *(c - absu(absu > d & absu <=c ))/2).^(-2)) ;

% 0,			              |u| >= c.

end
