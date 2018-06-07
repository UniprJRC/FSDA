function [bdp,eff,A,B,d] = HYPc(c,v,varargin)
%HYPc computes breakdown point and efficiency associated with constant chyperbolic tangent estimator (for a given value of k=sup CVC)
%
%
%
%<a href="matlab: docsearchFS('HYPc')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    c :     tuning constant c. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%    v :        number of response variables. Scalar. Number of variables of
%               the  dataset (for regression p=1)
%               UP TO NOW v=1 (JUST REGRESSION) TO DO FOR MULTIVARIATE
%               ANALYSIS
%
%  Optional input arguments:
%
%   k        : supremum of the change of variance curve
%              $supCVC(psi,x) x \in R$.
%              Default value is k=4.5
%               Example - 'k',4 
%               Data Types - double
%   shapeeff : location or shape efficiency. Scalar. If 1, the efficiency is referred to the shape else (default)
%              is referred to the location. TODO:Hac:shapeeff  
%               Example - 'shapeeff',1 
%               Data Types - double
%      param : parameters A, B and d. Vector with 3 elements.
%              Vector of length 3 specifying the parameters A, B and d of the
%              weight function of the hyperbolic tangent estimator.
%              param(1)=A param(2)=B param(3)=d
%              If these values are not supplied they will be automatically
%              calculated calling routine HYPck
%               Example - 'param',[0.0035;0.0317;0.0827] 
%               Data Types - double
%
% Output:
%
%     bdp      :  bdp. Scalar. Breakdown point associated to the supplied
%                 value of c
%     eff      :  eff. Scalar. Efficiency associated to the supplied
%                 value of c
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
% See also TBc, HAc, OPTc
%
% References:
%
%
% Hampel, F.R., Rousseeuw, P.J. and  Ronchetti E. (1981),
% The Change-of-Variance Curve and Optimal Redescending M-Estimators,
% "Journal of the American Statistical Association", Vol. 76,
% pp. 643-648 [HRR]

% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HYPc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %Analysis of breakdown point and asymptotic efficiency
    %at the normal distribution as a function of c in regression.
    cc=2:0.1:6;
    % BDPEFF = matrix which will contain
    % 1st column = value of c
    % 2nd column = breakdown point (bdp)
    % 3rd column = asympotic nominal efficiency (eff)
    % 4th column = value of parameter A
    % 5th column = value of parameter B
    % 6th column = value of parameter d
    BDPEFF=[cc' zeros(length(cc),5)];

    % Fixed value of k (sup CVC) which must be used
    kdef=4.5;

    jk=1;
    for c=cc
        [bdp,eff,A,B,d]=HYPc(c,1,'k',kdef);
        BDPEFF(jk,2:end)=[bdp, eff, A, B, d];
        jk=jk+1;
    end


    nr=2;
    nc=2;
    subplot(nr,nc,1)
    plot(BDPEFF(:,1),BDPEFF(:,2))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Breakdown point','Interpreter','none')

    subplot(nr,nc,2)
    plot(BDPEFF(:,1),BDPEFF(:,3))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Asymptotic efficiency','Interpreter','none')

    subplot(nr,nc,3)
    plot(BDPEFF(:,1),BDPEFF(:,4:5))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('A and B','Interpreter','none')

    subplot(nr,nc,4)
    plot(BDPEFF(:,1),BDPEFF(:,6))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('d','Interpreter','none')

    suplabel(['Constant k=' num2str(kdef)],'t');

%}

%% Beginning of code
options=struct('k',4.5,'shapeeff',0,'param','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:HYPc:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<2
    error('FSDA:HYPc:missingInputs','Initial value of p is missing');
end

if nargin > 3
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

k=options.k;

param=options.param;

if isempty(param)
    % Find value of A, B and d given c and k
    [A,B,d]=HYPck(c,k);
else
    A=param(1);
    B=param(2);
    d=param(3);
end

% Now compute the expectation of the rho function

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

% Convergence condition is E(\rho) = \rho(c) bdp
%  where \rho(c) for HYP is rhoc

bdp=Erho/rhoc;

eff=B^2/A;

end

%FScategory:UTISTAT
