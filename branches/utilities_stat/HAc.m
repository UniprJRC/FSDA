function [bdp,eff] = HAc(ctun,v,varargin)
%HAc computes breakdown point and efficiency associated with constant c 
%
%<a href="matlab: docsearchFS('HAc')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    ctun :     tuning constant c. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%    v :        number of response variables. Scalar. Number of variables of
%               the  dataset (for regression v=1)
%               UP TO NOW v=1 (JUST REGRESSION) TO DO FOR MULTIVARIATE
%               ANALYSIS
%
%  Optional input arguments:
%
%      param : tuning parameters. Vector. Vector of length 3 specifying the parameters a, b and c of the
%              weight function of the Hampel estimator.
%              param(1)=a param(2)=b param(3)=c
%              If these values are not supplied they will be automatically
%              set to a=2, b=4 c=8
%               Example - 'param',[1.5,3.5,8] 
%               Data Types - double
%   shapeeff : location or shape efficiency. Scalar. If 1, the efficiency is referred to the shape else (default)
%              is referred to the location. TODO:Hac:shapeeff  
%               Example - 'shapeeff',1 
%               Data Types - double
%
% Output:
%
%     bdp      :  bdp. Scalar. Breakdown point associated to the supplied
%                 value of c for Hampel rho function 
%     eff      :  eff. Scalar. Efficiency associated to the supplied
%                 value of c for Hampel rho function 
%
% More About:
%
% Function HApsi transforms vector u as follows.
%  \[
%  HApsi(u)  = \left\{   
%  \begin{array}{cc}
%    u & |u| <= a                                       \\
%    a \times sign(u) & a <= |u| < b                    \\
%    a \frac{c-|u|}{c-b} \times sign(u) & b <= |u| <  c \\
%    0 & |u| >= c 
%  \end{array} \right.
% \]
%
%             where $a$= ctun *param(1).
%                   $b$= ctun *param(2).
%                   $c$= ctun *param(3).
%
%             The default is
%                   $a$= 2*ctun. 
%                   $b$= 4*ctun. 
%                   $c$= 8*ctun. 
%
%	It is necessary to have 0 <= a <= b <= c
%
% See also: HYPc, TBc, OPTc
%
% References:
%
% Hoaglin, D.C., Mosteller, F., Tukey, J.W. (1982), "Understanding Robust and
% Exploratory Data Analysis", Wiley, New York.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('HAc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %% bdp and eff as function of c.
    %Analysis of breakdown point and asymptotic efficiency
    %at the normal distribution as a function of c in regression.
    cc=0.15:0.05:1.2;
    % BDPEFF = matrix which will contain
    % 1st column = value of c
    % 2nd column = breakdown point (bdp)
    % 3rd column = asympotic nominal efficiency (eff)
    BDPEFF=[cc' zeros(length(cc),2)];

    jk=1;
    for c=cc
        [bdp,eff]=HAc(c,1);
        BDPEFF(jk,2:end)=[bdp, eff];
        jk=jk+1;
    end


    nr=2;
    nc=1;
    subplot(nr,nc,1)
    plot(BDPEFF(:,1),BDPEFF(:,2))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Breakdown point','Interpreter','none')

    subplot(nr,nc,2)
    plot(BDPEFF(:,1),BDPEFF(:,3))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Asymptotic efficiency','Interpreter','none')

%}

%% Beginning of code
abcdef=[2 4 8];
options=struct('shapeeff',0,'param',abcdef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:HAc:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<2
        error('FSDA:HAc:missingInputs','Initial value of p is missing');
end

if nargin > 3
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end


param=options.param;

% Now compute the expectation of the rho function


a=(param(1)*ctun);
b=(param(2)*ctun);
c=(param(3)*ctun);

a2=a.^2/2;
b2=b.^2/2;
c2=c.^2/2;


phic=a*b-0.5*a^2+0.5*(c-b)*a;

% |u| <a
% Erhoa=  \int_-a^a u^2/2
Erhoa=0.5*v*gammainc(a2,(v+2)/2);
% rhoa = @(u,a,b,c)(0.5*u.^2.*(1/sqrt(2*pi)).*exp(-0.5*u.^2));
% Erhoack=integral(@(u)rhoa(u,a,b,c),-a,a);

% a< |u| <b
% Erhoab= 2 * \int_a^b (au-0.5a^2)
Erhoab=2*a*(normpdf(a)-normpdf(b))-a^2*(normcdf(b)-normcdf(a));
% rhoab = @(u,a,b,c)((a*u-0.5*a^2).*(1/sqrt(2*pi)).*exp(-0.5*u.^2));
% Erhoabck=2*integral(@(u)rhoab(u,a,b,c),a,b);

% b< |u| <c
% Erhobc = \int_b^c \rho(x) \Phi(x)
Erhobc1=2*(a*b-0.5*a^2+0.5*(c-b)*a*(1 -c^2/((c-b)^2)))*(normcdf(c)-normcdf(b));
Erhobc2=0.5*a*v*(gammainc(b2,(v+2)/2) -gammainc(c2,(v+2)/2)) /(c-b);
Erhobc3=2*a*c*(normpdf(b)-normpdf(c))/(c-b);
Erhobc=Erhobc1+Erhobc2+Erhobc3;



% |u| >c
Erhoc=phic*( 1-gammainc(c2,v/2) );


Erho= Erhoa+Erhoab+Erhobc+Erhoc;

% Convergence condition is E(\rho) = \rho(c) bdp
%  where \rho(c) for HYP is rhoc
bdp=Erho/phic;

% bet  = \int  \psi'(x) d \Phi(x)
% bet = \int_-a^a d \Phi(x) +2* \int_b^c -a/(c-b)
bet= gammainc(a2,v/2)+(gammainc(b2,v/2)-gammainc(c2,v/2))*a/(c-b);

% alph = \int \psi^2(x) d \Phi(x)
alph= v*gammainc(a2,(v+2)/2)...                                        % 2* \int_0^a x^2 f(x) dx
    +a.^2 .*(gammainc(b2,v/2)-gammainc(a2,v/2))...                     % 2* a^2 \int_a^b f(x) dx
    +(a./(c-b)).^2 .*(c.^2.*(gammainc(c2,v/2)-gammainc(b2,v/2)) ...    %(a./(c-b)).^2 (2 c^2 \int_b^c f(x) dx
    + v*(gammainc(c2,(v+2)/2)-gammainc(b2,(v+2)/2)) ...                %   + 2*  \int_b^c x^2 f(x) dx
    -2*c*v*sqrt(2/pi)*(gammainc(c2,(v+1)/2)-gammainc(b2,(v+1)/2)));        % +2 *2* \int_b^c |x| f(x)


% Remark: if v=1
% -2*c*v*sqrt(2/pi)*(gammainc(c2,(v+1)/2)-gammainc(b2,(v+1)/2)));
%     -4*c.*(normpdf(b)-normpdf(c))  );

% empeff = bet^2/alph = 1 / [var (robust estimator of location)]
eff=(bet^2)/alph;

end

%FScategory:UTISTAT
