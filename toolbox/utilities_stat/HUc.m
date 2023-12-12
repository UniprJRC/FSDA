function [bdp,eff] = HUc(c,v)
%HUc computes breakdown point and efficiency associated with constant c for Huber link
%
%<a href="matlab: docsearchFS('HUc')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    c :     tuning constant c. Scalar or vector.  greater than 0 which
%               controls the robustness/efficiency of the estimator
%    v :        number of response variables. Scalar. Number of variables of
%               the  dataset (for regression v=1)
%               Now it is implemented just for v=1
%
%  Optional input arguments:
%
%
% Output:
%
%     bdp      :  breakdown point. Breakdown point associated to the supplied
%                 value(s) of c. bdp has the same dimension of c.
%     eff      :  asymptocic efficiency. Efficiency associated to the supplied
%                 value(s) of c. eff has the same dimension of c.
%
%
% See also: TBc, HYPc, HAc, OPTc
%
% References:
%
% Huber, P.J. (1981), "Robust Statistics", Wiley.
%
% Huber, P.J. and Ronchetti, E.M. (2009), "Robust Statistics, 2nd Edition",
% Wiley.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('HUc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %HUc with just one output argument.
    [bdp]=HUc(2,1);
    disp('Break down point')
    disp(bdp)
%}

%{
    %HUc with 2 output arguments.
    [bdp,eff]=HUc(2,1)
    disp('Break down point and efficiency')
    disp(bdp)
    disp(eff)
%}



%{
    %Breakdown point and efficiency.
    %Analysis of breakdown point and asymptotic efficiency
    %at the normal distribution as a function of c in regression.
    c=0:0.01:6;
    [bdp,eff]=HUc(c,1);
    subplot(2,1,1)
    plot(c,bdp,'LineWidth',2)
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Breakdown point','Interpreter','none')
    subplot(2,1,2)
    plot(c,eff,'LineWidth',2)
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Asymptotic efficiency','Interpreter','none')

%}

%% Beginning of code


% bet  = \int_{-c}^c  \psi'(x) d \Phi(x)
% alph = \int_{-c}^c  \psi^2(x) d \Phi(x)
% eff = bet^2/alph = 1 / [var (robust estimator of location)]

bet= normcdf(c)-normcdf(-c);

alph= 2*(c.^2.*(1-normcdf(c))+normcdf(c)-0.5-c.*normpdf(c));
eff=(bet.^2)./alph;

varmedianinv=(2*normpdf(0))^2;
eff(c==0)=varmedianinv;
bdp=zeros(size(c));
bdp(c==0)=0.5;
end
%FScategory:UTISTAT