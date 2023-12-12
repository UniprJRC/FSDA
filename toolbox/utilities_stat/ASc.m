function [bdp,eff] = ASc(c, v)
%ASc computes breakdown point and efficiency associated with constant c for Andrew's rho function
%
%<a href="matlab: docsearchFS('ASc')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    c :     tuning constant c. Scalar. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%    v :        number of response variables. Scalar. Number of variables of
%               the  dataset (for regression v=1)
%               UP TO NOW v=1 (JUST REGRESSION) TO DO FOR MULTIVARIATE
%               ANALYSIS
%
%  Optional input arguments:
%
% Output:
%
%     bdp      :  bdp. Scalar. Breakdown point associated to the supplied
%                 value of c
%     eff      :  eff. Scalar. Efficiency associated to the supplied
%                 value of c
%
% More About:
%
%
%
% \[
% ASrho(u)= \left\{
%    \begin{array}{cc}
%   c (1-\cos (u / c))                   &  |u/c| \leq \pi  \\
%   2c                                   &  |u/c| > \pi   \\
% \end{array}
%    \right.
%  \]
%  
%
%
% See also: TBc, HYPc, HAc, PDc
%
%
% References:
%
% Andrews, D.F., Bickel, P.J., Hampel, F.R., Huber, P.J., Rogers, W.H., and
% Tukey, J.W. (1972), "Robust Estimates of Location: Survey and Advances",
% Princeton Univ. Press, Princeton, NJ. [p. 203]
%
% Andrews, D. F. (1974). A Robust Method for Multiple Linear Regression,
% "Technometrics", V. 16, pp. 523-531, https://doi.org/10.1080/00401706.1974.10489233
%
% Riani, M., Cerioli, A. and Torti, F. (2014), On consistency factors and
% efficiency of robust S-estimators, "TEST", Vol. 23, pp. 356-387.
% http://dx.doi.org/10.1007/s11749-014-0357-7
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('ASc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % bdp associated with a particular c.
    c=5;
    bdp = ASc(c,1)
%}

%{
    % bdp and eff associated with a particular c.
    c=5;
    [bdp, eff] = ASc(c,1)
%}

%{

    %% Breakdown vs efficiency.
    %Analysis of breakdown point and asymptotic efficiency
    %at the normal distribution as a function of c in regression.
    c=1:0.01:4;
    CC=[c' zeros(length(c),1)];
    jk=0;
    for j=c
        jk=jk+1;
         [bdp,eff]=ASc(j,1);
        CC(jk,2:3)=[bdp,eff];
    end

    
    subplot(2,1,1)
    plot(c',CC(:,2))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Breakdown point','Interpreter','none')
    subplot(2,1,2)
    plot(c',CC(:,3))
    xlabel('c','Interpreter','Latex','FontSize',16)
    ylabel('Asymptotic efficiency','Interpreter','none')
%}


%% Beginning of code

if nargin<=2 || v==1

    % Convergence condition is E(\rho) = \rho(c) bdp
    %  where \rho(c) for standardized rho function is 1
    Erho =integral(@(u)(ASrho(u,c)).*normpdf(u),-c*pi,c*pi);
    c2=2*c;
    bdp=(real(Erho)+2*c2*normcdf(c*pi,'upper'))/c2;

    % bet  = \int_{-c}^c  \psi'(x) d \Phi(x)
    % alph = \int_{-c}^c  \psi^2(x) d \Phi(x)
    % empeff = bet^2/alph = 1 / [var (robust estimator of location)]


    bet= integral(@(u)(ASpsider(u,c)).*normpdf(u),-c*pi,c*pi);

    alph= integral(@(u)((ASpsi(u,c)).^2).*normpdf(u),-c*pi,c*pi);
    eff=(real(bet)^2)/real(alph);

else
    error('FSDA:ASc:Wrongv','Not yet implemented for v>1')
end


end
%FScategory:UTISTAT