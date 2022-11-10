function [bdp,eff] = PDc(alpha)
%PDc computes breakdown point and efficiency associated with tuning constant alpha for minimum power divergence estimator
%
%
%<a href="matlab: docsearchFS('PDc')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    alpha :     tuning constant alpha. Scalar or Vector. Scalar greater than 0 which
%               controls the robustness/efficiency of the estimator
%
%  Optional input arguments:
%
%
% Output:
%
%     bdp      :  bdp. Scalar. Breakdown point associated to the supplied
%                 value of c
%     eff      :  eff. Scalar. Efficiency associated to the supplied
%                 value of c
%                 Remark: if alpha is a vector bdp and eff will also be vectors
%                 with the same size of alpha. For example bdp(3) and eff(3)
%                 are associated to alpha(3) ....
%
% See also: TBc, HYPc, HAc, OPTc
%
% References:
%
%
%  Riani, M. Atkinson, A.C., Corbellini A. and Perrotta A. (2020), Robust
%  Regression with Density Power Divergence: Theory, Comparisons and Data
%  Analysis, submitted.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('PDc')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %PDc with just one output argument.
    [bdp]=PDc(1)
    disp('Break down point')
    disp(bdp)
%}

%{
    %PDc with 2 output arguments.
    alpha=1;
    [bdp,eff]=PDc(alpha)
    disp('Break down point and efficienty')
    disp(['alpha=' num2str(alpha)])
    disp(['bdp=' num2str(bdp)])
    disp(['eff=' num2str(eff)])
%}



%{
    %Breakdown point and efficiency.
    %Analysis of breakdown point and asymptotic efficiency
    %at the normal distribution as a function of alpha in regression.
    c=0.01:0.01:1;
    [bdp,eff]=PDc(c);
    subplot(2,1,1)
    plot(c,bdp)
    xlabel('$\alpha$','Interpreter','Latex','FontSize',16)
    ylabel('Breakdown point','Interpreter','none')
    subplot(2,1,2)
    plot(c,eff)
    xlabel('$\alpha$','Interpreter','Latex','FontSize',16)
    ylabel('Asymptotic efficiency','Interpreter','none')

%}

%% Beginning of code
bdp=1-1./sqrt(1+alpha);
eff=(sqrt(1+2*alpha)./((1+alpha))).^3;

end
%FScategory:UTISTAT