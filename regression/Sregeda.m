function [out , varargout] = Sregeda(y,X,varargin)
%Sregeda computes S estimators in linear regression for a series of values of bdp
%
%<a href="matlab: docsearchFS('Sregeda')">Link to the help function</a>
%
%  Required input arguments:
%
%    y: Response variable. Vector. A vector with n elements that contains
%       the response variable. y can be either a row or a column vector.
%    X: Data matrix of explanatory variables (also called 'regressors') of
%       dimension (n x p-1). Rows of X represent observations, and columns
%       represent variables.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       automatically be excluded from the computations.
%
%  Optional input arguments:
%
%         bdp :  breakdown point. Scalar or vector.
%               It measures the fraction of outliers
%               the algorithm should resist. In this case any value greater
%               than 0 but smaller or equal than 0.5 will do fine.
%               The default for bdp is a sequence from 0.5 to 0.01
%               with step -0.01. The sequence is forced to be monotonically
%               decreasing.
%                 Example - 'bdp',[0.5 0.4 0.3 0.2 0.1]
%                 Data Types - double
%
%
%      bestr  : number of "best betas" to remember. Scalar. Scalar defining
%               number of "best betas" to remember from the subsamples.
%               These will be later iterated until convergence (default=5).
%                 Example - 'bestr',10
%                 Data Types - single | double
%
%
%     conflev :  Confidence level which is
%               used to declare units as outliers. Scalar.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%                 Example - 'conflev',0.99
%                 Data Types - double
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%
%     minsctol: tolerance for the iterative
%               procedure for finding the minimum value of the scale. Scalar.
%               Value of tolerance for the iterative
%               procedure for finding the minimum value of the scale
%               for each subset and each of the best subsets
%               (It is used by subroutine minscale.m)
%               The default value is 1e-7;
%                 Example - 'minsctol',1e-7
%                 Data Types - single | double
%
%
%        msg  : Level of output to display. Scalar. It controls whether
%                 to display or not messages on the screen.
%               If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the estimator
%               and the warnings about
%               'MATLAB:rankDeficientMatrix', 'MATLAB:singularMatrix' and
%               'MATLAB:nearlySingularMatrix' are set to off
%               else no message is displayed on the screen
%                 Example - 'msg',0
%                 Data Types - single | double
%
%       nocheck : Check input arguments. Boolean. If nocheck is equal to
%               true no check is performed on matrix y and matrix X. Notice
%               that y and X are left unchanged. In other words the
%               additional column of ones for the intercept is not added.
%               As default nocheck=false.
%               Example - 'nocheck',true
%               Data Types - boolean
%
%       nsamp   : Number of subsamples which will be extracted to find the
%                 robust estimator. Scalar. If nsamp=0 all subsets will be extracted.
%                 They will be (n choose p).
%                 If the number of all possible subset is <1000 the
%                 default is to extract all subsets otherwise just 1000.
%                 Example - 'nsamp',1000
%                 Data Types - single | double
%
%
%    refsteps : Number of refining iterations. Scalar. Number of refining
%               iterationsin each subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'refsteps',10
%                 Data Types - single | double
%
%refstepsbestr: number of refining iterations for each best subset. Scalar.
%               Scalar defining number of refining iterations for each
%               best subset (default = 50).
%                 Example - 'refstepsbestr',10
%                 Data Types - single | double
%
%     reftol  : tolerance for the refining steps. Scalar.
%               The default value is 1e-6;
%                 Example - 'reftol',1e-05
%                 Data Types - single | double
%
% reftolbestr : Tolerance for the refining steps. Scalar.
%               Tolerance for the refining steps
%               for each of the best subsets
%               The default value is 1e-8;
%                 Example - 'reftolbestr',1e-10
%                 Data Types - single | double
%
%
%     rhofunc : rho function. String. String which specifies the rho
%               function which must be used to weight the residuals.
%               Possible values are
%               'bisquare'
%               'optimal'
%               'hyperbolic'
%               'hampel'
%               'mdpd'.
%               'bisquare' uses Tukey's $\rho$ and $\psi$ functions.
%               See TBrho.m and TBpsi.m.
%               'optimal' uses optimal $\rho$ and $\psi$ functions.
%               See OPTrho.m and OPTpsi.m.
%               'hyperbolic' uses hyperbolic $\rho$ and $\psi$ functions.
%               See HYPrho.m and HYPpsi.m.
%               'hampel' uses Hampel $\rho$ and $\psi$ functions.
%               See HArho.m and HApsi.m.
%               'mdpd' uses Minimum Density Power Divergence $\rho$ and $\psi$ functions.
%               See PDrho.m and PDpsi.m.
%               The default is bisquare
%                 Example - 'rhofunc','optimal'
%                 Data Types - character
%
% rhofuncparam: Additional parameters for the specified rho function.
%               Scalar or vector.
%               For hyperbolic rho function it is possible to set up the
%               value of k = sup CVC (the default value of k is 4.5).
%               For Hampel rho function it is possible to define parameters
%               a, b and c (the default values are a=2, b=4, c=8)
%                 Example - 'rhofuncparam',5
%                 Data Types - single | double
%
%       plots : Plot on the screen. Scalar.
%               If plots = 1, generates a plot with the robust residuals
%               for each value of bdp. The confidence level used to draw the
%               confidence bands for the residuals is given by the input
%               option conflev. If conflev is not specified a nominal 0.975
%               confidence interval will be used.
%                 Example - 'plots',0
%                 Data Types - single | double
%
%  Output:
%
%  out :     A structure containing the following fields
%
%            out.Beta = matrix containing the S estimator of regression
%                       coefficients for each value of bdp
%            out.Scale= vector containing the estimate of the scale
%                       (sigma) for each value of bdp. This is the value of the objective function
%              out.BS = p x 1 vector containing the units forming best subset
%                       associated with S estimate of regression coefficient.
%              out.RES= n x length(bdp) matrix containing the robust
%                       scaled residuals for each value of bdp
%         out.Weights = n x length(bdp) vector containing the estimates of
%                       the weights for each value of bdp
%        out.Outliers = Boolean matrix containing the list of
%                       the units declared as outliers for each value of bdp using confidence
%                       level specified in input scalar conflev
%         out.conflev = confidence level which has been used to declare
%                       outliers.
%         out.Singsub = Number of subsets wihtout full rank. Notice that
%                       out.singsub(bdp(jj)) > 0.1*(number of subsamples) produces a
%                       warning
%         out.rhofunc = string identifying the rho function which has been
%                       used
%    out.rhofuncparam = vector which contains the additional parameters
%                       for the specified rho function which have been
%                       used. For hyperbolic rho function the value of
%                       k =sup CVC. For Hampel rho function the parameters
%                       a, b and c.  This field is present only if input
%                       argument 'rhofunc' is  'hyperbolic' or 'hampel'.
%           out.bdp   = vector which contains the values of bdp which have
%                       been used
%            out.y    = response vector y.
%            out.X    = data matrix X.
%           out.class = 'Sregeda'
%
%  Optional Output:
%
%            C        : matrix containing the indices of the subsamples
%                       extracted for computing the estimate (the so called
%                       elemental sets).
%
%
% See also: Sreg, MMreg, Taureg
%
% References:
%
% Riani, M., Cerioli, A., Atkinson, A.C. and Perrotta, D. (2014), Monitoring
% Robust Regression, "Electronic Journal of Statistics", Vol. 8, pp. 646-677.
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Sregeda')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    %% Sregeda with msg=0.
    % Run this code to see the output shown in the help file
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Sregeda(ycont,X,'msg',0);
    resfwdplot(out)
    ylabel('Scaled S residuals');
%}

%{
    % Sreg with optional input arguments.
    % Sreg with optimal rho function. Run this code to see the output shown in the help file
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Sregeda(ycont,X,'rhofunc','optimal');
%}


%{
    % Sreg with hyperbolic rho function.
    % Run this code to see the output shown in the help file
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Sregeda(ycont,X,'rhofunc','hyperbolic');
%}

%{
    %% Sreg on Stars data.
    % Run this code to see the Figure 2 of the article in the References
    load('stars');
    X=stars{:,1};
    y=stars{:,2};
    [out]=Sregeda(y,X,'rhofunc','bisquare');

    standard.Color={'b'}
    standard.xvalues=size(out.RES,1)-size(out.RES,2)+1:size(out.RES,1)
    fground.Color={'r'}
    resfwdplot(out,'standard',standard,'fground',fground)
    ylabel('Scaled S residuals');

    
    RHO = [];
    for i=1:49
        RHO(i,1) = corr(out.RES(:,i),out.RES(:,i+1),'type','Spearman');
        RHO(i,2) = corr(out.RES(:,i),out.RES(:,i+1),'type','Kendall');
        RHO(i,3) = corr(out.RES(:,i),out.RES(:,i+1),'type','Pearson');
    end
    minc = min(RHO);
    maxc = max(RHO);
    ylimits = [min(minc)*0.8,max(maxc)*1.1];
    figure;
    subplot(3,1,1);
    plot(out.bdp(1:49),RHO(:,1)');
    if strcmp(out.class,'Sregeda')
        set(gca,'XDir','reverse','ylim',ylimits);
        title('Spearman');
    end

    subplot(3,1,2);
    plot(out.bdp(1:49),RHO(:,2)');
    if strcmp(out.class,'Sregeda')
        set(gca,'XDir','reverse','ylim',ylimits);
        title('Kendall');
    end

    subplot(3,1,3);
    plot(out.bdp(1:49),RHO(:,3)');
    if strcmp(out.class,'Sregeda')
        set(gca,'XDir','reverse','ylim',ylimits);
        title('Pearson');
    end

%}

%% Beginning of code


nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% default value of break down point
bdpdef=0.5:-0.01:0.01;


% default values of subsamples to extract
ncomb=bc(n,p);
nsampdef=min(1000,ncomb);

% default value of number of refining iterations (C steps) for each extracted subset
refstepsdef=3;
% default value of tolerance for the refining steps convergence for  each extracted subset
reftoldef=1e-6;
% default value of number of best betas to remember
bestrdef=5;
% default value of number of refining iterations (C steps) for best subsets
refstepsbestrdef=50;
% default value of tolerance for the refining steps convergence for best subsets
reftolbestrdef=1e-8;
% default value of tolerance for finding the minimum value of the scale
% both for each extracted subset and each of the best subsets
minsctoldef=1e-7;

% rho (psi) function which has to be used to weight the residuals
rhofuncdef='bisquare';
%rhofuncdef='optimal';

if coder.target('MATLAB')
    
    % store default values in the structure options
    options=struct('intercept',true,'nsamp',nsampdef,'refsteps',refstepsdef,...
        'reftol',reftoldef,'refstepsbestr',refstepsbestrdef,'reftolbestr',reftolbestrdef,...
        'minsctol',minsctoldef,'bestr',bestrdef,'rhofunc',rhofuncdef,'rhofuncparam','','bdp',bdpdef,...
        'plots',0,'conflev',0.975,'nocheck',false,'msg',1);
    
    % check user options and update structure options
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:Sregeda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end


if coder.target('MATLAB')
    
    % Get user values of warnings
    warnrank=warning('query','MATLAB:rankDeficientMatrix');
    warnsing=warning('query','MATLAB:singularMatrix');
    warnnear=warning('query','MATLAB:nearlySingularMatrix');
    % Set them to off inside this function
    % At the end of the file they will be restored to previous values
    warning('off','MATLAB:rankDeficientMatrix');
    warning('off','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
    
end

bdp = options.bdp;              % break down point
refsteps = options.refsteps;    % refining steps
reftol = options.reftol;        % tolerance for refining steps
bestr = options.bestr;          % best betas for refining steps till convergence
nsamp = options.nsamp;          % subsamples to extract
minsctol = options.minsctol;    % tolerance for finding minimum value of the scale for each subset
refstepsbestr=options.refstepsbestr;  % refining steps for the best subsets
reftolbestr=options.reftolbestr;      % tolerance for refining steps for the best subsets
msg=options.msg;                % Scalar which controls the messages displayed on the screen
rhofunc=options.rhofunc;        % String which specifies the function to use to weight the residuals

if min(bdp)<0
    error('FSDA:Sregeda:Wrongbdp','elements of bdp must lie in the interval [0 0.5]')
elseif max(bdp)>0.5
    error('FSDA:Sregeda:Wrongbdp','elements of bdp must lie in the interval [0 0.5]')
else
    bdp=sort(bdp,'descend');
end

%% Extract in the rows of matrix C the indexes of all required subsets
[C,nselected] = subsets(nsamp,n,p,ncomb,msg);
% Store the indices in varargout
if nargout==2
    varargout={C};
end


% Store in output structure the outliers found with confidence level conflev
conflev = options.conflev;
out.conflev = conflev;

conflev = (conflev+1)/2;
seq = 1:n;
psifunc=struct;

% Define matrices which will store relevant quantities
lbdp=length(bdp);
% Beta= matrix which will contain beta coefficients
Beta=zeros(p,lbdp);
% Scale = vector which will contain the estimate of the scale
Scale=zeros(lbdp,1);
BS=zeros(p,lbdp);
Weights=zeros(n,lbdp);
Residuals=zeros(n,lbdp);
Singsub=zeros(lbdp,1);
Outliers=false(n,lbdp);
rhofuncparam=[];

for jj=1:length(bdp)
    
    if strcmp(rhofunc,'bisquare')
        % Tukey's biweight is strictly increasing on [0 c] and constant (equal to c^2/6) on [c \infty)
        % E(\rho) = kc = (c^2/6)*bdp = TBrho(c,c)*bdp, being kc the K of Rousseeuw
        % and Leroy (1987)
        
        % Compute tuning constant associated to the requested breakdown
        % point
        % For bdp =0.5 and Tukey biweight rho function c1=0.4046
        % Remark: given that in function OPTbdp rho function is defined in the interval 0---2c/3, 2c/3---3c/3, >3c/3
        % it is necessary to divide the output of OPTbdp by 3
        c=TBbdp(bdp(jj),1);
        % kc1 = E(rho) = sup(rho)*bdp
        kc=TBrho(c,c)*bdp(jj);
        
        
        psifunc.c1=c;
        psifunc.kc1=kc;
        psifunc.class='TB';
        
    elseif strcmp(rhofunc,'optimal')
        % Optimal rho function is strictly increasing on [0 3c] and constant (equal to 3.25c^2) on [3c \infty)
        % E(\rho) = kc = (3.25c^2)*bdp = TBrho(3*c,c)*bdp, being kc the K of
        % Rousseeuw and Leroy (1987)
        
        % Compute tuning constant associated to the requested breakdown
        % point
        % For bdp =0.5 and optimal rho function c= 1.2139
        c=OPTbdp(bdp(jj),1);
        % kc1 = E(rho) = sup(rho)*bdp
        kc=OPTrho(c,c)*bdp(jj);
        
        psifunc.c1=c;
        psifunc.kc1=kc;
        psifunc.class='OPT';
        
        
    elseif strcmp(rhofunc,'hyperbolic')
        
        if isempty(options.rhofuncparam)
            kdef=4.5;
        else
            kdef=options.rhofuncparam;
            kdef=kdef(1); % Instruction necessary for Ccoder
        end
        rhofuncparam=kdef;
        
        % Use (if possible) precalculated values of c,A,b,d and kc
        BDP=0.5:-0.01:0.01;
        KDEF=[4 4.5 5];
        [diffbdp,inddiffbdp]=min(abs(bdp(jj)-BDP));
        [diffk,inddiffk]=min(abs(kdef-KDEF));
        if diffbdp<1e-6 && diffk<1e-06
            % Load precalculated values of tuning constants
            Mat=coder.load('Hyp_BdpEff.mat','MatBDP');
            row=Mat.MatBDP(inddiffbdp,2:end,inddiffk);
            c=row(1); A=row(2); B=row(3); d=row(4); kc=row(5);
            
            %     % Use (if possible) precalculated values of c,A,b,d and kc
            %     if kdef == 4 && bdp(jj)==0.5
            %         c =2.158325031399727;
            %         A =1.627074124322223e-04;
            %         B =0.006991738279441;
            %         d =0.016982948780061;
            %         kc=0.010460153813287;
            %
            %     elseif kdef == 4.5 && bdp(jj)==0.5
            %         c =2.010311082005501;
            %         A =0.008931591866092;
            %         B =0.051928487236632;
            %         d =0.132017481327058;
            %         kc=0.074478627985759;
            %     elseif kdef == 5 && bdp(jj)==0.5
            %         c =1.900709968805313;
            %         A =0.023186529890225;
            %         B =0.083526860351552;
            %         d =0.221246910095216;
            %         kc=0.116380290077336;
            %     elseif kdef == 4.5 && bdp(jj)==0.25
            %         c =2.679452645778656;
            %         A =0.464174145115400;
            %         B =0.588821276233494;
            %         d =1.092639541625978;
            %         kc=0.410853066399912;
            
        else
            if coder.target('MATLAB')
                % Compute tuning constant associated to the requested breakdown
                % point
                [c,A,B,d]=HYPbdp(bdp(jj),1,kdef);
                % kc1 = E(rho) = sup(rho)*bdp
                kc=HYPrho(c,[c,kdef,A,B,d])*bdp(jj);
            else
                error('FSDA:Sregeda:WrongBdpHyp','Values of bdp for hyperbolic tangent estimator not supported for code generation')
            end
        end
        psifunc.c1=[c;kdef;A;B;d];
        psifunc.kc1=kc;
        
        psifunc.class='HYP';
        
        
    elseif strcmp(rhofunc,'hampel')
        
        if isempty(options.rhofuncparam)
            abc=[2;4;8];
        else
            abc=options.rhofuncparam;
        end
        rhofuncparam=abc(:);
        
        % Compute tuning constant associated to the requested breakdown
        % point
        c=HAbdp(bdp(jj),1,abc);
        % kc = E(rho) = sup(rho)*bdp
        kc=HArho(c*abc(3),[c; abc])*bdp(jj);
        
        
        
        psifunc.c1=[c;abc];
        psifunc.kc1=kc;
        
        psifunc.class='HA';
        
    elseif strcmp(rhofunc,'mdpd')
        % minimum density power divergence estimator
        % minimum density power divergence estimator
        
        c=PDbdp(bdp(jj));
        % kc1 = E(rho) = sup(rho)*bdp
        kc=bdp(jj);
        
        psifunc.c1=c;
        psifunc.kc1=kc;
        psifunc.class='PD';
        
    else
        error('FSDA:Sregeda:WrongRho','Specified rho function is not supported: possible values are ''bisquare'' , ''optimal'',  ''hyperbolic'', ''hampel'' ,''mdpd''')
    end
    
    if coder.target('MATLAB')
        XXrho=strcat(psifunc.class,'rho');
        hrho=str2func(XXrho);
    end
    
    bestbetas = zeros(bestr,p);
    bestsubset = bestbetas;
    bestscales = Inf * ones(bestr,1);
    sworst = Inf;
    
    % singsub = scalar which will contain the number of singular subsets which
    % are extracted (that is the subsets of size p which are not full rank)
    singsub=0;
    
    % ij is a scalar used to ensure that the best first bestr non singular
    % subsets are stored
    ij=1;
    
    % initialise and start timer.
    tsampling = ceil(min(nselected/100 , 1000));
    time=zeros(tsampling,1);
    
    
    for i = 1:nselected
        
        if i <= tsampling, tic; end
        
        % extract a subset of size p
        index = C(i,:);
        
        Xb = X(index,:);
        yb = y(index);
        
        % beta estimate
        beta = Xb\yb;
        
        if ~isnan(beta(1)) && ~isinf(beta(1))
            
            if refsteps > 0
                % do the refsteps refining steps
                % kc determines the breakdown point
                % c is linked to the biweight function
                tmp = IRWLSregS(y,X,beta,psifunc,refsteps,reftol);
                
                betarw = tmp.betarw;
                scalerw = tmp.scalerw;
                resrw = y - X * betarw;
            else
                % no refining steps
                betarw = beta;
                resrw = y - X * betarw;
                scalerw = median(abs(resrw))/.6745;
            end
            
            % to find s, save first the best bestr scales (deriving from non
            % singular subsets) and, from iteration bestr+1 (associated to
            % another non singular subset), replace the worst scale
            % with a better one as follows
            if ij > bestr
                % compute the objective function using current residuals and
                % the worst estimate of the scale among the bests previously
                % stored
                % scaletest = (1/n) \sum_i=1^n (u_i/(sworst*c))
                
                % Use function handle hrho. For example if
                % for optimal psi hrho=OPTrho
                
                
                if coder.target('MATLAB')
                    scaletest=mean(feval(hrho,resrw/sworst,psifunc.c1));
                    % OLD DELETED TO CHECK
                    % scaletest=mean(feval(hrho,resrw/sworst,c));
                    %scaletest = mean(TBrho(resrw/sworst,c));
                else
                    if strcmp(psifunc.class,'TB')
                        scaletest = mean(TBrho(resrw/sworst,psifunc.c1));
                        
                    elseif strcmp(psifunc.class,'OPT')
                        scaletest = mean(OPTrho(resrw/sworst,psifunc.c1));
                        
                    elseif strcmp(psifunc.class,'HA')
                        scaletest = mean(HArho(resrw/sworst,psifunc.c1));
                        
                    elseif strcmp(psifunc.class,'HYP')
                        scaletest = mean(HYPrho(resrw/sworst,psifunc.c1));
                        
                    elseif strcmp(psifunc.class,'PD')
                        scaletest = mean(PDrho(resrw/sworst,psifunc.c1));
                        
                    else
                        error('FSDA:Sreg:WrongRhoFunc','Wrong rho function supplied')
                    end
                    
                end
                
                
                if scaletest < kc
                    % Find position of the maximum value of previously stored
                    % best scales
                    [~,ind] = max(bestscales);
                    
                    
                    sbest = Mscale(resrw,psifunc,scalerw,minsctol);
                    %sbest1 = Mscale1(resrw,psifunc,scalerw,minsctol);
                    %sbest2 = minscale(resrw,psifunc.c1,psifunc.kc1,scalerw,minsctol);
                    %[sbest sbest1 sbest2]
                    
                    % Store sbest, betarw and indexes of the units forming the
                    % best subset associated with minimum value
                    % of the scale
                    bestscales(ind) = sbest;
                    bestbetas(ind,:) = betarw';
                    % best subset
                    bestsubset(ind,:)=index;
                    % sworst = the best scale among the bestr found up to now
                    sworst = max(bestscales);
                end
            else
                %bestscales(ij) = minscale(resrw,psifunc,scalerw,minsctol);
                
                bestscales(ij) = Mscale(resrw,psifunc,scalerw,minsctol);
                
                
                bestbetas(ij,:) = betarw';
                bestsubset(ij,:) = index;
                ij=ij+1;
            end
        else
            singsub=singsub+1;
        end
        
        % Write total estimate time to compute final estimate
        if i <= tsampling
            % sampling time until step tsampling
            time(i)=toc;
        elseif i==tsampling+1
            % stop sampling and print the estimated time
            if msg==1
                fprintf('Total estimated time to complete S estimate: %5.2f seconds \n', nselected*median(time));
            end
        end
        
    end
    
    % perform C-steps on best 'bestr' solutions, till convergence or for a
    % maximum of refstepsbestr steps using a convergence tolerance as specified
    % by scalar reftolbestr
    
    superbestscale = Inf;
    if ~coder.target('MATLAB')
        % Initializations necessary for MATLAB Ccoder
        superbestbeta=bestbetas(1,:)';
        superbestsubset=bestsubset(1,:);
        weights=y;
    end
    
    
    for i=1:bestr
        tmp = IRWLSregS(y,X,bestbetas(i,:)',psifunc,refstepsbestr,reftolbestr,bestscales(i));
        
        if tmp.scalerw < superbestscale
            superbestscale = tmp.scalerw;
            superbestbeta = tmp.betarw;
            superbestsubset = bestsubset(i,:);
            weights = tmp.weights;
        end
    end
    residuals=(y-X*superbestbeta)/superbestscale;
    outliers=seq( abs(residuals)>norminv(conflev) );
    
    if coder.target('MATLAB')
        if singsub/nselected>0.1
            disp('------------------------------')
            disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nselected) '%'])
        end
    end
    
    Residuals(:,jj)=residuals;
    Beta(:,jj)=superbestbeta;
    Scale(jj)=superbestscale;
    BS(:,jj)=superbestsubset;
    Weights(:,jj)=weights;
    Singsub(jj)=singsub;
    Outliers(outliers,jj)=true;
end

if coder.target('MATLAB')
    % Restore the previous state of the warnings
    warning(warnrank.state,'MATLAB:rankDeficientMatrix');
    warning(warnsing.state,'MATLAB:singularMatrix');
    warning(warnnear.state,'MATLAB:nearlySingularMatrix');
end

% Store in output structure monitored \beta, s, best subset, vector of
% S-weights and residuals
out.Beta = Beta;
out.Scale = Scale;
out.BS = BS;
out.Weights = Weights;
out.RES=Residuals;

% Store in output structure the number of singular subsets
out.Singsub=Singsub;
% Store in output structure the outliers
out.Outliers = Outliers;
% Store values of bdp which have been used
out.bdp=bdp;

out.class='Sregeda';

out.rhofunc=rhofunc;
% In case of Hampel or hyperbolic tangent estimator store the additional
% parameters which have been used
% For Hampel store a vector of length 3 containing parameters a, b and c
% For hyperbolic store the value of k= sup CVC
out.rhofuncparam=rhofuncparam;

if options.intercept==true
    % Store X (without the column of ones if there is an intercept)
    out.X=X(:,2:end);
else
    out.X=X;
end
% Store response
out.y=y;

if coder.target('MATLAB')
    % Plot residuals as function of the break down point
    if options.plots==1
        laby='Scaled S residuals';
        resfwdplot(out);
        ylabel(laby)
    end
end

end

% -------------------------------------------------------------------
% subfunction IRWLSreg
% -------------------------------------------------------------------

function outIRWLS = IRWLSregS(y,X,initialbeta,psifunc,refsteps,reftol,initialscale)
%IRWLSregS (iterative reweighted least squares) does refsteps refining steps from initialbeta for S estimator
%
%  Required input arguments:
%
%    y:         A vector with n elements that contains the response variable.
%               It can be both a row or column vector.
%    X :        Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p). Rows of X represent observations, and
%               columns represent variables.
% initialbeta : p x 1 vector containing initial estimate of beta
%     psifunc : a structure specifying the class of rho function to use, the
%               consistency factor, and the value associated with the
%               Expectation of rho in correspondence of the consistency
%               factor
%               psifunc must contain the following fields
%               c1 = consistency factor associated to required
%                    breakdown point
%               kc1= Expectation for rho associated with c1
%               class = string identyfing the rho (psi) function to use.
%                    Admissible values for class are 'bisquare', 'optimal'
%                    'hyperbolic' and 'hampel'
%               Remark: if class is 'hyperbolic' it is also necessary to
%                   specify parameters k (sup CVC), A, B and d
%               Remark: if class is 'hampel' it is also necessary to
%                   specify parameters a, b and c
%   refsteps  : scalar, number of refining (IRLS) steps
%   reftol    : relative convergence tolerance
%               Default value is 1e-7
%
%  Optional input arguments:
%
% initialscale: scalar, initial estimate of the scale. If not defined,
%               scaled MAD of residuals is used.
%
%  Output:
%
%  The output consists of a structure 'outIRWLS' containing the following fields:
%      betarw  : p x 1 vector. Estimate of beta after refsteps refining steps
%     scalerw  : scalar. Estimate of scale after refsteps refining step
%     weights  : n x 1 vector. Weights assigned to each observation
%
% In the IRWLS procedure the value of beta and the value of the scale are
% updated in each step

%% Beginning of code
c=psifunc.c1;
kc=psifunc.kc1;

% Residuals for the initialbeta
res = y - X * initialbeta;

% The scaled MAD of residuals is the initial scale estimate default value
if (nargin < 8)
    initialscale = median(abs(res))/.6745;
end

beta = initialbeta;
scale = initialscale;
newbeta =initialbeta;

if coder.target('MATLAB')
    
    XXrho=strcat(psifunc.class,'rho');
    hrho=str2func(XXrho);
    
    XXwei=strcat(psifunc.class,'wei');
    hwei=str2func(XXwei);
else
    weights=y; % initialization of weights necessary for MATLAB coder
end

iter = 0;
betadiff = 9999;

while ( (betadiff > reftol) && (iter < refsteps) )
    iter = iter + 1;
    
    if coder.target('MATLAB')
        
        % Solve for the scale
        meanrho=mean(feval(hrho,res/scale,c));
        scale = scale * sqrt(meanrho / kc );
        
        % Compute n x 1 vector of weights (using TB)
        weights = feval(hwei,res/scale,c);
        % weights = TBwei(res/scale,c);
    else
        if strcmp(psifunc.class,'TB')
            meanrho=mean(TBrho(res./scale,c));
            scale = scale * sqrt(meanrho / kc );
            weights = TBwei(res/scale,c);
            
        elseif strcmp(psifunc.class,'OPT')
            meanrho=mean(OPTrho(res/scale,c));
            scale = scale * sqrt(meanrho / kc );
            weights = OPTwei(res/scale,c);
            
        elseif strcmp(psifunc.class,'HA')
            meanrho=mean(HArho(res/scale,c));
            scale = scale * sqrt(meanrho / kc );
            weights = HAwei(res/scale,c);
            
        elseif strcmp(psifunc.class,'HYP')
            meanrho=mean(HYPrho(res/scale,c));
            scale = scale * sqrt(meanrho / kc );
            weights = HYPwei(res/scale,c);
            
        elseif strcmp(psifunc.class,'PD')
            meanrho=mean(PDrho(res/scale,c));
            scale = scale * sqrt(meanrho / kc );
            weights = PDwei(res/scale,c);
            
        else
            error('FSDA:Sreg:WrongRhoFunc','Wrong rho function supplied')
        end
        
    end
    
    sqweights = weights.^(1/2);
    
    % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
    Xw = bsxfun(@times, X, sqweights);
    yw = y .* sqweights;
    
    % estimate of beta from (re)weighted regression (RWLS)
    newbeta = Xw\yw;
    
    % exit from the loop if the new beta has singular values. In such a
    % case, any intermediate estimate is not reliable and we can just
    % keep the initialbeta and initial scale.
    if (any(isnan(newbeta)))
        newbeta = initialbeta;
        scale = initialscale;
        weights = NaN;
        break
    end
    
    % betadiff is linked to the tolerance (specified in scalar reftol)
    betadiff = norm(beta - newbeta,1) / norm(beta,1);
    
    % update residuals and beta
    res = y - X * newbeta;
    beta = newbeta;
    
end

% store final estimate of beta
outIRWLS.betarw = newbeta;
% store final estimate of scale
outIRWLS.scalerw = scale;
% store final estimate of the weights for each observation
outIRWLS.weights=weights;

end
%FScategory:REG-Regression

