function [out , varargout] = Taureg(y, X, varargin)
%Taureg computes Tau estimators in linear regression
%
%
%<a href="matlab: docsearchFS('Taureg')">Link to the help function</a>
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
%  intercept :  Indicator for constant term. Scalar. If 1, a model with
%               constant term will be fitted (default), else no constant
%               term will be included.
%               Example - 'intercept',1 
%               Data Types - double
%         bdp :  breakdown point. Scalar. 
%               It measures the fraction of outliers
%               the algorithm should resist. In this case any value greater
%               than 0 but smaller or equal than 0.5 will do fine (default=0.5).
%               Note that given bdp nominal
%               efficiency is automatically determined.
%                 Example - 'bdp',0.4
%                 Data Types - double      
%      eff     : nominal efficiency. Scalar.
%                Scalar defining nominal efficiency (i.e. a number between
%                 0.5 and 0.99). The default value is 0.95
%                 Asymptotic nominal efficiency is:
%                 Example - 'eff',0.99
%                 Data Types - double
%     rhofunc : rho function. String. String which specifies the rho function which must be used to
%               weight the residuals. Possible values are 'bisquare',
%               'optimal', 'hyperbolic', 'hampel'. 
%               'bisquare' uses Tukey's $\rho$ and $\psi$ functions, see
%               TBrho and TBpsi.
%               'optimal' uses optimal $\rho$ and $\psi$ functions, see
%               OPTrho and OPTpsi.
%               'hyperbolic' uses hyperbolic $\rho$ and $\psi$ functions,
%               see HYPrho and HYPpsi.
%               'hampel' uses Hampel $\rho$ and $\psi$ functions, see HArho
%               and HApsi. 
%               The default is bisquare
%                 Example - 'rhofunc','optimal' 
%                 Data Types - char
% rhofuncparam: Additional parameters for the specified rho function.
%               Scalar or vector.
%               For hyperbolic rho function it is possible to set up the
%               value of k = sup CVC (the default value of k is 4.5).
%               For Hampel rho function it is possible to define parameters
%               a, b and c (the default values are a=2, b=4, c=8)
%                 Example - 'rhofuncparam',5 
%                 Data Types - single | double
%       nsamp   : Number of subsamples which will be extracted to find the
%                 robust estimator. Scalar. If nsamp=0 all subsets will be extracted.
%                 They will be (n choose p).
%                 If the number of all possible subset is <1000 the
%                 default is to extract all subsets otherwise just 1000.
%                 Example - 'nsamp',1000 
%                 Data Types - single | double
%    refsteps : Number of refining iterations. Scalar. Number of refining iterationsin each
%               subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'nsamp',1000 
%                 Data Types - single | double
%     reftol  : scalar. Default value of tolerance for the refining steps.
%               The default value is 1e-6;
%                 Example - 'nsamp',1000 
%                 Data Types - single | double
%refstepsbestr: number of refining iterations for each best subset. Scalar.
%               Scalar defining number of refining iterations for each
%               best subset (default = 50).
%                 Example - 'refstepsbestr',10 
%                 Data Types - single | double
% reftolbestr : Tolerance for the refining steps. Scalar. 
%               Tolerance for the refining steps
%               for each of the best subsets
%               The default value is 1e-8;
%                 Example - 'reftolbestr',1e-10 
%                 Data Types - single | double
%     minsctol: tolerance for the iterative
%               procedure for finding the minimum value of the scale. Scalar. 
%               Value of tolerance for the iterative
%               procedure for finding the minimum value of the scale
%               for each subset and each of the best subsets
%               (It is used by subroutine minscale.m)
%               The default value is 1e-7;
%                 Example - 'minsctol',1e-7 
%                 Data Types - single | double
%      bestr  : number of "best betas" to remember. Scalar. Scalar defining number of "best betas" to remember from the
%               subsamples. These will be later iterated until convergence
%               (default=5)
%                 Example - 'bestr',10 
%                 Data Types - single | double
%     conflev :  Confidence level which is
%               used to declare units as outliers. Scalar. 
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%                 Example - 'conflev',0.99
%                 Data Types - double
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
%       nocheck : Check input arguments. Scalar. If nocheck is equal to 1 no check is performed on
%                 matrix y and matrix X. Notice that y and X are left
%                 unchanged. In other words the additional column of ones
%                 for the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1 
%               Data Types - double
%       plots : Plot on the screen. Scalar or structure.
%               If plots = 1, generates a plot with the robust residuals
%               against index number. The confidence level used to draw the
%               confidence bands for the residuals is given by the input
%               option conflev. If conflev is not specified a nominal 0.975
%               confidence interval will be used.
%                 Example - 'plots',0 
%                 Data Types - single | double
%       yxsave : Save matrices X and y. Scalar. If yxsave is equal to 1 the
%               response vector y and data matrix X are saved into the output
%                structure out. 
%               Default is 0, i.e. no saving is done.
%               Example - 'yxsave',1 
%               Data Types - double
%
%  Output:
%
%  out :     A structure containing the following fields
%
%            out.beta = vector containing the Tau estimator of regression
%                       coefficients
%            out.scale= scalar containing the estimate of the tau scale
%                       (sigma). This is the value of the objective function
%                        tau_scale = s_scale * average (\rho_2 (scaled residuals))
%              out.bs = p x 1 vector containing the units forming best subset
%                       associated with tau estimate of regression coefficient.
%        out.residuals= n x 1 vector containing the estimates of the robust
%                       scaled residuals
%        out.outliers = this output is present only if conflev has been
%                       specified. It is a vector containing the list of
%                       the units declared as outliers using confidence
%                       level specified in input scalar conflev
%         out.conflev = confidence level which is used to declare outliers.
%                       Remark: scalar out.conflev will be used to draw the
%                       horizontal line (confidence band) in the plot.
%         out.singsub = Number of subsets wihtout full rank. Notice that
%                       out.singsub > 0.1*(number of subsamples) produces a
%                       warning
%         out.weights = n x 1 vector containing the estimates of the weights
%         out.rhofunc = string identifying the rho function which has been
%                       used
%    out.rhofuncparam = vector which contains the additional parameters
%                       for the specified rho function which have been
%                       used. For hyperbolic rho function the value of
%                       k =sup CVC. For Hampel rho function the parameters
%                       a, b and c
%            out.y    = Response vector y. The field is present if option
%                       yxsave is set to 1.
%            out.X    = Data matrix X. The field is present if option
%                       yxsave is set to 1.
%           out.class = 'Taureg'
%
%  Optional Output:
%
%            C        : matrix containing the indices of the subsamples 
%                       extracted for computing the estimate (the so called
%                       elemental sets).
%
% See also: Sreg, MMreg
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
% Salibian-Barrera, M., Willems, G. and Zamar, R.H. (2008). The fast-tau
%   estimator for regression. Journal of Computational and Graphical
%   Statistics, 17, 659-682. (Referred below as SBWZ08)
% Yohai V.J. and Zamar R.H. (1988) High Breakdown-Point Estimates of
%   Regression by Means of the Minimization of an Efficient Scale,
%   Vol. 83, No. 402, pp. 406-413 (Referred below as YZ88)
%
% Acknowledgements:
%
% The kernel of the function is based on a MATLAB code downloaded from the
% web page of Dr. Matias Saliban-Barrera. However, all routines and
% subroutines have been completely redesigned, with considerable increase
% of the computational performance. Moreover in this function there is the
% possibility of choosing the rho (psi) function.
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Taureg')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % Taureg with all default options.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Taureg(ycont,X);
%}

%{
    %% Taureg with optional arguments.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Taureg(ycont,X,'plots',1);
%}

%{
    % Taureg with optimal rho function.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Taureg(ycont,X,'plots',1,'rhofunc','optimal');
%}

%{
    % Taureg with hyperbolic rho function.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    [out]=Taureg(ycont,X,'plots',1,'rhofunc','hyperbolic');
%}

%{
    % Taureg with Hampel rho function.
    % With parameters a=1.5 b=3.5 c=8. 
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    rhofuncparam=[1.5 3.5 8];
    [out]=Taureg(ycont,X,'plots',1,'rhofunc','hampel','rhofuncparam',rhofuncparam);
%}

%% Beginning of code
nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% default value of break down point
bdpdef=0.5;

% default value of asymptotic nominal efficiency
effdef=0.95;

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
% rhofuncdef='optimal';
% rhofuncdef='hampel';


% store default values in the structure options
options=struct('intercept',1,'nsamp',nsampdef,'refsteps',refstepsdef,...
    'reftol',reftoldef,'refstepsbestr',refstepsbestrdef,'reftolbestr',reftolbestrdef,...
    'minsctol',minsctoldef,'bestr',bestrdef,...
    'rhofunc',rhofuncdef,'bdp',bdpdef,'eff',effdef,'rhofuncparam','',...
    'plots',0,'conflev',0.975,'nocheck',0,'msg',1,'yxsave',0);

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:Taureg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options', the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end



% Get user values of warnings
warnrank=warning('query','MATLAB:rankDeficientMatrix');
warnsing=warning('query','MATLAB:singularMatrix');
warnnear=warning('query','MATLAB:nearlySingularMatrix');
% Set them to off inside this function
% At the end of the file they will be restored to previous values
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

bdp = options.bdp;              % break down point
eff = options.eff;              % break down point

refsteps = options.refsteps;    % refining steps
reftol = options.reftol;        % tolerance for refining steps
bestr = options.bestr;          % best betas for refining steps till convergence
nsamp = options.nsamp;          % subsamples to extract
minsctol = options.minsctol;    % tolerance for finding minimum value of the scale for each subset
refstepsbestr=options.refstepsbestr;  % refining steps for the best subsets
reftolbestr=options.reftolbestr;      % tolerance for refining steps for the best subsets
msg=options.msg;                % Scalar which controls the messages displayed on the screen
rhofunc=options.rhofunc;        % String which specifies the function to use to weight the residuals


% Find tuning constant c linked to rho function

% Optimal rho function is strictly increasing on [0 3c] and constant (equal to 3.25c^2) on [3c \infty)
% E(\rho) = kc = (3.25c^2)*bdp = TBrho(3*c,c)*bdp, being kc the K of
% Rousseeuw and Leroy (1987)


% Note that if \rho is standardized in such a way that (\rho(\infty))=1
% E(\rho) = kc = bdp

psifunc=struct;

if strcmp(rhofunc,'bisquare')
    
    % Tukey's biweight is strictly increasing on [0 c] and constant (equal to c^2/6) on [c \infty)
    % E(\rho) = kc = (c^2/6)*bdp = TBrho(c,c)*bdp, being kc the K of Rousseeuw
    % and Leroy (1987)
    
    
    % Compute tuning constant associated to the requested breakdown
    % point
    % For bdp =0.5 and Tukey biweight rho function c1=1.5476
    c1=TBbdp(bdp,1);
    % kc1 = E(rho) = sup(rho)*bdp
    kc1=TBrho(c1,c1)*bdp;
    
    % Compute tuning constant associated to the requested nominal efficiency
    % c2 = consistency factor for a given value of efficiency
    c2=TBeff(eff,1);
    % bdp2 = E(rho_2).
    % Note that if \rho is standardized in such a way that \rho(c)=1
    % E(rho_2) is nothing but the breakdown point associated to c2
    % else it is bdp2*sup(rho)
    kc2=TBc(c2,1)*TBrho(c2,c2);
    
    
    psifunc.c1=c1;
    psifunc.kc1=kc1;
    psifunc.c2=c2;
    psifunc.class='TB';
    
elseif strcmp(rhofunc,'optimal')
    
    % Compute tuning constant associated to the requested breakdown
    % point
    % For bdp =0.5 and optimal rho function c1=0.4046
    % Remark: given that in function OPTbdp rho function is defined in the interval 0---2c/3, 2c/3---3c/3, >3c/3
    % it is necessary to divide the output of OPTbdp by 3
    c1=OPTbdp(bdp,1)/3;
    % kc1 = E(rho) = sup(rho)*bdp
    kc1=OPTrho(3*c1,c1)*bdp;
    
    % Compute tuning constant associated to the requested nominal efficiency
    % c2 = consistency factor for a given value of efficiency
    % Remark: given that in function OPTeff rho function is defined in the interval 0---2c/3, 2c/3---3c/3, >3c/3
    % it is necessary to divide the output of OPTeff by 3
    c2=OPTeff(eff,1)/3;
    % b2 = E(rho_2).
    % Note that given that \rho is standardized in such a way that \rho(c)=1
    % E(rho_2) is nothing but the breakdown point associated to c2
    kc2=OPTc(3*c2,1);
    
    
    % Original precalculated values of c2 and b2 were
    %     c2 = 1.0900;
    %     bdp2 = 0.1278;
    
    psifunc.c1=c1;
    psifunc.kc1=kc1;
    psifunc.c2=c2;
    psifunc.class='OPT';
    
elseif strcmp(rhofunc,'hyperbolic')
    
    if isempty(options.rhofuncparam)
        kdef=4.5;
    else
        kdef=options.rhofuncparam;
    end
    rhofuncparam=kdef;
    
    % Use (if possible) precalculated values of c,A,b,d and kc
    if kdef == 4 && bdp==0.5
        c1 =2.158325031399727;
        A1 =1.627074124322223e-04;
        B1 =0.006991738279441;
        d1 =0.016982948780061;
        kc1=0.010460153813287;
        
    elseif kdef == 4.5 && bdp==0.5
        c1 =2.010311082005501;
        A1 =0.008931591866092;
        B1 =0.051928487236632;
        d1 =0.132017481327058;
        kc1=0.074478627985759;
    elseif kdef == 5 && bdp==0.5
        c1 =1.900709968805313;
        A1 =0.023186529890225;
        B1 =0.083526860351552;
        d1 =0.221246910095216;
        kc1=0.116380290077336;
        
    else
        
        % Compute tuning constant associated to the requested breakdown
        % point
        [c1,A1,B1,d1]=HYPbdp(bdp,1,kdef);
        % kc1 = E(rho) = sup(rho)*bdp
        kc1=HYPrho(c1,[c1,kdef,A1,B1,d1])*bdp;
    end
    
    if kdef == 4 && eff==0.85
        c2 =3.212800979614258;
        A2 =0.570183575755717;
        B2 =0.696172437281084;
        d2 =1.205900263786317;
        kc2=0.439232837342420;
    elseif kdef == 4.5 && eff==0.85
        c2 =3.032387733459473;
        A2 =0.615717108822885;
        B2 = 0.723435958485131;
        d2 =1.321987605094910;
        kc2=0.448833150947462;
    elseif kdef == 5 && eff==0.85
        c2 =2.911890029907227;
        A2 =0.650228046997054;
        B2 =0.743433840145084;
        d2 =1.419320821762087;
        kc2=0.455326016919854;
        
    elseif kdef == 4 && eff==0.95
        c2 =4.331634521484375;
        A2 =0.754327484845243;
        B2 =0.846528826589308;
        d2 =1.480099129676819;
        kc2=0.473872135913024;
    elseif kdef == 4.5 && eff==0.95
        c2 =3.866390228271484;
        A2 =0.791281464739131;
        B2 =0.867016329355630;
        d2 =1.610621500015260;
        kc2=0.479388475649576;
    elseif kdef == 5 && eff==0.95
        c2 =3.629499435424805;
        A2 =0.818876452066880;
        B2 =0.882004888111327;
        d2 =1.723768949508668;
        kc2=0.483053062139011;
        
    else
        
        % Compute tuning constant associated to the requested nominal efficiency
        % c2 = consistency factor for a given value of efficiency
        [c2,A2,B2,d2]=HYPeff(eff,1,kdef);
        % b2 = E(rho_2).
        % Note that given that if \rho is standardized in such a way that \rho(c)=1
        % E(rho_2) is nothing but the breakdown point associated to c2
        % else it is bdp2*sup(rho)
        kc2=HYPc(c2,1,'k',kdef,'param',[A2 B2 d2])*HYPrho(c2,[c2 kdef A2 B2 d2]);
    end
    
    psifunc.c1=[c1,kdef,A1,B1,d1];
    psifunc.kc1=kc1;
    
    psifunc.c2=[c2,kdef,A2,B2,d2];
    psifunc.class='HYP';
    
    c1=psifunc.c1;
    c2=psifunc.c2;
    
elseif strcmp(rhofunc,'hampel')
    
    if isempty(options.rhofuncparam)
        abc=[2,4,8];
    else
        abc=options.rhofuncparam;
    end
    rhofuncparam=abc;
    
    % Compute tuning constant associated to the requested breakdown
    % point
    c1=HAbdp(bdp,1,abc);
    % kc = E(rho) = sup(rho)*bdp
    kc1=HArho(c1*abc(3),[c1, abc])*bdp;
    
    
    % Compute tuning constant associated to the requested nominal efficiency
    % c2 = consistency factor for a given value of efficiency
    c2=HAeff(eff,1,abc);
    % b2 = E(rho_2).
    % Note that given that if \rho is standardized in such a way that \rho(c)=1
    % E(rho_2) is nothing but the breakdown point associated to c2
    % else it is bdp2*sup(rho)
    kc2=HAc(c2,1,'param',abc)* HArho(c2*abc(3),[c2 abc]);
    
    psifunc.c1=[c1,abc];
    psifunc.kc1=kc1;
    
    psifunc.c2=[c2,abc];
    psifunc.class='HA';
    
    c1=psifunc.c1;
    c2=psifunc.c2;
    
else
    error('FSDA:Taureg:WrongRho','Specified rho function is not supported: possible values are ''bisquare'' , ''optimal'',  ''hyperbolic'', ''hampel''')
    
end

XXrho=strcat(psifunc.class,'rho');
hrho=str2func(XXrho);


% bestbetas = matrix which contains the best betas (in the rows)
bestbetas = zeros(bestr,p);
bestsubset = bestbetas;
% worstind = scalar. Index associated to element of vector besttauscales containing
% the largest value of the tau scale among the bestr scales
worstind = 1;
% besttauscales = vector which contains the best tau scales
besttauscales = 1e20 * ones(bestr,1);
% bestscale = vector which contains the best S scales
bestscales = besttauscales;

worsts = Inf;
% the worst (largest estimate of the tau scale) among the bestr best ones
worsttau = Inf;
% worstres = vector of residuals associated with worsttau
worstres = y;

% singsub = scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size p which are not full rank)
singsub=0;

% ij is a scalar used to ensure that the best first bestr non singular
% subsets are stored
ij=1;

%% Extract in the rows of matrix C the indexes of all required subsets
[C,nselected] = subsets(nsamp,n,p,ncomb,msg);
% Store the indices in varargout
if nargout==2
    varargout={C};
end

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
        
        % do refsteps refining (concentration) steps
        if refsteps>0
            tmp = IRWLSregTau(y, X, beta, psifunc, refsteps, reftol);
            betarw = tmp.betarw;
            scalerw = tmp.scalerw;
            resrw = y - X * betarw;
        else
            % no refining steps
            betarw = beta;
            resrw = y - X * betarw;
            scalerw = median(abs(resrw))/.6745;
        end
        
        
        % To find the best estimate of the tau scale, save first the best bestr
        % scales (deriving from non singular subsets) and, from iteration
        % bestr+1 (associated to another non singular subset), replace the
        % worst scale with a better one as follows
        if ij > bestr
            % Condition 1 and 2 to understand whether this subset may
            % produce a smaller estimate ot the tau scale
            % Condition 1: equation (2.5) of SBWZ08
            % This condition is necessary and sufficient to obtain a
            % smaller value of the s scale
            % but it is just necessary (not sufficient) to obtain a
            % smaller estimate of the tau scale
            
            % Use function handle hrho. For example if
            % for optimal psi hrho=OPTrho
            % scaletest1 = mean(OPTrho(resrw/worsts,c1)) < kc1;
            scaletest1=mean(feval(hrho,resrw/worsts,c1))< kc1;
            
            % Condition 2: equation (2.6) of SBWZ08
            % Also this condition is necessary but not sufficient
            % scaletest2 = sum(OPTrho(resrw/worsts,c2)) < sum(OPTrho(worstres/worsts,c2));
            scaletest2=sum(feval(hrho,resrw/worsts,c2))< sum(feval(hrho,worstres/worsts,c2));
            
            if (scaletest1 || scaletest2)
                
                % Find M estimator of the scale using consistency factor
                % c1 (associated to break down point)
                news = Mscale(resrw, psifunc, scalerw, minsctol);
                % Find associated value of the tau scale
                % To be precise newtau is \tan(beta)*sqrt(bdp2). See equation
                % (1.3) of SBWZ08 therefore outside the loop it will be
                % necessary to divide newtau by sqrt bdp2.
                % Note that in this case rhoOptfun uses consistency factor c2
                % (associated to nominal asymptotic efficiency)
                
                % newtau = news * sqrt(mean(OPTrho(resrw/news,c2)));
                newtau = news * sqrt(mean(feval(hrho,resrw/news,c2)));
                
                % Given that the two previous conditions were just necessary but not
                % sufficient the following if is necessary
                if newtau < worsttau
                    bestscales(worstind) = news;
                    besttauscales(worstind) = newtau;
                    bestsubset(worstind,:)=index;
                    bestbetas(worstind,:) = betarw';
                    [worsttau,worstind] = max(besttauscales);
                    worsts = bestscales(worstind);
                    % vector worstres will be the input of condition 2
                    worstres = y - X * bestbetas(worstind,:)';
                end
            end
        else
            
            news = Mscale(resrw, psifunc, scalerw, minsctol);
            % news1 = minscale(resrw, c1, kc1, scalerw, minsctol);
            
            % newtau = news * sqrt(mean(OPTrho(resrw/news,c2)));
            newtau = news * sqrt(mean(feval(hrho,resrw/news,c2)));
            
            bestscales(ij) = news;
            besttauscales(ij) =newtau;
            bestsubset(ij,:) = index;
            bestbetas(ij,:) = betarw';
            if ij==bestr
                [worsttau,worstind] = max(besttauscales);
                worsts = bestscales(worstind);
                % vector worstres will be the input of condition 2
                worstres = y - X * bestbetas(worstind,:)';
            end
            ij=ij+1;
        end
    else
        singsub=singsub+1;
    end
    
    % Write total estimation time to compute final estimate
    if i <= tsampling
        % sampling time until step tsampling
        time(i)=toc;
    elseif i==tsampling+1
        % stop sampling and print the estimated time
        if msg==1
            fprintf('Total estimated time to complete tau estimate: %5.2f seconds \n', nselected*median(time));
        end
    end
    
end

% perform C-steps on best 'bestr' solutions, till convergence or for a
% maximum of refstepsbestr steps using a convergence tolerance as specified
% by scalar reftolbestr

superbesttauscale = Inf;

for i=1:bestr
    tmp = IRWLSregTau(y, X, bestbetas(i,:)', psifunc, refstepsbestr, reftolbestr, bestscales(i));
    resrw = y - X * tmp.betarw;
    
    % tauscalerw = tmp.scalerw * sqrt(mean(OPTrho(resrw/tmp.scalerw,c2)));
    tauscalerw = tmp.scalerw * sqrt(mean(feval(hrho,resrw/tmp.scalerw,c2)));
    
    
    if tauscalerw < superbesttauscale
        superbesttauscale = tauscalerw;
        superbestscale=tmp.scalerw;
        superbestbeta = tmp.betarw;
        superbestsubset = bestsubset(i,:);
        weights = tmp.weights;
    end
end

superbestscale = Mscale(y-X*superbestbeta, psifunc, superbestscale, minsctol);

% superbesttauscale = superbestscale * sqrt(mean(OPTrho((y-X*superbestbeta)/superbestscale,c2)));
superbesttauscale = superbestscale * sqrt(mean(feval(hrho,(y-X*superbestbeta)/superbestscale,c2)));


% Store in output structure \beta, s, best subset and vector of S-weights
out.beta = superbestbeta;

% Rescale the estimate of the scale with the breakdownpoint associated with
% consistency factor c2
superbesttauscale = superbesttauscale / sqrt(kc2);

out.scale = superbesttauscale;
out.bs = superbestsubset;
out.weights = weights;

% compute and store in output structure the Tau robust scaled residuals
out.residuals=(y-X*out.beta)/out.scale;

% Store in output structure the number of singular subsets
out.singsub=singsub;
if singsub/nselected>0.1
    disp('------------------------------')
    disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nselected) '%'])
end

% Restore the previous state of the warnings
warning(warnrank.state,'MATLAB:rankDeficientMatrix');
warning(warnsing.state,'MATLAB:singularMatrix');
warning(warnnear.state,'MATLAB:nearlySingularMatrix');


% Store in output structure the outliers found with confidence level conflev
conflev = options.conflev;
out.conflev = conflev;

conflev = (conflev+1)/2;
seq = 1:n;
out.outliers = seq( abs(out.residuals)>norminv(conflev) );

out.rhofunc=rhofunc;
% In case of Hampel or hyperbolic tangent estimator store the additional
% parameters which have been used
% For Hampel store a vector of length 3 containing parameters a, b and c
% For hyperbolic store the value of k= sup CVC
if exist('rhofuncparam','var')
    out.rhofuncparam=rhofuncparam;
end


if options.yxsave
    if options.intercept==1
        % Store X (without the column of ones if there is an intercept)
        out.X=X(:,2:end);
    else
        out.X=X;
    end
    % Store response
    out.y=y;
end

% Plot resindexplot with outliers highlighted
if options.plots==1
    laby='Scaled tau residuals';
    resindexplot(out.residuals,'conflev',out.conflev,'laby',laby,'numlab',out.outliers);
end

end


function outIRWLS = IRWLSregTau(y, X, initialbeta, psifunc, refsteps, reftol, initialscale)
%IRWLSregTau (iterative reweighted least squares) does refsteps refining steps from initialbeta in Tau estimator
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
%               c2 = consistency factor associated to required
%                    nominal efficiency
%               class = string identyfing the rho (psi) function to use.
%                    Admissible values for class are 'bisquare', 'optimal'
%                    'hyperbolic' and 'hampel'
%               Remark: if class is 'hyperbolic' it is also necessary to
%                   specify parameters k (sup CVC), A, B and d
%               Remark: if class is 'hampel' it is also necessary to
%                   specify parameters a, b and c
%
%
%  Optional input arguments:
%
%   refsteps  : scalar, number of refining (IRLS) steps
%   reftol    : relative convergence tolerance
%               Default value is 1e-7
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

% Remark:
% if psifunc = hyperbolic
%   c1 and c2 will be vectors of length 5
%   containing parameters ctuning,ktuning,A,B,d
% if psifunc = hampel
%   c1 and c2 will be vectors of length 4
%   containing parameters ctuning,a,b,c
% else (Tukey biweight and hyperbolic) c1 is a scalar
c1=psifunc.c1;


kc1=psifunc.kc1;
c2=psifunc.c2;

if (nargin < 5)
    refsteps=100;
end

if (nargin < 6)
    reftol=1e-7;
end

% Residuals for the initialbeta
res = y - X * initialbeta;

% The scaled MAD of residuals is the initial scale estimate default value
if (nargin < 7)
    initialscale = median(abs(res))/.6745;
end

beta = initialbeta;
scale = initialscale;

XXrho=strcat(psifunc.class,'rho');
hrho=str2func(XXrho);

XXpsix=strcat(psifunc.class,'psix');
hpsix=str2func(XXpsix);

XXwei=strcat(psifunc.class,'wei');
hwei=str2func(XXwei);


iter = 0;

betadiff=9999;
while (betadiff > reftol) && (iter < refsteps)
    
    iter = iter + 1;
    
    % One iteration for the scale
    % meanrho=mean( OPTrho(res/scale,c1) );
    meanrho=mean(feval(hrho,res/scale,c1));
    
    scale = scale*sqrt( meanrho / kc1 );
    
    scaledres = res/scale;
    
    % Compute n x 1 vector of weights (using requested weight function)
    
    % Wn_numer = \sum_{i=1}^n [ 2 \rho(x) -\psi(x)*x ]
    % \sum [ 2 \rho(scaledres) -\psi(scaledres)*scaledres ]
    % (using consistency factor c2)
    %
    % sum(2*OPTrho(scaledres,c2)-OPTpsix(scaledres,c2))
    % Wn_numer = sum(WnumerOptfun(scaledres,c2));
    % Wn_numer=sum(2*OPTrho(scaledres,c2)-OPTpsix(scaledres,c2));
    
    Wn_numer=sum(2*feval(hrho,scaledres,c2)-feval(hpsix,scaledres,c2));
    
    % Wn_denom = \sum_{i=1}^n psi(scaledres)*scaledres
    % See denominator of equation 2.8 of Yohai and Zamar
    % Wn_denom = sum(OPTpsix(scaledres,c1));
    Wn_denom = sum(feval(hpsix,scaledres,c1));
    
    Wn = Wn_numer/Wn_denom;
    % OPTwei = psi(x)/x
    % Vector weights is Wn*psi(x,c1)/x+psi(x,c2)/x
    % see page 409 YZ88 or page 664 below equation (2.1) of SBWZ08
    % weights = Wn*OPTwei(scaledres,c1)+OPTwei(scaledres,c2);
    weights = Wn*feval(hwei,scaledres,c1)+feval(hwei,scaledres,c2);
    
    
    sqweights = weights.^(1/2);
    
    % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
    Xw = bsxfun(@times, X, sqweights);
    yw = y .* sqweights;
    
    % New estimate of beta from (re)weighted regression (RWLS)
    newbeta = Xw\yw;
    
    
    if (any(isnan(newbeta)))
        newbeta = initialbeta;
        scale = initialscale;
        weights = NaN;
        break
    end
    
    %  betadiff = norm(beta - newbeta)/sqrt(p);
    
    % betadiff is linked to the tolerance (specified in scalar reftol)
    betadiff = norm(beta - newbeta,1) / norm(beta,1);
    
    
    res = y - X * newbeta;
    beta = newbeta;
end

outIRWLS.betarw = newbeta;
outIRWLS.scalerw = scale;
outIRWLS.iters = iter;
outIRWLS.weights=weights;
end

%FScategory:REG-Regression




