function [out , varargout] = MMreg(y,X,varargin)
%MMreg computes MM estimator of regression coefficients
%
%<a href="matlab: docsearchFS('MMreg')">Link to the help function</a>
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
%   intercept : Indicator for constant term. true (default) | false. 
%               Indicator for the constant term (intercept) in the fit,
%               specified as the comma-separated pair consisting of
%               'intercept' and either true to include or false to remove
%               the constant term from the model.
%               Example - 'intercept',false
%               Data Types - boolean
%
%  InitialEst : starting values of the MM-estimator. [] (default) or structure.
%               InitialEst must contain the following fields
%               InitialEst.beta =  v x 1 vector (estimate of the centroid)
%               InitialEst.scale = scalar (estimate of the scale parameter).
%               If InitialEst is empty (default)
%               program uses S estimators. In this last case it is
%               possible to specify the options given in function Sreg.
%               Example - 'InitialEst',[]
%               Data Types - struct
%  Soptions  :  options to pass to Sreg. Name value pairs.
%               Options if initial estimator is S and InitialEst is empty.
%               The options are: Srhofunc,Snsamp,Srefsteps, Sreftol, Srefstepsbestr,
%               Sreftolbestr, Sminsctol, Sbestr.
%               See function Sreg.m for more details on these options.
%               It is necessary to add to the S options the letter
%               S at the beginning. For example, if you want to use the
%               optimal rho function the supplied option is
%               'Srhofunc','optimal'. For example, if you want to use 3000
%               subsets, the supplied option is 'Snsamp',3000
%               Example - 'Snsamp',1000
%               Data Types - single | double
%
%
%               MM options
%
%      eff     : nominal efficiency. Scalar.
%                Scalar defining nominal efficiency (i.e. a number between
%                 0.5 and 0.99). The default value is 0.95
%                 Asymptotic nominal efficiency is:
%                 $(\int \psi' d\Phi)^2 / (\psi^2 d\Phi)$
%                 Example - 'eff',0.99
%                 Data Types - double
%     effshape : location or scale efficiency. dummy scalar.
%                If effshape=1 efficiency refers to shape
%                efficiency else (default) efficiency refers to location
%                 Example - 'effshape',1
%                 Data Types - double
%     rhofunc : rho function. String. String which specifies the rho function which must be used to
%               weight the residuals. Possible values are
%               'bisquare'
%               'optimal'
%               'hyperbolic'
%               'hampel'
%               'bisquare' uses Tukey's $\rho$ and $\psi$ functions.
%               See TBrho and TBpsi
%               'optimal' uses optimal $\rho$ and $\psi$ functions.
%               See OPTrho and OPTpsi
%               'hyperbolic' uses hyperbolic $\rho$ and $\psi$ functions.
%               See HYPrho and HYPpsi
%               'hampel' uses Hampel $\rho$ and $\psi$ functions.
%               See HArho and HApsi
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
%     refsteps  : Maximum iterations. Scalar.
%                 Scalar defining maximum number of iterations in the MM
%                 loop. Default value is 100.
%                 Example - 'refsteps',10
%                 Data Types - double
%       tol    : Tolerance. Scalar.
%                 Scalar controlling tolerance in the MM loop.
%                 Default value is 1e-7
%                 Example - 'tol',1e-10
%                 Data Types - double
%     conflev :  Confidence level which is
%               used to declare units as outliers. Scalar.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975
%                 Example - 'conflev',0.99
%                 Data Types - double
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
%       yxsave : the response vector y and data matrix X are saved into the output
%                structure out. Scalar.
%               Default is 0, i.e. no saving is done.
%               Example - 'yxsave',1
%               Data Types - double
%       msg    :  Level of output to display. Scalar. It controls whether
%                 to display or not messages on the screen
%                 If msg==1 (default) messages are displayed on the screen about
%                   step in which signal took place
%                 else no message is displayed on the screen.
%               Example - 'msg',1
%               Data Types - double
%
%  Output:
%
%
%  out :     A structure containing the following fields:
%       out.beta        =   p x 1 vector containing MM estimate of
%                           regression coefficients.
%       out.auxscale    =   scalar, S estimate of the scale (or supplied
%                           external estimate of scale, if option InitialEst
%                           is not empty).
%       out.residuals	=   n x 1 vector containing standardized MM
%                           residuals.
%      out.fittedvalues =   n x 1 vector containing the fitted values.
%                           out.residuals=(y-X*out.beta)/out.auxscale
%       out.weights     =   n x 1 vector. Weights assigned to each observation
%       out.Sbeta       :   p x 1 vector containing S estimate of regression
%                           coefficients (or supplied initial external
%                           estimate of regression coefficients, if option
%                           InitialEst is not empty)
%       out.Ssingsub    =   Number of subsets without full rank in the S
%                           preliminary part. Notice that
%                           out.singsub > 0.1*(number of subsamples)
%                           produces a warning
%       out.outliers    :   1 x k vectors containing the outliers which
%                           have been found
%       out.conflev     =   Confidence level that was used to declare outliers
%           out.rhofunc =   string identifying the rho function which has been
%                           used. If a different rho function is specified
%                           for S and MM loop then insted of out.rhofunc we
%                           will have out.rhofuncS and out.rhofuncMM.
%      out.rhofuncparam =   vector which contains the additional parameters
%                           for the specified rho function which have been
%                           used. For hyperbolic rho function the value of
%                           k =sup CVC. For Hampel rho function the parameters
%                           a, b and c. If a different rho function is
%                           specified for S and MM loop then insted of
%                           out.rhofuncparam we will have out.rhofuncparamS
%                           and out.rhofuncparamMM.
%            out.y      =   response vector Y. The field is present if option
%                           yxsave is set to 1.
%            out.X      =   data matrix X. The field is present if option
%                           yxsave is set to 1.
%       out.class       =   'MMreg'
%
%
%  Optional Output:
%
%            C        : matrix containing the indices of the subsamples
%                       extracted for computing the estimate (the so called
%                       elemental sets).
%
%
% See also: Sreg
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Acknowledgements:
%
% This function follows the lines of MATLAB/R code developed during the
% years by many authors.
% For more details see http://www.econ.kuleuven.be/public/NDBAE06/programs/
% and the R library robustbase http://robustbase.r-forge.r-project.org/
% The core of these routines, e.g. the resampling approach, however, has
% been completely redesigned, with considerable increase of the
% computational performance.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('MMreg')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % MMreg with all default options.
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
    [out]=MMreg(ycont,X);
%}

%{
    % MMreg with optional input arguments.
    % MMreg using the hyperbolic rho function.
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
    [out]=MMreg(ycont,X,'Srhofunc','optimal');
%}

%{
    % MMreg with optional input arguments.
    % MMreg using the OLS estimates ac InitialEst.
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
    % OLS estimates
    bols=[ones(n,1) X]\y;
    res=y-[ones(n,1) X]*bols;
    sols=sqrt((res'*res)/(n-p-1));
    InitialEst.beta=bols;
    InitialEst.scale=sols;
    [out]=MMreg(ycont,X,'InitialEst',InitialEst);
%}

%{
    %% Comparing the output of different MMreg runs.
    state=100;
    randn('state', state);
    n=100;
    X=randn(n,3);
    bet=[3;4;5];
    y=3*randn(n,1)+X*bet;
    y(1:20)=y(1:20)+13;

    %For outlier detection we consider both the nominal individual 1%
    %significance level and the simultaneous Bonferroni confidence level.

    % Define nominal confidence level
    conflev=[0.99,1-0.01/length(y)];
    % Define number of subsets
    nsamp=3000;
    % Define the main title of the plots
    titl='';

    % MM  estimators
    [outMM]=MMreg(y,X,'conflev',conflev(1));
    laby='Scaled MM residuals';
    resindexplot(outMM.residuals,'title',titl,'laby',laby,'numlab','','conflev',conflev)
    % In this example MM estimator seems to detect half of the outlier with a Bonferroni significance level.
    % By simply changing the seed to 543 (state=543), using a Bonferroni size
    %of 1%, no unit is declared as outlier and just half of them using the 99%
    %band.
%}

%{
    % Comparison between direct call to MMreg and call to Sreg and MMregcore.
    n=30;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    % Two different rho functions are used for S and MM
    rhofuncS='hyperbolic';
    rhofuncMM='hampel';
    % Direct call to MMreg
    [out]=MMreg(ycont,X,'Srhofunc',rhofuncS,'rhofunc',rhofuncMM,'Snsamp',0);

    % Call to Sreg and then to MMregcore
    [outS]=Sreg(ycont,X,'rhofunc',rhofuncS,'nsamp',0);
    outMM=MMregcore(ycont,X,outS.beta,outS.scale,'rhofunc',rhofuncMM);
    disp('Difference between direct call to S and the calls to Sreg and MMregcore')
    max(abs([out.beta-outMM.beta]))
%}

%% Beginning of code

% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% default values for the initial S estimate:

% default value of break down point
Sbdpdef=0.5;

% default values of subsamples to extract
ncomb=bc(n,p);
Snsampdef=min(1000,ncomb);

% default value of number of refining iterations (C steps) for each extracted subset
Srefstepsdef=3;
% default value of tolerance for the refining steps convergence for  each extracted subset
Sreftoldef=1e-6;
% default value of number of best locs to remember
Sbestrdef=5;
% default value of number of refining iterations (C steps) for best subsets
Srefstepsbestrdef=50;
% default value of tolerance for the refining steps convergence for best subsets
Sreftolbestrdef=1e-8;
% default value of tolerance for finding the minimum value of the scale
% both for each extracted subset and each of the best subsets
Sminsctoldef=1e-7;

% rho (psi) function which has to be used to weight the residuals
Srhofuncdef='bisquare';
rhofuncdef=Srhofuncdef;

options=struct('intercept',1,'InitialEst','','Snsamp',Snsampdef,'Srefsteps',Srefstepsdef,...
    'Sbestr',Sbestrdef,'Sreftol',Sreftoldef,'Sminsctol',Sminsctoldef,...
    'Srefstepsbestr',Srefstepsbestrdef,'Sreftolbestr',Sreftolbestrdef,...
    'Sbdp',Sbdpdef,'Srhofunc',Srhofuncdef,'Srhofuncparam','','nocheck',0,'eff',0.95,'effshape',0,...
    'rhofunc',rhofuncdef,'rhofuncparam','',...
    'refsteps',100,'tol',1e-7,'conflev',0.975,'plots',0,'yxsave',0,'msg',1);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:MMreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:MMreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if nargin > 2
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

msg=options.msg;

% intercept=options.intercept;

% InitialEst = structure which contains initial estimate of beta and sigma
% If InitialEst is empty then initial estimates of beta and sigma come from
% S-estimation
InitialEst=options.InitialEst;

if isempty(InitialEst)
    
    bdp = options.Sbdp;              % break down point
    refsteps = options.Srefsteps;    % refining steps
    bestr = options.Sbestr;          % best locs for refining steps till convergence
    nsamp = options.Snsamp;          % subsamples to extract
    reftol = options.Sreftol;        % tolerance for refining steps
    minsctol = options.Sminsctol;    % tolerance for finding minimum value of the scale for each subset
    refstepsbestr=options.Srefstepsbestr;  % refining steps for the best subsets
    reftolbestr=options.Sreftolbestr;      % tolerance for refining steps for the best subsets
    
    rhofuncS=options.Srhofunc;              % rho function which must be used for S estimator
    rhofuncparamS=options.Srhofuncparam;    % eventual additional parameters associated to the rho function for S estimator
    
    
    
    % first compute S-estimator with a fixed breakdown point
    
    % SR is the routine which computes S estimates of beta and sigma in regression
    % Note that intercept is taken care of by chkinputR call.
    if nargout==2
        [Sresult , C] = Sreg(y,X,'nsamp',nsamp,'bdp',bdp,'refsteps',refsteps,'bestr',bestr,...
            'reftol',reftol,'minsctol',minsctol,'refstepsbestr',refstepsbestr,...
            'reftolbestr',reftolbestr,'rhofunc',rhofuncS,'rhofuncparam',rhofuncparamS,...
            'nocheck',1,'msg',msg);
        
        varargout = {C};
    else
        Sresult = Sreg(y,X,'nsamp',nsamp,'bdp',bdp,'refsteps',refsteps,'bestr',bestr,...
            'reftol',reftol,'minsctol',minsctol,'refstepsbestr',refstepsbestr,...
            'reftolbestr',reftolbestr,'rhofunc',rhofuncS,'rhofuncparam',rhofuncparamS,...
            'nocheck',1,'msg',msg);
    end
    
    bs = Sresult.beta;
    ss = Sresult.scale;
    singsub=Sresult.singsub;
    
    
else
    bs = InitialEst.beta;
    ss = InitialEst.scale;
    singsub=0;
    % In this case there is no preliminary S estimator
    rhofuncS=[];
    % Sresult is an empty structure
    Sresult=struct;
end


rhofuncMM=options.rhofunc;              % rho function which must be used for MM loop
rhofuncparamMM=options.rhofuncparam;    % eventual additional parameters associated to the rho function for MM loop


% Asymptotic nominal efficiency (for location or shape)
eff = options.eff;

% effshape = scalar which specifies whether nominal efficiency refers to location or scale
effshape = options.effshape;

% refsteps = maximum number of iteration in the MM step
refsteps = options.refsteps;

% tol = tolerance to declare convergence in the MM step
tol = options.tol;

% MMregcore = function which does IRWLS steps from initialbeta (bs) and sigma (ss)
% Notice that the estimate of sigma (scale) remains fixed
plots=options.plots;
conflev=options.conflev;

outIRW = MMregcore(y,X,bs,ss,'eff',eff,'effshape',effshape,...
    'rhofunc',rhofuncMM,'rhofuncparam',rhofuncparamMM,...
    'refsteps',refsteps,'reftol',tol,'conflev',conflev,'plots',plots,'nocheck',1,'msg',msg);


out = struct;
out.beta = outIRW.beta;
out.auxscale = ss;
out.fittedvalues = X*outIRW.beta;
out.residuals = (y-out.fittedvalues)/ss; % MM scaled residuals
out.Sbeta = bs;
out.Ssingsub=singsub;
out.weights=outIRW.weights;
out.outliers=outIRW.outliers;
out.conflev=conflev;
out.class='MMreg';
if isequal(rhofuncS,rhofuncMM)
    out.rhofunc=rhofuncS;
    % In case of Hampel or hyperbolic tangent estimator store the additional
    % parameters which have been used
    % For Hampel store a vector of length 3 containing parameters a, b and c
    % For hyperbolic store the value of k= sup CVC
    if isfield(Sresult, 'rhofuncparam')
        out.rhofuncparam=Sresult.rhofuncparam;
    end
else
    out.rhofuncS=rhofuncS;
    out.rhofuncMM=rhofuncMM;
    % In case of Hampel or hyperbolic tangent estimator store the additional
    % parameters which have been used
    % For Hampel store a vector of length 3 containing parameters a, b and c
    % For hyperbolic store the value of k= sup CVC
    if isfield(Sresult, 'rhofuncparam')
        out.rhofuncparamS=Sresult.rhofuncparam;
    end
    
    if isfield(outIRW,'rhofuncparam')
        out.rhofuncparamMM=outIRW.rhofuncparam;
    end
    
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

end
%FScategory:REG-Regression