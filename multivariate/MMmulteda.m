function [out , varargout] = MMmulteda(Y,varargin)
%MMmulteda computes MM estimators in multivariate analysis for a series of values of eff
%
%
%<a href="matlab: docsearchFS('MMmulteda')">Link to the help function</a>
%
%  Required input arguments:
%
% Y :           Input data. Matrix.
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
%  Optional input arguments:
%
%      InitialEst : starting values of the MM-estimator. [] (default) or structure.
%                  InitialEst must contain the following fields:
%                   InitialEst.loc0 =  1 x v vector (estimate of the centroid)
%                   InitialEst.shape0 = v x v matrix (estimate of the shape matrix)
%                   InitialEst.auxscale = scalar (estimate of the scale parameter).
%                  If InitialEst is empty (default)
%                  program uses S estimators. In this last case it is
%                   possible to specify the options given in function Smult.
%               Example - 'InitialEst',[]
%               Data Types - struct
%  Soptions  :  options if initial estimator is S and InitialEst is empty.
%               Srhofunc,Snsamp,Srefsteps, Sreftol, Srefstepsbestr,
%               Sreftolbestr, Sminsctol, Sbestr.
%               See function Smult.m for more details on these options.
%               It is necessary to add to the S options the letter
%               S at the beginning. For example, if you want to use the
%               optimal rho function the supplied option is
%               'Srhofunc','optimal'. For example, if you want to use 3000
%               subsets, the supplied option is 'Snsamp',3000
%               Example - 'Snsamp',1000
%               Data Types - single | double
%         Sbdp :  breakdown point. Scalar.
%               It measures the fraction of outliers
%               the algorithm should resist. In this case any value greater
%               than 0 but smaller or equal than 0.5 will do fine (default=0.5).
%               Note that given bdp nominal
%               efficiency is automatically determined.
%                 Example - 'Sbdp',0.4
%                 Data Types - double
%      Sbestr  : number of "best betas" to remember. Scalar. Scalar defining number of "best betas" to remember from the
%               subsamples. These will be later iterated until convergence
%               (default=5)
%                 Example - 'Sbestr',10
%                 Data Types - single | double
%     Sminsctol: tolerance for the iterative
%               procedure for finding the minimum value of the scale. Scalar.
%               Value of tolerance for the iterative
%               procedure for finding the minimum value of the scale
%               for each subset and each of the best subsets
%               (It is used by subroutine minscale.m)
%               The default value is 1e-7;
%                 Example - 'Sminsctol',1e-7
%                 Data Types - single | double
%       Snsamp   : Number of subsamples which will be extracted to find the
%                 robust estimator. Scalar. If nsamp=0 all subsets will be extracted.
%                 They will be (n choose p).
%                 If the number of all possible subset is <1000 the
%                 default is to extract all subsets otherwise just 1000.
%                 Example - 'Snsamp',1000
%                 Data Types - single | double
%    Srefsteps : Number of refining iterations. Scalar. Number of refining iterationsin each
%               subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%                 Example - 'Srefsteps',0
%                 Data Types - single | double
%     Sreftol  : scalar. Default value of tolerance for the refining steps.
%               The default value is 1e-6;
%                 Example - 'Sreftol',1e-8
%                 Data Types - single | double
%Srefstepsbestr: number of refining iterations for each best subset. Scalar.
%               Scalar defining number of refining iterations for each
%               best subset (default = 50).
%                 Example - 'Srefstepsbestr',10
%                 Data Types - single | double
% Sreftolbestr : Tolerance for the refining steps. Scalar.
%               Tolerance for the refining steps
%               for each of the best subsets
%               The default value is 1e-8;
%                 Example - 'Sreftolbestr',1e-10
%                 Data Types - single | double
%
%               MM options
%
%
%      eff     : nominal efficiency. Scalar or vector.
%                Vector defining nominal efficiency (i.e. a series of numbers between
%                 0.5 and 0.99). The default value is the sequence 0.5:0.01:0.99
%                 Asymptotic nominal efficiency is:
%                 $(\int \psi' d\Phi)^2 / (\psi^2 d\Phi)$
%                 Example - 'eff',[0.85 0.90 0.95 0.99]
%                 Data Types - double
%     effshape : location or scale effiicency. dummy scalar. 
%                If effshape=1 efficiency refers to shape 
%                efficiency else (default) efficiency refers to location
%                 Example - 'effshape',1
%                 Data Types - double
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
%       plots : Plot on the screen. Scalar or structure.
%               If plots is a structure or scalar equal to 1, generates:
%               (1) a plot of Mahalanobis distances against index number. The
%               confidence level used to draw the confidence bands for
%               the MD is given by the input option conflev. If conflev is
%               not specified a nominal 0.975 confidence interval will be
%               used.
%               (2) a scatter plot matrix with the outliers highlighted.
%               If plots is a structure it may contain the following fields
%                   plots.labeladd = if this option is '1', the outliers in the
%                       spm are labelled with their unit row index. The
%                       default value is labeladd='', i.e. no label is
%                       added.
%                   plots.nameY = cell array of strings containing the labels of
%                       the variables. As default value, the labels which
%                       are added are Y1, ...Yv.
%                 Example - 'plots',0
%                 Data Types - single | double
%       nocheck : Check input arguments. Scalar. If nocheck is equal to 1
%                 no check is performed on
%                 matrix Y. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%
% Output:
%
%  out :     A structure containing the following fields
%
%         out.Loc  = length(bdp)-by-v  matrix containing MM estimate of
%                   location for each value of eff
%         out.Shape= v-by-v-by-length(eff) 3D array  containing robust
%                   estimate of the shape for each value of eff. Remark det|shape|=1
%         out.Scale= length(eff) vector. Robust estimate of the scale for
%                    each value of eff
%         out.Cov  = v-by-v-by-length(eff) 3D array  containing robust estimate of
%                    Note that  out.scale(i)^2 * out.shape(:,:,i) = robust estimate of
%                   covariance matrix
%           out.bs = (v+1) x 1 vector containing the units forming best subset
%                    associated with MM estimate of location.
%           out.MAL= n x length(eff) matrix containing the estimates of the robust
%                       Mahalanobis distances (in squared units) for each
%                       value of eff
%     out.Outliers = n-by-length(bdp) matrix containing true for the outliers.
%                    It is a Boolean matrix containing the list of the
%                    units declared as outliers for each value of bdp using
%                    confidence level specified in input scalar conflev
%      out.Weights = n x length(eff) matrix containing the weights for each
%                       value of eff
%      out.conflev = Confidence level that was used to declare outliers
%      out.singsub = Number of subsets without full rank. Notice that
%                    out.singsub > 0.1*(number of subsamples) produces a
%                    warning
%           out.eff= vector which contains the values of efficiency which have
%                       been used
%            out.Y = Data matrix Y.
%        out.class = 'MMmulteda'
%
%  Optional Output:
%
%            C        : matrix containing the indices of the subsamples
%                       extracted for computing the estimate (the so called
%                       elemental sets).
%
% See also: Smulteda, MMmult
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), Robust Statistics, Theory
% and Methods, Wiley, New York.
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('MMmulteda')">Link to the help page for this function</a>
% Last modified 31-05-2016
%
%
% Examples:

%{
    %% MMmult with all default options.
    load('swiss_banknotes');
    Y=swiss_banknotes.data;
    Y=Y(1:100,:);
    [outMM]=MMmulteda(Y);
    malfwdplot(outMM);
%}

%{
    %% MMmult with optional arguments.
    Y = load('geyser.txt');
    [out1]=MMmulteda(Y,'conflev',0.99,'plots',1);
%}

%{
    % MMmulteda with exctracted subsamples.
    load('swiss_banknotes');
    Y=swiss_banknotes.data;
    Y=Y(1:100,:);
    [outMM,C]=MMmulteda(Y);
%}

%% Beginning of code

%chkinputM does not do any check if option nocheck=1
nnargin=nargin;
vvarargin=varargin;
[Y,n,v] = chkinputM(Y,nnargin,vvarargin);

seq=1:n;
% default values for the initial S estimate:

% default value of break down point
Sbdpdef=0.5;
% default values of subsamples to extract
Snsampdef=1000;
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

% default values of nominal efficiency which are used
eff=0.5:0.01:0.99;


options=struct('InitialEst','','Snsamp',Snsampdef,'Srefsteps',Srefstepsdef,...
    'Sbestr',Sbestrdef,'Sreftol',Sreftoldef,'Sminsctol',Sminsctoldef,...
    'Srefstepsbestr',Srefstepsbestrdef,'Sreftolbestr',Sreftolbestrdef,...
    'Sbdp',Sbdpdef,...
    'nocheck',0,'eff',eff,'effshape',0,'refsteps',100,'tol',1e-7,...
    'conflev',0.975,'plots',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:MMmulteda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

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
    
    % first compute S-estimator with a fixed breakdown point
    if nargout==2
        [Sresult , C] = Smult(Y,'nsamp',nsamp,'bdp',bdp,'refsteps',refsteps,'bestr',bestr,...
            'reftol',reftol,'minsctol',minsctol,'refstepsbestr',refstepsbestr,...
            'reftolbestr',reftolbestr,...
            'nocheck',1);
        varargout = {C};
    else
        Sresult = Smult(Y,'nsamp',nsamp,'bdp',bdp,'refsteps',refsteps,'bestr',bestr,...
            'reftol',reftol,'minsctol',minsctol,'refstepsbestr',refstepsbestr,...
            'reftolbestr',reftolbestr,...
            'nocheck',1);
    end
    
    auxscale = Sresult.scale;
    shape0 = Sresult.shape;
    loc0 = Sresult.loc;
else
    auxscale = InitialEst.auxscale;
    shape0 = InitialEst.shape0;
    loc0 = InitialEst.loc0;
end


%% MM part

% Asymptotic nominal efficiency (for location or shape)
eff=options.eff;

% effshape = scalar which specifies whether nominal efficiency refers to
% location or scale
effshape=options.effshape;

% refsteps = maximum number of iterations in the MM step
refsteps=options.refsteps;

% tol = tolerance to declare convergence in the MM step
tol = options.tol;
conflev=options.conflev;

% Initialize quantities to store for each value of eff
leff=length(eff);
Weights=zeros(n,leff);
MAL=zeros(n,leff);
Outliers=false(n,leff);
% Loc= matrix which will contain location coefficients
Loc=zeros(leff,v);
% Covar = 3D array which will contain the estimate of the covariance matrix
Covar=zeros(v,v,leff);
% Shape = 3D array which will contain the estimate of the shape matrix
Shape=zeros(v,v,leff);


for jj=1:leff
    % MMmultcore = function which does IRWLS steps from initial loc (loc0) and
    % initial shape matrix (Shape0). The estimate of sigma (auxscale) remains
    % fixed inside this routine
    outIRW = MMmultcore(Y,loc0,shape0,auxscale,'eff',eff(jj),'effshape',effshape,'refsteps',refsteps,'reftol',tol,'conflev',conflev,'nocheck',1);
    
    
    md = mahalFS(Y,outIRW.loc,outIRW.cov);
    
    outliers = seq(md > chi2inv(conflev,v) );
    
    MAL(:,jj)=md;
    Loc(jj,:)=outIRW.loc;
    Shape(:,:,jj)=outIRW.shape;
    Covar(:,:,jj)=outIRW.cov;
    Weights(:,jj)=outIRW.weights;
    Outliers(outliers,jj)=true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.Loc     = Loc;         % robust estimate of location
out.Shape   = Shape;       % robust estimate of shape matrix
out.cov     = Covar;       % robust estimate of covariance matrix
out.Weights = Weights;
out.MAL=MAL;
out.Y = Y;
out.eff=eff;
out.scale=auxscale;
out.class   = 'MMmulteda';

% Store in output structure the outliers found with confidence level conflev
out.conflev = conflev;

plo=options.plots;

% Plot residuals as function of the break down point
if plo==1
    laby='MM Mahalanobis distances';
    malfwdplot(out);
    ylabel(laby)
end

end
%FScategory:MULT-Multivariate




