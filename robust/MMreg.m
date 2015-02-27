function [out , varargout] = MMreg(y,X,varargin)
%MMreg computes MM estimator of regression coefficients
%
%<a href="matlab: docsearchFS('mmreg')">Link to the help function</a>
%
%  Required input arguments:
%
%    y:         A vector with n elements that contains the response variable.
%               It can be both a row of column vector.
%    X :        Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations,
%               and columns represent variables. 
%
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%
%  Optional input arguments:
%
%   intercept : If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%  InitialEst : a structure containing starting values of the MM-estimator.
%               The structure must contain
%               - loc =  v x 1 vector (estimate of the centroid)
%               - scale = scalar (estimate of the scale parameter).
%               If InitialEst is empty (default)
%               program uses S estimators. In this last case it is
%               possible to specify the options given in function Sreg.
%
%               Soptions (if InitialEst is empty): see function Sreg.
%               Remark: it is necessary to add to the S options the letter
%               S at the beginning. For example, if you want to use the
%               optimal rho function the supplied option is
%               'Srhofunc','optimal'. For example, if you want to use 3000
%               subsets, the supplied option is 'Snsamp',3000
%                   
%
%               MM options
%
%      eff     : scalar defining nominal efficiency (i.e. a number between
%                 0.5 and 0.99). The default value is 0.95
%                 Asymptotic nominal efficiency is:
%                 (\int \psi' d\Phi)^2 / (\psi^2 d\Phi)
%     effshape : dummy scalar. If effshape=1 efficiency refers to shape 
%                efficiency else (default) efficiency refers to location
%     refsteps  : scalar defining maximum number of iterations in the MM
%                 loop. Default value is 100.
%       tol    : scalar controlling tolerance in the MM loop.
%                 Default value is 1e-7
%     conflev:  Scalar between 0 and 1 containing the confidence level
%               used to declare units as outliers. Usually conflev = 0.95,
%               0.975, 0.99 (individual alpha) or 1-0.05/n, 1-0.025/n,
%               1-0.01/n (simultaneous alpha). Default value is 0.975.
%      nocheck: Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0. 
%       plots : Scalar or structure.
%               If plots = 1, generates a plot of scaled residuals against
%               index number. The confidence level used to draw the
%               confidence bands for the scaled residuals is given by the
%               input option conflev. If conflev is not specified a nominal
%               0.975 confidence interval will be used.%
%       yxsave : scalar that is set to 1 to request that the response 
%                vector y and data matrix X are saved into the output
%                structure out. Default is 0, i.e. no saving is done.
%
%  Output:
%
%
%  The output consists of a structure 'out' containing the following fields:
%       out.beta        :   p x 1 vector containing MM estimate of 
%                           regression coefficients
%       out.auxscale    :   scalar, S estimate of the scale (or supplied
%                           external estimate of scale, if option InitialEst  
%                           is not empty)
%       out.residuals	:   n x 1 vector containing standardized MM
%                           residuals
%                           out.residuals=(y-X*out.beta)/out.auxscale
%       out.weights     :   n x 1 vector. Weights assigned to each observation
%       out.Sbeta       :   p x 1 vector containing S estimate of regression
%                           coefficients (or supplied initial external
%                           estimate of regression coefficients, if option
%                           InitialEst is not empty)
%       out.Ssingsub    :   Number of subsets without full rank in the S 
%                           preliminary part. Notice that 
%                           out.singsub > 0.1*(number of subsamples) 
%                           produces a warning
%       out.outliers    :   1 x k vectors containing the outliers which
%                           have been found
%       out.conflev     :   Confidence level that was used to declare outliers
%       out.class       :   'MM'
%           out.rhofunc :   string identifying the rho function which has been
%                           used
%      out.rhofuncparam :   vector which contains the additional parameters
%                           for the specified rho function which have been
%                           used. For hyperbolic rho function the value of
%                           k =sup CVC. For Hampel rho function the parameters
%                           a, b and c
%            out.y      :   response vector Y. The field is present if option 
%                           yxsave is set to 1.
%            out.X      :   data matrix X. The field is present if option 
%                           yxsave is set to 1.
%  Optional Output:
%
%            C     : matrix of the indices of the samples extracted for
%                    computing the estimate
%
% References:
%
% ``Robust Statistics, Theory and Methods'' by Maronna, Martin and Yohai;
% Wiley 2006.
%
% Acknowledgements
%
% This function follows the lines of MATLAB/R code developed during the
% years by many authors.
% For more details see http://www.econ.kuleuven.be/public/NDBAE06/programs/
% and the R library robustbase http://robustbase.r-forge.r-project.org/
% The core of these routines, e.g. the resampling approach, however, has
% been completely redesigned, with considerable increase of the
% computational performance.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mmreg')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % MMreg with all default options
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
    % MMreg using the hyperbolic rho function 
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

%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[y,X,~,~] = chkinputR(y,X,nnargin,vvarargin);

% default values for the initial S estimate:

% default value of break down point
Sbdpdef=0.5;
% default values of subsamples to extract
Snsampdef=20;
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


options=struct('intercept',1,'InitialEst','','Snsamp',Snsampdef,'Srefsteps',Srefstepsdef,...
    'Sbestr',Sbestrdef,'Sreftol',Sreftoldef,'Sminsctol',Sminsctoldef,...
    'Srefstepsbestr',Srefstepsbestrdef,'Sreftolbestr',Sreftolbestrdef,...
    'Sbdp',Sbdpdef,'Srhofunc',Srhofuncdef,'Srhofuncparam','','nocheck',0,'eff',0.95,'effshape',0,...
    'refsteps',100,'tol',1e-7,'conflev',0.975,'plots',0,'yxsave',0);

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
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end


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
    
    rhofunc=options.Srhofunc;           % rho function which must be used
    rhofuncparam=options.Srhofuncparam;    % eventual additional parameters associated to the rho function
    
    
    % first compute S-estimator with a fixed breakdown point
    
    % SR is the routine which computes S estimates of beta and sigma in regression
    % Note that intercept is taken care of by chkinputR call.
    if nargout==2
        [Sresult , C] = Sreg(y,X,'nsamp',nsamp,'bdp',bdp,'refsteps',refsteps,'bestr',bestr,...
            'reftol',reftol,'minsctol',minsctol,'refstepsbestr',refstepsbestr,...
            'reftolbestr',reftolbestr,'rhofunc',rhofunc,'rhofuncparam',rhofuncparam,...
            'nocheck',1);

        varargout = {C};
    else
        Sresult = Sreg(y,X,'nsamp',nsamp,'bdp',bdp,'refsteps',refsteps,'bestr',bestr,...
            'reftol',reftol,'minsctol',minsctol,'refstepsbestr',refstepsbestr,...
            'reftolbestr',reftolbestr,'rhofunc',rhofunc,'rhofuncparam',rhofuncparam,...
            'nocheck',1);
    end
    
    bs = Sresult.beta;
    ss = Sresult.scale;
    singsub=Sresult.singsub;
else
    bs = InitialEst.beta;
    ss = InitialEst.scale;
    singsub=0;
end

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
    'refsteps',refsteps,'reftol',tol,'conflev',conflev,'plots',plots,'nocheck',1);


out = struct;
out.beta = outIRW.beta;
out.auxscale = ss;
out.residuals = (y-X*outIRW.beta)/ss; % MM scaled residuals
out.Sbeta = bs;
out.Ssingsub=singsub;
out.weights=outIRW.weights;
out.outliers=outIRW.outliers;
out.conflev=conflev;
out.class='MM';
out.rhofunc=rhofunc;
% In case of Hampel or hyperbolic tangent estimator store the additional
% parameters which have been used
% For Hampel store a vector of length 3 containing parameters a, b and c
% For hyperbolic store the value of k= sup CVC
if exist('rhofuncparam','var')
    out.rhofuncparam=rhofuncparam;
end

if options.yxsave
    if options.intercept==1;
        % Store X (without the column of ones if there is an intercept)
        out.X=X(:,2:end);
    else
        out.X=X;
    end
    % Store response
    out.y=y;
end

end
