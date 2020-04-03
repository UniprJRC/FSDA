function [reduced_est, reduced_model, msgstr] = LTStsVarSel(y,varargin)
%LTStsVarSel does variable selection in the robust time series model LTSts
%
%<a href="matlab: docsearchFS('LTStsVarSel')">Link to the help function</a>
%
% LTSts requires variable selection when the optimal model parameters are
% not known in advance. This happens in particular when the function has to
% be applied to many heterogeneous datasets in an automatic way, possibly
% on a regular basis (that is the model parameters are expected to change
% over time even for datasets associated to the same phenomenon).
%
% The approach consists in iteratively eliminating the less significant
% estimated model parameter, starting from an over-parametrized model. The
% model parameters are re-estimated with LTSts at each step, until all the
% p-values are below a given threshold. Then, the output is a reduced time
% series model with significant parameters only.
%
%
%  Required input arguments:
%
%    y:         Time series to analyze. Vector. A row or a column vector
%               with T elements, which contains the time series.
%
%
%  Optional input arguments:
%
%    model:     model type. Structure. A structure which specifies the
%               (over-parametrized) model which will be used to initialise
%               the variable selection process. The model structure is
%               identical to the one defined for function LTSts: for
%               convenience, we list the fields also here:
%
%               model.s = scalar (length of seasonal period). For monthly
%                         data s=12 (default), for quartely data s=4, ...
%               model.trend = scalar (order of the trend component).
%                       trend = 0 implies no trend
%                       trend = 1 implies linear trend with intercept (default),
%                       trend = 2 implies quadratic trend
%                       trend = 3 implies cubic trend
%                       Admissible values for trend are, 0, 1, 2 and 3.
%                       In the paper RPRH to denote the order of the trend
%                       symbol A is used. If this field is not present into
%                       input structure model, model.trend=2 is used.
%               model.seasonal = scalar (integer specifying number of
%                        frequencies, i.e. harmonics, in the seasonal
%                        component. Possible values for seasonal are
%                        $1, 2, ..., [s/2]$, where $[s/2]=floor(s/2)$.
%                        For example:
%                        if seasonal =1 (default) we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 sin ( 2 \pi t/s)$;
%                        if seasonal =2 we have:
%                        $\beta_1 \cos( 2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s)
%                        + \beta_3 \cos(4 \pi t/s) + \beta_4 \sin (4 \pi t/s)$.
%                        Note that when $s$ is even the sine term disappears
%                        for $j=s/2$ and so the maximum number of
%                        trigonometric parameters is $s-1$.
%                        If seasonal is a number greater than 100 then it
%                        is possible to specify how the seasonal component
%                        grows over time.
%                        For example, seasonal = 101 implies a seasonal
%                        component which just uses one frequency
%                        which grows linearly over time as follows:
%                        $(1+\beta_3 t)\times ( \beta_1 cos( 2 \pi t/s) +
%                        \beta_2 \sin ( 2 \pi t/s))$.
%                        For example, seasonal =201 implies a seasonal
%                        component which just uses one frequency
%                        which grows in a quadratic way over time as
%                        follows:
%                        $(1+\beta_3 t + \beta_4  t^2)\times( \beta_1 \cos(
%                        2 \pi t/s) + \beta_2 \sin ( 2 \pi t/s))$.
%                        seasonal =0 implies a non seasonal model.
%                       In the paper RPRH to denote the number of
%                       frequencies of the seasonal component
%                       symbol B is used, while symbol G is used to denote
%                       the order of the trend of the seasonal component.
%                       Therefore, for example, model.seasonal=201
%                       corresponds to B=1 and G=2, while model.seasonal=3
%                       corresponds to B=3 and G=0;
%               model.X  =  matrix of size T-by-nexpl containing the
%                         values of nexpl extra covariates which are likely
%                         to affect y.
%               model.lshift = scalar greater or equal than 0 which
%                         specifies whether it is necessary to include a
%                         level shift component. lshift = 0 (default)
%                         implies no level shift component. If lshift is an
%                         interger greater then 0 then it is possible to
%                         specify the moment to start considering level
%                         shifts. For example if lshift =13 then the
%                         following additional parameters are estimated
%                          $\beta_{LS1}* I(t \geq beta_{LS2})$ where $\beta_{LS1}$
%                          is a real number and $\beta_{LS2}$ is an integer
%                          which assumes values 14, 14, ..., T-13.
%                         In general, the level shift which are considered
%                         are referred to times (lshift+1):(T-lshift).
%                       In the paper RPRH $\beta_{LS1}$ is denoted with
%                       symbol $\delta_1$, while, $\beta_{LS2}$ is denoted
%                       with symbol $\delta_2$.
%               model.ARp = scalar greater or equal than 0 which
%                         specifies the length of the autoregressive
%                         component. The default value of model.ARp is 0,
%                         that is there is no autoregressive component.
%                 Example - 'model', model
%                 Data Types - struct
%               Remark: the default overparametrized model is for monthly
%               data with a quadratic
%               trend (3 parameters) + seasonal component with just one
%               harmonic (2 parameters), no additional explanatory
%               variables, no level shift and no AR component that is
%                               model=struct;
%                               model.s=12;
%                               model.trend=1;
%                               model.seasonal=1;
%                               model.X='';
%                               model.lshift=0;
%                               model.ARp=0;
%               Using the notation of the paper RPRH we have A=1, B=1; and
%               $\delta_1=0$.
%
%    thPval:    threshold for pvalues. Scalar. A value between 0 and 1.
%               An estimated parameter/variable is eliminated if the
%               associated pvalue is below thPval. Default is thPval=0.01.
%                 Example - 'thPval',0.05
%                 Data Types - double
%
%    plots:     Plots on the screen. Scalar.
%               If plots = 1, the typical LTSts plots will be shown on the
%               screen. The default value of plot is 0 i.e. no plot is
%               shown on the screen.
%                 Example - 'plots',1
%                 Data Types - double
%
%    msg:       Messages on the screen. Scalar.
%               Scalar which controls whether LTSts will display or not
%               messages on the screen. Deafault is msg=0, that is no
%               messages are displayed on the screen. If msg==1 messages
%               displayed on the screen are about estimated time to compute
%               the estimator and the warnings about
%               'MATLAB:rankDeficientMatrix', 'MATLAB:singularMatrix' and
%               'MATLAB:nearlySingularMatrix'.
%               Example - 'msg',1
%               Data Types - double
%
%  dispresults : Display results of final fit. Boolean. If dispresults is
%               true, labels of coefficients, estimated coefficients,
%               standard errors, tstat and p-values are shown on the
%               screen in a fully formatted way. The default value of
%               dispresults is false.
%               Example - 'dispresults',true
%               Data Types - logical
%
%
%  Output:
%
%  reduced_est:  A reduced model structure obtained by eliminating
%                    parameters that are non-significant. It is a structure
%                    containing the typical input model fields for function
%                    LTSts (refer to LTSts for details):
%                    model.s = the optimal length of seasonal period.
%                    model.trend = the optimal order of the trend.
%                    model.seasonal = the optimal number of frequencies in
%                      the seasonal component.
%                    model.lshift = the optimal level shift position.
%                    model.X = a matrix containing the values of the extra
%                      covariates which are likely to affect y. If the
%                      imput model specifies autoregressive components
%                      in model.ARp, then the selected ones will be also
%                      included in model.X.
%
% reduced_model:  Structure containing the output fields of the optimal model.
%                    The fields are those of function LTSts (refer to
%                    LTSts for details):
%                    out.B = matrix of estimated beta coefficients.
%                    out.h = number of observations that have determined
%                      the initial LTS estimator.
%                    out.bs = vector of the units with the smallest
%                      squared residuals before the reweighting step.
%                    out.Hsubset = matrix of the units forming best H
%                      subset for each tentative level shift considered.
%                    out.numscale2 = matrix of the values of the lts.bestr
%                      smallest values of the target function.
%                    out.BestIndexes = matrix of indexes associated with
%                      the best nbestindexes solutions.
%                    out.Likloc = matrix containing local sum of squares
%                      of residuals determining the best position of level
%                      shift.
%                    out.RES = matrix containing scaled residuals for all
%                      the units of the original time series monitored in
%                      steps lshift+1, lshift+2, ....
%                    out.yhat = vector of fitted values after final step.
%                    out.residuals = vector of scaled residuals from
%                      after final NLS step.
%                    out.weights = vector of weights after adaptive
%                      reweighting.
%                    out.scale = final scale estimate of the residuals
%                      using final weights.
%                    out.conflev = confidence level used to declare outliers.
%                    out.outliers = vector of the units declared outliers.
%                    out.singsub = number of subsets wihtout full rank.
%                    out.y = response vector y.
%                    out.X = data matrix X containing trend, seasonal, expl
%                       (with autoregressive component) and lshift.
%                    out.class = 'LTSts'.
%
%
%  Optional Output:
%
%         msgstr     : String containing the last warning message.
%                      This relates to the execution of the LTS.
%
% See also LTSts
%
% References:
%
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('LTStsVarSel')">Link to the help function</a>
%
%$LastChangedDate:: 2019-08-31 00:40:12 #$: Date of the last commit
%
% Examples:
%
%
%{
    % run LTStsVarSel with all default options.
    % data model
    model=struct;
    model.trend=1;                  % linear trend
    model.trendb=[0 1];             % parameters of the linear trend
    model.s=12;                     % monthly time series
    model.seasonal=1;               % 1 harmonic with linear trend
    model.seasonalb=[10 10];        % parameter for one harmonic with linear trend
    model.lshiftb=100;              % level shift amplitude
    model.lshift= 30;               % level shift amplitude
    model.signal2noiseratio = 100;  % signal to noise
    
    rng('default')
    n = 100;                        % sample size
    tmp = rand(n,1);
    model.X = tmp.*[1:n]';          % a extra covariate
    model.Xb = 1;                   % beta coefficient of the covariate
    % generate data
    out_sim=simulateTS(n,'plots',1,'model',model);

    %run LTStsVarSel with all default options
    rng(1);
    [out_model_0, out_reduced_0] = LTStsVarSel(out_sim.y);
 
    % optional: add a FS step to the LTSts estimator
    % outFS = FSRts(out_sim.y,'model',out_model_0);
    % To be fixed: 'Non existent user option found-> '    'ARp'
%}

%{
    % run LTStsVarSel starting from a specific over-parametrized model.

    % sample size
    n = 100;                       
    tmp = rand(n,1);
    model.X = tmp.*[1:n]';          % a extra covariate
    model.Xb = 1;                   % beta coefficient of the covariate
    out_sim=simulateTS(n,'plots',1,'model',model);
    % complete model to be tested.
    overmodel=struct;
    overmodel.trend=2;              % quadratic trend
    overmodel.s=12;                 % monthly time series
    overmodel.seasonal=303;         % number of harmonics
    overmodel.lshift=4;             % position where to start monitoring level shift
    overmodel.X=tmp.*[1:n]';

    % pval threshold
    thPval=0.01;

    [out_model_1, out_reduced_1] = LTStsVarSel(out_sim.y,'model',overmodel,'thPval',thPval,'plots',1);

%}

%{
    % run LTStsVarSel starting from over-parametrized model with autoregressive components.
    % add three autoregressive components to the complete model.
     n = 100;                        % sample size
    tmp = rand(n,1);
    model.X = tmp.*[1:n]';          % a extra covariate
    model.Xb = 1;                   % beta coefficient of the covariate
    out_sim=simulateTS(n,'plots',1,'model',model);
    % complete model to be tested.
    overmodel=struct;
    overmodel.trend=2;              % quadratic trend
    overmodel.s=12;                 % monthly time series
    overmodel.seasonal=303;         % number of harmonics
    overmodel.lshift=4;             % position where to start monitoring level shift
    overmodel.X=tmp.*[1:n]';
    overmodel.ARp=3;

    % pval threshold
    thPval=0.01;
     
    [out_model_2, out_reduced_2] = LTStsVarSel(out_sim.y,'model',overmodel,'thPval',thPval);
%}

%{
    % run LTStsVarSel with default options and return warning messages.
        % data model
    model=struct;
    model.trend=1;                  % linear trend
    model.trendb=[0 1];             % parameters of the linear trend
    model.s=12;                     % monthly time series
    model.seasonal=1;               % 1 harmonic with linear trend
    model.seasonalb=[10 10];        % parameter for one harmonic with linear trend
    model.lshiftb=100;              % level shift amplitude
    model.lshift= 30;               % level shift amplitude
    model.signal2noiseratio = 100;  % signal to noise
    
    n = 100;                        % sample size
    tmp = rand(n,1);
    model.X = tmp.*[1:n]';          % a extra covariate
    model.Xb = 1;                   % beta coefficient of the covariate
    % generate data
    out_sim=simulateTS(n,'plots',1,'model',model);
    [out_model_3, out_reduced_3, messages] = LTStsVarSel(out_sim.y);
%}


%% Beginning of code 

% Input parameters checking

warning('off','all');

if nargin<1
    error('FSDA:LTStsVarSel:MissingInputs','Input time series is missing');
end
% if nargin<2
%     error('FSDA:LTStsVarSel:MissingInputs','Provide an initial (over-parametrised) model');
% end

% Set up defaults for the over-parametrized model
modeldef          = struct;
modeldef.trend    = 2;        % quadratic trend
modeldef.s        = 12;       % monthly time series
modeldef.seasonal = 303;      % three harmonics growing cubically (B=3, G=3)
modeldef.X        = [];       % no extra explanatory variable
modeldef.lshift   = 0;        % no level shift
modeldef.ARp      = 0;        % no autoregressive component

options=struct('model',modeldef, 'thPval', 0.01, ...
    'plots',0,'msg',0,'dispresults',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:LTSts3:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all optional arguments were present in structure options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:LTSts:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% Put User options inside modeldef
if ~isequal(options.model,modeldef)
    fld=fieldnames(options.model);
    chkoptions(modeldef,fld)
    for i=1:length(fld)
        modeldef.(fld{i})=options.model.(fld{i});
    end
end
% and finally set the over-parametrized model to start with
model = modeldef;

thPval = options.thPval;
plots  = options.plots;
msg    = options.msg;
dispresults = options.dispresults;

%% Step 1: estimate model parameters with LTSts for the over-parametrized input model

n   = size(y,1);     % number of observations in the input dataset
h1  = round(n*0.9);  % default for h (num. obs. for the LTS estimator)

% Estimate the parameters based on initial full model
out_LTSts = LTSts(y,'model',model,'nsamp',500,'h',h1,...
    'plots',plots,'msg',msg,'dispresults',dispresults,'SmallSampleCor',1);
if plots
    a=gcf;
    if isfield(model,'X')
        title(a.Children(end),['trend =' num2str(model.trend) ', seas = ' num2str(model.seasonal) ' every ' num2str(model.s) ' months , X = ' num2str(size(model.X,2)) ]);
    else
        title(a.Children(end),['trend =' num2str(model.trend) ', seas = ' num2str(model.seasonal) ' every ' num2str(model.s) ' months']);
    end
end

%% Step 2: iterate and reduce the model

% Step (2a) tests the model parameters. It consists in identifying the
% largest p-value among the parameters of:
% - level shift,
% - harmonics,
% - covariates,
% - the largest degree component of:
%   * trend,
%   * amplitude of harmonics.
%
% Step (2b) re-estimates the model. It remove the less significant parameter
% and re-estimates the model with LTSts. Remark: the level shift position
% is not re-estimated, but the original estimate is kept as is.


% Initializations

% AllPvalSig = dycotomic variable that becomes 1 if all variables are
% significant. This means that the variables are to be kept and the
% iterative procedure should stop.
AllPvalSig=0;


% iniloop associated with level shift, which has to be tested only once
iniloop=1;

% fwd search index
lLSH = n-model.lshift*2;

% initialize a flag to check the initial presence of level shift
lshift_present = 0;

% Iterative model reduction.
while AllPvalSig == 0
    % The loop terminates when all p-values are smaller than thPval. In
    % this case AllPvalSig will become equal to 1
    
    rownam=out_LTSts.Btable.Properties.RowNames;
    seqp=1:length(rownam);
    
    % Position of the last element of the trend component
    posLastTrend=max(seqp(contains(rownam,'b_trend')));
    if posLastTrend>1
        LastTrendPval=out_LTSts.Btable{posLastTrend,'pval'};
    else
        LastTrendPval=0;
    end
    
    posX=seqp(contains(rownam,'b_X'));
    if ~isempty(posX)
        % if iniloop is 0 the pval of the last expl variable is in reality
        % the pval of the level shift component and therefore it has to be
        % recalibrated
        PvalX=out_LTSts.Btable{posX,'pval'};
        if iniloop==0
            tstatX=out_LTSts.Btable{posX,'t'};
            %lLSH: fwd search index
            %abs(tstatX(end)): minimum deletion residuals
            %size(out_LTSts.Btable,1)-1: number of explanatory variables
            lsdet=FSRinvmdr([lLSH abs(tstatX(end))],size(out_LTSts.Btable,1)-1);
            %lsdet(1,2) = confidence level of each value of mdr.
            PvalX(end)=1-lsdet(1,2);
        end
        
        [maxPvalX,posmaxPvalX]=max(PvalX);
        %posmaxPvalX=posX(posmaxPvalX);
    else
        maxPvalX=0;
        posmaxPvalX=[];
    end
    
    % tre=cellfun(@isempty,strfind(rownam,'b_varamp'));
    posLastVarAmpl=max(seqp(contains(rownam,'b_varamp')));
    
    if ~isempty(posLastVarAmpl)
        LastVarAmplPval=out_LTSts.Btable{posLastVarAmpl,'pval'};
    else
        LastVarAmplPval=0;
    end
    
    % delete first (if necessary the time varying harmonic rather than the
    % unique harmonic)
    if LastVarAmplPval>0
        LastHarmonicPval=out_LTSts.LastHarmonicPval;
        sea=(num2str(model.seasonal));
        if strcmp(sea(end),'1')
            LastHarmonicPval=0;
        else
        end
    else
        LastHarmonicPval=out_LTSts.LastHarmonicPval;
    end
    
    if model.lshift>0
        LevelShiftPval=out_LTSts.LevelShiftPval;
        posLS=out_LTSts.posLS;
    else
        LevelShiftPval=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%

    posAR=seqp(contains(rownam,'b_AR'));
    % Initialize pvalue of last AR component to be zero.
    LastARPval=0;
    if ~isempty(posAR)
        % Position of the last element of the AR component
        posLastAR=max(seqp(contains(rownam,'b_AR')));
        if posLastAR>1
            LastARPval=out_LTSts.Btable{posLastAR,'pval'};
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%
    
    % Group all p-values into a vector
    Pvalall = [LastTrendPval;...
        LastHarmonicPval;...
        maxPvalX; LastVarAmplPval;...
        LevelShiftPval; LastARPval];
    
    % Start of model reduction (step 2b)
    [maxPvalall,indmaxPvalall]=max(Pvalall);
    
    if maxPvalall>thPval
        switch indmaxPvalall
            case 1
                % if indmaxPvalall == 1
                % Remove from model the last term of the trend component
                if msg==1 || plots==1
                    removed =['Removing trend of order ' num2str(model.trend)];
                end
                model.trend=model.trend-1;
            case 2
                % elseif indmaxPvalall ==2
                % Remove from model the last term of the seaonal component that
                % is remove last harmonic
                if msg==1 || plots==1
                    tmp = num2str(model.seasonal);
                    tmp = tmp(end);
                    removed =['Removing harmonic number ' tmp];
                end
                model.seasonal= model.seasonal-1;
            case 3
                % elseif indmaxPvalall ==3
                % Remove from model the non signif expl var
                if posmaxPvalX==size(model.X,2) && iniloop==0
                    if msg==1 || plots==1
                        removed ='Removing level shift component';
                    end
                else
                    if msg==1 || plots==1
                        removed =['Removing expl. variable number ' num2str(posmaxPvalX)];
                    end
                end
                model.X(:,posmaxPvalX)= [];
            case 4
                % elseif indmaxPvalall ==4
                % Remove from model the high order term of non linear seasonality
                if msg==1 || plots==1
                    strseaso=num2str(model.seasonal);
                    removed = ['Removing amplitude of order ' strseaso(1) ' of seas. comp.'];
                end
                model.seasonal= model.seasonal-100;
            case 5
                % elseif indmaxPvalall ==5
                % Remove from model the level shift component
                model.lshift=0;
                if msg==1 || plots==1
                    removed ='Remove level shift component';
                end
            case 6
                if msg==1 || plots==1
                    strAR=num2str(model.ARp);
                    removed = ['Removing AR component of order ' strAR];
                end
                model.ARp=model.ARp-1;
            otherwise
                %else
        end
        
        if msg==1 || plots==1
            disp(removed)
        end
        % keep the level shift component and transfer it to X part
        if model.lshift>0 && iniloop==1
            Xls=[zeros(posLS-1,1); ones(n-posLS+1,1)];
            if ~isfield(model,'X')
                model.X=[];
            end
            model.X=[model.X Xls];
            lshift_present = posLS;
            model.lshift=0;
            iniloop=0;
        end
        
        % Re-run the model but do not re-estimate the position of the
        % level shift
        [out_LTSts]=LTSts(out_LTSts.y,'model',model,'nsamp',100,...
            'plots',plots,'msg',msg,'dispresults',dispresults,'h',h1,'SmallSampleCor',1);
        
        if plots==1
            a=gcf;
            title(a.Children(end),{['trend = ' num2str(model.trend) ...
                ', seas = ' num2str(model.seasonal) ...
                ' every ' num2str(model.s) ...
                ' months', ', LS = ' num2str(lshift_present), ...
                ', X = ' num2str(size(model.X,2)-1) ],num2str(removed)});
            if lshift_present>0
                hold on;
                line([lshift_present lshift_present],[min(out_LTSts.y),max(out_LTSts.y)],...
                    'Color','black','LineStyle','--','Linewidth',1);
            end
        end
    else
        % All the variables are significant, variable selection procedure
        % stops
        AllPvalSig=1;
    end
    
end

% If level shift is present and model.X is not empty, it means that the last
% column of model.X is the level shift explanatory variable
if lshift_present > 0 && ~isempty(model.X)
    tmp = find(model.X(:,end)>0);
    model.lshift=tmp(1);
    out_LTSts.posLS=tmp(1);
    model.X(:,end) = [];
    out_LTSts.Btable.Properties.RowNames(end)={'b_lshift'};
end

reduced_est   = model;
reduced_model = out_LTSts;

if msg==1
    disp('The final select model has these parameters:')
    disp(reduced_model);
end

[msgstr, ~] = lastwarn;

end
%FScategory:REG-Regression