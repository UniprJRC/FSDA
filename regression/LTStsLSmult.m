function out = LTStsLSmult(y,varargin)
%LTStsLSmult extends LTSts to the detection of multiple Level Shifts in time series
%
%<a href="matlab: docsearchFS('LTStsLSmult')">Link to the help function</a>
%
%  Required input arguments:
%
%        y  :   Time series to analyze. Vector. A row or a column vector
%               with T elements, which contains the time series.
%               Data Types - double.
%
%       Optional input arguments:
%
%     maxLS    : maximum number of Level Shifts. Scalar. Maximum number of Level Shifts that the function looks for. Default is 5.
%                Example - 'maxLS', 8
%                Data Types - double
%    startLS   : first position in which Level Shift should be tested.
%                   Scalar. It is the low extreme of the interval for
%                   the parameter model.lshift for the LTSts function.
%                   Default is 3. 
%                Example - 'startLS', 15
%                Data Types - double
%   sampleLS   : vector of the positions in which Level Shift should be
%               tested. Vector. It is be the vector for the parameter 
%                   model.lshift of the LTSts function.
%                   Default is 3. 
%                Example - 'sampleLS', [5 16 27 35]
%                Data Types - double
%   alphaLS    : maximum pvalue to consider a Level Shift significant. 
%               Scalar. It is a threshold for detecting significant Level Shifts. 
%               Default is 0.01. 
%                Example - 'alphaLS', 0.1
%                Data Types - double
%   thresLS    : minimum length of the time series to look for Level Shift.
%               Scalar. In order to detecting Level Shifts the series should 
%               have more than thresLS observations. 
%               Default is 20. 
%               Example - 'thresLS', 50
%   alphaLTS   : complement of the onfidence level. Scalar. Scalar between 
%               0 and 1 containing the complement (to 1) of the confidence 
%               level which is used to declare units as outliers. Usually 
%               alpha = 0.05/n, 0.025/n, 1-0.01/n (simultaneous alpha). 
%               Default value is 0.01/n.
%                 Example - 'alphaLTS',0.1
%                 Data Types - double
%   bdpLTS     : breakdown point. Scalar. It measures the fraction of outliers
%               the algorithm should resist. In this case any value greater
%               than 0 but smaller or equal than 0.5 will do fine. Default
%               value 0.1.
%                 Example - 'bdp',0.4
%                 Data Types - double
%        msg  : Messages on the screen. Boolean.
%               Scalar which controls whether to display or not messages
%               on the screen. If msg==true (default) messages are displayed on
%               the screen about estimated time to compute the estimator
%               and the warnings about 'MATLAB:rankDeficientMatrix',
%               'MATLAB:singularMatrix' and 'MATLAB:nearlySingularMatrix'
%               are set to off else no message is displayed on the screen
%               Example - 'msg',true
%               Data Types - logical
%       plots : Plots on the screen. Scalar.
%               If plots = 1, a two panel plot will be shown on the screen.
%               The upper panel contains the orginal time series with
%               fitted values. The bottom panel will contain the plot
%               of robust residuals against index number. The confidence
%               level which is used to draw the horizontal lines associated
%               with the bands for the residuals is specified in input
%               option conflev. 
%               The default value of plot is 0 i.e. no plot is shown on the
%               screen.
%                 Example - 'plots',1
%                 Data Types - double
%
%  Output:
%
%  out :     A structure containing the following fields
%            out.R2 = Scalar containing the goodness of fit.
%            out.yhat = Vector of the fitted values
%            out.LSpos = Vector of the positions of the detected Level
%                       Shifts
%            out.LSpval = Vector of the pvalues of the detected Level
%                       Shifts
%            out.LSud =  Vector of the sign (positive or negative) of the detected Level
%                       Shifts residuals
%            out.outliers = Vector of the positions of the detected
%                       outliers
%            out.outX = Structure containing the final output of LTSts. It
%               is useful later to run LTStsVarSel with in the parameter
%               model also the item model.X, containing the dummy variables with
%               Level Shifts positions. 
%
% See also LXS, wedgeplot
%
% References:
%
% Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
% Monitoring of Many Time Series with Application to Fraud Detection,
% "Econometrics and Statistics". [RPRH]
%
%
% Copyright 2008-2023.
% Written by Marco Riani, Domenico Perrotta, Peter
% Rousseeuw and Mia Hubert
%
%
%<a href="matlab: docsearchFS('LTSts')">Link to the help function</a>
%
%$LastChangedDate:: 2019-12-15 21:09:21 #$: Date of the last commit

% Examples:


%{
    %% A synthetic example.
    Y=rand(50,1);
    Y(35:end)=Y(35:end)+5;
    Y(7:20)=Y(7:20)-3;
    out  = LTStsLSmult(Y,'msg',true,'plots',1);
%}

%{
    % Trade data examples.
    % Two examples taken from the (extended version of the) series in:
    % Rousseeuw, P.J., Perrotta D., Riani M. and Hubert, M. (2018), Robust
    % Monitoring of Many Time Series with Application to Fraud Detection,
    % "Econometrics and Statistics".

    load('P_12119085_TQV_full.txt'); % KE-GB
    load('P_17049075_TQV_full.txt'); % UA-LT
    
    dates1 = P_12119085_TQV_full(:,1);
    yin1   = P_12119085_TQV_full(:,2);
    
    dates2 = P_17049075_TQV_full(:,1);
    yin2   = P_17049075_TQV_full(:,2);
    
    out1  = LTStsLSmult(yin1,'msg',true,'plots',1);
    title('P_12119085_KE_GB','interpreter','none','Fontsize',20);
    pause(1);
    out2  = LTStsLSmult(yin2,'msg',true,'plots',1);
    title('P_17049075_UA_LT','interpreter','none','Fontsize',20);
%}


%{
% Multiple level shift and variable selection. Example 1.
% Detection of multiple Level Shifts followed by variable selection on the
% dataset of example before
  load('P_17049075_TQV_full.txt'); % UA-LT
  dates2 = P_17049075_TQV_full(:,1);
  yin2   = P_17049075_TQV_full(:,2);
  out = LTStsLSmult(yin2,'maxLS',4,'alphaLTS',0.01,...
    'alphaLS',0.01,'thresLS',0.01,'plots',1,'msg',1);
  outX = out.outX;
model.trend = 2;
model.lshift = 0;
model.seasonal=303;
model.X = outX(:,3:end);

[out_model_1, out_reduced_1] = LTStsVarSel(yin2,'model',model,'plots',1);

[out_LTSts]=LTSts(yin2,'model',out_model_1,...
    'plots',1,'msg',0,'dispresults',1,'SmallSampleCor',1,'conflev',1-0.01/length(yin2));

%}

%{
% Multiple level shift and variable selection. Example 2.
% Detection of multiple Level Shifts followed by variable selection on the
% dataset of example before
  load('P_12119085_TQV_full.txt'); % UA-LT
  dates2 = P_12119085_TQV_full(:,1);
  yin2   = P_12119085_TQV_full(:,2);
  out = LTStsLSmult(yin2,'maxLS',4,'alphaLTS',0.01,...
    'alphaLS',0.01,'thresLS',0.01,'plots',1,'msg',1);
  outX = out.outX;
model.trend = 2;
model.lshift = 0;
model.seasonal=303;
model.X = outX(:,3:end);

[out_model_1, out_reduced_1] = LTStsVarSel(yin2,'model',model,'plots',1);

[out_LTSts]=LTSts(yin2,'model',out_model_1,...
    'plots',1,'msg',0,'dispresults',1,'SmallSampleCor',1,'conflev',1-0.01/length(yin2));

%}

%% Beginning of code

options=struct('msg',false,'plots',0,'maxLS',5,'alphaLTS',0.05,'bdpLTS',0.1,...
    'alphaLS',0.01,'thresLS',20,'startLS',3,'sampleLS',[]);
for i=1:2:length(varargin)
    options.(varargin{i})=varargin{i+1};
end

msg         = options.msg;
plots       = options.plots;
maxLS       = options.maxLS;
alphaLTS    = options.alphaLTS;
bdpLTS      = options.bdpLTS;
alphaLS     = options.alphaLS;
thresLS     = options.thresLS;
startLS     = options.startLS;
sampleLS    = options.sampleLS;

T = length(y);

model       = struct;
model.trend = 0;
LSpos       = nan(maxLS,1);
LSud        = nan(maxLS,1);
LSpval      = nan(maxLS,1);
LSsignif    = true;
ij=1;

model.X=[];
if T>thresLS
    if isempty(sampleLS) 
        model.lshift = startLS:T-startLS;
    else
        sampleLS(sampleLS<startLS)   = [];
        sampleLS(sampleLS>T-startLS) = [];
        model.lshift = sampleLS;
    end
end
while LSsignif==true &&  ij<=maxLS
    outTent=LTSts(y,'model',model,'conflev',1-alphaLTS/T,'plots',0,'msg',0,'bdp',bdpLTS);
    %in case there is no LS or only 1 save some results
    if T <= thresLS
        outTent.LevelShiftPval = NaN;
        outTent.posLS = NaN;
    end
    if ij==1
        out_LTSts_tmp=outTent; %#ok<NASGU>
    end
    % if the pvalue of the current LS is significant, look for another one
    if outTent.LevelShiftPval < alphaLS
        if msg
            disp(['significant LS at position ' num2str(outTent.posLS) ...
                ' (pval = ' num2str(outTent.LevelShiftPval) ')']);
        end
        LSpval(ij) = outTent.LevelShiftPval;
        LSpos(ij)  = outTent.posLS;
        if T > thresLS
            LSud(ij) = sign(outTent.yhat(outTent.posLS)-outTent.yhat(outTent.posLS-1));
        else
            LSud(ij) = NaN;
        end
        if T > thresLS
            model.X(outTent.posLS:T,end+1)=1;
        else
            model.X = [];
        end
        ij=ij+1;
        % Structure out_LTSts will contain the final model
        out_LTSts_tmp=outTent; %#ok<NASGU>
    else
        LSsignif=false;
        wedgeplot(outTent,'extradata',[y outTent.yhat]); 
    end
end

%% after having detected LS, call LTSts to identify outliers, etc

model.trend  = 0;
model.lshift = 0;
%model.seasonal = 2;           % two harmonics
outLTSts = LTSts(y,'model',model,'conflev',1-alphaLTS/T,'plots',plots,'msg',0,'bdp',bdpLTS,'yxsave',true);
yhat = outLTSts.yhat;
outliers = outLTSts.outliers;
goodobs = 1:T;  goodobs(outliers)=[];
res = (y-yhat);
resgood = res(goodobs);
% numerator of estimate of the error variance
numS2 = (resgood'*resgood);
% deviation from the mean
ytilde = y(goodobs)-mean(y(goodobs));
% total deviance
devtot = ytilde'*ytilde;
% R squared ;
R2=1-numS2/devtot;
outX=outLTSts.X;
%LSpos(isnan(LSpos)) =[];
%LSpval(isnan(LSpos))=[];

out.R2 = R2;
out.yhat = yhat; 
out.LSpos = LSpos;
out.LSpval = LSpval;
out.LSud = LSud; 
out.outliers = outliers;
out.outX = outX;
end
%FScategory:REG-Regression
