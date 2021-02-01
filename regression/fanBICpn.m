function [out]=fanBICpn(outFSRfan, varargin)
%fanBICpn uses the output of FSRfan called with input option family 'YJpn' to choose la_P and la_N
%
% This function finds the best values of transformation parameter for positive and
% negative observations in linear regression using BIC and agreement
% index (AGI).
%
%<a href="matlab: docsearchFS('fanBICpn')">Link to the help function</a>
%
% Required input arguments:
%
%  outFSRfan :  Structure created with function FSRfan. Structure.
%               Structure containing the following fields
%outFSRfan.Score  =  (n-init) x length(la)+1 matrix:
%               1st col = fwd search index;
%               2nd col = value of the score test in each step
%               of the fwd search for la;
%  outFSRfan.Scorep = (n-init) x 2 matrix containing the values of the
%               score test for positive observations for each value of the
%               transformation parameter.
%               1st col = fwd search index;
%               2nd col = value of the (positive) score test in each step
%               of the fwd search for la;
% outFSRfan.Scoren  = (n-init) x 2 matrix containing the values of the
%               score test for positive observations for each value of the
%               transformation parameter:
%               1st col = fwd search index;
%               2nd col = value of the (negative) score test in each step
%               of the fwd search for la;
% outFSRfan.Scoreb  = (n-init) x 2+1 matrix containing the values of the
%               score test for the joint presence of both constructed
%               variables (associated with positive and negative
%               observations) for each value of the transformation
%               parameter.  In this case the reference distribution is the
%               $F$ with 2 and subset_size-p degrees of freedom.
%               1st col = fwd search index (subset_size);
%               2nd col = value of the score test in each step
%               of the fwd search for la;
% outFSRfan.Scoremle  = (n-init) x 2 matrix containing the values of the
%               (score) likelihood ratio test for the joint presence of
%               both constructed variables (associated with positive and
%               negative observations) for each value of the transformation
%               parameter.  In this case the reference distribution is the
%               $F$ with 2 and subset_size-p degrees of freedom.
%               1st col = fwd search index (subset_size);
%               2nd col = value of the score test in each step of the fwd search for la;
%               outFSRfan.Scoremle is present only if FSRfan has been
%               called with input option scoremle set to true.
%       outFSRfan.la   =  scalar containing the value of lambda for which FSRfan
%               was computed.
%       outFSRfan.bs   =  matrix of size p x 1 containing the units forming
%               the initial subset for lambda.
%         outFSRfan.y  = a vector containing the response
%         outFSRfan.X  = a matrix containing the explanatory variables.
%                        Note that this matrix includes the column of ones.
%      Data Types - struct
%
% Optional input arguments:
%
%   laRangeAndStep :  values of laP and laN to explore. Vector of length 3.
%                  The grid search for laP and laN is done in the grid
%                  outFSRfan.la(1)-laRangeAndStep(1):laRangeAndStep(2):outFSRfan.la(1)+laRangeAndStep(3).
%                  For example, If laRangeAndStep(1)=0.8,  laRangeAndStep(2) is 0.1 and
%                  laRangeAndStep(3)=0.6 and input parameter outFSRfan.la(1)=0.5, the
%                  grid search for laP and laN is done in the grid
%                  (0.5-0.8):0.1:(0.5+0.6). The default value for
%                  laRangeAndStep is [0.75 0.25 0.75].
%                   Example - 'laRangeAndStep',[1 0.2 1.2]
%                   Data Types - double
%
%
%       conflev :   Confidence level. Scalar. Confidence level
%                   to evaluate the exceedances in hte fanplot.
%                   Default confidence level is 0.9999 that is signals are
%                   considered when there is an exceedance for confidence
%                   level for at least 3 consecutive times.
%                   Example - 'conflev',[0.999]
%                   Data Types - double
%    fraciniFSR :   fraction of observations to initialize search for outlier detection.
%                   Scalar. After exceedance procedure based on the score
%                   test a subset of obverations in agreement with a
%                   transformation is found. On this subset we perform
%                   outlier detection using FSR. fraciniFSR specifies the
%                   fraction of observations to start monitoring
%                   exceedances of the minimum deletion residuals. The
%                   default value of fraciniFSR is 0.8.
%                   Example - 'fraciniFSR',0.85
%                   Data Types - double
%
%       init    : Step to start monitoring exceedances. Scalar. It specifies the initial
%                 subset size to start monitoring exceedances of the
%                 fanplot. If init is not specified it set equal
%                 to round(n*0.6):
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
%
%
%      bonflev  : Signal to use to identify outliers. Scalar. Option to be
%                used if the distribution of the data is
%                 strongly non normal and, thus, the general signal
%                 detection rule based on consecutive exceedances cannot be
%                 used. In this case bonflev can be:
%                 - a scalar smaller than 1 which specifies the confidence
%                   level for a signal and a stopping rule based on the
%                   comparison of the minimum MD with a
%                   Bonferroni bound. For example if bonflev=0.99 the
%                   procedure stops when the trajectory exceeds for the
%                   first time the 99% bonferroni bound.
%                 - A scalar value greater than 1. In this case the
%                   procedure stops when the residual trajectory exceeds
%                   for the first time this value.
%                 Default value is '', which means to rely on general rules
%                 based on consecutive exceedances.
%               Example - 'bonflev',0.99
%               Data Types - double
%
%      scoremle: likelihood ratio test for the two different transformation
%                parameters $\lambda_P$ and $\lambda_N$. Boolean. If
%                scoremle is true and field Scoremle is present in input
%                structure outFSRfan is present we check exceedance of the
%               threshold according to according to likelihood
%               ratio test else we check exceedance of the threshold
%               according to outFSRfan.Scoreb.
%               Example - 'scoremle',true
%               Data Types - logical
%
%       nsamp   :   Number of subsamples to extract. Scalar.
%                   Number of subsamples which will be extracted in order to find
%                   initial subset for each candidate value of lambda.
%                   If nsamp=0 all subsets will be
%                   extracted. They will be (n choose p). Remark: if the
%                   number of all possible subset is <1000 the default is
%                   to extract all subsets otherwise just 1000.
%                   Example - 'nsamp',1000
%                   Data Types - double
%
%       msg    :  Level of output to display. Scalar. It controls whether
%                 to display or not messages on the screen
%                 If msg==1 (default) messages are displayed on the screen about
%                   values or la_P and la_N which are being analyzed.
%                 else no message is displayed on the screen.
%               Example - 'msg',1
%               Data Types - double
%
% plots    :    Plot on the screen. Scalar structure.
%
%               Case 1: plots option used as scalar.
%               - If plots=0, plots are not generated.
%               - If plots=1 (default), 4 heatmaps are shown on the screen.
%                 The first plot ("BIC") shows the values of BIC, the
%                 second ("AGI") shows the values of the agreement index,
%                 the third ('Obs') the number of observations in agreement
%                 with the transformation excluding the outliers and the
%                 fourth ('R2c') the final value of R2 (corrected for truncation).
%
%               Case 2: plots option used as struct.
%                 If plots is a structure it may contain the following fields:
%                 plots.name = cell array of strings which enables to
%                   specify which plot to display.
%                   plots.name={'Obs'; 'BIC'; 'AGI'; 'R2c'};
%                   is exactly equivalent to plots=1
%                   For the explanation of the above plots see plots=1.
%                   If plots.name=={ 'Obs'; 'BIC'; 'AGI'; 'R2c';...
%                                   'ObsWithOut'; 'AGIW'; 'R2'}; or
%                   plots.name={'all'};
%                   it is also possible to view the heatmap referred to the
%                   number of obserations in agreement with transformation
%                   before outlier detection ('ObsWithOut') the weighted
%                   version of the agreement index ('AGIW') and the
%                   orginal value of R2 before correction for truncation.
%                   Example - 'plots', 1
%                   Data Types - single | double | struct
%
%
%
%  Output:
%
%
%         out:   structure which contains the following fields
%
% out.Summary = k-by-9 table where k is the number of values of laPosxlaNeg which have been considered.
%              out.Summary contains the following information:
%              1st column= value of laPos (transformation for positive
%              values of y);
%              2nd column= value of laNeg (transformation for negative
%              values of y);
%              3rd col = number of observations in agreement with the
%              transformation before outlier detection;
%              4th col = number of observations in agreement with the
%              transformation after outlier detection;
%              5th col = value of BIC;
%              6th col = value of the agreement index;
%              7th col = value of the agreement index weighted;
%              8th col = value of R2 based on observations in agreement
%                   with transformation after outlier detection.
%              9th col = value of R2 corrected for elliptical trunction.
%  out.labestBIC = vector of length 2 containing best values of laP and laN
%              according to BIC.
%  out.labestAGI = vector of length 2 containing best values of laP and laN
%              according to agreement index.
%  out.ty =    transformed response accordint to out.labestBIC.
%     out.rsq  = the multiple R-squared value for the transformed values
%      out.y  = n x 1 vector containing the original y values.
%      out.X  = n x p matrix containing the original X matrix.
%
% See also: FSRfan, fanBIC, fanplot
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
% Atkinson, A.C. and Riani, M. (2002a), Tests in the fan plot for robust,
% diagnostic transformations in regression, "Chemometrics and Intelligent
% Laboratory Systems", Vol. 60, pp. 87-100.
% Atkinson, A.C. Riani, M., Corbellini A. (2019), The analysis of
% transformations for profit-and-loss data, Journal of the Royal
% Statistical Society, Series C, "Applied Statistics",
% https://doi.org/10.1111/rssc.12389
% Atkinson, A.C. Riani, M. and Corbellini A. (2020), The Box-Cox
% Transformation: Review and Extensions, "Statistical Science", in press.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('fanBICpn')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of use of FSRfanBICpn with all default options.
    % Load the Investment funds data.
    YY=load('fondi_large.txt');
    y=YY(:,2);
    X=YY(:,[1 3]);
    yXplot(y,X);
    n=length(y);
    [outFSRfan]=FSRfan(y,X,'plots',0,'init',round(n*0.3),'nsamp',10000,'la',[0 0.25 0.5 0.75 1 1.25],'msg',0,'family','YJ');
    [outini]=fanBIC(outFSRfan,'plots',0);
    % labest is the best value imposing the constraint that positive and
    % negative observations must have the same tramsformation parameter.
    labest=outini.labest;
    % Compute test for positive and test for negative using labest
    [outFSRfanpn]=FSRfan(y,X,'msg',0,'family','YJpn','la',labest,'plots',0);
    % Check if two different transformations are needed for positive and negative values
    out=fanBICpn(outFSRfanpn);
%}

%{
    %% Example of the use of option laRangeAndStep.
    % Use simulated data from Atkinson Riani and Corbellini (2020)
    rng('default')
    rng(10)
    n=1000;
    p=3;
    kk=200;
    X=randn(n,p);
    beta=[ 1; 1; 1]*0.3;
    sig=0.5;
    eta=X*beta;

    init=6;
    lapos=1.5;
    laneg=0;

    y=eta+sig*randn(n,1);
    % Data contamination
    y(1:kk)=y(1:kk)-1.9;
    ypos=y>0;
    ytra=y;
    ytra(ypos)=normYJ(y(ypos),[],lapos,'inverse',true,'Jacobian',false);
    ytra(~ypos)=normYJ(y(~ypos),[],laneg,'inverse',true,'Jacobian',false);
    y=ytra;

    % Initial fan plot
    outFSRfan=FSRfan(y,X,'la',[0.5 0.75 1 1.25 1.5],'family','YJ','plots',0,'init',init,'msg',0);

    % Find best value of lambda according to BIC
    % (same value of lambda for positive and negative observations).
    [outUniqueLambda]=fanBIC(outFSRfan,'plots',0);
    BIC=outUniqueLambda.BIC;
    labest=outUniqueLambda.labest;

    % Check if two different transformations are needed for positive and negative values
    [outFSRfanpn]=FSRfan(y,X,'msg',0,'family','YJpn','la',labest);
    % option laRangeAndStep
    laRangeAndStep=[1.5 0.25 0.5];
    out=fanBICpn(outFSRfanpn,'laRangeAndStep',laRangeAndStep);
%}

%{
    %% Example of the use of options fraciniFSR and plots.
    % Balance sheets data.
    XX=load('BalanceSheets.txt');
    % Define X and y
    y=XX(:,6);
    X=XX(:,1:5);
    n=length(y);
    la=[0 0.25 0.5 0.75 1 1.25];
    [outFSRfan]=FSRfan(y,X,'plots',1,'init',round(n*0.3),'nsamp',5000,'la',la,'msg',0,'family','YJ');
    [outini]=fanBIC(outFSRfan,'plots',0);
    % labest is the best value imposing the constraint that positive and
    % negative observations must have the same tramsformation parameter.
    labest=outini.labest;
    % Compute test for positive and test for negative using labest
    indexlabest=find(labest==la);
    % Find initial subset to initialize the search.
    lms=outFSRfan.bs(:,indexlabest);
    [outFSRfanpn]=FSRfan(y,X,'msg',0,'family','YJpn','la',labest,'plots',0,'lms',lms);
    % Check if two different transformations are needed for positive and negative values
    % Start monitoring the exceedances in the subset in agreement with a
    % transformation from 90 per cent.
    fraciniFSR=0.90;
    % option plots (just show the BIC and the agreement index plot).
    plots=struct;
    plots.name={'BIC','AGI'};
    nsamp=2000;
    out=fanBICpn(outFSRfanpn,'fraciniFSR',fraciniFSR,'plots',plots,'nsamp',nsamp);
%}



%% Beginning of code

% fopen() not being supported by localthreads prevent us
% from using this code
% 9.8 is MATLAB 2020a where  parpool('threads') was first introduced
% numbertotest = 9.8;
% MLver=verLessThanFS(numbertotest);
%
%
% pp = gcp('nocreate');
% if MLver == true &&  isempty(pp)
%     parpool('local');
% elseif MLver == false && isempty(pp)
%     parpool('threads');
% end

Xwithintercept=outFSRfan.X;
[n,pwithintercept]=size(Xwithintercept);

y=outFSRfan.y;
logn=log(n);

init=floor(n*0.6);
conflev=0.999;
bonflev='';
msg=1;
plots=1;
laRangeAndStep=[0.75 0.25 0.75];
labest=outFSRfan.la(1);
nsamp=2000;
% fraciniFSR = fraction of units to initialize routine FSR for outlier
% detection on the subset of units in agreement with a particular
% transformation
fraciniFSR=0.80;

if length(outFSRfan.la)>1
    warning('FSDA:fanBICpn:labestNotScalar',['Input structure outFSRfan contains the values ' ...
        'of the score test for more than one value of lambda. Here we assume that labest is the first ' ...
        'that is we assume that labest is ' num2str(labest)])
end

scoremle=false;
UserOptions=varargin(1:2:length(varargin));


if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:fanBIC:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    options=struct('plots',plots,'init',init,'conflev',conflev,...
        'bonflev',bonflev,'msg',msg,'laRangeAndStep',laRangeAndStep, ...
        'scoremle',scoremle, 'nsamp', nsamp,'fraciniFSR',fraciniFSR);
    
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    
    plots=options.plots;
    scoremle=options.scoremle;
    conflev=options.conflev;
    laRangeAndStep=options.laRangeAndStep;
    msg=options.msg;
    bonflev=options.bonflev;
    nsamp=options.nsamp;
    fraciniFSR=options.fraciniFSR;
end


% Extract the indexes of the subsets once and for all
nsamp = subsets(nsamp,n,pwithintercept);


if length(laRangeAndStep)~=3
    error('FSDA:fanBICpn:WrongInputOpt', 'Input option laRangeAndStep must be a vector of length 3.');
end

nonnegs=y>=0;
negs=~nonnegs;
ynonnegs=y(nonnegs);
ynegs=y(negs);
logynonnegsp1=log(ynonnegs+1);
log1mynegs=log(1-ynegs);
SumLogYp=sum(logynonnegsp1);
SumLogYn=sum(-log1mynegs);

hhseq=(init:n)';
vt = norminv(0.5*(1+hhseq/n));
corrfactor = (1-2*(n./hhseq).*vt.*normpdf(vt));
corrfactor(end)=1;

% If this condition is true it means that the trajectory for positive or
% negative have gone beyond the bands along the search


laposcand=labest-laRangeAndStep(1):laRangeAndStep(2):labest+laRangeAndStep(3);
lanegcand=laposcand;


booInit=outFSRfan.Scorep(:,1)>=init;
Scorep=outFSRfan.Scorep(booInit,:);
Scoren=outFSRfan.Scoren(booInit,:);


% Scorep=outFSRfan.Scorep;
% outFSRfan.Scoren
if sum(Scorep(:,2)>Scoren(:,2))> size(Scorep,1)/2
    % in this case lapos must be greater than labest and laneg must be smaller than labest
    laposcand=laposcand(laposcand>=labest);
    lanegcand=lanegcand(lanegcand<=labest);
else
    laposcand=laposcand(laposcand<=labest);
    lanegcand=lanegcand(lanegcand>=labest);
end

ytra=y;
Exc=zeros(length(laposcand)*length(lanegcand),9);
BBlacell=cell(length(laposcand)*length(lanegcand),1);

Excnegj=zeros(length(lanegcand),9);
BBlacellj=cell(length(lanegcand),1);
ynegs=y(negs);
ynonnegs=y(~negs);
ij=1;
llanegcand=length(lanegcand);
for jjpos=1:length(laposcand)
    laposj=laposcand(jjpos);
    ytra(~negs)=normYJ(ynonnegs,[],laposj,'inverse',false,'Jacobian',false);
    
    parfor jjneg=1:llanegcand
        corrfactor1=corrfactor;
        lanegj=lanegcand(jjneg);
        ytra1=ytra;
        ytra1(negs)=normYJ(ynegs,[],lanegj,'inverse',false,'Jacobian',false);
        
        [outtest]=FSRfan(ytra1,Xwithintercept,'plots',0,'init',init,'nsamp',nsamp,'la',1,...
            'msg',0,'family','YJall','conflev',conflev,'scoremle',scoremle,'nocheck',1);
        
        if msg==1
            % title([num2str(laposj) ' -- ' num2str(lanegj)])
            disp(['Analyzing la_P=' num2str(laposj)  ' and la_N='  num2str(lanegj)])
        end
        
        if scoremle==true
            outtest.Scoreb=outtest.Scoremle;
        end
        Sco=[outtest.Score outtest.Scorep(:,2) outtest.Scoren(:,2) outtest.Scoreb(:,2)];
        % Exceedance on the global Score test
        [mmTest,BBla] = fanBICpncore(ytra1,Xwithintercept,Sco,bonflev,outtest.bs,init,conflev,fraciniFSR);
        
        % Store BBla
        BBlacellj{jjneg}=BBla;
        
        % vt = variance of the truncated normal distribution
        % 1-hh/n is the trimming percentage
        hh=sum(BBla==2);
        
        vt = norminv(0.5*(1+hh/n));
        if hh<n
            factor = 1/(1-2*(n/hh)*vt.*normpdf(vt));
        else
            factor=1;
        end
        
        
        % Store value of R2 in the step before the entry of the
        % first outlier
        yb=ytra1(BBla==2);
        Xb=Xwithintercept(BBla==2,:);
        b=Xb\yb;
        e=yb-Xb*b;
        sumSqres=(e'*e);
        sigma2hat=sumSqres/hh;
        ym=yb-sum(yb)/hh;
        R2=1-sumSqres/(ym'*ym);
        
        % logJ = log of the Jacobian for all the observations
        logJ=(laposj-1)*SumLogYp +(lanegj-1)*SumLogYn;
        
        %  BIC (the larger the better)
        BIC=2*(-0.5*n*log(factor*sigma2hat)+logJ)-(pwithintercept+2+n-hh)*logn;
        
        Scopn=[outtest.Scorep(:,1:2) outtest.Scoren(:,2) outtest.Score(:,2) outtest.Scoreb(:,2)];
        
        R2c=R2/factor;
        
        boo=Scopn(:,1)<=hh;
        ScopnSel=Scopn(boo,:);
        corrfactorSel=corrfactor1(boo);
        tPNdiff=abs(ScopnSel(:,2)-ScopnSel(:,3));
        
        tPNdiffMean=mean(tPNdiff);
        sumwei=sum(corrfactorSel);
        tPNdiffMeanW=sum(tPNdiff.*corrfactorSel)/sumwei;
        
        
        % disp(['pnMean' num2str(pnMean)])
        tdiff=abs(ScopnSel(:,4));
        t0diffMean=mean(tdiff);
        t0diffMeanW=sum(tdiff.*corrfactorSel)/sumwei;
        
        
        AGI=1./(tPNdiffMean*t0diffMean*(factor.^2));
        AGIW=1./(tPNdiffMeanW*t0diffMeanW*(factor.^2));
        
        
        Excnegj(jjneg,:)=[laposj lanegj mmTest(1:2)' BIC AGI AGIW R2c  R2];
        % ij=ij+1;
    end
    Exc(ij:ij+llanegcand-1,:)=Excnegj;
    BBlacell(ij:ij+llanegcand-1)=BBlacellj;
    ij=ij+llanegcand;
end

% Set to NaN the combination of values of laP and laN which are not
% admissible
booAgiNaN=isnan(Exc(:,6));
Exc(booAgiNaN,3:end)=NaN;

varnames={'laP';'laN';'ObsWithOut';'Obs';'BIC'; 'AGI';'AGIW';'R2c'; 'R2'};
Exctable=array2table(Exc,'VariableNames',varnames);

% poscell= intersect(find(Exctable.laP==laposj),find(Exctable.laN==lanegj));
% BBla=BBlacell{poscell};
% fitlm(X(BBla==2,:),ytra(BBla==2))
% yXplot(ytra,X,'group',BBla)
[~, indmaxBIC]=max(Exctable{:,'BIC'});
% disp('Best values of transformation parameters')
labestBIC=Exctable{indmaxBIC,1:2};
[~, indmaxAGI]=max(Exctable{:,'AGI'});
labestAGI=Exctable{indmaxAGI,1:2};

% plotdef = list of the plots which are produced by default (is plots=1)
plotdef={'Obs','BIC', 'AGI','R2c'};
% plotall = list of all available plots
plotall=varnames(3:end);

if isstruct(plots)
    fplots=fieldnames(plots);
    
    d=find(strcmp('name',fplots));
    if d>0
        name=plots.name;
        if ~iscell(name)
            error('FSDA:fanBICpn:Wronginput','plots.name must be a cell')
        end
        if strcmp(name,{'all'})
            name=plotall;
        else
            % Check that the specified names is in the list of available names.
            chkoptions(cell2struct(plotall,plotall),name)
        end
    else
        name=plotdef;
    end
    
elseif plots==1
    name=plotdef;
else
    name='';
end


namej = 'Obs';
d=strcmp(namej,name);
if max(d)>0
    figure('Name',namej,'Visible','on');
    hold('off')
    if ~verLessThanFS(9.2) % >2016b
        h=heatmap(Exctable,'laP','laN','ColorVariable',namej,'MissingDataColor','w');
        title('Numb obs. in agreement with the different transformations')
        h.XLabel= '\lambda_P';
        h.YLabel= '\lambda_N';
    else
        text(0.1,0.5,'Heatmap cannot be shown in this version of MATLAB','Units','normalized')
        text(0.1,0.3,'You need at least MATLAB 2017a','Units','normalized')
    end
end


namej = 'R2';
d=strcmp(namej,name);
if max(d)>0
    figure('Name',namej,'Visible','on');
    if ~verLessThanFS(9.2) % >2016b
        h=heatmap(Exctable,'laP','laN','ColorVariable',namej,'MissingDataColor','w');
        title('R2 not corrected')
        h.XLabel= '\lambda_P';
        h.YLabel= '\lambda_N';
    else
        text(0.1,0.5,'Heatmap cannot be shown in this version of MATLAB','Units','normalized')
        text(0.1,0.3,'You need at least MATLAB 2017a','Units','normalized')
    end
end

namej = 'R2c';
d=strcmp(namej,name);
if max(d)>0
    figure('Name',namej,'Visible','on');
    if ~verLessThanFS(9.2) % >2016b
        h=heatmap(Exctable,'laP','laN','ColorVariable',namej,'MissingDataColor','w');
        title('R2 corrected for truncation')
        h.XLabel= '\lambda_P';
        h.YLabel= '\lambda_N';
    else
        text(0.1,0.5,'Heatmap cannot be shown in this version of MATLAB','Units','normalized')
        text(0.1,0.3,'You need at least MATLAB 2017a','Units','normalized')
    end
    
end

namej = 'AGI';
d=strcmp(namej,name);
if max(d)>0
    figure('Name',namej,'Visible','on');
    if ~verLessThanFS(9.2) % >2016b
        h=heatmap(Exctable,'laP','laN','ColorVariable',namej,'MissingDataColor','w');
        % title('1/(|ScoP-ScoN|_c*|Sco|_c)')
        title(['Agreement index: best \lambda_P=' num2str(labestAGI(1)) ' best  \lambda_N=' num2str(labestAGI(2))])
        h.XLabel= '\lambda_P';
        h.YLabel= '\lambda_N';
    else
        text(0.1,0.5,'Heatmap cannot be shown in this version of MATLAB','Units','normalized')
        text(0.1,0.3,'You need at least MATLAB 2017a','Units','normalized')
    end
end

namej = 'AGIW';
d=strcmp(namej,name);
if max(d)>0
    [~, indmax]=max(Exctable{:,'AGIW'});
    labestAGIw=Exctable{indmax,1:2};
    
    figure('Name',namej,'Visible','on');
    if ~verLessThanFS(9.2) % >2016b
        h=heatmap(Exctable,'laP','laN','ColorVariable',namej,'MissingDataColor','w');
        title(['Agreement index weighted: best \lambda_P=' num2str(labestAGIw(1)) ' best  \lambda_N=' num2str(labestAGIw(2))])
        h.XLabel= '\lambda_P';
        h.YLabel= '\lambda_N';
    else
        text(0.1,0.5,'Heatmap cannot be shown in this version of MATLAB','Units','normalized')
        text(0.1,0.3,'You need at least MATLAB 2017a','Units','normalized')
    end
end

namej = 'BIC';
d=strcmp(namej,name);
if max(d)>0
    figure('Name',namej,'Visible','on');
    if ~verLessThanFS(9.2) % >2016b
        h=heatmap(Exctable,'laP','laN','ColorVariable',namej,'MissingDataColor','w');
        title(['BIC: best \lambda_P=' num2str(labestBIC(1)) ' best  \lambda_N=' num2str(labestBIC(2))])
        h.XLabel= '\lambda_P';
        h.YLabel= '\lambda_N';
    else
        text(0.1,0.5,'Heatmap cannot be shown in this version of MATLAB','Units','normalized')
        text(0.1,0.3,'You need at least MATLAB 2017a','Units','normalized')
    end
    
end

out=struct;
out.Summary=Exctable;
% out.mmstop=mmstop;
% out.BBla=BBla;
out.labestBIC=labestBIC;
out.labestAGI=labestAGI;
out.ty=normYJpn(y,1,labestBIC);
out.rsq=Exctable{indmaxBIC,'R2'};
out.y=y;
out.X=outFSRfan.X;

%      out.ty  = n x 1 vector containing the transformed y values.
%     out.rsq  = the multiple R-squared value for the transformed values
%      out.y  = n x 1 vector containing the original y values.
%      out.X  = n x p matrix containing the original X matrix.

end


function [mmstop, BBla] = fanBICpncore(y,Xwithintercept,Sco,bonflev,bs,init,conflev,fraciniFSR)
% Check exceedance of Tcentral Tpos Tneg or Tboth
% Sco has 5 columns fwd index, global test, pos, neg and both
n=length(y);
seq=1:n;
% Divide central part from final part of the search
istep = n-floor(13*sqrt(n/200));

nmdr=size(Sco,1);
upperbound=conflev+(1-conflev)/2;

th=norminv(upperbound);


mmstop=n*ones(2,1);
BBla=zeros(n,1);

Scolaj=abs(Sco(:,2));
ScolajP=abs(Sco(:,3));
ScolajN=abs(Sco(:,4));
ScolajB=abs(Sco(:,5));

%quant = 0.999
quant=conflev;
thB=finv(quant,2,Sco(:,1)-size(Xwithintercept,2)-2);

signal=false;
for i=3:nmdr
    
    if i<istep-init+1 % CENTRAL PART OF THE SEARCH
        % order interchange between Tp and Tn for a quintuplet and max
        % interchange greater than 1
        PosMinNeg=sum(Sco(i-2:i+2,3)<Sco(i-2:i+2,4))==5 && max(Sco(i-2:i+2,4)-Sco(i-2:i+2,3))>1;
        % FandT=(ScolajB(i)>thB(i) && ScolajB(i-1)>thB(i) && ScolajB(i-2)>thB(i) && ScolajB(i+1)>thB(i) && ScolajB(i+2)>thB(i));
        
        % F test is greater than threhold for 5 consecutive times (quintuplet) and both tP and tN
        % are greater then threshold
        FandT=sum(ScolajB(i-2:i+2)>thB(i))==5 &&  (ScolajP(i)>th ||  ScolajN(i)>th);
        
        % Extreme triplet out of the envelope for general test or
        % quintuplet of Tp or Tn out of the threshold or
        if (Scolaj(i)>th && Scolaj(i-1)>th && Scolaj(i+1)>th) || ...
                (ScolajP(i)>th && ScolajP(i-1)>th && ScolajP(i-2)>th && ScolajP(i+1)>th && ScolajP(i+2)>th) || ...
                (ScolajN(i)>th && ScolajN(i-1)>th && ScolajN(i-2)>th && ScolajN(i+1)>th && ScolajN(i+2)>th) || ...
                FandT == true || ...
                PosMinNeg == true
            mdagger=Sco(i-2,1);
            signal=true;
            break
        end
    else
        if  Scolaj(i)>th || ScolajP(i)>th || ScolajN (i)>th
            % if  Scolaj(i)>th || ScolajB(i)>thB(i)
            mdagger=Sco(i-1,1);
            signal=true;
            break
        end
    end
    % Now get out of the loop and do outlier detection
end

if signal==true
    mmstop(1)=mdagger;
    % Find units belonging to subset for laj in step mdagger
    [~,BB] = FSRbsb(y,Xwithintercept,bs,'init',init,'nocheck',1,'msg',0,'bsbsteps',unique([init mdagger]));
    % good = good units for score test
    good=seq(~isnan(BB(:,end)));
else
    good=seq;
end
BBla(good,1)=1;

n1=length(good);
initlarger=round(n1*fraciniFSR);

% lms contains the new indexs of bs inside good. lms will be used to
% initialize the search without calling LXS
lms=zeros(length(bs),1);
for i=1:length(bs)
    newindexbsi=find(good==bs(i),1);
    if isempty(newindexbsi)
        lms=1;
        break
    else
        lms(i)=newindexbsi;
    end
end

%  Outlier detection
outl=FSR(y(good),Xwithintercept(good,:),'bonflev',bonflev,...
    'lms',lms,'nocheck',1,'plots',0,'msg',0,'init',initlarger);

if isnan(outl.ListOut)
    
    LowerEnv=FSRenvmdr(n1,size(Xwithintercept,2),'init',initlarger);
    
    % Find last time in which trajectory went below the 1% lower envelope
    % and remained below it
    % outl.mdr(find(outl.mdr(:,2)>LowerEnv(:,2),1,'last'),1)
    if outl.mdr(end,2)<LowerEnv(end,2)
        for jj=1:size(outl.mdr,1)-1
            if outl.mdr(end-jj,2)>LowerEnv(end-jj,2)
                break
            else
            end
        end
        
        Sizenewgood=outl.mdr(end-jj,1);
        
        bsStepInit=seq(~isnan(BB(:,1)));
        if Sizenewgood<init
            [~,BBn] = FSRbsb(y,Xwithintercept,bs,'init',Sizenewgood,'nocheck',1,'msg',0,...
                'bsbsteps',Sizenewgood);
        else
            [~,BBn] = FSRbsb(y,Xwithintercept,bsStepInit,'init',length(bsStepInit),'nocheck',1,'msg',0,...
                'bsbsteps',Sizenewgood);
        end
        
        goodf=seq(~isnan(BBn(:,1)));
        BBla(goodf)=2;
    else
        goodf=1:length(good);
        BBla(good(goodf))=2;
    end
else
    
    goodf=setdiff(1:length(good),outl.ListOut);
    BBla(good(goodf))=2;
end
mmstop(2,1)=length(goodf);
end

%FScategory:VIS-Reg