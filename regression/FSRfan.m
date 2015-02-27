function [out]=FSRfan(y,X,varargin)
%FSRfan monitors the values of the score test statistic for each lambda
%
%<a href="matlab: docsearchFS('FSRfan')">Link to the help function</a>
%
%  Required input arguments:
%
%    y: A vector with n elements that contains the response variable. y can
%       be both a row of column vector.
%    X: Data matrix of explanatory variables (also called 'regressors') of
%       dimension (n x p-1). Rows of X represent observations, and columns
%       represent variables. Missing values (NaN's) and infinite values
%       (Inf's) are allowed, since observations (rows) with missing or
%       infinite values will automatically be excluded from the
%       computations. NOTICE THAT THE INTERCEPT MUST ALWAYS BE INCLUDED.
%
%  Optional input arguments:
%
%   intercept   :   If 1, a model with constant term will be fitted
%                   (default), if 0, no constant term will be included.
%       nocheck :   Scalar. If nocheck is equal to 1 no check is performed
%                   on matrix y and matrix X. Notice that y and X are left
%                   unchanged. In other words the additional column of ones
%                   for the intercept is not added. As default nocheck=0.
%           la  :   vector which specifies for which values of the
%                   transformation parameter it is necessary to compute the
%                   score test.
%                   Default value of lambda is la=[-1 -0.5 0 0.5 1]; that
%                   is the five most common values of lambda
%           h   :   The number of observations that have determined the
%                   least trimmed (median of) squares estimator.
%                   Generally h is an integer greater or equal than
%                   [(n+size(X,2)+1)/2] but smaller then n
%       nsamp   :   Number of subsamples which will be extracted to find
%                   the robust estimator. If nsamp=0 all subsets will be
%                   extracted. They will be (n choose p). Remark: if the
%                   number of all possible subset is <1000 the default is
%                   to extract all subsets otherwise just 1000.
%       lms     :   Scalar. If lms=1 (default) Least Median of Squares is
%                   computed, else Least trimmed of Squares is computed.
%       init    :   scalar which specifies the initial subset size to start
%                   monitoring the value of the score test, if init is not
%                   specified it will be set equal to:
%                    p+1, if the sample size is smaller than 40;
%                    min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%       plots   :   Scalar. If plots=1 the fan plot is produced
%                   else (default) no plot is produced
%
%                   REMARK: all the following options work only if plots=1
%
%       conflev :   confidence level for the bands (default is 0.99 that is
%                   we plot two horizontal lines in correspondence of value
%                   -2.58 and 2.58
%       titl    :   a label for the title (default: 'Fan plot')
%       labx    :   a label for the x-axis (default: 'Subset size m')
%       laby    :   a label for the y-axis (default:'Score test statistic')
%       xlimx   :   vector with two elements controlling minimum and maximum
%                   of the x axis. Default value is [init n]
%       ylimy   :   vector with two elements controlling minimum and
%                   maximum of the y axis. Default value for
%                   ylimy(1)=max(min(score_test),-20). Default value for
%                   ylimy(2)=min(max(score_test),20).
%       lwd     :   Scalar which controls linewidth of the curves which
%                   contain the score test. Default line width=2.
%       lwdenv  :   Scalar which controls the width of the lines associated
%                   with the envelopes. Default is lwdenv=1.
%       FontSize:   Scalar which controls the font size of the labels of
%                   the axes. Default value is 12.
%    SizeAxesNum:   Scalar which controls the size of the numbers of the
%                   axes. Default value is 10.
%         msg   :   scalar which controls whether to display or not
%                   messages on the screen If msg==1 (default) messages are
%                   displayed on the screen about estimated time to compute
%                   the LMS (LTS) for each value of lamabda else no message
%                   is displayed on the screen
%       tag     :   string which identifies the handle of the plot which
%                   is about to be created. The default is to use tag
%                   'pl_fan'. Notice that if the program finds a plot which
%                   has a tag equal to the one specified by the user, then
%                   the output of the new plot overwrites the existing one
%                   in the same window else a new window is created
%
%  Output:
%
%    The output consists of a structure 'out' containing the following fields:
%  out.Score  = (n-init) x length(la)+1 matrix containing the values of the
%               score test for each value of the transformation parameter
%               1st col = fwd search index
%               2nd col = value of the score test in each step of the
%               fwd search for la(1)
%               ...........
%               end col = value of the score test in each step of the fwd
%               search for la(end)
%  out.la     = vector containing the values of lambda for which fan plot
%               is constructed
%  out.bs     = matrix of size p x length(la) containing the units forming
%               the initial subset for each value of lambda
%  out.Un     = cell of size length(la).
%               out.Un{i} is a n-init) x 11 matrix which contains the unit(s) included in
%               the subset at each step in the search associated with la(i)
%               REMARK: in every step the new subset is compared with the old subset. Un
%               contains the unit(s) present in the new subset but not in
%               the old one Un(1,:) for example contains the unit included
%               in step init+1 ... Un(end,2) contains the units included in the
%               final step of the search
%  out.y      = A vector with n elements that contains the response
%               variable which has been used
%  out.X      = Data matrix of explanatory variables
%               which has been used (it also contains the column of ones if
%               input option intercept was missing or equal to 1)
%
%
% See also
%
% References:
%
%   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
%   Springer Verlag, New York.
%   Atkinson, A.C. and Riani, M. (2002a). Tests in the fan plot for robust,
%   diagnostic transformations in regression, Chemometrics and Intelligent
%   Laboratory Systems, Vol. 60, pp. 87–100.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRfan')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % Store values of the score test statistic
    % for the five most common values of lambda
    [out]=FSRfan(y,X);
    % Produce a fan plot and display it on the screen
    fanplot(out);
%}

%{
    % Produce a personalized fan plot with required font sizes
    % for labels and axes
    [out]=FSRfan(y,X,'plots',1,'FontSize',16,'SizeAxesNum',16);
%}

%{
    % specify lambda
    la=[-1 -0.5 0 0.5];
    % Produce a fan plot for each value of \lambda inside vector la
    [out]=FSRfan(y,X,'la',la,'plots',1);
    % Extract in matrix Un the units which entered the search in each step
    Unsel=cell2mat(out.Un);
    lla=length(la);
    nr=size(Unsel,1)/lla;
    Un=[Unsel(1:nr,1) reshape(Unsel(:,2),nr,lla)];
%}

%{
 % construct fan plot specifying the confidence level and the initial
 % starting point for monitoring
    [out]=FSRfan(y,X,'init',size(X,2)+2,'nsamp',0,'conflev',0.95,'plots',1);
%}

%{
    % starting point based on LTS, extraction of all subsamples, construct
    % fan plot specifying the confidence level and the initial starting
    % point for monitoring based on p+2 observations strong line width for
    % lines associated with the confidence bands
    [out]=FSRfan(y,X,'init',size(X,2)+2,'nsamp',0,'lms',0,'lwdenv',3,'plots',1);
%}

%{
    [out]=FSRfan(y,X,'init',4,'nsamp',0);
%}

%{
    % Fan plot using fidelity cards data
    XX=load('loyalty.txt');
    namey='Sales'
    y=XX(:,end);
    nameX={'Number of visits', 'Age', 'Number of persons in the family'};
    X=XX(:,1:end-1);
    % la = vector contanining the most common values of the transformation
    % parameter
    la=[0 1/3 0.4 0.5];
    % Store the score test statistics for the specified values of lambda
    % and automatically produce the fan plot
    [out]=FSRfan(y,X,'la',la,'init',size(X,2)+2,'plots',1,'lwd',3);
%}


%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

% If the number of all possible subsets is <1000 the default is to extract
% all subsets, otherwise just 1000.
ncomb=bc(n,p);
nsampdef=min(1000,ncomb);

% REMARK: a fast approximation of the bc computed above is:
% ncomb=floor(exp( gammaln(n+1) - gammaln(n-p+1) - gammaln(p+1) ) + .5);

hdef=floor(0.5*(n+p+1));
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end

options=struct('la',[-1 -0.5 0 0.5 1],'h',hdef,...
    'nsamp',nsampdef,'lms',1,'plots',0,'init',init,'conflev',0.99,'titl','Fan plot','labx','Subset size m',...
    'laby','Score test statistic','xlimx','','ylimy','','lwd',2,'lwdenv',1,'FontSize',12,'SizeAxesNum',10,...
    'tag','pl_fan','intercept',1,'msg',1,'nocheck',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRfan:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

h=options.h;
lms=options.lms;
plo=options.plots;
nsamp=options.nsamp;
msg=options.msg;

% Specify where to send the output of the current procedure if options plot
% =1
if plo==1
    h1=findobj('-depth',1,'tag',options.tag);
    if (~isempty(h1))
        clf(h1);
        figure(h1)
        axes;
    else
        figure;
        % include specified tag in the current plot
        set(gcf,'tag',options.tag);
    end
end

la=options.la;
init=options.init;
if  init <p+1;
    mess=sprintf(['Attention : init should be larger than p+1. \n',...
        'It is set to p+2.']);
    disp(mess);
    init=p+2;
end


%% Start of the forward search

seq=(1:n)';

%  Unlai is a Matrix whose 2:11th col contains the unit(s) just included.
Unlai = cat(2 , (init+1:n)' , NaN(n-init,10));

% Un = cell which will contain the matrices Unlai for each value of lambda
Un=cell(length(la),1);

lla=length(la);

% Initialize matrix which will contain the score test
Sco=[((init):n)'  NaN(n-init+1,lla)];

% The second column of matrix r will contain the OLS residuals at each step
% of the forward search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100=100*(1:1:ceil(n/100));

binit=zeros(p,lla);

% loop over the values of \lambda
for i=1:lla;
    
    % Construct transformed z according to power tansformation
    if abs(la(i))<1e-8;
        z=log(y);
    else
        z=y.^la(i);
    end
    
    % Find initial subset to initialize the search using as y transformed
    % vector z
    [out]=LXS(z,X,'lms',lms,'h',h,'nsamp',nsamp,'nocheck',1,'msg',msg);
    bsb=out.bs;
    
    % Store information about the units forming subset for each value of
    % lambda
    binit(:,i)=out.bs';
    
    % bsb=[1 8 12 15];
    %ini0 = initial value for forward search loop
    ini0=length(bsb);
    
    % FS loop for a particular value of vector la
    zb=z(bsb);
    yb=z(bsb);
    Xb=X(bsb,:);
    
    % last correctly computed beta oefficients
    blast=NaN(p,1);
    
    if (rank(Xb)~=p)
        warning('FSRfan:message','The provided initial subset does not form full rank matrix');
        % FS loop will not be performed
    else
        for mm=ini0:n;
            % if n>200 show every 100 steps the fwd search index
            if msg==1 && n>200;
                if length(intersect(mm,seq100))==1;
                    disp(['m=' int2str(mm)]);
                end
            end
            
            
            if (mm>=init)
                % Compute and store the value of the score test
                [outSC]=Score(yb,Xb,'la',la(i),'nocheck',1);
                
                % Store S2 for the units belonging to subset
                Sco(mm-init+1,i+1)=outSC.Score;
            end
            
            % Compute b using transformed vector zb
            NoRankProblem=(rank(Xb) == p);
            if NoRankProblem  % rank is ok
                b=Xb\zb;
                blast=b;   % Store correctly computed b for the case of rank problem
            else   % number of independent columns is smaller than number of parameters
                warning('FSR:FSRfan','Rank problem in step %d: Beta coefficients are used from the most recent correctly computed step',mm);
                b=blast;
            end
            
            
            % e= (n x 1) vector of residual for all units using b estimated
            % using subset
            e=z-X*b;
            
            % r_i =e_i^2
            r(:,2)=e.^2;
            
            
            if mm<n;
                
                % store units forming old subset in vector oldbsb
                oldbsb=bsb;
                
                % order the r_i and include the smallest among the units
                %  forming the group of potential outliers
                ord=sortrows(r,2);
                
                % bsb= units forming the new  subset
                bsb=ord(1:(mm+1),1);
                
                Xb=X(bsb,:);  % subset of X
                yb=y(bsb);    % subset of y
                zb=z(bsb);    % subset of z
                
                if mm>=init;
                    unit=setdiff(bsb,oldbsb);
                    if length(unit)<=10
                        Unlai(mm-init+1,2:(length(unit)+1))=unit;
                    else
                        % ALSO INCLUDE VALUE OF LAMBDA
                        disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        Unlai(mm-init+1,2:end)=unit(1:10);
                    end;
                end
            end
            
        end
    end  % rank check
    % Store in cell Un matrix Unlai
    Un{i}=Unlai;
end

%% Structure returned by function FSRfan
out=struct;

out.Score=Sco;
out.la=la;
out.bs=binit;
out.Un=Un;
out.y=y;
out.X=X;

if plo==1
    
    % plot the lines associated with the score test lwd = line width of the
    % trajectories which contain the score test
    lwd=options.lwd;
    plot1=plot(Sco(:,1),Sco(:,2:end),'LineWidth',lwd);
    
    % Specify the line type for the units inside vector units
    slin={'-';'--';':';'-.'};
    slin=repmat(slin,ceil(lla/4),1);
    
        set(plot1,{'LineStyle'},slin(1:lla));

        % set the x and y axis
    xlimx=options.xlimx;
    ylimy=options.ylimy;
    
    if ~isempty(xlimx)
        xlim(xlimx);
    end
    
    if isempty(ylimy)
        
        % Use default limits for y axis
        ylim1=max(-20,min(min(Sco(:,2:end))));
        ylim2=min(20,max(max(Sco(:,2:end))));
        ylim([ylim1 ylim2]);
    else
        
        % Use limits specified by the user
        ylim(ylimy);
    end
    
    % Confidence bands lwdenv = line width of the curves associated with
    % the envelopes
    lwdenv=options.lwdenv;
    conflev=options.conflev;
    v=axis;
    quant=sqrt(chi2inv(conflev,1));
    line([v(1),v(2)],[quant,quant],'color','r','LineWidth',lwdenv);
    line([v(1),v(2)],[-quant,-quant],'color','r','LineWidth',lwdenv);
    
    if size(la,2)>1
        la=la';
    end
    text(n*ones(lla,1),Sco(end,2:end)',num2str(la));
    
    % Main title of the plot and labels for the axes
    labx=options.labx;
    laby=options.laby;
    titl=options.titl;
    
    title(titl);
    
    % FontSize = font size of the axes labels
    FontSize =options.FontSize;
    
    % Add to the plot the labels for values of la Add the horizontal lines
    % representing asymptotic confidence bands
    xlabel(labx,'Fontsize',FontSize);
    ylabel(laby,'Fontsize',FontSize);
    
    % FontSizeAxesNum = font size for the axes numbers
    SizeAxesNum=options.SizeAxesNum;
    set(gca,'FontSize',SizeAxesNum)
    box on
    
end

end

