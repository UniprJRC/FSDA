function [outCp] = FSRcp(y,X,smallp,varargin)
%FSRcp monitors Cp and AIC for all models of interest of size smallp
%
%<a href="matlab: docsearchFS('FSRcp')">Link to the help function</a>
%
% Required input arguments:
%
% y:            Response variable. Vector. A vector with n elements that contains the response variable.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
% X :           Predictor variables. Matrix. Data matrix of explanatory variables (also called
%               'regressors') of dimension (n x (bigP-1)).
%               The intercept will be added in automatic way, so that the
%               dimension of the full model is bigP
%               Rows of X represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
% smallp:       number of variables in the reduced models. Scalar. Scalar
%               which specifies the number of variables in the
%               reduced models which will be considered. For example if
%               smallp=3, all possible subsets containing 2 columns of
%               matrix X will be considered. Notice that the dimension of
%               each submodel is 3 because to each submodel the column of
%               ones is added automatically.
%
% Optional input arguments:
%
%    intercept :  Indicator for constant term. true (default) | false. 
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%     nocheck   : Check input arguments. Scalar.
%                        If nocheck is equal to 1 no check is performed on
%                       matrix y and matrix X. Note that y and X are left
%                       unchanged. In other words the additioanl column of ones
%                       for the intercept is not added. As default nocheck=1.
%                       Example - 'nocheck',1
%                       Data Types - double
%           h   :       number of observations that have determined the least
%                       trimmed squares estimator. Integer.
%                       h is an integer greater than
%                       smallp+1 but smaller then n. The default value of h is
%                       [(n+smallp+1)/2]
%                       Example - 'h',3
%                       Data Types - double
%           lms :    Criterion to use to find the initlal  subset to
%                       initialize the search. Scalar. If lms=1 (default) Least Median of Squares is
%                       computed, else Least Trimmed of Squares is computed.
%                       Example - 'lms',1
%                       Data Types - double
%          nomes:       Displaying time message. Scalar. If nomes is equal to 1 (default) no message about
%                       estimated time to compute LMS (LTS) for each considered
%                       model is displayed, else a message about estimated time
%                       is displayed.
%                       Example - 'lms',1
%                       Data Types - double
%         nsamp : Number of subsamples which will be extracted to find the
%                       robust estimator. Scalar.
%                       If nsamp=0 all subsets will be extracted.
%                       They will be (n choose smallp).
%                       Example - 'nsamp',1000
%                       Data Types - double
%                       Remark: if the number of all possible subset is <1000 the
%                       default is to extract all subsets otherwise just 1000.
%          init :       Search initialization. Scalar.
%                       It specifies the initial subset size to start
%                       monitoring the required quantities, if init is not
%                       specified it will be set equal to
%                       smallp+1, if the sample size is smaller than 40;
%                       min(3*smallp+1,floor(0.5*(n+smallp+1))), otherwise.
%                       Example - 'init',100 starts monitoring from step m=100
%                       Data Types - double
%           aic :      Akaike's information criterion. Scalar.
%                       If aic=1 the value of AIC is also stored in each
%                       step of the search else (default) only Mallows Cp is stored
%                       Example - 'aic',1
%                       Data Types - double
%        labels :   names of the explanatory variables. Cell array of strings.
%                       Cell array of strings of length bigP-1 containing the  names of the explanatory variables.
%                       If labels is a missing  value the following sequence of strings will be
%                       automatically created for X
%                      (1,2,3,4,5,6,7,8,9,A,B,C,D,E,E,G,H,I,J,K,...,Z)
%                       Example - 'labels',{'Time','1','2','3','4','5','6','7','8'}
%                       Data Types - cell
%      fin_step :  portion of the search which has to be
%                 monitored to choose the best models. Scalar.
%                 If fin_step is an integer greater
%                 or equal 1, it refers to the number of steps.
%                 For example if fin_step=10 the program considers the last
%                 10 steps to choose the best models.
%                 If fin_step is a real number alpha (0<alpha<0.5) in the
%                 interval (0 0.5] than the program considers the last
%                 round(n*alpha) steps. As default fin_step=round(n*0.2)
%                 that is the last 20% of the steps are considered.
%                       Example - 'fin_step',1
%                       Data Types - double
%                 Remark1: the number of best models to consider is
%                 controlled by scalar first_k (see below).
%                 Remark2: if fin_step is an empty value, no selection is
%                 done and all trajectories of Cp are displayed (in this
%                 case the value of first_k below is ignored, all models are
%                 considered of interest and output matrix outCp.Ajout is
%                 equal to an empty value).
%       first_k :  number of best models to
%                      consider in each of the last fin_step. Scalar.
%                       For example if first_k=5 in each of the last fin_step, the models which had
%                       the 5 smallest values of Cp are considered. As default
%                       first_k=3
%                       Example - 'first_k',5
%                       Data Types - double
%          Excl : Matrix which contains the models which surely do not have
%                 to be considered. Matrix.
%                   As default Excl=''
%                   For example if smallp=3, bigP=6 and
%                   Excl = [23; 24; 27]; the three models 23, 24, and 27 are skipped
%                       Example - 'Excl',[23; 24]
%                       Data Types - double
%    ExclThresh : Exclusion threshold. Scalar.
%                 Exclusion threshold associated to the upper
%                 percentage point of the F distribution of Cp which
%                 defines the threshold for declaring models as unuseful.
%                 The default value of ExclThresh is 0.99999 that is the
%                 models whose minimum value of Cp in the part of the
%                 search defined by fin_step is above ExclThresh are stored
%                 in output matrix outCp.Ajout. Notice that ExclThresh must
%                 be smaller than 1
%                 Example - 'ExclThresh',0.6
%                 Data Types - double
%         plots : Plot on the screen. Scalar.
%                 If plots==1 a plot is created on the screen which
%                 contains the trajectories of Cp monitored along the
%                 search with confidence bands
%                 If plots==2 two plots are generated. The first contains
%                 the trajectories of Cp monitored along the search with
%                 confidence bands. The second contains the trajectories of
%                 AIC monitored along the search
%                 else (default) no plot is shown on the screen
%                 Example - 'plots',1
%                 Data Types - double
%        labout :If labout=1 the output LABOUT contains the list of models
%                 whose Cp values are inacceptable. Scalar. Default: no
%                 model is created.
%                 Example - 'labout',1
%                 Data Types - double
%                 Remark: the options below only work if plots is equal 1.
%
%         quant : It specifies the quantiles which are used to
%                 produce Cp envelopes. Vector.
%                The elements of quant are numbers
%                 between 0 and 1. The default value of quant is
%                 quant=[0.025 0.5 0.975];
%                 Example - 'quant',0.1
%                 Data Types - double
%         steps : Steps to add labels. Vector. It specifies in which steps of the plot which
%                 monitors Cp it is necessary to include the labels of the
%                 models which have been previously chosen. 
%                 The default is to write the labels of the models in steps
%                 round([n*0.6  n*0.8  n]);
%                 Example - 'steps',[4 8]
%                 Data Types - double
%       titl    : a label for the title. Character.
%               default is ['Forward Cp' p= num2str(smallp)]
%                 Example - 'titl','my title'
%                 Data Types - char
%       labx    : a label for the x-axis. Character.
%                   default: 'Subset size m'
%                 Example - 'labx','my label'
%                 Data Types - double
%       laby    : a label for the y-axis. Character.
%                   default:''
%                 Example - 'laby','my label'
%                 Data Types - char
%       xlimx   :  minimum and maximum on the x axis. Vector.
%                 Default value is '' (automatic scale)
%                 Example - 'xlimx',[0 1]
%                 Data Types - double
%       ylimy   : minimum and maximum on the y axis. Vector.
%                 Default value is '' (automatic scale)
%                 Example - 'ylimx',[0 1]
%                 Data Types - double
%       lwd     : linewidth of the curves which contain the score test.
%                   Scalar.
%                 Default line width=2
%                 Example - 'linewidth',6
%                 Data Types - double
%       lwdenv  :  width of the lines associated
%                 with the envelopes. Scalar.
%                  Default is lwdenv=1
%                 Example - 'lwdenv',6
%                 Data Types - double
%       FontSize: font size of the labels of
%                 the axes and of the labels inside the plot. Scalar.
%                 Default value is 12
%                 Example - 'FontSize',20
%                 Data Types - double
%    SizeAxesNum: size of the numbers of the axes. Scalar.
%                 Default value is 10
%                 Example - 'SizeAxesNum',30
%                 Data Types - double
%   selunitcolor: colors to be used for the Cp trajectories. Cell array of strings.
%                   If selunittype is not specified or if
%                 it is an empty value default Matlab colors are used.
%                 Example - 'selunitcolor',{'b';'g';'r'}
%                 Data Types - cell
%   selunittype : line types of the Cp trajectories. Cell array of strings.
%                 If selunittype is not specified or if
%                 it is an empty value all possible line styles are used.
%                 Example - 'selunittype',{'-';'--';':';'-.'}
%                 Data Types - cell
% Output:
%
%  The output consists of a structure 'outCp' containing the following fields:
%
%         outCp.MAL = (n-init+1) x (k+1) matrix.
%                 Mallows Cp monitored along the search:
%                   1st col is fwd search index;
%                   2nd col is associated with first selected model;
%                   3rd col is associated with second selected model;
%                   ...;
%                   (k+1)th col is associated with k-th selected model.
%                   Notice that k<=(n choose smallp) and that all
%                   models contain the constant.
%         outCp.AIC = (n-init+1) x (k+1) matrix.
%                AIC monitored along the search:
%                   1st col is fwd search index;
%                   2nd col is associated with first selected model;
%                   3rd col is associated with second selected model;
%                   ...;
%                   (k+1)th col is associated with k-th selected model.
%                   Remark 1: k<=(n choose smallp).
%                   Remark 2: all models contain the constant.
%                   Remark 3: matrix AIC is produced only if input option
%                   aic=1.
%      outCp.UnAll =    cell of dimension k. Each element of the cell is a
%                   (n-init) x 11 matrix containing the unit(s) included
%                   in the subset at each step of the search.
%                   REMARK: in every step the new subset is compared with the old
%                   subset. Un contains the unit(s) present in the new
%                   subset but not in the old one.
%      outCp.LAB    =    cell array of strings of length k containing the labels of the
%                   models which have been extracted. First element of LAB
%                   is associated with second column of matrix MAL...
%     outCp.Ajout  =    numeric matrix which contains the list of the models whose Cp
%                   values are inacceptable.
%                   The number of columns of matrix Ajout is equal to
%                   smallp-1
%                   This information is useful because in this way it is
%                   possible to skip the computation of the submodels of
%                   the rows of matrix Ajout.
%                   For example if smallp=3, bigP=6 and
%                   Ajout = [ 23; 24; 27 ]
%                   the three models 23, 24, and 27 always have Cp values
%                   much greater than the threshold (that is variables
%                   2,3,4,7 are considered unimportant).
%    outCp.LABout  = cell array of strings which contains as
%                   strings the list of models which are unacceptable.
%                   LABout is created only if input option labout=1.
%
% See also: FSR, FSReda
%
% References:
%
% Atkinson, A.C. and Riani, M. (2008), A robust and diagnostic information
% criterion for selecting regression models, "Journal of the Japan
% Statistical Society", Vol. 38, pp. 3-14.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRcp')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % FSRcp with all default options.
    % Extract the best models of size 4.
    % Common part to all examples: load Ozone dataset.
    X=load('ozone.txt');
    % Transform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    smallp=4;
    [Cpmon]=FSRcp(y,X,smallp);
%}

%{
    %% FSRcp with optional arguments.
    % Extract the best models of size 4, and show the plot
    % of forward Cp.
    X=load('ozone.txt');
    % Transform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    smallp=4;
    [outCp]=FSRcp(y,X,smallp,'plots',1);
%}

%{
    % Use labels defined by the user.
    % Extract the best models of size 4 and show the plot of Cp. All the
    % default options are used, apart from labels, therefore the plot of Cp
    % and the output matrix Cpmon.MAL only contains the searches associated
    % with the smallest 3 values of Cp in the last 16 steps of the search.
    X=load('ozone.txt');
    % Transform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    smallp=4;
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpmon]=FSRcp(y,X,smallp,'plots',1,'labels',labels);
%}


%{
    % Extract and show the trajectories of all models of size 4 of Cp.
    % Notice that in this last case the forward plot becomes
    % unreadable.
    X=load('ozone.txt');
    % Transform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    smallp=4;
    [Cpmon]=FSRcp(y,X,smallp,'plots',1,'fin_step','');
%}

%{
    % Extract the best models of size 5 and plot monitoring of Cp.
    % Extract 1000 subsets to initialize the search and use labels defined
    % by the user.
    X=load('ozone.txt');
    % Transform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    smallp=5;
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpmon]=FSRcp(y,X,smallp,'nsamp',1000,'plots',1,'labels',labels);
%}

%{
    %% Extract the best models of size 6 and 5 and plot monitoring of Cp.
    % Extract 1000 subsets to initialize the search andse labels defined by
    % the user. Exclude the searches of the models which were unacceptable
    % for smallp=5.
    X=load('ozone.txt');
    % Transform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    smallp=6;
    labels={'Time','1','2','3','4','5','6','7','8'};
    [Cpmon6]=FSRcp(y,X,smallp,'nsamp',1000,'plots',1,'labels',labels);
    smallp=5;
    [Cpmon5]=FSRcp(y,X,smallp,'nsamp',1000,'Excl',Cpmon6.Ajout,'plots',1,'labels',labels);
%}

%{
    % Customizing the graphical options.
    % In the following example we play with the graphical options
    X=load('ozone.txt');
    % Transform the response using logs
    X(:,end)=log(X(:,end));
    % Add a time trend
    X=[(-40:39)' X];
    % Define y
    y=X(:,end);
    % Define X
    X=X(:,1:end-1);
    labels={'Time','1','2','3','4','5','6','7','8'};
   [Cpmon]=FSRcp(y,X,smallp,'plots',1,'labels',labels,'xlimx',[40 80],'lwdenv',5,'lwd',4,'FontSize',25,'SizeAxesNum',20);
%}


%% Beginning of code

% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% If the number of all possible subsets is <1000 the default is to extract
% all subsets otherwise just 1000.
ncomb=bc(n,smallp);
nsampdef=min(1000,ncomb);
% Note that a fast approximation of the bc computed above is:
% ncomb=floor(exp( gammaln(n+1) - gammaln(n-p+1) - gammaln(p+1) ) + .5);

% specify the steps of the search in which it is necessary to write labels
% of the models
hdef=floor(0.5*(n+smallp+1));

if n<40
    init=smallp+1;
else
    init=min(3*smallp+1,floor(0.5*(n+smallp+1)));
end
steplabdef=round([hdef n*0.8 n]);

options=struct('intercept',1,'h',hdef,...
    'nsamp',nsampdef,'lms',1,'init',init,'nomes',1,...
    'labels','','fin_step',round(n*0.2),'first_k',3,'aic',0,...
    'selunitcolor','','selunittype','',...
    'plots',0,'quant',[0.025 0.5 0.975],...
    'steps',steplabdef,'Excl','','labout',0,'ExclThresh',0.99999,...
    'titl','','labx','Subset size m','laby','',...
    'xlimx','','ylimy','','lwd',2,'lwdenv',1,'FontSize',12,'SizeAxesNum',10,'nocheck',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRcp:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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


init=options.init;
h=options.h;
lms=options.lms;
nsamp=options.nsamp;
% Information about plotting or not the trajectories
plots=options.plots;

% Information about storing AIC or not
aic=options.aic;

nomes=options.nomes;
%% Beginning of procedure

% Information about the labels of the variables
labels=options.labels;

% Create the labels of the different models which have been considered
Aj=combsFS(1:(p-1),smallp-1);
LAB=cellstr(num2str(Aj,'%d,'));

if isempty(labels)
    % Default labels for the explanatory variables
    labels={'1','2','3','4','5','6','7','8','9','A','B','C',...
        'D','E','F','G','H','I','J','K','L','M','N','O',...
        'P','Q','R','S','T','U','V','X','Y','Z'};
end

lLAB=length(LAB{1,1});
for j=1:size(LAB,1)
    LAB{j,1}(lLAB)=[];
end

for i=1:p-1
    LAB=regexprep(LAB,num2str(i),labels(i));
end

Excl=options.Excl;

[rAj,cAj]=size(Aj);
if ~isempty(Excl)
    
    ajsel=1:rAj;
    for i=1:rAj
        for j=1:size(Excl,1)
            % If the following condition is fulfilled means that the model
            % including the previous set of explanatory variables had
            % already been analyzed and has been found unuseful
            if length(intersect(Aj(i,:),Excl(j,:)))==cAj
                ajsel(i)=0;
                break;
            else
            end
        end
    end
    % Redefine matrix Aj and vector of strings LAB to include just the
    % searches for the models which are relevant
    ajsel=ajsel(ajsel>0);
    Aj=Aj(ajsel,:);
    
    if ~isempty(Aj)
        rAj=size(Aj,1);
        LAB=LAB(ajsel);
    end
end

% No model of interest for this size of smallp
if ~isempty(Aj)
    Aj=Aj+1;
    
    %% Initialise key matrices
    
    % sequence from 1 to n.
    seq=(1:n)';
    
    % The second column of matrix R will contain the OLS residuals at each step
    % of the search
    r=[seq zeros(n,1)];
    
    % If n is very large, the step of the search is printed every 100 step
    % seq100 is linked to printing
    seq100=100*(1:1:ceil(n/100));
    
    
    % Matrix MAL will contain the monitoring of Mallows Cp in each step of the fwd
    % search.
    MAL=[(init:n)' zeros(n-init+1,rAj)];
    
    
    % Matrix AIC will contain the monitoring of AIC in each step of the fwd
    % search.
    if aic
        AIC=MAL;
    end
    
    UnAll=cell(rAj,1);
    
    %  Un is a Matrix whose 2nd column:11th cols contain the unit(s) just
    %  included.
    Un = cat(2 , (init+1:n)' , NaN(n-init,10));
    
    
    % Information about the number of last steps which have to be chosen to extract the
    % models with best Cp
    fin_step=options.fin_step;
    
    if fin_step<0.5
        fin_step=round(n*fin_step);
    end
    
    % Variable subset loop
    for jj=1:rAj
        
        sel=[1 (Aj(jj,:))];
        Xred=X(:,sel);
        
        % Find initial subset to initialize the search using reduced subset of
        % variables
        [out]=LXS(y,Xred,'lms',lms,'h',h,'nsamp',nsamp,'nocheck',1,'nomes',nomes);
        bsb=out.bs;
        Xb=Xred(bsb,:); % Subset of X using reduced set of expl. variables
        Xball=X(bsb,:); % Subset of X using full set of expl. variables
        yb=y(bsb);
        
        
        %% Forward search loop
        
        for mm=smallp:n
            % if n>200 show every 100 steps the fwd search index
            if n>500
                if length(intersect(mm,seq100))==1
                    disp(['m=' int2str(mm)]);
                end
            end
            
            
            bred=Xb\yb;
            
            
            rred=yb-Xb*bred;
            sred=rred'*rred;
            
            % ols on the same
            % subset but on the full set of explanatory variables
            ball=Xball\yb;
            rall=yb-Xball*ball;
            sall=rall'*rall/(mm-p);
            
            if mm>=init
                % Now compute AIC criterion for that particular subset of data
                if aic==1
                    AIC(mm-init+1,jj+1)=mm*log(2*pi) +mm*sall  +sred/sall+2*smallp;
                end
                MAL(mm-init+1,jj+1)=sred/sall-mm+2*smallp;
            end
            
            % e = vector of residuals for all units using b estimated using subset
            e=y-Xred*bred;
            
            r(:,2)=e.^2;
            
            if mm<n
                
                % store units forming old subset in vector oldbsb
                oldbsb=bsb;
                
                % order the r_i and include the smallest among the units forming
                % the group of potential outliers
                ord=sortrows(r,2);
                
                % bsb= units forming the new  subset
                bsb=ord(1:(mm+1),1);
                
                Xb=Xred(bsb,:);  % subset of X using reduced set of variables
                Xball=X(bsb,:);
                yb=y(bsb);    % subset of y
                
                if mm>=init
                    unit=setdiff(bsb,oldbsb);
                    % If the interchange involves more than 10 units, store only the
                    % first 10.
                    if length(unit)<=10
                        Un(mm-init+1,2:(length(unit)+1))=unit;
                    else
                        disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        disp(['Number of units which entered=' int2str(length(unit))]);
                        
                        Un(mm-init+1,2:end)=unit(1:10);
                    end
                end
                
                
            end
        end
        
        UnAll{jj}=Un;
    end
    
    
    
    
    if ~isempty(fin_step)
        % Find in the last fin_steps what are the searches with the smallest
        % first_k Mallows Cp
        [~,IX]=sort(MAL(end-fin_step(1)+1:end,2:end),2);
        
        % Information about the number of best searches to consider
        first_k=options.first_k;
        
        first=min(first_k,size(IX,2));
        sto=reshape(IX(:,1:first),fin_step(1)*first,1);
        % units vector which contains the information about the
        % searches which satisfy the previous criteria
        units=unique(sto);
        % lunits = number of selected searches
        lunits=length(units);
    else
        % all trajectories are considered
        units=(1:(size(MAL,2)-1))';
        
        % lunits = number of searches which must be labelled refers to all
        % searches
        lunits=length(LAB);
    end
    
    outCp=struct;
    outCp.MAL=MAL(:,[1;(units+1)]);
    Unsel=UnAll(units);
    outCp.UnAll=Unsel;
    outCp.LAB=LAB(units);
    
    % Plotting part
    if plots==1 || plots==2
        
        % Quantiles which are used for the envelopes of Cp
        quant=options.quant;
        
        % Compute and superimpose envelopes based on the F distribution
        EnvF=[MAL(:,1) zeros(size(MAL,1),length(quant))];
        
        for i=1:size(MAL,1)
            ast=finv(quant,p-smallp,MAL(i,1)-p);
            ast=ast*(p-smallp)+2*smallp-p;
            EnvF(i,2:end)=ast;
        end
        
        
        numtext=LAB;
        
        % Information about the steps in which labels have to be added
        steps=options.steps;
        
        % Check that the values of steps are all greater than the smallest
        % value of the subset size (m0) and smaller or equal than n
        
        residuals=MAL(:,2:end)';
        x=MAL(:,1);
        
        if max(steps)>n
            fprintf(['One of the steps which has beeen chosen is greater than n. \n',...
                'It is deleted.']);
            steps=steps(steps<=n);
        end
        if min(steps)<x(1)
            fprintf(['One of the steps which has beeen chosen is smaller than m0. \n',...
                'It is deleted.']);
            steps=steps(steps>=x(1));
        end
        
        
        % lsteps = number of steps for which it is necessary to add the labels
        lsteps=length(steps);
        % lall = number of steps * number of models to label
        lall=lunits*lsteps;
        
        figure;
        
        % lwd = line width of the trajectories of Cp
        lwd=options.lwd;
        
        % plot the trajectories for models specified in units
        plot1=plot(x,residuals(units,:),'tag','Cp','LineWidth',lwd);
        
        % SET SOME FIGURE PROPERTIES OF THE PLOT WHICH MONITORS Cp
        % FontSize = font size of the axes labels
        FontSize =options.FontSize;
        
        
        
        % Labels of the models in particular steps of the search specified
        % in vector steps
        text(reshape(repmat(steps,lunits,1),lall,1),reshape(residuals(units,steps-x(1)+1),lall,1),...
            reshape(repmat(numtext(units),1,lsteps),lall,1),'FontSize',FontSize)
        
        hold('on');
        % Set line width of the curves which represent the envelopes
        lwdenv=options.lwdenv;
        line(EnvF(:,1),EnvF(:,2:end),'LineStyle','--','Color','r','LineWidth',lwdenv);
        
        % Add text associated with the envelopes which have been used
        text(repmat(n,length(quant),1),EnvF(end,2:end)',num2str(quant'));
        
        
        % Main title of the plot and labels for the axes
        labx=options.labx;
        laby=options.laby;
        titl=options.titl;
        
        if ~isempty(titl)
            title(titl);
        else
            % Title
            title(['Forward Cp, p='  num2str(smallp)]);
        end
        
        
        % Add to the plot the labels for values of la
        % Add the horizontal lines representing asymptotic confidence bands
        xlabel(labx,'Fontsize',FontSize);
        ylabel(laby,'Fontsize',FontSize);
        
        % SizeAxesNum = font size for the axes numbers
        SizeAxesNum=options.SizeAxesNum;
        set(gca,'FontSize',SizeAxesNum)
        
        % set the x and y axis
        xlimx=options.xlimx;
        ylimy=options.ylimy;
        
        if ~isempty(xlimx)
            xlim(xlimx);
        end
        
        if ~isempty(ylimy)
            %     % Use default limits for y axis
            %     ylim1=max(-20,min(min(Sco(:,2:end))));
            %     ylim2=min(20,max(max(Sco(:,2:end))));
            %     ylim([ylim1 ylim2]);
            % else
            %     % Use limits specified by the user
            ylim(ylimy);
        end
        
        % Specify the line type for the units inside vector units
        slintyp=options.selunittype;
        
        if ~isempty(slintyp)
            slintyp=repmat(slintyp,ceil(length(units)/length(slintyp)),1);
            set(plot1(units),{'LineStyle'},slintyp(units));
        else
            slintyp={'-';'--';':';'-.'};
            % slintyp={'-'};
            slintyp=repmat(slintyp,ceil(length(units)/length(slintyp)),1);
            % NOTE NOTE set(plot1,{'Line'},slintyp(1:length(units)));
            set(plot1,{'LineStyle'},slintyp(1:length(units)));
        end
        
        % if requested, set the color of the selected trajectories
        % note that if selunitcolor contains more than one color, e.g.
        % options.selunitcolor = {'b';'g';'r'},
        % then the colors of the trajectories alternate.
        scol=options.selunitcolor;
        if ~isempty(scol)
            scol=repmat(scol,ceil(lunits/length(scol)),1);
            set(plot1,{'Color'},scol(1:lunits));
        end
        
        
        if plots==2
            if ~isempty(fin_step)
                % Find in the last fin_steps which are the searches with the smallest
                % first_k Mallows Cp
                [~,IX]=sort(AIC(end-fin_step(1)+1:end,2:end),2);
                
                % Information abou the number of best searches to consider
                first_k=options.first_k;
                
                first=min(first_k,size(IX,2));
                sto=reshape(IX(:,1:first),fin_step(1)*first,1);
                % units vector which contains the information abou the
                % searches which satisfy the previous criteria
                units=unique(sto);
                % lunits = number of selected searches
                lunits=length(units);
            else
                % all trajectories are considered
                units=1:(size(AIC,2)-1);
                
                % lunits = number of searches which must be labelled refers to all
                % searches
                lunits=length(LAB);
            end
            
            numtext=LAB;
            
            % lsteps = number of steps for which it is necessary to add the labels
            lsteps=length(steps);
            % lall = number of steps * number of models to label
            lall=lunits*lsteps;
            
            figure;
            residuals=AIC(:,2:end)';
            x=AIC(:,1);
            
            % plot the trajectories for models specified in units
            plot(x,residuals(units,:),'tag','Cp');
            
            title(['Forward AIC, p='  num2str(smallp)]);
            text(reshape(repmat(steps,lunits,1),lall,1),reshape(residuals(units,steps-x(1)+1),lall,1),reshape(repmat(numtext(units),1,lsteps),lall,1))
            
        end
    end
    
    if aic==1
        outCp.AIC=AIC;
    end
    
    
    if ~isempty(fin_step)
        % Construct the list of the models which surely contain no information
        % that is those which in the final part of the search always have a value
        % of Cp greater than the ExclThresh quantile
        
        cor=((p-smallp)+2*smallp-p);
        ast=finv(options.ExclThresh,p-smallp,n-p)*cor;
        
        Ajout=Aj(min(MAL(end-fin_step(1)+1:end,2:end))>ast,:)-1;
        LABout=LAB(min(MAL(end-fin_step(1)+1:end,2:end))>ast);
        
        % Now select models which are clearly good and therefore their
        % nested models must be considered for sure
        %astmin=cor*finv(0.975,p-smallp,n-p);
        %Ajsel=Aj(max(MAL(end-fin_step(1)+1:end-5,2:end))<astmin,:)-1;
        
        outCp.Ajout=Ajout;
        if options.labout
            outCp.LABout=LABout;
        end
    else
        % all models are considered of interst
    end
else
    % No model of interest
    % In this case the output is an empty structure with no fields
    
    outCp=struct([]);
    % In this case the models to exclude are all
    % outCp.Ajout=Excl;
end
end
%FScategory:REG-ModelSelection