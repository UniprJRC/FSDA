function [out]=FSRfanCens(y,X,varargin)
%FSRfanCens monitors the values of the signed sqrt LR test for each lambda in the tobit model
%
%<a href="matlab: docsearchFS('FSRfanCens')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that
%               contains the response
%               variable. It can be either a row or a column vector.
%    X :        Predictor variables. Matrix. Data matrix of explanatory
%               variables (also called 'regressors') of dimension (n x
%               p-1). Rows of X represent observations, and columns
%               represent variables. Missing values (NaN's) and infinite
%               values (Inf's) are allowed, since observations (rows) with
%               missing or infinite values will automatically be excluded
%               from the computations. NOTICE THAT THE INTERCEPT MUST
%               ALWAYS BE INCLUDED.
%
%  Optional input arguments:
%
%
% balancedSearch:   Balanced search. Scalar logical.
%                   If Balanced search the proportion of observations in
%                   the subsets equals (as much as possible) the proportion
%                   of units in the sample. The default value of
%                   balancedSearch is true.
%                   Example - 'balancedSearch',false
%                   Data Types - logical
%
%        family :   string which identifies the family of transformations which
%                   must be used. Character. Possible values are 'BoxCox'
%                   or 'YJ' (default)
%                   The Box-Cox family of power transformations equals
%                   $(y^{\lambda}-1)/\lambda$ for $\lambda$ not equal to zero,
%                   and $\log(y)$ if $\lambda = 0$.
%                   The Yeo-Johnson (YJ) transformation is the Box-Cox
%                   transformation of $y+1$ for nonnegative values, and of
%                   $|y|+1$ with parameter $2-\lambda$ for $y$ negative.
%                   Remember that BoxCox can be used just
%                   if input left is positive. Yeo-Johnson family of
%                   transformations does not have this limitation.
%                   Example - 'family','YJ'
%                   Data Types - char
%
%
%       init    :   Search initialization. Scalar.
%                   It specifies the initial subset size to start
%                   monitoring the value of the score test, if init is not
%                   specified it will be set equal to:
%                    p+1, if the sample size is smaller than 40;
%                    min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%                    Example - 'init',100 starts monitoring from step m=100
%                    Data Types - double
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%           la  :   values of the transformation parameter for which it is
%                   necessary to compute the LR test. Vector.
%                   Default value of lambda is la=[-1 -0.5 0 0.5 1]; that
%                   is the five most common values of lambda
%                   Example - 'la',[-1 -0.5]
%                   Data Types - double
%
%    left :     left limit for the censored dependent variable. Scalar.
%               If set to -Inf, the dependent variable is assumed to be not
%               left-censored; default value of left is zero (classical
%               Tobit model).
%               Example - 'left',1
%               Data Types - double
%
%       lms     :   Criterion to use to find the initial
%                 subset to initialize the search. Scalar.
%                 If lms=1 (default) Least Median of Squares is
%                 computed, else Least Trimmed Squares is computed.
%                   Example - 'lms',0
%                   Data Types - double
%
%         msg   : Level of output to display. Boolean.
%                   Boolean scalar which controls whether to display or not
%                   messages on the screen.
%                   If msg==true (default) messages are
%                   displayed on the screen about estimated time to compute
%                   the LMS (LTS) for each value of lambda else no message
%                   is displayed on the screen
%                  Example - 'msg',false
%                  Data Types - logical
%
%       nocheck :   Check input arguments. Boolean.
%                   If nocheck is equal to true no check is performed
%                   on matrix y and matrix X. Notice that y and X are left
%                   unchanged. In other words the additional column of ones
%                   for the intercept is not added. As default nocheck=false.
%                   Example - 'nocheck',true
%                   Data Types - boolean
%
%        right :    right limit for the censored dependent variable. Scalar.
%                   If set to Inf, the dependent variable is assumed to be not
%                   right-censored; default value of left is Inf (classical
%                   Tobit model).
%                   Example - 'right',800
%                   Data Types - double
%
%       nsamp   :   Number of subsamples which will be extracted to find
%                   the initiasl robust estimator. Scalar.
%                   If nsamp=0 all subsets will be extracted. They will be
%                   (n choose p). Remark: if the number of all possible
%                   subset is <1000, the default is to extract all subsets
%                   otherwise just 1000. If nsamp is a matrix of size
%                   r-by-p, it contains in the rows the subsets which sill
%                   have to be extracted. For example, if p=3 and nsamp=[ 2
%                   4 9; 23 45 49; 90 34 1]; the first subset is made up of
%                   units [2 4 9], the second subset of units [23 45 49]
%                   and the third subset of units [90 34 1];
%                   Example - 'nsamp',1000
%                   Data Types - double
%
%
%       plots   :  Plot on the screen. Scalar.
%                   If plots=1 the fan plot is produced
%                   else (default) no plot is produced.
%                   Example - 'plots',1
%                   Data Types - double
%                   REMARK: all the following options work only if plots=1
%
%
%       conflev :   Confidence level. Scalar or vector. Confidence level
%                   for the bands (default is 0.99, that is we plot two
%                   horizontal lines in correspondence of value -2.58 and
%                   2.58).
%                   Example - 'conflev',[0.9 0.95 0.99]
%                   Data Types - double
%
%       FontSize:   font size of the labels of the axes. Scalar.
%                   Default value is 12.
%                   Example - 'FontSize',20
%                   Data Types - double
%
%       labx    :   a label for the x-axis. Character.
%                   default: 'Subset size m'
%                   Example - 'labx','my labx'
%                   Data Types - char
%
%       laby    :   a label for the y-axis. Character.
%                   default:'Score test statistic'
%                   Example - 'laby','my laby'
%                   Data Types - char
%
%       lwd     :   linewidth of the curves which
%                   contain the score test. Scalar.
%                   Default line width=2.
%                   Example - 'lwd',5
%                   Data Types - double
%
%       lwdenv  :   width of the lines associated
%                   with the envelopes. Scalar.
%                   Default is lwdenv=1.
%                   Example - 'lwdenv',5
%                   Data Types - double
%
%    SizeAxesNum:   Scalar which controls the size of the numbers of the
%                   axes. Scalar.
%                   Default value is 10.
%                  Example - 'SizeAxesNum',12
%                  Data Types - double
%
%       tag     :   handle of the plot which is about to be created.
%                   Character. The default is to use tag 'pl_fan'. Notice
%                   that if the program finds a plot which has a tag equal
%                   to the one specified by the user, then the output of
%                   the new plot overwrites the existing one in the same
%                   window else a new window is created Example -
%                   'tag','mytag' Data Types - char
%
%       titl    :   a label for the title. Character.
%                   default: 'Fan plot'
%                   Example - 'titl','my title'
%                   Data Types - char
%
%       xlimx   :   Minimum and maximum of the x axis. Vector.
%                   Default value is [init n]
%                   Example - 'xlimx',[0 1]
%                   Data Types - double
%
%       ylimy   :  Minimum and maximum of the y axis. Vector.
%                   Default value for ylimy(1)=max(min(score_test),-20).
%                   Default value for ylimy(2)=min(max(score_test),20).
%                   Example - 'ylimy',[0 1]
%                   Data Types - double
%
%  Output:
%
%         out:   structure which contains the following fields
%
%  out.Score  = (n-init+1) x length(la)+1 matrix containing the values of the
%               signed sqrt LR test for each value of the transformation parameter:
%               1st col = fwd search index;
%               2nd col = value of the signed sqrt LR test in each step of the
%               fwd search for la(1);
%               ...........
%               end col = value of the signed sqrt LR test test in each step of the fwd
%               search for la(end).
%   out.Exflag = (n-init+1) x *length(la)+1 matrix containing the
%               reason fminunc stopped in the maximization of the
%               unconstrained likelihood. Matrix of integers.
%               Column 2 is associated with the search which has
%               ordered the data using la(1);
%               .........
%               Column *length(la)+1 is associated with
%               the search which has ordered the data using
%               la(length(la)).
%    out.laMLE = (n-init+1) x *length(la)+1 matrix containing the values of the
%               maximum likelihood estimate of lambda.
%               Column 2 is associated with the search which has
%               ordered the data using la(1);
%               .........
%               Column *length(la)+1 is associated with
%               the search which has ordered the data using
%               la(length(la)).
%  out.la     = vector containing the values of lambda for which fan plot
%               is constructed
%  out.bs     = matrix of size p x length(la) containing the units forming
%               the initial subset for each value of lambda
%  out.Un     = cell of size length(la).
%               out.Un{i} is a n-init) x 11 matrix which contains the
%               unit(s) included in the subset at each step in the search
%               associated with la(i).
%               REMARK: in every step the new subset is compared with the old subset. Un
%               contains the unit(s) present in the new subset but not in
%               the old one. Un(1,:), for example, contains the unit included
%               in step init+1 ... ; Un(end,2) contains the units included in the
%               final step of the search
%  out.y      = A vector with n elements that contains the response
%               variable which has been used
%  out.X      = Data matrix of explanatory variables
%               which has been used (it also contains the column of ones if
%               input option intercept was missing or equal to 1)
%   out.class = 'FSRfanCens'.
%
%
% See also: regressCensTra, regressCens, FSRfan, Score, ScoreYJ, ScoreYJpn, fanBIC, fanBICpn
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRfanCens')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% FSRfanCens with simulated data.
    % Generate Tobit regression data
    rng('default')
    rng(200)
    n=200;
    lambda=0.5;
    p=3;
    sigma=0.5;
    beta=ones(p,1);
    X=randn(n,p)+1;
    epsilon=randn(n,1);
    
    y=X*beta+sigma*epsilon;
    y=normYJ(y,1,lambda,'inverse',true,'Jacobian',true);
    qq=quantile(y,0.3);
    y(y<=qq)=qq;
    left=min(y);
    % Compute fanplot for 3 values of lambda
    [out]=FSRfanCens(y,X,'left',left,'la',[0 0.5 1],'init',round(n*0.5));
%}

%{
    %% Double censored simulated data.
    rng('default')
    rng(4)
    n=200;
    lambda=0;
    p=3;
    sigma=0.5;
    beta=1*ones(p,1);
    X=randn(n,p)+1;
    epsilon=randn(n,1);
    
    y=X*beta+sigma*epsilon;
    y=normYJ(y,1,lambda,'inverse',true);
    
    % Create left censored observations
    qq=quantile(y,0.3);
    y(y<=qq)=qq;
    
    % Create right censored observations
     qq=quantile(y,0.9);
    y(y>=qq)=qq;
    
    % Specify the values of the transformation parameter
    la=[-0.5 0 0.5];
    left=min(y);
    right=max(y);
    % Call FSRfanCens and produce the fanplot
    [out]=FSRfanCens(y,X,'left',left,'right',right,'la',la,'init',round(n*0.5));
%}

%{

%}

%{
%}

%{

%}

%{
  
%}

%{

%}




%% Beginning of code

% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = aux.chkinputR(y,X,nnargin,vvarargin);

%% User options

% If the number of all possible subsets is <1000, the default is to extract
% all subsets, otherwise just 1000.
ncomb=bc(n,p);
nsamp=min(1000,ncomb);

if n<40
    init=round(n/2);
else
    init=(floor(0.5*(n+p+1)));
end
family='YJ';
plo=1;
lms=1;
conflev=0.99;
msg=true;
tag='pl_fan';
la=[-1 -0.5 0 0.5 1];
balancedSearch=true;

nocheck=false;
lwd=2;
lwdenv=1;
FontSize=12;
xlimx=[];
ylimy=[];
SizeAxesNum=10;
labx='Subset size m';
laby='Signed sqrt LR test statistic';
intercept=true;
titl='Fan plot';
left=0;
right=Inf;

if coder.target('MATLAB')

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)

        options=struct('la',la,...
            'nsamp',nsamp,'lms',lms,'plots',plo,'init',init,'conflev',conflev,...
            'titl',titl,'labx',labx,...
            'laby',laby,'xlimx',xlimx,'ylimy',ylimy,'lwd',lwd,...
            'lwdenv',lwdenv,'FontSize',FontSize,'SizeAxesNum',SizeAxesNum,...
            'tag',tag,'intercept',intercept,'msg',msg,'nocheck',nocheck,'family',family,...
            'left',left,'right',right,'balancedSearch',balancedSearch);

        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSRfanCens:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        aux.chkoptions(options,UserOptions)
    end
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    lms=options.lms;
    nsamp=options.nsamp;
    msg=options.msg;
    nocheck=options.nocheck;
    family=options.family;
    la=options.la;
    init=options.init;
    left=options.left;
    right=options.right;
    balancedSearch=options.balancedSearch;

    if coder.target('MATLAB')
        plo=options.plots;
        tag=options.tag;
        lwd=options.lwd;
        lwdenv=options.lwdenv;
        conflev=options.conflev;
        labx=options.labx;
        laby=options.laby;
        titl=options.titl;
        ylimy=options.ylimy;
        xlimx=options.xlimx;
    end
end


if strcmp(family,'BoxCox')
    if left<=0
        error('FSDA:FSRfanCens:WrongFamily','BoxCox has been chosen but left truncation point is <=0J')
    end
    BoxCox=1;
elseif strcmp(family,'YJ')
    BoxCox=0;
else
    BoxCox=NaN;
    if coder.target('MATLAB')
        warning('FSDA:FSRfanCens:WrongFamily','Transformation family which has been chosen is not supported')
        error('FSDA:FSRfanCens:WrongFamily','Supported values are BoxCox or YJ')
    end
end

% Specify where to send the output of the current procedure if options plot
% =1
if coder.target('MATLAB')
    if plo==1
        h1=findobj('-depth',1,'tag',tag);
        if (~isempty(h1))
            clf(h1);
            figure(h1)
            axes;
        else
            figure;
            % include specified tag in the current plot
            set(gcf,'tag',tag);
        end
    end
end

if  init <p+1
    fprintf(['Attention : init should be larger than p+1. \n',...
        'It is set to p+2.']);
    init=p+2;
end


warnrank=warning('query','MATLAB:rankDeficientMatrix');
warnsing=warning('query','MATLAB:singularMatrix');
warnnear=warning('query','MATLAB:nearlySingularMatrix');
% Set them to off inside this function at the end of the file they will be
% restored to previous values
warning('off','MATLAB:rankDeficientMatrix');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

%% Start of the forward search
zeron1=false(n,1);

% Initialization of the n x 1 Boolean vector which contains a true in
% correspondence of the units belonging to subset in each step
bsbT=zeron1;

seq=(1:n)';

%  Unlai is a Matrix whose 2:11th col contains the unit(s) just included.
Unlai = cat(2 , (init+1:n)' , NaN(n-init,10));

% Un = cell which will contain the matrices Unlai for each value of lambda
lla=length(la);
Un=cell(lla,1);
if coder.target('MATLAB')
    % For code generation, before you use a cell array element, you must assign
    % a value to it. When you use cell to create a variable-size cell array,
    % for example, cell(1,n), MATLAB® assigns an empty matrix to each element.
    % However, for code generation, the elements are unassigned. For code
    % generation, after you use cell to create a variable-size cell array, you
    % must assign all elements of the cell array before any use of the cell
    % array.
    for i=1:lla
        Un{i,1}=Unlai;
    end
end

% Initialize matrix which will contain the signed sqrt LR test
Sco=[((init):n)'  NaN(n-init+1,lla)];
Exflag=Sco;
laMLE=Sco;

% The second column of matrix r will contain the OLS residuals at each step
% of the forward search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100 = 1000*(1:1:ceil(n/1000));
seq100(seq100>n)=[];
seq100boo=false(n,1);
seq100boo(seq100)=true;


binit=zeros(p,lla);

[nsampArray] = subsets(nsamp,n,p);

if balancedSearch==true

    yleft=seq(y==left);
    yright=seq(y==right);
    inleft=Inf(length(yleft),1);
    inright=Inf(length(yright),1);
    propleft=length(yleft)/n;
    propright=length(yright)/n;
end


% loop over the values of \lambda
for i=1:lla

    if BoxCox==1
        z=normBoxCox(y,1,la(i),'Jacobian',false,'inverse',false);
    else
        z=normYJ(y,1,la(i),'Jacobian',false,'inverse',false);
    end

    % Find initial subset to initialize the search using as y transformed
    % vector z
    [out]=LXS(z,X,'lms',lms,'nsamp',nsampArray,'nocheck',true,'msg',msg,...
        'nomes',0,'bonflevoutX','','conflev',0.99,'rew',0,'yxsave',false,'intercept',true);
    [~,indres]=sort(abs(out.residuals));
    bsb=indres(1:(init-3));

    bsbT(bsb)=true;

    %ini0 = initial value for forward search loop
    ini0=length(bsb);


    yb=y(bsb);
    Xb=X(bsb,:);

    if nocheck==false && (rank(Xb)~=p)
        if coder.target('MATLAB')
            warning('FSRfanCens:Wrngbsb','Initial subset does not form full rank matrix');
        else
            disp('Initial subset does not form full rank matrix')
        end
        % FS loop will not be performed
    else
        for mm=ini0:n
            % if n>1000 show every 100 steps the fwd search index
            if  msg==true && seq100boo(mm) == true
                disp(['m=' int2str(mm)]);
            end

            % Compute and store the value of signed sqrt lik ratio test
            outSC=regressCensTra(yb,Xb,"dispresults", ...
                false,"left",left,"right",right,'nocheck',true,'la0',la(i));

            if (mm>=init)

                % Store value of the signed sqrt LR test
                Sco(mm-init+1,i+1)=outSC.signLR;
                % Store information about convergence
                Exflag(mm-init+1,i+1)=outSC.Exflag;
                laMLE(mm-init+1,i+1)=outSC.Beta(end,1);
            end

            b=outSC.beta0(1:end-2,1);
            jacla0=outSC.jacla0;

            % e= (n x 1) vector of residuals for all units using b estimated
            % using subset and transformed response (normalized with the
            % Jacobian based on uncensored observations belonging to subset)
            e=z*jacla0-X*b;

            % r_i =e_i^2
            r(:,2)=e.^2;

            if mm<n

                % store units forming old subset in vector oldbsbT
                oldbsbT=bsbT;

                if balancedSearch==true
                    % Make sure that the proportion of truncated observations in
                    % the subset is as much as possible equal to the original one.
                    nleftsubset=round(propleft*(mm+1));
                    nrightsubset=round(propright*(mm+1));

                    yinleft=inleft;
                    [~,yleftind]=sort(r(yleft,2));
                    yinleft(yleftind(1:nleftsubset))=0;

                    yinright=inright;
                    [~,yrightind]=sort(r(yright,2));
                    yinright(yrightind(1:nrightsubset))=0;

                    r(yleft,2)=yinleft;
                    r(yright,2)=yinright;
                end


                % order the r_i and include the smallest among the units
                % forming the group of potential outliers
                [~,ord]=sort(r(:,2));


                % bsb= units forming the new subset
                bsb=ord(1:(mm+1),1);

                bsbT=zeron1;
                bsbT(bsb)=true;


                Xb=X(bsb,:);  % subset of X
                yb=y(bsb);    % subset of y

                if mm>=init

                    % unit = vector containing units which just entered subset;
                    % unit=setdiff(bsb,oldbsbT);
                    % new instruction to find unit
                    unit=seq(bsbT & ~oldbsbT);

                    if length(unit)<=10
                        Unlai(mm-init+1,2:(length(unit)+1))=unit;
                    else
                        if msg==true
                            disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        end
                        Unlai(mm-init+1,2:end)=unit(1:10);
                    end
                end
            end

        end
    end  % rank check
    % Store in cell Un matrix Unlai
    Un{i}=Unlai;
end

%% Structure returned by function FSRfan
Sco=real(Sco);
out=struct;

out.Score=Sco;
out.Exflag=Exflag;
out.laMLE=laMLE;
out.la=la;
out.bs=binit;
out.Un=Un;
out.y=y;
out.X=X;


out.class='FSRfanCens';

% Restore the previous state of the warnings
warning(warnrank.state,'MATLAB:rankDeficientMatrix');
warning(warnsing.state,'MATLAB:singularMatrix');
warning(warnnear.state,'MATLAB:nearlySingularMatrix');


if coder.target('MATLAB')

    if plo==1

        % Specify the line type for the units inside vector units
        slin={'-';'--';':';'-.'};
        slin=repmat(slin,ceil(lla/4),1);

        % Specify the color for the trajectories
        ColorOrd=[{[0 0 1]}; {[1 0 0]}; {[1 0 1]}; {[1 1 0]}; {[0 1 0]}; {[0 1 1]}];
        ColorOrd=repmat(ColorOrd,4,1);


        % plot the lines associated with the score test lwd = line width of the
        % trajectories which contain the score test
        plot1=plot(Sco(:,1),Sco(:,2:end),'LineWidth',lwd);

        set(plot1,{'Color'}, ColorOrd(1:lla,:));

        set(plot1,{'LineStyle'},slin(1:lla));


        rangeaxis=axis;
        quant = sqrt(chi2inv(conflev,1));
        numconflev=length(conflev);
        V=repmat([rangeaxis(1);rangeaxis(2)],1,2*numconflev);
        QUANT=[[quant;quant],[ -quant;-quant]];
        line(V, QUANT,'LineWidth',lwdenv,'color','r','LineWidth',lwdenv);



        if size(la,2)>1
            la=la';
        end
        text(n*ones(lla,1),Sco(end,2:end)',num2str(la));

        % set the x and y axis

        if ~isempty(xlimx)
            xlim(xlimx);
        end

        if isempty(ylimy)

            % Use default limits for y axis
            ylim1=max(-20,min([min(Sco(:,2:end)) -quant]));
            ylim2=min(20,max([max(Sco(:,2:end)) quant]));
            ylim([ylim1 ylim2]);
        else

            % Use limits specified by the user
            ylim(ylimy);
        end


        % Main title of the plot and labels for the axes

        title(titl);

        % FontSize = font size of the axes labels

        % Add to the plot the labels for values of la. Add the horizontal lines
        % representing asymptotic confidence bands
        xlabel(labx,'Fontsize',FontSize);
        ylabel(laby,'Fontsize',FontSize);

        % FontSizeAxesNum = font size for the axes numbers
        set(gca,'FontSize',SizeAxesNum)
        box on

    end
end

end
%FScategory:REG-Transformations
