function [out]=FSRaddt(y,X,varargin)
%FSRaddt produces t deletion tests for each explanatory variable.
%
%<a href="matlab: docsearchFS('FSRaddt')">Link to the help function</a>
%
% Required input arguments:
%
%    y: A vector with n elements that contains the response variable. y can
%       be both a row of column vector.
%    X: Data matrix of explanatory variables (also called 'regressors') of
%       dimension (n x p-1). Rows of X represent observations, and columns
%       represent variables.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       automatically be excluded from the computations.
%
% Optional input arguments:
%
%   intercept   :  Indicator for constant term. Scalar. 
%                       If 1, a model with constant term will be fitted
%                       (default), if 0, no constant term will be included.
%                        Example - 'intercept',1 
%                       Data Types - double
%           h   :      The number of observations that have determined the
%                       least trimmed squares estimator. Scalar.
%                       h is an integer greater or
%                       equal than [(n+size(X,2)+1)/2] but smaller then n
%                       Example - 'h',round(n*0,75) 
%                       Data Types - double
%       nsamp   : Number of subsamplse which will be extracted to find the
%                       robust estimator. Scalar.
%                       If nsamp=0 all subsets will be
%                       extracted. They will be (n choose p). Remark: if the
%                       number of all possible subset is <1000 the default is to
%                       extract all subsets otherwise just 1000.
%                       Example - 'nsamp',1000 
%                        Data Types - double
%       lms     :    Criterion to use to find the initlal subset to
%                       initialize the search. Scalar,  vector or structure.
%                       If lms=1 (default) Least Median of Squares is
%                       computed, else Least Trimmed of Squares is computed. else
%                       (default) no plot is produced
%                       Example - 'lms',1 
%                       Data Types - double
%       init    :       Search initialization. Scalar. 
%                       Scalar which specifies the initial subset size to start
%                       monitoring exceedances of minimum deletion residual, if
%                       init is not specified it will be set equal to:
%                       p+1, if the sample size is smaller than 40;
%                       min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%                       Example - 'init',100 starts monitoring from step m=100 
%                       Data Types - double
%       plots   :    Plot on the screen. Scalar.
%                       If plots=1 a plot with forward deletion
%                        t-statistics is produced
%                        Example - 'plots',1 
%                        Data Types - double
%        nameX  : Add variable labels in plot. Cell array of strings.
%                       Cell array of strings of length p containing the labels of
%                       the varibles of the regression dataset. If it is empty
%                       (default) the sequence X1, ..., Xp will be created
%                       automatically
%                       Example - 'nameX',{'NameVar1','NameVar2'} 
%                       Data Types - cell
%       lwdenv  : Line width for envelopes. Scalar.
%                       Line width for envelopes based on student T (default is 2)
%                        Example - 'lwdenv',1 
%                        Data Types - double
%        quant  :  Confidence quantiles for the envelopes. Vector.
%                       Confidence quantiles for the envelopes of deletion t
%                        stat. Default is [0.005 0.995] (i.e. a 99% pointwise
%                       confidence interval)
%                        Example - 'quant',[0.025 0.975]
%                        Data Types - double
%       lwdt       : Line width for deletion T stat. Scalar.
%                       (default is 2)
%                        Example - 'lwdt',1 
%                        Data Types - double
%       nocheck : Check input arguments. Scalar.
%                       If nocheck is equal to 1 no check is performed on
%                       matrix y and matrix X. Notice that y and X are left
%                       unchanged. In other words the additional column of ones
%                       for the intercept is not added. As default nocheck=0.
%                       Example - 'nocheck',1 
%                       Data Types - double
%       titl    :       a label for the title. Character.
%                       (default: '')
%                       Example - 'titl','Example' 
%                       Data Types - char
%       labx    :   a label for the x-axis. Character.
%                       (default: 'Subset size m')
%                       Example - 'labx','Subset' 
%                       Data Types - char
%       laby    :     a label for the y-axis. Character.
%                       (default: 'Deletion t statistics')
%                       Example - 'laby','statistics' 
%                       Data Types - char
%     FontSize:  the font size of the labels of
%                       the axes and of the labels inside the plot. Scalar.
%                       Default value is 12
%                       Example - 'FontSize',11 
%                       Data Types - double
% SizeAxesNum: size of the numbers of the axes. Scalar.
%                       Default value is 10
%                       Example - 'SizeAxesNum',11 
%                       Data Types - double
%          ylimy:    minimum and maximum of the y axis. Vector.
%                        Default value is '' (automatic scale)
%                       Example - 'ylimy',[0 1] 
%                       Data Types - double
%          xlimx:   minimum and maximum of the x axis. Vector.
%                       Default value is '' (automatic scale)
%                       Example - 'xlimy',[0 1] 
%                       Data Types - double
% Output:
%
%         out:   structure which contains the following fields
%
%
%  out.Tdel=    (n-init+1) x (p+1) matrix containing the monitoring of
%               deletion t stat in each step of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = deletion t stat for first explanatory variable
%               3rd col = deletion t stat for second explanatory variable
%               ...
%               (p+1)th col = deletion t stat for pth explanatory variable
%
%  out.S2del=   (n-init+1) x (p+1) matrix containing the monitoring of
%               deletion t stat in each step of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = deletion t stat for first explanatory variable
%               3rd col = deletion t stat for second explanatory variable
%               ...
%               (p+1)th col = deletion t stat for pth explanatory variable
%    out.Una=   cell of size p.
%               out.Una{i} (i=1, ..., p) is a (n-init) x 11 matrix which
%               contains the unit(s) included in the subset at each step in
%               the search which excludes the ith explanatory variable. 
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,:) for example contains
%               the unit included in step init+1 ... Un(end,2) contains the
%               units included in the final step of the search
%
% See also addt
%
% References: 
%
%   Atkinson, A.C. and Riani, M. (2002b). Forward search added variable t
%   tests and the effect of masked outliers on model selection. Biometrika,
%   Vol. 89, pp. 939-946.
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRaddt')">Link to the help function</a>
% Last modified 11-06-2016

% Examples

%{
    % FSRaddt with all default options.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    [out]=FSRaddt(y,X);
%}

%{
    %%FSRaddt with optional arguments.
    % We perform a variable selection on the 'famous' stack loss data using different transformation scales for the response.
    load('stack_loss');
    y=stack_loss.data(:,4);
    X=stack_loss.data(:,1:3);
    %We start with a fan plot based on first-order model and the five most common values of ? (Figure below).
    [out]=FSRfan(y,X,'plots',1);
    %The fan plot shows that the square root transformation, ?= 0.5, is supported by all the data, with the absolute value of the statistic always less than 1.5. The evidence for all other transformations depends on which observations have been deleted: the log transformation is rejected when some of the suspected outliers are introduced into the data although it is acceptable for all the data: ?= 1 is rejected as soon as any of the suspected outliers are present.
    %Given that the transformation for the response which is chosen depends on the number of units declared as outliers we perform a variable selection using the original scale, the square root and the log transformation.
    %Robust variable selection using original untransformed values of the response
    % Monitoring of deletion t stat in the original scale
    [out]=FSRaddt(y,X,'plots',1,'quant',[0.025 0.975]);
    %Robust variable selection using square root values
    % Monitoring of deletion t stat using transformed response based on the square root
    [out]=FSRaddt(y.^0.5,X,'plots',1,'quant',[0.025 0.975]);
    %Robust variable selection using log transformed values of the response
    % Monitoring of deletion t stat using log transformed values
    [out]=FSRaddt(log(y),X,'plots',1,'quant',[0.025 0.975]);
    %Conclusion: the forward analysis based on the deletion t statistics clearly reveals that variable X3, independently from the transformation which is chosen and the number of outliers which are declared, is NOT significant.
%}

%{
    % FSRaddt with optional arguments.
    % Example of use of FSRaddt with plot of deletion t with personalized line
    % width for the envelopes and personalized confidence interval.
    n=200;
    p=3;
    X=randn(n,p);
    y=randn(n,1);
    kk=9;
    y(1:kk)=y(1:kk)+6;
    X(1:kk,:)=X(1:kk,:)+3;
    [out]=FSRaddt(y,X,'plots',1,'quant',[0.025 0.975]);
%}

%{
    % FSRaddt with plots (transformed wool data).
    load('wool');
    y=log(wool.data(:,end));
    X=wool.data(:,1:end-1);
    [out]=FSRaddt(y,X,'plots',1);
%}

%{
    %% FSRaddt with labels for the columns of matrix X.
    % Line width equal to 3 for the curves representing envelopes; 
    % line width equal to 4 for the curves associated with deletion t stat.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    [out]=FSRaddt(y,X,'plots',1,'nameX',{'F1','F2','F3'},'lwdenv',3,'lwdt',4);
%}


%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

% Reduce the number of variables by 1 becuse we delete each variable in
% turn
p1=p-1;

%% User options

% If the number of all possible subsets is <1000 the default is to extract
% all subsets, otherwise just 1000.
ncomb=bc(n,p1);
nsampdef=min(1000,ncomb);

% The default value of h is floor(0.5*(n+p+1))
hdef=floor(0.5*(n+p+1));
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));  
end

options=struct('h',hdef,...
    'nsamp',nsampdef,'lms',1,'plots',0,...
    'init',init,'nameX','','lwdenv',2,'quant',[0.005 0.995],'lwdt',2,'xlimx','','ylimy','',...
    'titl','','labx','Subset size m','laby','Deletion t statistics','FontSize',12,'SizeAxesNum',10,'nocheck',0,'intercept',1);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRaddt:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
plo=options.plots;
nsamp=options.nsamp;

% intcolumn = the index of the first constant column found in X, or empty.
% Used here to check if X includes the constant term for the intercept.
intcolumn = find(max(X,[],1)-min(X,[],1) == 0, 1);

%vars: list of the variables
vars=1:p; vars(intcolumn)=[];
ini0=length(vars);


% sequence from 1 to n
seq=(1:n)';

% The second column of matrix R will contain the OLS residuals
% at each step of the forward search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100=100*(1:1:ceil(n/100));


%% Start of the forward search
%
% Matrix Tdel will contain the deletion T-stat for each explanatory variable in each step of the fwd
% search.
% 1st col = fwd search index
% 2nd col = added T stat for first explanatory variable
% ...
% ini0+1 col = added T stat for final explanatory variable
Tdel=[(init:n)' zeros(n-init+1,ini0)];


% Matrix S2del will contain the deletion \sigma^2 estimate when corresponding
% explanatory variable is deleted from the forward search
% search.
% 1st col = fwd search index
% 2nd col = estimate of \sigma^2 when first explanatory variable is deleted
% ...
% ini0+1 col = estimate of \sigma^2 when first explanatory variable is
% deleted
S2del=[(init:n)' zeros(n-init+1,ini0)];


% Una = cell which will contain the matrices Un for each deleted explanatory
% variable
Una=cell(ini0,1);

%  Un is a Matrix whose 2:11th col contains the unit(s) just
%  included.
Uni = cat(2 , (init+1:n)' , NaN(n-init,10));


j=0;
% A different forward search for each variable which is excluded
for i=vars

    j=j+1;

    % Initialize matrix Un
    Un=Uni;

    %w: the variable to test;
    %Xred: the remaining variables.
    w=X(:,i);
    Xred=X;
    Xred(:,i)=[];

    % Find initial subset to initialize the search (in the search which
    % excludes variable w
    [outLXS]=LXS(y,Xred,'lms',lms,'h',h,'nsamp',nsamp,'nocheck',1);
    bsb=outLXS.bs;

    Xb=Xred(bsb,:); % Subset of reduced X matrix (in the search which excludes variable w)
    yb=y(bsb); % Subset of y  (in the search which excludes variable w)
    wb=w(bsb); % Subset of vector w (excluded explanatory variable)

    % Initialize the fwd search (excluding variable w)
    for mm=ini0:n
        % if n>200 show every 100 steps the fwd search index
        if n>200
            if length(intersect(mm,seq100))==1
                disp(['m=' int2str(mm)]);
            end
        end

        b=Xb\yb;

        if (mm>=init)
            % compute added t test
            [outADDT]=addt(yb,Xb,wb,'intercept',0);
            % Store added tstat
            Tdel(mm-init+1,i)=outADDT.Tadd;
            % store added estimate of S2
            S2del(mm-init+1,i)=outADDT.S2add;
        end

        % e = vector of residual for all units using b estimated using subset
        e=y-Xred*b;

        % Second column of matrix r contains squared residuals
        r(:,2)=e.^2;

        if mm<n

            % store units forming old subset in vector oldbsb
            oldbsb=bsb;

            % order the r_i and select the smallest m+1
            ord=sortrows(r,2);

            % bsb= units forming the new  subset
            bsb=ord(1:(mm+1),1);

            Xb=Xred(bsb,:);  % subset of X
            yb=y(bsb);    % subset of y
            wb=w(bsb);    % subset of w


            % if mm>=init it is necessary to store the required quantities
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

    % Store in cell Una matrix Un
    Una{j}=Un;
end

out.Tdel=Tdel;
out.S2del=S2del;
out.Una=Una;

if plo==1

    lwdt=options.lwdt;
    plot1=plot(Tdel(:,1),Tdel(:,2:end),'LineWidth',lwdt);

    % Specify the line type for the units inside vector units
    slin={'-';'--';':';'-.'};
    slin=repmat(slin,ceil(p/4),1);
    set(plot1,{'LineStyle'},slin(1:p-1));


    % Main title of the plot and labels for the axes
    labx=options.labx;
    laby=options.laby;
    titl=options.titl;

    title(titl);

    % FontSize = font size of the axes labels
    FontSize =options.FontSize;

    % Add to the plot the labels for values of la
    % Add the horizontal lines representing asymptotic confidence bands
    xlabel(labx,'Fontsize',FontSize);
    ylabel(laby,'Fontsize',FontSize);

    % SizeAxesNum = font size for the axes numbers
    SizeAxesNum=options.SizeAxesNum;
    set(gca,'FontSize',SizeAxesNum)


    lwdenv=options.lwdenv;

    FontSize=12;
    kx=0;
    ky=0;
    % x coordinates where to put the messages about envelopes

    if ~isempty(options.xlimx)
        xlimx=options.xlimx;
        xlim(xlimx);

        % x coordinates where to put the messages about envelopes
        xcoord=max([xlimx(1) init]);
    else
        xcoord=init;
    end

    if ~isempty(options.ylimy)
        ylimy=options.ylimy;
        ylim(ylimy);
    end


    quant=options.quant;

    % The instruction 'repmat(quant,length(Tdel),1)' replicates vector quant length(Tdel) times
    % The instruction 'repmat(Tdel(:,1),1,length(quant))' replicates vector
    % Tdel(:,1) length(quant) times
    Tdelenv=tinv(repmat(quant,length(Tdel(:,1)),1),repmat(Tdel(:,1),1,length(quant))-p);

    for i=1:length(quant)

        % Superimpose chosen envelopes
        if quant(i)>=0.005 && quant(i) <=0.995
            line(Tdel(:,1),Tdelenv(:,i),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        else
            line(Tdel(:,1),tinv(quant(i),Tdel(:,1)-p),'LineWidth',lwdenv,'LineStyle','--','Color',[0  0 0],'tag','env');
        end

        [figx, figy] = dsxy2figxy(gca, xcoord,Tdelenv(Tdel(:,1)==xcoord,i));
        if isempty(figy) || figy<0
            figy=0;
        end
        if isempty(figx) || figx<0
            figx=0;
        end

        if isempty(figy) || figy>1
            figy=1;
        end
        if isempty(figx) || figx>1
            figx=1;
        end
        
        annotation(gcf,'textbox',[figx figy kx ky],'String',{[num2str(100*quant(i)) '%']},...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'EdgeColor','none',...
            'BackgroundColor','none',...
            'FitBoxToText','off',...
            'FontSize',FontSize);
    end

    if isempty(options.nameX)
        nameX=cellstr(num2str((vars-1)','X%d'));
    else
        nameX=options.nameX;
    end

    % Add labels at the end of the search
    text(n*ones(ini0,1),Tdel(end,2:end)',nameX,'FontSize',FontSize);

end

end
%FScategory:REG-ModelSelection