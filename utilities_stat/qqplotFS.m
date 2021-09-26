function Y=qqplotFS(res,varargin)
%qqplot of studentized residuals with envelopes
%
%<a href="matlab: docsearchFS('qqplotFS')">Link to the help function</a>
%
% Displays a quantile-quantile plot of the quantiles of the sample
% studentized residuals versus the theoretical quantile values from a
% normal distribution. If the distribution of residuals is normal, then the
% data plot appears linear. A confidence level is added to the band.
%
%
%  Required input arguments:
%
%    res: vector containing studentizd residuals. Vector.
%          Vector with n elements containing the studentized residuals from
%          a regression model
%
%  Optional input arguments:
%
%     conflev :  Confidence level which is
%               used to compute confidence bands of studentized residuals. Scalar
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.90
%                 Example - 'conflev',0.99
%                 Data Types - double
%
%    intercept :  Indicator for constant term. true (default) | false. 
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%     nsimul :  number of simulations to compute the envelopes. Scalar. The
%               default value is 1000.
%               Example - 'nsimul',300
%               Data Types - double
%
%    X :        Predictor variable. Matrix. Data matrix of explanatory
%               variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables. This is the matrix which has
%               been used to produce stdentized residuals.
%               If this optional argument is
%               missing we take as matrix X the column of ones.
%                   Example - 'X',randn(n,3)
%                   Data Types - double
%
%       plots : Plot on the screen. Scalar.
%               If plots = 1, a plot which shows the
%               robust the qqplot of residuals with envelopes is shown on the
%               screen. The confidence level which is used to draw the
%               horizontal lines associated with the bands for the
%               residuals is specified in input option conflev. If
%               conflev is missing a nominal 0.90 confidence interval will
%               be used.
%                 Example - 'plots',1
%                 Data Types - double
%
%               h : the axis handle of a figure where to send the qqplot.
%                   This can be used to host the qqplot in a subplot of a
%                   complex figure formed by different panels (for example a panel
%                   with qqplot from a classical ols estimator and another
%                   with qqplot from a robust regression).
%                   Example -'h',h1 where h1=subplot(2,1,1)
%                   Data Types - Axes object (supplied as a scalar)
%
%       tag     :   handle of the plot which is about to be created.
%                   Character. The default is to use tag 'pl_qq'. Notice
%                   that if the program finds a plot which has a tag equal
%                   to the one specified by the user, then the output of
%                   the new plot overwrites the existing one in the same
%                   window else a new window is created Example -
%                   'tag','mytag' Data Types - char
%
% Output:
%
%       Y:  matrix with raws data on which the plot is based. n-by-3 matrix.
%            1st col = standard normal quantiles.
%            2nd col = quantiles of input sample of studentized residual.
%            3rd col = lower confidence band of quantiles of studentized residuals.
%            4th col = upper confidence band of quantiles of studentized residuals.
%
% See also: qqplot, fitlm
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('qqplotFS')">Link to the help function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit

% Examples:

%{
    %% qqplot for the multiple regression data.
    % This is an example of the use of options X and plots
    load('multiple_regression.txt');
    y=multiple_regression(:,4);
    X=multiple_regression(:,1:3);

    outLM=fitlm(X,y,'exclude','');
    res=outLM.Residuals{:,3};

    qqplotFS(res,'X',X,'plots',1);
    title('qqplot of stud. res.')
    % No outlier appears
%}

%{
    %% qqplot with envelopes for the Wool data.
    % Compare the results using untransformed and transformed data.
    % This is an example of the use of option h
    XX=load('wool.txt');
    y=(XX(:,end));
    lny=log(y);
    X=XX(:,1:end-1);
    outLM=fitlm(X,y,'exclude','');
    res=outLM.Residuals{:,3};

    outLMtra=fitlm(X,lny,'exclude','');
    restra=outLMtra.Residuals{:,3};

    h1=subplot(1,2,1);
    qqplotFS(res,'X',X,'plots',1,'h',h1);
    title('QQplot using untransformed data')
    h2=subplot(1,2,2);
    qqplotFS(restra,'X',X,'plots',1,'h',h2);
    title('QQplot using transformed data')
%}


%% Beginning of code

n=length(res);
X=ones(n,1);

options=struct('intercept',true,'X',X,'plots',0,'conflev',0.90,'nsimul',1000,...
    'tag','pl_qq','h','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:qqplotFS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


% Write in structure 'options' the options chosen by the user
if nargin > 1
    for j=1:2:length(varargin)
        options.(varargin{j})=varargin{j+1};
    end
end

X=options.X;
plo=options.plots;
conflev=options.conflev;
nsimul=options.nsimul;
intercept=options.intercept;
h=options.h;
tag=options.tag;
conflow=(1-conflev)/2;

quantinf=round(nsimul*conflow);
quantsup=round(nsimul*(1-conflow));

nsimul=options.nsimul;
if size(res,1)==1
    res = res';
end

sorres  =  sort(res);
xcoo = plotpos(sorres);
xcoo = norminv(xcoo);

% hdat = line(xcoo,sorres,'LineStyle','none','Marker','+');
Y=[xcoo sorres];

if intercept==true || intercept == 1
    X=[ones(n,1) X];
end

% Compute matrix A = Matrix I - H
[~, R] = qr(X,0);
E = X/R;
% Notice that:
% -E*E' = matrix -H = -X*inv(X'X)*X' computed through qr decomposition
A = -E*E';
seq=1:n;
linind=sub2ind(size(A),seq,seq);
A(linind)=1+A(linind);
onesminushii=sqrt(diag(A));


% each column of matrix epsi contains the n raw residuals for simulation j=1,
% ..., nsimul
% e=M\epsilon
epsi=randn(n,nsimul);
for j=1:1000
    epsi(:,j)=sort(A*epsi(:,j)./onesminushii);
end

ConfBand=zeros(n,2);
for i=1:n
    eo=sort(epsi(i,:));
    ConfBand(i,:)=[eo(quantinf) eo(quantsup)];
end

% Specify where to send the output of the current procedure if options plot
% =1
if plo==1
    
    % Create the figure that will host the qqplot
    
    h1=findobj('-depth',1,'tag',tag);
    if (~isempty(h1)) && isempty(h)
        clf(h1);
        hfig=figure('Name', 'qq plot', 'NumberTitle', 'off', 'Tag',tag);
        clf reset
    else
        hfig=figure('Name', 'qq plot', 'NumberTitle', 'off', 'Tag',tag);
        % include specified tag in the current plot
        % set(gcf,'tag',options.tag);
    end
    
    % hfig = figure('Name', 'qq plot', 'NumberTitle', 'off',...
    %     'Tag',tag);
    
    % Get figure's axis
    afig = axes('Parent',hfig);
    
    % Plot the resindexplot and add relevant labels
    plot(afig,Y(:,1),Y(:,2),'o')
    line(afig,Y(:,1),ConfBand(:,1),'Color','r')
    line(afig,Y(:,1),ConfBand(:,2),'Color','r')
    % xlabel('Standard normal quantiles')
    % ylabel('Studentized residuals')
    labx='Standard normal quantiles';
    laby='Studentized residuals';
    titl='';
    FontSize=12;
    SizeAxesNum=12;
    
    if ~isempty(h)
        % Eventually send the resindexplot into a different figure/subplot
        hfigh = get(h,'Parent');
        
        set(hfigh,'Name','QQ plot','NumberTitle','off');
        set(h,'Tag','qq_subplot');
        copyobj(allchild(afig),h);
        pause(0.0000001);
        delete(hfig);
        
        % Fix the y-axis
        set(h,'YLimMode', 'manual');
        % Add title and axis labels for the figure with subplots, and set their FontSize
        title(gca,titl);
        xlabel(gca,labx,'Fontsize',FontSize);
        ylabel(gca,laby,'Fontsize',FontSize);
        % Set the font size for the axes numbers
        set(gca,'FontSize',SizeAxesNum);
        
    else
        % If the resindexplot has not to be sent in a different figure/subplot
        % add the figure title and axis labels, and set their FontSize
        title(afig,titl);
        xlabel(afig,labx,'Fontsize',FontSize);
        ylabel(afig,laby,'Fontsize',FontSize);
        
        
    end
    
end

Y=[Y ConfBand];
end

function pp = plotpos(sx)
%PLOTPOS Compute plotting positions for a probability plot
%   PP = PLOTPOS(SX) compute the plotting positions for a probability
%   plot of the columns of SX (or for SX itself if it is a vector).
%   SX must be sorted before being passed into PLOTPOS.  The ith
%   value of SX has plotting position (i-0.5)/n, where n is
%   the number of rows of SX.  NaN values are removed before
%   computing the plotting positions.

[n, m] = size(sx);
if n == 1
    sx = sx';
    n = m;
    m = 1;
end

nvec = sum(~isnan(sx));
pp = repmat((1:n)', 1, m);
pp = (pp-.5) ./ repmat(nvec, n, 1);
pp(isnan(sx)) = NaN;
end

%FScategory:VIS-Reg
