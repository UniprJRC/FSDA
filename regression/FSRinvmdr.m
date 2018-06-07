function [MDRinv] = FSRinvmdr(mdr,p,varargin)
%FSRinvmdr converts values of minimum deletion residual into confidence levels
%
%<a href="matlab: docsearchFS('FSRinvmdr')">Link to the help function</a>
%
%  Required input arguments:
%
%    mdr : Minimum deletion residuals. Matrix. n-m0 x 2 matrix containing: 
%          1st col = fwd search index; 
%          2nd col = minimum deletion residual . 
%    p : Number of explanatory variables. Scalar. Number of explanatory variables of the underlying dataset
%           (including the intercept if present)
%
%  Optional input arguments:
%
%        n:     size of the sample. Scalar.
%               If it is not specified
%               it is set equal to mdr(end,1)+1
%               Example - 'n',10
%               Data Types - double
%   plots :  Plot on the screen. Scalar or structure.
%               It specify whether it is necessary to
%               plot in normal coordinates the value of mdr
%               If plots = 1, a plot which shows the
%               confidence level of mdr in each step is shown on the
%               screen.
%               Remark. three horizontal lines associated respectively with
%               values  0.01 0.5 and 0.99  are added to the plot
%               If plots is a structure the user can specify the following options
%                   conflev = vector containing horizontal lines associated
%                       with confidence levels
%                   conflevlab = scalar if it is equal 1 labels associated with
%                       horizontal lines are shown on the screen
%                   xlim = minimum and maximum on the x axis
%                   ylim = minimum and maximum on the y axis
%                   LineWidth = Line width of the trajectory of mdr in
%                   normal coordinates
%                   LineStyle = Line style of the
%                   trajectory of mle of transformation parameters
%                   LineWidthEnv = Line width of the horizontal lines
%                   Tag = tag of the plot (default is pl_mdrinv)
%                   FontSize = font size of the text labels which identify
%                   the trajectories
%                 Example - 'plots',1 
%                 Data Types - double
%  Output:
%
%   MDRinv:     Confidence levels. Matrix. Matrix with n-m0 rows (same rows
%               of input matrix mdr) and 3 columns: 
%               1st col = fwd search index from m0 to n-1. 
%               2nd col = confidence level of each value of mdr. 
%               3rd col = confidence level in normal coordinates: 50% conf
%               level becomes norminv(0.50)=0; 99% conf level becomes norminv(0.99)=2.33. 
%
%
% See also FSRenvmdr, LXS.m, FSREDA.m
%
% References:
%
% Atkinson, A.C. and Riani, M. (2006), Distribution theory and
% simulations for tests of outliers in regression, "Journal of
% Computational and Graphical Statistics", Vol. 15, pp. 460-476.
% Riani, M. and Atkinson, A.C. (2007), Fast calibrations of the forward
% search for testing multiple outliers in regression, "Advances in Data
% Analysis and Classification", Vol. 1, pp. 123-141.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRinvmdr')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{  
    % FSRinvmdr with all default options.
    % Example of finding the confidence level of MDRenv, where MDRenv is
    % the matrix of 99 per cent confidence envelopes based on 1000
    % observations and 5 explanatory variables.
    % MDinv is a matrix which in the second column contains
    % all values equal to 0.99
    p=5;
    MDRenv=FSRenvmdr(1000,p,'prob',0.99);
    MDRinv=FSRinvmdr(MDRenv,p);      
%}

%{
    %% FSRinvmdr with optional arguments.
    % Example of finding confidence level of mdr for untransformed wool
    % data.
    % In the example, the values of mdr are plotted and then transformed 
    % into observed confidence levels.
    % The output is plotted in normal coordinates.

    load('wool.txt','wool');
    y=wool(:,4);
    X=wool(:,1:3);
    % The line below shows the plot of mdr
    [out]=FSR(y,X,'nsamp',0,'plots',1);
    MDRinv=FSRinvmdr(out.mdr,size(X,2)+1,'plots',1);
%}

%{
    % Resuperimposing envelopes.
    % Comparison of resuperimposing envelopes using mdr coordinates and normal
    % coordinates again on wool data.
    load('wool.txt','wool');
    y=wool(:,4);
    X=wool(:,1:3);
    % The line below shows the plot of mdr
    [out]=FSR(y,X,'nsamp',0,'plots',2);

    n0=16:19;
    quantplo=[0.01 0.5 0.99 0.999 0.9999 0.99999];
    ninv=norminv(quantplo);
    lwdenv=2;
    ij=0;
    supn0=max(n0);

    for jn0=n0;
        ij=ij+1;
        MDRinv = FSRinvmdr(out.mdr,4,'n',jn0);
        % Resuperimposed envelope in normal coordinates
        subplot(2,2,ij)
        plot(MDRinv(:,1),norminv(MDRinv(:,2)),'LineWidth',2)
        xlim([0 supn0])
        v=axis;
        line(v(1:2)',[ninv;ninv],'color','g','LineWidth',lwdenv,'LineStyle','--','Tag','env');
        text(v(1)*ones(length(quantplo),1),ninv',strcat(num2str(100*quantplo'),'%'));
        % line(MDRinv(:,1),norminv(MDRinv(:,2)),'LineWidth',2)
        title(['Resuperimposed envelope n=' num2str(jn0)]);
    end
%}

%{
    % Example with Hospital Data.
    % Comparison of resuperimposing envelopes using mdr coordinates and normal
    % coordinates at particular steps.
    load('hospital.txt');
    y=hospital(:,5);
    X=hospital(:,1:4);
    % exploratory analysis through the yXplot
    out=FSR(y,X,'nsamp',20000,'plots',2,'lms',0);

    n0=[54 58 62 63];
    quantplo=[0.01 0.5 0.99 0.999 0.9999 0.99999];
    ninv=norminv(quantplo);
    lwdenv=2;
    supn0=max(n0);

    figure;
    ij=0;
    for jn0=n0;
        ij=ij+1;
        [MDRinv] = FSRinvmdr(out.mdr,5,'n',jn0);
        % Plot for each step of the fwd search the values of mdr translated in
        % Plot for each step of the fwd search the values of mdr translated in
        % terms of normal quantiles
        subplot(2,2,ij)
        plot(MDRinv(:,1),norminv(MDRinv(:,2)),'LineWidth',2)
        xlim([0 supn0])
        v=axis;
        line(v(1:2)',[ninv;ninv],'color','g','LineWidth',lwdenv,'LineStyle','--','Tag','env');
        text(v(1)*ones(length(quantplo),1),ninv',strcat(num2str(100*quantplo'),'%'));
        line(MDRinv(:,1),norminv(MDRinv(:,2)),'LineWidth',2)
        title(['Resuperimposed envelope n=' num2str(jn0)]);
    end

%}



%% Input parameters checks
n=mdr(end,1)+1;
options=struct('n',n,'plots','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRinvmdr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

if nargin>2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

n=options.n;

if ~isscalar(n) || isempty(n) || isnan(n)
error('FSDA:FSRinvmdr:missingInputs','n must be scalar non empty and non missing!!');end

if ~isscalar(p) || isempty(n) || isnan(p)
    error('FSDA:FSRinvmdr:missingInputs','p must be scalar non empty and non missing!!!');
end


%% Find confidence level of each value of mdr

% mm = column vector which contains fwd search index.
mm =mdr(:,1);

% Compute variance of the truncated normal distribution.
% mm/n is the percentage of observations inside subset.
a=norminv(0.5*(1+mm/n));
corr=1-2*(n./mm).*a.*normpdf(a);

qq=-mm-1+(mm+1)./(2*tcdf(mdr(:,2).*sqrt(corr),mm-p)-1);
qq=qq./(n-mm);
% mdr= 1-confidence level
mdr=fcdf(qq,2*(n-mm),2*(mm+1));
% mdrt = confidence level
mdrt=1-mdr;

% Confidence level in normal coordinates
% One may wonder why mmdncoo cannot be computed as follows
% mmdncoo=norminv(mmt);
% The reason is that while norminv(1-\epsilon)=inf but
% -norminv(\epsilon) is exact where \epsilon is a very very small number
MDRncoord=-norminv(mdr);


MDRinv=[mm mdrt MDRncoord];


%% Plotting part

if ~isempty(options.plots)
   
    plots=options.plots;
    
    if isstruct(plots)
        
        fplots=fieldnames(plots);
        
        d=find(strcmp('xlim',fplots));
        if d>0
            xlimx=plots.xlim;
        else
            xlimx='';
        end
        
        
        d=find(strcmp('LineWidth',fplots));
        if d>0
            LineWidth=plots.LineWidth;
        else
            LineWidth=2;
        end
        
        % LineWidthEnv = line width of the horizontal lines associated
        % with required confidence envelopes
        % the plot of monitoring of MLE of transformation parameters or the
        % horizontal lines associated with the asymptotic confidence levels
        % for the likelihood ratio test
        d=find(strcmp('LineWidthEnv',fplots));
        if d>0
            LineWidthEnv=plots.LineWidthEnv;
        else
            LineWidthEnv=1;
        end
        
        % Specify the line type for the trajectory of mdr
        % in normal coordinate
        d=find(strcmp('LineStyle',fplots));
        if d>0
            LineStyle=plots.LineStyle;
        else
            LineStyle={'-'};
        end
        
        % Horizontal lines associated with confidence level
        d=find(strcmp('conflev',fplots));
        if d>0
            conflev=plots.conflev;
        else
            conflev=[0.01 0.5 0.99];
        end
        
        d=find(strcmp('ylim',fplots));
        if d>0
            ylimy=plots.ylim;
        else
            ylimy=[min([MDRncoord;norminv(min(conflev))]) max([MDRncoord;norminv(max(conflev))])];
        end
        
        d=find(strcmp('conflevlab',fplots));
        if d>0
            conflevlab=plots.conflevlab;
        else
            conflevlab=1;
        end
        
        d=find(strcmp('Tag',fplots));
        if d>0
            tag=plots.Tag;
        else
            tag='pl_mdrinv';
        end
        
        % Font size for the number of the axes and the titles of the
        % axes
        d=find(strcmp('FontSize',fplots));
        if d>0
            FontSize=plots.FontSize;
        else
            FontSize=12;
        end
        
        % Font size for text messages associated with the
        % labels of the horizontal lines
        d=find(strcmp('FontSizeLab',fplots));
        if d>0
            FontSizeLab=plots.FontSizeLab;
        else
            FontSizeLab=12;
        end
        
    else
        
        xlimx=[MDRinv(1,1) MDRinv(end,1)+1];
        ylimy=[min([MDRncoord;norminv(0.009)]) max([MDRncoord;norminv(0.991)])];
        LineWidth=2;
        LineWidthEnv=1;
        tag='pl_mdrinv';
        FontSize=12;
        FontSizeLab=12;
        LineStyle={'-'};
        conflev=[0.01 0.5 0.99];
        conflevlab=1;
    end
    
    % Specify where to send the output of the monitoring of mdr
    % in normal coordinates
    hmle=findobj('-depth',1,'tag',tag);
    if (~isempty(hmle))
        clf(hmle);
        figure(hmle)
        axes;
    else
        figure;
        set(gcf,'Name','MDR in normal coordinates');
    end
    
    plot1=plot(MDRinv(:,1),MDRncoord,'LineWidth',LineWidth);
    
    % Specify the line style of the trajectory of mdr in normal coordinates
    set(plot1,{'LineStyle'},LineStyle);
    
    
    set(gcf,'Tag',tag)
    
    if ~isempty(xlimx)
        xlim(xlimx);
    end
    
    if ~isempty(ylimy)
        ylim(ylimy);
    end
    
    % Add horizontal lines associated with confidenve levels
    v=axis;
    ninv=norminv(conflev);
    for i=1:length(conflev)
        if conflev(i)==0.5
            col=[1 0.69 0.39];
        elseif conflev(i)>=0.01 &&  conflev(i)<=0.99
            col='b';
        else
            col='r';
        end
        line(v(1:2)',[ninv(i);ninv(i)],'color',col,'LineWidth',LineWidthEnv,'LineStyle','--','Tag','env');
    end
    
    % add text label associated to horizontal confidence levels
    if conflevlab==1
        text(v(1)*ones(length(ninv),1),ninv+0.2,strcat(num2str(100*conflev'),'%'),...
            'FontSize',FontSizeLab,'HorizontalAlignment','Left');
    end
    
    xlabel('Subset size m','FontSize',FontSize);
    ylabel('MDR in normal coordinates','FontSize',FontSize);
    set(gca,'FontSize',FontSize)
    
end

end

%FScategory:REG-Regression

