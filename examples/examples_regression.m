%% EXAMPLES OF ROBUST REGRESSION
% examples_regression shows a series of analysis of regression datasets
% Copyright 2008-2023.
% Written by FSDA team
%
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

%% FD (Forbes data) -- Forward EDA (Exploratory Data Analysis with FS)
% Interactive_example
clearvars;close all;
% scatterplot of data: one point looks outlying
load('forbes.txt');
y=forbes(:,2);
X=forbes(:,1);
plot(X,y,'o');
xlabel('Boiling point')
ylabel('100 x log(pressure)')
set(gcf,'Name', 'Plot of y versus X','NumberTitle', 'off')
% running the search
[out]=LXS(y,X);
[out]=FSReda(y,X,out.bs);
% Plot minimum deletion residual
mdrplot(out,'xlimx',[6 17],'ylimy',[0 13]);
% Now, some interactive brushing starting from the monitoring residuals
% plot. Once a set of trajectories is highlighted in the monitoring residual plot,
% the corresponding units are highlighted in the other plots
databrush=struct;
databrush.bivarfit='i1';
databrush.selectionmode='Rect'; % Rectangular selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
databrush.RemoveTool='on'; %remove yellow

cascade;
resfwdplot(out,'databrush',databrush);


%% FD (Forbes data) -- EDA and Analysis using S estimators
% Interactive_example
% scatterplot of data: one point looks outlying
clearvars;close all;
load('forbes.txt');
y=forbes(:,2);
X=forbes(:,1);
plot(X,y,'o');
xlabel('Boiling point')
ylabel('100 x log(pressure)')
set(gcf,'Name', 'Plot of y versus X','NumberTitle', 'off')
[out]=Sreg(y,X);

%Remark: if you want to use MMestimators, simply replace
% [out]=Sreg(y,X);
% with
% [out]=MMreg(y,X);

% Now, some interactive brushing starting from the index plot of residuals
% Once a set of trajectories is highlighted in the index plot of residuals
% the corresponding units are highlighted in the scatter plot
databrush=struct;databrush.selectionmode='Rect'; % Rectangular selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
databrush.RemoveTool='on'; %remove yellow

cascade;
out.X=X;
out.y=y;
resindexplot(out,'databrush',databrush);


%% MR: (Multiple regression data): Forward EDA with default options for resfwdplot
clearvars;close all;

load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);

% yX plot: called just using as input y and X, it reduces to a standard gplotmatrix
yXplot(y,X);

% LMS using 1000 subsamples
[out]=LXS(y,X,'nsamp',10000);
% Forward Search
[out]=FSReda(y,X,out.bs);
out1=out;
% Create scaled squared residuals
out1.RES=out.RES.^2;
% resfwdplot with all default options
resfwdplot(out1);
cascade;

%% MR: (Multiple regression data): S estimators with 2 values of breakdown point
clearvars;close all;

load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);

% yX plot: called just using as input y and X, it reduces to a standard gplotmatrix
yXplot(y,X);

conflev=[0.95 0.99];
% Note that the pattern of residuals changes completely
% Using bdp=0.5 the outliers are correctly found, on the other hand using
% bdp=0.25 the masking effect is clear
figure;
h1=subplot(2,1,1);
bdp=0.25;
[out]=Sreg(y,X,'nsamp',3000,'bdp',bdp);
resindexplot(out,'h',h1,'conflev',conflev);
ylabel(['Breakdown point =' num2str(bdp)])
h2=subplot(2,1,2);
bdp=0.5;
[out]=Sreg(y,X,'nsamp',3000,'bdp',bdp);
resindexplot(out,'h',h2,'conflev',conflev,'numlab',{6});
ylabel(['Breakdown point =' num2str(bdp)])
cascade;



%% MR: (Multiple regression data): MM estimators with 2 values of efficiency
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);

% yX plot: called just using as input y and X, it reduces to a standard gplotmatrix
yXplot(y,X);

% MMreg using two different level of efficiency
conflev=[0.95 0.99];
% Note that the pattern of residuals changes completely
% Using eff=0.90 the outliers are correctly found, on the other hand using
% eff=0.95 the masking effect is clear
figure;
h1=subplot(2,1,1);
eff=0.90;
[out]=MMreg(y,X,'Snsamp',3000,'eff',eff);
resindexplot(out,'h',h1,'conflev',conflev);
ylabel(['Eff.=' num2str(eff)])
h2=subplot(2,1,2);
eff=0.95;
[out]=MMreg(y,X,'Snsamp',3000,'eff',eff);
resindexplot(out,'h',h2,'conflev',conflev);
ylabel(['Eff.=' num2str(eff)])
cascade;

%% MR: Forward EDA with personalized options for resfwdplot
% Options for the "standard trajectories"
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);
% LMS using 10000 subsamples
[out]=LXS(y,X,'nsamp',10000);
% Forward Search
[out]=FSReda(y,X,out.bs);
out1=out;
% Create scaled squared residuals
out1.RES=out.RES.^2;

standard=struct;
standard.Color={'g'};
standard.LineWidth=1;
standard.xlim=[12 60];
standard.ylim=[0 30];
standard.SizeAxesNum=9;
standard.SizeAxesLab=10;

% Options for the trajectories in foreground
fground=struct;
fground.Color={'b'};
fground.LineWidth=2;
fground.funit=[21 47 38 31 30 9 43];

% Options for the trajectories in background
bgroud=struct;
bground.bthresh=4;
bground.bstyle='greysh';

resfwdplot(out1,'standard',standard,'fground',fground,'bground',bground);


%% MR (Multiple regression data): Forward EDA using personalized datatooltip
% Interactive_example
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);
% LMS using 1000 subsamples
[out]=LXS(y,X,'nsamp',10000);
% Forward Search
[out]=FSReda(y,X,out.bs);
out1=out;
% Create scaled squared residuals
out1.RES=out.RES.^2;

datatooltip=struct;
datatooltip.SubsetLinesColor=[1 0 0];
resfwdplot(out,'datatooltip',datatooltip)

%% MR (Multiple regression data): Forward EDA using persistent brushing
% Interactive_example
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);
% LMS using 1000 subsamples
[out]=LXS(y,X,'nsamp',10000);
% Forward Search
[out]=FSReda(y,X,out.bs);
out1=out;
% Create scaled squared residuals
out1.RES=out.RES.^2;

% plot minimum deletion residual with personalized options
mdrplot(out,'ylimy',[1 4.2],'xlimx',[10 60],'FontSize',14,'SizeAxesNum',14,'lwdenv',2);

% Persistent brushing on the plot of the scaled residuals. The plot is:
fground.flabstep='';        % without labels at steps 0 and n
fground.fthresh=3.5^2;      % threshold which defines the trajectories in foreground
fground.LineWidth=1.5;      % personalised linewidth for trajectories in foreground
fground.Color={'r'};        % personalised color (red lines) for trajectories in foreground

databrush=struct;
databrush.bivarfit='';
databrush.selectionmode='Rect'; % Rectangular selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
databrush.Pointer='hand'; % Hand cursor point while selecting
databrush.FlagSize='8'; % Size of the brushed points
databrush.RemoveTool='on'; % Remove yellow selection after finishing brushing
resfwdplot(out1,'fground',fground,'databrush',databrush);

%% MR: Forward EDA persistent brushing with other options
% Interactive_example
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);
% LMS using 1000 subsamples
[out]=LXS(y,X,'nsamp',10000);
% Forward Search
[out]=FSReda(y,X,out.bs);
out1=out;
% Create scaled squared residuals
out1.RES=out.RES.^2;

fground=struct;
fground.fthresh=3.1^2;      % threshold which defines the trajectories in foreground
fground.LineStyle={'--' '-.' ':'};    % different line styles for different foreground trajectories
fground.Color={'b';'g';'c';'m';'y';'k'};  % different colors for different foreground trajectories

databrush=struct;
databrush.bivarfit='';
databrush.selectionmode='Rect'; % Rectangular selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection

resfwdplot(out1,'fground',fground,'databrush',databrush);

%% MR: Forward EDA persistent brushing with labels at specific steps (e.g. 15 and 20).
% Interactive_example
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);
% LMS using 1000 subsamples
[out]=LXS(y,X,'nsamp',10000);
% Forward Search
[out]=FSReda(y,X,out.bs);
out1=out;
% Create scaled squared residuals
out1.RES=out.RES.^2;

fground.flabstep=[15 20];
databrush=struct;
databrush.bivarfit='';
databrush.selectionmode='Lasso'; % Lasso selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
resfwdplot(out1,'fground',fground,'databrush',databrush);


%% MR (Multiple regression data): S estimators with persistent brushing
% Interactive_example
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);
% S regression using 10000 subsamples
[out]=Sreg(y,X,'nsamp',10000,'yxsave',1);

% index plot of residuals with option databrush

databrush=struct;
databrush.bivarfit='';
databrush.selectionmode='Rect'; % Rectangular selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
databrush.Pointer='hand'; % Hand cursor point while selecting
databrush.FlagSize='8'; % Size of the brushed points
databrush.RemoveTool='on'; % Remove yellow selection after finishing brushing
resindexplot(out,'databrush',databrush);


%% MR: Added variable plot to show the importance of units 9 21 30 31 38 47
clearvars;close all;
load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);

% Computes in a new plot the added variable plot with and without the
% outliers
figure;
% Set Font Size for the title
fsiztitl=12;
% Set Font Size for the labels on the axes
SizeAxesNum=12;
outADD=addt(y,X(:,2:3),X(:,1),'plots',1,'units',[9 21 30 31 38 47]','textlab','y','FontSize',fsiztitl,'SizeAxesNum',SizeAxesNum);

% Computes in a new plot the added variable plot with and without unit 43
figure;
out43=addt(y,X(:,2:3),X(:,1),'plots',1,'units',43','textlab','y','FontSize',fsiztitl,'SizeAxesNum',SizeAxesNum);


%% HD: Hawkins data: Forward EDA
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);
yXplot(y,X);
[outLXS]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,outLXS.bs);
% Set graphical parameters for standard trajectories.
fsiztitl=16;
SizeAxesNum=15;

% Plot minimum deletion residual using personalized graphical options
mdrplot(out,'ylimy',[1 8],'xlimx',[25 128],'FontSize',fsiztitl,'SizeAxesNum',SizeAxesNum,'lwdenv',2,'quant',[0.01 0.5 0.99 0.9999]);

% Plot of resfwdplot with some options for trajectories in foreground
clear fground
fground.LineWidth=2;
fground.flabstep='';
fground.funit='';
resfwdplot(out,'fground',fground);

%% HD: example of resfwdplot with personalized options
% Explore different options for coloring the trajectories forming the
% monitoring residual plot
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

[outLXS]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,outLXS.bs);

standard=struct;
standard.xlim=[21 128];
standard.ylim=[-4.5 4.5];
standard.SizeAxesNum=12;
standard.SizeAxesLab=14;

fground=struct;
fground.Color={'b';'g';'c';'m';'y';'k'};
fground.flabstep=[28 80 ];
fground.fthresh=[-3 3];

bground=struct;
bground.bthresh=2;
bground.bstyle='greysh';
resfwdplot(out,'standard',standard,'fground',fground,'bground',bground);

%% HD: monitoring of S2 and beta coefficients
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

[outLXS]=LXS(y,X,'nsamp',10000);
[~,~,~,Bols,S2] = FSRmdr(y,X,outLXS.bs,'bsbmfullrank',0);

% Plot of the monitoring of S2 and R2
figure;
subplot(1,2,1)
plot(S2(:,1),S2(:,2))
xlabel('Subset size m');
ylabel('S2');
subplot(1,2,2)
plot(S2(:,1),S2(:,3))
xlabel('Subset size m');
ylabel('R2');

% Plots of the monitoring of "Estimates of beta coefficients"
figure;
for j=3:size(Bols,2)
    subplot(2,4,j-2)
    plot(Bols(:,1),Bols(:,j))
    xlabel('Subset size m');
    ylabel(['b' num2str(j-2)]);
end

%%  HD: resfwdplot persistent brushing
% Interactive_example
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

[outLXS]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,outLXS.bs);

databrush=struct;
databrush.bivarfit='';
databrush.selectionmode='Brush'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
% Note that if you wish to have the labels for the brushed units in the
% yXplot it is necessary to add the instruction in the line below
% databrush.labeladd='1'
resfwdplot(out,'databrush',databrush);

%%  HD: resfwdplot with persistent brushing and line fit resuperimposing
% Interactive_example
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

[outLXS]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,outLXS.bs);

databrush=struct;
databrush.bivarfit='2';
databrush.selectionmode='Brush'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
resfwdplot(out,'databrush',databrush);

%% HD: manual envelope resuperimposition in the plot of MDR
% (minimum deletion residual)
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

[outLXS]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,outLXS.bs);

mdrplot(out,'envm',86,'ylimy',[0 6],'tag','min86');
mdrplot(out,'envm',87,'ylimy',[0 6],'tag','min87');


%% HD: automatic procedure for outlier detection
% Use of FSR starting with 1000 subsamples
% focusing in the output plots to the interval 1-6 on the y axis and
% to steps 30-100.
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);
FSR(y,X,'nsamp',1000,'init',20,'ylim',[1 6],'xlim',[30 100]);

%% HD: analysis with S estimators (using two values of breakdown point)
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

% Simultaneous confidence level at 1%
conflev=1-0.01/length(y);
figure;
h1=subplot(2,1,1);
bdp=0.25;
[out]=Sreg(y,X,'nsamp',3000,'bdp',bdp);
resindexplot(out,'h',h1,'conflev',conflev);
ylabel(['Breakdown point =' num2str(bdp)])
h2=subplot(2,1,2);
bdp=0.5;
[out]=Sreg(y,X,'nsamp',3000,'bdp',bdp);
resindexplot(out,'h',h2,'conflev',conflev);
ylabel(['Breakdown point =' num2str(bdp)])

%% HD: analysis with MM estimators (using two values of efficiency)
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

% Simultaneous confidence level at 1%
conflev=1-0.01/length(y);

% MMreg using two different level of efficiency
% In this case the pattern of residuals reamins stable
figure;
h1=subplot(2,1,1);
eff=0.85;
[out]=MMreg(y,X,'Snsamp',3000,'eff',eff);
resindexplot(out,'h',h1,'conflev',conflev);
ylabel(['Efficiency =' num2str(eff)])
h2=subplot(2,1,2);
eff=0.95;
[out]=MMreg(y,X,'Snsamp',3000,'eff',eff);
resindexplot(out,'h',h2,'conflev',conflev);
ylabel(['Efficiency =' num2str(eff)])


%% HD: plot of minimum deletion residual with datatooltip
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);

[outLXS]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,outLXS.bs);
mdrplot(out,'ylimy',[1 8],'datatooltip',1);

%% HD: Example of dynamic brushing starting highlighting from the mdrplot
% Interactive_example
% Example of dynamic brushing starting highlighting from the mdrplot
% (that is the plot of minimum deletion residual)
clearvars;close all;
load('hawkins.txt');
y=hawkins(:,9);
X=hawkins(:,1:8);
[out]=LXS(y,X,'lms',0,'nsamp',10000);
[out]=FSReda(y,X,out.bs);
fground=struct;
fground.fthresh=2;
resfwdplot(out,'datatooltip','','fground',fground);

databrush=struct;
databrush.selectionmode='Brush'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of selected steps while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selecting
mdrplot(out,'ylimy',[1 8],'databrush',databrush);

%% WD (Wool data): forward EDA with untransformed data
% Interactive_example
clearvars;close all;
load('wool.txt','wool');
y=wool(:,4);
X=wool(:,1:3);
[out]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,out.bs);
mdrplot(out,'ylimy',[0.5 7]);
databrush=struct;
databrush.bivarfit='0';
databrush.selectionmode='Rect'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
resfwdplot(out,'databrush',databrush);


%% WD (Wool data): analysis using S estimators and brushing
% Interactive_example
clearvars;close all;
load('wool.txt','wool');
y=wool(:,4);
X=wool(:,1:3);

[out]=Sreg(y,X,'nsamp',3000,'yxsave',1);
resindexplot(out,'databrush',1);

%% WD: fan plot
% Log transformation is strongly suggested
clearvars;close all;
load('wool.txt','wool');
y=wool(:,4);
X=wool(:,1:3);
[outfan]=FSRfan(y,X,'plots',1,'init',5);

%% WD: automatic outlier detection procedure using logged observations
clearvars;close all;
load('wool.txt','wool');
y=wool(:,4);
X=wool(:,1:3);
[outFSR]=FSR(log(y),X);

%% WD: analysis using LTS
% Interactive_example
clearvars;close all;
load('wool.txt','wool');
y=wool(:,4);
X=wool(:,1:3);
[out]=LXS(log(y),X,'lms',0,'yxsave',1);
% Simultaneous confidence level
conflev = 1- 0.01/length(y);
resindexplot(out,'databrush',1,'conflev',conflev);


%% WD: analysis using MM estimators and residual brushing
% Interactive_example
clearvars;close all;
load('wool.txt','wool');
y=wool(:,4);
X=wool(:,1:3);
[out]=MMreg(log(y),X,'yxsave',1);
% Simultaneous confidence level
conflev = 1- 0.01/length(y);
resindexplot(out,'databrush',1,'conflev',conflev);

%% SD (Stack loss data): forward EDA
clearvars;close all;
load('stack_loss.txt');
y=stack_loss(:,4);
X=stack_loss(:,1:3);
nameX={'x1=air flow', 'x2=cooling water inlet temperature' 'x3=10 ? (acid concentration ?50)'};
namey={'y=Stack loss'};
yXplot(y,X,'nameX',nameX,'namey',namey);

[out]=LXS(y,X,'nsamp',0);
[out1]=FSReda(y,X,out.bs,'init',5);
% Create scaled squared residuals
%out1.RES=out1.RES.^2;
% Plot of minimum deletion residual
mdrplot(out1,'ylimy',[0.5 5],'xlimx',[5 21]);

%lab={'r1' 'r2' 'r3' 'r4' 'r5' 'r6' 'r7' 'r8' 'r9' 'r10' 'r11' 'r12' 'r13' 'r14' 'r15' 'r16' 'r17'};

% Plot of monitoring of residuals with datatooltip
resfwdplot(out1,'datatooltip',1)


%% SD: Brush starting from the monitoring residuals plot
% Interactive_example
clearvars;close all;
load('stack_loss.txt');
y=stack_loss(:,4);
X=stack_loss(:,1:3);
[out]=LXS(y,X,'nsamp',0);
[out1]=FSReda(y,X,out.bs,'init',5);

databrush=struct;
databrush.bivarfit='2';
databrush.selectionmode='Rect'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
resfwdplot(out1,'databrush',databrush);


%% SD: analysis using LTS
clearvars;close all;
load('stack_loss.txt');
y=stack_loss(:,4);
X=stack_loss(:,1:3);
[out]=LXS(y,X,'lms',0);
% Simultaneous confidence level
conflev = 1- 0.01/length(y);
resindexplot(out,'conflev',conflev);


%% SD: analysis using MM estimators
clearvars;close all;
load('stack_loss.txt');
y=stack_loss(:,4);
X=stack_loss(:,1:3);
[out]=MMreg(y,X,'eff',0.95);
% Simultaneous confidence level
conflev = 1- 0.01/length(y);
resindexplot(out,'conflev',conflev);



%% SD: Fan plot
clearvars;close all;
load('stack_loss.txt');
y=stack_loss(:,4);
X=stack_loss(:,1:3);
[out]=FSRfan(y,X,'plots',1);
fieldnames(out)

%% SD: Variable selection
% Monitoring of deletion t stat (untransformed values)
clearvars;close all;
load('stack_loss.txt');
y=stack_loss(:,4);
X=stack_loss(:,1:3);

figure;
FSRaddt(y,X,'plots',1,'quant',[0.025 0.975],'titl','Original scale');
% Monitoring of deletion t stat (square root scale)
figure;
FSRaddt(y.^0.5,X,'plots',1,'quant',[0.025 0.975],'titl','Square root scale');
% Monitoring of deletion t stat (log scale)
figure;
FSRaddt(log(y),X,'plots',1,'quant',[0.025 0.975],'titl','Log scale');
% Please notice that variable X3 is never significant in any reasonable
% scale

%% SP (Hospital data): fan plot
% Load logged hospital data
clearvars;close all;
load('hospitalFS.txt');
y=exp(hospitalFS(:,5));
X=hospitalFS(:,1:4);
% Fan plot
[outs]=FSRfan(y,X,'nsamp',10000,'plots',1);

%% SP: yX plot for the 2 hospitals
clearvars;close all;
load('hospitalFS.txt');
y=exp(hospitalFS(:,5));
X=hospitalFS(:,1:4);

y1=log(y);
n=length(y);
% exploratory analysis through the yXplot
group=ones(n,1);
group(55:108)=2;
unigroup=unique(group);
styp={'+';'o';'s';'d';'^';'v';'*';'x';'>';'<';'p';'h';'.'};

clr='brcmykgbrcmykgbrcmykg';
p=size(X,2);
p1=1:p;
namey=char('Logged survival time');
nameX={'x1=blood clotting score' 'x2=prognostic index' 'x3=enzyme test' 'x4=score for liver function'};

[H,AX,BigAx] = gplotmatrix(X,y1,group,clr(unigroup),char(styp{unigroup}),8,'on',[],nameX,namey);
set(H(:,:,1),'DisplayName','First hospital');
set(H(:,:,2),'DisplayName','Second hospital');


%% SP: Fwd search with EDA purposes
clearvars;close all;
load('hospitalFS.txt');
y1=hospitalFS(:,5);
X=hospitalFS(:,1:4);
[out]=LXS(y1,X,'nsamp',10000,'lms',0);
[out1]=FSReda(y1,X,out.bs);
% The plot of minimum deletion residual shows a peak out of the envelope in
% the central part of the search
mdrplot(out1,'quant',[0.01 0.5 0.99 0.9999 0.99999],'ylimy',[1 5],'lwdenv',2,'xlimx',[10 110]);

%% SP: persistent brushing starting from the plot of minimum deletion residual
% Interactive_example
clearvars;close all;
load('hospitalFS.txt');
y1=hospitalFS(:,5);
X=hospitalFS(:,1:4);
[out]=LXS(y1,X,'nsamp',10000,'lms',0);
[out1]=FSReda(y1,X,out.bs);

databrush=struct;
databrush.selectionmode='Rect'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.labeladd='1'; % Write labels of selected units in the yX matrix
databrush.RemoveLabels='off'; % Do not remove labels after selection
mdrplot(out1,'databrush',databrush);

%% SP: automatic outlier detection procedure based on FS
clearvars;close all;
load('hospitalFS.txt');
y1=hospitalFS(:,5);
X=hospitalFS(:,1:4);
[out]=LXS(y1,X,'nsamp',10000,'lms',0);
fieldnames(out)
% [out1]=FSReda(y1,X,out.bs);
[outFS]=FSR(y1,X,'init',20,'lms',0,'plots',2);
fieldnames(outFS)

%% SP: analysis using LMS and LTS
clearvars;close all;
load('hospitalFS.txt');
y1=hospitalFS(:,5);
X=hospitalFS(:,1:4);

% Define nominal confidence level
conflev=[0.99 1-0.01/length(y1)];
% Define number of subsets
nsamp=10000;
% Define the main title of the plots
titl='';
% Define upper xlim
cc=108;

% LMS with no rewighting
[outLMS]=LXS(y1,X,'nsamp',nsamp,'conflev',conflev(1));
h1=subplot(2,2,1);
laby='Scaled LMS residuals';
resindexplot(outLMS.residuals,'h',h1,'title',titl,'laby',laby,'numlab','','conflev',conflev,'xlimx',[0 cc])

% LTS with no rewighting
h2=subplot(2,2,2);
[outLTS]=LXS(y1,X,'nsamp',nsamp,'conflev',conflev(1),'lms',0);
laby='Scaled LTS residuals';
resindexplot(outLTS.residuals,'h',h2,'title',titl,'laby',laby,'numlab','','conflev',conflev,'xlimx',[0 cc]);

% LMS with reweighting
[outLMSr]=LXS(y1,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1);
h3=subplot(2,2,3);
laby='Scaled reweighted LMS residuals';
resindexplot(outLMSr.residuals,'h',h3,'title',titl,'laby',laby,'numlab','','conflev',conflev,'xlimx',[0 cc])

% LTS with reweighting
[outLTSr]=LXS(y1,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1,'lms',0);
h4=subplot(2,2,4);
laby='Scaled reweighted LTS residuals';
resindexplot(outLTSr.residuals,'h',h4,'title',titl,'laby',laby,'numlab','','conflev',conflev,'xlimx',[0 cc]);


%% SD: analysis using S estimators with 2 values of breakdown point
clearvars;close all;
load('hospitalFS.txt');
y=hospitalFS(:,5);
X=hospitalFS(:,1:4);
% Simulataneous confidence level
conflev = 1- 0.01/length(y);

% Sreg using two different level of breakdown point
% Using bdp=0.5 it is clear that the first 54 units have a pattern of residuals
% which is different from the remaining 54
figure;
h1=subplot(2,1,1);
bdp=0.25;
[out]=Sreg(y,X,'nsamp',3000,'bdp',bdp);
resindexplot(out,'h',h1,'conflev',conflev);
ylabel(['Breakdown point =' num2str(bdp)])
h2=subplot(2,1,2);
bdp=0.5;
[out]=Sreg(y,X,'nsamp',3000,'bdp',bdp);
resindexplot(out,'h',h2,'conflev',conflev);
ylabel(['Breakdown point =' num2str(bdp)])


%% SP: analysis using MM estimators
clearvars;close all;
load('hospitalFS.txt');
y=hospitalFS(:,5);
X=hospitalFS(:,1:4);
[out]=MMreg(y,X,'eff',0.85);
% Simulataneous confidence level
conflev = 1- 0.01/length(y);
resindexplot(out,'conflev',conflev);


%% SP: variable selection using added t-tests
clearvars;close all;
load('hospitalFS.txt');
y1=hospitalFS(:,5);
X=hospitalFS(:,1:4);

[outFS]=FSRaddt(y1,X,'init',20,'lms',0,'plots',1);
fieldnames(outFS)

%% LD (Loyalty cards data): fan plot
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4); %#ok<SUSENS>
X=loyalty(:,1:3);
namey='Sales';
nameX={'Number of visits', 'Age', 'Number of persons in the family'};
% yXplot
yXplot(y,X,'nameX',nameX,'namey',namey);

% Compute fan plot to find best value of transformation parameter
[out]=FSRfan(y,X,'plots',1,'la',[-1 -0.5  0 1/4 1/3 0.4 0.5 1]);
% Dynamic Brushing starting from the fan plot
%Example of the use of FlagSize, namey, namex, lwd,FontSize, SizeAxesNum.
%FlagSize controls how large must be the highlighted points. It is a
%parameter of selectdataFS.
fanplot(out,'xlimx',[10 520],'lwd',1.5,'FontSize',11,'SizeAxesNum',11)

%% LD: dynamic brushing from the fan plot with dynamic brushing
% Interactive_example
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);
namey='Sales';
nameX={'Number of visits', 'Age', 'Number of persons in the family'};

% Compute fan plot to find best value of transformation parameter
[out]=FSRfan(y,X,'plots',1,'la',[-1 -0.5  0 1/4 1/3 0.4 0.5 1]);
%FlagSize controls how large must be the highlighted points. It is a
%parameter of selectdataFS.
fanplot(out,'xlimx',[10 520],'lwd',1.5,'FontSize',11,'SizeAxesNum',11,'nameX',nameX,'namey',namey,'databrush',{'selectionmode' 'Brush'...
    'multivarfit' '2' 'FlagSize' '5'})
% If you wish to do persistent brushing from the fan plot
% uncomment the following line. Notice that multiple trajectories can be selected
% fanplot(out,'databrush',{'selectionmode' 'Rect' 'persist' 'on' 'selectionmode','Brush'})


%% LD: forward EDA on transformed data
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);

y1=y.^(0.4);
[out]=LXS(y1,X,'nsamp',10000);
[out]=FSReda(y1,X,out.bs);
% Monitoring scaled residuals with all default parameters
resfwdplot(out);
% Monitoring of minimum deletion residual
mdrplot(out);

%% LD: Interactive monitoring of the trajectories of scaled residuals
% Interactive_example
% using persistent brushing
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);

y1=y.^(0.4);
[out]=LXS(y1,X,'nsamp',10000);
[out]=FSReda(y1,X,out.bs);

databrush=struct;
databrush.bivarfit='2';
databrush.selectionmode='Rect'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='off'; % Write labels of trajectories while selecting
databrush.RemoveLabels='on'; % Do not remove labels after selection
resfwdplot(out,'databrush',databrush);

%% LD: monitoring of modified Cook distance
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);

y1=y.^(0.4);
[out]=LXS(y1,X,'nsamp',10000);
[out]=FSReda(y1,X,out.bs);

plot(out.coo(:,1),out.coo(:,3))
ylabel('Modified Cook distance');
xlabel('Subset size m');
xlim([20 510]);


%% LD: monitoring of "Estimates of beta coefficients"
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);
nameX={'Number of visits', 'Age', 'Number of persons in the family'};

y1=y.^(0.4);
[out]=LXS(y1,X,'nsamp',10000);
[out]=FSReda(y1,X,out.bs);

figure;
for j=3:size(out.Bols,2)
    subplot(2,2,j-2)
    plot(out.Bols(:,1),out.Bols(:,j))
    xlim([10 510])
    xlabel('Subset size m');
    ylabel(nameX(j-2));
end
% text(-13,1.5,'Monitoring of the elements of estimated beta coefficient');

%% LD: Monitoring of "Normality test"
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);

y1=y.^(0.4);
[out]=LXS(y1,X,'nsamp',10000);
[out]=FSReda(y1,X,out.bs);

figure;
xlimx=[100 510];
subplot(2,2,1);
plot(out.nor(:,1),out.nor(:,2));
title('Asymmetry test');
xlim(xlimx);

subplot(2,2,2)
plot(out.nor(:,1),out.nor(:,3))
title('Kurtosis test');
xlim(xlimx);

subplot(2,2,3:4)
plot(out.nor(:,1),out.nor(:,4))
xlim(xlimx);
title('Normality test');
xlabel('Subset size m');


%% LD: Monitoring of deletion t statistics
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);
y1=y.^(0.4);

[out]=FSRaddt(y1,X,'plots',1,'quant',[0.025 0.975]);
fieldnames(out)

%% LD: Automatic outlier detection procedure on transformed data
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);
y1=y.^(0.4);
nameX={'Number of visits', 'Age', 'Number of persons in the family'};

namey1='Sales^{0.4}';
[outFS]=FSR(y1,X,'namey',namey1,'nameX',nameX);

%% LD: analysis using S estimators and brushing
% Interactive_example
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);
y1=y.^(0.4);
nameX={'Number of visits', 'Age', 'Number of persons in the family'};

namey1='Sales^{0.4}';

[out]=Sreg(y1,X,'nsamp',3000,'yxsave',1);
resindexplot(out,'databrush',1,'nameX',nameX,'namey',namey1);

%% LD: analysis using MM estimators and brushing
% Interactive_example
clearvars;close all;
load('loyalty.txt');
y=loyalty(:,4);
X=loyalty(:,1:3);
y1=y.^(0.4);
nameX={'Number of visits', 'Age', 'Number of persons in the family'};

namey1='Sales^{0.4}';

[out]=MMreg(y1,X,'Snsamp',3000,'yxsave',1);
resindexplot(out,'databrush',1,'nameX',nameX,'namey',namey1);



%% SL Salinity data: fan plot
clearvars;close all;
load('salinity.txt');
y=salinity(:,4);
X=salinity(:,1:3);
namey='y: biweekly average salinity';
nameX={'x1: salinity lagged two weeks', 'x2: trend', 'x3: water flow'};
% yXplot
yXplot(y,X,'nameX',nameX,'namey',namey);

% Compute fan plot using 5 most common values of lambda
% This dataset is an example of data which do not have to be transformed
[out]=FSRfan(y,X,'plots',1);
fieldnames(out)

%% SL: forward EDA
clearvars;close all;
load('salinity.txt');
y=salinity(:,4);
X=salinity(:,1:3);
[out]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,out.bs);
% Monitoring of residuals with datatooltip
resfwdplot(out,'datatooltip',1,'tag','resfwdplot')
% levfwdplot(out,'datatooltip',1,'tag','levfwdplot','selunit','0.25')


%% SL: Interactive monitoring of the trajectories of scaled residuals using
% Interactive_example

% Interactive monitoring of the trajectories of scaled residuals using
% persistent brushing
clearvars;close all;
load('salinity.txt');
y=salinity(:,4);
X=salinity(:,1:3);
[out]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,out.bs);

databrush=struct;
databrush.bivarfit='0';
databrush.selectionmode='Rect'; % Brush selection
databrush.persist='off'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='on'; % Do not remove labels after selection

% Highlight trajectories for units 15, 16  and 17
% Note that units 15 and 17 show a large residual only in the final steps
% of the search when unit 16 comes in
fground=struct;
fground.funit=15:17;
fground.LineWidth=3;
resfwdplot(out,'databrush',databrush,'fground',fground);


%% SL: analysis using S and MM estimators
clearvars;close all;
load('salinity.txt');
y=salinity(:,4);
X=salinity(:,1:3);

% Notice that the pattern of residuals  which comes out from the use of
% MMestimators is very different from those of the fwd search and S
% residuals
conflev=1-0.01/length(y);
figure;
h1=subplot(2,1,1);
[out]=MMreg(y,X,'Snsamp',3000);
resindexplot(out,'h',h1,'conflev',conflev);
ylabel('MM residuals')
h2=subplot(2,1,2);
[out]=Sreg(y,X,'nsamp',3000);
resindexplot(out,'h',h2,'conflev',conflev);


%% OD (Ozone data): fan plot
clearvars;close all;
load('ozone.txt','ozone');
y=ozone(:,9);
% Add a time trend to design matrix X
X=[(1:length(y))' ozone(:,1:8)];
% Compute fan plot using 5 most common values of lambda
% Log transformation is clearly the best
[out]=FSRfan(y,X,'plots',1);
fieldnames(out)

%% OD: forward analysis with residuals brushing
% Interactive_example
clearvars;close all;
load('ozone.txt','ozone');
y=ozone(:,9);
% Add a time trend to design matrix X
X=[(1:length(y))' ozone(:,1:8)];

[out]=LXS(log(y),X,'nsamp',10000);
[out]=FSReda(log(y),X,out.bs);

% Interactive monitoring of the trajectories of scaled residuals using
% persistent brushing
databrush=struct;
databrush.bivarfit='0';
databrush.selectionmode='Rect'; % Brush selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
resfwdplot(out,'databrush',databrush);

%% OD: analysis using S and MM estimators
clearvars;close all;
load('ozone.txt','ozone');
y=ozone(:,9);
% Add a time trend to design matrix X
X=[(1:length(y))' ozone(:,1:8)];

% In this dataset the pattern of residuals  which comes out from the use of
% MMestimator is very similar to that which uses S estimator
conflev=1-0.01/length(y);
figure;
h1=subplot(2,1,1);
[outS]=Sreg(y,X,'nsamp',3000);
resindexplot(outS,'h',h1,'conflev',conflev);
ylabel('S residuals')
h2=subplot(2,1,2);
[outMM]=MMregcore(y,X,outS.beta,outS.scale);
resindexplot(outMM,'h',h2,'conflev',conflev);
ylabel('MM residuals')


%% FP (Fishery product): preliminary analysis
clearvars;close all;
load('fishery.mat');
y=fishery{:,2};
X=fishery{:,1};
% Plot of the original data
plot(X,y,'*');
xlabel('Quantity (Tons)');
ylabel('Values (Thousands of Euros)');

%% FP: Dynamic brushing from the fan plot without persistent option
% Interactive_example
clearvars;close all;
% Multiple trajectories can be selected
load('fishery.mat');
y=fishery{:,2};
X=fishery{:,1};

[out]=FSRfan(y,X,'plots',1,'la',[0 0.5 1]);
fanplot(out,'ylimy',[-40,20],'databrush',{'selectionmode' 'Rect' 'persist' '' 'selectionmode','Brush'},'conflev',1-0.001/length(y))


%% PD (Poison data): Fan plot
clearvars;close all;
load('poison.txt');
y=poison(:,end); %#ok<SUSENS>
X=poison(:,1:6);

[out]=FSRfan(y,X,'plots',1,'intercept',0,'ylimy',[-14 3]);
la=[-1 -0.5 0 0.5 1];
Unsel=cell2mat(out.Un);
lla=length(la);
nr=size(Unsel,1)/lla;
Un=[Unsel(1:nr,1) reshape(Unsel(:,2),nr,lla)];

%% PDM1 (Singly modified poison data): Fan plot
clearvars;close all;
load('poison.txt');
y=poison(:,end);
X=poison(:,1:6);
y(8)=0.13;
[out]=FSRfan(y,X,'plots',1,'intercept',0,'ylimy',[-11 9]);
fieldnames(out)

%% PDM2 (Doubly modified poison data): Fan plot
clearvars;close all;
load('poison.txt');
y=poison(:,end);
X=poison(:,1:6);
y(8)=0.13;
y(38)=0.14;
[out]=FSRfan(y,X,'plots',1,'intercept',0,'ylimy',[-8 11]);
fieldnames(out)

%% PDM4 (Multiply modified poison data): Fan plot
clearvars;close all;
load('poison.txt');
y=poison(:,end);
X=poison(:,1:6);
y(6)=0.14;
y(9)=0.08;
y(10)=0.07;
y(11)=0.06;
[out]=FSRfan(y,X,'plots',1,'intercept',0,'ylimy',[-10 22],'nsamp',30000,'init',10);
fieldnames(out);

%% PDM4: Forward exploratory analysis
clearvars;close all;
load('poison.txt');
y=poison(:,end);
X=poison(:,1:6);
y(6)=0.14;
y(9)=0.08;
y(10)=0.07;
y(11)=0.06;

y1=y.^(-1);
[out]=LXS(y1,X,'intercept',0);
[out]=FSReda(y1,X,out.bs,'intercept',0);
% Monitoring scaled residuals: a label is written for the residuals greater
% than 2
fground=struct;
fground.fthresh=2;
resfwdplot(out,'fground',fground);

%% PDM4: Automatic outlier detection procedure
clearvars;close all;
load('poison.txt');
y=poison(:,end);
X=poison(:,1:6);
y(6)=0.14;
y(9)=0.08;
y(10)=0.07;
y(11)=0.06;
y1=y.^(-1);

[out]=FSR(y1,X,'intercept',0);
fieldnames(out)

%% CD (Simulated contaminated data): EDA
clearvars;close all;
n=100;
p=5;

rng(1,'shr3cong');
% Remark if you run this example with a version of MATLAB older than 20011b
% you have to comment the previous four lines and uncomment the following
% lines
% state=1;
% mtstream = RandStream('shr3cong','Seed',state);
% RandStream.setDefaultStream(mtstream);
% defaultStream = RandStream.getDefaultStream();
% reset(defaultStream)

% Remark if you run this example with a version of MATLAB older than 7.9
% you have to comment the previous four lines and uncomment the following
% line
% state=1;
% randn('state', state);rng(1,'shr3cong');
% Remark if you run this example with a version of MATLAB older than 7.9
% you have to comment the previous four lines and uncomment the following
% line
% randn('state', state);
X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;


% LMS using 10000 subsamples
[outLXS]=LXS(y,X,'nsamp',10000,'lms',0);
% Forward Search
[out]=FSReda(y,X,outLXS.bs);
%out.RES=out.RES.^2;
% Monitoring of minimum deletion residual
mdrplot(out);

%% CD: Monitoring of Cook distance
n=100;
p=5;
rng(1,'shr3cong');

X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;


% LMS using 10000 subsamples
[outLXS]=LXS(y,X,'nsamp',10000,'lms',0);
% Forward Search
[out]=FSReda(y,X,outLXS.bs);

figure;
plot(out.coo(:,1),out.coo(:,2))
ylabel('Cook distance');
xlabel('Subset size m');
xlim([20 100]);

%% CD: Monitoring of modified Cook distance
clearvars;close all;
figure;
n=100;
p=5;
rng(1,'shr3cong');


X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;


% LMS using 10000 subsamples
[outLXS]=LXS(y,X,'nsamp',10000,'lms',0);
% Forward Search
[out]=FSReda(y,X,outLXS.bs);

plot(out.coo(:,1),out.coo(:,3))
ylabel('Modified Cook distance');
xlabel('Subset size m');
xlim([20 100]);

%% CD: Monitoring of "Estimates of beta coefficients"
clearvars;close all;
n=100;
p=5;
rng(1,'shr3cong');


X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;


% LMS using 10000 subsamples
[outLXS]=LXS(y,X,'nsamp',10000,'lms',0);
% Forward Search
[out]=FSReda(y,X,outLXS.bs);

figure;
for j=2:size(out.Bols,2)
    subplot(2,4,j-1)
    plot(out.Bols(:,1),out.Bols(:,j))
    % xlim([10 510])
    xlabel('Subset size m');
    %ylabel(nameX(j-1));
end


%% CD: Monitoring of S2 and R2
clearvars;close all;
n=100;
p=5;
rng(1,'shr3cong');

X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;


% LMS using 10000 subsamples
[outLXS]=LXS(y,X,'nsamp',10000,'lms',0);
% Forward Search
[out]=FSReda(y,X,outLXS.bs);

figure;
subplot(1,2,1)
plot(out.S2(:,1),out.S2(:,2))
xlabel('Subset size m');
ylabel('S2');
subplot(1,2,2)
plot(out.S2(:,1),out.S2(:,3))
xlabel('Subset size m');
ylabel('R2');

%% CD: Monitoring of "Normality test"
clearvars;close all;
n=100;
p=5;

rng(1,'shr3cong');

X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;


% LMS using 10000 subsamples
[outLXS]=LXS(y,X,'nsamp',10000,'lms',0);
% Forward Search
[out]=FSReda(y,X,outLXS.bs);

figure;
lwdenv=2;
xlimx=[10 100];
subplot(2,2,1);
plot(out.nor(:,1),out.nor(:,2));
title('Asymmetry test');
xlim(xlimx);
quant=chi2inv(0.99,1);
v=axis;
line([v(1),v(2)],[quant,quant],'color','r','LineWidth',lwdenv);

subplot(2,2,2)
plot(out.nor(:,1),out.nor(:,3))
title('Kurtosis test');
xlim(xlimx);
v=axis;
line([v(1),v(2)],[quant,quant],'color','r','LineWidth',lwdenv);

subplot(2,2,3:4)
plot(out.nor(:,1),out.nor(:,4));
xlim(xlimx);
quant=chi2inv(0.99,2);
v=axis;
line([v(1),v(2)],[quant,quant],'color','r','LineWidth',lwdenv);

title('Normality test');
xlabel('Subset size m');

%% CD: Automatic outlier detection procedure
% Forward search
clearvars;close all;
n=100;
p=5;

rng(1,'shr3cong');

X = randn(n,p);
y = randn(n,1);
y(1:30) = y(1:30)+5;
out     = FSR(y,X); %#ok<NASGU>

%% CD: LTS and LMS
clearvars;close all;
n=100;
p=5;

rng(1,'shr3cong');
X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;


conflev=[0.99, 1-0.01/length(y)];
% Define the main title of the plots
titl='';

% Define number of subsamples
nsamp=10000;

% LMS with no rewighting
[outLMS]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1));
h1=subplot(2,2,1);
laby='Scaled LMS residuals';
resindexplot(outLMS.residuals,'h',h1,'title',titl,'laby',laby,'numlab','','conflev',conflev)

% LTS with no rewighting
h2=subplot(2,2,2);
[outLTS]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'lms',0);
laby='Scaled LTS residuals';
resindexplot(outLTS.residuals,'h',h2,'title',titl,'laby',laby,'numlab','','conflev',conflev);

% LMS with reweighting
[outLMSr]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1);
h3=subplot(2,2,3);
laby='Scaled reweighted LMS residuals';
resindexplot(outLMSr.residuals,'h',h3,'title',titl,'laby',laby,'numlab','','conflev',conflev)

% LTS with reweighting
[outLTSr]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1,'lms',0);
h4=subplot(2,2,4);
laby='Scaled reweighted LTS residuals';
resindexplot(outLTSr.residuals,'h',h4,'title',titl,'laby',laby,'numlab','','conflev',conflev);


%% CD: LTS and LMS (permuting the order of the contaminated data)
clearvars;close all;
n=100;
p=5;

rng(1,'shr3cong');


X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;

% Permute values in X and y
s=randsample(n,n);
y=y(s);
X=X(s,:);

conflev=[0.99, 1-0.01/length(y)];
% Define the main title of the plots
titl='';

% Define number of subsamples
nsamp=10000;

% LMS with no rewighting
[outLMS]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1));
h1=subplot(2,2,1);
laby='Scaled LMS residuals';
resindexplot(outLMS.residuals.^2,'h',h1,'title',titl,'laby',laby,'numlab','','conflev',conflev)


% LTS with no rewighting
h2=subplot(2,2,2);
[outLTS]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'lms',0);
laby='Scaled LTS residuals';
resindexplot(outLTS.residuals,'h',h2,'title',titl,'laby',laby,'numlab','','conflev',conflev);

% LMS with reweighting
[outLMSr]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1);
h3=subplot(2,2,3);
laby='Scaled reweighted LMS residuals';
resindexplot(outLMSr.residuals,'h',h3,'title',titl,'laby',laby,'numlab','','conflev',conflev)

% LTS with reweighting
[outLTSr]=LXS(y,X,'nsamp',nsamp,'conflev',conflev(1),'rew',1,'lms',0);
h4=subplot(2,2,4);
laby='Scaled reweighted LTS residuals';
resindexplot(outLTSr.residuals,'h',h4,'title',titl,'laby',laby,'numlab','','conflev',conflev);


%% CD: S and MM (permuting the order of the contaminated data)
clearvars;close all;
n=100;
p=5;

rng(1,'shr3cong');


X=randn(n,p);
y=randn(n,1);
y(1:30)=y(1:30)+5;

% Permute values in X and y
s=randsample(n,n);
y=y(s);
X=X(s,:);

conflev=[0.99, 1-0.01/length(y)];
% Define the main title of the plots
titl='';

% Define number of subsamples
nsamp=10000;

% S estimators with bdp=0.5
bdp=0.5;
[outS050]=Sreg(y,X,'nsamp',nsamp,'bdp',bdp);
h1=subplot(2,2,1);
laby=['S residuals bdp=' num2str(bdp)];
resindexplot(outS050,'h',h1,'title',titl,'laby',laby,'numlab','','conflev',conflev)


% S estimators with bdp=0.25
h2=subplot(2,2,2);
bdp=0.25;
[outS]=Sreg(y,X,'nsamp',nsamp,'bdp',bdp);
laby=['S residuals bdp=' num2str(bdp)];
resindexplot(outS,'h',h2,'title',titl,'laby',laby,'numlab','','conflev',conflev);

% MM estimators with eff=0.85
eff=0.85;
[outMM85]=MMregcore(y,X,outS050.beta,outS050.scale,'eff',eff);
h3=subplot(2,2,3);
laby=['MM residuals eff=' num2str(eff)];
resindexplot(outMM85,'h',h3,'title',titl,'laby',laby,'numlab','','conflev',conflev)

% MM estimators with eff=0.95
eff=0.95;
[outMM95]=MMregcore(y,X,outS050.beta,outS050.scale,'eff',eff);
h4=subplot(2,2,4);
laby=['MM residuals eff=' num2str(eff)];
resindexplot(outMM95,'h',h4,'title',titl,'laby',laby,'numlab','','conflev',conflev)


%% ST: Stars dataset (analysis using FS, exploratory data analysis)
clearvars;close all;
stars=load('stars.txt');
y=stars(:,2);
X=stars(:,1);
plot(X,y,'o');
xlabel('Log effective surface temperature')
ylabel('Log light intensity')
set(gca,'XDir','reverse');

% LMS using 1000 subsamples
[out]=LXS(y,X,'nsamp',10000);
% Forward Search
[out]=FSReda(y,X,out.bs);
out1=out;
% Create scaled squared residuals
out1.RES=out.RES.^2;
% resfwdplot with all default options
resfwdplot(out1);

%% ST: Stars dataset (analysis using FS, automatic outlier detection)
clearvars;close all;
stars=load('stars.txt');
y=stars(:,2);
X=stars(:,1);
plot(X,y,'o');
xlabel('Log effective surface temperature')
ylabel('Log light intensity')
set(gca,'XDir','reverse');
[out]=FSR(y,X,'plots',2);

%% ST: Stars dataset (analysis using S and MM)
clearvars;close all;
stars=load('stars.txt');
y=stars(:,2);
X=stars(:,1);

conflev=[0.99, 1-0.01/length(y)];
% Define the main title of the plots
titl='';

% Define number of subsamples
nsamp=1000;

% S estimators with bdp=0.5
bdp=0.5;
[outS050]=Sreg(y,X,'nsamp',nsamp,'bdp',bdp);
h1=subplot(2,2,1);
laby=['S residuals bdp=' num2str(bdp)];
resindexplot(outS050,'h',h1,'title',titl,'laby',laby,'numlab','','conflev',conflev)


% S estimators with bdp=0.25
h2=subplot(2,2,2);
bdp=0.25;
[outS]=Sreg(y,X,'nsamp',nsamp,'bdp',bdp);
laby=['S residuals bdp=' num2str(bdp)];
resindexplot(outS,'h',h2,'title',titl,'laby',laby,'numlab','','conflev',conflev);

% MM estimators with eff=0.85
eff=0.85;
[outMM85]=MMregcore(y,X,outS050.beta,outS050.scale,'eff',eff);
h3=subplot(2,2,3);
laby=['MM residuals eff=' num2str(eff)];
resindexplot(outMM85,'h',h3,'title',titl,'laby',laby,'numlab','','conflev',conflev)

% MM estimators with eff=0.95
eff=0.95;
[outMM95]=MMregcore(y,X,outS050.beta,outS050.scale,'eff',eff);
h4=subplot(2,2,4);
laby=['MM residuals eff=' num2str(eff)];
resindexplot(outMM95,'h',h4,'title',titl,'laby',laby,'numlab','','conflev',conflev)

%% TCLUST-REG: choice of the hyperparameters through monitoring
%
% This example shows how to use in an automatic way the monitoring
% functions of FSDA associated to clustering methods (TCLUST family). In
% the example:
% - We start using tclustregIC to monitor several choices for the number of
%   groups and restriction factor values, fixing a reasonable trimming
%   level.
% - Then we extract a set of relevant solutions with tclustICsol and
%   we visualize them with carbikeplot.
% - We use the information extracted by tclustICsol to identify "the best"
%   solution. Intuitively, the best solution corresponds to the "car" of
%   bigger area among those represented in the carbikeplot.
% - Given the best solution for k and c, we run tclustregeda to monitor what
%   happens for different possible trimming levels, and we choose the best
%   solution based on the ARI indexes of the various clusterings. If the
%   ARI indexes are similar, we choose the solution that preserves data and
%   gives more efficiency to the estimator.
%
% do not plot un-necessary graphics
doPlots = false;
% generate data with mixsim
rng(372,'twister');
p=3;
k=3;
Q=MixSimreg(k,p,'MaxOmega',0.00001,'restrfactor',2);
n=300;
cont=30;
[y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib,'noiseunits',cont);
% plot the generated data
y(id==-1) = y(id==-1) + 2*randn(cont,1);
yXplot(y,X,id);

% Use tclustregIC to monitor the effect of k and c, for alpha fixed (ssume
% it is overestimated)
alpha = 0.15;
kvec  = 2:1:7;
cvec  = [1,2,4,8,16,32,64];
outIC = tclustregIC(y,X,'whichIC','CLACLA','kk',kvec,'cc',cvec,'alphaLik',alpha,'plots',0);

% Extracts a set of best relevant solutions ...
[outICsol] = tclustICsol(outIC,'whichIC','CLACLA','plots',doPlots,'NumberOfBestSolutions',5,'ThreshRandIndex',0.7);

% ... and visualise them with the carbike plot, which highlights the most
% relevant one in intuitive way.
[hcb,areas] = carbikeplot(outICsol,'SpuriousSolutions',true);

% Use the information extracted by tclustICsol to identify the best
% solution.
[truesol,~]  = ismember(outICsol.CLACLAbs(:,end),'true');
truesoli     = find(truesol);
[amax,iamax] = max(areas(truesoli,2));
if numel(truesoli) > 0
    if amax>0
        % take the true solution with larger area
        kopt = outICsol.CLACLAbs{truesoli(iamax),1};  % optimal number of groups
        copt = outICsol.CLACLAbs{truesoli(iamax),2};  % optimal nrestriction factor
    else
        % if areas are all zero, take the true solution with larger k
        kopt = outICsol.CLACLAbs{truesoli(1),1};
        copt = outICsol.CLACLAbs{truesoli(1),2};
    end
else
    % if there are no true solutions, take the one with larger k
    kopt = outICsol.CLACLAbs{truesol(1),1};
    copt = outICsol.CLACLAbs{truesol(1),2};
end

% Finally, use tclustregeda to monitor alpha, with k and c estimated by tclustregIC
alphaLikvec   = 0.15:-0.05:0;
outEDA        = tclustregeda(y,X,kopt,copt,alphaLikvec,0,'plots',doPlots);

% retrieve the optimal alpha
[ARImax]   = max(outEDA.Amon(:,2));
iARImax    = find(outEDA.Amon(:,2) == ARImax);
alphaopt   = outEDA.Amon(iARImax(end),1);

% clustering on the first PC scores
% With the parameters sugggested by the monitoring the main structure of the data emerges
outTCLUST    = tclustreg(y,X,kopt,copt,alphaopt,0,'plots',1);
idxTCLUST    = outTCLUST.idx;

%% As above, for the Fishery2003 data

% Monitoring of fishery2003 with or without intercept.
clear all; close all;
rng(123)

intercept = 0;

load Fishery2003.mat;
X = Fishery2003{:,2};
X = X + 10^(-5) * abs(randn(size(X)));
y = Fishery2003{:,3};
y = y + 10^(-5) * abs(randn(size(y)));
id = Fishery2003{:,1};

% yXplot(yid,Xid);

% Use tclustregIC to monitor the effect of k and c, for alpha reasonably fixed
alpha = 0.01;
kvec  = 1:1:4;
cvec  = [1,2,4,8];
outIC = tclustregIC(y,X,'whichIC','CLACLA','kk',kvec,'cc',cvec,'alphaLik',alpha,'intercept',intercept,'plots',0);

% Extracts a set of best relevant solutions ...
[outICsol] = tclustICsol(outIC,'whichIC','CLACLA','plots',0,'NumberOfBestSolutions',5,'ThreshRandIndex',0.7);

% ... and visualise them with the carbike plot, which highlights the most
% relevant one in intuitive way.
[hcb,areas] = carbikeplot(outICsol,'SpuriousSolutions',false);

% Use the information extracted by tclustICsol to identify the best
% solution.
[truesol,~]  = ismember(outICsol.CLACLAbs(:,end),'true');
truesoli     = find(truesol);
[amax,iamax] = max(areas(truesoli,2));
if numel(truesoli) > 0
    if amax>0
        % take the true solution with larger area
        kopt = outICsol.CLACLAbs{truesoli(iamax),1};  % optimal number of groups
        copt = outICsol.CLACLAbs{truesoli(iamax),2};  % optimal nrestriction factor
    else
        % if areas are all zero, take the true solution with larger k
        kopt = outICsol.CLACLAbs{truesoli(1),1};
        copt = outICsol.CLACLAbs{truesoli(1),2};
    end
else
    % if there are no true solutions, take the one with larger k
    kopt = outICsol.CLACLAbs{truesol(1),1};
    copt = outICsol.CLACLAbs{truesol(1),2};
end

% Finally, use tclustregeda to monitor alpha, with k and c estimated by tclustregIC
alphaLikvec   = 0.05:-0.01:0;
outEDA        = tclustregeda(y,X,kopt,copt,alphaLikvec,0,'intercept',intercept,'plots',0);

% retrieve the optimal alpha
[ARImax]   = max(outEDA.Amon(:,2));
iARImax    = find(outEDA.Amon(:,2) == ARImax);
alphaopt   = outEDA.Amon(iARImax(end),1);

% final clustering % with the parameters sugggested by the monitoring
outTCLUST    = tclustreg(y,X,kopt,copt,alphaopt,0,'intercept',intercept,'plots',1);
idxTCLUST    = outTCLUST.idx;

