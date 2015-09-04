%examples_multivariate shows a series of analysis of multivariate datasets
% Copyright 2008-2015.
% Written by FSDA team

% Last modified 06-Feb-2015

%% HD (Heads data) analysis using univariate boxplots
clearvars;close all;
load('head.mat');
Y=head.data;
cnames=head.colnames;
% Compare the output with Figure 3.9 p. 97
boxplot(Y,'labels',cnames,'LabelOrientation','inline');
% Label the outliers with the unit number
%Find in current plot the handles associated with the univariate outliers
a=findobj(gca,'tag','Outliers');
[n,v]=size(Y);
for j=1:v;
    % Get the X and Y coordinates of the univariate outliers
    aYdata=get(a(j),'Ydata');
    aXdata=get(a(j),'Xdata');
    
    % Loop over the outliers inside variable v-j+1
    for i=1:length(aYdata);
        % Find the row number of the outliers
        ind=find(aYdata(i)==Y(:,end-j+1),n);
        % Add the label to the current plot
        text(aXdata(1),aYdata(i),num2str(ind'))
    end
end

%% HD (Heads data) analysis using S and MM estimators
clearvars;close all;
load('head.mat');
Y=head.data;
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')


%% HD (Heads data) -- Forward EDA (Exploratory Data Analysis):
clearvars;close all;
% scatterplot of data: two points look outlying
load('head.mat');
Y=head.data;
cnames=head.colnames;
% Scatter plot matrix
gplotmatrix(Y,[],[],[],[],[],[],[],cnames);
set(findall(gcf,'Type','text'),'HorizontalAlignment','Right','rotation',45,'Interpreter','None','FontName','Arial')


%% HD qqplot based on quantiles of beta distribution
% compare the output with Figure 1.2, p.9 or ARC(2004)
clearvars;close all;
figure;
load('head.mat');
Y=head.data;
[n,v]=size(Y);

cor=(n/((n-1)^2));

% malas= vector of scaled Mahalanobis distances (MD)
malas=mahal(Y,Y)*cor;
% Set the beta
pd = ProbDistUnivParam('beta',[v/2 0.5*(n-v-1)]);
% qqplot of scaled MD against beta distribution
h=qqplot(malas,pd);
xplo=get(h(1),'XData')';

nsimul=100;
malasim=zeros(n,nsimul);
for j=1:nsimul;
    Ysim=randn(n,v);
    % In each simulation store scaled MD
    malasim(:,j)=sort(cor*mahal(Ysim,Ysim));
end

Yband=sort(malasim,2);
quant=[0.05 0.95];
% add lower and upper confidence band
line(xplo,Yband(:,nsimul*quant),'color','g');

%% HD Preliminary analysis based on robust bivariate ellipses and minMD
% Figure 3.2
clearvars;close all;
load('head.mat');
Y=head.data;
[fre]=unibiv(Y,'plots',1,'textlab',1,'rf',0.999);

fre=sortrows(fre,4);

bs=fre(1:size(Y,2),1);
% Forward search with EDA purposes and plot of minimum MD
% Figure 1.3
[out]=FSMeda(Y,bs,'plots',1,'init',60); %#ok<*NASGU>

%% HD Plot of max MD and gap
% Interactive_example
clearvars;close all;
load('head.mat');
Y=head.data;

[fre]=unibiv(Y,'rf',0.75);
bs=fre((fre(:,2)==0 & fre(:,3)==0),1);
[out]=FSMeda(Y,bs,'plots',1,'init',68);

figure;
% Figure 3.2, p. 92
msr=out.msr;
% Top panel: max and mth
subplot(2,1,1);
plot1=plot(msr(:,1),msr(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')

h=legend('Maximum MD','m-th ordered MD','Location','NorthWest');
set(h,'FontSize',14);

% Bottom panel: gap and (m+1)th - mth
% Figure 3.4 p. 93
subplot(2,1,2);
gap=out.gap;
plot1=plot(gap(:,1),gap(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')
xlabel('Subset size m')
h=legend('gap','(m+1)th ordered MD - (m)th ordered MD','Location','SouthEast');
set(h,'FontSize',14);


%% HD persistent brushing from the malfwdplot (scaled distances)
% Interactive_example
clearvars;close all;
load('head.mat');
Y=head.data;
[fre]=unibiv(Y,'rf',0.75);
bs=fre(1:size(Y,2)+1,1);
[out]=FSMeda(Y,bs,'init',60,'scaled',1);

% Plot of minimum scaled Mahalanobis distance
mmdplot(out,'scaled',1);
% Now, some interactive brushing starting from the monitoring residuals
% plot. Once a set of trajectories is highlighted in the monitoring residual plot,
% the corresponding units are highlighted in the other plots
databrush=struct;
databrush.selectionmode='Rect'; % Rectangular selection
databrush.persist='on'; % Enable repeated mouse selections
databrush.Label='on'; % Write labels of trajectories while selecting
databrush.RemoveLabels='off'; % Do not remove labels after selection
% Compare with Figure 3.6 p. 95
malfwdplot(out,'databrush',databrush);

%% HD: EDA monitoring of the estimated correlation matrix
clearvars;close all;
load('head.mat');
Y=head.data;
Y1=zscore(Y);
[fre]=unibiv(Y1,'rf',0.75);
bs=fre(1:size(Y,2)+2,1);
[out]=FSMeda(Y1,bs,'init',60);

S2cov=out.S2cov;
plot1=plot(S2cov(:,1),S2cov(:,2:end));
xlabel('Subset size m')
ylabel('Elements of correlation matrix')
v=size(Y,2);
v1=v*(v+1)/2;
slintyp={'--' '-' '-.' ':'}';
slintyp=repmat(slintyp,ceil(v1/length(slintyp)),1);
set(plot1,{'LineStyle'},slintyp(1:v1));
%% HD: analysis of transformations
% FS based on untransformed data H_0:\lambda=1 for all variables
clearvars;close all;
load('head.mat');
Y=head.data;
% Monitoring of likelihood ratio test
% Compare the output with Figure 4.13 p. 172 of ARC (2004)
[out]=FSMtra(Y,'plotslrt',1);

%% HD: analysis of transformations
clearvars;close all;
load('head.mat');
Y=head.data;
% FS based on untransformed data H_0:\lambda=1 for variable 4
% Monitoring of likelihood ratio test
% Compare the output with Figure 4.14 p. 173 of ARC (2004)
[out]=FSMtra(Y,'ColToTra',4,'plotslrt',1);


%% HD: confirmation of transformation
clearvars;close all;
% Compare the output with Figure 4.15 p. 174 of ARC (2004)
load('head.mat');
Y=head.data;
v=size(Y,2);
plotslrt=struct;
plotslrt.ylim=[-3.2 3.2];
plotslrt.xlim=[110 200];
[out]=FSMfan(Y,ones(v,1),'init',110,'plotslrt',plotslrt);

%% HD: confirmation of transformation
% Compare the output with Figure 4.65 p. 222 of ARC (2004)
clearvars;close all;
load('head.mat');
Y=head.data;
n=size(Y,1);
[out]=FSRfan(Y(:,4),ones(n,1),'plots',1,'intercept',0,'xlimx',[60 210],'ylimy',[-5 4.4]);


%% HD: random starts
clearvars;close all;
load('head.mat');
Y=head.data;
[n,v]=size(Y);

init=20;
nsimul=200;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmdeasy(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));

% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');


%% TR: (Track records): spm
clearvars;close all;
close all;
load('recordfg');
Y=recordfg.data;
% Scatter plot matrix (Figure 1.4 of ARC)
gplotmatrix(Y);

%% TR: qqplot of MD based on the beta distribution
% compare the output with Figure 1.2, p.9 or ARC(2004)
clearvars;close all;
load('recordfg');
Y=recordfg.data;
[n,v]=size(Y);

cor=(n/((n-1)^2));

% malas= vector of scaled Mahalanobis distances (MD)
malas=mahal(Y,Y)*cor;
% Set the beta
pd = ProbDistUnivParam('beta',[v/2 0.5*(n-v-1)]);
% qqplot of scaled MD against beta distribution
h=qqplot(malas,pd);
xplo=get(h(1),'XData')';

nsimul=1000;
malasim=zeros(n,nsimul);
for j=1:nsimul;
    Ysim=randn(n,v);
    % In each simulation store scaled MD
    malasim(:,j)=sort(cor*mahal(Ysim,Ysim));
end

Yband=sort(malasim,2);
quant=[0.05 0.95];
% add lower and upper confidence band
line(xplo,Yband(:,nsimul*quant),'color','g');


%% TR: Forward EDA
% Preliminary analysis based on robust bivariate ellipses
clearvars;close all;
load('recordfg');
Y=recordfg.data;
[fre]=unibiv(Y);
fre=sortrows(fre,4);

bs=fre(1:size(Y,2),1);
% Forward search with EDA purposes and plot of minimum MD
% Create Figure 1.6 p. 14
[out]=FSMeda(Y,bs,'plots',1,'init',10);

%% TR: plot of scaled Mahalanobis distances
% compare with Figure 3.12 p. 100
clearvars;close all;
load('recordfg');
Y=recordfg.data;
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'scaled',1,'init',15);
% [out]=FSMeda(Y,bs,'init',10);
standard=struct;
standard.ylim=[0 22];
standard.LineStyle={'--' '-' '-.' ':'};    % different line styles for different standard trajectories
standard.Color={'b';'g';'c';'m';'y';'k'};  % different colors for different standard trajectories
malfwdplot(out,'standard',standard);

%% TR: plot of max and mth distance
clearvars;close all;
load('recordfg');
Y=recordfg.data;
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',15);
% Figure 3.16, p. 104
msr=out.msr;
% Top panel: max and mth
subplot(2,1,1);
plot1=plot(msr(:,1),msr(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')

h=legend('Maximum MD','m-th ordered MD','Location','NorthWest');
set(h,'FontSize',14);

% Bottom panel: gap and (m+1)th - mth
% Figure 3.17 p. 104
subplot(2,1,2);
gap=out.gap;
plot1=plot(gap(:,1),gap(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')
xlabel('Subset size m')
h=legend('gap','(m+1)th ordered MD - (m)th ordered MD','Location','NorthWest');
set(h,'FontSize',14);

%% TR: automatic outlier detection procedure on original data
clearvars;close all;
load('recordfg');
Y=recordfg.data;
FSM(Y,'plots',2)


%% TR: analysis using S and MM estimators
clearvars;close all;
load('recordfg');
Y=recordfg.data;
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')


%% TR: analysis of transformations
% Track records
clearvars;close all;
load('recordfg');
Y=recordfg.data;
n=size(Y,1);
Y1=repmat(max(Y),n,1);
Y=Y./Y1;
la0=[-1 -3];
tagsmle={'mle-1' 'mle-3'};
plotsmle=struct;
plotsmle.ylim=[-4 0];

ii=1;
for la=la0;
    plotsmle.Tag=tagsmle{ii};
    FSMtra(Y,'plotsmle',plotsmle,'onelambda',1,'la0',la,'init',20);
    ii=ii+1;
end
% Compare these 2 plots with Figure 4.70 p. 225 of ARC (2004)

%% TR: analysis of transformations
% Track records
clearvars;close all;
load('recordfg');
Y=recordfg.data;
n=size(Y,1);
Y1=repmat(max(Y),n,1);
Y=Y./Y1;
la0=[-1 -2 -3 -4];
% If lik=1 it is also possible to monitor the plots of mle of transformation parameters
lik=0;

tagslrt={'lrt-1' 'lrt-2' 'lrt-3' 'lrt-4'};
tagsmle={'mle-1' 'mle-2' 'mle-3' 'mle-4'};
plotsmle=struct;
plotsmle.ylim=[-4 0];

plotslrt=struct;
plotslrt.ylim=[0 21];

ii=1;
for la=la0;
    plotsmle.Tag=tagsmle{ii};
    plotslrt.Tag=tagslrt{ii};
    if lik==1;
        FSMtra(Y,'plotsmle',plotsmle,'plotslrt',plotslrt,'onelambda',1,'la0',la,'init',20);
    else
        FSMtra(Y,'plotslrt',plotslrt,'onelambda',1,'la0',la,'init',20);
    end
    ii=ii+1;
end
% Compare these 4 plots with Figure 4.50 p. 207 of ARC (2004)



%% TR: confirmation of transformations
% confirmatory search
% Compare the plot with Figure 4.51 p. 208 of ARC (2004)
clearvars;close all;
load('recordfg');
Y=recordfg.data;
[n , v]=size(Y);
Y1=repmat(max(Y),n,1);
Y=Y./Y1;
% FS based on with H_0:\lambda=-3 for all variables
plotslrt=struct;
plotslrt.ylim=[-6.2 3.2];
plotslrt.xlim=[28 59];
laAround=-4:0;
[out]=FSMfan(Y,-3*ones(v,1),'ColToComp',[2 3 5 7],'laAround',laAround,'init',28,'plotslrt',plotslrt);


%% TR transformed: parallel coordinates plot
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;
n=size(Y,1);
group=ones(n,1);
group([14 33 55 36])=2;
parallelcoords(Y,'LineWidth',1.5,'Standardize','on','Group',group)

%% TR transformed: malfwdplot
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;

[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'scaled',1,'init',15);
% [out]=FSMeda(Y,bs,'init',10);
standard=struct;
standard.ylim=[0 10];
standard.LineStyle={'--' '-' '-.' ':'};    % different line styles for different standard trajectories
standard.Color={'b';'g';'c';'m';'y';'k'};  % different colors for different standard trajectories
malfwdplot(out,'standard',standard);

%% TR transformed: spm
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;
cnames=recordfg.colnames;
% Scatter plot matrix
gplotmatrix(Y,[],[],[],[],[],[],[],cnames);

%% TR transformed: plot of max inside and min outside
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;
figure;
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',15);
% Figure 4.71 left panel
msr=out.msr;
subplot(2,1,1);
plot1=plot(msr(:,1),msr(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')

h=legend('Maximum MD','m-th ordered MD','Location','NorthWest');
set(h,'FontSize',14);

% Bottom panel: min and (m+1)th
% Figure 4.71 (right panel) p. 226
subplot(2,1,2);
mmd=out.mmd;
plot1=plot(mmd(:,1),mmd(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')
xlabel('Subset size m')
h=legend('min','(m+1)th ordered MD','Location','NorthWest');
set(h,'FontSize',14);


%% TR transformed: automatic outlier detection procedure (part I)
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;
FSM(Y,'init',40,'plots',2)

%% TR transformed: automatic outlier detection procedure (part II)
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;
FSM(Y,'plots',2)


%% TR transformed: random starts
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;
[n,v]=size(Y);

init=20;
nsimul=200;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));


% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');

%% TR transformed: analysis using S and MM estimators
clearvars;close all;
load('recordfg');
Y=recordfg.data;
Y=Y.^-3;
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')

%% SB: (Swiss banknotes) Forward EDA minMD
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
group=ones(200,1);
group(101:200)=2;
spmplot(Y,group)


%% SB: (Swiss banknotes) Forward EDA minMD
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',30);
mmdplot(out,'mplus1',1);


%% SB: Forward EDA malfwdplot
% Create figure 3.28 p.117 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',30,'scaled',1);
fground=struct;
fground.flabstep='';
malfwdplot(out,'fground',fground);

%% SB: Forward EDA gapplot
% Create figure 3.29 p.117 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',30);
gap=out.gap;
plot1=plot(gap(:,1),gap(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')
h=legend('min outside - max inside','(m+1)th ordered MD - mth ordered MD','Location','NorthWest');
set(h,'FontSize',14);

%% SB:  analysis using S and MM estimators
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')


%% SB: Forward EDA gapplot starting with the first 20 obs (geniune notes)
% Create figure 3.31 p.118 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
% start in the group of geniune notes
bs=1:20;
[out]=FSMeda(Y,bs,'init',30);
gap=out.gap;
plot1=plot(gap(:,1),gap(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')
h=legend('min outside - max inside','(m+1)th ordered MD - mth ordered MD','Location','NorthWest');
set(h,'FontSize',14);

%% SB: Forward EDA gapplot starting with the obs 101:120 (forged notes)
% Create figure 3.36 p.124 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
% start in the group of forgeries
bs=(1:20)+100;
[out]=FSMeda(Y,bs,'init',30);
gap=out.gap;
plot1=plot(gap(:,1),gap(:,2:end),'k','LineWidth',2);
set(plot1(1),'LineStyle',':')
h=legend('min outside - max inside','(m+1)th ordered MD - mth ordered MD','Location','NorthWest');
set(h,'FontSize',14);

%% SB: monitoring of MD starting in the genuine group
% Create figure 3.30 p.118 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
% start in the group of geniune notes
bs=1:20;
[out]=FSMeda(Y,bs,'init',30,'scaled',1);
malfwdplot(out);

%% SB: monitoring of MD starting in the group of forged notes
% Create figure 3.35 p.123 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
% start in the group of geniune notes
bs=101:120;
[out]=FSMeda(Y,bs,'init',30,'scaled',1);
fground=struct;
fground.flabstep='';
malfwdplot(out,'fground',fground);

%% SB: Forward EDA compare scaled and unscaled MD
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
fground=struct;
fground.flabstep='';
% Start the search with the first 20 observations of the forged notes
bs=1:20;
% Create upper panel of figure 2.5 p. 67 of ARC 2004
[outsc]=FSMeda(Y,bs,'init',20,'scaled',1);
malfwdplot(outsc,'tag','scaled','fground',fground);

% Create lower panel of figure 2.5 p. 67 of ARC 2004
[out]=FSMeda(Y,bs,'init',20);
malfwdplot(out,'tag','unscaled','fground',fground);

%% SB: random starts
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
[n,v]=size(Y);

init=20;
nsimul=200;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));


% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');

%% SB (genuine notes): Forward EDA malfwdplot
% Create figure 3.41 p.129 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=Y(1:100,:);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',30,'scaled',1);
malfwdplot(out);

%% SB (genuine notes): Forward EDA maxMDplot
% Create figure 3.41 p.129 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=Y(1:100,:);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',30);
msr=out.msr;
plot1=plot(msr(:,1),msr(:,2:end),'k','LineWidth',2);
set(plot1(2),'LineStyle',':')
h=legend('mth ordered MD','max inside','Location','NorthWest');
set(h,'FontSize',14);


%% SB (genuine notes): elements of cov
% Create figure 3.43 p. 130 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
% if corr=1 it is possible to monitor the elements of the correlation
% matrix
corr=0;
if corr==1;
    Y1=zscore(Y(1:100,:));
else
    Y1=Y(1:100,:);
end

[fre]=unibiv(Y1,'rf',0.75);
bs=fre(1:size(Y,2)+2,1);
[out]=FSMeda(Y1,bs,'init',60);

S2cov=out.S2cov;
plot1=plot(S2cov(:,1),S2cov(:,2:end));
xlabel('Subset size m')
if corr==1;
    ylabel('Elements of correlation matrix')
else
    ylabel('Elements of covariance matrix')
end
v=size(Y,2);
v1=v*(v+1)/2;
slintyp={'--' '-' '-.' ':'}';
slintyp=repmat(slintyp,ceil(v1/length(slintyp)),1);


set(plot1,{'LineStyle'},slintyp(1:v1));


aco=triu(ones(v,v));
ind = find(abs(aco)>0);
[I,J]=ind2sub(v,ind);
lab=cellstr([num2str(I) num2str(J)]);
text(S2cov(end,1)*ones(v1,1),S2cov(end,2:end)',lab)

%% SB (genuine notes): univariate boxplots
% Create figure 3.44 p.129 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=zscore(Y(1:100,:));
cnames=swiss_banknotes.colnames;
% Compare the output with Figure 3.9 p. 97
boxplot(Y,'labels',cnames);
% Label the outliers with the unit number
%Find in current plot the handles associated with the univariate outliers
a=findobj(gca,'tag','Outliers');
[n,v]=size(Y);
for j=1:v;
    % Get the X and Y coordinates of the univariate outliers
    aYdata=get(a(j),'Ydata');
    aXdata=get(a(j),'Xdata');
    
    % Loop over the outliers inside variable v-j+1
    for i=1:length(aYdata);
        % Find the row number of the outliers
        ind=find(aYdata(i)==Y(:,end-j+1),n);
        % Add the label to the current plot
        text(aXdata(1),aYdata(i),num2str(ind'))
    end
end

%% SB (genuine notes): automatic outlier detection procedure
% Create figure 3.41 p.129 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=Y(1:100,:);
[out]=FSM(Y,'plots',2);

%% SB (genuine notes): brushing from mmdplot
% Interactive_example
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=Y(1:100,:);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',30,'scaled',0);
databrush=struct;
databrush.labeladd='1';
databrush.Label='on';
databrush.RemoveLabels  = 'off';
databrush.BrushShape='rect';
% Control the characteristics of the units in foreground in the resfwdplot
fground=struct;
% Labels will be included for units (these are the units entering the last
% six steps of the search
fground.funit=[1 13 40 70 71 5];
% FontSize of the Labels
fground.FontSize=14;
% The Labels will be put when m=30
fground.flabstep=30;
standard.xlim=[28 100];
malfwdplot(out,'standard',standard,'fground',fground);
mmdplot(out,'databrush',databrush);

%% SB:  analysis using S and MM estimators
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=Y(1:100,:);
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')



%% SB (forged notes): elements of cov
% Create figure 3.43 p. 130 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
% if corr=1 it is possible to monitor the elements of the correlation
% matrix
corr=0;
if corr==1;
    Y1=zscore(Y(101:200,:));
else
    Y1=Y(101:200,:);
end

[fre]=unibiv(Y1,'rf',0.75);
bs=fre(1:size(Y,2)+2,1);
[out]=FSMeda(Y1,bs,'init',60);

S2cov=out.S2cov;
plot1=plot(S2cov(:,1),S2cov(:,2:end));
xlabel('Subset size m')
if corr==1;
    ylabel('Elements of correlation matrix')
else
    ylabel('Elements of covariance matrix')
end
v=size(Y,2);
v1=v*(v+1)/2;
slintyp={'--' '-' '-.' ':'}';
slintyp=repmat(slintyp,ceil(v1/length(slintyp)),1);
set(plot1,{'LineStyle'},slintyp(1:v1));



%% SB (forged notes): Forward EDA malfwdplot
% Create figure 3.46 p.132 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=Y(101:200,:);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',30,'scaled',1);
malfwdplot(out);

%% SB (forged notes): automatic outlier detection procedure
% Create figure 3.41 p.129 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=Y(101:200,:);
[out]=FSM(Y);

%% SB (forged notes): univariate boxplots
% Create figure 3.50 p.135 of ARC 2004
clearvars;close all;
load('swiss_banknotes');
Y=swiss_banknotes.data;
Y=zscore(Y(101:200,:));
cnames=swiss_banknotes.colnames;
% Compare the output with Figure 3.9 p. 97
boxplot(Y,'labels',cnames);
% Label the outliers with the unit number
%Find in current plot the handles associated with the univariate outliers
a=findobj(gca,'tag','Outliers');
[n,v]=size(Y);
for j=1:v;
    % Get the X and Y coordinates of the univariate outliers
    aYdata=get(a(j),'Ydata');
    aXdata=get(a(j),'Xdata');
    
    % Loop over the outliers inside variable v-j+1
    for i=1:length(aYdata);
        % Find the row number of the outliers
        ind=find(aYdata(i)==Y(:,end-j+1),n);
        % Add the label to the current plot
        text(aXdata(1),aYdata(i),num2str(100+ind'))
    end
end



%% C2: two overlapping clusters
% This dataset has been used in Atkinson and Riani (2007)
clearvars;close all;
Y=load('clus2over.txt');
n=size(Y,1);
group=ones(n,1);
group(501:end)=2;
spmplot(Y,group);
% compare the output with Figure 4 of Atkinson and Riani (2007)

%%  C2:  analysis using S and MM estimators
clearvars;close all;
Y=load('clus2over.txt');
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')




%% C2: random starts
clearvars;close all;
Y=load('clus2over.txt');
[n,v]=size(Y);

init=20;
nsimul=100;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));

% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');
% compare the output with Figure 5 of Atkinson and Riani (2007)

%% OF: Old Faithful data
% Plot of the data
% This dataset has been used in Atkinson and Riani (2007)
clearvars;close all;
Y = load('geyser.txt');
plot(Y(:,1),Y(:,2),'o')

%% OF: Old Faithful data
clearvars;close all;
Y = load('geyser.txt');
plot(Y(:,1),Y(:,2),'o')

[n,v]=size(Y);

init=20;
nsimul=100;

mmdStore=nan(n-init,nsimul);
ij=0;
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    if ~isnan(mmd)
        mmdStore(:,j)=mmd(:,2);
    else
        ij=ij+1;
    end
end
disp(['Number of random starts with no full rank subsets =' num2str(ij)])
% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
% REMARK: in the instruction below it is possible to use instruction
% plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
% provided mmd is not NaN (that is provided that the last value of the random 
% starts loop produces a no full rank subset)
plot1=plot((init:n-1)',mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));

% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');
% compare the output with Figure 5 of Atkinson and Riani (2007)

%% OF: analysis using S and MM estimators
clearvars;close all;
Y = load('geyser.txt');
plot(Y(:,1),Y(:,2),'o')
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')



%% 6080 data
% Interactive_example
clearvars;close all;
Y=load('sixty_eighty.txt');
[fre]=unibiv(Y,'plots',1,'textlab',1,'rf',0.5);
fre=sortrows(fre,4);
init=20;
bs=fre(1:init,1);
% Forward search with EDA purposes
[out]=FSMeda(Y,bs,'plots',0,'init',init,'scaled',1);
malfwdplot(out,'databrush',1);

%% 6080 data (random starts)
clearvars;close all;
Y=load('sixty_eighty.txt');

[n,v]=size(Y);

init=20;
nsimul=100;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));


% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');

%% 6080 data: analysis using S and MM estimators
clearvars;close all;
Y=load('sixty_eighty.txt');
plot(Y(:,1),Y(:,2),'o')
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')



%% 3C 3 clusters 2 outliers
% Interactive_example
clearvars;close all;
Y=load('three_clust_2outl.txt');
[fre]=unibiv(Y,'plots',0,'textlab',1,'rf',0.5);
fre=sortrows(fre,4);
init=20;
bs=fre(1:init,1);
% Forward search with EDA purposes
[out]=FSMeda(Y,bs,'plots',0,'init',init,'scaled',1);
databrush=struct;
databrush.persist='on';
malfwdplot(out,'databrush',databrush);

%% 3C 3 clusters 2 outliers (random starts)
% Plot of the data
% This dataset has been used in Atkinson and Riani (2007)
clearvars;close all;
Y=load('three_clust_2outl.txt');
% plot(Y(:,1),Y(:,2),'o')

[n,v]=size(Y);

init=20;
nsimul=300;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));


% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');

%% BD bridge data
% Interactive_example
clearvars;close all;
Y=load('databri.txt');
[fre]=unibiv(Y,'plots',0,'textlab',1,'rf',0.5);
fre=sortrows(fre,4);
init=20;
bs=fre(1:init,1);
% Forward search with EDA purposes
[out]=FSMeda(Y,bs,'plots',0,'init',init,'scaled',1);
databrush=struct;
databrush.persist='on';
malfwdplot(out,'databrush',databrush);

%% BD bridge data (random starts)
% Plot of the data
clearvars;close all;
Y=load('databri.txt');
% plot(Y(:,1),Y(:,2),'o')

[n,v]=size(Y);

init=20;
nsimul=300;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));

% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');


%% FD Financial data
% Interactive_example
clearvars;close all;
Y=load('fondi.txt');
[fre]=unibiv(Y,'plots',0,'textlab',1,'rf',0.5);
fre=sortrows(fre,4);
init=20;
bs=fre(1:init,1);
% Forward search with EDA purposes
[out]=FSMeda(Y,bs,'plots',0,'init',init,'scaled',1);
databrush=struct;
databrush.persist='on';
malfwdplot(out,'databrush',databrush);

%% FD Financial data (random starts)
% Plot of the data
clearvars;close all;
Y=load('fondi.txt');
% plot(Y(:,1),Y(:,2),'o')

[n,v]=size(Y);

init=20;
nsimul=300;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);


set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));


% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');


%% FD: analysis using S and MM estimators
clearvars;close all;
% Plot of the data
Y=load('fondi.txt');
plot(Y(:,1),Y(:,2),'o')
[n,v]=size(Y);

conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')



%% DD Diabetes data
% Interactive_example
clearvars;close all;
Y=load('diabetes.txt');
[fre]=unibiv(Y,'plots',0,'textlab',1,'rf',0.5);
fre=sortrows(fre,4);
init=20;
bs=fre(1:init,1);
% Forward search with EDA purposes
[out]=FSMeda(Y,bs,'plots',0,'init',init,'scaled',1);
databrush=struct;
databrush.persist='on';
malfwdplot(out,'databrush',databrush);

%% DD Diabetes data (random starts)
clearvars;close all;
% Plot of the data
Y=load('diabetes.txt');

[n,v]=size(Y);

init=20;
nsimul=300;

mmdStore=NaN(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    if ~isnan(mmd)
        mmdStore(:,j)=mmd(:,2);
    end
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));

% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m');

%% MU Mussels data (Untransformed) -- Forward EDA (Exploratory Data Analysis):
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
v=size(Y,2);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',15,'scaled',1);
standard=struct;
standard.ylim=[0 8];
malfwdplot(out,'standard',standard);

%% MU Mussels data (Untransformed) -- Forward EDA (Exploratory Data Analysis):
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
v=size(Y,2);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',15);
standard=struct;
standard.ylim=[0 8];
mmdplot(out);

%% MU: automatic outlier detection procedure on original data
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
FSM(Y)

%% MU Mussels data analysis of transformations
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
% FS based on with H_0:\lambda=[1 0.5 1 0 1/3]
% Compare plot of mle with Figure 4.21 p. 178 of ARC (2004)
% Compare plot of lrt with Figure 4.20 p. 178 of ARC (2004)
[out]=FSMtra(Y,'la0',[0.5 0 0.5 0 0],'plotsmle',1,'plotslrt',1);


%% MU transformed: plot of MD
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
la=[0.5 0 0.5 0 0];
v=size(Y,2);
Y=normBoxCox(Y,1:v,la);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',15,'scaled',1);
malfwdplot(out);


%% MU transformed: plot of mmd
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
la=[0.5 0 0.5 0 0];
v=size(Y,2);
Y=normBoxCox(Y,1:v,la);
[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:size(Y,2)+5,1);
[out]=FSMeda(Y,bs,'init',15);
mmdplot(out);



%% MU comparison of mmd for untransformed and transformed data
clearvars;close all;
close all;
load('mussels.mat');
Y=mussels.data;
la=[0.5 0 0.5 0 0];
v=size(Y,2);
Y1=normBoxCox(Y,1:v,la);
[fre]=unibiv(Y1);
fre=sortrows(fre,[3 4]);
bs=fre(1:v+5,1);
[out1]=FSMeda(Y1,bs,'init',15);
% subplot(2,1,1);
mmdplot(out1,'tag','pl1');

[fre]=unibiv(Y);
fre=sortrows(fre,[3 4]);
bs=fre(1:v+5,1);
[out]=FSMeda(Y,bs,'init',15);
% subplot(2,1,2);
mmdplot(out);


%% MU transformed: automatic outlier detection procedure
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
la=[0.5 0 0.5 0 0];
v=size(Y,2);
Y=normBoxCox(Y,1:v,la);
FSM(Y,'plots',2)

%% MU transformed with la=[1/3 1/3 1/3 0 0] automatic outlier detection procedure
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
la=[1/3 1/3 1/3 0 0];
v=size(Y,2);
Y=normBoxCox(Y,1:v,la);
FSM(Y,'plots',2)

%% MU transformed random starts
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
la=[0.5 0 0.5 0 0];
[n,v]=size(Y);
Y=normBoxCox(Y,1:v,la);

init=20;
nsimul=300;

mmdStore=zeros(n-init,nsimul);
for j=1:nsimul
    mmd = FSMmmd(Y,0,'init',init);
    mmdStore(:,j)=mmd(:,2);
end

% Plot minMD with random starts
figure;
hold('on');
% Plot lines of empirical quantiles
LineStyle={'-';'--';':';'-.'};
plot1=plot(mmd(:,1),mmdStore,'LineWidth',2);
slintyp=repmat(LineStyle,ceil(nsimul/length(LineStyle)),1);
fcol={'b';'g';'r';'c';'m';'y';'k'};
fcol=repmat(fcol,ceil(nsimul/length(fcol)),1);

set(plot1,{'LineStyle'},slintyp(1:nsimul));
set(plot1,{'Color'},fcol(1:nsimul));

% Plots lines of theoretical quantiles using order statistics
mmdT=FSMenvmmd(n,v,'exact',1,'init',init);
line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
xlabel('Subset size m')



%% MU transformed: S and MM estimators
clearvars;close all;
load('mussels.mat');
Y=mussels.data;
la=[1/3 1/3 1/3 0 0];
[n,v]=size(Y);
Y=normBoxCox(Y,1:v,la);


conflev=[0.95 0.99 1-0.01/n];

[outS]=Smult(Y);
figure;
h1=subplot(2,1,1);
malindexplot(outS,v,'h',h1,'conflev',conflev);
ylabel('S estimator')

h2=subplot(2,1,2);

[outMM]=MMmultcore(Y,outS.loc,outS.shape,outS.scale);
malindexplot(outMM,v,'h',h2,'conflev',conflev);
ylabel('MM estimator')
