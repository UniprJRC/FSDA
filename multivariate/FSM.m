function [out]=FSM(Y,varargin)
%FSM gives an automatic outlier detection procedure in multivariate analysis
%
%<a href="matlab: docsearchFS('FSM')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Input data. Matrix.
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
% Optional input arguments:
%
%          m0   : Initial subset size or vector which contains the list of the units forming
%                 initial subset. Scalar or vector.
%                 The default is to start the search with v+1 units which
%                 consisting of those observations which are not outlying
%                 on any scatterplot, found as the intersection of all
%                 points lying within a robust contour containing a
%                 specified portion of the data (Riani and Zani 1997) and
%                 inside the univariate boxplot. Remark: if m0 is a vector
%                 option below crit is ignored.
%                 Example - 'm0',5
%                 Data Types - double
%       crit    : It specified the criterion to be used to
%                 initialize the search. Character.
%                 if crit='md' the units which form initial subset are
%                  those which have the smallest m0 pseudo Mahalanobis
%                  distances computed using procedure unibiv (bivariate
%                  robust ellipses).
%                 if crit='biv' sorting is done first in
%                  terms of times units fell outside robust bivariate
%                  ellipses and then in terms of pseudoMD. In other words,
%                  the units forming initial subset are chosen first among
%                  the set of those which never fell outside robust
%                  bivariate ellipses then among those which fell only once
%                  outside bivariate ellipses ... up to reach m0.
%                 if crit='uni' sorting is done first in
%                  terms of times units fell outside univariate boxplots
%                  and then in terms of pseudoMD. In other words,
%                  the units forming initial subset are chosen first among
%                  the set of those which never fell outside
%                  univariate boxplots then among those which fell only once
%                  outside univariate boxplots... up to reach m0.
%               Example - 'crit','md'
%               Data Types - char
%                 Remark: as the user can see the starting point of the
%                 search is not going to affect at all the results of the
%                 analysis. The user can explore this point with his own
%                 datasets.
%                 Remark: if crit='biv' the user can also supply in scalar rf
%                 (see below) the confidence level of the bivariate
%                 ellipses.
%        rf     : confidence level for bivariate ellipses. Scalar. The default is
%                 0.95. This option is useful only if crit='biv'.
%                 Example - 'rf',0.9
%                 Data Types - double
%       init    : Point where to start monitoring required diagnostics. Scalar.
%                 Note that if bsb is suppliedinit>=length(bsb). If init is not
%                 specified it will be set equal to floor(n*0.6).
%                 Example - 'init',50
%                 Data Types - double
%       plots   : plot of minimum Mahalanobis distance.
%                  Scalar or structure. If plots is a missing value or is a scalar equal to 0 no
%                 plot is produced.
%                 If plots is a scalar equal to 1 (default) the plot
%                 of minimum MD with envelopes based on n observations and
%                 the scatterplot matrix with the outliers highlighted is
%                 produced.
%                 If plots is a scalar equal to 2 the additional plots of
%                 envelope resuperimposition are
%                 produced.
%                 If plots is a structure it may contain the following fields:
%                   plots.ylim = vector with two elements controlling minimum and maximum
%                       on the y axis. Default value is '' (automatic
%                       scale);
%                   plots.xlim = vector with two elements controlling minimum and maximum
%                       on the x axis. Default value is '' (automatic
%                       scale);
%                   plots.resuper = vector which specifies for which steps it is
%                       necessary to show the plots of resuperimposed envelopes
%                       if resuper is not supplied a plot of each step in which the
%                       envelope is resuperimposed is shown. Example if resuper =[85 87]
%                       plots of resuperimposedenvelopes are shown at steps
%                       m=85 and m=87;
%                   plots.ncoord = scalar. If ncoord=1 plots are shown in normal
%                       coordinates else (default) plots are shown in
%                       traditional mmd coordinates;
%                   plots.labeladd = If this option is '1', the outliers in the
%                       spm are labelled with the unit row index. The
%                       default value is labeladd='', i.e. no label is
%                       added;
%                   plots.nameY = cell array of strings containing the labels of
%                       the variables. As default value, the labels which are
%                       added are Y1, ...Yv;
%                   plots.lwd =  Scalar which controls line width of the curve which
%                       contains the monitoring of minimum Mahalanobis
%                       distance. Default line of lwd=2.
%                   plots.lwdenv = Scalar which controls linewidth of the
%                       envelopes. Default value of lwdenv=2.
%               Example - 'plots',2
%               Data Types - double
%      bonflev  : option that might be used to identify extreme outliers
%                 when the distribution of the data is strongly non normal.
%                 Scalar.
%                 In these circumstances, the general signal detection rule based on
%                 consecutive exceedances cannot be used. In this case
%                 bonflev can be:
%                 - a scalar smaller than 1, which specifies the confidence
%                   level for a signal and a stopping rule based on the
%                   comparison of the minimum deletion residual with a
%                   Bonferroni bound. For example if bonflev=0.99 the
%                   procedure stops when the trajectory exceeds for the
%                   first time the 99% bonferroni bound.
%                 - A scalar value greater than 1. In this case the
%                   procedure stops when the residual trajectory exceeds
%                   for the first time this value.
%                 Default value is ' ', which means to rely on general rules
%                 based on consecutive exceedances.
%                 Example - 'bonflev',0.7
%                 Data Types - double
%       msg     : It controls whether to display or not messages
%                 on the screen. Scalar.
%                 If msg==1 (default) messages about the progression of the
%                 search are displayed on the screen otherwise only error
%                 messages will be displayed.
%                 Example - 'msg',0
%                 Data Types - double
%   nocheck     : It controls whether to perform checks on matrix Y.Scalar.
%                 If nocheck is equal to 1 no check is performed.
%                 As default nocheck=0.
%                 Example - 'nocheck',1
%                 Data Types - double
%
%
% Output:
%
%         out:   structure which contains the following fields
%
%out.outliers=  k x 1 vector containing the list of the units declared as
%               outliers or NaN if the sample is homogeneous
% out.mmd    =  (n-init) x 2 matrix.
%               1st col = fwd search index;
%               2nd col = value of minimum Mahalanobis Distance in each step
%               of the fwd search.
% out.Un     =  (n-init) x 11 Matrix which contains the unit(s) included
%               in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one. Un(1,2) for example
%               contains the unit included in step init+1. Un(end,2)
%               contains the units included in the final step of the search.
% out.nout    = 2 x 5 matrix containing the number of times mmd went out
%               of particular quantiles.
%               First row contains quantiles 1 99 99.9 99.99 99.999 per cent;
%               Second row contains the frequency distribution. It is NaN
%               if bonflev threshold is used.
% out.loc     = 1 x v  vector containing location of the data.
% out.cov     = v x v robust estimate of covariance matrix.
% out.md      = n x 1 vector containing the estimates of the robust
%               Mahalanobis distances (in squared units). This vector
%               contains the distances of each observation from the
%               location of the data, relative to the scatter matrix cov.
% out.class  =  'FSM'.
%
%
% See also: FSMeda, unibiv.m, FSMmmd.m
%
% References:
%
% Riani, M., Atkinson, A.C. and Cerioli, A. (2009), Finding an unknown
% number of multivariate outliers, "Journal of the Royal Statistical
% Society Series B", Vol. 71, pp. 201-221.
% Cerioli, A., Farcomeni, A. and Riani M. (2014), Strong consistency and
% robustness of the Forward Search estimator of multivariate location
% and scatter, "Journal of Multivariate Analysis", Vol. 126,
% pp. 167-183, http://dx.doi.org/10.1016/j.jmva.2013.12.010
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('FSM')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% FSM with all default options.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y; Ycont(1:5,[1,3]) = Ycont(1:5,[1,3])+sign(randn(5,2))*4.5;
    [out]=FSM(Ycont);
    title('Outliers detected by FSM','Fontsize',24,'Interpreter','LaTex');
%}

%{
    %% FSM with optional arguments.
    % FSM with plots showing envelope superimposition.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=FSM(Ycont,'plots',2);
%}

%{
    %% FSM with plots showing envelope superimposition in normal
    % coordinates.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    plots=struct;
    plots.ncoord=1;
    [out]=FSM(Ycont,'plots',plots);
%}


%{
    % Monitor the exceedances from m=200 without showing plots.
    n=1000;
    v=10;
    Y=randn(n,v);
    [out]=FSM(Y,'init',200,'plots',0);
%}

%{
    % Choosing an initial subset formed by the three observations with
    % the smallest Mahalanobis Distance.
    n=100;
    v=3;
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=FSM(Ycont,'m0',5,'crit','md');
%}

%{
    % Forgery Swiss banknotes examples.

    load('swiss_banknotes');

    % Monitor the exceedances of Minimum Mahalanobis Distance
    [out]=FSM(swiss_banknotes.data(101:200,:),'plots',1);

    % Control minimum and maximum on the x axis
    plots=struct;
    plots.xlim=[60 90];
    [out]=FSM(swiss_banknotes.data(101:200,:),'plots',plots);

    % Monitor the exceedances of Minimum Mahalanobis Distance using normal coordinates for mmd.
    plots.ncoord=1;
    [out]=FSM(swiss_banknotes.data(101:200,:),'plots',plots);
%}

%% Input parameters checking
%chkinputM does not do any check if option nocheck=1
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);

[n,v]=size(Y);
seq=(1:n)';

hdef=floor(n*0.6);

if v>1
    critdef='md';
else
    critdef ='uni';
end

options=struct('m0',v+1,'init',hdef,'crit',critdef,'rf',0.95,...
    'plots',1,'msg',1,'bonflev','','nocheck',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
plo=options.plots;


msg=options.msg;
crit=options.crit;
m0=options.m0;

% fsizeannot is a scalar which Font Size of the annotations which are
% shown on the screen
fsizeannot=11;

%% Start of the forward search

if length(m0)>1
    bs=m0;
else
    % Confidence level for robust bivariate ellipses
    rf=options.rf;
    
    % Find initial subset to initialize the search
    [fre]=unibiv(Y,'rf',rf);
    
    if strcmp(crit,'md')==1
        % The user has chosen to select the intial subset according to the
        % smallest m0 pseudo MD Select only the potential bivariate outliers
        fre=sortrows(fre,4);
    elseif strcmp(crit,'biv')==1
        fre=sortrows(fre,[3 4]);
    elseif strcmp(crit,'uni')==1
        fre=sortrows(fre,[2 4]);
    else
        error('FSDA:FSM:WrongInputOpt','Supplied options to initialize the search does not exist. crit must be ''md'' ''biv'' or ''uni''');
    end
    
    % initial subset
    bs=fre(1:m0,1);
    
    % the subset need to be incremented if it is not full rank. We also
    % treat the unfortunate case when the rank of the matrix is v but a
    % column is constant.
    incre = 1;
    %the second condition is added to treat subset with a constant
    %variable. This situation does not decrease the rank of Y, but it
    %decreases the rank of ym (i.e. Y-mean(Y)) inside FSMmmd.
    while (rank(Y(bs,:))<v) || min(max(Y(bs,:)) - min(Y(bs,:))) == 0
        bs=fre(1:m0+incre,1);
        incre = incre+1;
    end
    
    % To make sure that new value of init is minimum lenght of bs for which
    % the Y matrix is full rank
    if init<length(bs)
        init=length(bs);
    end
    
end



% Compute Minimum Mahalanobis Distance for each step of the search
if n<5000
    [mmd,Un,bb] = FSMmmd(Y,bs,'init',init,'nocheck',1,'msg',msg);
else
    [mmd,Un] = FSMmmd(Y,bs,'init',init,'nocheck',1,'msg',msg);
end

if isnan(mmd)
    out = nan;
    return
end


bonflev=options.bonflev;



%% Part 1. Signal detection and validation
signal=0;
sto=0;
extram3='';
extram2='';
strplot='';
resup=2;
if msg
    disp('-------------------------')
    disp('Signal detection loop');
end

nmmd=size(mmd,1);


if ~isempty(bonflev)
    if bonflev<1
        [gbonf] = FSMbonfbound(n,v,'prob',bonflev,'init',init);
        bonfthresh=gbonf;
    else
        bonfthresh=bonflev*ones(n-init,1);
    end
else
    if nmmd<4
        error('FSDA:FSM:WrongRationv','ratio n/v too small; modify init (i.e. decrease initial subset size)')
    end
    
    
    quant=[0.99;0.999;0.9999;0.99999;0.01;0.5];
    % Compute theoretical envelops for minimum Mahalanobis Distance based on all
    % the observations for the above quantiles.
    [gmin] = FSMenvmmd(n,v,'prob',quant,'init',init);
    % gmin = the matrix which contains envelopes based on all observations.
    % 1st col of gmin = fwd search index
    % 2nd col of gmin = 99% envelope
    % 3rd col of gmin = 99.9% envelope
    % 4th col of gmin = 99.99% envelope
    % 5th col of gmin = 99.999% envelope
    % 6th col of gmin = 1% envelope
    % 7th col of gmin = 50% envelope
    % Thus, set the columns of gmin where the theoretical quantiles are located.
    [c99 , c999 , c9999 , c99999 , c001 , c50] = deal(2,3,4,5,6,7);
    
    % Store in nout the number of times the observed mmd (d_min) lies above:
    [out99 , out999 , out9999 , out99999 , out001] = deal( ...
        mmd(mmd(:,2)>gmin(:,c99),:) , ...       % the 99% envelope
        mmd(mmd(:,2)>gmin(:,c999),:) , ...      % the 99.9% envelope
        mmd(mmd(:,2)>gmin(:,c9999),:) , ...     % the 99.99% envelope
        mmd(mmd(:,2)>gmin(:,c99999),:) , ...    % the 99.999% envelope
        mmd(mmd(:,2)<gmin(:,c001),:) );         % the 1% envelope
    
    nout = [[1 99 999 9999 99999]; ...
        [size(out001,1) size(out99,1) size(out999,1) size(out9999,1) size(out99999,1)]];
    
    % NoFalseSig = boolean linked to the fact that the signal is good or not
    NoFalseSig=0;
    
    % NoFalseSig is set to 1 if the condition for an INCONTROVERTIBLE SIGNAL is
    % fulfilled.
    n9999 = nout(2,nout(1,:)==9999);
    if (n9999>=10)
        NoFalseSig=1;
        if msg
            disp('Observed curve of d_min is at least 10 times greater than 99.99% envelope'); % exact number is int2str(n9999)
            disp('--------------------------------------------------');
        end
    end
    
    % Divide central part from final part of the search
    istep = n-floor(13*sqrt(n/200));
end


%% Stage 1a: signal detection
% Signal detection is based on monitoring consecutive triplets or single
% extreme values
% Signal detection loop
for i=3:nmmd
    
    if isempty(bonflev)
        
        if i<istep-init+1 % CENTRAL PART OF THE SEARCH
            % Extreme triplet or an extreme single value
            % Three consecutive values of d_min above the 99.99% threshold or 1
            % above 99.999% envelope
            if ((mmd(i,2)>gmin(i,c9999) && mmd(i+1,2)>gmin(i+1,c9999) && mmd(i-1,2)>gmin(i-1,c9999)) || mmd(i,2)>gmin(end,c99) || mmd(i,2)>gmin(i,c99999))
                if msg
                    disp(['Tentative signal in central part of the search: step m=' int2str(mmd(i,1)) ' because']);
                end
                if (mmd(i,2)>gmin(i,c9999) && mmd(i+1,2)>gmin(i+1,c9999) && mmd(i-1,2)>gmin(i-1,c9999))
                    if msg
                        disp(['dmin('  int2str(mmd(i,1)) ',' int2str(n) ')>99.99% and dmin(' int2str(mmd(i-1,1)) ',' int2str(n) ')>99.99% and rmin(' int2str(mmd(i+1,1)) ',' int2str(n) ')>99.99%']);
                    end
                    strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99.99\%$ and $d_{min}(' int2str(mmd(i-1,1)) ',' int2str(n) ')>99.99\%$ and $d_{min}(' int2str(mmd(i+1,1)) ',' int2str(n) ')>99.99\%$'];
                    mmdsel=mmd(i-1:i+1,1:2);
                end
                
                if (mmd(i,2)>gmin(i,c99999))
                    if msg
                        disp(['dmin(' int2str(mmd(i,1)) ',' int2str(n) ')>99.999%']);
                    end
                    strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99.999\%$'];
                    mmdsel=mmd(i-1:i+1,1:2);
                end
                
                if (mmd(i,2)>gmin(end,c99))
                    if msg
                        disp(['dmin(' int2str(mmd(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the central part of the search.']);
                    end
                    strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
                    mmdsel=mmd(i:i,1:2);
                    NoFalseSig=1; % i.e., no need of further validation
                end
                '------------------------------------------------';
                
                signal=1;
            end
        elseif i<size(mmd,1)-1 % FINAL PART OF THE SEARCH
            % Extreme couple adjacent to an exceedance
            % Two consecutive values of mmd above the 99.99% envelope and 1 above 99%
            if ((mmd(i,2)>gmin(i,c999) && mmd(i+1,2)>gmin(i+1,c999) && mmd(i-1,2)>gmin(i-1,c99)) || (mmd(i-1,2)>gmin(i-1,c999) && mmd(i,2)>gmin(i,c999) && mmd(i+1,2)>gmin(i+1,c99)) || mmd(i,2)>gmin(end,c99) || mmd(i,2)>gmin(i,c99999))
                'Signal in final part of the search: step '; mmd(i,1); 'because';
                if (mmd(i,2)>gmin(i,c999) && mmd(i+1,2)>gmin(i+1,c999) && mmd(i-1,2)>gmin(i-1,c99))
                    if msg
                        disp(['dmin('  int2str(mmd(i,1)) ',' int2str(n) ')>99.9% and dmin('  int2str(mmd(i+1,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mmd(i-1,1)) ',' int2str(n) ')>99%']);
                    end
                    strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99.9\%$ and $d_{min}(' int2str(mmd(i-1,1)) ',' int2str(n) ')>99\%$ and $d_{min}(' int2str(mmd(i+1,1)) ',' int2str(n) ')>99.9\%$'];
                    mmdsel=mmd(i-1:i+1,1:2);
                end
                
                if (mmd(i-1,2)>gmin(i-1,c999) && mmd(i,2)>gmin(i,c999) && mmd(i+1,2)>gmin(i+1,c99))
                    if msg
                        disp(['drmin('  int2str(mmd(i-1,1)) ',' int2str(n) ')>99.9% and dmin('  int2str(mmd(i,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mmd(i+1,1)) ',' int2str(n) ')>99%']);
                    end
                    strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99.9\%$ and $d_{min}(' int2str(mmd(i-1,1)) ',' int2str(n) ')>99.9\%$ and $d_{min}(' int2str(mmd(i+1,1)) ',' int2str(n) ')>99\%$'];
                    mmdsel=mmd(i-1:i+1,1:2);
                end
                
                if (mmd(i,2)>gmin(end,c99))
                    if msg
                        disp(['dmin('  int2str(mmd(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the final part of the search.']);
                    end
                    strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
                    mmdsel=mmd(i:i,1:2);
                end
                
                % Extreme single value
                if mmd(i,2)>gmin(i,c99999)
                    if msg
                        disp(['dmin('  int2str(mmd(i,1)) ',' int2str(n) ')>99.999%']);
                    end
                    strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99.999\%$'];
                    mmdsel=mmd(i:i,1:2);
                end
                
                '------------------------------------------------';
                % Signal is always considered true if it takes place in the
                % final part of the search
                NoFalseSig=1;
                signal=1;
            end
        elseif (mmd(i,2)>gmin(i,c999) || mmd(i,2)>gmin(end,c99)) && i==size(mmd,1)-1
            % potential couple of outliers
            signal=1;
            if msg
                disp('Signal is in penultimate step of the search');
            end
            
            if (mmd(i,2)>gmin(i,c999))
                if msg
                    disp(['dmin(' int2str(mmd(i,1)) ',' int2str(n) ')>99.9%']);
                end
                strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99.9\%$'];
            end
            
            if (mmd(i,2)>gmin(end,c99))
                if msg
                    disp(['dmin('  int2str(mmd(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the final part of the search.']);
                end
                strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
            end
            mmdsel=mmd(i:i,1:2);
        elseif  mmd(i,2)>gmin(i,c99) && i==size(mmd,1)
            % a single outlier
            signal=1;
            if msg
                disp('Signal is in final step of the search');
            end
            strplot=['$d_{min}(' int2str(mmd(i,1)) ',' int2str(n) ')>99\%$ at final step'];
            mmdsel=mmd(i:i,1:2);
        end
        
        %% Stage 1b: signal validation
        if (signal==1)
            if msg
                disp('-------------------')
                disp('Signal validation');
            end
            % mdag is $m^\dagger$
            mdag=mmd(i,1);
            
            if mmd(i,1)<n-2
                % Check if the signal is incontrovertible
                % Incontrovertible signal = 3 consecutive values of d_min >
                % 99.999% threshold
                if mmd(i,2)>gmin(i,c99999) && mmd(i-1,2)>gmin(i-1,c99999) &&  mmd(i+1,2)>gmin(i+1,c99999)
                    if msg
                        disp(['3 consecutive values of d_min greater than 99.999% envelope in step mdag= ' int2str(mmd(i,1))]);
                    end
                    NoFalseSig=1;
                    extram3='Extreme signal';
                end
            else
                NoFalseSig=1;
            end
            
            % if the following statement is true, observed curve of d_min is
            % above 99.99% and later is below 1%: peak followed by dip
            if size(mmd,1)>mdag-mmd(1,1)+31
                if sum(mmd(i+1:i+31,2)<gmin(i+1:i+31,c001))>=2
                    NoFalseSig=1;  % Peak followed by dip
                    extram2='Peak followed by dip (d_min is above 99.99% threshold and in the sucessive 30 steps goes below 1% envelope';
                end
            else
                if sum(mmd(i+1:end,2) < gmin(i+1:end,c001))>=2
                    NoFalseSig=1;  %Peak followed by dip in the final part of the search';
                    extram2='Peak followed by dip (d_min is above 99.99% threshold and in the sucessive 30 steps goes below 1% envelope)';
                end
            end
            
            % if at this point NoFalseSig==0 it means that:
            % 1) n9999<10
            % 2) signal tool place in the central part of the search
            % 3) signal was not incontrovertible
            % 4) there was not a peak followed by dip
            if NoFalseSig==0
                % Compute the final value of the envelope based on
                % mmd(i+1,1)=mdagger+1 observations
                [gval]=FSMenvmmd(mdag+1,v,'prob',0.01,'init',mdag);
                if mmd(i,2)<gval(1,2)
                    if msg
                        disp('false signal in step');
                        disp(['mdag='  int2str(mdag)]);
                    end
                    % increase mdag of the search by one unit
                    mdag=0;
                else
                    NoFalseSig=1;
                end
            end
            
            % If the signal has been validated get out of the signal detection
            % loop and move to stage 2: superimposition of the envelopes
            if (NoFalseSig==1)
                if msg
                    disp('Validated signal');
                end
                break;
            end
        end
    else
        % Outlier detection based on Bonferroni threshold
        if (mmd(i,2)>bonfthresh(i,end))
            if msg
                disp(['$d_min$(' int2str(mmd(i,1)) ',' int2str(n) ')>99% Bonferroni level']);
            end
            strplot=['$d_min$(' int2str(mmd(i,1)) ',' int2str(n) ')>99\%$ (Bonferroni level)'];
            mmdsel=mmd(i:i,1:2);
            
            signal=1;
            break
        end
    end
end

%% Create figure containing mmd + envelopes based on all the observations.
% if plo is a structure or if it is a scalar different from 0
% then produce a plot
if isstruct(plo) || (~isstruct(plo) && plo~=0)
    
    % get screen size [left, bottom, width, height]
    scrsz = get(0,'ScreenSize');
    
    figure1 = figure('Position',[1 scrsz(4)/2.5 scrsz(3)/3 scrsz(4)/2],'PaperSize',[20.98 29.68],'Name','Envelopes based on all the observations');
    axes1 = axes('Parent',figure1);
    
    % Create a pan-handle for the figure ...
    hpan_figure1 = pan(figure1);
    % ... and listen to pan events using callback functions
    set(hpan_figure1,'ActionPreCallback',@figure1_precallback);
    set(hpan_figure1,'ActionPostCallback',@figure1_postcallback);
    
    if isstruct(plo)
        
        fplo=fieldnames(plo);
        
        d=find(strcmp('xlim',fplo));
        if d>0
            xlimx=plo.xlim;
        else
            xlimx='';
        end
        
        d=find(strcmp('ylim',fplo));
        if d>0
            ylimy=plo.ylim;
        else
            ylimy='';
        end
        
        
        d=find(strcmp('ncoord',fplo));
        if d>0
            ncoord=plo.ncoord;
        else
            ncoord=0;
        end
        
        d=find(strcmp('lwd',fplo));
        if d>0
            lwd=plo.lwd;
        else
            lwd=2;
        end
        
        d=find(strcmp('lwdenv',fplo));
        if d>0
            lwdenv=plo.lwdenv;
        else
            lwdenv=2;
        end
        
        
    else
        xlimx='';
        ylimy='';
        ncoord=0;
        lwd=2;
        lwdenv=2;
    end
    
    if isempty(xlimx)
        xl1=init-3; xl2=mmd(end,1);
    else
        xl1=xlimx(1);
        xl2=xlimx(2);
    end
    
    
    % set the xlimits
    % ylimits are set later according to option ncoord (plot in normal
    % coordinates)
    xlim([xl1 xl2]);
    
    kx=0; ky=0;
    
    box('on'); hold('all');
    
    
    
    if isempty(bonflev)
        
        if ncoord ~=1
            
            % Set the ylimits in mdr coordinates
            if isempty(ylimy)
                yl1=min([gmin(:,c001);mmd(:,2)]);
                yl2=max([gmin(:,c999);mmd(:,2)]);
                
            else
                yl1=ylimy(1);
                yl2=ylimy(2);
            end
            
            ylim([yl1 yl2]);
            
            plot(mmd(:,1),mmd(:,2),'LineWidth',lwd);
            
            
            % Superimpose 1%, 99%, 99.9% envelopes based on all the observations
            line(gmin(:,1),gmin(:,[c001 c99 c999]),'Parent',axes1,'LineWidth',lwdenv,'LineStyle','--','Color',[0 0 1]);
            % Superimpose 99.99% and 99.999% envelopes based on all the observations
            line(gmin(:,1),gmin(:,[c9999 c99999]),'Parent',axes1,'LineWidth',lwdenv,'LineStyle','--','Color',[1 0 0]);
            % Superimpose 50% envelope based on all the observations
            line(gmin(:,1),gmin(:,c50),'Parent',axes1,'LineWidth',lwdenv,'LineStyle','--','Color',[1 0.69 0.39]);
            
            % Property-value pairs which are common to all quantile-labels
            PrVaCell{1,1} = 'HorizontalAlignment'; PrVaCell{2,1} = 'center';
            PrVaCell{1,2} = 'EdgeColor'; PrVaCell{2,2} = 'none';
            PrVaCell{1,3} = 'BackgroundColor'; PrVaCell{2,3} = 'none';
            PrVaCell{1,4} = 'FitBoxToText'; PrVaCell{2,4} = 'off';
            PrVaCell{1,5} = 'Tag'; PrVaCell{2,5} = 'quantile_label';
            
            % Create textbox with 1% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c001));
            if figy>=0 && figy<=1 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'1%'},...
                    'UserData',[gmin(:,1) gmin(:,c001)],...
                    PrVaCell{:},'FontSize',fsizeannot);
            end
            
            % Create textbox with 99% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c99));
            
            if figy>=0 && figy<=1 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'99%'},...
                    'UserData',[gmin(:,1) gmin(:,c99)],...
                    PrVaCell{:},'FontSize',fsizeannot);
            end
            
            % Create textbox with 50% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c50));
            
            if figy>=0 && figy<=1 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'50%'},...
                    'UserData',[gmin(:,1) gmin(:,c50)],...
                    PrVaCell{:},'FontSize',fsizeannot);
            end
            
            % Create textbox with 99.9% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c999));
            if figy>=0 && figy<=1 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'99.9%'},...
                    'UserData',[gmin(:,1) gmin(:,c999)],...
                    PrVaCell{:},'FontSize',fsizeannot);
            end
            
            % Create textbox with 99.99% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c9999));
            if figy<=1 && figy>=0 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'99.99%'},...
                    'UserData',[gmin(:,1) gmin(:,c9999)],...
                    PrVaCell{:},'FontSize',fsizeannot);
            end
            
            if gmin(1,c99999)<=yl2
                % Create textbox with 99.999% label
                [figx, figy] = dsxy2figxy(gca, init, gmin(1,c99999));
                if figy<=1 && figy>=0 && figx>=0 && figx<=1
                    annotation(figure1,'textbox',[figx figy kx ky],...
                        'String',{'99.999%'},...
                        'UserData',[gmin(:,1) gmin(:,c99999)],...
                        PrVaCell{:},'FontSize',fsizeannot);
                end
            end
            
            % Add string which informs about the step where signal took place
            if signal==1
                strsig=['Signal is in step $m=' int2str(mdag) '$ because'];
                stem(mmdsel(:,1),mmdsel(:,2),'LineWidth',1,...
                    'Color',[0.4784 0.06275 0.8941], 'DisplayName','Signal');
            else
                strsig='No signal during the search';
            end
            
            % Property-value pairs which are common to the next latex annotations
            PrVaCell{1,1} = 'Interpreter'; PrVaCell{2,1} = 'latex';
            PrVaCell{1,2} = 'HorizontalAlignment'; PrVaCell{2,2} = 'center';
            PrVaCell{1,3} = 'FitBoxToText'; PrVaCell{2,3} = 'on';
            PrVaCell{1,4} = 'EdgeColor'; PrVaCell{2,4} = 'none';
            PrVaCell{1,5} = 'BackgroundColor'; PrVaCell{2,5} = 'none';
            
            % latex annotations informing that the envelopes are based on
            % all the observations
            strmin=['$d_{min}(m,' int2str(n) ')$. '];
            annotation(figure1,'textbox',[0.5 0.90 kx ky],'String',{[strmin strsig]},...
                PrVaCell{:},'FontSize',fsizeannot);
            
            annotation(figure1,'textbox',[0.5 0.85 kx ky],'String',{strplot},...
                PrVaCell{:},'FontSize',fsizeannot);
            
            if istep<=xl2
                % Add vertical line which divides central part from final part of the
                % search
                [figx, figy] = dsxy2figxy(gca, istep, yl1);
                [figx, figy2] = dsxy2figxy(gca, istep, yl2);
                if figy2>=1
                    figy2=1;
                else
                end
                if figx>=0
                    annotation(figure1,'line',[figx figx],[figy figy2],...
                        'UserData',[istep yl1 yl2],...
                        'Tag','FinalPartLine')
                end
            end
        else % show the above plot in normal coordinates
            
            
            % Transform mmd in normal coordinates
            mmdinv = FSMinvmmd(mmd,v);
            
            % Set ylim in normal coordinates
            if isempty(ylimy)
                yl1=min([-2.33;mmdinv(:,3)]);
                yl2=max([3.09;mmdinv(:,3)]);
                
            else
                yl1=ylimy(1);
                yl2=ylimy(2);
            end
            
            
            % Plot in normal coordinates
            plot(mmdinv(:,1),mmdinv(:,3),'LineWidth',lwd);
            
            ylim([yl1 yl2]);
            
            % Add horizontal lines associated with confidence levels
            vaxis=axis;
            ninv=norminv(quant);
            for i=1:length(quant)
                if quant(i)==0.5
                    col=[1 0.69 0.39];
                elseif quant(i)>=0.01 &&  quant(i)<=0.999
                    col='b';
                else
                    col='r';
                end
                line(vaxis(1:2)',[ninv(i);ninv(i)],'color',col,'LineWidth',lwdenv,'LineStyle','--','Tag','env');
            end
            text(vaxis(1)*ones(length(ninv),1),ninv+0.2,strcat(num2str(100*quant),'%'),...
                'FontSize',12,'HorizontalAlignment','Left');
            
            
            % Add string which informs about the step where signal took place
            if signal==1
                strsig=['Signal is in step $m=' int2str(mdag) '$ because'];
                
                seli=find(mmdinv(:,1)>=mmdsel(1,1),size(mmdsel,1));
                invmmdsel=mmdinv(seli,:);
                stem(invmmdsel(:,1),invmmdsel(:,3),'LineWidth',1,...
                    'Color',[0.4784 0.06275 0.8941], 'DisplayName','Signal');
                
            else
                strsig='No signal during the search';
            end
            
            % Property-value pairs which are common to the next latex annotations
            PrVaCell{1,1} = 'Interpreter'; PrVaCell{2,1} = 'latex';
            PrVaCell{1,2} = 'HorizontalAlignment'; PrVaCell{2,2} = 'center';
            PrVaCell{1,3} = 'FitBoxToText'; PrVaCell{2,3} = 'on';
            PrVaCell{1,4} = 'EdgeColor'; PrVaCell{2,4} = 'none';
            PrVaCell{1,5} = 'BackgroundColor'; PrVaCell{2,5} = 'none';
            
            
            % latex annotations informing that the envelopes are based on
            % all the observations
            % ycoordinates of the messages displayed on the screen
            ycoordannott=0.90;
            ycoordannotb=0.85;
            
            strmin=['$d_{min}(m,' int2str(n) ')$. '];
            annotation(figure1,'textbox',[0.5 ycoordannott kx ky],'String',{[strmin strsig]},...
                PrVaCell{:},'FontSize',fsizeannot);
            
            annotation(figure1,'textbox',[0.5 ycoordannotb kx ky],'String',{strplot},...
                PrVaCell{:},'FontSize',fsizeannot);
            
            if istep<=xl2
                % Add vertical line which divides central part from final part of the
                % search
                [figx, figy] = dsxy2figxy(gca, istep, yl1);
                [figx, figy2] = dsxy2figxy(gca, istep, yl2);
                if figy2>=1
                    figy2=1;
                else
                end
                if figx>=0
                    annotation(figure1,'line',[figx figx],[figy figy2],...
                        'UserData',[istep yl1 yl2],...
                        'Tag','FinalPartLine');
                end
            end
            
        end
    else
        % Set the ylimits in mdr coordinates
        if isempty(ylimy)
            yl1=min(mmd(:,2));
            yl2=max([bonfthresh(:,2);mmd(:,2)]);
            
        else
            yl1=ylimy(1);
            yl2=ylimy(2);
        end
        
        ylim([yl1 yl2]);
        
        
        plot(mmd(:,1),mmd(:,2),'LineWidth',lwd);
        
        % Superimpose Bonferroni line to the plot
        line(bonfthresh(:,1),bonfthresh(:,end),'Parent',axes1,'LineWidth',lwdenv,'LineStyle','--','Color',[0 0 1]);
        % Property-value pairs which are common to the next latex annotations
        
        PrVaCell{1,1} = 'Interpreter'; PrVaCell{2,1} = 'latex';
        PrVaCell{1,2} = 'HorizontalAlignment'; PrVaCell{2,2} = 'center';
        PrVaCell{1,3} = 'FitBoxToText'; PrVaCell{2,3} = 'on';
        PrVaCell{1,4} = 'EdgeColor'; PrVaCell{2,4} = 'none';
        PrVaCell{1,5} = 'BackgroundColor'; PrVaCell{2,5} = 'none';
        
        
        if size(bonfthresh,2)>1
            % latex annotations informing that the envelopes are based on
            % all the observations
            strmin='Exceedance based on Bonferroni threshold';
            annotation(figure1,'textbox',[0.5 0.90 kx ky],'String',strmin,...
                PrVaCell{:});
            msgtxt=['$d_{min}(' num2str(mmd(i,1)) ',' int2str(n) ')>' num2str(100*bonflev) '$\% envelope'];
            annotation(figure1,'textbox',[0.5 0.80 kx ky],'String',msgtxt,PrVaCell{:},'FontSize',fsizeannot);
        else
            strmin='Exceedance based on user supplied threshold';
            annotation(figure1,'textbox',[0.5 0.90 kx ky],'String',strmin,...
                PrVaCell{:},'FontSize',fsizeannot);
            msgtxt=['$d_{min}(' num2str(mmd(i,1)) ',' int2str(n) ')>$' num2str(bonflev)];
            annotation(figure1,'textbox',[0.5 0.80 kx ky],'String',msgtxt,PrVaCell{:},'FontSize',fsizeannot);
        end
        if signal==1
            stem(mmd(i,1),mmd(i,2),'LineWidth',1,...
                'Color',[0.4784 0.06275 0.8941], 'DisplayName','Signal');
        end
    end
    
end

%% Part 2: envelope resuperimposition
% if a validated signal tool place, superimposition of the envelopes starts
% from m^\dagger-1

if (signal==1)
    if isempty(bonflev)
        if msg
            disp('-------------------------------');
            disp(['Start resuperimposing envelopes from step m=' int2str(mdag-1)]);
        end
        
        if isstruct(plo)
            d=find(strcmp('resuper',fplo));
            if d>0
                resuper=plo.resuper;
            else
                resuper=seq;
            end
        else
            if plo==2
                resuper=seq;
            else
                resuper='';
            end
        end
        
        
        if ~isempty(resuper)
            % jwind is associated with subplot window number
            % nr is the number of row panes in the plot
            % nc is the number of columns panes in the plot
            jwind=1;
            nc=2;
            if mmd(i,1)>=n-2
                nr=1;
            else
                nr=2;
            end
            
            figure2 = figure('PaperSize',[20.98 29.68],'Name','Resuperimposed envelopes #1');
            % Create axes
            axes('Parent',figure2);
        end
        
        % First resuperimposed envelope is based on mdag-1 observations
        % Notice that mmd(i,1) = m dagger
        for tr=(mdag-1):(n)
            % Compute theoretical envelopes based on tr observations
            gmin1=FSMenvmmd(tr,v,'prob',[0.99; 0.999; 0.01; 0.5],'init',init);
            
            for ii=(i-1):size(gmin1,1)
                
                % CHECK IF STOPPING RULE IS FULFILLED
                % ii>=size(gmin1,1)-2 = final, penultimate or antepenultimate value
                % of the resuperimposed envelope based on tr observations
                if mmd(ii,2)>gmin1(ii,c99) && ii>=size(gmin1,1)-2
                    % Condition S1
                    mes=['$d_{min}('   int2str(mmd(ii,1)) ',' int2str(tr) ')>99$\% envelope'];
                    if msg
                        disp(['Superimposition stopped because d_{min}(' int2str(mmd(ii,1)) ',' int2str(tr) ')>99% envelope']);
                        disp(mes);
                    end
                    sto=1;
                    break;
                elseif ii<size(gmin1,1)-2 &&  mmd(ii,2)>gmin1(ii,c999)
                    % Condition S2
                    mes=['$d_{min}('   int2str(mmd(ii,1)) ',' int2str(tr) ')>99.9$\% envelope'];
                    if msg
                        disp(['Superimposition stopped because d_{min}(' int2str(mmd(ii,1)) ',' int2str(tr) ')>99.9% envelope']);
                    end
                    sto=1;
                    break;
                else
                    % mmd is inside the envelopes, so keep resuperimposing
                end
            end
            
            if ~isempty(resuper) && ~isempty(intersect(resuper,tr))
                % Plotting part
                subplot(nr,nc,jwind,'Parent',figure2);
                ylim([yl1 yl2]); xlim([xl1 xl2]);
                box('on'); hold('on');
                
                if ncoord~=1  % superimposition in original coordinates
                    % Show curve of mmd up to step tr-1 (notice that the envelope is
                    % based on tr observations. Step tr-1 in matrix mmd is
                    % (tr-1)-mmd(1,1)+1=tr-mmd(1,1)
                    plot(mmd(1:(tr-mmd(1,1)),1),mmd(1:(tr-mmd(1,1)),2),'LineWidth',lwd);
                    
                    
                    % Display the lines associated with 1%, 50% 99% and 99.9% envelopes
                    line(gmin1(:,1),gmin1(:,[2 3 4]),'LineWidth',lwdenv,'LineStyle','--','Color',[0 0 1]);
                    line(gmin1(:,1),gmin1(:,5),'LineWidth',lwdenv,'LineStyle','--','Color',[0.3 0.3 0.2]);
                    
                else % superimposition in normal coordinates
                    
                    % Transform mmd in normal coordinates using
                    % tr observations
                    mmdinvr = FSMinvmmd(mmd,v,'n',tr);
                    
                    plot(mmdinvr(1:(tr-mmd(1,1)),1),mmdinvr(1:(tr-mmd(1,1)),3),'LineWidth',lwd);
                    
                    
                    for iq=[1 2 5 6]
                        if quant(iq)==0.5
                            col=[1 0.69 0.39];
                        elseif quant(iq)>=0.01 &&  quant(iq)<=0.999
                            col='b';
                        else
                            col='r';
                        end
                        line(vaxis(1:2)',[ninv(iq);ninv(iq)],'color',col,'LineWidth',lwdenv,'LineStyle','--','Tag','env');
                    end
                    
                end
                
                
                strtemp=['$d_{min}(m,' int2str(tr) ')$'];
                % get the position of the current pane
                gposcurax=get(gca,'position');
                % set the width and the height of the current pane but keep
                % unaltered the distance from bottom left corner
                set(gca,'position',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3)*1.3,gposcurax(4)*1.21]);
                % Subplots located on the left hand side
                plleft=1:nc:(nr*(nc-1)+1);
                % Subplots located on the right hand side
                plright=nc:nc:nr*nc;
                % For all plots not located on the left and the right hand side
                % delete numbers on the y axis (YTickLabels)
                getYTickLab=get(gca,'YTickLabel');
                if isempty(intersect(jwind,[plleft plright])) && nr*nc>1
                    set(gca,'YTickLabel',[]);
                end
                % For all plots not located on the right hand side put the
                % yaxis location on the right
                if ~isempty(intersect(jwind,plright)) && nr*nc>1  && nc>1
                    set(gca,'YAxisLocation','right');
                end
                % For all plots which are not on the bottom hand side delete
                % numbers on the x axis (XTickLabels)
                if isempty(intersect(jwind,(nr-1)*nc+1:(nr*nc)))
                    set(gca,'XTickLabel',[]);
                else
                    % For the plots on the bottom side add the xlabel
                    xlabel('Subset size m');
                end
                annotation(figure2,...
                    'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)+0.05],...
                    'String',{strtemp},'Tag','mes1',...
                    PrVaCell{:},'FontSize',fsizeannot);
                
                hold('off');
                jwind=jwind+1;
            end
            
            if sto==1
                if ~isempty(resuper) && ~isempty(intersect(resuper,tr))
                    % Write on the plot the reason why the procedure stopped
                    annotation(figure2,...
                        'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)],...
                        'String',{mes},'Tag','mes2',...
                        PrVaCell{:},'FontSize',fsizeannot);
                    
                    % Unless for all plots not located on the right hand side
                    % For the final plot put the yaxis location on the right
                    % Unless it is the first on the left hand side
                    if isempty(intersect(jwind-1,plleft)) && nr*nc>1  && nc>1
                        set(gca,'YTickLabel',getYTickLab);
                        set(gca,'YAxisLocation','right');
                    end
                    
                    % add X label again to the last plot
                    % set(gca,'XTickLabel',get(get(figure1,'Children'),'XTick'))
                    xlabel('Subset size m');
                    
                    % if jwind=2 than the width of the last plot is enlarged to
                    % full screen
                    if jwind==2 && nr*nc>1
                        set(gca,'Position',[0.1 0.1 0.85 0.85])
                        % Relocate the two messages
                        set(findall(gcf,'Tag','mes1'),'Units','normalized','Position',[0.3 0.7 0.5 0.1])
                        set(findall(gcf,'Tag','mes2'),'Units','normalized','Position',[0.3 0.6 0.5 0.1])
                    end
                end
                
                break;
            end
            if ~isempty(resuper) && ~isempty(intersect(resuper,tr))
                if jwind==nr*nc+1
                    jwind=1;
                    figure2 = figure('PaperSize',[20.98 29.68],'Name',['Resuperimposed envelopes #' int2str(resup)]);
                    resup=resup+1;
                    % Create axes (Following line should not be necessary
                    % axes('Parent',figure2);
                end
            end
            
        end
        
        %% Stage 2a: subset validation
        % In this part we check whether the subset is homogeneous. In other
        % words we verify conditions H1 and H2
        % tr= m^\dagger+k+1
        % m^\dagger+k=tr-1
        % m*=mmd(ii,1)
        % Condition H2
        % Check if stopping rule takes place at m* <m^\dagger+k
        if (mmd(ii,1)<tr-1)
            % Condition H2b and H2a
            if sum(gmin1(ii+1:end,4)>mmd(ii+1:size(gmin1,1),2))>0
                if msg
                    disp(['Subsample of ' int2str(tr-1) ' units is not homogeneous because the curve was above 99.99% and later it was below 1%']);
                    disp('----------------------------------------');
                end
                % Find m^{1%} that is the step where mmd goes below the 1%
                % threshold for the first time
                % ginfd = concatenate all the steps from m^*+1 to m^\dagger+k-1
                gfind=[gmin1(i+1:end,1) gmin1(i+1:end,4)>mmd(i+1:size(gmin1,1),2)];
                % select from gfind the steps where mmd was below 1% threshold
                % gfind(1,1) contains the first step where mmd was below 1%
                gfind=gfind(gfind(:,2)>0,1);
                % find maximum in the interval m^\dagger=mmd(i,1) to the step
                % prior to the one where mmd goes below 1% envelope
                tr=sortrows(mmd(i:gfind(1,1)-mmd(1,1),1:2),2);
                tr=tr(end,1);
                if msg
                    disp('Probably there are two overlapping groups');
                    disp(['Using the criterion of the maximum, the group of homogenous obs. is=' int2str(tr)]);
                end
                % tr is redefined and is associated with the step associated to
                % the maximum value of d_min
                % try=sormcl[rows(sormcl),1]+1;
                tr=tr+1;
            else
                if msg
                    disp(['Subsample of ' int2str(tr-1) ' units is homogeneous']);
                end
            end
        else
        end
        ndecl=n-tr+1;
    else
        ndecl=n-mmd(i,1);
    end
    
    if msg
        disp('----------------------------');
        disp('Final output');
        disp(['Number of units declared as outliers=' int2str(ndecl)]);
    end
else
    if msg
        disp('Sample seems homgeneous, no outlier has been found');
    end
    ndecl=0;
end

%% End of the forward search

if msg
    if isempty(bonflev)
        disp('Summary of the exceedances');
        disp(nout);
    end
    
    if ~isempty(extram3)
        disp(extram3);
    else
    end
    
    if ~isempty(extram2)
        disp(extram2);
    else
    end
end

group=ones(n,1);

if ndecl>0
    % Now find the list of the units declared as outliers
    
    % bsel=~isnan(bb(:,tr-init+1));
    % ListOut=setdiff(1:n,bsel,1);
    if n<5000
        outliers=seq(isnan(bb(:,end-ndecl)));
    else
        [~,bb] = FSMbsb(Y,bs,'bsbsteps',n-ndecl,'init',n-ndecl,'nocheck',1);
        outliers=seq(isnan(bb));
    end
    group(outliers)=2;
    
else
    outliers=NaN;
end

%compute locatione and covariance matrix
goodobs=setdiff(1:n,outliers);
loc=mean(Y(goodobs,:));
cova=cov(Y(goodobs,:));
md=mahalFS(Y,loc,cova);
%% Scatter plot matrix with the outliers shown with a different symbol
if v<=15
    if isstruct(plo) || (~isstruct(plo) && plo~=0)
        figure('Tag','pl_spm_outliers');
        if ~isstruct(plo)
            spmplot(Y,group,1);
        else
            spmplot(Y,group,plo);
        end
        set(gcf,'Name','FSM: scatter plot matrix with outliers highlighted');
    end
else
    if msg
        disp('There are more than 15 variables, spmplot is not shown')
    end
end
%% Structure returned by function FSM
% If you wish that the output also contains the list of units not declared
% as outliers, please uncomment the two following lines.
%ListIn=seq(~isnan(bb(:,end-ndecl)));
%out.ListIn=ListIn;

out=struct;
out.outliers=outliers;
out.loc=loc;
out.cov=cova;
out.md=md;
out.mmd=mmd;
out.Un=Un;
if isempty(bonflev)
    out.nout=nout;
else
    out.nout =NaN;
end
out.class  =  'FSM';

%% Callback functions used to "pin" quantile labels and vertical line to axes.

    function figure1_precallback(obj,~)
        % When the panning action starts, the quantile labels and the
        % vertical line which identifies the final part of the search are
        % made invisible.
        
        hanno = findall(obj,'Tag', 'quantile_label');
        set(hanno,'Visible','off');
        hanno2 = findall(obj, 'Tag' , 'FinalPartLine');
        set(hanno2,'Visible','off');
        
    end

    function figure1_postcallback(obj,evd)
        %When the panning action terminates, the new positions of the
        %labels and of the vertical line are set.
        
        % axis limits
        xxlim = get(evd.Axes,'XLim');
        xmin = xxlim(1); xmax=xxlim(2);
        
        % QUANTILES ANNOTATION: the handles
        hanno1 = findall(obj, 'Tag' , 'quantile_label');
        % the quantiles values at each step of the search
        xy = get(hanno1,'UserData');
        % the step 'init'
        init = xy{1,1}(1);
        % the max between init and the lower axis limit
        x = floor(max(xmin , init));
        % the position in the quantiles matrix of x
        xi = x - init + 1;
        % xp is the actual position of the annotation in the plot: normally
        % xp = x, but we move it towards right a bit (by 3 steps) when x is
        % so that the annotation would overlap with the quantiles values
        % labels.
        if x-xmin>0
            xp = x;
        else
            xp = x-3;
        end
        % the figure coordinates which correspond to the annotations positions.
        [y99999, y9999 , y999 , y50 , y99 , y1] = ...
            deal(xy{1,1}(xi,2) , xy{2,1}(xi,2) , xy{3,1}(xi,2) , xy{4,1}(xi,2) , xy{5,1}(xi,2) , xy{6,1}(xi,2));
        [figx, figy1] = dsxy2figxy(evd.Axes, xp, y1);
        [figx, figy99] = dsxy2figxy(evd.Axes, xp, y99);
        [figx, figy50] = dsxy2figxy(evd.Axes, xp, y50);
        [figx, figy999] = dsxy2figxy(evd.Axes, xp, y999);
        [figx, figy9999] = dsxy2figxy(evd.Axes, xp, y9999);
        [figx, figy99999] = dsxy2figxy(evd.Axes, xp, y99999);
        positionValCell(1,1) = {[figx figy99999 0 0]};
        positionValCell(2,1) = {[figx figy9999 0 0]};
        positionValCell(3,1) = {[figx figy999 0 0]};
        positionValCell(4,1) = {[figx figy50 0 0]};
        positionValCell(5,1) = {[figx figy99 0 0]};
        positionValCell(6,1) = {[figx figy1 0 0]};
        % set the new figure coordinates and make the annotations visible again
        set(hanno1 , {'Position'} , positionValCell, 'Visible','on');
        
        % Uncomment the next lines to display in the command window the new
        % positions
        % disp(['Min axis value:' int2str(xmin)]);
        % disp(['x value:' int2str(x)]);
        % disp(['xi value:' int2str(xi)]);
        % disp(['xp value:' int2str(xp)]);
        
        % VERTICAL LINE which identifies the final part of the search.
        % For that line we only need to recompute the X position.
        hanno2 = findall(obj, 'Tag' , 'FinalPartLine');
        LineData = get(hanno2,'UserData');
        istep = LineData(1); yl1 = LineData(2); yl2 = LineData(3);
        if istep > xmax, istep = xmax; end
        if istep < xmin, istep = xmin; end
        [figx, figy]  = dsxy2figxy(evd.Axes, istep, yl1);
        
        positionLineOri=get(hanno2,'Position');
        positionLineOri(1)=figx;
        
        positionLineCell = {positionLineOri};
        xCell = {[figx figx]};
        
        set(hanno2 , {'Position'} , positionLineCell, {'X'} , xCell, 'Visible','on');
        
    end

end
%FScategory:MULT-Multivariate
