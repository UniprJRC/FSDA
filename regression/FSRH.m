function [out]=FSRH(y,X,varargin)
%FSRH gives an automatic outlier detection procedure in heteroskedastic linear regression
%
%
%<a href="matlab: docsearch('fsrh')">Link to the help function</a>
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
%   intercept   : If 1, a model with constant term will be fitted (default),
%                 if 0, no constant term will be included.
%           h   : The number of observations that have determined the least
%                 trimmed squares estimator. h is an integer greater or
%                 equal than p but smaller then n. Generally if the purpose
%                 is outlier detection h=[0.5*(n+p+1)] (default value). h
%                 can be smaller than this threshold if the purpose is to find
%                 subgroups of homogeneous observations.
%       nsamp   : Number of subsamples which will be extracted to find the
%                 robust estimator. If nsamp=0 all subsets will be extracted.
%                 They will be (n choose p).
%                 Remark: if the number of all possible subset is <1000 the
%                 default is to extract all subsets otherwise just 1000.
%       lms     : Scalar,  vector or structure.
%                 lms specifies the criterion to use to find the initlal
%                 subset to initialize the search (LMS, LTS with
%                 concentration steps, LTS without concentration steps
%                 or subset supplied directly by the user).
%                 The default value is 1 (Least Median of Squares
%                 is computed to initialize the search). On the other hand,
%                 if the user wants to initialze the search with LTS with
%                 all the default options for concentration steps then
%                 lms=2. If the user wants to use LTS without
%                 concentration steps, lms can be a scalar different from 1
%                 or 2. If lms is a struct it is possible to control a
%                 series of options for concentration steps (for more
%                 details see option lms inside LXS.m)
%                 LXS.m. 
%                 If, on the other hand, the user wants to initialize the
%                 search with a prespecified set of units there are two
%                 possibilities
%                 1) lms can be a vector with length greater than 1 which 
%                 contains the list of units forming the initial subset.
%                 For example, if the user wants to initialize the search
%                 with units 4, 6 and 10 then lms=[4 6 10];
%                 2) lms is a struct which contains a field named bsb which
%                 contains the list of units to initialize the search. For
%                 example, in the case of simple regression through the
%                 origin with just one explanatory variable, if the user
%                 wants to initialize the search with unit 3 then
%                 lms=struct; lms.bsb=3;
%       plots   : Scalar.
%                 If plots=1 (default) the plot of minimum deletion
%                 residual with envelopes based on n observations and the
%                 scatterplot matrix with the outliers highlighted is
%                 produced.
%                 If plots=2 the user can also monitor the intermediate
%                 plots based on envelope superimposition.
%                 else no plot is produced.
%       init    : scalar which specifies the initial subset size to start
%                 monitoring exceedances of minimum deletion residual, if
%                 init is not specified it set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%        scoring: if scoring=1 a scoring algorithm is executed, if scoring=0,
%                 standard algorithm is executed.
%       exact   : scalar, if it is equal to 1 the calculation of the quantiles
%                 of the T and F distribution is based on functions finv
%                 and tinv from the Matlab statistics toolbox, else the
%                 calculations of the former quantiles is based on
%                 functions invcdff and invcdft.
%                 The solution has a tolerance of 1e-8 (change variable tol
%                 in files invcdff.m and invcdft.m if required
%                 Remark: the use of functions tinv and finv is more precise
%                 but requires more time.
%       nocheck : Scalar. If nocheck is equal to 1 no check is performed on
%                 matrix y and matrix X. Notice that y and X are left
%                 unchanged. In other words the additional column of ones
%                 for the intercept is not added. As default nocheck=0.
%    bivarfit : This option adds one or more least square lines, based on
%                 SIMPLE REGRESSION of y on Xi, to the plots of y|Xi.
%                 bivarfit = ''
%                   is the default: no line is fitted.
%                 bivarfit = '1'
%                   fits a single ols line to all points of each bivariate
%                   plot in the scatter matrix y|X.
%                 bivarfit = '2'
%                   fits two ols lines: one to all points and another to
%                   the group of the genuine observations. The group of the
%                   potential outliers is not fitted.
%                 bivarfit = '0'
%                   fits one ols line to each group. This is useful for the
%                   purpose of fitting mixtures of regression lines.
%                 bivarfit = 'i1' or 'i2' or 'i3' etc.
%                   fits an ols line to a specific group, the one with
%                   index 'i' equal to 1, 2, 3 etc. Again, useful in case
%                   of mixtures.
%       multivarfit : This option adds one or more least square lines, based on
%                 MULTIVARIATE REGRESSION of y on X, to the plots of y|Xi.
%                 multivarfit = ''
%                   is the default: no line is fitted.
%                 multivarfit = '1'
%                   fits a single ols line to all points of each bivariate
%                   plot in the scatter matrix y|X. The line added to the
%                   scatter plot y|Xi is avconst + Ci*Xi, where Ci is the
%                   coefficient of Xi in the multivariate regression and
%                   avconst is the effect of all the other explanatory
%                   variables different from Xi evaluated at their centroid
%                   (that is overline{y}'C))
%                 multivarfit = '2'
%                   equal to multivarfit ='1' but this time we also add the
%                   line based on the group of unselected observations
%                   (i.e. the normal units).
%      labeladd : If this option is '1',  we label the outliers with the
%                 unit row index in matrices X and y. The default value is
%                 labeladd='', i.e. no label is added.
%       nameX  :  cell array of strings of length p containing the labels of
%                 the variables of the regression dataset. If it is empty
%                 (default) the sequence X1, ..., Xp will be created
%                 automatically
%       namey  :  character containing the label of the response
%       ylim   :  vector with two elements controlling minimum and maximum
%                 on the y axis. Default value is '' (automatic scale)
%       xlim   :  vector with two elements controlling minimum and maximum
%                 on the x axis. Default value is '' (automatic scale)
%      bonflev  : option to be used if the distribution of the data is
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
%       msg    :  scalar which controls whether to display or not messages
%                 on the screen
%                 If msg==1 (default) messages are displayed on the screen about
%                   step in which signal took place
%                 else no message is displayed on the screen
% bsbmfullrank :  scalar which tells how to behave in case subset at step m
%                 (say bsbm) produces a non singular X. In other words,
%                 this options controls what to do when rank(X(bsbm,:)) is
%                 smaller then number of explanatory variables. If
%                 bsbmfullrank =1 (default) these units (whose number is
%                 say mnofullrank) are constrained to enter the search in
%                 the final n-mnofullrank steps else the search continues
%                 using as estimate of beta at step m the estimate of beta
%                 found in the previous step.
%
%
% Output:
%
%  The output consists of a structure 'out' containing the following fields:
% out.ListOut=  k x 1 vector containing the list of the units declared as
%               outliers or NaN if the sample is homogeneous
% out.mdr    =  (n-init) x 2 matrix
%               1st col = fwd search index
%               2nd col = value of minimum deletion residual in each step
%               of the fwd search
% out.Un     =  (n-init) x 11 Matrix which contains the unit(s) included
%               in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
% out.nout    = 2 x 5 matrix containing the number of times mdr went out
%               of particular quantiles.
%               First row contains quantiles 1 99 99.9 99.99 99.999.
%               Second row contains the frequency distribution.
% out.constr  = This output is produced only if the search found at a
%               certain step a non singular matrix X. In this case the
%               search run in a constrained mode, that is including the
%               units which produced a singular matrix in the last n-constr
%               steps. out.constr is a vector which contains the list of
%               units which produced a singular X matrix
%
% See also: FSReda, LXS.m
%
% References:
%
%   Atkinson A.C., Riani M. and Torti F. (2015), Robust methods for
%   heteroskedastic regression, submitted (ART)
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearch('fsrh')">Link to the help page for this function</a>
%
% Last modified 06-Feb-2015

% Examples:

%{
    XX=load('tradeH.txt')
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    %% Plot the data
    plot(X,y,'o')
    fs=14;
    ylabel('Value','FontSize',fs)
    xlabel('Quantity','FontSize',fs)

    % Run the automatic outlier detection procedure 
    out=FSRH(y,X,'init',round(length(y)/2),'plots',1,'ylim',[1.6 3]);

    % Create figure 3 of ART
    figure
    subplot(2,2,1)
    n=length(y);
    seq=1:n;
    sel=setdiff(seq,out.ListOut);
    hold('on')
    plot(X(sel),y(sel),'o')
    plot(X(out.ListOut),y(out.ListOut),'rx','MarkerSize',12,'LineWidth',2)
    % fs=12;
    % ylabel('Value','FontSize',fs)
    % xlabel('Quantity','FontSize',fs)
    % set(gca,'FontSize',fs)
    %
    subplot(2,2,2)
    plot(out.Hetero(:,1),out.Hetero(:,2))
    xlabel('Subset size m')
    kk=20;
    xlim([out.Hetero(1,1) out.Hetero(end,1)+kk])
    ylim([1.7 2.7])
    title('\alpha')
    subplot(2,2,3)
    plot(out.Hetero(:,1),log(out.Hetero(:,3)))
    title('log(\theta)')
    xlim([out.Hetero(1,1) out.Hetero(end,1)+kk])
    %ylim([5 7.5])
    xlabel('Subset size m')
    subplot(2,2,4)
    plot(out.S2(:,1),out.S2(:,2))
    xlim([out.Hetero(1,1) out.Hetero(end,1)+kk])
    ylim([0 300000])
    title('\sigma^2')
    xlabel('Subset size m')
%}




%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options
lowexceed=0;

% intcolumn = the index of the first constant column found in X, or empty.
% Used here to check if X includes the constant term for the intercept.
% The variable 'intercept' will be used later for plotting.
intcolumn = find(max(X,[],1)-min(X,[],1) == 0,1);
if any(intcolumn) && p>1;
    intercept=1;
else
    intercept=0;
end

% If the number of all possible subsets is <1000 the default is to extract
% all subsets, otherwise just 1000.
ncomb=bc(n,p);
nsampdef=min(1000,ncomb);
% Notice that a fast approximation of the bc computed above is:
% ncomb=floor(exp( gammaln(n+1) - gammaln(n-p+1) - gammaln(p+1) ) + .5);

% The default value of h is floor(0.5*(n+p+1))
hdef=floor(0.5*(n+p+1));
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end
% ini0=init;

options=struct('h',hdef,...
    'nsamp',nsampdef,'lms',1,'plots',1,...
    'init',init,'exact',1,...
    'labeladd','','bivarfit','','multivarfit','',...
    'xlim','','ylim','','nameX','','namey','',...
    'msg',1,'nocheck',0,'intercept',1,'bonflev','',...
    'bsbmfullrank',1,'scoring',1);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('Error:: number of supplied options is invalid. Probably values for some parameters are missing.');
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

exact=options.exact;
init=options.init;
h=options.h;
lms=options.lms;
plo=options.plots;
nsamp=options.nsamp;
labeladd=options.labeladd;
bivarfit=options.bivarfit;
multivarfit=options.multivarfit;
xlimx=options.xlim;
ylimy=options.ylim;
msg=options.msg;
bsbmfullrank=options.bsbmfullrank;
scoring=options.scoring;

bonflev=options.bonflev;
if ~isempty(bonflev)
    if bonflev<1;
        [gbonf] = FSRbonfbound(n,p,'prob',bonflev,'init',init);
        bonfthresh=gbonf;
    else
        bonfthresh=bonflev*ones(n-init,1);
    end
end

%% Start of the forward search

seq=1:n;

% Compute theoretical envelopes based on all observations
quant=[0.99;0.999;0.9999;0.99999;0.01;0.5;0.00001];
% Compute theoretical envelops for minimum deletion residual based on all
% the observations for the above quantiles.
[gmin] = FSRenvmdr(n,p,'prob',quant,'init',init,'exact',exact);
% gmin = the matrix which contains envelopes based on all observations.
% 1st col of gmin = fwd search index
% 2nd col of gmin = 99% envelope
% 3rd col of gmin = 99.9% envelope
% 4th col of gmin = 99.99% envelope
% 5th col of gmin = 99.999% envelope
% 6th col of gmin = 1% envelope
% 7th col of gmin = 50% envelope
% Thus, Set the columns of gmin where the theoretical quantiles are located.
[c99 , c999 , c9999 , c99999 , c001 , c50] = deal(2,3,4,5,6,7);

iter=0;

% Use as initial subset the one supplied by the user or the best according
% to LMS or LTS

if length(lms)>1 || (isstruct(lms) && isfield(lms,'bsb'));
    if length(lms)>1
        bs=lms;
    else
        bs=lms.bsb;
    end
    if init<length(bs)
        init=length(bs);
    end
    
    % Compute Minimum Deletion Residual for each step of the search
    [mdr,Un,bb,~,~,Hetero,WEI] = FSRHmdr(y,X,bs,'init',init,'plots',0,'nocheck',1,'msg',msg,'scoring',scoring);
    
    if size(mdr,2)<2
        if length(mdr)>=n/2;
            disp('More than half of the observations produce a singular X matrix')
            disp('X is badly defined')
            disp('If you wish to run the procedure using for updating the values of beta of the last step in which there was fll rank use option bsbmfullrank=0')
            out.ListOut=setdiff(seq,mdr);
            
        else
            disp('Bad starting point which produced a singular matrix, please restart the search from a different starting point or use option bsbmfullrank=0 ')
            
        end
        
        out.mdr = NaN;
        out.Un  = NaN;
        out.nout= NaN;
        return
    end
else % initial subset is not supplied by the user
    % Find initial subset to initialize the search
    [out]=LXS(y,X,'lms',lms,'h',h,'nsamp',nsamp,'nocheck',1,'msg',msg);
    
    if out.s0==0
        disp('More than half of the observations produce a linear model with a perfect fit')
        % Just return the outliers found by LXS
        %out.ListOut=out.outliers;
        %return
    end
    
    bs=out.bs;
    mdr=0;
    constr='';
    
    
    while size(mdr,2)<2 && iter <6
        % Compute Minimum Deletion Residual for each step of the search
        % The instruction below is surely executed once.
        [mdr,Un,bb,~,S2,Hetero,WEI] = FSRHmdr(y,X,bs,'init',init,'plots',0,'nocheck',1,'msg',msg,'constr',constr,'bsbmfullrank',bsbmfullrank,'intercept',intercept,'scoring',scoring);
        
        % If FSRmdr run without problems mdr has two columns. In the second
        % column it contains the value of the minimum deletion residual
        % monitored in each step of the search
        
        % If mdr has just one columns then one of the following two cases took place:
        % isnan(mdr)=1 ==> in this case initial subset was not full rank
        % mdr has just one column ==> in this case, even if the initial
        %    subset was full rank, the search has found at a certain step
        %    m<n/2 a list of units which produce a singular matrix. In this
        %    case LXS is rerun excluding these units which gave rise to a
        %    non singular matrix
        
        if size(mdr,2)<2
            if length(mdr)>=n/2;
                disp('More than half of the observations produce a singular X matrix')
                disp('If you wish to run the procedure using for updating the values of beta of the last step in which there was fll rank use option bsbmfullrank=0')
                
                out.ListOut=setdiff(seq,mdr);
                
                return
            elseif isnan(mdr(1,1))
                % INITIAL SUBSET WAS NOT FULL RANK
                % restart LXS without the units forming
                % initial subset
                bsb=setdiff(seq,out.bs);
                [out]=LXS(y(bsb),X(bsb,:),'lms',lms,'nsamp',nsamp,'nocheck',1,'msg',msg);
                bs=bsb(out.bs);
                
                
            else
                % INITIAL SUBSET WAS FULL RANK BUT THE SEARCH HAS FOUND A
                % SET OF OBSERVATIONS CONSTR <n/2  WHICH PRODUCED A SINGULAR
                % MATRIX. IN THIS CASE NEW LXS IS BASED ON  n-constr OBSERVATIONS
                iter=iter+1;
                bsb=setdiff(seq,mdr);
                constr=mdr;
                [out]=LXS(y(bsb),X(bsb,:),'lms',lms,'nsamp',nsamp,'nocheck',1,'msg',msg);
                bs=bsb(out.bs);
            end
        end
    end
    
    
end


if iter >=5
    disp('No convergence')
    disp('No convergence')
    disp('No convergence')
    %     out.mdr = NaN;
    %     out.Un  = NaN;
    %     out.nout= NaN;
    %     error('no convergence')
end
bool=mdr(:,1)>=init;
mdr=mdr(bool,:);
gmin=gmin(gmin(:,1)>=init,:);
Hetero=Hetero(bool,:);

% Store in nout the number of times the observed mdr (d_min) lies above:
[out99 , out999 , out9999 , out99999 , out001] = deal( ...
    mdr(mdr(:,2)>gmin(:,c99),:) , ...       % the 99% envelope
    mdr(mdr(:,2)>gmin(:,c999),:) , ...      % the 99.9% envelope
    mdr(mdr(:,2)>gmin(:,c9999),:) , ...     % the 99.99% envelope
    mdr(mdr(:,2)>gmin(:,c99999),:) , ...    % the 99.999% envelope
    mdr(mdr(:,2)<gmin(:,c001),:) );         % the 1% envelope

nout = [[1 99 999 9999 99999]; ...
    [size(out001,1) size(out99,1) size(out999,1) size(out9999,1) size(out99999,1)]];

% NoFalseSig = boolean linked to the fact that the signal is good or not
NoFalseSig=0;

% NoFalseSig is set to 1 if the condition for an INCONTROVERTIBLE SIGNAL is
% fulfilled.
n9999 = nout(2,nout(1,:)==9999);
if (n9999>=10);
    NoFalseSig=1;
    if msg;
        disp('Observed curve of r_min is at least 10 times greater than 99.99% envelope'); % exact number is int2str(n9999)
        disp('--------------------------------------------------');
    end
end

% Divide central part from final part of the search
istep = n-floor(13*sqrt(n/200));

%% Part 1. Signal detection and validation
nmdr=size(mdr,1);
if nmdr<4
    error('ratio n/p too small; modify init (i.e. decrease initial subset size)')
end
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

%% Stage 1a: signal detection
% Signal dection is based on monitoring consecutive triplets or single
% extreme values

% Signal detection loop
for i=3:nmdr;
    % Check if signal must be based on consecutive exceedances of envelope
    % of mdr or on exceedance of global Bonferroni level
    if isempty(bonflev)
        
        if i<istep-init+1; % CENTRAL PART OF THE SEARCH
            % Extreme triplet or an extreme single value
            % Three consecutive values of d_min above the 99.99% threshold or 1
            % above 99.999% envelope
            if ((mdr(i,2)>gmin(i,c9999) && mdr(i+1,2)>gmin(i+1,c9999) && mdr(i-1,2)>gmin(i-1,c9999)) || mdr(i,2)>gmin(end,c99) || mdr(i,2)>gmin(i,c99999));
                if msg
                    disp(['Tentative signal in central part of the search: step m=' int2str(mdr(i,1)) ' because']);
                end
                if (mdr(i,2)>gmin(i,c9999) && mdr(i+1,2)>gmin(i+1,c9999) && mdr(i-1,2)>gmin(i-1,c9999));
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.99% and rmin(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99.99% and rmin(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99.99%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.99\%$ and $r_{min}(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99.99\%$ and $r_{min}(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99.99\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end
                
                if (mdr(i,2)>gmin(i,c99999));
                    if msg
                        disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')>99.999%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.999\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end;
                
                if (mdr(i,2)>gmin(end,c99));
                    if msg
                        disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the central part of the search.']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
                    mdrsel=mdr(i,1:2);
                    NoFalseSig=1; % i.e., no need of further validation
                end;
                
                '------------------------------------------------';
                
                signal=1;
                
            elseif (mdr(i,2)<gmin(i,end)) && lowexceed==1  % && mdr(i,1)>round(n/2); % exceedance of the lower band
                if msg
                    disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')<0.00001% Bonferroni signal in the central part of the search.']);
                end
                strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')<0.00001\%$  (Bonferroni signal)'];
                mdrsel=mdr(i,1:2);
                NoFalseSig=1; % i.e., no need of further validation
                signal=2;
                mdag=mdr(i,1);
            else
            end
        elseif i<size(mdr,1)-1; % FINAL PART OF THE SEARCH
            % Extreme couple adjacent to an exceedance
            % Two consecutive values of mdr above the 99.99% envelope and 1 above 99%
            if ((mdr(i,2)>gmin(i,c999) && mdr(i+1,2)>gmin(i+1,c999) && mdr(i-1,2)>gmin(i-1,c99)) || (mdr(i-1,2)>gmin(i-1,c999) && mdr(i,2)>gmin(i,c999) && mdr(i+1,2)>gmin(i+1,c99)) || mdr(i,2)>gmin(end,c99) || mdr(i,2)>gmin(i,c99999));
                if msg
                    disp(['Signal in final part of the search: step ' num2str(mdr(i,1)) ' because']);
                end
                if (mdr(i,2)>gmin(i,c999) && mdr(i+1,2)>gmin(i+1,c999) && mdr(i-1,2)>gmin(i-1,c99));
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i+1,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i-1,1)) ',' int2str(n) ')>99%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9\%$ and $r_{min}(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99\%$ and $r_{min}(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99.9\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end
                
                if (mdr(i-1,2)>gmin(i-1,c999) && mdr(i,2)>gmin(i,c999) && mdr(i+1,2)>gmin(i+1,c99));
                    if msg
                        disp(['rmin('  int2str(mdr(i-1,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.9% and rmin('  int2str(mdr(i+1,1)) ',' int2str(n) ')>99%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9\%$ and $r_{min}(' int2str(mdr(i-1,1)) ',' int2str(n) ')>99.9\%$ and $r_{min}(' int2str(mdr(i+1,1)) ',' int2str(n) ')>99\%$'];
                    mdrsel=mdr(i-1:i+1,1:2);
                end
                
                if (mdr(i,2)>gmin(end,c99));
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the final part of the search.']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
                    mdrsel=mdr(i:i,1:2);
                end
                
                % Extreme single value above the upper threshold
                if mdr(i,2)>gmin(i,c99999);
                    if msg
                        disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99.999%']);
                    end
                    strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.999\%$'];
                    mdrsel=mdr(i:i,1:2);
                end
                
                '------------------------------------------------';
                % Signal is always considered true if it takes place in the
                % final part of the search
                NoFalseSig=1;
                signal=1;
            end
        elseif (mdr(i,2)>gmin(i,c999) || mdr(i,2)>gmin(end,c99)) && i==size(mdr,1)-1;
            % potential couple of outliers
            signal=1;
            if msg
                disp('Signal is in penultimate step of the search');
            end
            
            if (mdr(i,2)>gmin(i,c999));
                if msg
                    disp(['rmin(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9%']);
                end
                strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99.9\%$'];
            end;
            
            if (mdr(i,2)>gmin(end,c99));
                if msg
                    disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')>99% at final step: Bonferroni signal in the final part of the search.']);
                end
                strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step (Bonferroni signal)'];
            end
            mdrsel=mdr(i:i,1:2);
        elseif  mdr(i,2)>gmin(i,c99) && i==size(mdr,1);
            % a single outlier
            signal=1;
            if msg
                disp('Signal is in final step of the search');
            end
            strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ at final step'];
            mdrsel=mdr(i:i,1:2);
            
        elseif mdr(i,2)<gmin(i,end) && lowexceed==1; % exceedance of the lower extreme envelope
            signal=2;
            if msg
                disp(['rmin('  int2str(mdr(i,1)) ',' int2str(n) ')<0.0001%']);
            end
            strplot=['$r_{min}(' int2str(mdr(i,1)) ',' int2str(n) ')<0.0001\%$'];
            mdrsel=mdr(i:i,1:2);
            mdag=mdr(i,1);
        end
        
        %% Stage 1b: signal validation
        if signal==1;
            if msg
                disp('-------------------')
                disp('Signal validation exceedance of upper envelopes');
            end
            % mdag is $m^\dagger$
            mdag=mdr(i,1);
            
            if mdr(i,1)<n-2;
                % Check if the signal is incontrovertible
                % Incontrovertible signal = 3 consecutive values of d_min >
                % 99.999% threshold
                if mdr(i,2)>gmin(i,c99999) && mdr(i-1,2)>gmin(i-1,c99999) &&  mdr(i+1,2)>gmin(i+1,c99999);
                    if msg
                        disp(['3 consecutive values of r_min greater than 99.999% envelope in step mdag= ' int2str(mdr(i,1))]);
                    end
                    NoFalseSig=1;
                    extram3='Extreme signal';
                end;
            else
                NoFalseSig=1;
            end;
            
            % if the following statement is true, observed curve of r_min is
            % above 99.99% and later is below 1%: peak followed by dip
            if size(mdr,1)>mdag-mdr(1,1)+31;
                if sum(mdr(i+1:i+31,2)<gmin(i+1:i+31,c001))>=2;
                    NoFalseSig=1;  % Peak followed by dip
                    extram2='Peak followed by dip (d_min is above 99.99% threshold and in the sucessive steps goes below 1% envelope';
                end;
            else
                if sum(mdr(i+1:end,2) < gmin(i+1:end,c001))>=2;
                    NoFalseSig=1;  %Peak followed by dip in the final part of the search';
                    extram2='Peak followed by dip (r_min is above 99.99% threshold and in the sucessive steps goes below 1% envelope)';
                end;
            end;
            
            % if at this point NoFalseSig==0 it means that:
            % 1) n9999<10
            % 2) signal tool place in the central part of the search
            % 3) signal was not incontrovertible
            % 4) there was not a peak followed by dip
            if NoFalseSig==0;
                % Compute the final value of the envelope based on
                % mdr(i+1,1)=mdagger+1 observations
                %[gval]=FSRenvmdr(mdag+1,p,'prob',0.01,'m0',mdag);
                [gval]=FSRenvmdr(mdag+1,p,'prob',0.01,'init',mdag);
                if mdr(i,2)<gval(1,2);
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
            if (NoFalseSig==1);
                if msg
                    disp('Validated signal');
                end
                break
            end
        elseif signal==2
            if msg
                disp('-------------------')
                disp('Signal validation exceedance of extreme lower envelope');
            end
            break
        end
    else
        % Outlier detection based on Bonferroni threshold
        if (mdr(i,2)>bonfthresh(i,end)) && mdr(i,1)>floor(0.5*n);
            if msg
                disp(['mdr(' int2str(mdr(i,1)) ',' int2str(n) ')>99% Bonferroni level']);
            end
            strplot=['$mdr(' int2str(mdr(i,1)) ',' int2str(n) ')>99\%$ (Bonferroni level)'];
            mdrsel=mdr(i:i,1:2);
            
            signal=1;
            break
        end
    end
    
end

%% Create figure containing mdr + envelopes based on all the observations.
if plo==1 || plo ==2
    % get screen size [left, bottom, width, height]
    scrsz = get(0,'ScreenSize');
    
    figure1 = figure('Position',[1 scrsz(4)/2.5 scrsz(3)/3 scrsz(4)/2],'PaperSize',[20.98 29.68],'Name','Envelopes based on all the observations');
    axes1 = axes('Parent',figure1);
    
    % Create a pan-handle for the figure ...
    hpan_figure1 = pan(figure1);
    % ... and listen to pan events using callback functions
    set(hpan_figure1,'ActionPreCallback',@fprecallback);
    set(hpan_figure1,'ActionPostCallback',@fpostcallback);
    
    %     % Create a zoom-handle for the figure ...
    %     hzoom_figure1 = zoom(figure1);
    %     % ... and listen to pan events using callback functions
    %     set(hzoom_figure1,'ActionPreCallback',@figure1_precallback);
    %     set(hzoom_figure1,'ActionPostCallback',@figure1_postcallback_zoom);
    
    if isempty(xlimx)
        xl1=init-3; xl2=mdr(end,1);
    else
        xl1=xlimx(1);
        xl2=xlimx(2);
    end
    
    if isempty(ylimy)
        
        yl1=min([gmin(:,c001);mdr(:,2)]);
        yl2=max([gmin(:,c999);mdr(:,2)]);
    else
        yl1=ylimy(1);
        yl2=ylimy(2);
    end
    
    xlim([xl1 xl2]);
    ylim([yl1 yl2]);
    
    
    kx=0; ky=0;
    
    box('on'); hold('all');
    
    plot(mdr(:,1),mdr(:,2));
    
    
    if isempty(bonflev)
        
        % 1%, 99%, 99.9% envelopes based on all the observations
        line(gmin(:,1),gmin(:,[c001 c99 c999]),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
        % 99.99% and 99.999% envelopes based on all the observations
        line(gmin(:,1),gmin(:,[c9999 c99999]),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
        if signal==2
            line(gmin(:,1),gmin(:,end),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
        end
        
        % 50% envelope based on all the observations
        line(gmin(:,1),gmin(:,c50),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[1 0.69 0.39]);
        
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
                PrVaCell{:});
        end
        
        % Create textbox with 99% label
        [figx, figy] = dsxy2figxy(gca, init, gmin(1,c99));
        
        if figy>=0 && figy<=1 && figx>=0 && figx<=1
            annotation(figure1,'textbox',[figx figy kx ky],...
                'String',{'99%'},...
                'UserData',[gmin(:,1) gmin(:,c99)],...
                PrVaCell{:});
        end
        
        % Create textbox with 50% label
        [figx, figy] = dsxy2figxy(gca, init, gmin(1,c50));
        
        if figy>=0 && figy<=1 && figx>=0 && figx<=1
            annotation(figure1,'textbox',[figx figy kx ky],...
                'String',{'50%'},...
                'UserData',[gmin(:,1) gmin(:,c50)],...
                PrVaCell{:});
        end
        
        if signal==1
            % Create textbox with 99.9% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c999));
            if figy>=0 && figy<=1 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'99.9%'},...
                    'UserData',[gmin(:,1) gmin(:,c999)],...
                    PrVaCell{:});
            end
            
            % Create textbox with 99.99% label
            [figx, figy] = dsxy2figxy(gca, init, gmin(1,c9999));
            if figy<=1 && figy>=0 && figx>=0 && figx<=1
                annotation(figure1,'textbox',[figx figy kx ky],...
                    'String',{'99.99%'},...
                    'UserData',[gmin(:,1) gmin(:,c9999)],...
                    PrVaCell{:});
            end
            
            if gmin(1,c99999)<=yl2
                % Create textbox with 99.999% label
                [figx, figy] = dsxy2figxy(gca, init, gmin(1,c99999));
                if figy<=1 && figy>=0 && figx>=0 && figx<=1
                    annotation(figure1,'textbox',[figx figy kx ky],...
                        'String',{'99.999%'},...
                        'UserData',[gmin(:,1) gmin(:,c99999)],...
                        PrVaCell{:});
                end
            end
        else
        end
        
        % Add string which informs about the step where signal took place
        if signal==1 || signal==2
            strsig=['Signal is in step $m=' int2str(mdag) '$ because'];
            stem(mdrsel(:,1),mdrsel(:,2),'LineWidth',1,...
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
        strmin=['$r_{min}(m,' int2str(n) ')$. '];
        annotation(figure1,'textbox',[0.5 0.9 kx ky],'String',{[strmin strsig]},...
            PrVaCell{:});
        
        annotation(figure1,'textbox',[0.5 0.85 kx ky],'String',{strplot},...
            PrVaCell{:});
        
        if istep<=xl2
            % Add vertical line which divides central part from final part of the
            % search
            [figx, figy] = dsxy2figxy(gca, istep, yl1);
            [figx, figy2] = dsxy2figxy(gca, istep, yl2);
            if figy2>=1
                figy2=1;
            else
            end
            if figx>=0;
                annotation(figure1,'line',[figx figx],[figy figy2],...
                    'UserData',[istep yl1 yl2],...
                    'Tag','FinalPartLine');
            end
        end
    else
        % Superimpose Bonferroni line to the plot
        line(gmin(:,1),bonfthresh(:,end),'Parent',axes1,'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
        
        % Property-value pairs which are common to the next latex annotations
        PrVaCell=cell(2,5);
        PrVaCell{1,1} = 'Interpreter'; PrVaCell{2,1} = 'latex';
        PrVaCell{1,2} = 'HorizontalAlignment'; PrVaCell{2,2} = 'center';
        PrVaCell{1,3} = 'FitBoxToText'; PrVaCell{2,3} = 'on';
        PrVaCell{1,4} = 'EdgeColor'; PrVaCell{2,4} = 'none';
        PrVaCell{1,5} = 'BackgroundColor'; PrVaCell{2,5} = 'none';
        
        
        
        if size(bonfthresh,2)>1
            % latex annotations informing that the envelopes are based on
            % all the observations
            strmin='Exceedance based on Bonferroni threshold';
            annotation(figure1,'textbox',[0.5 0.9 kx ky],'String',strmin,...
                PrVaCell{:});
            msg=['$r_{min}(' num2str(mdr(i,1)) ',' int2str(n) ')>' num2str(100*bonflev) '$\% envelope'];
            annotation(figure1,'textbox',[0.5 0.8 kx ky],'String',msg,PrVaCell{:});
        else
            strmin='Exceedance based on user supplied threshold';
            annotation(figure1,'textbox',[0.5 0.9 kx ky],'String',strmin,...
                PrVaCell{:});
            msg=['$r_{min}(' num2str(mdr(i,1)) ',' int2str(n) ')>$' num2str(bonflev)];
            annotation(figure1,'textbox',[0.5 0.8 kx ky],'String',msg,PrVaCell{:});
        end
    end
    if signal==1
        stem(mdr(i,1),mdr(i,2),'LineWidth',1,...
            'Color',[0.4784 0.06275 0.8941], 'DisplayName','Signal');
    end
end

%% Part 2: envelope resuperimposition
% if a validated signal tool place, superimposition of the envelopes starts
% from m^\dagger-1

if signal==1 || signal==2;
    if isempty(bonflev) && signal==1
        
        if msg
            disp('-------------------------------');
            disp(['Start resuperimposing envelopes from step m=' int2str(mdag-1)]);
        end
        
        if plo==2;
            % jwind is associated with subplot window number
            % nr is the number of row panes in the plot
            % nc is the number of columns panes in the plot
            jwind=1;
            nc=2;
            if mdr(i,1)>=n-2;
                nr=1;
            else
                nr=2;
            end
            
            figure2 = figure('PaperSize',[20.98 29.68],'Name','Resuperimposed envelopes #1');
            % Create axes
            axes('Parent',figure2);
        end
        
        % First resuperimposed envelope is based on mdag-1 observations
        % Notice that mdr(i,1) = m dagger
        for tr=(mdag-1):(n);
            % Compute theoretical envelopes based on tr observations
            gmin1=FSRenvmdr(tr,p,'prob',[0.99; 0.999; 0.01; 0.5],'init',init);
            
            for ii=(i-1):size(gmin1,1);
                
                % CHECK IF STOPPING RULE IS FULFILLED
                % ii>=size(gmin1,1)-2 = final, penultimate or antepenultimate value
                % of the resuperimposed envelope based on tr observations
                if mdr(ii,2)>gmin1(ii,c99) && ii>=size(gmin1,1)-2;
                    % Condition S1
                    mes=['$r_{min}('   int2str(mdr(ii,1)) ',' int2str(tr) ')>99$\% envelope'];
                    if msg
                        disp(['Superimposition stopped because r_{min}(' int2str(mdr(ii,1)) ',' int2str(tr) ')>99% envelope']);
                        disp(mes);
                    end
                    sto=1;
                    break
                elseif ii<size(gmin1,1)-2 &&  mdr(ii,2)>gmin1(ii,c999);
                    % Condition S2
                    mes=['$r_{min}('   int2str(mdr(ii,1)) ',' int2str(tr) ')>99.9$\% envelope'];
                    if msg
                        disp(['Superimposition stopped because r_{min}(' int2str(mdr(ii,1)) ',' int2str(tr) ')>99.9% envelope']);
                    end
                    sto=1;
                    break
                else
                    % mdr is inside the envelopes, so keep resuperimposing
                end
            end
            
            if plo==2
                % Plotting part
                subplot(nr,nc,jwind,'Parent',figure2);
                ylim([yl1 yl2]); xlim([xl1 xl2]);
                box('on'); hold('on');
                
                % Show curve of mdr up to step tr-1 (notice that the envelope is
                % based on tr observations. Step tr-1 in matrix mdr is
                % (tr-1)-mdr(1,1)+1=tr-mdr(1,1)
                plot(mdr(1:(tr-mdr(1,1)),1),mdr(1:(tr-mdr(1,1)),2));
                
                % Display the lines associated with 1%, 99% and 99.9% envelopes
                line(gmin1(:,1),gmin1(:,[2 3 4]),'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
                line(gmin1(:,1),gmin1(:,5),'LineWidth',2,'LineStyle','--','Color',[0.3 0.3 0.2]);
                
                strtemp=['$r_{min}(m,' int2str(tr) ')$'];
                % get the position of the current pane
                gposcurax=get(gca,'position');
                % set the width and the height of the current pane but keep
                % unaltered the distance from bottom left corner
                set(gca,'position',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3)*1.3,gposcurax(4)*1.3]);
                % Subplots located on the left hand side
                plleft=1:nc:(nr*(nc-1)+1);
                % Subplots located on the right hand side
                plright=nc:nc:nr*nc;
                % For all plots not located on the left and the right hand side
                % delete numbers on the y axis (YTickLabels)
                getYTickLab=get(gca,'YTickLabel');
                if isempty(intersect(jwind,[plleft plright])) && nr*nc>1;
                    set(gca,'YTickLabel',[]);
                end
                % For all plots not located on the right hand side put the
                % yaxis location on the right
                if ~isempty(intersect(jwind,plright)) && nr*nc>1  && nc>1;
                    set(gca,'YAxisLocation','right');
                end
                % For all plots which are not on the bottom hand side delete
                % numbers on the x axis (XTickLabels)
                if isempty(intersect(jwind,(nr-1)*nc+1:(nr*nc)));
                    set(gca,'XTickLabel',[]);
                else
                    % For the plots on the bottom side add the xlabel
                    xlabel('Subset size m');
                end
                annotation(figure2,...
                    'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)+0.05],...
                    'String',{strtemp},'Tag','mes1',...
                    PrVaCell{:});
                
                hold('off');
                jwind=jwind+1;
            end
            
            if sto==1;
                if plo==2
                    % Write on the plot the reason why the procedure stopped
                    annotation(figure2,...
                        'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)],...
                        'String',{mes},'Tag','mes2',...
                        PrVaCell{:});
                    
                    % Unless for all plots not located on the right hand side
                    % For the final plot put the yaxis location on the right
                    % Unless it is the first on the left hand side
                    if isempty(intersect(jwind-1,plleft)) && nr*nc>1  && nc>1;
                        set(gca,'YTickLabel',getYTickLab);
                        set(gca,'YAxisLocation','right');
                    end
                    
                    % add X label again to the last plot
                    set(gca,'XTickLabel',get(get(figure1,'Children'),'XTick'))
                    xlabel('Subset size m');
                    
                    % if jwind=2 than the width of the last plot is enlarged to
                    % full screen
                    if jwind==2 && nr*nc>1;
                        set(gca,'Position',[0.1 0.1 0.85 0.85])
                        % Relocate the two messages
                        set(findall(gcf,'Tag','mes1'),'Units','normalized','Position',[0.3 0.7 0.5 0.1])
                        set(findall(gcf,'Tag','mes2'),'Units','normalized','Position',[0.3 0.6 0.5 0.1])
                    end
                end
                
                break
            end
            if plo==2
                if jwind==nr*nc+1;
                    jwind=1;
                    figure2 = figure('PaperSize',[20.98 29.68],'Name',['Resuperimposed envelopes #' int2str(resup)]);
                    resup=resup+1;
                    % Create axes (Following line should not be necessary
                    %axes('Parent',figure2);
                end
            end
            
        end
        
        %% Stage 2a: subset validation
        % In this part we check whether the subset is homogeneous. In other
        % words we verify conditions H1 and H2
        % tr= m^\dagger+k+1
        % m^\dagger+k=tr-1
        % m*=mdr(ii,1)
        % Condition H2
        % Check if stopping rule takes place at m* <m^\dagger+k
        if (mdr(ii,1)<tr-1);
            % Condition H2b and H2a
            if sum(gmin1(ii+1:end,4)>mdr(ii+1:size(gmin1,1),2))>0;
                if msg
                    disp(['Subsample of ' int2str(tr-1) ' units is not homogeneous because the curve was above 99.99% and later it was below 1%']);
                    disp('----------------------------------------');
                end
                % Find m^{1%} that is the step where mdr goes below the 1%
                % threshold for the first time
                % ginfd = concatenate all the steps from m^*+1 to m^\dagger+k-1
                gfind=[gmin1(i+1:end,1) gmin1(i+1:end,4)>mdr(i+1:size(gmin1,1),2)];
                % select from gfind the steps where mdr was below 1% threshold
                % gfind(1,1) contains the first step where mdr was below 1%
                gfind=gfind(gfind(:,2)>0,1);
                % find maximum in the interval m^\dagger=mdr(i,1) to the step
                % prior to the one where mdr goes below 1% envelope
                if length(gfind)==1
                    tr=gfind;
                else
                    tr=sortrows(mdr(i:gfind(1,1)-mdr(1,1),1:2),2);
                    tr=tr(end,1);
                end
                if msg
                    disp('Probably there are two overlapping groups');
                    disp(['Using the criterion of the maximum, the group of homogenous obs. is=' int2str(tr)]);
                end
                % tr is redefined and is associated with the step associated to
                % the maximum value of r_min
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
        
    elseif  isempty(bonflev) && signal==2 % exceedance of the lower threshold
        
        
        if plo==2
            nr=3; nc=3;
            resup=1;
            jwind=1;
            figure2 = figure('PaperSize',[20.98 29.68],'Name',['Resuperimposed envelopes #' int2str(resup)]);
            
        end
        
        
        
        for tr=(mdag-1):-1:max(mdag-80,init+1);
            % Compute theoretical envelopes based on tr observations
            gmin1=FSRenvmdr(tr,p,'prob',[0.99; 0.01; 0.5; 0.001],'init',init);
            
            ii=size(gmin1,1);
            
            % CHECK IF STOPPING RULE IS FULFILLED
            % resuperimposed envelope based on tr observations is greater
            % than 0.1% threshold
            if mdr(ii,2)>gmin1(ii,5);
                % Condition N1
                mes=['$r_{min}('   int2str(mdr(ii,1)) ',' int2str(tr) ')>0.1$\% envelope'];
                if msg
                    disp(['Superimposition stopped because r_{min}(' int2str(mdr(ii,1)) ',' int2str(tr) ')>1% envelope']);
                    disp(mes);
                end
                sto=1;
                
            elseif tr==max(mdag-80,init+1);
                disp('Exceedance of the lower envelope in the intial part of the search')
                disp('Please decrease the value of init')
                % mdr is still below the 0.1% lower envelope, so keep resuperimposing
            else
            end
            
            
            if plo==2
                % Plotting part
                subplot(nr,nc,jwind,'Parent',figure2);
                ylim([yl1 yl2]); xlim([xl1 xl2]);
                box('on'); hold('on');
                
                % Show curve of mdr up to step tr-1 (notice that the envelope is
                % based on tr observations. Step tr-1 in matrix mdr is
                % (tr-1)-mdr(1,1)+1=tr-mdr(1,1)
                plot(mdr(1:(tr-mdr(1,1)),1),mdr(1:(tr-mdr(1,1)),2));
                
                % Display the lines associated with 1%, 50% and 99% envelopes
                line(gmin1(:,1),gmin1(:,2:4),'LineWidth',2,'LineStyle','--','Color',[0 0 1]);
                
                % Display the lines associated with 0.1% envelopes
                line(gmin1(:,1),gmin1(:,5),'LineWidth',2,'LineStyle','--','Color',[1 0 0]);
                
                
                strtemp=['$r_{min}(m,' int2str(tr) ')$'];
                % get the position of the current pane
                gposcurax=get(gca,'position');
                % set the width and the height of the current pane but keep
                % unaltered the distance from bottom left corner
                set(gca,'position',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3)*1.3,gposcurax(4)*1.3]);
                % Subplots located on the left hand side
                plleft=1:nc:(nr*(nc-1)+1);
                % Subplots located on the right hand side
                plright=nc:nc:nr*nc;
                % For all plots not located on the left and the right hand side
                % delete numbers on the y axis (YTickLabels)
                getYTickLab=get(gca,'YTickLabel');
                if isempty(intersect(jwind,[plleft plright])) && nr*nc>1;
                    set(gca,'YTickLabel',[]);
                end
                % For all plots not located on the right hand side put the
                % yaxis location on the right
                if ~isempty(intersect(jwind,plright)) && nr*nc>1  && nc>1;
                    set(gca,'YAxisLocation','right');
                end
                % For all plots which are not on the bottom hand side delete
                % numbers on the x axis (XTickLabels)
                if isempty(intersect(jwind,(nr-1)*nc+1:(nr*nc)));
                    set(gca,'XTickLabel',[]);
                else
                    % For the plots on the bottom side add the xlabel
                    xlabel('Subset size m');
                end
                annotation(figure2,...
                    'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)+0.05],...
                    'String',{strtemp},'Tag','mes1',...
                    PrVaCell{:});
                
                hold('off');
                jwind=jwind+1;
            end
            
            if sto==1;
                if plo==2
                    % Write on the plot the reason why the procedure stopped
                    annotation(figure2,...
                        'textbox',[gposcurax(1)-0.1,gposcurax(2),gposcurax(3),gposcurax(4)],...
                        'String',{mes},'Tag','mes2',...
                        PrVaCell{:});
                    
                    % Unless for all plots not located on the right hand side
                    % For the final plot put the yaxis location on the right
                    % Unless it is the first on the left hand side
                    if isempty(intersect(jwind-1,plleft)) && nr*nc>1  && nc>1;
                        set(gca,'YTickLabel',getYTickLab);
                        set(gca,'YAxisLocation','right');
                    end
                    
                    % add X label again to the last plot
                    set(gca,'XTickLabel',get(get(figure1,'Children'),'XTick'))
                    xlabel('Subset size m');
                    
                    % if jwind=2 than the width of the last plot is enlarged to
                    % full screen
                    if jwind==2 && nr*nc>1;
                        set(gca,'Position',[0.1 0.1 0.85 0.85])
                        % Relocate the two messages
                        set(findall(gcf,'Tag','mes1'),'Units','normalized','Position',[0.3 0.7 0.5 0.1])
                        set(findall(gcf,'Tag','mes2'),'Units','normalized','Position',[0.3 0.6 0.5 0.1])
                    end
                end
                
                break
            end
            if plo==2
                if jwind==nr*nc+1;
                    jwind=1;
                    resup=resup+1;
                    figure2 = figure('PaperSize',[20.98 29.68],'Name',['Resuperimposed envelopes #' int2str(resup)]);
                    % Create axes (Following line should not be necessary
                    %axes('Parent',figure2);
                end
            end
            
        end
        
        ndecl=n-tr;
        
        
    else
        ndecl=n-mdr(i,1);
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
    disp('Summary of the exceedances');
    disp(nout);
end

if msg
    if ~isempty(extram3);
        disp(extram3);
    else
    end
    
    if ~isempty(extram2);
        disp(extram2);
    else
    end
end

group=ones(n,1);

% Plot entry order of the untis
% plot([Un(1,1)-1;Un(:,1)],bb','x')

if ndecl>0;
    % Now find the list of the units declared as outliers
    % bsel=~isnan(bb(:,tr-init+1));
    % ListOut=setdiff(1:n,bsel,1);
    ListOut=seq(isnan(bb(:,end-ndecl)));
    
    % Add to ListOut all the units which have equal values in terms of X
    % and to y to those declared as outliers
    add=zeros(round(n*5),1);
    good=setdiff(seq,ListOut);
    Xy=[X y];
    ij=0;
    for i=1:length(ListOut)
        for j=1:length(good)
            if isequal(Xy(good(j),:),Xy(ListOut(i),:))
                ij=ij+1;
                add(ij)=good(j);
            end
            %   disp(['i' num2str(i) 'j' num2str(j)])
        end
    end
    if ij>0
        add=add(1:ij);
        add=unique(add);
        ListOut=[ListOut,add'];
    end
    ndecl=length(ListOut);
    group(ListOut)=2;
else
    ListOut=NaN;
end

%% Scatter plot matrix with the outliers shown with a different symbol

if plo==1 || plo==2
    figure;
    if isempty(options.namey)
        namey=char('y');
    else
        namey=options.namey;
    end
    
    if intercept==1
        if isempty(options.nameX)
            nameX=cellstr(num2str((1:p-1)','X%d'));
        else
            nameX=options.nameX;
        end
        [H,~,~]=gplotmatrix(X(:,2:end),y,group,'br','+o',[],[],[],nameX,namey);
    else
        if isempty(options.nameX)
            nameX=cellstr(num2str((1:p)','X%d'));
        else
            nameX=options.nameX;
        end
        [H,~,~]=gplotmatrix(X,y,group,'br','+o',[],[],[],nameX,namey);
    end
    set(gcf,'Name','Scatter plot matrix y|X with outliers highlighted');
    
    if ndecl>0
        set(H(:,:,1),'DisplayName','Good units');
        set(H(:,:,2),'DisplayName','Outliers');
    else
        set(H,'DisplayName','Units');
    end
    
    % save the indices of the outliers (ListOut) to the
    % 'UserData' field of the second group of H(:,:,2)
    if ~isnan(ListOut)
        set(H(:,:,2), 'UserData' , ListOut');
    end
    
    % The following line adds objects to the panels of the yX
    % add2yX(H,AX,BigAx,outadd,group,ListOut,bivarfit,multivarfit,labeladd)
    add2yX('intercept',intercept,'bivarfit',bivarfit,'multivarfit',multivarfit,'labeladd',labeladd);
end

%% Structure returned by function FSR
out=struct;
out.ListOut=ListOut;

% If you wish that the output also contains the list of units not declared
% as outliers, please uncomment the two following lines.
% ListIn=seq(~isnan(bb(:,end-ndecl)));
% out.ListIn=ListIn;

out.mdr=mdr;
out.Un=Un;
out.nout=nout;
out.Hetero=Hetero;
out.S2=S2;
out.WEI=WEI;

%% Callback functions used to "pin" quantile labels and vertical line to axes.

    function fprecallback(obj,~)
        % When the panning action starts, the quantile labels and the
        % vertical line which identifies the final part of the search are
        % made invisible.
        
        hanno = findall(obj,'Tag', 'quantile_label');
        set(hanno,'Visible','off');
        hanno2 = findall(obj, 'Tag' , 'FinalPartLine');
        set(hanno2,'Visible','off');
        
    end

    function fpostcallback(obj,evd)
        %When the panning action terminates, the new positions of the
        %labels and of the vertical line are set.
        
        % axis limits
        xlim = get(evd.Axes,'XLim'); xmin = xlim(1); xmax=xlim(2);
        
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
        % the figure coordinates which correspond to the anotations positions.
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

