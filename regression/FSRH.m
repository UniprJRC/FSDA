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
%       scoring:  if scoring=1 a scoring algorithm is executed, if scoring=0,
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
    fs=12;
    ylabel('Value','FontSize',fs)
    xlabel('Quantity','FontSize',fs)
    set(gca,'FontSize',fs)
    
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
        error('FSDA:FSRH:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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

init=options.init;
h=options.h;
lms=options.lms;
nsamp=options.nsamp;
msg=options.msg;
bsbmfullrank=options.bsbmfullrank;
scoring=options.scoring;

%% Start of the forward search

seq=1:n;

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
Hetero=Hetero(bool,:);

%% Call core function which computes exceedances to thresholds of mdr
[out]=FSRcore(y,X,n,p,mdr,init,Un,bb,options);

out.mdr=mdr;
out.Un=Un;
out.Hetero=Hetero;
out.S2=S2;
out.WEI=WEI;

end

