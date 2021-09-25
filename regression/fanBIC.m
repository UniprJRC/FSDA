function [out]=fanBIC(outFSRfan,varargin)
%fanBIC uses the output of FSRfan to choose the best value of the transformation parameter in linear regression
%
%<a href="matlab: docsearchFS('fanBIC')">Link to the help function</a>
%
% Required input arguments:
%
%  outFSRfan :  Structure created with function FSRfan. Structure.
%               Structure containing the following fields
%outFSRfan.Score  =  (n-init) x length(la)+1 matrix:
%               1st col = fwd search index;
%               2nd col = value of the score test in each step
%               of the fwd search for la(1);
%               ...;
%               last col  =  value of the score test in each step
%               of the fwd search for la(end).
%       outFSRfan.la   =  vector containing the values of lambda for which FSRfan
%               was computed
%       outFSRfan.bs   =  matrix of size p x length(la) containing the units forming
%               the initial subset for each value of lambda.
%         outFSRfan.y  = a vector containing the response
%         outFSRfan.X  = a matrix containing the explanatory variables
%      Data Types - struct
%
% Optional input arguments:
%
%
%       conflev :   Confidence level. Scalar. Confidence level
%                   to evaluate the exceedances in hte fanplot.
%                   Default confidence level is 0.9999 that is signals are
%                   considered when there is an exceedance for confidence
%                   level for at least 3 consecutive times.
%                   Example - 'conflev',[0.999]
%                   Data Types - double
%
%       init    : Step to start monitoring exceedances. Scalar.
%                 It specifies the initial
%                 subset size to start monitoring exceedances of the
%                 fanplot. If init is not specified it set equal
%                 to round(n*0.6).
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
%
%        family :   string which identifies the family of transformations which
%                   must be used. Character. Possible values are 'BoxCox'
%                   (default), 'YJ', 'YJpn' or 'YJall'.
%                   The Box-Cox family of power transformations equals
%                   $(y^{\lambda}-1)/\lambda$ for $\lambda$ not equal to zero,
%                   and $\log(y)$ if $\lambda = 0$.
%                   The Yeo-Johnson (YJ) transformation is the Box-Cox
%                   transformation of $y+1$ for nonnegative values, and of
%                   $|y|+1$ with parameter $2-\lambda$ for $y$ negative.
%                   Remember that BoxCox can be used just
%                   if input y is positive. Yeo-Johnson family of
%                   transformations does not have this limitation.
%                   Example - 'family','YJ'
%                   Data Types - char
%
%      bonflev  : Signal to use to identify outliers. Scalar. Option to be
%                used if the distribution of the data is
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
%               Example - 'bonflev',0.99
%               Data Types - double
%
%       plots   :  Plot on the screen. Scalar.
%                   If plots=1 a three panel plot will be produced. The
%                   left panel contains the BIC for the various values of
%                   lambda, the right panel the index of agreement with
%                   MLE, while the bottom panel the fraction of
%                   observations in agreement with the different values of
%                   lambda.
%                   Example - 'plots',1
%                   Data Types - double
%
%       tag     :   Handle of the plot. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag pl_fanBIC. Notice that if the
%                   program finds a plot which has a tag equal to the one
%                   specified by the user, then the output of the new plot
%                   overwrites the existing one in the same window else a
%                   new window is created.
%                   Example - 'tag','pl_myfanBIC'
%                   Data Types - char
%
%
%
%  Output:
%
%
%         out:   structure which contains the following fields
%
% out.BIC  = length(la)-by-3 matrix containing in the first column the
%             values of lambda, in the second column the values of BIC and
%             in the third column the values of the "Agreement index". The
%             agreement index is the reciprocal of the mean of the absolute
%             values of the score test computed in the interval init:h. The
%             default value of init is n*0.6 (see input option init) and h
%             is the number of clean observations in agreement with a
%             particular transformation. h is contained in the third column
%             of out.mmstop. The value of the index is rescaled with the
%             variance of the truncated normal distribution, in order to
%             give more weight to the searches with larger values of h.
% out.mmstop  = length(la)-by-3 matrix containing in the first column the
%             values of lambda, in the second column the number of units in
%             agreement with the different values of lambda and in the
%             third column the number of units not declared as outliers in
%             the subsequent outlier detection procedure.
% out.BBla = n-by-length(la) matrix containing information about the
%               outlier(s) for each value of lambda.
%               If out.BBla(i,j)=0 means that unit i (i=1, 2, ...n) is not
%               in agrement with la(j) j=1, 2, ..., length(la).
%               If out.BBla(i,j)=1 means that unit i (i=1, 2, ...n) is
%               in agrement with la(j) j=1, 2, ..., length(la) but has been
%               declared as outlier in the subsequent outlier detection
%               procedure.
%               If out.BBla(i,j)=2 means that unit i (i=1, 2, ...n) is
%               in agrement with la(j) j=1, 2, ..., length(la) and has not been
%               declared as outlier in the subsequent outlier detection
%               procedure.
%    out.labest= scalar. Value of lambda associated with the largest BIC
%                value.
%
% See also: FSRfan, fanplot
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
% Atkinson, A.C. and Riani, M. (2002a), Tests in the fan plot for robust,
% diagnostic transformations in regression, "Chemometrics and Intelligent
% Laboratory Systems", Vol. 60, pp. 87-100.
% Atkinson, A.C. Riani, M., Corbellini A. (2019), The analysis of
% transformations for profit-and-loss data, Journal of the Royal
% Statistical Society, Series C, "Applied Statistics",
% https://doi.org/10.1111/rssc.12389
% Atkinson, A.C. Riani, M. and Corbellini A. (2020), The Box-Cox
% Transformation: Review and Extensions, "Statistical Science", in press.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('fanBIC')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% fanBIC with all default options.
    % load the wool data.
    XX=load('wool.txt');
    y=XX(:,end);
    X=XX(:,1:end-1);
    % FSRfan and fanplot with all default options
    [outFSR]=FSRfan(y,X,'msg',0);
    out=fanBIC(outFSR);
%}
%
%{
    %% BIC plot with optional arguments.
    % FSRfan and fanBIC with specified lambda.
    load('loyalty.txt');
    y=loyalty(:,4);
    X=loyalty(:,1:3);
    % la = vector containing the grid of values to use for the
    % transformation parameter
    la=-1:0.1:1;
    [outFSRfan]=FSRfan(y,X,'la',la,'msg',0,'plots',0);
    out=fanBIC(outFSRfan);
%}
%


%% Beginning of code

la=outFSRfan.la;
nla=length(la);
X=outFSRfan.X;
y=outFSRfan.y;
bs=outFSRfan.bs;
[n,pwithintercept]=size(X);

seq=1:n;
logn=log(n);
% Divide central part from final part of the search
istep = n-floor(13*sqrt(n/200));

init=floor(n*0.6);
conflev=0.9999;
bonflev='';
family='YJ';
plots=1;
tag='pl_fanBIC';

UserOptions=varargin(1:2:length(varargin));


if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:fanBIC:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    
    options=struct('plots',plots,'init',init,'conflev',conflev,...
        'bonflev',bonflev,...
        'tag',tag,'family',family);
    
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    
    plots=options.plots;
    family=options.family;
    conflev=options.conflev;
    tag=options.tag;
end

if strcmp(family,'BoxCox')
    BoxCox=1;
elseif strcmp(family,'YJ')
    BoxCox=0;
    
else
    warning('FSDA:fanBIC:WrongFamily','Transformation family which has been chosen is not supported')
    error('FSDA:fanBIC:WrongFamily','Supported values are BoxCox or YJ ')
end

Sco=outFSRfan.Score;
bool=Sco(:,1)>=init;
Sco=Sco(bool,:);

nmdr=size(Sco,1);
upperbound=conflev+(1-conflev)/2;

th=norminv(upperbound);

mmstop=[la(:) n*ones(length(la),2)];
BBla=zeros(n,nla);

% SumLogY is a necessary ingredient to find the Jacobian
nonnegs=y>=0;
SumLogY=sum( log(   (1 + abs(y)).^(2 * nonnegs - 1)) );

% Initialize matrix BIC
% 2nd col = value of BIC for each value of la
% 3rd col = agreement index (the larger the better)
BIC=[la(:) zeros(nla,2)];

for j=1:nla
    Scolaj=abs(Sco(:,j+1));
    signal=false;
    for i=3:nmdr
        
        if i<istep-init+1 % CENTRAL PART OF THE SEARCH
            % Extreme triplet out of the envelope
            if Scolaj(i)>th && Scolaj(i-1)>th && Scolaj(i+1)>th
                mdagger=Sco(i-2,1);
                signal=true;
                break
            end
        else
            if  Scolaj(i)>th
                mdagger=Sco(i-1,1);
                signal=true;
                break
            end
        end
        % Now get out of the loop and do outlier detection
    end
    
    if BoxCox == 1
        ytraj=normBoxCox(y,1,la(j));
    else
        % Transform the data using Yeo and Johnson family or BoxCox family
        ytraj=normYJ(y,1,la(j),'Jacobian',false);
    end
    
    if signal==true
        mmstop(j,2)=mdagger;
        % Find units belonging to subset for laj in step mdagger
        [~,BB] = FSRbsb(ytraj,X,bs(:,j),'intercept',true,'init',mdagger,'nocheck',1,'msg',0);
        % good = good units for score test
        good=seq(~isnan(BB(:,1)));
    else
        good=seq;
    end
    n1=length(good);
    ini=min([round(n1*0.8) n1-6]);
    
    % Subsequent outlier detection
    [~,newbsbj]=intersect(good,bs(:,j));

    if length(newbsbj)<size(bs,1)
       newbsbj=1;
    end
    outl=FSR(ytraj(good),X(good,:),'nocheck',1,'bonflev',bonflev,'plots',0,'msg',0,'init',ini,'lms',newbsbj);
    if isscalar(outl.mdr)
         outl=FSR(ytraj(good),X(good,:),'nocheck',1,'bonflev',bonflev,'plots',0,'msg',0,'init',ini);
    end

    LowerEnv=FSRenvmdr(n1,size(X,2)+1,'init',ini);
    % If the value of mdr is always below the lower envelope
    if sum(outl.mdr(:,2)>LowerEnv(:,2)) ==0
        goodf=1:ini;
    else
        
        goodf=setdiff(1:length(good),outl.ListOut);
    end
    BBla(good,j)=1;
    BBla(good(goodf),j)=2;
    hh=length(goodf);
    mmstop(j,3)=hh;
    
    
    % logJ = log of the Jacobian
    logJ=(la(j)-1)*SumLogY;
    yb=ytraj(BBla(:,j)==2);
    Xb=X(BBla(:,j)==2,:);
    b=Xb\yb;
    e=yb-Xb*b;
    sigma2hat=(e'*e)/hh;
    
    vt = norminv(0.5*(1+hh/n));
    if hh<n
        factor = 1/(1-2*(n/hh)*vt.*normpdf(vt));
    else
        factor=1;
    end
    
    BIC(j,2)=2*(-0.5*n*log(factor*sigma2hat)+logJ)-(pwithintercept+1+n-hh)*logn;
    
    % Compute the "Agreement index"
    boo=Sco(:,1)<=hh;
    BIC(j,3)=1/(mean(Scolaj(boo))*factor);
end


 
% Find best value of lambda according to BIC
[~,imax]=max(BIC(:,2));
labest=la(imax);

% Find best value of lambda according to "Agreement index"
[~,imaxSI]=max(BIC(:,3));

% If the two indexes produce different answers delete the values of BIC and
% agreement index for which h (number of units in agreement with
% transformation and without outliers) is equal to init (m_M)
if imax~=imaxSI
    BIC(mmstop(:,2)<=init,2:3)=NaN;
    [~,imax]=max(BIC(:,2));
    labest=la(imax);
end

if plots == 1
    % Specify where to send the output of the current procedure
    h=findobj('-depth',1,'tag',tag);
    if (~isempty(h))
        clf(h);
        figure(h)
        axes;
    else
        figure;
    end
    set(gcf,'Name',['BIC for lambda=' mat2str(outFSRfan.la) ]);
    
    nr=4;
    nc=2;
    
    % Panel on the left
    subplot(nr,nc,[1 3 5])
    plot(BIC(:,1),BIC(:,2),'-ok')
    if length(la)<=6
        set(gca,'XTick',BIC(:,1));
    else
        set(gca,'XTick',BIC(1:2:length(la),1));
    end
    ylabel('BIC')
    xlabel('\lambda')
    hold('on')
    plot(BIC(imax,1),BIC(imax,2),'o','MarkerFaceColor','b')
    xlim([-0.1+BIC(1,1),BIC(end,1)+0.1])
    set(gca,'Xgrid','on')
    if length(la)>6
        set(gca,'XTickLabelRotation',90)
    end
    
    % Panel on the right
    subplot(nr,nc,[2 4 6])
    plot(BIC(:,1),BIC(:,3),'-ok')
    if length(la)<=6
        set(gca,'XTick',BIC(:,1));
    else
        set(gca,'XTick',BIC(1:2:length(la),1));
    end
    ylabel('Agreement  index')
    xlabel('\lambda')
    hold('on')
    
    plot(BIC(imaxSI,1),BIC(imaxSI,3),'o','MarkerFaceColor','b')
    xlim([-0.1+BIC(1,1),BIC(end,1)+0.1])
    set(gca,'Xgrid','on')
    if length(la)>6
        set(gca,'XTickLabelRotation',90)
    end
    
    % Bottom panel
    subplot(nr,nc,[7 8])
    hold('on')
    bar(la, mmstop(:,2)/n,'w')
    bar(la, mmstop(:,3)/n,'r')
    ylim([min(mmstop(:,3))/n-0.02 1])
    set(gca,'XTick',BIC(:,1));
    xlabel('\lambda')
    if length(la)>6
        set(gca,'XTickLabelRotation',90)
    end
    hold('off')
    
    
    % tag the figure
    set(gcf,'Tag',tag)
    
    titl=['Best \lambda='  num2str(la(imax)) '. Number of cleaned obs.='  num2str(mmstop(imax,3))];
    if verLessThanFS(9.5)
        suplabel(titl,'t');
    else
        sgtitle(titl)
    end
end

out=struct;
out.BIC=BIC;
out.mmstop=mmstop;
out.BBla=BBla;
out.labest=labest;
end
%FScategory:VIS-Reg