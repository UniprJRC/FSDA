function [out] = FSRedaCens(y,X,bsb,varargin)
%FSRedaCens enables to monitor several quantities in each step of the forward search
%
%<a href="matlab: docsearchFS('FSRedaCens')">Link to the help function</a>
%
% Required input arguments:
%
%    y:         Response variable. Vector. Response variable, specified as
%               a vector of length n, where n is the number of
%               observations. Each entry in y is the response for the
%               corresponding row of X.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%   X :         Data matrix of explanatory variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables. Missing values (NaN's) and
%               infinite values (Inf's) are allowed, since observations
%               (rows) with missing or infinite values will automatically
%               be excluded from the computations.
%   bsb :       list of units forming the initial
%               subset. Vector or scalar. If bsb=0 (default), then the procedure starts with p
%               units randomly chosen, else if bsb is not 0, the search will
%               start with m0=length(bsb).
%
% Optional input arguments:
%
% balancedSearch:   Balanced search. Scalar logical.
%                   If Balanced search the proportion of observations in
%                   the subsets equals (as much as possible) the proportion
%                   of units in the sample. The default value of
%                   balancedSearch is true.
%                   Example - 'balancedSearch',false
%                   Data Types - logical
%
%  conflev:   confidence levels to be used to compute confidence interval
%             for the elements of $\beta$ and for $\sigma^2$. Vector.
%             The default value of conflev is [0.95 0.99] that
%             is 95% and 99% confidence intervals are computed.
%             Example - 'conflev',[0.90 0.93]
%             Data Types - double
%
%
% init :      Search initialization. Scalar.
%             It specifies the point where to initialize the search
%             and start monitoring required diagnostics. If init is not
%             specified it will be set equal to :
%             p+1, if the sample size is smaller than 40;
%             min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%             Example - 'init',100 starts monitoring from step m=100
%             Data Types - double
%
%    left :     left limit for the censored dependent variable. Scalar.
%               If set to -Inf, the dependent variable is assumed to be not
%               left-censored; default value of left is zero (classical
%               Tobit model).
%               Example - 'left',1
%               Data Types - double
%
% intercept : Indicator for constant term. true (default) | false.
%             Indicator for the constant term (intercept) in the fit,
%             specified as the comma-separated pair consisting of
%             'Intercept' and either true to include or false to remove
%             the constant term from the model.
%             Example - 'intercept',false
%             Data Types - boolean
%
%  nocheck:  Check input arguments. Boolean.
%            If nocheck is equal to true, no check is performed on
%            matrix y and matrix X. Notice that y and X are left
%            unchanged. In other words the additional column of ones for
%            the intercept is not added. As default nocheck=false. The
%            controls on h, alpha and nsamp still remain
%            Example - 'nocheck',true
%            Data Types - boolean
%
%    right :    right limit for the censored dependent variable. Scalar.
%               If set to Inf, the dependent variable is assumed to be not
%               right-censored; default value of right is Inf (classical
%               Tobit model).
%               Example - 'right',800
%               Data Types - double
%
%  tstat:    the kind of t-statistics which have to be monitored.
%            Character.
%            tstat = 'trad' implies  monitoring of traditional t
%            statistics (out.Tols). In this case the estimate of $\sigma^2$ at step m
%            is based on $s^2_m$ (notice that $s^2_m<<\sigma^2$ when m/n is
%            small) tstat = 'scal' (default) implies monitoring of
%            rescaled t statistics In this case the estimate of
%            $\sigma^2$ at step m is based on $s^2_m / var_{truncnorm(m/n)}$
%            where $var_{truncnorm(m/n)}$ is the variance of the truncated
%            normal distribution.
%            Example - 'tstat','trad'
%            Data Types - char
%
% Remark:       The user should only give the input arguments that have to
%               change their default value. The name of the input arguments
%               needs to be followed by their value. The order of the input
%               arguments is of no importance.
%
%               Missing values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations. y can
%               be both a row of column vector.
%
% Output:
%
%         out:   structure which contains the following fields
%
%   out.RES=        n x (n-init+1) = matrix containing the monitoring of
%               scaled residuals:
%               1st row = residual for first unit;
%               ...;
%               nth row = residual for nth unit.
%   out.LEV=        (n+1) x (n-init+1) = matrix containing the monitoring of
%               leverage:
%               1st row = leverage for first unit:
%               ...;
%               nth row = leverage for nth unit.
%    out.BB=        n x (n-init+1) matrix containing the information about the units belonging
%               to the subset at each step of the forward search:
%               1st col = indexes of the units forming subset in the
%               initial step;
%               ...;
%               last column = units forming subset in the final step (all
%               units).
%   out.mdr=        n-init x 3 matrix which contains the monitoring of minimum
%               deletion residual or (m+1)ordered residual at each step of
%               the forward search:
%               1st col = fwd search index (from init to n-1);
%               2nd col = minimum deletion residual;
%               3rd col = (m+1)-ordered residual.
%               Remark: these quantities are stored with sign, that is the
%               min deletion residual is stored with negative sign if
%               it corresponds to a negative residual.
%   out.msr=        n-init+1 x 3 = matrix which contains the monitoring of
%               maximum studentized residual or m-th ordered residual:
%               1st col = fwd search index (from init to n);
%               2nd col = maximum studentized residual;
%               3rd col = (m)-ordered studentized residual.
%   out.nor=        (n-init+1) x 4 matrix containing the monitoring of
%               normality test in each step of the forward search:
%               1st col = fwd search index (from init to n);
%               2nd col = Asymmetry test;
%               3rd col = Kurtosis test;
%               4th col = Normality test.
%  out.Bols=        (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward
%               search.
%    out.S2=       (n-init+1) x 5 matrix containing the monitoring of S2 or
%                   R2, F test, in each step of the forward search:
%               1st col = fwd search index (from init to n);
%               2nd col = monitoring of S2;
%               3rd col = monitoring of R2;
%               4th col = monitoring of rescaled S2.
%               In this case the
%               estimated of $\sigma^2$ at step m is divided by the
%               consistency factor (to make the estimate asymptotically
%               unbiased).
%               5th col = monitoring of F test. Note that an asymptotic
%               unbiased estimate of sigma2 is used.
%   out.coo=    (n-init+1) x 3 matrix containing the monitoring of Cook or
%               modified Cook distance in each step of the forward search:
%               1st col = fwd search index (from init to n);
%               2nd col = monitoring of Cook distance;
%               3rd col = monitoring of modified Cook distance.
%  out.Tols=    (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated t-statistics (as specified in option input
%               'tstat'.
%               in each step of the forward search
%   out.Un=        (n-init) x 11 Matrix which contains the unit(s)
%               included in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one. Un(1,2), for example, contains
%               the unit included in step init+1. Un(end,2) contains the
%               units included in the final step of the search.
%  out.betaINT = Confidence intervals for the elements of $\beta$.
%                 betaINT is a (n-init+1)-by-2*length(confint)-by-p 3D array.
%                 Each third dimension refers to an element of beta:
%                 betaINT(:,:,1) is associated with first element of beta;
%                 ...;
%                 betaINT(:,:,p) is associated with last element of beta.
%                 The first two columns contain the lower
%                 and upper confidence limits associated with conflev(1).
%                 Columns three and four contain the lower
%                 and upper confidence limits associated with conflev(2);
%                 ...;
%                 The last two columns contain the lower
%                 and upper confidence limits associated with conflev(end).
%
%                 For example, betaint(:,3:4,5) contain the lower and upper
%                 confidence limits for the fifth element of beta using
%                 confidence level specified in the second element of input
%                 option conflev.
%out.sigma2INT = confidence interval for $\sigma^2$.
%                1st col = fwd search index;
%                2nd col = lower confidence limit based on conflev(1);
%                3rd col = upper confidence limit based on conflev(1);
%                4th col = lower confidence limit based on conflev(2);
%                5th col = upper confidence limit based on conflev(2);
%                ...
%                penultimate col = lower confidence limit based on conflev(end);
%                last col = upper confidence limit based on conflev(end);
%     out.y=     A vector with n elements that contains the response
%               variable which has been used
%     out.X=    Data matrix of explanatory variables
%               which has been used (it also contains the column of ones if
%               input option intercept was missing or equal to 1)
%  out.class =  'FSReda'.
%
%
%
% See also LXS.m, FSReda.m, FSRfanCens
%
% References:
%
% Atkinson, A.C. and Riani, M. (2000), "Robust Diagnostic Regression
% Analysis", Springer Verlag, New York.
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRedaCens')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % FSRedaCens with all default options.
    % Example of use of FSRedaCens based on a starting point coming
    % from LMS.
    n=200;
    p=3;
    rng default
    rng(123456)
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1)+1;
    % Contaminated data
    ycont=y;
    cont=1:5;
    ycont(cont)=ycont(cont)+5;
    ycont(ycont<=0)=0;
    [out]=LXS(ycont,X,'nsamp',1000);
    out=FSRedaCens(ycont,X,out.bs);
    fground=struct;
    fground.funit=cont;
    resfwdplot(out,'fground',fground)
%}

%{
    % FSRedaCens with optional argument.
    % Example of use of function FSReda using a random start and traditional
    % t-stat monitoring.
    n=200;
    p=3;
    rng default
    rng(123456)
    X=randn(n,p);
    % Uncontaminated data
    y=randn(n,1);
    % Contaminated data
    ycont=y;
    ycont(1:5)=ycont(1:5)+6;
    ycont(ycont<=0)=0;
    out=FSRedaCens(ycont,X,0,'tstat','trad');
%}

%{
    %% Monitoring of residuals using the affairs dataset.
    % In the example of Kleiber and Zeileis (2008, p. 142), the number of a
    % person's extramarital sexual inter-courses ("affairs") in the past year
    % is regressed on the person's age, number of years married, religiousness,
    % occupation, and won rating of the marriage. The dependent variable is
    % left-censored at zero and not right-censored. Hence this is a standard
    % Tobit model which can be estimated by the following lines
    load affairs.mat
    X=affairs{:,["age" "yearsmarried" "religiousness" "occupation" "rating"]};
    y=affairs{:,"affairs"};
    outLXS=LXS(y,X);
    [~,sor]=sort(abs(outLXS.residuals))
    out=FSRedaCens(y,X,sor(1:100));
    resfwdplot(out)
%}

%{
    %%  Outliers and a Lower Threshold example.
    rng default
    rng(2)
    n=300;
    lambda=-0.5;
    p=5;
    sigma=0.1;
    beta=1*ones(p,1);
    X=0.2*randn(n,p);
    epsilon=randn(n,1);
    
    y=X*beta+sigma*epsilon;
    y=normYJ(y,1,lambda,'inverse',true,'Jacobian',false);
    
    sel=1:30;
    y(sel)=y(sel)+1.2;
    
    qq=quantile(y,0.3);
    y(y<=qq)=qq;
    left=min(y);
    right=Inf;
    
    % See function FSRfanCens on the procedure to find the correct
    % transformation
    yf=normYJ(y,1,lambda,'inverse',false,'Jacobian',false);
    leftf=normYJ(left,1,lambda,'inverse',false,'Jacobian',false);
    rightf=normYJ(right,1,lambda,'inverse',false,'Jacobian',false);
    
    zlimits=[leftf rightf];
    % Call to FSRedaCens
    outLXS=LXS(yf,X);
    out=FSRedaCens(yf,X,outLXS.bs,'left',leftf,'right',rightf,'init',100);
    fground.funit=1:30;
    resfwdplot(out,'fground',fground);
%}

%% Beginning of code

if nargin<2
    error('FSDA:FSRedaCens:missingInputs','not enough input args')
end


% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = aux.chkinputR(y,X,nnargin,vvarargin);

%% User options
init=round(n*0.75);

conflevdef=[0.95 0.99];
left=0;
right=Inf;
balancedSearch=true;

options=struct('intercept',true,'init',init,'tstat','scal',...
    'nocheck',false,'conflev',conflevdef, ...
    'balancedSearch',balancedSearch, 'left',left,'right',right);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRedaCens:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    aux.chkoptions(options,UserOptions)
end

if nargin > 3

    % We now overwrite inside structure options the default values with
    % those chosen by the user
    % Notice that in order to do this we use dynamic field names
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

if bsb==0
    Ra=1; nwhile=0;
    while and(Ra,nwhile<100)
        bsb=randsample(n,p);
        Xb=X(bsb,:);
        Ra=~(rank(Xb)==p);
        nwhile=nwhile+1;
    end
    if nwhile==100
        warning('FSDA:FSRedaCens:NoFullRank','Unable to randomly sample full rank matrix');
    end
    yb=y(bsb);
else
    Xb=X(bsb,:);
    yb=y(bsb);
end

ini0=length(bsb);

nocheck=options.nocheck;
left=options.left;
right=options.right;
balancedSearch=options.balancedSearch;

% check init
init=options.init;
if  init <p
    fprintf(['Attention : init should be larger than p+1. \n',...
        'It is set to p.']);
    init=p;
elseif init<ini0
    fprintf(['Attention : init should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    init=ini0;
elseif init>=n
    fprintf(['Attention : init should be smaller than n. \n',...
        'It is set to n-1.']);
    init=n-1;
end

%% Declare matrices to store quantities

% sequence from 1 to n
seq=(1:n)';

% complementary of bsb
ncl=setdiff(seq,bsb);

% The second column of matrix R will contain the OLS residuals
% at each step of the forward search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100=100*(1:1:ceil(n/100));

zer=NaN(n-init,2);
zer1=NaN(n-init+1,2);

% Matrix Bols will contain the beta coeff. in each step of the fwd search
% The first row will contain the units forming initial subset
Bols=[(init:n)' NaN(n-init+1,p)];

% Vector of the beta coefficients from the last correctly calculated step
% Used in case the rank of Xb is less than p
blast=NaN(p,1);

% S2=(n-init+1) x 3  matrix which will contain
% 1st col = fwd search index
% 2nd col = S2= \sum e_i^2 / (m-p)
% 3rd col = R^2
% 4th col = (\sum e_i^2 / (m-p)) / (consistency factor) to make the
% estimate asymptotically unbiased
% 5th col = F test (based on unbiased estimate of sigma2)
S2=[(init:n)' NaN(n-init+1,4)];

% mdr= (n-init) x 3 matrix
% 1st column = fwd search index
% 2nd column deletion residual among observations non belonging to the
% subset
% 3rd column (m+1)-th ordered residual
% They are stored with sign, that is the min deletion residual
% is stored with negative sign if it corresponds to a negative residual
mdr=[(init:n-1)'  zer];

% mdr= (n-init+1) x 3 matrix which will contain max studentized residual
%  among bsb and m-th studentized residual
msr=[(init:n)'  zer1];

% Coo= (n-init) x 3 matrix which will contain Cook distances
%  (2nd col) and modified Cook distance (3rd col)
coo=[((init+1):n)'  NaN(n-init,6)];

% nor= (n-init+1) x 3 matrix which will contain asymmetry (2nd col)
% kurtosis (3rd col) and normality test (4th col)
nor=[(init:n)'  zer1];

% Matrix RES will contain the residuals for each unit in each step of the forward search
% The first row refers to the residuals of the first unit
RES=NaN(n,n-init+1);
RES(:)=NaN;

% Matrix BB will contain the units forming subset in each step of the forward search
% The first column contains the units forming subset at step init
% The first row is associated with the first unit
BB=RES;

% Matrix LEV will contain the leverage of the units forming subset in each step of the forward search
% The first column contains the leverage associated with the units forming subset at step init
% The first row is associated with the first unit
LEV=RES;

%  Un= Matrix whose 2nd column:11th col contain the unit(s) just included
Un=NaN(n-init,10);
Un=[(init+1:n)' Un];

%  Tols = Matrix whose columns contain t statistics specified in option
%  tstat
Tols=Bols;



% vector conflev
conflev=options.conflev;
conflev=1-(1-conflev)/2;
lconflev=length(conflev);

% betaINT will contain the confidence intervals for the elements of $\beta$
% betaINT is a (n-init+1)-by-2*length(confint)-by-p 3D array.
% Each third dimension refers to an element of beta
% betaINT(:,:,1) is associated with first element of beta
% .....
% betaINT(:,:,p) is associated with last element of beta
% The first two columns contain the lower
% and upper confidence limits associated with conflev(1).
% Columns three and four contain the lower
% and upper confidence limits associated with conflev(2)
% ....
% The last two columns contain the lower
% and upper confidence limits associated with conflev(end)
betaINT=NaN(n-init+1,2*lconflev,p);
% sigma2INT confidence interval for $\sigma^2$
sigma2INT=[(init:n)' zeros(n-init+1,2*lconflev)];

% opts is a structure which contains the options to use in linsolve
opts=struct;
opts.RECT = true;
opts.LT =false;
opts.UT =false;

if balancedSearch==true
    yleft=seq(y==left);
    yright=seq(y==right);
    inleft=Inf(length(yleft),1);
    inright=Inf(length(yright),1);
    propleft=length(yleft)/n;
    propright=length(yright)/n;
end

%% Start of the forward search
if nocheck==false && rank(Xb)~=p
    warning('FSDA:FSReda:NoFullRank','Initial subset does not form full rank matrix');
    % FS loop will not be performed
else
    for mm=ini0:n

        % if n>200 show every 100 steps the fwd search index
        if n>200
            if isscalar(intersect(mm,seq100))
                disp(['m=' int2str(mm)]);
            end

        end

        if  mm< p+5

            % Implicitly control the rank of Xb checking the condition number
            % for inversion (which in the case of a rectangular matrix is
            % nothing but the rank)
            % Old instruction was b=Xb\yb;
            [b,condNumber]=linsolve(Xb,yb,opts);
            % disp([mm condNumber])
            if condNumber<p
                NoRankProblem =false;
            else
                NoRankProblem =true;
            end


            if NoRankProblem  % rank is ok
                resBSB=yb-Xb*b;
                blast=b;   % Store correctly computed b for the case of rank problem
            else   % number of independent columns is smaller than number of parameters
                warning('FSDA:FSReda','Rank problem in step %d: Beta coefficients are used from the most recent correctly computed step',mm);
                b=blast;
            end
        else
            NoRankProblem =true;
            outTOB=regressCens(yb,Xb,'left',left,'right',right,'intercept',false);
            % Store values of beta coefficients
            b=outTOB.Beta(1:end-1,1);
            resBSB=yb-Xb*b;
        end


        if (mm>=init)

            % Store Units belonging to the subset
            BB(bsb,mm-init+1)=bsb;

            if NoRankProblem

                % Store beta coefficients
                Bols(mm-init+1,2:p+1)=b';


                % Measure of asymmetry
                sqb1=(sum(resBSB.^3)/mm) / (sum(resBSB.^2)/mm)^(3/2);

                % Measure of Kurtosis  */
                b2=(sum(resBSB.^4)/mm) / (sum(resBSB.^2)/mm)^2;

                % Asymmetry test
                nor(mm-init+1,2)=  (mm/6)*  sqb1  ^2  ;

                % Kurtosis test
                nor(mm-init+1,3)=(mm/24)*((b2 -3)^2);

                % Normality test
                nor(mm-init+1,4)=nor(mm-init+1,2)+nor(mm-init+1,3);

                % Store leverage for the units belonging to subset
                % hi contains leverage for all units
                % It is a proper leverage for the units belonging to subset
                % It is a pseudo leverage for the units not belonging to the subset
                mAm=Xb'*Xb;

                mmX=inv(mAm);
                dmmX=diag(mmX);
                % Notice that we could replace the following line with
                % hi=sum((X/mAm).*X,2); but there is no gain since we need
                % to compute dmmX=diag(mmX);
                hi=sum((X*mmX).*X,2); %#ok<MINV>

                LEV(bsb,mm-init+1)=hi(bsb);
            end % no rank problem
        end

        if (mm>p)

            % store res. sum of squares/(mm-k)
            % Store estimate of \sigma^2 using units forming subset
            if NoRankProblem
                Sb=(resBSB)'*(resBSB)/(mm-p);
            end

        else
            Sb=0;
        end

        % e= vector of residual for all units using b estimated using subset
        yhat=X*b;
        e=y-yhat;

        if (mm>=init)
            % Store all residuals
            RES(:,mm-init+1)=e;

            if NoRankProblem
                % Store S2 for the units belonging to subset
                S2(mm-init+1,2)=Sb;

                % Store rescaled version of S2 in the fourth column
                % Compute the variance of the truncated normal distribution
                if mm<n
                    a=norminv(0.5*(1+mm/n));
                    corr=1-2*(n./mm).*a.*normpdf(a);
                else
                    corr=1;
                end
                Sbrescaled=Sb/corr;
                S2(mm-init+1,4)=Sbrescaled;

                % Store F test
                % Compute regression sum of squares (divided by p-1)
                yhatb=yhat(bsb);
                if options.intercept==1
                    yhatm=yhatb-sum(yhatb)/mm;
                else
                    yhatm=yhatb;
                end
                SSRadd=yhatm'*yhatm/(p-1);
                S2(mm-init+1,5)=SSRadd/S2(mm-init+1,4);


                % Store maximum studentized residual
                % among the units belonging to the subset
                msrsel=sort(abs(resBSB)./sqrt(Sb*(1-hi(bsb))));
                msr(mm-init+1,2)=msrsel(mm);

                % Store R2
                S2(mm-init+1,3)=1-var(resBSB)/var(yb);
            end

        end

        r(:,2)=e.^2;

        if mm>init

            if NoRankProblem

                % Store in the second column of matrix coo the Cook
                % distance
                bib=Bols(mm-init+1,2:p+1)-Bols(mm-init,2:p+1);
                if S2(mm-init+1,2)>0
                    coo(mm-init,2)=bib*mAm*(bib')/(p*S2(mm-init+1,2));
                end

                if length(unit)>5
                    unit=unit(1:5);
                end
                if S2(mm-init,2)>0
                    coo(mm-init,3:length(unit)+2)= 1./(1-hi(unit)).* sqrt(((mm-p)/p)*hi(unit).*r(unit,2)./S2(mm-init,2));
                end
            end % NoRankProblem
        end

        if mm<n
            if mm>=init
                if NoRankProblem
                    % ord = matrix whose first col (divided by S2(i)) contains the deletion residuals
                    % for all units. For the units belonging to the subset these are proper deletion residuals
                    ord = [(r(:,2)./(1+hi)) e];

                    % Store minimum deletion residual in 2nd col of matrix mdr
                    selmdr=sortrows(ord(ncl,:),1);
                    mdr(mm-init+1,2)=sign(selmdr(1,2))*sqrt(selmdr(1,1)/S2(mm-init+1,2));

                    % Store (m+1) ordered pseudodeletion residual in 3rd col of matrix
                    % mdr
                    selmdr=sortrows(ord,1);
                    mdr(mm-init+1,3)=sign(selmdr(mm+1,2))*sqrt(selmdr(mm+1,1)/S2(mm-init+1,2));
                end % NoRankProblem
            end

            % store units forming old subset in vector oldbsb
            oldbsb=bsb;

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
            %  forming the group of potential outliers
            % ord=sortrows(r,2);
            [~,ord]=sort(r(:,2));

            % bsb= units forming the new  subset
            bsb=ord(1:(mm+1),1);

            Xb=X(bsb,:);  % subset of X
            yb=y(bsb);    % subset of y

            if mm>=init
                unit=setdiff(bsb,oldbsb);
                if length(unit)<=10
                    Un(mm-init+1,2:(length(unit)+1))=unit;
                else
                    disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                    Un(mm-init+1,2:end)=unit(1:10);
                end
            end


            if mm < n-1
                % ncl= units forming the new noclean
                ncl=ord(mm+2:n,1);
            end
        end

        if mm >= init
            if NoRankProblem
                if strcmp(options.tstat,'scal')

                    Tols(mm-init+1,2:end)=sqrt(corr)*Bols(mm-init+1,2:end)./sqrt(Sb*dmmX');

                elseif strcmp(options.tstat,'trad')
                    Tols(mm-init+1,2:end)=Bols(mm-init+1,2:end)./sqrt(Sb*dmmX');
                end

                % Compute highest posterior density interval for each value of
                % Tinvcdf = required quantiles of T distribution
                % consider just upper quantiles due to symmetry
                Tinvcdf=tinv(conflev,mm-p);
                % IGinvcdf = required quantiles of Inverse Gamma distribution
                Chi2invcdf=chi2inv([conflev 1-conflev],mm-p);

                c=sqrt(Sbrescaled*dmmX);
                for j=1:lconflev
                    betaINT(mm-init+1,j*2-1:j*2,:)=[ b-Tinvcdf(j)*c  b+Tinvcdf(j)*c]';
                    sigma2INT(mm-init+1,(j*2):(j*2+1))=[Sbrescaled*(mm-p)/Chi2invcdf(j) Sbrescaled*(mm-p)/Chi2invcdf(j+lconflev)];
                end

            end % NoRankProblem
        end
    end
end   %Rank check

RES=RES/sqrt(S2(end,2));

%% Structure returned by function FSReda
out=struct;
out.RES=RES;
out.LEV=LEV;
out.BB=BB;
out.mdr=mdr;
out.msr=msr;
out.nor=nor;
out.Bols=Bols;
out.S2=S2;
out.coo=coo;
out.Tols=Tols;
out.Un=Un;
out.betaINT=betaINT;
out.sigma2INT=sigma2INT;
out.y=y;
out.X=X;
out.class='FSRedaCens';
end
%FScategory:REG-Regression
