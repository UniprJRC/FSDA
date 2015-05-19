function [out] = FSReda(y,X,bsb,varargin)
%FSReda enables to monitor several quantities in each step of the Bayeisan forward search
%
%<a href="matlab: docsearchFS('FSReda')">Link to the help function</a>
%
% Required input arguments:
%
%   y:          A vector with n elements that contains the response variables.
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
%   bsb :       vector containing the list of units forming the initial
%               subset, if bsb=0 (default) then the procedure starts with p
%               units randomly chosen else if bsb is not 0 the search will
%               start with m0=length(bsb)
%
% Optional input arguments:
%
%   intercept : Indicator for constant term. Scalar.
%                     If 1, a model with constant term will be fitted (default),
%                     if 0, no constant term will be included.
%                       Example - 'intercept',1 
%                       Data Types - double
%        init :      Search initialization. Scalar.
%                      It specifies the point where to initialize the search
%                       and start monitoring required diagnostics. if init is not
%                       specified it will be set equal to :
%                       p+1, if the sample size is smaller than 40;
%                       min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%                       Example - 'init',100 starts monitoring from step m=100 
%                       Data Types - double
%      nocheck:  Check input arguments. Scalar. 
%                       If nocheck is equal to 1 no check is performed on
%                       matrix y and matrix X. Notice that y and X are left
%                       unchanged. In other words the additional column of ones for
%                       the intercept is not added. As default nocheck=0. The
%                       controls on h, alpha and nsamp still remain
%                       Example - 'nocheck',1 
%                       Data Types - double
%        tstat:      the kind of t-statistics which have to be monitored.
%               Character.
%               tstat = 'trad' implies  monitoring of traditional t
%               statistics (out.Tols). In this case the estimate of \sigma^2 at step m
%               is based on s^2_m (notice that s^2_m<<\sigma^2 when m/n is
%               small) tstat = 'resc' (default) implies monitoring of
%               rescaled t statistics In this scale the estimate of
%               \sigma^2 at step m is based on s^_m / var_truncnorm(m/n)
%               where var_truncnorm(m/n) is the variance of the truncated
%               normal distribution.
%               Example - 'tstat','trad'
%               Data Types - char
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
%   The output consists of a structure 'out' containing the following fields:
%   RES:        n x (n-init+1) = matrix containing the monitoring of
%               scaled residuals
%               1st row = residual for first unit ......
%               nth row = residual for nth unit.
%   LEV:        (n+1) x (n-init+1) = matrix containing the monitoring of
%               leverage
%               1st row = leverage for first unit ......
%               nth row = leverage for nth unit.
%    BB:        n x (n-init+1) matrix containing the information about the units belonging
%               to the subset at each step of the forward search.
%               1st col = indexes of the units forming subset in the initial step
%               ...
%               last column = units forming subset in the final step (all units)
%   mdr:        n-init x 3 matrix which contains the monitoring of minimum
%               deletion residual or (m+1)ordered residual  at each step of
%               the forward search.
%               1st col = fwd search index (from init to n-1)
%               2nd col = minimum deletion residual
%               3rd col = (m+1)-ordered residual
%               Remark: these quantities are stored with sign, that is the
%               min deletion residual is stored with negative sign if
%               it corresponds to a negative residual
%   msr:        n-init+1 x 3 = matrix which contains the monitoring of
%               maximum studentized residual or m-th ordered residual
%               1st col = fwd search index (from init to n)
%               2nd col = maximum studentized residual
%               3rd col = (m)-ordered studentized residual
%   nor:        (n-init+1) x 4 matrix containing the monitoring of
%               normality test in each step of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = Asymmetry test
%               3rd col = Kurtosis test
%               4th col = Normality test
%  Bols:        (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward search
%    S2:        (n-init+1) x 3 matrix containing the monitoring of S2 or R2
%               in each step of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = monitoring of S2
%               3rd col = monitoring of R2
%   Coo:        (n-init+1) x 3 matrix containing the monitoring of Cook or
%               modified Cook distance in each step of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = monitoring of Cook distance
%               3rd col = monitoring of modified Cook distance
%  Tols:        (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated t-statistics (as specified in option input 'tstat'
%               in each step of the forward search
%    Un:        (n-init) x 11 Matrix which contains the unit(s)
%               included in the subset at each step of the fwd search
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,2) for example contains
%               the unit included in step init+1 Un(end,2) contains the
%               units included in the final step of the search
%     y:        A vector with n elements that contains the response
%               variable which has been used
%     X:        Data matrix of explanatory variables
%               which has been used (it also contains the column of ones if
%               input option intercept was missing or equal to 1)
%class :        string FSReda.
%
%
% See also LXS.m, FSRbsb.m
%
% References:
%
%   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
%   Springer Verlag, New York.
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('fsreda')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
       %Example of use of FSReda based on a starting point coming
       %from LMS
        n=200;
        p=3;
        randn('state', 123456);
        X=randn(n,p);
        % Uncontaminated data
        y=randn(n,1);
        % Contaminated data
        ycont=y;
        ycont(1:5)=ycont(1:5)+6;
        [out]=LXS(y,X,'nsamp',1000);
        out=FSReda(y,X,out.bs);
%}

%{
    %Example of use of function FSReda using a random start and traditional
    %t-stat monitoring
    out=FSReda(y,X,0,'tstat','trad');
%}

%{
    %Examples with real data: wool data
    xx=load('wool.txt');
    X=xx(:,1:3);
    y=log(xx(:,4));
    [out]=LXS(y,X,'nsamp',0);
    [out]=FSReda(y,X,out.bs,'tstat','scal');
%}


%% Input parameters checking
if nargin<1
    stack_loss=load('stack_loss.txt');
    y=stack_loss(:,4);
    X=stack_loss(:,1:3);
    % LMS using all subsets
    [out]=LXS(y,X,'nsamp',0);
    % Forward search with EDA purposes
    [out1]=FSReda(y,X,out.bs,'init',5);
    % Create scaled squared residuals
    out1.RES=out1.RES.^2;
    % Specify the characteristics of the highlighted trajectories 
    fground=struct;
    fground.LineStyle={'-';'--';':';'-.'};
    fground.LineWidth=3;
    fground.Color={'r'};
    % Plot of monitoring of scaled squared residuals 
    resfwdplot(out1,'fground',fground)
    title('Monitoring squared scaled residuals for Stack loss dataset')
    
    noargmsg = sprintf(['To use FSReda with your own data y (response vector) and X (design matrix), type:\n\n' ...
                     '    [out]=LXS(y,X) \n' ...
                     '    [out]=FSReda(y,X,out.bs)\n' ...
                     '    resfwdplot(out) \n\n' ...
                     'Type "help FSReda" for more information\n\n' ...
                     'Click OK to view the monitoring residual plot for the stack loss data.']);
                 titl='Example of the use of the Forward search with EDA purposes';
     [cdata, colormap] = imread('logo.png','BackgroundColor',(240/255)*[1 1 1]); 
     msgbox(noargmsg,titl,'custom',cdata,colormap,'modal');
    return
end

if nargin>3
    
    % chklist = vector containing the names of optional arguments
    chklist=varargin(1:2:length(varargin));
    
    % chkint gives the position of the option nocheck in vector chklist
    % It is empty if the option nocheck is not specified by the user
    chkint=find(strcmp('nocheck',chklist));
else
    chkint='';
end

% If nocheck=1 skip checks on y and X
if ~isempty(chkint) && cell2mat(varargin(2*chkint))==1;
    [n,p]=size(X);
else
    
    % The first argument which is passed is y
    if nargin<1 || isempty(y)==1
        error('FSDA:FSReda:missingInputs','vector y is missing');
    end
    
    [m,q]=size(y);
    
    % If y is a row vector it is transformed in a column vector
    if q~=1
        if m==1
            y=y';
        else
            error('FSDA:FSReda:Wrongy','y is not one-dimensional.');
        end
    end
    
    
    % The second argument which is passed is X
    if nargin<2 || isempty(X)==1
        error('FSDA:FSReda:missingInputs','Input matrix X not specified.');
    else
        % X must be a 2-dimensional array
        % (line below correspods to if ndims(X)>2)
        if ~ismatrix(X)
            error('FSDA:FSReda:WrongX','X must be a matrix (2D array)');
        end
    end
    
    % Missing values are removed from X and y
    na.X=~isfinite(X*ones(size(X,2),1)); % na.X is a logical vector
    na.y=~isfinite(y);
    if size(na.X,1)~=size(na.y,1)
        error('FSDA:FSReda:WrongInputs','Number of observations in X and y not equal.');
    end
    
    % Observations with missing or infinite values are ommitted.
    ok=~(na.X|na.y); % | = Element-wise logical OR
    X=X(ok,:);
    y=y(ok,:);
    
    % Now n is the new number of non missing observations
    n=length(y);
    
    if nargin<3
        error('FSDA:FSReda:missingInputs','Initial subset is missing');
    end
    
    % Now check if the intercept has to be added
    if nargin > 3
        
        % chklist = vector containing the names of optional arguments
        chklist=varargin(1:2:length(varargin));
        
        % chkint is non empty if the user has specified the option intercept
        % chkint gives the position of the option intercept in vector chklist
        % It is empty if the option intercept is not specified by the user
        chkint=find(strcmp('intercept',chklist));
        
        
        % if a value for the intercept has been specified
        % and this value is equal to 1
        % then add the colum of ones to matrix X
        if isempty(chkint) || cell2mat(varargin(2*chkint))==1
            X=cat(2,ones(n,1),X); % add column of ones
        end
    else
        % If the user has not specified a value for the intercept than the
        % column of ones is automatically attached
        X=cat(2,ones(n,1),X); % add column of ones
    end;
    
    % p is the number of parameters to be estimated
    p=size(X,2);
    
    if n < p
        error('FSDA:FSReda:NsmallerP',['Need more observations than variables: n=' int2str(size(X,1)) ' and p=' int2str(size(X,2)) ]);
    end
    
    rk=rank(X);
    if rk < p
        error('FSDA:FSReda:NoFullRank','Matrix X is singular');
    end
end

%% User options
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end
options=struct('intercept',1,'init',init,'tstat','scal','nocheck',0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSReda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

if nargin > 3
    
    % We now overwrite inside structure options the default values with
    % those chosen by the user
    % Notice that in order to do this we use dynamic field names
    for i=1:2:length(varargin);
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
        warning('FSDA:FSReda:NoFullRank','Unable to randomly sample full rank matrix');
    end
    yb=y(bsb);
else
    Xb=X(bsb,:);
    yb=y(bsb);
end

ini0=length(bsb);

% check init
init=options.init;
if  init <p;
    mess=sprintf(['Attention : init should be larger than p+1. \n',...
        'It is set to p.']);
    disp(mess);
    init=p;
elseif init<ini0;
    mess=sprintf(['Attention : init should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    disp(mess);
    init=ini0;
elseif init>=n;
    mess=sprintf(['Attention : init should be smaller than n. \n',...
        'It is set to n-1.']);
    disp(mess);
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
S2=[(init:n)' zer1];

% mdr= (n-init) x 3 matrix
% 1st column = fwd search index
% 2nd col min deletion residual among observerations non belonging to the
% subset
% 3rd column (m+1)-th ordered residual
% They are stored with sign, that is the min deletion residual
% is stored with negative sign if it corresponds to a negative residual
mdr=[(init:n-1)'  zer];

% mdr= (n-init+1) x 3 matrix which will contain max studentized residual
%  among bsb and m-th studentized residual
msr=[(init:n)'  zer1];

% coo= (n-init) x 3 matrix which will contain Cook distances
%  (2nd col) and modified Cook distance (3rd col)
coo=[((init+1):n)'  NaN(n-init,6)];

% nor= (n-init+1) x 3 matrix which will contain asymmetry (2nd col)
% kurtosis (3rd col) and normality test (4th col)
nor=[(init:n)'  zer1];

% Matrix RES will contain the resisuals for each unit in each step of the forward search
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


%% Start of the forward search
if (rank(Xb)~=p)
    warning('FSDA:FSReda:NoFullRank','The provided initial subset does not form full rank matrix');
    % FS loop will not be performed
else
    for mm=ini0:n;
        
        % if n>200 show every 100 steps the fwd search index
        if n>200;
            if length(intersect(mm,seq100))==1;
                disp(['m=' int2str(mm)]);
            end
        end
        
        NoRankProblem=(rank(Xb) == p);
        if NoRankProblem  % rank is ok
            b=Xb\yb;
            resBSB=yb-Xb*b;
            blast=b;   % Store correctly computed b for the case of rank problem
        else   % number of independent columns is smaller than number of parameters
            warning('FSDA:FSReda','Rank problem in step %d: Beta coefficients are used from the most recent correctly computed step',mm);
            b=blast;
        end
        
        if (mm>=init);
            
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
                % It is a proper leverage for the units belonging to susbet
                % It is a pseudo leverage for the unit not belonging to the subset
                mAm=Xb'*Xb;
                
                mmX=inv(mAm);
                dmmX=diag(mmX);
                % Notice that we could replace the lowwing line with 
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
            end;
            
        else
            Sb=0;
        end
        
        % e= vector of residual for all units using b estimated using subset
        e=y-X*b;
        
        if (mm>=init)
            % Store all residuals
            RES(:,mm-init+1)=e;
            
            if NoRankProblem
                % Store S2 for the units belonging to subset
                S2(mm-init+1,2)=Sb;
                
                % Store maximum studentized residual
                % among the units belonging to the subset
                msrsel=sort(abs(resBSB)./sqrt(Sb*(1-hi(bsb))));
                msr(mm-init+1,2)=msrsel(mm);
                
                % Store R2
                S2(mm-init+1,3)=1-var(resBSB)/var(yb);
            end
            
        end
        
        r(:,2)=e.^2;
        
        if mm>init;
            
            if NoRankProblem
                
                % Store in the second column of matrix coo the Cook
                % distance
                bib=Bols(mm-init+1,2:p+1)-Bols(mm-init,2:p+1);
                if S2(mm-init+1,2)>0;
                    coo(mm-init,2)=bib*mAm*(bib')/(p*S2(mm-init+1,2));
                end
                
                if length(unit)>5;
                    unit=unit(1:5);
                end
                if S2(mm-init,2)>0;
                    coo(mm-init,3:length(unit)+2)= 1./(1-hi(unit)).* sqrt(((mm-p)/p)*hi(unit).*r(unit,2)./S2(mm-init,2));
                end
            end % NoRankProblem
        end
        
        if mm<n;
            if mm>=init;
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
            
            % order the r_i and include the smallest among the units
            %  forming the group of potential outliers
            ord=sortrows(r,2);
            
            % bsb= units forming the new  subset
            bsb=ord(1:(mm+1),1);
            
            Xb=X(bsb,:);  % subset of X
            yb=y(bsb);    % subset of y
            
            if mm>=init;
                unit=setdiff(bsb,oldbsb);
                if length(unit)<=10
                    Un(mm-init+1,2:(length(unit)+1))=unit;
                else
                    disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                    Un(mm-init+1,2:end)=unit(1:10);
                end;
            end
            
            
            if mm < n-1;
                
                % ncl= units forming the new noclean
                ncl=ord(mm+2:n,1);
                
            end;
        end
        
        if mm >= init;
            if NoRankProblem
                if strcmp(options.tstat,'scal')
                    % Compute the variance of the truncated normal distribution
                    if mm<n
                        a=norminv(0.5*(1+mm/n));
                        corr=1-2*(n./mm).*a.*normpdf(a);
                    else
                        corr=1;
                    end
                    Tols(mm-init+1,2:end)=sqrt(corr)*Bols(mm-init+1,2:end)./sqrt(Sb*dmmX');
                elseif strcmp(options.tstat,'trad')
                    Tols(mm-init+1,2:end)=Bols(mm-init+1,2:end)./sqrt(Sb*dmmX');
                end
            end % NoRankProblem
        end
    end
end   %Rank check

%% Structure returned by function FSReda
out=struct;
out.RES=RES/sqrt(S2(end,2));
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
out.y=y;
out.X=X;
out.class='FSReda';
end

