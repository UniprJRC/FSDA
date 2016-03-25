function [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,bsb,varargin)
%FSRmdr computes minimum deletion residual and other basic linear regression quantities in each step of the search
%
%<a href="matlab: docsearchFS('fsrmdr')">Link to the help function</a>
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
%  X :          Predictor variables. Matrix. Matrix of explanatory
%               variables (also called 'regressors') of dimension n x (p-1)
%               where p denotes the number of explanatory variables
%               including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless you explicitly remove it using input option
%               intercept, so do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%  bsb     :    list of units forming the initial subset. Vector. If bsb=0
%               (default) then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
%
% Optional input arguments:
%
%  init :       Search initialization. Scalar.
%               It specifies the point where to initialize the search and
%               start monitoring required diagnostics. If it is not
%               specified it is set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
%  intercept :  Indicator for constant term. Scalar. If 1, a model with
%               constant term will be fitted (default), if 0, no constant
%               term will be included.
%               Example - 'intercept',1
%               Data Types - double
%  plots :      Plot on the screen. Scalar. If equal to one a plot of
%               minimum deletion residual appears  on the screen with 1%,
%               50% and 99% confidence bands else (default) no plot is
%               shown.
%               Example - 'plots',1
%               Data Types - double
%               Remark: the plot which is produced is very simple. In order
%               to control a series of options in this plot and in order to
%               connect it dynamically to the other forward plots it is
%               necessary to use function mdrplot.
%  nocheck:     Check input arguments. Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additioanl column of ones for
%               the intercept is not added. As default nocheck=0. The
%               controls on h, alpha and nsamp still remain
%               Example - 'nocheck',1
%               Data Types - double
%  msg  :       Level of output to display. Scalar. It controls whether to
%               display or not messages about great interchange on the
%               screen If msg==1 (default)
%               messages are displayed on the screen
%               else no message is displayed on the screen
%               Example - 'msg',1
%               Data Types - double
%  constr :     Constrained search. Vector. r x 1 vector which contains the list of units which are
%               forced to join the search in the last r steps. The default
%               is constr=''.  No constraint is imposed
%               Example - 'constr',[1:10] forces the first 10 units to join
%               the subset in the last 10 steps
%               Data Types - double
% bsbmfullrank :What to do in case subset at step m (say bsbm) produces a
%               non singular X. Scalar.
%               This options controls what to do when rank(X(bsbm,:)) is
%               smaller then number of explanatory variables.
%               If bsbmfullrank = 1 (default is 1) these units (whose number
%               is say mnofullrank) are constrained to enter the search in
%               the final n-mnofullrank steps else the search continues
%               using as estimate of beta at step m the estimate of beta
%               found in the previous step.
%               Example - 'bsbmfullrank',1
%               Data Types - double
%   bsbsteps :  Save the units forming subsets. Vector. It specifies for
%               which steps of the fwd search it
%               is necessary to save the units forming subsets. If bsbsteps
%               is 0 we store the units forming subset in all steps. The
%               default is store the units forming subset in all steps if
%               n<=5000, else to store the units forming subset at steps
%               init and steps which are multiple of 100. For example, as
%               default, if n=753 and init=6,
%               units forming subset are stored for
%               m=init, 100, 200, 300, 400, 500 and 600.
%               Example - 'bsbsteps',[100 200] stores the unis forming
%               subset in steps 100 and 200.
%               Data Types - double
%  Remark:      The user should only give the input arguments that have to
%               change their default value.
%               The name of the input arguments needs to be followed by
%               their value. The order of the input arguments is of no
%               importance.
%
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations. y can be both a row of column vector.
%
% Output:
%
%  mdr:          n -init x 2 matrix which contains the monitoring of minimum
%               deletion residual at each step of the forward search.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum deletion residual.
%               REMARK: if in a certain step of the search matrix is
%               singular, this procedure checks ohw many observations
%               produce a singular matrix. In this case mdr is a column
%               vector which contains the list of units for which matrix X
%               is non singular.
%  Un:          Units included in each step. Matrix.
%               (n-init) x 11 Matrix which contains the unit(s) included
%               in the subset at each step of the search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
%  BB:          Units belonging to subset in each step. Matrix.
%               n x (n-init+1) or n-by-length(bsbsteps) matrix (depending on input
%               option bsbsteps) which contains information about the units
%               belonging to the subset at each step of the forward search.
%               BB has the following structure
%               1-st row has number 1 in correspondence of the steps in
%                   which unit 1 is included inside subset and a missing
%                   value for the other steps
%               ......
%               (n-1)-th row has number n-1 in correspondence of the steps in
%                   which unit n-1 is included inside subset and a missing
%                   value for the other steps
%               n-th row has number n in correspondence of the steps in
%                   which unit n is included inside subset and a missing
%                   value for the other steps
%               The size of matrix BB is n x (n-init+1) if option input
%               bsbsteps is 0 else the size is n-by-length(bsbsteps).
%  Bols:        OLS coefficents. Matrix.
%               (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward
%               search.
%  S2:          S2 and R2. Matrix.
%               (n-init+1) x 3 matrix containing the monitoring of S2 (2nd
%               column)and R2 (third column) in each step of the forward
%               search.
%
% See also: FSR, FSReda
%
% References:
%
%   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
%   Springer Verlag, New York.
%   Atkinson, A.C. and Riani, M. (2006). Distribution theory and
%   simulations for tests of outliers in regression. Journal of
%   Computational and Graphical Statistics, Vol. 15, pp. 460-476
%   Riani, M. and Atkinson, A.C. (2007). Fast calibrations of the forward
%   search for testing multiple outliers in regression, Advances in Data
%   Analysis and Classification, Vol. 1, pp. 123-141.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('fsrmdr')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % FSRmdr with all default options.
    % Compute minimum deletion residual.
    % Monitor minimum deletion residual in each step of the forward search.
    % Common part to all examples: load fishery dataset.
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
     [mdr] = FSRmdr(y,X,out.bs);
     plot(mdr(:,1),mdr(:,2))
     title('Monitoring of minimum deletion residual')
%}

%{
    % FSRmdr with optional arguments.
    % Choose step to start monitoring.
    % Compute minimum deletion residual and start monitoring it from step
    % 60.
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr] = FSRmdr(y,X,out.bs,'init',60);
%}

%{
    % Analyze units entering the search in the final steps.
    % Compute minimum deletion residual and analyze the units entering
    % subset in each step of the fwd search (matrix Un).  As is well known,
    % the FS provides an ordering of the data from those most in agreement
    % with a suggested model (which enter the first steps) to those least in
    % agreement with it (which are included in the final steps).
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un] = FSRmdr(y,X,out.bs);
%}


%{
    % Units forming subset in each step.
    % Obtain detailed information about the units forming subset in each
    % step of the forward search (matrix BB).
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un,BB] = FSRmdr(y,X,out.bs);
%}

%{
    % Monitor \( \hat  \beta \).
    % Monitor how the estimates of beta coefficients changes as the subset
    % size increases (matrix Bols).
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un,BB,Bols] = FSRmdr(y,X,out.bs);
%}

%{
    % Monitor $s^2$.
    % Monitor the estimate of $\sigma^2$ in each step of the fwd search
    % (matrix S2).
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs);
    plot(S2(:,1),S2(:,2))
    title('Monitoring of s2')
%}

%{
    % Specify a regression model without intercept.
    % FSRmdr using a regression model without intercept.
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs,'intercept','0');
%}

%{
    % Example of the use of option nocheck.
    %FSRmdr applied without doing any checks on y and X variables.
     load('fishery');
     y=fishery.data(:,2);
     X=fishery.data(:,1);
     % Find starting subset
     [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs,'nocheck',1);
%}

%{
%% Monitoring of s2, and the evolution of beta coefficients for the Hawkins dataset
    load('hawkins.txt');
    y=hawkins(:,9);
    X=hawkins(:,1:8);
    [out]=LXS(y,X,'nsamp',10000);
    [~,~,~,Bols,S2] = FSRmdr(y,X,out.bs);
    %The forward plot of s2 shows that initially the estimate is virtually
    %zero. The four line segments comprising the plot correspond to the four
    %groups of observations

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

    %The forward plots of the estimates of the beta coefficients show that they are virtually constant until m = 86, after which they start fluctuating in different directions.

    % Plots of the monitoring of "Estimates of beta coefficients"
    figure;
    for j=3:size(Bols,2);
        subplot(2,4,j-2)
        plot(Bols(:,1),Bols(:,j))
        xlabel('Subset size m');
        ylabel(['b' num2str(j-2)]);
    end
%}

%{
    %% Store units forming subsets in selected steps.
    % In this example the units forming subset are stored just for
    % selected steps.
    load('hawkins.txt');
    y=hawkins(:,9);
    X=hawkins(:,1:8);
    [out]=LXS(y,X,'nsamp',10000);
    [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs,'bsbsteps',[30 60]);
    % BB has just two columns
    % First column contains information about units forming subset at step m=30
    % sum(~isnan(BB(:,1))) is 30
    % Second column contains information about units forming subset at step m=60
    % sum(~isnan(BB(:,2))) is 60
    disp(sum(~isnan(BB(:,1))))
    disp(sum(~isnan(BB(:,2))))
%}


%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

if n<40
    initdef=p+1;
else
    initdef=min(3*p+1,floor(0.5*(n+p+1)));
end

% Default for vector bsbsteps which indicates for which steps of the fwd
% search units forming subset have to be saved
if n<=5000
    bsbstepdef = 0;
else
    bsbstepdef = [initdef 100:100:100*floor(n/100)];
end

options=struct('intercept',1,'init',initdef,'plots',0,'nocheck',0,'msg',1,...
    'constr','','bsbmfullrank',1,'bsbsteps',bsbstepdef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRmdr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<3
    error('FSDA:FSRmdr:missingInputs','Initial subset is missing');
end
%init1=options.init;
if nargin > 3
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
% And check if the optional user parameters are reasonable.

if bsb==0;
    Ra=1; nwhile=1;
    while and(Ra,nwhile<100)
        bsb=randsample(n,p);
        Xb=X(bsb,:);
        Ra=(rank(Xb)<p);
        nwhile=nwhile+1;
    end
    if nwhile==100
        warning('FSDA:FSRmdr:NoFullRank','Unable to randomly sample full rank matrix');
    end
    yb=y(bsb);
else
    Xb=X(bsb,:);
    yb=y(bsb);
end

ini0=length(bsb);

% check init
init1=options.init;
if  init1 <p+1;
    mess=sprintf(['Attention : init1 should be larger than p. \n',...
        'It is set to p+1.']);
    disp(mess);
    init1=p+1;
elseif init1<ini0;
    mess=sprintf(['Attention : init1 should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    disp(mess);
    init1=ini0;
elseif init1>=n;
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    disp(mess);
    init1=n-1;
end


msg=options.msg;
constr=options.constr;
bsbmfullrank=options.bsbmfullrank;
bsbsteps=options.bsbsteps;


%% Initialise key matrices

% sequence from 1 to n.
seq=(1:n)';

% the set complementary to bsb.
ncl=setdiff(seq,bsb);

% The second column of matrix R will contain the OLS residuals at each step
% of the search
r=[seq zeros(n,1)];

% If n is very large (>500), the step of the search is printed every 500 step
% seq500 is linked to printing
seq500=500*(1:1:ceil(n/500));

% Matrix Bols will contain the beta coefficients in each step of the fwd
% search. The first column of Bols contains the fwd search index
Bols=[(init1:n)' NaN(n-init1+1,p)];     %initial value of beta coefficients is set to NaN

% S2 = (n-init1+1) x 3 matrix which will contain:
% 1st col = fwd search index
% 2nd col = S2= \sum e_i^2 / (m-p)
% 3rd col = R^2
S2=[(init1:n)' NaN(n-init1+1,2)];        %initial value of S2 (R2) is set to NaN

% mdr = (n-init1-1) x 2 matrix which will contain min deletion residual
% among nobsb r_i^*
mdr=[(init1:n-1)'  NaN(n-init1,1)];      %initial value of mdr is set to NaN

% Matrix BB will contain the units forming subset in each step (or in
% selected steps) of the forward search. The first column contains
% information about units forming subset at step init1.
if bsbsteps == 0
    bsbsteps=init1:n;
    BB = NaN(n,n-init1+1);
else
    % The number of columns of matrix BB is equal to the number of steps
    % for which bsbsteps is greater or equal than init1
    bsbsteps=bsbsteps(bsbsteps>=init1);
    BB = NaN(n,length(bsbsteps));
    %   BB = NaN(n, sum(bsbsteps>=init1));
end

% ij = index which is linked with the columns of matrix BB. During the
% search every time a subset is stored inside matrix BB ij increases by one
ij=1;


%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init1+1:n)' , NaN(n-init1,10));


%% Start of the forward search
if (rank(Xb)~=p)
    warning('FSDA:FSRmdr:NoFullRank','Supplied initial subset does not produce full rank matrix');
    warning('FSDA:FSRmdr:NoFullRank','FS loop will not be performed');
    mdr=NaN;
    % FS loop will not be performed
else
    for mm=ini0:n;
        
        
        
        % if n>5000 show every 500 steps the fwd search index
        if msg==1 && n>5000;
            if length(intersect(mm,seq500))==1;
                disp(['m=' int2str(mm)]);
            end
        end
        
        
        NoRankProblem=(rank(Xb) == p);
        if NoRankProblem  % rank is ok
            b=Xb\yb;
            resBSB=yb-Xb*b;
            blast=b;   % Store correctly computed b for the case of rank problem
        else   % number of independent columns is smaller than number of parameters
            if bsbmfullrank
                Xbx=Xb;
                nclx=ncl;
                bsbx=zeros(n,1);
                bsbx(1:mm)=bsb;
                norank=1;
                while norank ==1
                    
                    norank=1;
                    % Increase the size of the subset by one unit iteratively until you
                    % obtain a full rank matrix
                    for i=1:length(nclx)
                        Xbb=[Xbx;X(nclx(i),:)];
                        if rank(Xbb)==p
                            norank=0;
                        else
                            bsbx(1:size(Xbb,1))=[bsbx(1:size(Xbb,1)-1);nclx(i)];
                            Xbx=X(bsbx(1:size(Xbb,1)),:);
                            nclx=setdiff(seq,bsbx(1:size(Xbb,1)));
                            norank=1;
                            break
                        end
                    end
                end
                % check how many observations produce a singular X matrix
                bsbsing=bsbx(1:size(Xbb,1)-1);
                
                if msg==1
                    warning('FSDA:FSRmdr','Rank problem in step %d:',mm);
                    disp('Observations')
                    disp(bsbsing')
                    disp('produce a singular matrix')
                end
                mdr=bsbsing;
                Un=NaN;
                BB=NaN;
                Bols=NaN;
                S2=NaN;
                return
                
            else
                disp(['Matrix without full rank at step m=' num2str(mm)])
                disp('Estimate of \beta which is used is based on previous step with full rank')
                b=blast;
                % disp([mm b'])
            end
        end
        e=y-X*b;  % e = vector of residual for all units using b estimated using subset
        r(:,2)=e.^2;
        
        if (mm>=init1);
            % Store units belonging to the subset
            if intersect(mm,bsbsteps)==mm
                BB(bsb,ij)=bsb;
                ij=ij+1;
            end
            
            if NoRankProblem
                % Store beta coefficients if there is no rank problem
                Bols(mm-init1+1,2:p+1)=b';
                % Compute and store estimate of sigma^2
                Sb=(resBSB)'*(resBSB)/(mm-p);
                S2(mm-init1+1,2)=Sb;
                % Store R2
                S2(mm-init1+1,3)=1-var(resBSB)/var(yb);
                
                if mm<n
                    mAm=Xb'*Xb;
                    
                    % Take minimum deletion residual among noBSB
                    % hi (the leverage for the units not belonging to the
                    % subset) is defined as follows
                    % hi=diag(X(ncl,:)*inv(Xb'*Xb)*(X(ncl,:))');
                    
                    % Take units not belonging to bsb
                    Xncl = X(ncl,:);
                    
                    % mmX=inv(mAm);
                    % hi = sum((Xncl*mmX).*Xncl,2);
                    hi=sum((Xncl/mAm).*Xncl,2);
                    
                    ord = [(r(ncl,2)./(1+hi)) e(ncl)];
                    
                    % Store minimum deletion residual in matrix mdr
                    selmdr=sortrows(ord,1);
                    if S2(mm-init1+1,2)==0
                        warning('FSDA:FSRmdr:ZeroS2','Value of S2 at step %d is zero, mdr is NaN',mm-init1+1);
                    else
                        mdr(mm-init1+1,2)=sqrt(selmdr(1,1)/S2(mm-init1+1,2));
                    end
                end  %if mm<n
            end   %~RankProblem
        end     %mm>=init1
        
        if mm<n;
            
            % store units forming old subset in vector oldbsb
            oldbsb=bsb;
            
            % order the r_i
            
            % units inside vector constr are forced to join the search in
            % the final k steps
            if ~isempty(constr) && mm<n-length(constr)
                r(constr,2)=Inf;
            end
            ord=sortrows(r,2);
            
            % bsb= units forming the new  subset
            bsb=ord(1:(mm+1),1);
            
            Xb=X(bsb,:);  % subset of X
            yb=y(bsb);    % subset of y
            
            if mm>=init1;
                unit=setdiff(bsb,oldbsb);
                
                % If the interchange involves more than 10 units, store only the
                % first 10.
                if length(unit)<=10
                    Un(mm-init1+1,2:(length(unit)+1))=unit;
                else
                    if msg==1
                        disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        disp(['Number of units which entered=' int2str(length(unit))]);
                        Un(mm-init1+1,2:end)=unit(1:10);
                    end
                end
            end
            
            if mm < n-1;
                if ~isempty(constr) && mm<n-length(constr)-1
                    % disp(mm)
                    ncl=ord(mm+2:n,1);    % ncl= units forming the new noclean
                    ncl=setdiff(ncl,constr);
                else
                    ncl=ord(mm+2:n,1);    % ncl= units forming the new noclean
                end
                
            end
        end   % if mm<n
    end  % for mm=ini0:n loop
    
    % Plot minimum deletion residual with 1%, 50% and 99% envelopes
    if options.plots==1
        quant=[0.01;0.5;0.99];
        % Compute theoretical envelops for minimum deletion residual based on all
        % the observations for the above quantiles.
        [gmin] = FSRenvmdr(n,p,'prob',quant,'init',init1,'exact',1);
        plot(mdr(:,1),mdr(:,2));
        
        % Superimpose 1%, 99%, 99.9% envelopes based on all the observations
        lwdenv=2;
        % Superimpose 50% envelope
        line(gmin(:,1),gmin(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
        % Superimpose 1% and 99% envelope
        line(gmin(:,1),gmin(:,2),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        line(gmin(:,1),gmin(:,4),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        
        xlabel('Subset size m');
        ylabel('Monitoring of minimum deletion residual');
    end
    
end % rank check

end