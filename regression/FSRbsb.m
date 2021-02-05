function [Un,BB] = FSRbsb(y,X,bsb,varargin)
%FSRbsb returns the units belonging to the subset in each step of the forward search
%
%<a href="matlab: docsearchFS('FSRbsb')">Link to the help function</a>
%
% Required input arguments:
%
%    y    :     Response variable. Vector. A vector with n elements that contains
%               the response variable. y can be either a row or a column vector.
%  X :          Predictor variables. Matrix.
%               Matrix of explanatory variables (also called 'regressors')
%               of dimension n x (p-1) where p denotes the number of
%               explanatory variables including the intercept.
%               Rows of X represent observations, and columns represent
%               variables. By default, there is a constant term in the
%               model, unless you explicitly remove it using input option
%               intercept, so do not include a column of 1s in X. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%  bsb :        list of units forming the initial subset. Vector | 0. If
%               bsb=0 then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
%
% Optional input arguments:
%
%       init  :     Search initialization. Scalar.
%                   It specifies the initial subset size to start
%                   monitoring units forming subset
%                   Example - 'init',100 starts the search from step m=100
%                   Data Types - double
%
%    intercept :  Indicator for constant term. true (default) | false.
%                 Indicator for the constant term (intercept) in the fit,
%                 specified as the comma-separated pair consisting of
%                 'Intercept' and either true to include or false to remove
%                 the constant term from the model.
%                 Example - 'intercept',false
%                 Data Types - boolean
%
%    nocheck  :    Check input arguments. Boolean.
%                  If nocheck is equal to 1 no check is performed on
%                  matrix y and matrix X. Notice that y and X are left
%                  unchanged. In other words the additional column of ones for
%                  the intercept is not added. As default nocheck=false.
%                  Example - 'nocheck',true
%                  Data Types - boolean
%
%   bsbsteps :  Save the units forming subsets in selected steps. Vector.
%               It specifies for which steps of the fwd search it is
%               necessary to save the units forming subset. If bsbsteps is
%               0 we store the units forming subset in all steps. The
%               default is store the units forming subset in all steps if
%               n<=5000, else to store the units forming subset at steps
%               init and steps which are multiple of 100. For example, as
%               default, if n=7530 and init=6, units forming subset are
%               stored for
%               m=init, 100, 200, ..., 7500.
%               Example - 'bsbsteps',[100 200] stores the unis forming
%               subset in steps 100 and 200.
%               Data Types - double
%
%       plots   : Plot on the screen. Scalar.
%                 If plots=1 the monitoring units plot is displayed on the
%                 screen. The default value of plots is 0 (that is no plot
%                 is produced on the screen).
%                 Example - 'plots',1
%                 Data Types - double
%
%       msg    :  Level of output to display. Boolean. It controls whether
%                 to display or not messages on the screen
%                 If msg==true (default) messages are displayed on the screen about
%                   step of the fwd search
%                 else no message is displayed on the screen.
%               Example - 'msg',true
%               Data Types - boolean
%
% Output:
%
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
%  BB:          Units belonging to subset in each step or selected steps. Matrix.
%               n-by-(n-init+1) or n-by-length(bsbsteps) matrix which
%               contains the units belonging to the subset at each step (or
%               in selected steps as specified by optional vector bsbsteps)
%               of the forward search.
%               More precisely:
%               BB(:,1) contains the units forming subset in step bsbsteps(1);
%               ....;
%               BB(:,end) contains the units forming subset in step  bsbsteps(end);
%               Row 1 of matrix BB is referred to unit 1;
%               ......;
%               Row n of matrix BB is referred to unit n;
%               Units not belonging to subset are denoted with NaN.
%
% See also FSRBbsb, FSRHbsb
%
% See also: FSReda
%
% References:
%
% Chaloner, K. and Brant, R. (1988), A Bayesian Approach to Outlier
% Detection and Residual Analysis, "Biometrika", Vol. 75, pp. 651-659.
% Riani, M., Corbellini, A. and Atkinson, A.C. (2018), Very Robust Bayesian
% Regression for Fraud Detection, "International Statistical Review",
% http://dx.doi.org/10.1111/insr.12247
% Atkinson, A.C., Corbellini, A. and Riani, M., (2017), Robust Bayesian
% Regression with the Forward Search: Theory and Data Analysis, "Test",
% Vol. 26, pp. 869-886, DOI 10.1007/s11749-017-0542-6
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRbsb')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % FSRbsb with all default options.
    load('fishery');
    y=fishery{:,2};
    X=fishery{:,1};
    [out]=LXS(y,X,'nsamp',10000);
    Un = FSRbsb(y,X,out.bs);
%}

%{
    %% FSRbsb with optional arguments.
    % Monitoring units plot for fishery dataset
    load('fishery');
    y=fishery{:,2};
    X=fishery{:,1};
    [out]=LXS(y,X,'nsamp',10000);
    Un = FSRbsb(y,X,out.bs,'plots',1);
%}

%{
    %% Monitoring the units belonging to subset.
    state=1000;
    randn('state', state);
    n=100;
    X=randn(n,3);
    bet=[3;4;5];
    y=3*randn(n,1)+X*bet;
    y(1:20)=y(1:20)+15;
    [outLMS]=LXS(y,X);
    bsb=outLMS.bs;
    % Store in matrix BB the units belonging to subset in each step of the forward search
    [Un,BB] = FSRbsb(y,X,bsb);
    % Create the 'monitoring units plot'
    figure;
    seqr=[Un(1,1)-1; Un(:,1)];
    plot(seqr,BB','bx');
    xlabel('Subset size m');
    ylabel('Monitoring units plot');
    % The plot, which monitors the units belonging to subset in each step of
    % the forward search shows that apart from unit 11 which enters the
    % search in step m=78 all the other contaminated units enter the search
    % in the last 19 steps.

    % if we consider the seed state=500, we obtain a plot showing that the
    % 20 contaminated units enter the search in the final 20 steps.
    state=500;
    randn('state', state);
    X=randn(n,3);
    y=3*randn(n,1)+X*bet;
    y(1:20)=y(1:20)+15;
    [outLMS]=LXS(y,X);
    bsb=outLMS.bs;
    % Store in matrix BB the units belonging to subset in each step of the forward search
    [Un,BB] = FSRbsb(y,X,bsb);
    % Create the 'monitoring units plot'
    figure;
    seqr=[Un(1,1)-1; Un(:,1)];
    plot(seqr,BB','bx');
    xlabel('Subset size m');
    ylabel('Monitoring units plot');
%}

%{
    % Example with monitoring from step 60.
    load('fishery');
    y=fishery{:,2};
    X=fishery{:,1};
    [out]=LXS(y,X,'nsamp',10000);
    Un = FSRbsb(y,X,out.bs,'init',60);
%}

%{
    % FSR using a regression model without intercept.
    load('fishery');
    y=fishery{:,2};
    X=fishery{:,1};
    [out]=LXS(y,X);
    bsb=out.bs;
    [Un,BB] = FSRbsb(y,X,out.bs,'intercept','0');
%}

%{
    %FSR applied without doing any checks on y and X variables.
    load('fishery');
    y=fishery{:,2};
    X=fishery{:,1};
    [out]=LXS(y,X,'nsamp',10000);
    [Un,BB] = FSRbsb(y,X,out.bs,'nocheck','1');
%}

%% Beginning of code

% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options
if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end

if init<length(bsb)
    init=length(bsb);
end

bsbstepdef='';

if coder.target('MATLAB')
    
    options=struct('intercept',1,'init',init,'nocheck',false,'plots',0,...
        'bsbsteps',bsbstepdef,'msg',true);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSRbsb:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
end

if nargin<3
    error('FSDA:FSRbsb:missingInputs','Initial subset is missing');
end

if nargin > 3
    % We now overwrite inside structure options the default values with
    % those chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% Initialize Xb (necessary for MATLAB coder)
% otherwise the following message appears
% Variable 'Xb' is not fully defined on some execution paths.
Xb=X(1:p,:);

if bsb==0
    Ra=true; nwhile=1;
    while and(Ra,nwhile<100)
        bsb=randsample(n,p);
        Xb=X(bsb,:);
        Ra=(rank(Xb)<p);
        nwhile=nwhile+1;
    end
    if nwhile==100
        if coder.target('MATLAB')
            warning('FSDA:FSRbsb:NoFullRank','Unable to randomly sample full rank matrix');
        end
    end
    yb=y(bsb);
else
    Xb=X(bsb,:);
    yb=y(bsb);
end

ini0=length(bsb);

% check init
init=options.init;
if init <p
    if coder.target('MATLAB')
    fprintf(['Attention : init should be larger than p-1. \n',...
        'It is set to p.']);
    end
    init=p;
elseif init<ini0
    if coder.target('MATLAB')
    fprintf(['Attention : init should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    end
    init=ini0;
elseif init>=n
    if coder.target('MATLAB')
    fprintf(['Attention : init should be smaller than n. \n',...
        'It is set to n-1.']);
    end
    init=n-1;
end
init=init(1);
msg=options.msg;
nocheck=options.nocheck;

%% Initialise key matrices

zeron1=false(n,1);

% Initialization of the n x 1 Boolean vector which contains a true in
% correspondence of the units belonging to subset in each step
bsbT=zeron1;
bsbT(bsb)=true;

% sequence from 1 to n.
seq = (1:n)';

% The second column of matrix R will contain the OLS residuals at each step
% of the forward search
r = [seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100 = 100*(1:1:ceil(n/100));
seq100=seq100(seq100<=n);
seq100boo=false(n,1);
seq100boo(seq100)=true;


bsbsteps=options.bsbsteps;
% Matrix BB will contain the units forming subset in each step (or in
% selected steps) of the forward search. The first column contains
% information about units forming subset at step init1.
if isempty(bsbsteps)
    % Default for vector bsbsteps which indicates for which steps of the fwd
    % search units forming subset have to be saved
    if n<=5000
        bsbsteps = (init:1:n)';
    else
        bsbsteps = [init init+100-mod(init,100):100:100*floor(n/100)]';
    end
    BB = NaN(n,length(bsbsteps),'single');
elseif bsbsteps==0
    bsbsteps=(init:n)';
    BB = NaN(n,n-init+1,'single');
else
    if min(bsbsteps)<init && coder.target('MATLAB')
        warning('FSDA:FSMbsb:WrongInit','It is impossible to monitor the subset for values smaller than init');
    end
    boo=(bsbsteps>=init);
    bsbsteps=bsbsteps(boo);
    
    BB = NaN(n,length(bsbsteps),'single');
end

boobsbsteps=false(n,1);
boobsbsteps(bsbsteps)=true;


%  UN is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init+1:n)' , NaN(n-init,10));

% The last correctly computed beta oefficients
blast=NaN(p,1);

%% Forward search loop
if nocheck==false && rank(Xb)~=p
    if coder.target('MATLAB')
        warning('FSDA:FSRbsb:NoFullRank','The provided initial subset does not form a full rank matrix');
    end
    % FS loop will not be performed
    Un=NaN;
else
    % ij = index which is linked with the columns of matrix BB. During the
    % search every time a subset is stored inside matrix BB ij icreases by one
    ij=1;
    
    for mm = ini0:n
        % if n>200 show every 100 steps the fwd search index
        if n>200 && msg== true
            
            if  seq100boo(mm) == true
                % OLD CODE if length(intersect(mm,seq100))==1
                disp(['m=' int2str(mm)]);
            end
        end
        
        % Store units belonging to the subset
        if (mm>=init)
            if boobsbsteps(mm)==true
                % OLD CODE
                % intersect(mm,bsbsteps)==mm
                BB(bsb,ij)=bsb;
                ij=ij+1;
            end
        end
        
        % Compute beta coefficients using subset
        if nocheck==1
            NoRankProblem=true;
        else
            rankXb=rank(Xb);
            NoRankProblem=(rankXb == p);
        end
        
        if NoRankProblem  % rank is ok
            b = Xb\yb;
            blast=b;
        else
            b=blast;    % in case of rank problem, the last orrectly computed coefficients are used
            if coder.target('MATLAB')
                warning('FSDA:FSRbsb:NoFullRank','Rank problem in step %d: Beta coefficients are used from the most recent correctly computed step',mm);
            end
        end
        
        % e= vector of residual for all units using b estimated using subset
        e=y-X*b;
        
        r(:,2)=e.^2;
        
        if mm<n
            
            % store units forming old subset in vector oldbsb
            oldbsbT=bsbT;
            
            % order the r_i and include the smallest among the units forming
            % the group of potential outliers
            % ord=sortrows(r,2);
            [~,ord]=sort(r(:,2));
            
            % bsb= units forming the new subset
            bsb=ord(1:(mm+1),1);
            
            bsbT=zeron1;
            bsbT(bsb)=true;
            
            
            Xb=X(bsb,:);  % subset of X
            yb=y(bsb);    % subset of y
            
            if mm>=init
                
                % unit = vector containing units which just entered subset;
                % unit=setdiff(bsb,oldbsb);
                % new instruction to find unit
                unit=seq(bsbT & ~oldbsbT);
                
                % If the interchange involves more than 10 units, store only the
                % first 10.
                if (size(unit,2)<=10)
                    Un(mm-init+1,2:(size(unit,1)+1)) = unit;
                else
                    disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                    Un(mm-init+1,2:end) = unit(1:10)';
                end
            end
        end
    end
    
    plots=options.plots;
    if coder.target('MATLAB') &&  plots==1
        % Create the 'monitoring units plot'
        figure;
        plot(bsbsteps,BB','bx')
        xlabel('Subset size m');
        ylabel('Monitoring units plot');
    end
end
%FScategory:REG-Regression