function [mdr,Un,BB,Bgls,S2,Hetero,WEI] = FSRHmdr(y,X,Z,bsb,varargin)
%FSRHmdr computes minimum deletion residual and other basic linear regression quantities in each step of the heteroskedastic search
%
%<a href="matlab: docsearch('fsrhmdr')">Link to the help function</a>
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
%  X :          Predictor variables in the regression equation. Matrix.
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
%     Z :       Predictor variables in the scedastic equation. Matrix.
%               n x r matrix or vector of length r.
%               If Z is a n x r matrix it contains the r variables which
%               form the scedastic function as follows (if input option art==1)
%               \[
%               \omega_i = 1 + exp(\gamma_0 + \gamma_1 Z(i,1) + ...+ \gamma_{r} Z(i,r))
%               \]
%               If Z is a vector of length r it contains the indexes of the
%               columns of matrix X which form the scedastic function as
%               follows
%               \[
%               \omega_i = 1 +  exp(\gamma_0 + \gamma_1 X(i,Z(1)) + ...+
%               \gamma_{r} X(i,Z(r)))
%               \]
%
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedasticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax
%                    FSRHmdr(y,X,X(:,[3 5]),0)
%               or the sintax
%                    FSRHmdr(y,X,[3 5],0)
%  bsb :        list of units forming the initial subset. Vector | 0. If
%               bsb=0 then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
%
% Optional input arguments:
%
%  init :       Search initialization. Scalar.
%               It specifies the point where to start monitoring
%               required diagnostics. If it is not specified it is set
%               equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%               The minimum value of init is 0. In this case in the first
%               step we just use prior information
%               Example - 'init',100 starts monitoring from step m=100
%               Data Types - double
%  intercept :   Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double
% modeltype:    Parametric function to be used in the skedastic equation.
%               String.
%               If modeltype is 'arc' (default) than the skedastic function is
%               modelled as follows
%               \[
%               \sigma^2_i = \sigma^2 (1 + \exp(\gamma_0 + \gamma_1 Z(i,1) +
%                           \cdots + \gamma_{r} Z(i,r)))
%               \]
%               on the other hand, if modeltype is 'har' then traditional
%               formulation due to Harvey is used as follows
%               \[
%               \sigma^2_i = \exp(\gamma_0 + \gamma_1 Z(i,1) + \cdots +
%                           \gamma_{r} Z(i,r)) =\sigma^2 (\exp(\gamma_1
%                           Z(i,1) + \cdots + \gamma_{r} Z(i,r))
%               \]
%               Example - 'modeltype','har'
%               Data Types - string
%  plots :    Plot on the screen. Scalar.
%               If equal to one a plot of Bayesian minimum deletion residual
%               appears  on the screen with 1 per cent, 50 per cent and 99 per cent confidence
%               bands else (default) no plot is shown.
%               Remark. the plot which is produced is very simple. In order
%               to control a series of options in this plot and in order to
%               connect it dynamically to the other forward plots it is necessary to use
%               function mdrplot
%                 Example - 'plots',1
%                 Data Types - double
%  nocheck:   Check input arguments. Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additional column of ones for
%               the intercept is not added. As default nocheck=0.
%               Example - 'nocheck',1
%               Data Types - double
%  msg  :    Level of output to display. Scalar.
%               It controls whether to display or not messages
%               about great interchange on the screen
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen
%               Example - 'msg',1
%               Data Types - double
% gridsearch:   Algorithm to be used. Scalar.
%               If gridsearch ==1 grid search will be used else the
%               scoring algorith will be used.
%               Example - 'gridsearch',0
%               Data Types - double
%               REMARK: the grid search has only been implemented when
%               there is just one explantory variable which controls
%               heteroskedasticity
%  constr :    units which are forced to join the search in the last r steps. Vector.
%               r x 1 vector. The default is constr=''.  No constraint is imposed
%               Example - 'constr',[1 6 3]
%               Data Types - double
% bsbmfullrank :It tells how to behave in case subset at step m
%               (say bsbm) produces a non singular X. Scalar.
%               In other words, this options controls what to do when rank(X(bsbm,:)) is
%               smaller then number of explanatory variables. If
%               bsbmfullrank = 1 (default is 1) these units (whose number is
%               say mnofullrank) are constrained to enter the search in
%               the final n-mnofullrank steps else the search continues
%               using as estimate of beta at step m the estimate of beta
%               found in the previous step.
%               Example - 'bsbmfullrank',0
%               Data Types - double
%   bsbsteps :  Save the units forming subsets (and weights vector) in
%               selected steps. Vector.
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
%  mdr:         Minimum deletion residual. Matrix.
%               n -init x 2 matrix which contains the monitoring of minimum
%               deletion residual at each step of the forward search.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum deletion residual.
%               REMARK: if in a certain step of the search matrix is
%               singular, this procedure checks how many observations
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
%  Bgls:        GLS beta coefficents. Matrix.
%               (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward
%               search.
%  S2:          S2 and R2. Matrix.
%               (n-init+1) x 3 matrix containing the monitoring of S2 (2nd
%               column)and R2 (third column) in each step of the forward
%               search.
%  Hetero :     Coefficients in the heteroskedastic equation. Matrix.
%               (n-init1+1) x (r+1) matrix containing:
%                  1st col = fwd search index;
%                  2nd col = estimate of first coeff in the scedastic
%                  equation;
%                  ...
%                  (r+1)-th col = estimate of last coeff in the scedastic equation.
%   WEI:     Weights. Matrix.
%               n x (n-init+1) or n-by-length(bsbsteps) matrix (depending on input
%               option bsbsteps) which contains information about the
%               weights assigned to each unit to make the regression equation
%               skedastic.
%            More precisely, if $var (\epsilon)= \sigma^2
%            Omega=diag(omegahat)$ the weights which are stored are
%            $omegahat.^(-0.5)$;
%
% See also:   FSRmdr
%
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
%<a href="matlab: docsearch('FSRHmdr')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    %% FSRHmdr with all default options.
    % Common part to all examples: load tradeH dataset (used in the paper ART).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    mdr=FSRHmdr(y,X,Z,0);
%}

%{
    % FSRHmdr with optional arguments.
    % Specifying the search initialization.
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    mdr=FSRHmdr(y,X,Z,0,'init',round(length(y)/2));
%}

%{
    % Analyze units entering the search in the final steps.
    % Compute minimum deletion residual and analyze the units entering
    % subset in each step of the fwd search (matrix Un).  As is well known,
    % the FS provides an ordering of the data from those most in agreement
    % with a suggested model (which enter the first steps) to those least in
    % agreement with it (which are included in the final steps).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    [mdr,Un]=FSRHmdr(y,X,Z,0,'init',round(length(y)/2));
%}

%{
    % Units forming subset in each step.
    % Obtain detailed information about the units forming subset in each
    % step of the forward search (matrix BB).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    [mdr,Un,BB]=FSRHmdr(y,X,Z,0,'init',round(length(y)/2));
%}

%{
    % Monitor $\hat  \beta$.
    % Monitor how the estimates of beta coefficients changes as the subset
    % size increases (matrix Bols).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    [mdr,Un,BB,Bols]=FSRHmdr(y,X,Z,0,'init',round(length(y)/2));
%}

%{
    % Monitor $s^2$.
    % Monitor the estimate of sigma^2 in each step of the fwd search
    % (matrix S2).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    [mdr,Un,BB,Bols,S2]=FSRHmdr(y,X,Z,0,'init',round(length(y)/2));
%}

%{
    %% Monitoring the estimates of the scedastic equation.
    % With plot of the \alpha parameter.
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    [mdr,Un,BB,Bols,S2,Hetero]=FSRHmdr(y,X,Z,0,'init',round(length(y)/2));
    plot(Hetero(:,1),Hetero(:,2))
%}

%{
    % Monitoring the estimates of the weights.
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    [mdr,Un,BB,Bols,S2,Hetero,WEI]=FSRHmdr(y,X,Z,0,'init',round(length(y)/2));
    plot(S2(:,1),WEI')
    title('Monitoring of the weights')
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

bsbstepdef='';

options=struct('intercept',1,'init',initdef,'plots',0,'nocheck',0,'msg',1,...
    'constr','','bsbmfullrank',1,'modeltype','art','gridsearch',0,'bsbsteps',bsbstepdef);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSRHmdr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<4
    error('FSDA:FSRHmdr:missingInputs','Initial subset is missing');
end

if nargin > 4
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end


% Z = n-by-r matrix which contains the explanatory variables for
% heteroskedasticity
if size(Z,1)~=n
    % Check if intercept was true
    intercept=options.intercept;
    if intercept==1
        Z=X(:,Z+1);
    else
        Z=X(:,Z);
    end
end


if bsb==0;
    Ra=1; nwhile=1;
    while and(Ra,nwhile<100)
        bsb=randsample(n,p);
        Xb=X(bsb,:);
        Ra=(rank(Xb)<p);
        nwhile=nwhile+1;
    end
    if nwhile==100
        warning('FSDA:FSRHmdr:NoFullRankMatrix','Unable to randomly sample full rank matrix');
    end
    Zb=Z(bsb,:);
    yb=y(bsb);
    
else
    Xb=X(bsb,:);
    Zb=Z(bsb,:);
    yb=y(bsb);
end

gridsearch=options.gridsearch;
if gridsearch==1 && size(Z,2)>1
    warning('FSDA:FSRHmdr:WrongInputOpts','To perform a grid search you cannot have more than one varaible responsible for heteroskedasticity');
    warning('FSDA:FSRHmdr:WrongInputOpts','Scoring algorith is used');
    gridsearch=0;
end

modeltype=options.modeltype;

if strcmp(modeltype,'art') ==1
    art=1;
else
    art=0;
end

ini0=length(bsb);

% check init
init=options.init;
if  init <p+1;
    mess=sprintf(['Attention : init should be larger than p. \n',...
        'It is set to p+1.']);
    disp(mess);
    init=p+1;
elseif init<ini0;
    mess=sprintf(['Attention : init should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    disp(mess);
    init=ini0;
elseif init>=n;
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    disp(mess);
    init=n-1;
end

msg=options.msg;
constr=options.constr;
intercept=options.intercept;

%% Initialise key matrices


% sequence from 1 to n.
seq=(1:n)';

% the set complementary to bsb.
ncl=setdiff(seq,bsb);

% The second column of matrix R will contain the OLS residuals at each step
% of the search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100=100*(1:1:ceil(n/100));

% Matrix Bgls will contain the GLS beta coefficients in each step of the fwd
% search. The first column of Bgls contains the fwd search index
Bgls=[(init:n)' NaN(n-init+1,p)];     %initial value of beta coefficients is set to NaN

% S2 = (n-init1+1) x 3 matrix which will contain:
% 1st col = fwd search index
% 2nd col = S2= \sum e_i^2 / (m-p)
% 3rd col = R^2
S2=[(init:n)' NaN(n-init+1,2)];        %initial value of S2 (R2) is set to NaN

% mdr = (n-init1-1) x 2 matrix which will contain min deletion residual
% among nobsb r_i^*
mdr=[(init:n-1)'  NaN(n-init,1)];      %initial value of mdr is set to NaN


bsbsteps=options.bsbsteps;
% Matrix BB will contain the units forming subset in each step (or in
% selected steps) of the forward search. The first column contains
% information about units forming subset at step init1.
% WEI= matrix which contains in each column the estimate of the weights
if isempty(bsbsteps)
    % Default for vector bsbsteps which indicates for which steps of the fwd
    % search units forming subset have to be saved
    if n<=5000
        bsbsteps = init:1:n;
    else
        bsbsteps = [init init+100-mod(init,100):100:100*floor(n/100)];
    end
    BB = NaN(n,length(bsbsteps),'single');
    WEI = NaN(n,length(bsbsteps));
    
elseif bsbsteps==0
    bsbsteps=init:n;
    BB = NaN(n,n-init+1,'single');
    WEI = NaN(n,n-init+1);
else
    if min(bsbsteps)<init
        warning('FSDA:FSMbsb:WrongInit','It is impossible to monitor the subset for values smaller than init');
    end
    bsbsteps=bsbsteps(bsbsteps>=init);
    
    BB = NaN(n,length(bsbsteps),'single');
    WEI = NaN(n,length(bsbsteps));
end




% ij = index which is linked with the columns of matrix BB. During the
% search every time a subset is stored inside matrix BB ij icreases by one
ij=1;


% Hetero = (n-init1+1) x (r+1) matrix which will contain:
% 1st col = fwd search index
% 2nd col = estimate of first coeff of scedastic equation
%...
% (r+1) col = estimate of last coeff of scedastic equation
Hetero=[(init:n)' NaN(n-init+1,size(Z,2)+1)];


%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init+1:n)' , NaN(n-init,10));

hhh=1;
%% Start of the forward search
if (rank(Xb)~=p)
    warning('FSDA:FSRHmdr:message','Supplied initial subset does not produce full rank matrix');
    warning('FSDA:FSRHmdr:message','FS loop will not be performed');
    mdr=NaN;
    % FS loop will not be performed
else
    for mm=ini0:n;
        
        % if n>200 show every 100 steps the fwd search index
        if msg==1 && n>5000;
            if length(intersect(mm,seq100))==1;
                disp(['m=' int2str(mm)]);
            end
        end
        
        NoRankProblem=(rank(Xb) == p);
        if NoRankProblem  % rank is ok
            if art==1
                if  mm > 5  && gridsearch ~=1
                    % Use scoring
                    HET=regressHart(yb,Xb,Zb,'nocheck',1,'intercept',intercept);
                else
                    if size(Zb,2)==1
                        % Use grid search algorithm if Z has just one column
                        HET=regressHart_grid(yb,Xb,exp(Zb),'nocheck',1,'intercept',intercept);
                    else
                        HET=regressHart(yb,Xb,Zb,'nocheck',1,'intercept',intercept);
                    end
                end
                
                % gam=HET.GammaOLD;
                % alp=HET.alphaOLD;
                % omegahat=1+real(X(:,end).^alp)*gam;%equaz 6 di paper
                
                omegahat=1+exp(HET.Gamma(1,1))*exp(Z*HET.Gamma(2:end,1));
            else
                if  mm > 5  && gridsearch ~=1
                    % Use scoring
                    HET=regressHhar(yb,Xb,Zb,'intercept',intercept,'nocheck',1);
                else
                    if size(Zb,2)==1
                        % Use grid search algorithm if Z has just one column
                        HET=regressHhar_grid(yb,Xb,exp(Zb),'intercept',intercept,'nocheck',1);
                    else
                        HET=regressHhar(yb,Xb,Zb,'intercept',intercept,'nocheck',1);
                    end
                end
                
                
                omegahat=exp(Z*HET.Gamma(2:end,1));
                
            end
            
            sqweights = omegahat.^(-0.5);
            
            % Xw and yw are referred to all the observations
            % They contains transformed values of X and y using estimates
            % at step m
            % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
            Xw = bsxfun(@times, X, sqweights);
            yw = y .* sqweights;
            Xb=Xw(bsb,:);
            yb=yw(bsb);
            % The instruction below should not be necessary
            % Zb=Z(bsb,:);
            
            % b=Xb\yb;   % HHH
            b=HET.Beta(:,1);
            resBSB=yb-Xb*b;
        else   % number of independent columns is smaller than number of parameters
            error('FSDA:FSRHmdr:NoFullRank','Not full rank stop')
        end
        % HHH
        if hhh==1
            e=yw-Xw*b;  % e = vector of residual for all units using b estimated using subset
        else
            e=y-X*b;  % e = vector of residual for all units using b estimated using subset
        end
        
        r(:,2)=e.^2;
        
        if (mm>=init);
            
            % Store units belonging to the subset and weights
            if intersect(mm,bsbsteps)==mm
                % Store units belonging to subset
                BB(bsb,ij)=bsb;
                
                % Store weights
                WEI(:,ij)=sqweights;
                
                ij=ij+1;
            end
            
            
            if NoRankProblem
                % Store beta coefficients if there is no rank problem
                Bgls(mm-init+1,2:p+1)=b';
                % Store parameters of the scedastic equation
                Hetero(mm-init+1,2:end)=HET.Gamma(:,1)';
                
                % Compute and store estimate of sigma^2
                % using y and X transformed
                Sb=(resBSB)'*(resBSB)/(mm-p);
                S2(mm-init+1,2)=Sb;
                % Store R2
                S2(mm-init+1,3)=1-var(resBSB)/var(yb);
                
                if mm<n
                    mAm=Xb'*Xb; % HHH
                    
                    % Take minimum deletion residual among noBSB
                    % hi (the leverage for the units not belonging to the
                    % subset) is defined as follows
                    % hi=diag(X(ncl,:)*inv(Xb'*Xb)*(X(ncl,:))');
                    
                    % Take units not belonging to bsb
                    if hhh==1
                        Xncl = Xw(ncl, :);
                    else
                        Xncl = X(ncl,:); % HHH
                    end
                    % mmX=inv(mAm);
                    % hi = sum((Xncl*mmX).*Xncl,2);
                    hi=sum((Xncl/mAm).*Xncl,2);
                    
                    if hhh==1
                        ord = [(r(ncl,2)./(1+hi)) e(ncl)];
                    else
                        ord = [(r(ncl,2)./(omegahat(ncl)+hi)) e(ncl)];
                    end
                    
                    % Store minimum deletion residual in matrix mdr
                    selmdr=sortrows(ord,1);
                    if S2(mm-init+1,2)==0
                        warning('FSDA:FSRHmdr:ZeroS2','Value of S2 at step %d is zero, mdr is NaN',mm-init+1);
                    else
                        mdr(mm-init+1,2)=sqrt(selmdr(1,1)/HET.sigma2);
                    end
                end  %if mm<n
            end   %~RankProblem
        end     %mm>=init1
        
        if mm<n;
            
            % store units forming old subset in vector oldbsb
            oldbsb=bsb;
            
            % order the r_i
            % r(:,2)=(y-X*b).^2;  % e = vector of residual for all units using b estimated using subset
            
            
            % units inside vector constr are forced to join the search in
            % the final k steps
            if ~isempty(constr) && mm<n-length(constr)
                r(constr,2)=Inf;
            end
            ord=sortrows(r,2);
            
            % bsb= units forming the new  subset
            bsb=ord(1:(mm+1),1);
            
            Xb=X(bsb,:);  % subset of X     HHH  OK
            yb=y(bsb);    % subset of y     HHH   OK
            Zb=Z(bsb,:);  % subset of Z
            
            if mm>=init;
                unit=setdiff(bsb,oldbsb);
                
                % If the interchange involves more than 10 units, store only the
                % first 10.
                if length(unit)<=10
                    Un(mm-init+1,2:(length(unit)+1))=unit;
                else
                    if msg==1
                        disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        disp(['Number of units which entered=' int2str(length(unit))]);
                        Un(mm-init+1,2:end)=unit(1:10);
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
        [gmin] = FSRenvmdr(n,p,'prob',quant,'init',init,'exact',1);
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
%FScategory:REG-Hetero