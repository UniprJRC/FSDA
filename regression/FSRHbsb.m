function [Un,BB] = FSRHbsb(y,X,Z,bsb,varargin)
%FSRHbsb returns the units belonging to the subset in each step of the heteroskedastic forward search
%
%<a href="matlab: docsearchFS('FSRHbsb')">Link to the help function</a>
%
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
%               Therefore, if for example the explanatory variables
%               responsible for heteroscedasticity are columns 3 and 5
%               of matrix X, it is possible to use both the sintax:
%                    FSRHbsb(y,X,X(:,[3 5]),0)
%               or the sintax:
%                    FSRHbsb(y,X,[3 5],0)
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
%  intercept :  Indicator for constant term. true (default) | false. 
%               Indicator for the constant term (intercept) in the fit,
%               specified as the comma-separated pair consisting of
%               'intercept' and either true to include or false to remove
%               the constant term from the model.
%               Example - 'intercept',false
%               Data Types - boolean
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
%               If msg==1 (default) messages are displayed on the screen
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
%  constr   :   units which are forced to join the search in the last r steps. Vector.
%               r x 1 vector. The default is constr=''.  No constraint is imposed
%               Example - 'constr',[1 6 3]
%               Data Types - double
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
%   plots   :   Plot on the screen. Scalar.
%               If plots=1 the monitoring units plot is displayed on the
%               screen. The default value of plots is 0 (that is no plot
%               is produced on the screen).
%               Example - 'plots',1
%               Data Types - double
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
%  BB:          Units belonging to subset in each step. Matrix.
%               n x (n-init+1) matrix which contains the units belonging to the
%               subset at each step of the forward search.
%               1st col = index forming subset in the initial step
%               ...
%               last column = units forming subset in the final step (i.e.
%               all units).
%
%
% See also:   FSRbsb, FSRBbsb
%
%
% References:
%
% Atkinson, A.C., Riani, M. and Torti, F. (2016), Robust methods for
% heteroskedastic regression, "Computational Statistics and Data Analysis",
% Vol. 104, pp. 209-222, http://dx.doi.org/10.1016/j.csda.2016.07.002 [ART]

%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSRHbsb')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % FSRHbsb with all default options.
    % Common part to all examples: load tradeH dataset (used in the paper ART).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    Un=FSRHbsb(y,X,Z,[1:10]);
%}

%{
    %% FSRHbsb with optional arguments.
    % Suppress all messages about interchange with option msg.
    % Common part to all examples: load tradeH dataset (used in the paper ART).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    Un=FSRHbsb(y,X,Z,[1:10],'plots',1,'msg',0);
%}


%{
    % Monitoring the units belonging to subset in each step.
    % Common part to all examples: load tradeH dataset (used in the paper ART).
    XX=load('tradeH.txt');
    y=XX(:,2);
    X=XX(:,1);
    X=X./max(X);
    Z=log(X);
    [~,Un,BB]=FSRHmdr(y,X,Z,[1:10]);
    [Unchk,BBchk]=FSRHbsb(y,X,Z,[1:10]);
    % Test for equality BB and BBchk
    disp(isequaln(BB,BBchk))
    % Test for equality Un and Unchk
    disp(isequaln(Un,Unchk))
%}


%% Beginning of code

% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

if n<40
    initdef=p+1;
else
    initdef=min(3*p+1,floor(0.5*(n+p+1)));
end

if initdef<length(bsb)
    initdef=length(bsb);
end

bsbstepdef='';

options=struct('intercept',1,'init',initdef,'plots',0,'nocheck',0,'msg',1,...
    'constr','','modeltype','art','gridsearch',0,'bsbsteps',bsbstepdef);

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
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end


% Z = n-by-r matrix which contains the explanatory variables for
% heteroskedasticity
if size(Z,1)~=n
    % Check if interecept was true
    intercept=options.intercept;
    if intercept==1
        Z=X(:,Z+1);
    else
        Z=X(:,Z);
    end
end


if bsb==0
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
if  init <p+1
    fprintf(['Attention : init should be larger than p. \n',...
        'It is set to p+1.']);
    init=p+1;
elseif init<ini0
    fprintf(['Attention : init should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    init=ini0;
elseif init>=n
    fprintf(['Attention : init should be smaller than n. \n',...
        'It is set to n-1.']);
    init=n-1;
end

msg=options.msg;
constr=options.constr;
intercept=options.intercept;

%% Initialise key matrices


% sequence from 1 to n.
seq=(1:n)';

% The second column of matrix R will contain the OLS residuals at each step
% of the search
r=[seq zeros(n,1)];

% If n is very large, the step of the search is printed every 100 step
% seq100 is linked to printing
seq100=100*(1:1:ceil(n/100));

bsbsteps=options.bsbsteps;
% Matrix BB will contain the units forming subset in each step (or in
% selected steps) of the forward search. The first column contains
% information about units forming subset at step init1.
if isempty(bsbsteps)
    % Default for vector bsbsteps which indicates for which steps of the fwd
    % search units forming subset have to be saved
    if n<=5000
        bsbsteps = init:1:n;
    else
        bsbsteps = [init init+100-mod(init,100):100:100*floor(n/100)];
    end
    BB = NaN(n,length(bsbsteps),'single');
elseif bsbsteps==0
    bsbsteps=init:n;
    BB = NaN(n,n-init+1,'single');
else
    if min(bsbsteps)<init
        warning('FSDA:FSMbsb:WrongInit','It is impossible to monitor the subset for values smaller than init');
    end
    bsbsteps=bsbsteps(bsbsteps>=init);
    
    BB = NaN(n,length(bsbsteps),'single');
end

%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init+1:n)' , NaN(n-init,10));

hhh=1;
%% Start of the forward search
if (rank(Xb)~=p)
    warning('FSDA:FSRHmdr:message','Supplied initial subset does not produce full rank matrix');
    warning('FSDA:FSRHmdr:message','FS loop will not be performed');
    % FS loop will not be performed
else
    % ij = index which is linked with the columns of matrix BB. During the
    % search every time a subset is stored inside matrix BB ij icreases by one
    ij=1;
    
    for mm=ini0:n
        
        % if n>200 show every 100 steps the fwd search index
        if msg==1 && n>5000
            if length(intersect(mm,seq100))==1
                disp(['m=' int2str(mm)]);
            end
        end
        
        NoRankProblem=(rank(Xb) == p);
        if NoRankProblem  % rank is ok
            if art==1
                if  mm > 5  && gridsearch ~=1
                    % Use scoring
                    HET=regressHart(yb,Xb(:,2:end),Zb,'intercept',intercept);
                else
                    if size(Zb,2)==1
                        % Use grid search algorithm if Z has just one column
                        HET=regressHart_grid(yb,Xb(:,2:end),exp(Zb),'intercept',intercept);
                    else
                        HET=regressHart(yb,Xb(:,2:end),Zb,'intercept',intercept);
                    end
                end
                
                % gam=HET.GammaOLD;
                % alp=HET.alphaOLD;
                % omegahat=1+real(X(:,end).^alp)*gam;%equaz 6 di paper
                
                omegahat=1+exp(HET.Gamma(1,1))*exp(Z*HET.Gamma(2:end,1));
            else
                if  mm > 5  && gridsearch ~=1
                    % Use scoring
                    HET=regressHhar(yb,Xb(:,2:end),Zb,'intercept',intercept);
                else
                    if size(Zb,2)==1
                        % Use grid search algorithm if Z has just one column
                        HET=regressHhar_grid(yb,Xb(:,end),exp(Zb),'intercept',intercept);
                    else
                        HET=regressHhar(yb,Xb(:,2:end),Zb,'intercept',intercept);
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
            b=HET.Beta(:,1);
        else   % number of independent columns is smaller than number of parameters
            error('FSDA:FSRHbsb:NoFullRank','Not full rank stop')
        end
        % HHH
        if hhh==1
            e=yw-Xw*b;  % e = vector of residual for all units using b estimated using subset
        else
            e=y-X*b;  % e = vector of residual for all units using b estimated using subset
        end
        
        r(:,2)=e.^2;
        
        % Store units belonging to the subset
        if (mm>=init)
            if intersect(mm,bsbsteps)==mm
                BB(bsb,ij)=bsb;
                ij=ij+1;
            end
        end
        
        
        
        if mm<n
            
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
            
            Xb=X(bsb,:);  % subset of X
            yb=y(bsb);    % subset of y
            Zb=Z(bsb,:);  % subset of Z
            
            if mm>=init
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
        end   % if mm<n
    end  % for mm=ini0:n loop
    
end % rank check

plots=options.plots;
if plots==1
    % Create the 'monitoring units plot'
    figure;
    plot(bsbsteps,BB','bx')
    xlabel('Subset size m');
    ylabel('Monitoring units plot');
end

end
%FScategory:REG-Hetero