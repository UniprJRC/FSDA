function [Bmdr,Un,BB,BBayes,S2Bayes] = FSRBmdr(y, X, beta0, R, tau0, n0, bsb, varargin)
%FSRBmdr computes minimum deletion residual and other basic linear regression 
%quantities in each step of the search.
%
%
%<a href="matlab: docsearch('FSRBmdr')">Link to the help function</a>
%
% Required input arguments2
%
%  y:            A vector with n elements that contains the response variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%  X :           Data matrix of explanatory variables (also called
%               'regressors') of dimension (n x p-1).
%               Rows of X represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
%  bsb :         list of units forming the initial subset, if bsb=0
%               (default) then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
%
% Optional input arguments:
%
%  init :        scalar, specifies the point where to initialize the search
%               and start monitoring required diagnostics. If it is not
%               specified it is set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%  intercept :   If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%  plots :       If equal to one a plot of minimum deletion residual
%               appears  on the screen with 1%, 50% and 99% confidence
%               bands else (default) no plot is shown.
%               Remark: the plot which is produced is very simple. In order
%               to control a series of options in this plot and in order to
%               connect it dynamically to the other forward plots it is necessary to use
%               function mdrplot
%  nocheck:      Scalar. If nocheck is equal to 1 no check is performed on
%               matrix y and matrix X. Notice that y and X are left
%               unchanged. In other words the additioanl column of ones for
%               the intercept is not added. As default nocheck=0. The
%               controls on h, alpha and nsamp still remain
%  msg  :        scalar which controls whether to display or not messages
%               about great interchange on the screen
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen
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
%  Bmdr:        n -init x 2 matrix which contains the monitoring of minimum
%               deletion residual at each step of the forward search.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum deletion residual.
%  Un:           (n-init) x 11 Matrix which contains the unit(s) included
%               in the subset at each step of the search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
%  BB:           n x (n-init+1) matrix which the units belonging to the
%               subset at each step of the forward search.
%               1st col = index forming subset in the initial step
%               ...
%               last column = units forming subset in the final step (i.e.
%               all units).
%  BBayes:      (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward
%               search.
%  S2Bayes:     (n-init+1) x 3 matrix containing the monitoring of S2 (2nd
%               column)and R2 (third column) in each step of the forward
%               search.
% See also
%
% References:
%
%   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
%   Springer Verlag, New York.
%   Atkinson, A.C. and Riani, M. (2006). Distribution theory and
%   simulations for tests of outliers in regression. Journal of
%   Computational and Graphical Statistics, Vol. 15, pp. 460–476
%   Riani, M. and Atkinson, A.C. (2007). Fast calibrations of the forward
%   search for testing multiple outliers in regression, Advances in Data
%   Analysis and Classification, Vol. 1, pp. 123–141.
%
% Copyright 2008-2015.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti
%            and Vytis Kopustinskas (2009-2010)
%
%<a href="matlab: docsearch('FSRBmdr')">Link to the help function</a>
% Last modified 15-Nov-2011

% Examples:

%{
%Common part to all examples: load fishery dataset
%
 load('fishery');
 y=fishery.data(:,1);
 X=fishery.data(:,2);
 [out]=LXS(y,X,'nsamp',10000);
%}

%{
% FSR with all default options.
[mdr,Un,BB,Bols,S2] = FSRBmdr(y,X,out.bs);
%}

%{
% FSRmdr monitoring from step 60.
[mdr,Un,BB,Bols,S2] = FSRBmdr(y,X,out.bs,'init',60);
%}

%{
% FSRmdr using a regression model without intercept.
[mdr,Un,BB,Bols,S2] = FSRBmdr(y,X,out.bs,'intercept','0');
%}



%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

init=0;
options=struct('intercept',1,'init',init,'plots',0,'nocheck',0,'msg',1);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('Error:: number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin<3
    error('Initial subset is missing');
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
    bsb=randsample(n,p);
    Xb=X(bsb,:);
    yb=y(bsb);
else
    Xb=X(bsb,:);
    yb=y(bsb,:);
end

ini0=numel(bsb);

% check init
init1=options.init;
if  init1 <0;
    mess=sprintf(['Attention : init1 should be larger than 0. \n',...
        'It is set to 0.']);
    disp(mess);
    init1=0;
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

% Matrix BB will contain the beta coefficients in each step of the fwd
% search. The first row will contain the units forming initial subset
BBayes=[(init1:n)' NaN(n-init1+1,p)];     %initial value of beta coefficients is set to NaN
% initialize the space on the SE array with NaNs
% S2 = (n-init1+1) x 3 matrix which will contain:
% 1st col = fwd search index
% 2nd col = S2= \sum e_i^2 / (m-p)
% 3rd col = R^2
S2Bayes=[(init1:n)' NaN(n-init1+1,2)];        %initial value of S2 (R2) is set to NaN

% mdr = (n-init1-1) x 2 matrix which will contain min deletion residual
% among nobsb r_i^*
Bmdr=[(init1:n-1)'  NaN(n-init1,1)];      %initial value of mdr is set to NaN

% Matrix BB will contain the units forming subset in each step of the
% forward search. The first column contains the units forming subset at
% step init1
BB = NaN(n,n-init1+1);


%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init1+1:n)' , NaN(n-init1,10));

%% Start of the forward search

%mj=1;
for mm=ini0:n;
    
    % if n>200 show every 100 steps the fwd search index
%     if msg==1 && n>2000;
%         if length(intersect(mm,seq100))==1;
%             disp(['m==' int2str(mm)]);
%             % figure;
%             offset1=1;
%         end
%     end
    
    % call bayesian procedure
    [bayes]=regressB(y, X(:,2:end), beta0, R, tau0, n0, 'bsb', bsb);
    
    % bayesian beta
    b=bayes.beta1;
    
    % hpdi
    pos=mm-init+1;
    
    % please note:
    % FS uses standard residuals to advance but
    % at each iteration saves the bayesian residuals
    
    resBSB=yb-Xb*b;
    
    e=y-X*b;    % e = vector of residual for all units using b estimated using subset
    % precalc mAm=X0'*X0 + Xb'*Xb , where R=X0'*X0
    
    r(:,2)=e.^2;
    
    if (mm>=init1);
        % store Units belonging to the subset
        BB(bsb,mm-init1+1)=bsb;
        
        % Stores beta coefficients if there is no rank problem
        BBayes(mm-init1+1,2:p+1)=b';
        % Compute and store estimate of sigma^2
        
        % Sb=(resBSB)'*(resBSB)/(mm-p);
        Sb=1/bayes.tau1;
        
        S2Bayes(mm-init1+1,2)=Sb;
        % Store R2
        S2Bayes(mm-init1+1,3)=1-var(resBSB)/var(yb);
        
        if mm<n
            % Store Bayesian mdr
            Bmdr(mm-init1+1,2)= min(abs(bayes.res(ncl,2)));
        end
        
    end
    
    if mm<n;
        
        % store units forming old subset in vector oldbsb
        oldbsb=bsb;
        
        % order the r_i
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
            ncl=ord(mm+2:n,1);    % ncl= units forming the new noclean
        end
    end   % if mm<n
    %mj=mj+n0;
end  % for mm=ini0:n loop

end