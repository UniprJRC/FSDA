function [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,bsb,varargin)
%FSRmdr computes minimum deletion residual and other basic linear regression quantities in each step of the search
%
%<a href="matlab: docsearch('FSRmdr')">Link to the help function</a>
%
% Required input arguments:
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
%  init :       scalar, specifies the point where to initialize the search
%               and start monitoring required diagnostics. If it is not
%               specified it is set equal to:
%                   p+1, if the sample size is smaller than 40;
%                   min(3*p+1,floor(0.5*(n+p+1))), otherwise.
%  intercept :  If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%  plots :      If equal to one a plot of minimum deletion residual
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
%  msg  :       scalar which controls whether to display or not messages
%               about great interchange on the screen
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen
%  constr :     r x 1 vector which contains the list of units which are
%               forced to join the search in the last r steps. The default
%               is constr=''.  No constraint is imposed
% bsbmfullrank :scalar which tells how to behave in case subset at step m
%               (say bsbm) produces a non singular X. In other words,
%               this options controls what to do when rank(X(bsbm,:)) is
%               smaller then number of explanatory variables. If
%               bsbmfullrank = 1 (default is 1) these units (whose number is
%               say mnofullrank) are constrained to enter the search in
%               the final n-mnofullrank steps else the search continues
%               using as estimate of beta at step m the estimate of beta
%               found in the previous step.
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
%  mdr:          n -init x 2 matrix which contains the monitoring of minimum
%               deletion residual at each step of the forward search.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum deletion residual.
%               REMARK: if in a certain step of the search matrix is
%               singular, this procedure checks ohw many observations
%               produce a singular matrix. In this case mdr is a column
%               vector which contains the list of units for which matrix X
%               is non singular.
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
%  Bols:         (n-init+1) x (p+1) matrix containing the monitoring of
%               estimated beta coefficients in each step of the forward
%               search.
%  S2:           (n-init+1) x 3 matrix containing the monitoring of S2 (2nd
%               column)and R2 (third column) in each step of the forward
%               search.
%
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
% Copyright 2008-2014.
% Written by FSDA team
%
%
%<a href="matlab: docsearch('FSRmdr')">Link to the help function</a>
% Last modified 08-Dec-2013

% Examples:

%{
%Common part to all examples: load fishery dataset
%
 load('fishery');
 y=fishery.data(:,2);
 X=fishery.data(:,1);
 [out]=LXS(y,X,'nsamp',10000);
%}

%{
% FSR with all default options.
[mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs);
%}

%{
% FSRmdr monitoring from step 60.
[mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs,'init',60);
%}

%{
% FSRmdr using a regression model without intercept.
[mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs,'intercept','0');
%}

%{
%FSRmdr applied without doing any checks on y and X variables.
[mdr,Un,BB,Bols,S2] = FSRmdr(y,X,out.bs);
%}



%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

if n<40
    init=p+1;
else
    init=min(3*p+1,floor(0.5*(n+p+1)));
end
options=struct('intercept',1,'init',init,'plots',0,'nocheck',0,'msg',1,...
    'constr','','bsbmfullrank',1);

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
    Ra=1; nwhile=1;
    while and(Ra,nwhile<100)
        bsb=randsample(n,p);
        Xb=X(bsb,:);
        Ra=(rank(Xb)<p);
        nwhile=nwhile+1;
    end
    if nwhile==100
        warning('FSRmdr:message','Unable to randomly sample full rank matrix');
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
Bols=[(init1:n)' NaN(n-init1+1,p)];     %initial value of beta coefficients is set to NaN

% S2 = (n-init1+1) x 3 matrix which will contain:
% 1st col = fwd search index
% 2nd col = S2= \sum e_i^2 / (m-p)
% 3rd col = R^2
S2=[(init1:n)' NaN(n-init1+1,2)];        %initial value of S2 (R2) is set to NaN

% mdr = (n-init1-1) x 2 matrix which will contain min deletion residual
% among nobsb r_i^*
mdr=[(init1:n-1)'  NaN(n-init1,1)];      %initial value of mdr is set to NaN

% Matrix BB will contain the units forming subset in each step of the
% forward search. The first column contains the units forming subset at
% step init1
BB = NaN(n,n-init1+1);

%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init1+1:n)' , NaN(n-init1,10));

%% Start of the forward search
if (rank(Xb)~=p)
    warning('FSRmdr:message','Supplied initial subset does not produce full rank matrix');
    warning('FSRmdr:message','FS loop will not be performed');
    mdr=NaN;
    % FS loop will not be performed
else
    for mm=ini0:n;
        
        
        
        % if n>200 show every 100 steps the fwd search index
        if msg==1 && n>500;
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
            BB(bsb,mm-init1+1)=bsb;
            
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
                        warning('FSRmdr:ZeroS2','Value of S2 at step %d is zero, mdr is NaN',mm-init1+1);
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