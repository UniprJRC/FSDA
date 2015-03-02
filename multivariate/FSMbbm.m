function BBmsel = FSMbbm(Y,bsb,msel,varargin)
%FSMbbm gives the units belonging to subset at step(s) msel of the forward search
%
%<a href="matlab: docsearchFS('FSMbbm')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Y = n x p data matrix; n observations
%               and p variables
%               Rows of Y represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
% bsb :         list of units forming the initial subset, if bsb=0
%               (default) then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
% msel :        vector which specifies for which steps of the search the
%               units forming subset must be saved
%
%
% Optional input arguments:
%
% init :        scalar, specifies the point where to start monitoring
%               required diagnostics. Note that if bsb is supplied
%               init>=length(bsb). If init is not specified it will
%               be set equal to floor(n*0.6).
% plots :       scalar. If plots=1, a plot of the monitoring of minMD among
%               the units not belonging to the subset is produced on the
%               screen with 1% 50% and 99% confidence bands
%               else (default) no plot is produced.
%  msg  :       scalar which controls whether to display or not messages
%               about great interchange on the screen
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen
% nocheck :     Scalar. If nocheck is equal to 1 no check is performed on
%               matrix Y. As default nocheck=0.
%
% Remark:       The user should only give the input arguments that have to
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
%
% BBmsel  :     Matrix of size n-times-length(msel) containing units
%               belonging to subset in the steps of
%               the search defined by input vector msel
%               More precisely
%               BBmsel(:,1) contains the units forming subset in step mmsel(1)
%               ....
%               BBmsel(:,end) contains the units forming subset in step mmsel(end)
%               Row 1 of matrix BBmsel is referred to unit 1
%               ......
%               Row n of matrix BBmsel is referred to unit n
%               Units not belonging to subset are denoted with NaN
%
% See also FSMeda, FSM.m, FSMmmd
%
% References:
%
%   Atkinson Riani and Cerioli (2004), Exploring multivariate data with the
%   forward search Springer Verlag, New York.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSMbbm')">Link to the help function</a>
% Last modified 06-Feb-2015


% Examples:

%{
    % Run this code to see the output shown in the help file
    n=200;
    v=3;
    randn('state',123456);
    Y=randn(n,v);
    %Contaminated data
    Ycont=Y;
    seqcont=11:10:51;
    disp('Contaminated units')
    disp(seqcont)
    Ycont(seqcont,:)=Ycont(seqcont,:)+2.5;
    % Analyse the units forming subset in step msel=195
    msel=195;
    mmsel=FSMbbm(Ycont,0,msel);
    seq=1:n;
    disp(['Units outside subset at step m=' num2str(msel)])
    disp(setdiff(seq,mmsel))
%}

%% Beginning of code
% Input parameters checking
%chkinputM does not do any check if option nocheck=1
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);

%% Input parameters checking
[n,v]=size(Y);

hdef=floor(n*0.6);

options=struct('init',hdef,'plots',0,'msg',1,'nocheck',0);

if nargin<2
        error('FSDA:FSMbbm:missingInputs','Initial subset is missing')
end

if nargin<3
    error('FSDA:FSMbbm:missingInputs','Vector which contains the selected steps of the search is missing');
end

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSMbbm:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

%init1=options.init;
if nargin > 2
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

% And check if the optional user parameters are reasonable.

if bsb==0;
    Ra=1; nwhile=1;
    while and(Ra,nwhile<100)
        % Extract a random sample of size v+1
        bsb=randsample(n,v+1);
        % Check if the var-cov matrix of the random sample is full (i.e =v)
        Ra=(rank(cov(Y(bsb,:)))<v);
        nwhile=nwhile+1;
    end
    if nwhile==100
        warning('FSDA:FSMbbm:NoFullRank','Unable to randomly sample full rank matrix');
    end
else
end


ini0=length(bsb);


if  min(msel) <ini0+1
    mess=sprintf(['Attention : min(mmsel) should be larger than  \n',...
        'length(bsb)']);
    disp(mess);
    error('FSDA:FSMbbm:Wrongmmsel','Redefine input argument mmsel')
elseif msel>=n;
    mess=sprintf(['Attention : max(mmsel) should be smaller than n. \n',...
        'It is set to n-1.']);
    disp(mess)
    error('FSDA:FSMbbm:Wrongmmsel','Redefine input argument mmsel')
end


% Initialize matrix which will contain the units forming subset in each
% step of the fwd search
BBmsel=NaN(n,length(msel),'single');

if (rank(Y(bsb,:))<v)
    warning('FSDA:FSMbbm:NoFullRank','The supplied initial subset is not full rank matrix');
    % FS loop will not be performed
else
    
    jj=1;
    for mm = ini0:max(msel)
        
        % select subset
        Yb=Y(bsb,:);
        
        if max(mm==msel);
            BBmsel(bsb,jj)=bsb;
            jj=jj+1;
        end
        
        % Find vector of means inside subset
        % Note that ym is a row vector
        ym=mean(Yb);
        
        % Squared Mahalanobis distances computed using QR decomposition
        % Ym=Y-one*ym;
        Ym = bsxfun(@minus,Y, ym);
        
        [~,R]=qr(Ym(bsb,:),0);
        
        
        % Remark: u=(Ym/R)' should be much faster than u=inv(R')*Ym';
        u=(Ym/R);
        % Compute square Mahalanobis distances
        mala=((mm-1)*sum(u.^2,2));
        
        [~,zs]= sort(mala);
        
        if mm<n
            
            % the dimension of subset increases by one unit.
            % vector bsb contains the indexes corresponding to the units of
            % the new subset
            bsb=zs(1:mm+1,1);
            
        else
        end
    end
end


