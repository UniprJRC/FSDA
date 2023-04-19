function out = FSMeda(Y,bsb,varargin)
%FSMeda performs forward search in multivariate analysis with exploratory data analysis purposes
%
%
%<a href="matlab: docsearchFS('FSMeda')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Input data. Matrix. n x v data matrix; n observations and v
%               variables. Rows of Y represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
%               Data Types - single | double
% bsb :         Units forming subset. Vector. List of units forming the initial subset.
%               If bsb=0 (default) then the procedure starts with v units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb).
%               Data Types - single | double
%
% Optional input arguments:
%
% init :       Point where to start monitoring required diagnostics. Scalar.
%              Note that if bsb is supplied, init>=length(bsb). If init is not
%              specified it will be set equal to floor(n*0.6).
%                 Example - 'init',50
%                 Data Types - double
%
% plots :    It specify whether it is necessary to produce the plots of the
%               monitoring of minMD.
%                 Scalar. If plots=1, a plot of the monitoring of minMD among
%               the units not belonging to the subset is produced on the
%               screen with 1 per cent, 50 per cent and 99 per cent confidence bands
%               else (default), all plots are suppressed.
%               Example - 'plots',0
%               Data Types - double
%
%  msg  :       It controls whether to display or not messages
%               about great interchange on the screen. Scalar.
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen.
%                 Example - 'msg',0
%                 Data Types - double
%
% scaled:     It controls whether to monitor scaled Mahalanobis distances.
%               Scalar.
%               If scaled=1  Mahalanobis distances
%               monitored during the search are scaled using ratio of determinant.
%               If scaled=2  Mahalanobis distances
%               monitored during the search are scaled using asymptotic consistency factor.
%               The default value is 0 that is Mahalanobis distances are
%               not scaled.
%                 Example - 'scaled',0
%                 Data Types - double
%
% nocheck     : It controls whether to perform checks on matrix Y.Scalar.
%               If nocheck is equal to 1 no check is performed on
%               matrix Y. As default nocheck=0.
%                 Example - 'nocheck',1
%                 Data Types - double
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
%               computations.
%
% Output:
%
%         out:   structure which contains the following fields
%
%   out.MAL=        n x (n-init+1) = matrix containing the monitoring of
%               Mahalanobis distances.
%               1st row = distance for first unit;
%               ...;
%               nth row = distance for nth unit.
%    out.BB=        n x (n-init+1) matrix containing the information about the units belonging
%               to the subset at each step of the forward search.
%               1st col = indexes of the units forming subset in the
%               initial step;
%               ...;
%               last column = units forming subset in the final step (all
%               units).
%   out.mmd=        n-init x 3 matrix which contains the monitoring of minimum
%               MD or (m+1)th ordered MD  at each step of
%               the forward search.
%               1st col = fwd search index (from init to n-1);
%               2nd col = minimum MD;
%               3rd col = (m+1)th-ordered MD.
%   out.msr=        n-init+1 x 3 = matrix which contains the monitoring of
%               maximum MD or mth ordered MD.
%               1st col = fwd search index (from init to n);
%               2nd col = maximum MD;
%               3rd col = mth-ordered MD.
%    out.gap=       n-init+1 x 3 = matrix which contains the monitoring of
%               the gap (difference between minMD outside subset and max.
%               inside).
%               1st col = fwd search index (from init to n);
%               2nd col = min MD - max MD;
%               3rd col = (m+1)th ordered MD - mth ordered distance.
%   out.Loc=        (n-init+1) x (v+1) matrix containing the monitoring of
%               estimated of the means for each variable in each step of
%               the forward search.
%  out.S2cov=       (n-init+1) x (v*(v+1)/2+1) matrix containing the monitoring
%               of the elements of the covariance matrix in each step
%               of the forward search.
%               1st col = fwd search index (from init to n);
%               2nd col = monitoring of S(1,1);
%               3rd col = monitoring of S(1,2);
%               ...;
%               end col = monitoring of S(v,v).
%  out.detS=        (n-init+1) x (2) matrix containing the monitoring of
%               the determinant of the covariance matrix
%               in each step of the forward search.
%    out.Un=        (n-init) x 11 Matrix which contains the unit(s)
%               included in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,2) for example contains
%               the unit included in step init+1 Un(end,2) contains the
%               units included in the final step of the search
%     out.Y=        Original data input matrix
% out.class=    'FSMeda'
%
% See also FSMmmd.m, FSM.m
%
% References:
%
% Atkinson, A.C., Riani, M. and Cerioli, A. (2004), "Exploring multivariate
% data with the forward search", Springer Verlag, New York.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSMeda')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    % FSMeda with all default options.
    % Run the FS on a simulated dataset by choosing an initial subset
    % formed by the three observations with the smallest Mahalanobis
    % Distance.
    n=100;
    v=3;
    m0=4;
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [fre]=unibiv(Y);
    %create an initial subset with the 3 observations with the lowest
    %Mahalanobis Distance
    fre=sortrows(fre,4);
    bs=fre(1:m0,1);
    [out]=FSMeda(Ycont,bs);
%}

%{
    %% FSMeda with optional arguments.
    % Monitoring the evolution of minimum Mahalanobis distance.
    n=100;
    v=3;
    m0=3;
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [fre]=unibiv(Y);
    %create an initial subset with the 3 observations with the lowest
    %Mahalanobis Distance
    fre=sortrows(fre,4);
    bs=fre(1:m0,1);
    [out]=FSMeda(Ycont,bs,'plots',1);
%}

%{
    %% Example with the Swiss bank notes data.
    load('swiss_banknotes');
    Y=swiss_banknotes{:,:};
    [fre]=unibiv(Y);
    %create an initial subset with the 3 observations with the lowest
    %Mahalanobis Distance
    fre=sortrows(fre,4);
    m0=20;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'plots',1,'init',30);
%}

%{
    % Example with the Emilia Romagna data.
    load('emilia2001')
    Y=emilia2001{:,:};
    [fre]=unibiv(Y);
    %create an initial subset with the 30 observations with the lowest
    %Mahalanobis Distance
    fre=sortrows(fre,4);
    m0=30;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y,bs,'init',100);
    % Minimum Mahalanobis distance
    % Compare the plot with Figure 1.12 p. 21, ARC (2004)
    mmdplot(out,'ylimy',[6 14])
    % Analysis of the last 16 units to enter the forward search
    % Compare the results with Table 1.3 p. 21
    disp(out.Un(end-15:end,:));
%}


%{
    % Example with the Emilia Romagna data (all variables).
    load('emilia2001')
    Y=emilia2001{:,:};
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify variables y16 y23 y25 y26
    sel=[16 23 25 26];
    sel=[25 26];
    Y1=Y;
    Y1(:,sel)=100-Y1(:,sel);
                        
    la0demo=[0.5,0.25,0,1,0.25,0,0,0.25,0.5];
    la0weal=[0.25,0.5,0.5,1,1,0.5,-1/3,0.25,0.25,-1];
    la0work=[0.25,0,1,0,0,0.25,1,1,1];
    la0C2=[la0demo(1:5) la0work(1:4) la0demo(6:9) la0weal la0work(5:9)];
    Y1tr=normBoxCox(Y1,1:28,la0C2);
    [fre]=unibiv(Y1tr);
    %create an initial subset with the 30 observations with the lowest
    %Mahalanobis Distance
    fre=sortrows(fre,4);
    m0=30;
    bs=fre(1:m0,1);
    [out]=FSMeda(Y1tr,bs,'init',100,'scaled',1);
    % Minimum Mahalanobis distance
    [out]=FSMeda(Y1tr,bs,'init',100);
    mmdplot(out,'ylimy',[5 26])
    
    standard=struct;
    standard.ylim=[4 17];
    malfwdplot(out,'standard',standard);
%}


%% Beginning of code 

% Input parameters checking
%chkinputM does not do any check if option nocheck=1
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);

if nargin<2
    error('FSDA:FSMeda:missingInputs','Initial subset is missing')
end

[n,v]=size(Y);

hdef=floor(n*0.6);
options=struct('init',hdef,'plots',0,'msg',1,'scaled',0,'nocheck',0);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSMeda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

%init1=options.init;
if nargin > 2
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

scaled=options.scaled;

verLessThan2016b=verLessThan('matlab','9.1');

% And check if the optional user parameters are reasonable.

if bsb==0
    Ra=1; nwhile=1;
    while and(Ra,nwhile<100)
        % Extract a random sample of size v+1
        bsb=randsample(n,v+1);
        % Check if the var-cov matrix of the random sample is full (i.e =v)
        Ra=(rank(cov(Y(bsb,:)))<v);
        nwhile=nwhile+1;
    end
    if nwhile==100
        warning('FSDA:FSMeda:NoFullRank','Unable to randomly sample full rank matrix');
    end
else
end

ini0=length(bsb);

% check init
init1=options.init;
msg=options.msg;

if  init1 <v+1
    mess=sprintf(['Attention : init1 should be larger than v. \n',...
        'It is set to v+1.']);
    fprintf('%s\n',mess);
    init1=v+1;
elseif init1<ini0
    mess=sprintf(['Attention : init1 should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(length(bsb)) ]);
    fprintf('%s\n',mess);
    init1=ini0;
elseif init1>=n
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    fprintf('%s\n',mess);
    init1=n-1;
end


% Initialize matrix which will contain the units forming subset in each
% step of the fwd search
BB=NaN(n,n-init1+1);

% Initialize matrix which will contain the MD monitored in each
% step of the fwd search
MAL=BB;

% Initialize matrix which will contain the means of the variables monitored in
% each step
loc=cat(2,(init1:n)',NaN(n-init1+1,v));

%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2,(init1+1:n)',NaN(n-init1,10));

%  mmd has three columns
%  1st col = dimension of the subset
%  2nd col min. MD among the units ont belonging to the subset
%  3rd col (m+1) ordered MD
mmd=[(init1:n-1)' zeros(n-init1,2)];

%  msr has three columns
%  1st col = dimension of the subset
%  2nd col max. MD among the units  belonging to the subset
%  3rd col m-th ordered MD
msr=[(init1:n)' zeros(n-init1+1,2)];

%  msr has three columns
%  1st col = dimension of the subset
%  2nd col min. MD among the units not belonging to the subset - max inside
%  3rd col (m+1)-th ordered MD - m-th ordered MD
gap=mmd;

% Initialize matrix which contains the monitoring of the elements of the
% covariance matrix in each step
S2cov=[(init1:n)' zeros(n-init1+1,0.5*v*(v+1))];

% Initialize matrix which contains the monitoring of the determinant and of the
% trace of the covariance matrix in each step
% 2nd col = determinant
% 3rd col = trace
detS=msr;

% percn = scalar which controls up to which point of the search it is
% better to use linear indexing to extract the units forming subset. For
% example percn=0.85*n means that units belonging to susbet are found using
% linear indexing up to step m=0.85*n. After m=0.85*n units belonging to
% subset are found using a n-by-1 logical vector
percn=.85*n;
% nrepmin = scalar which controls the maximum number of repeated minima
% which must be taken in order to find new subset
nrepmin=10;


if n<32768
    bsb=int16(bsb);
    seq=int16((1:n)');
    minMDindex=zeros(1,1,'int16');
    
    % unitadd and bsbradd will be vectors which contain the k units which
    % are added using k repeated minima. More precisely the units which
    % were not in the previous subset are included in vector unitadd, those
    % which were in the previous subset are included in vector bsbradd
    unitadd=zeros(nrepmin,1,'int16');
    bsbradd=unitadd;
    
    
else
    bsb=int32(bsb);
    seq=int32((1:n)');
    minMDindex=zeros(1,1,'int32');
    
    % unitadd and bsbradd will be vectors which contain the k units which
    % are added using k repeated minima. More precisely the units which
    % were not in the previous subset are included in vector unitadd, those
    % which were in the previous subset are included in vector bsbradd
    unitadd=zeros(nrepmin,1,'int32');
    bsbradd=unitadd;
    
end

zeron1=false(n,1);

% Initialization of the n x 1 Boolean vector which contains a true in
% correspondence of the units belonging to subset in each step
bsbT=zeron1;
bsbT(bsb)=true;

% unit is the vector which will contain the units which enter subset at each
% step. It is initialized as a vector of zeros
unit=zeros(ini0,1,'int16');
lunit=length(unit);


if (rank(Y(bsb,:))<v)
    warning('FSDA:FSMeda:NoFullRank','The supplied initial subset is not full rank matrix');
    disp('FSMeda:message: FS loop will not be performed');
    out=struct;
else
    
    for mm = ini0:n
        
        % Extract units forming subset
        if mm<=percn
            Yb=Y(bsb,:);
        else
            Yb=Y(bsbT,:);
        end
        
        
        % Find vector of means inside subset
        % Note that ym is a row vector
        ym=sum(Yb,1)/mm;
        
        
        % Ym = n-by-v matrix containing deviations from the means computed
        % using units forming subset
        % Ym=Y-one*ym;
        if verLessThan2016b
            Ym = bsxfun(@minus,Y, ym);
        else
            Ym = Y - ym;
        end
        
        
        
        if mm-lunit>v+1
            
            % Find new S
            if lunit>1
                % S0=S;
                % Find units which left subset
                % Inefficient code is
                % unitout=setdiff(oldbsb,bsb);
                
                % unitoutT = Boolean for units which left subset
                % ~oldbsbF = units which were in previous subset
                % ~bsbT = units which are not in the current subset
                % unitoutT=~oldbsbF & ~bsbT;
                % Given that \not A intersect \not B = \not (A U B)
                
                if mm>percn || rankgap>nrepmin
                    % unitoutT=~(~oldbsbT | bsbT);
                    unitoutT = oldbsbT & ~bsbT;
                    unitout=seq(unitoutT);
                end
                
                lunitout=length(unitout);
                mi=sum(Y(unitout,:),1)/lunitout;
                
                % bsbr units which remained in subset
                % old inefficient code
                % bsbr=setdiff(oldbsb,unitout);
                
                % If mm>percn or if rankgap is greater than nrepmin, the units
                % which remained in subset are found using Boolean
                % operations
                % else they were immediately stored when repeated minima
                % were taken
                if mm>percn || rankgap>nrepmin
                    % oldbsbT = units which were in previous subset
                    % bsbT = units which are in current subset
                    bsbr=oldbsbT & bsbT;
                end
                mibsbr=sum(Y(bsbr,:),1)/(mm-1-lunitout);
                
                zi=sqrt(lunitout*(mm-1-lunitout)/(mm-1))*(mi-mibsbr);
                Szi=S*zi';
                % S=S+(S*(zi')*zi*S)/(1-zi*S*(zi'));
                S=S+Szi*(Szi')/(1-zi*Szi);
                if lunitout>1
                    for i=1:lunitout
                        zi=Y(unitout(i),:)-mi;
                        Szi=S*zi';
                        % S=S+(S*(zi')*zi*S)/(1-zi*S*(zi'));
                        S=S+Szi*(Szi')/(1-zi*Szi);
                    end
                end
            else
                lunitout=0;
                mibsbr=meoldbsb;
            end
            
            % mi = mean of units entering subset
            mi=sum(Y(unit,:),1)/lunit;
            % zi=sqrt(kin*(mm-1-k)/(mm-1-k+kin))*(mi-mean(Y(bsbr,:),1));
            zi=sqrt(lunit*(mm-1-lunitout)/(mm-1-lunitout+lunit))*(mi-mibsbr);
            Szi=S*zi';
            % S=S+(S*(zi')*zi*S)/(1-zi*S*(zi'));
            S=S-Szi*(Szi')/(1+zi*Szi);
            if lunit>1
                %mi=mean(Y(unit,:),1);
                for i=1:lunit
                    zi=Y(unit(i),:)-mi;
                    Szi=S*zi';
                    % S=S-(S*(zi')*zi*S)/(1+zi*S*(zi'));
                    S=S-Szi*(Szi')/(1+zi*Szi);
                end
            end
            
            % Compute Mahalanobis distance using updating formulae
            % Note that up for n>30000 it seems faster to use bsxfun rather
            % than .*
            if n<30000
                MD=(mm-1)*sum((Ym*S).*Ym,2);
            else
                MD=(mm-1)*sum(bsxfun(@times,mtimes(Ym,S),Ym),2);
            end
            
            
        else % In the initial step of the search the inverse is computed directly
            if mm>percn
                S=inv(Ym(bsbT,:)'*Ym(bsbT,:));
                [~,R]=qr(Ym(bsbT,:),0);
            else
                S=inv(Ym(bsb,:)'*Ym(bsb,:));
                [~,R]=qr(Ym(bsb,:),0);
            end
            if sum(isinf(S(:)))>0
                warning('FSDA:FSMmmd:NoFullRank',['Subset at step mm= ' num2str(mm) ' is not full rank matrix']);
                disp('FS loop will not be performed')
                out = nan;
                return
            end
            
            u=(Ym/R);
            % Compute squared Mahalanobis distances
            MD=(mm-1)*sum(u.^2,2);
        end
        
        MD=sqrt(MD);
        
        if mm>=init1
            
            if mm<=percn
                Ymbsb=Ym(bsb,:);
            else
                Ymbsb=Ym(bsbT,:);
            end
            % covYb=cov(Yb);
            covYb=Ymbsb'*Ymbsb/(mm-1);
            detcovYb=det(covYb);
            
            if scaled==1
                MD=MD*(detcovYb^(1/(2*v)));
            end
            
            if mm<=percn
                BB(bsb,mm-init1+1)=bsb;
            else
                BB(bsbT,mm-init1+1)=seq(bsbT);
            end
            % Store the means
            loc(mm-init1+1,2:end)=ym;
            
            % Store the trace and the determinant
            detS(mm-init1+1,2:end)=[detcovYb sum(diag(covYb))];
            
            % Store the elements of the covariance matrix
            aco=triu(covYb);
            aco=aco(abs(aco(:))>1e-15);
            S2cov(mm-init1+1,2:end)=aco';
            
            % Store MD inside matrix MAL
            MAL(:,mm-init1+1)=MD;
        end
        
        %         if verLessThan2017a
        MDsor=sort(MD);
        %         else
        %             MDsor=sort(MD,'ComparisonMethod','real');
        %         end
        
        if mm>=init1
            if mm<=percn
                % Store max and mmth ordered MD
                msr(mm-init1+1,2:3)= [max(MD(bsb)) MDsor(mm)];
            else
                % Store max and mmth ordered MD
                msr(mm-init1+1,2:3)= [max(MD(bsbT)) MDsor(mm)];
            end
        end
        
        if mm<n
            
            % MDmod contains modified Mahalanobis distances. The
            % Mahalanobis distance of the units belonging to subset are set
            % to inf because we need to consider the minimum of the units
            % outside subset
            MDmod=MD;
            
            if mm>percn
                MDmod(bsbT)=Inf;
            else
                MDmod(bsb)=Inf;
            end
            
            
            % oldbsbF=bsbF;
            oldbsb=bsb;
            oldbsbT=bsbT;
            % Take minimum distance of the units not belonging to subset
            [minMD,minMDindex(:)]=min(MDmod);
            
            % MDltminT = n x 1 Boolean vector which is true if corresponding MD is
            % smaller or equal minMD
            MDltminT=MD<=minMD;
            
            % MDltminbsb = n x 1 Boolean vector (if m>percn) or
            % int32 vector containing the units which certainly remain inside subset
            % i.e. those which have a true in MDltminT and belong to previous subset
            if mm>percn
                MDltminbsb=MDltminT & oldbsbT;
            else
                MDltminbsb=MDltminT(oldbsb);
            end
            
            
            % Find number of units of old subset which have a MD <= minMD
            mmtry=sum(MDltminbsb);
            
            % rankgap is the difference between m+1 (size of new size) and
            % the number of units of old subset which have a distance <=
            % minMD. For example if rankgap is =3, three more units must be
            % added.
            rankgap=mm+1-mmtry;
            
            if rankgap==1
                % Just one new unit entered subset
                unit= minMDindex;
                
                % Compute new bsbT and new bsb
                if mm<=percn
                    % new bsb is equal to oldbsb plus unit which just entered
                    bsb=[oldbsb;unit];
                end
                % bsbT is equal to old bsbT after adding a single true in
                % correspondence of the unit which entered subset
                bsbT(minMDindex)=true;
                
            elseif rankgap>1 && rankgap <=nrepmin
                
                % MDmod is the vector of Mahalanobis distance which will have
                % a Inf in correspondence of the units of old subset which
                % had a MD smaller than minMD
                MDmod=MD;
                
                
                % Find bsbrini, i.e. the vector which will contain the
                % units which remain in the subset in the next step
                % Note that bsbrini is defined using Boolean vector bsbT
                % when mm is greater than percn otherwise it uses numerical
                % vector bsb
                if mm<=percn
                    % bsbrini = vector containing the list of the units
                    % which certainly remain inside subset (i.e. those
                    % which have a MD smaller than minMD). In order to find
                    % bsbr we must check whether the k units which will be
                    % included were or not in the previous subset
                    bsbini=MDltminT(bsb);
                    bsbrini=bsb(bsbini);
                    % unitout = list of the units which potentially left
                    % subset. We say potentially because there are still k
                    % units to be included
                    % unitout=bsb(~bsbini);
                    unitout=bsb(~bsbini);
                    MDmod(bsbrini)=Inf;
                    
                else
                    % bsbriniT = Boolean vector which is true if the
                    % corresponding unit belonged to previous subset and
                    % has a MD smaller than minMD. This vector is nothing
                    % but the Boolean version of bsbrini.
                    % bsbriniT=MDltminT & bsbT;
                    bsbriniT= MDltminbsb;
                    MDmod(bsbriniT)=Inf;
                    
                end
                
                
                kk=1; zz=1;
                % In the following loop we add k units to form the new
                % subset of m+1 units Note that if the difference between
                % m+1 and the rank of the min outside subset is equal to rankgap,
                % than at most rankgap minima must be calculated to find
                % the the (m+1)-th order statistic
                for jj=1:rankgap
                    
                    [~,minMDindex(:)]=min(MDmod);
                    % minMDindex = index of the unit which is about to
                    % enter subset. We check whether unit minMDindex
                    % belonged or not to previous subset If unit minMDindex
                    % belonged to previous subset than a true is added into
                    % vector bsbriniT and the unit is included in vector
                    % bsbradd If unit minMDindex did not belong to previous
                    % subset, than minMDindex is included in vector unitadd
                    if bsbT(minMDindex)
                        if mm<=percn
                            bsbradd(zz)=minMDindex;
                            zz=zz+1;
                            % Delete from vector unitout (containing the
                            % list of the units which went out of the
                            % subset) element minMDindex
                            unitout=unitout(unitout~=minMDindex);
                        else
                            bsbriniT(minMDindex)=true;
                        end
                    else
                        unitadd(kk)=minMDindex;
                        kk=kk+1;
                    end
                    % disp(posunit(posncl1))
                    MDmod(minMDindex)=Inf;
                end
                
                % unit = vector containing all units which enter the new subset
                % but did not belong to previous subset
                unit=unitadd(1:kk-1);
                % bsbr = vector containing all units which enter the new
                % subset and were also in the previous subset
                % bsb = units forming new subset.
                if mm<=percn
                    
                    bsbr=[bsbrini;bsbradd(1:zz-1)];
                    bsb=[bsbr;unit];
                else
                    % After the instruction which follows bsbriniT
                    % will be exactly equal to bsbT
                    % Note that bsbT has been computed through the 3 following instructions
                    % -----------    bsbriniT=MDltminT & bsbT;
                    % -----------    bsbriniT(minMDindex)=true;
                    % -----------    bsbriniT(unit)=true;
                    bsbriniT(unit)=true;
                end
                
                
                
                % Compute bsbT (Boolean vector which identifies new subset)
                if mm<=percn
                    bsbT=zeron1;
                    bsbT(bsb)=true;
                else
                    bsbT=bsbriniT;
                end
                
            else %  rankgap>nrepmin
                
                % Old sorting
                %                 [~,MDsor]=sort(MD);
                %                 bsb=MDsor(1:mm+1);
                %                 bsbT=zeron1;
                %                 bsbT(bsb)=true;
                
                
                % New sorting based on quickselectFS
                [ksor]=quickselectFS(MD,mm+1,minMDindex);
                bsbT=MD<=ksor;
                
                if sum(bsbT)==mm+1
                    if mm<=percn
                        bsb=seq(bsbT);
                    end
                else
                    bsbmin=seq(MD<ksor);
                    bsbeq=seq(MD==ksor);
                    
                    bsb=[bsbmin;bsbeq(1:mm+1-length(bsbmin))];
                    
                    bsbT=zeron1;
                    bsbT(bsb)=true;
                end
                
                % unit = vector containing units which just entered subset;
                unit=find(bsbT & ~oldbsbT);
                
            end
            
            
            
            if mm>=init1
                
                % mmd contains minimum of Mahalanobis distances among
                % the units which are not in the subset at step m
                mmd(mm-init1+1,2:3)=[minMD MDsor(mm+1)];
                
                % store gap
                gap(mm-init1+1,2:3)=mmd(mm-init1+1,2:3)-msr(mm-init1+1,2:3);
                
            end
            
            % store mean of units forming old subset
            meoldbsb=ym;
            
            lunit=length(unit);
            
            if (mm>=init1)
                if (lunit<=10)
                    Un(mm-init1+1,2:(lunit+1))=unit;
                else
                    if msg==1
                        disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        disp(['Number of units which entered=' int2str(lunit)]);
                    end
                    Un(mm-init1+1,2:end)=unit(1:10);
                end
            end
        else
        end
    end
    
    
    
    if scaled==1
        % Scale distances with ratio of determinants
        detcovY=(detcovYb^(1/(2*v)));
        % Rescale MD
        MAL=MAL/detcovY;
        % Rescale minimum MD
        mmd(:,2:3)=mmd(:,2:3)/detcovY;
        % Rescale max MD
        msr(:,2:3)=msr(:,2:3)/detcovY;
    elseif scaled==2
        % Scale distances with asymptotic consistency factor
        mm=init1:n;
        a=chi2inv(mm/n,v);
        corr=sqrt((n./mm).*(chi2cdf(a,v+2)));
        MAL=bsxfun(@times,MAL,corr);
        mmd(:,2:3)=bsxfun(@times,mmd(:,2:3),corr(2:end)');
        % Rescale max MD
        msr(:,2:3)=bsxfun(@times,msr(:,2:3),corr');
        
    end
    
    
    % Plot minimum Mahalanobis distance with 1%, 50% and 99% envelopes
    if options.plots==1
        figure;
        quant=[0.01;0.5;0.99];
        % Compute theoretical envelops for minimum Mahalanobis distance based on all
        % the observations for the above quantiles.
        [gmin] = FSMenvmmd(n,v,'prob',quant,'init',init1,'scaled',scaled);
        plot(mmd(:,1),mmd(:,2),'tag','data_mmd');
        
        % include specified tag in the current plot
        set(gcf,'tag','pl_mmd');
        set(gcf,'Name', 'Monitoring of Minimum Mahalnobis distance', 'NumberTitle', 'off');
        
        % Superimpose 1%, 99%, 99.9% envelopes based on all the observations
        lwdenv=2;
        % Superimpose 50% envelope
        line(gmin(:,1),gmin(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
        % Superimpose 1% and 99% envelope
        line(gmin(:,1),gmin(:,2),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        line(gmin(:,1),gmin(:,4),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        
        xlabel('Subset size m');
        if scaled==1
            ylabel('Monitoring of scaled minimum Mahalanobis distance');
        else
            ylabel('Monitoring of minimum Mahalanobis distance');
        end
        
    end
    
    
    
    % Divide each column of detS by the final value at the end of the search
    detS(:,2)=detS(:,2)/detS(end,2);
    detS(:,3)=detS(:,3)/detS(end,3);
    
    out.MAL=MAL;
    out.BB=BB;
    out.msr=msr;
    out.mmd=mmd;
    out.gap=gap;
    out.S2cov=S2cov;
    out.detS=detS;
    out.Un=Un;
    out.Y=Y;
    out.Loc=loc;
    out.class='FSMeda';
end

end
%FScategory:MULT-Multivariate
