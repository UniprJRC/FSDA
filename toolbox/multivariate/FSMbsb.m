function [Un,BB] = FSMbsb(Y,bsb,varargin)
%FSMbsb gives the units belonging to subset at step(s) msel of the forward search
%
%<a href="matlab: docsearchFS('FSMbsb')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Input data. Matrix.
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
% bsb :         Units forming subset. Vector. List of units forming the initial subset.
%               If bsb=0 (default) then the procedure starts with p units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb)
%                Data Types - single|double
%
%
% Optional input arguments:
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
% init :       Point where to start monitoring required diagnostics. Scalar.
%              Note that if bsb is supplied, init>=length(bsb). If init is not
%              specified it will be set equal to floor(n*0.6).
%                 Example - 'init',50
%                 Data Types - double
%
%  msg  :   It controls whether to display or not messages
%               about great interchange on the screen. Boolean.
%               If msg==true (default) messages are displyed on the screen
%               else no message is displayed on the screen
%                 Example - 'msg',false
%                 Data Types - logical
%
% nocheck :   It controls whether to perform checks on matrix Y. Scalar.
%             If nocheck is equal to 1 no check is performed on matrix Y.
%             As default nocheck=0.
%                 Example - 'nocheck',1
%                 Data Types - double
%
% plots :     Plot on the screen. Scalar.
%               If plots=1, a plot of the monitoring of minMD among
%               the units not belonging to the subset is produced on the
%               screen with 1 per cent, 50 per cent and 99 per cent confidence bands
%               else (default) no plot is produced.
%               Example - 'plots',0
%               Data Types - double
%
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
% See also FSMeda, FSM.m, FSMmmd, FSRbsb, FSRHbsb, FSRBbsb
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
%<a href="matlab: docsearchFS('FSMbsb')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

%{
    % FSMbsb with all default options.
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
    Un=FSMbsb(Ycont,0);
%}

%{
    % FSMbsb with optional argument.
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
    [~,BBsel]=FSMbsb(Ycont,0,'bsbsteps',msel);
    disp(['Units outside subset at step m=' num2str(msel)])
    disp(setdiff(1:n,BBsel))
%}

%{
    %% Monitoring the units belonging to subset in each step.
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
    [~,BBsel]=FSMbsb(Ycont,0,'bsbsteps',msel);
    seq=1:n;
    disp(['Units outside subset at step m=' num2str(msel)])
    disp(setdiff(seq,BBsel))
%}

%{
    % Specifying the point where to start monitoring.
    % Specifying the point where to start monitoring units belongng to subset.
    n=200;
    v=3;
    randn('state',123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    seqcont=11:10:51;
    disp('Contaminated units')
    disp(seqcont)
    Ycont(seqcont,:)=Ycont(seqcont,:)+2.5;
    % Analyse the units forming subset in step msel=195
    msel=195;
    [Un,BBsel]=FSMbsb(Ycont,0,'bsbsteps',msel,'init',100);
    seq=1:n;
    disp(['Units outside subset at step m=' num2str(msel)])
    disp(setdiff(seq,BBsel))
%}


%% Beginning of code

% Input parameters checking
%chkinputM does not do any check if option nocheck=1
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);

%% Input parameters checking
[n,v]=size(Y);

initdef=floor(n*0.6);

bsbstepdef='';


if nargin<2
    error('FSDA:FSMbsb:missingInputs','Initial subset is missing')
end

if coder.target('MATLAB')
    options=struct('init',initdef,'plots',0,'msg',true,'nocheck',0,'bsbsteps',bsbstepdef);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSMbsb:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
end

if nargin>2
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% And check if the optional user parameters are reasonable.
if bsb==0
    Ra=true; nwhile=1;
    while and(Ra,nwhile<100)
        % Extract a random sample of size v+1
        bsb=randsample(n,v+1);
        % Check if the var-cov matrix of the random sample is full (i.e =v)
        Ra=(rank(cov(Y(bsb,:)))<v);
        nwhile=nwhile+1;
    end
    if nwhile==100
        if coder.target('MATLAB')
            warning('FSDA:FSMbsb:NoFullRank','Unable to randomly sample full rank matrix');
        else
            disp('FSDA:FSMbsb:NoFullRank','Unable to randomly sample full rank matrix');
        end
    end
else
end

% percn = scalar which controls up to which point of the search it is
% better to use linear indexing to extract the units forming subset. For
% example percn=0.85*n means that units belonging to susbet are found using
% linear indexing up to step m=0.85*n. After m=0.85*n units belonging to
% subset are found using a n-by-1 logical vector
percn=.85*n;
% nrepmin = scalar which controls the maximum number of repeated minima
% which must be taken in order to find new subset
nrepmin=10;

if coder.target('MATLAB')
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
else
    seq=(1:n)';
    minMDindex=zeros(1,1);
    unitadd=zeros(nrepmin,1);
    bsbradd=unitadd;
end
zeron1=false(n,1);
bsb=bsb(:);

% Initialization of the n x 1 Boolean vector which contains a true in
% correspondence of the units belonging to subset in each step
bsbT=zeron1;
bsbT(bsb)=true;


% Initialization for Matlab coder
if ~coder.target('MATLAB')
    rankgap=0;
    S=zeros(v,v);
    meoldbsb=zeros(1,v);
    oldbsbT=bsbT;
    bsbr=zeros(n,1);
    unitout=bsbr;
    bsbriniT=bsbT;
    bsbrini=bsbr;
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
    if coder.target('MATLAB')
        mess=sprintf(['Attention : init1 should be >= length of supplied subset. \n',...
            'It is set equal to ' num2str(length(bsb)) ]);
        fprintf('%s\n',mess);
    else
        disp('Attention : init1 should be >= length of supplied subset. It is set equal to size of supplied subset')
    end
    init1=ini0;
elseif init1>=n
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    fprintf('%s\n',mess);
    init1=n-1;
end

bsbsteps=options.bsbsteps;
% Matrix BB will contain the units forming subset in each step (or in
% selected steps) of the forward search. The first column contains
% information about units forming subset at step init1.
if isempty(bsbsteps)
    % Default for vector bsbsteps which indicates for which steps of the fwd
    % search units forming subset have to be saved
    if n<=5000
        bsbsteps = (init1:1:n)';
    else
        bsbsteps = [init1 init1+100-mod(init1,100):100:100*floor(n/100)]';
    end
    if coder.target('MATLAB')
        BB = NaN(n,length(bsbsteps),'single');
    else
        BB = NaN(n,length(bsbsteps));
    end
elseif bsbsteps==0
    bsbsteps=(init1:n)';
    if coder.target('MATLAB')
        BB = NaN(n,n-init1+1,'single');
    else
        BB = NaN(n,n-init1+1);
    end
else
    if min(bsbsteps)<init1
        if coder.target('MATLAB')
            warning('FSDA:FSMbsb:WrongInit','It is impossible to monitor the subset for values smaller than init');
        else
            disp('FSDA:FSMbsb:WrongInit','It is impossible to monitor the subset for values smaller than init');
        end
    end
    bsbsteps=bsbsteps(bsbsteps>=init1);

    if coder.target('MATLAB')
        BB = NaN(n,length(bsbsteps),'single');
    else
        BB = NaN(n,length(bsbsteps));
    end
end

%  Un is a Matrix whose 2nd column:11th col contain the unit(s) just
%  included.
Un = cat(2 , (init1+1:n)' , NaN(n-init1,10));

% unit is the vector which will contain the units which enter subset at each
% step. It is initialized as a vector of zeros
if coder.target('MATLAB')
    unit=zeros(ini0,1,'int16');
else
    unit=zeros(ini0,1);
end
lunit=length(unit);

% If the subset Y(bsb,:) is not full rank or a column is constant, then we
% return as output an empty structure.
if (rank(Y(bsb,:))<v) || min(max(Y(bsb,:)) - min(Y(bsb,:))) == 0
    if coder.target('MATLAB')
        warning('FSDA:FSMbsb:NoFullRank','The supplied initial subset does not produce a full rank matrix');
    else
        disp('FSDA:FSMmmd:NoFullRank','The supplied initial subset does not produce a full rank matrix');
    end
    disp('FS loop will not be performed')
    Un=NaN;
    BB=NaN;
else
    % ij = index which is linked with the columns of matrix BB. During the
    % search every time a subset is stored inside matrix BB ij icreases by one
    ij=1;

    for mm = ini0:n

        % Extract units forming subset
        if mm<=percn
            Yb=Y(bsb,:);
        else
            Yb=Y(bsbT,:);
        end

        % If required, store units forming subset at each step
        if (mm>=init1)

            if intersect(mm,bsbsteps(:))==mm
                if mm<=percn
                    BB(bsb,ij)=bsb;
                else
                    BB(bsbT,ij)=seq(bsbT);
                end
                ij=ij+1;
            end
        end

        % Find vector of means inside subset
        % Note that ym is a row vector
        ym=sum(Yb,1)/mm;



        % Ym = n-by-v matrix containing deviations from the means computed
        % using units forming subset
        % Ym=Y-one*ym;
        Ym = bsxfun(@minus,Y, ym);

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
                    bsbrT=oldbsbT & bsbT;
                    mibsbr=sum(Y(bsbrT,:),1)/(mm-1-lunitout);
                else
                    mibsbr=sum(Y(bsbr,:),1)/(mm-1-lunitout);
                end

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
                if coder.target('MATLAB')
                    warning('FSDA:FSMbsb:NoFullRank',['Subset at step mm= ' num2str(mm) ' is not full rank matrix']);
                else
                    disp('FSDA:FSMbsb:NoFullRank','Subset at step mm is not full rank matrix');
                end
                disp('FS loop will not be performed')

                Un=NaN;
                BB=NaN;
                return
            end

            u=(Ym/R);
            % Compute squared Mahalanobis distances
            MD=(mm-1)*sum(u.^2,2);
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
end

if coder.target('MATLAB')
    plots=options.plots;
    if plots==1
        % Create the 'monitoring units plot'
        figure;
        plot(bsbsteps,BB','bx')
        xlabel('Subset size m');
        ylabel('Monitoring units plot');
    end
end
end
%FScategory:MULT-Multivariate
