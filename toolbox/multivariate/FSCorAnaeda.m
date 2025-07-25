function out = FSCorAnaeda(N,varargin)
%FSCorAnaeda performs forward search in correspondence analysis with exploratory data analysis purposes
%
%
%<a href="matlab: docsearchFS('FSCorAnaeda')">Link to the help function</a>
%
% Required input arguments:
%
%  N :  contingency table or structure.
%               Array or table of size I-by-J or structure. If N is a
%               structure it contains the following fields:
%               N.N = contingency table in array format of size I-by-J.
%               N.Ntable = contingency table in table format of size  I-by-J.
%               N.loc = initial location estimate for the matrix of Profile
%                   rows of the contingency table (row vector or length J).
%               N.weights= I x 1 vector containing the proportion of the
%                    mass of each rows of matrix N in the computation of
%                    the MCD estimate of location. If N.weigths(2)=0.1 it
%                    means that row 2 of the contingency table contributes
%                    with 10 per cent of its mass. The initial subset is
%                    based on N.weights.
%               N.NsimStore= array of size I-by-J times nsimul containing
%                   in each column the nsimul simulated contingency tables.
%               Note that input structure N can be conveniently created by
%               function mcdCorAna.
%               If N is not a struct it is possible to specify the rows
%               of the contingency table forming initial subset with input
%               option bsb.
%                Data Types - single|double
%
% Optional input arguments:
%
% bsb  :       Initial subset. Vector of positive integers containing the
%              indexes of the rows of the contingency table which have to
%              be used to initialize the forward search. If bsb is empty
%              and required input argument is a struct N.loc will be used.
%              If bsb is supplied and N is a struct N.loc is ignored.
%              The default value of bsb is empty, and if N is not a struct
%              a random subset containing round(n/5) units will be used.
%                 Example - 'bsb',[3 6 8 10 12 14]
%                 Data Types - double
%
%   conflev :  simultaneous confidence interval to declare units as
%              outliers. Scalar.
%              The default value of conflev is 0.99, that is a 99 per
%              cent simultaneous confidence level.
%              Confidence level are based on simulated contingency
%              tables. This input argument is ignored if optional input
%              argument mmdEnv is not missing
%              Example - 'conflev',0.99
%              Data Types - numeric
%
% init :       Point where to start monitoring required diagnostics. Scalar.
%              Note that if init is not
%              specified it will be set equal to floor(n*0.6).
%              where the total number of units in the contingency table.
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
% resc :        Rescale or not the envelopes. Boolean. It controls whether to rescale or not the envelopes of min
%               MD when if in the initial part of the search  is steadily
%               above or below the 5 and 95 per cent confidence bands. The
%               default value of resc is true.
%                 Example - 'resc',false
%                 Data Types - logical
%
%  msg  :       It controls whether to display or not messages
%               about great interchange on the screen. Scalar.
%               If msg==1 (default) messages are displayed on the screen
%               else no message is displayed on the screen.
%                 Example - 'msg',0
%                 Data Types - double
%
%
% Output:
%
%         out:   structure which contains the following fields
%
%   out.MAL=    I x (n-init+1) = matrix containing the monitoring of
%               Mahalanobis distances.
%               1st row = distance for first row;
%               ...;
%               Ith row = distance for Ith row.
%    out.BB=    I-by-(n-init+1) matrix containing the information about the units belonging
%               to the subset at each step of the forward search.
%               1st col = indexes of the units forming subset in the
%               initial step;
%               ...;
%               last column = units forming subset in the final step (all
%               units).
%               Note that the numbers inside out.BB vary in the interval [0
%               1] and represent the proportions in which each unit is
%               represented in the subset. 0.1 means that the associated
%               row is represented in the subset with 10 per cent of its
%               mass.
%   out.mmd=    n-init-by-2 matrix which contains the monitoring of minimum
%               MD or (m+1)th ordered MD  at each step of
%               the forward search.
%               1st col = fwd search index (from init to n-1);
%               2nd col = minimum MD;
%   out.ine=    n-init-by-2 matrix which contains the monitoring of inertia
%               at each step of the forward search.
%               1st col = fwd search index (from init to n);
%               2nd col = inertia;
%   out.Loc=     (n-init+1)-by-J matrix containing the monitoring of
%               estimated  means for each variable in each step of
%               the forward search.
%    out.Un=        (n-init) x 11 Matrix which contains the unit(s)
%               included in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,2) for example contains
%               the unit included in step init+1 Un(end,2) contains the
%               units included in the final step of the search
%     out.N=    Original contingency table, in array format.
%     out.Ntable = Original contingency table in table format (if initially
%               supplied).
%
%     out.Y = array of size I-by-J containing row profiles.
% out.class=    'FSCorAnaeda'
%
% See also mcdCorAna.m, FSMeda.m, malfwdplot.m
%
% References:
%
% Atkinson, A.C., Riani, M. and Cerioli, A. (2004), "Exploring multivariate
% data with the forward search", Springer Verlag, New York.
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSCorAnaeda')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    % FSCorAnaeda with all default options.
    % Generate contingency table of size 50-by-5 with total sum of n_ij=2000.
    I=50;
    J=5;
    n=2000;
    % nrowt = column vector containing row marginal totals
    nrowt=(n/I)*ones(I,1);
    % ncolt = row vector containing column marginal totals
    ncolt=(n/J)*ones(1,J);
    out1=rcontFS(I,J,nrowt,ncolt);
    N=out1.m144;
    RAW=mcdCorAna(N,'plots',0);
    ini=round(sum(sum(RAW.N))/4);
    out=FSCorAnaeda(RAW);
%}

%{
    %% FSCorAnaeda with optional arguments.
    % Generate contingency table of size 50-by-5 with total sum of n_ij=2000.
    I=50;
    J=5;
    n=2000;
    % nrowt = column vector containing row marginal totals
    nrowt=(n/I)*ones(I,1);
    % ncolt = row vector containing column marginal totals
    ncolt=(n/J)*ones(1,J);
    out1=rcontFS(I,J,nrowt,ncolt);
    N=out1.m144;
    RAW=mcdCorAna(N,'plots',0);
    ini=round(sum(sum(RAW.N))/4);
    out=FSCorAnaeda(RAW,'plots',1);
%}

%{
    % FSCorAnaeda starting from a random initial subset.
    % Generate contingency table of size 50-by-5 with total sum of n_ij=2000.
    I=50;
    J=5;
    n=2000;
    % nrowt = column vector containing row marginal totals
    nrowt=(n/I)*ones(I,1);
    % ncolt = row vector containing column marginal totals
    ncolt=(n/J)*ones(1,J);
    out1=rcontFS(I,J,nrowt,ncolt);
    N=out1.m144;
    % The first input argument is a contingency table and no initial subset
    % and no initial location is supplied
    out=FSCorAnaeda(N,'plots',1);
%}

%% Beginning of code

% Initialization of bsb as an empty vector.
bsb=[];
label=[];

if isstruct(N)

    RAW=N;
    loc=N.loc;

    if isfield(N,'Ntable')
        Ntable=N.Ntable;
        if istable(N.Ntable)
            label=string(N.Ntable.Properties.RowNames);
        elseif istimetable(N.Ntable)
            label=string(N.Ntable.Properties.RowTimes);
        else
        end
    else
        Ntable=N.N;
    end



    % Find the contingency table associated to N.loc
    % weights is the vector which tells us how each row of the contingency
    % table is represented inside subset.
    weights=N.weights;
    N=N.N;
    % Note that some rows of Niter are equal to 0 and must be deleted.
    Niter=N.*weights;
    seqI=1:size(N,1);
    % bsb contains the rows of N which had weights strictly greater than 0
    bsbini=seqI(weights>0);

    rowstodel=sum(Niter,2)==0;
    Niter(rowstodel,:)=[];
    Nisstruct=true;

else
    Nisstruct=false;
    RAW=N;
    if istable(N)
        Ntable=N;
        N=table2array(N);
    else
        Ntable=array2table(N);
    end
end
% n= total sample size
n=round(sum(N,'all'));

[I,J]=size(N);

% P = correspondence matrix  containing relative frequencies
P = (1/n) * N;

% r= vector which contains row masses = centroids of the column profiles
r  = sum(P,2) ;

rtimesn=r*n;

% Row profiles (equation 4.14)
% Dr=diag(r);
% ProfilesRows = Dr^(-1) * P;
% Every row of matrix P is divided by the row total
ProfileRows=P./r;


init1=floor(n*0.5);
plots=0;
% Simultaneous confidence envelope to declare the outliers
conflev=0.99;
resc=true;

if nargin > 1

    options=struct('init',init1,'plots',plots,'msg',1,'bsb',[],'conflev',conflev,'resc',resc);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSCorAnaeda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        aux.chkoptions(options,UserOptions)
    end

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    % check init
    init1=options.init;
    msg=options.msg;
    plots=options.plots;
    bsb=options.bsb;
    resc=options.resc;
end

if isempty(bsb)
    ini0=round(n/5);
    if Nisstruct == true
        % In this case initial estimate of location is supplied in input structure RAW.loc;
        % Find vector of means inside subset
        ym=loc;
        bsb=bsbini;
    else
        % A random subset of ini0 units will be extracted.
        indsamp=randsample(I,I);

        cumsumnjdot=cumsum(rtimesn(indsamp));
        indexesCR=find(cumsumnjdot<ini0,1,'last');



        if isempty(indexesCR)
            indexesCR=0;
            unitstoADD=ini0;
        else
            % Find how many units must be included from
            % row indexesCR+1 of the initial contingency table;
            unitstoADD=ini0-cumsumnjdot(indexesCR);
        end

        bsb=indsamp(1:indexesCR+1);

        Niter=N(bsb,:);
        Niter(end,:)=unitstoADD*Niter(end,:)/rtimesn(indsamp(indexesCR+1));

        % Compute centroid based on the random extracted ini0 units
        ym  = sum(Niter,1)/ini0;
    end
else
    % Elements of N(bsb,:) are used to find initial estimate of location
    Nini=N(bsb,:);
    Niter=N(bsb,:);
    % ini0 is the sum of elements in the contingency table Nini
    ini0=sum(Nini,'all');
    ym=sum(Nini,1)/ini0;
end


if  init1 <J
    mess=sprintf(['Attention : init1 should not be smaller  than J. \n',...
        'It is set to J.']);
    fprintf('%s\n',mess);
    init1=J;
elseif init1<ini0
    mess=sprintf(['Attention : init1 should be >= length of supplied subset. \n',...
        'It is set equal to ' num2str(ini0) '=sum(N(bsb,:),''all'')']);
    fprintf('%s\n',mess);
    init1=ini0;
elseif init1>=n
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n-1.']);
    fprintf('%s\n',mess);
    init1=n-1;
end


% Initialize matrix which will contain the rows forming subset in each
% step of the fwd search
numcol=(n-init1+1);

BB=NaN(I,numcol);

% Initialize matrix which will contain the MD monitored in each
% step of the fwd search
MAL=BB;

% Initialize matrix which will contain the means of the variables monitored in
% each step
Loc=cat(2,(init1:n)',NaN(numcol,J));

%  Un is a Matrix whose 2nd column:11th col contains the row(s) just
%  included.
Un = cat(2,(init1+1:n)',NaN(n-init1,10));

%  mmd has two columns
%  1st col = dimension of the subset
%  2nd col min. MD among the units not belonging to the subset
mmd=[(init1:n-1)' zeros(n-init1,1)];

% 2nd col ine=chi2/mm statistic using the units belonging to subset
ine=[(init1:n)' zeros(n-init1+1,1)];

% funz2chi2 = anonymous function with 2 input args x and Ntheovect
% which computes the Chi2 statistic
% x is the vec operator of the contingency table in array format
% Ntheovec is the vec operator applied to the matrix of theoretical
% frequencies
% Ntheo=nrowt*ncolt/n;
% Ntheovec=Ntheo(:);
funzchi2=@(x,Ntheovec) sum(((x-Ntheovec).^2)./Ntheovec);

zeron1=false(I,1);
bsbT=zeron1;
bsbT(bsb)=true;

for mm = ini0:n

    % Mahalanobis distances from initialloc and Initialshape
    % mahaldist = r.*mahalFS(ProfileRows, initialloc, initialcov);
    MD=mahalCorAna(ProfileRows,ym);

    mahaldist=sqrt(MD);

    if mm>=init1

        nrowt=sum(Niter,2);
        % Find the proportion of units belonging to each row in the subset
        bsbProp=nrowt./sum(N(bsb,:),2);

        % Store value of Inertia (Chi2/m) at step mm
        ncolt=sum(Niter,1);
        Ntheo=nrowt*ncolt/mm;
        chi2mm=funzchi2(Niter(:),Ntheo(:));
        ine(mm-init1+1,2)=chi2mm/mm;

        % BB(bsb,mm-init1+1)=bsb;
        BB(bsb,mm-init1+1)=bsbProp;

        % Store the means
        Loc(mm-init1+1,2:end)=ym;

        % Store the trace and the determinant
        % detS(mm-init1+1,2:end)=[detcovYb sum(diag(covYb))];


        % Store weighted MD inside matrix MAL
        MAL(:,mm-init1+1)=mahaldist .*r;
    end


    if mm<n

        % oldbsb=bsb;
        oldbsbT=bsbT;
        % sort squared MD distances multiplied by row masses and by n
        [mahaldistsor,indsortdist]=sort(mahaldist.*rtimesn);


        % The rows sortdist(1:indexesCR) will be completely
        % represented;
        cumsumnjdot=cumsum(rtimesn(indsortdist))  +0.000001;
        indexesCR=find(cumsumnjdot<mm+1,1,'last');

        if isempty(indexesCR)
            indexesCR=0;
            unitstoADD=mm+1;
        else
            % Find how many units must be included from
            % row indexesCR+1;
            unitstoADD=mm+1-cumsumnjdot(indexesCR);
        end


        bsb=indsortdist(1:indexesCR+1);

        bsbT=zeron1;
        bsbT(bsb)=true;

        Niter=N(bsb,:);
        Niter(end,:)=unitstoADD*Niter(end,:)/rtimesn(indsortdist(indexesCR+1));

        % Compute new centroid based on mm+1 units
        ym  = sum(Niter,1)/(mm+1);


        % call to setdiff is much slower
        % unit=setdiff(bsb,oldbsb);
        unit=find(bsbT & ~oldbsbT);
        lunit=length(unit);


        if mm>=init1



            % mmd contains minimum of Mahalanobis distances among
            % the units which are not in the subset at step m
            % store minMD 
            IndexminMD=find(cumsumnjdot>=mm+1,1,'first');

            mmd(mm-init1+1,2)= mahaldistsor(IndexminMD);


            if (lunit<=2)
                if isempty(unit) && (mm-init1)>1
                    Un(mm-init1+1,2:end-1)=Un(mm-init1,2:end-1);

                else
                    Un(mm-init1+1,2:(lunit+1))=unit;
                    if ~isempty(unit)
                        if length(unit)>3
                            if msg==true
                                disp(['Warning: more than one unit entered in step' num2str(mm)])
                            end
                        end
                        Un(mm-init1+1,end)=unit(1);
                    end
                end


            else
                if msg==1
                    disp(['Warning: interchange greater than 1 when m=' int2str(mm)]);
                    disp(['Number of units which entered=' int2str(lunit)]);
                end
                Un(mm-init1+1,2:length(unit)+1)=unit(1:length(unit));
            end
        end
    end
end

% Plot minimum Mahalanobis distance with 1%, 50% and 99% envelopes
if plots==1
    figure;
    quant=[0.05;0.5;0.95;0.01;conflev];
    % Compute theoretical envelops for minimum Mahalanobis distance based on all
    % the observations for the above quantiles.
    disp('Creating empirical confidence band for minimum (weighted) Mahalanobis distance')
    [gmin] = FSCorAnaenv(RAW,'prob',quant,'init',init1);

    % Trajectory of mmd is rescaled it in the initial part it is outside
    % the envelopes
    if resc==true
        warmup=500;
        warmup=min([find(mmd(:,1)>round(n/2),1),warmup]);
        if sum(mmd(1:warmup,2)<gmin(1:warmup,2))>warmup/2   || sum(mmd(1:warmup,2)>gmin(1:warmup,4))>warmup/2
            coeff=mean(gmin(1:warmup,3))-mean(mmd(1:warmup,2));
            gmin(:,2:end)=gmin(:,2:end)-coeff;
        end

    end

    figure;
    plot(mmd(:,1),mmd(:,2),'tag','data_mmd');

    % include specified tag in the current plot
    set(gcf,'tag','pl_mmd');
    set(gcf,'Name', 'Monitoring of Minimum (weighted) Mahalanobis distance', 'NumberTitle', 'off');

    lwdenv=2;
    % Superimpose 50% envelope
    line(gmin(:,1),gmin(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
    % Superimpose 1% and conflev% envelope
    line(gmin(:,1),gmin(:,5),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    line(gmin(:,1),gmin(:,end),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');

    xlabel('Subset size m');
    ylabel('Monitoring of minimum (weighted) Mahalanobis distance');
    yline(max(gmin(:,end)))

    sel=~ismissing(Un(:,end));

    if isempty(label)
        text(mmd(sel,1),mmd(sel,end)*1.05,num2str(Un(sel,end)));
    else
        text(mmd(sel,1),mmd(sel,end)*1.05,label(Un(sel,end)));
    end

end

out.MAL=MAL;
out.BB=BB;
out.mmd=mmd;
out.ine=ine;
out.Un=Un;
out.N=N;
out.Ntable=Ntable;
out.Loc=Loc;
out.class='FSCorAnaeda';
out.Y=ProfileRows;

end
%FScategory:MULT-Multivariate
