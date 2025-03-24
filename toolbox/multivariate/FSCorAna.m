function out = FSCorAna(N,varargin)
%FSCorAna performs automatic outlier based on the forward search in correspondence analysis
%
%
%<a href="matlab: docsearchFS('FSCorAna')">Link to the help function</a>
%
% Required input arguments:
%
%  N :  contingency table or structure.
%               Array or table of size I-by-J or structure. If N is a
%               structure it contains the following fields:
%               N.N = contingency table in array format of size I-by-J.
%               N.Ntable = this field is not compulsory where Ntable is a table
%                   or a timetable. If this field is present the
%                   label of the rows which are used are taken from
%                   RAW.Ntable.Properties.RowTimes (in presence of a timetable)
%                   RAW.Ntable.Properties.RowNames (in presence of a table).
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
%               N.simulateUnderH0 = boolean. If it is true the simulated
%                   contingency tables  have been specified under H0.
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
%              outliers. Scalar inside (0, 1).
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
%      mmdEnv : Matrix which contains the precalculated empirical envelopes of
%               minimum Mahalanobis distance. Matrix or scalar missing
%               value (default). If this optional input argument is not missing the empirical
%               envelopes are taken from this optional argument and are not
%               calculated. First column is subset size Second column is 1
%               per cent simultaneous empirical envelope. Third column is 50
%               per cent simultaneous empirical envelope. Fourth column is
%               conflev per cent simultaneous empirical envelope which is
%               used to detect the outliers. The default value of
%               mmdStoreSim is a missing value, that is the envelopes are
%               based on the N.NsimStore pregenerated contingency tables or
%               if N.NsimStore is not present are generated assuming
%               independence between rows and columns
%               Example - 'mmdEnv',[]
%               Data Types - double
%
%    StoreSim  :  Store minimum Mahalanobis distance quantiles.
%                 Boolean. Boolean which specifies whether to store or not
%                 as fields named mmdStore the simulated envelopes of the
%                 minimum Mahalanobis distance monitored along the search.
%               Example - 'plots',0
%               Data Types - double
%
%
%  msg  :       It controls whether to display or not messages
%               about envelope creation on the screen. Logical.
%               If msg==1 (default) messages are displayed on the screen
%               else no message is displayed on the screen.
%                 Example - 'msg',false
%                 Data Types - logical
%
% resc :        Rescale or not the envelopes. Boolean. It controls whether to rescale or not the envelopes of min
%               MD when if in the initial part of the search  is steadily
%               above or below the 5 and 95 per cent confidence bands. The
%               default value of resc is true.
%                 Example - 'resc',false
%                 Data Types - logical
%
%      label : row labels. Cell or vector of strings.
%               Cell or vector of strings of length n containing the labels
%               of the rows. If input is a table or a timetable the row
%               labels are automatically taken from the row names.
%                   Example - 'label',{'UK' ...  'IT'}
%                   Data Types - cell or characters or vector of strings
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
% addRowNames  : add or not names of the rows to the plot of min MD.
%               Boolean. If this option is equal to true (default) the
%               first time a row is included in the subset is shown in the
%               plot with the corresponding row label or row number.
%               Example - 'addRowNames',false
%               Data Types - logical
%
% Output:
%
%         out:   structure which contains the following fields
%
%out.outliers=  k x 1 vector containing the list of the units declared as
%               outliers or empty value if the sample is homogeneous
%   out.mmd=    n-init-by-5 matrix which contains the monitoring of minimum
%               MD  at each step of
%               the forward search.
%               1st col = fwd search index (from init to n-1);
%               2nd col = minimum MD weighted by row mass;
%               3rd col = 1 per cent envelope;
%               4th col = 50 per cent envelope;
%               5th col = conflev per cent envelope;
%
%   out.ine=    n-init-by-2 matrix which contains the monitoring of inertia
%               at each step of the forward search.
%               1st col = fwd search index (from init to n);
%               2nd col = inertia;
%
%    out.Un=        (n-init) x 11 Matrix which contains the unit(s)
%               included in the subset at each step of the fwd search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,2) for example contains
%               the unit included in step init+1 Un(end,2) contains the
%               units included in the final step of the search
%     out.N=    Original contingency table, in array format.
% out.loc     = 1 x v  vector containing location of the data.
% out.md      = n x 1 vector containing the estimates of the robust
%               Mahalanobis distances (in squared units) mutiplied by the row masses.
%               This vector
%               contains the distances of each observation from the
%               location of the data, relative to the scatter matrix cov.
% out.thresh  = threshold for minMD with which outliers have been declared
% out.conflev = simultaneous confidence level which has been used to
%               declare the outliers.
% out.simulateUnderH0 = boolean. If it is true the simulated
%                   contingency tables have been specified under H0.
% out.nsimul       = number of simulations which have been used to create
%                   the envelopes. This information is taken from
%                   size(N.NsimStore,2)
%     out.Y = array of size I-by-J containing row profiles.
% out.class=    'FSCorAna'
%
% See also mcdCorAna.m, FSCorAnaeda.m, FSR.m
%
% References:
%
% Atkinson, A.C., Riani, M. and Cerioli, A. (2004), "Exploring multivariate
% data with the forward search", Springer Verlag, New York.
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSCorAna')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    % FSCorAna with all default options (input is output from mcdCorAna).
    % Generate contingency table of size 50-by-5 with total sum of n_ij=2000.
    I=20;
    J=5;
    n=5000;
    % nrowt = column vector containing row marginal totals
    nrowt=(n/I)*ones(I,1);
    % ncolt = row vector containing column marginal totals
    ncolt=(n/J)*ones(1,J);
    out1=rcontFS(I,J,nrowt,ncolt);
    N=out1.m144;
    RAW=mcdCorAna(N,'plots',0);
    ini=round(sum(sum(RAW.N))/4);
    out=FSCorAna(RAW);
%}

%{
    % FSCorAna with input contingency table.
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
    out=FSCorAna(N,'plots',1);
%}


%{
    %% Use pregenerated contingency tables to find envelopes for mmd.
    load clothes.mat
    % Now FSCorAna uses the pregenerated tables coming from mcdCorAna.
    % Example of findEmpiricalEnvelope a struct
    findEmp=struct;
    % Generate nsimul contingency tables
    findEmp.nsimul=1000;
    % Under the null hypothesis of independence
    findEmp.underH0=true;
    % Store the nsimul robust distance sorted (for each row)
    findEmp.StoreSim=true;
    RAW=mcdCorAna(clothes,'plots',0,'findEmpiricalEnvelope',findEmp);
    out=FSCorAna(RAW);
%}

%% Beginning of code

% Input parameters checking
% chkinputM does not do any check if option nocheck=1
% Initialization of bsb as an empty vector.
bsb=[];
label=[];
if isstruct(N)

    if isfield(N,'Ntable')
        if istable(N.Ntable)
            label=string(N.Ntable.Properties.RowNames);
        elseif istimetable(N.Ntable)
            label=string(N.Ntable.Properties.RowTimes);
        else
        end
    end


    if isfield(N,'simulateUnderH0')
        simulateUnderH0=N.simulateUnderH0;
    else
        simulateUnderH0=true;
    end


    RAW=N;
    loc=N.loc;

    % Find the contingency table associated to N.loc
    % weights is the vector which tells us how each row of the contingency
    % table is represented inside subset.
    weights=N.weights;
    N=N.N;
    seqI=1:size(N,1);
    % bsb contains the rows of N which had weights strictly greater than 0
    bsbini=seqI(weights>0);

    % Note that some rows of Niter are equal to 0 and must be deleted.
    Niter=N.*weights;
    rowstodel=sum(Niter,2)==0;
    Niter(rowstodel,:)=[];
    Nisstruct=true;
else
    Nisstruct=false;
    RAW=N;
    simulateUnderH0=true;
end

if istable(N) || istimetable(N)
    if istable(N)
        label=string(N.Properties.RowNames);
    elseif istimetable(N)
        label=string(N.Ntable.Properties.RowTimes);
    else
    end

    N=table2array(N);



end

% n= total sample size
n=round(sum(N,'all'));

[I,J]=size(N);

% P = correspondence matrix  containing relative frequencies
P = (1/n) * N;

% r= vector which contains row masses = centroids of the column profiles
r  = sum(P,2) ;

rtimesn=round(r*n);

% Dr=diag(r);
% ProfilesRows = Dr^(-1) * P;
% Every row of matrix P is divided by the row total
ProfileRows=P./r;

init1=floor(n*0.25);
plots=1;
% Simultaneous confidence envelope to declare the outliers
conflev=0.99;
StoreSim=false;
mmdEnv=[];
msg=true;

% funz2chi2 = anonymous function with 2 input args x and Ntheovect
% which computes the Chi2 statistic
% x is the vec operator of the contingency table in array format
% Ntheovec is the vec operator applied to the matrix of theoretical
% frequencies
% Ntheo=nrowt*ncolt/n;
% Ntheovec=Ntheo(:);
funzchi2=@(x,Ntheovec) sum(((x-Ntheovec).^2)./Ntheovec);

resc=true;
addRowNames=true;
if nargin > 1

    options=struct('init',init1,'plots',plots,'msg',msg,'bsb',[],'conflev',conflev, ...
        'StoreSim',StoreSim,'mmdEnv',mmdEnv,'label',label,'resc',resc,'addRowNames',addRowNames);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSCorAna:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    conflev=options.conflev;
    StoreSim=options.StoreSim;
    mmdEnv=options.mmdEnv;
    label=options.label;
    resc=options.resc;
    addRowNames=options.addRowNames;
end

if conflev<=0 || conflev>=1
    error('FSDA:FSCorAna:WrongConfInt','Confidence level must be a number in the interval (0 1)');
end

if isempty(bsb)
    ini0=round(n/5);
    if Nisstruct == true
        % In this case initial estimate of location is supplied in input structure RAW.loc;
        % Find vector of means inside subset
        ym=loc;
        bsb=bsbini;
        % Note that Niter has already been defined
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
    % Niter=N(bsb,:);
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


%  Un is a Matrix whose 2nd column:11th col contains the row(s) just
%  included.
Un = cat(2,(init1+1:n)',NaN(n-init1,3));

%  mmd has two columns
%  1st col = dimension of the subset
%  2nd col min. MD among the units not belonging to the subset
mmd=[(init1:n-1)' zeros(n-init1,1)];

% 2nd col ine=chi2/mm statistic using the units belonging to subset
ine=[(init1:n)' zeros(n-init1+1,1)];


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
        % Store value of Inertia (Chi2/m) at step mm
        ncolt=sum(Niter,1);
        Ntheo=nrowt*ncolt/mm;
        chi2mm=funzchi2(Niter(:),Ntheo(:));
        ine(mm-init1+1,2)=chi2mm/mm;
    end

    if mm<n

        % oldbsb=bsb;
        oldbsbT=bsbT;
        % sort MD distances (not squared) multiplied by row masses and by n
        [mahaldistsor,indsortdist]=sort(mahaldist.*rtimesn);


        % The rows sortdist(1:indexesCR) will be completely
        % represented;
        cumsumnjdot=cumsum(rtimesn(indsortdist));
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

            % mmd contains minimum of Mahalanobis distances *masses *n
            % among the units which are not in the subset at step m store
            % minMD
            IndexminMD=find(cumsumnjdot>=mm+1,1,'first');

            mmd(mm-init1+1,2)= mahaldistsor(IndexminMD);


            if (lunit<=2)
                if isempty(unit) && (mm-init1)>1
                    Un(mm-init1+1,2:end-1)=Un(mm-init1,2:end-1);

                else
                    Un(mm-init1+1,2:(lunit+1))=unit;
                    if ~isempty(unit)
                        if length(unit)>1
                            if msg==true
                                disp(['Warning: more than one unit entered in step' num2str(mm)])
                            end
                        end
                        Un(mm-init1+1,end)=unit(1);
                    end
                end


            else
                %if msg==1
                %    disp(['Warning: interchange greater than 1 when m=' int2str(mm)]);
                %    disp(['Number of units which entered=' int2str(lunit)]);
                %end
                Un(mm-init1+1,2:length(unit)+1)=unit(1:length(unit));
            end
        end
    end
end

out=struct;

if isempty(mmdEnv)
    % Plot mmd with envelopes
    % quant=[0.01;0.5;conflev];
    quant=[0.05;0.5;0.95;0.01;conflev];
    % Compute theoretical envelops for minimum Mahalanobis distance based on all
    % the observations for the above quantiles.
    if msg==true
        disp('Creating empirical confidence band for minimum (weighted) Mahalanobis distance')
    end
    [gmin,gine,nsimul] = FSCorAnaenv(RAW,'prob',quant,'init',init1);
    if StoreSim ==true
        out.mmdEnv=gmin;
        out.ineEnv=gine;
    end

else
    % Use precalculated empirical confidence envelope of min Mahalanobis
    % distance
    gmin=mmdEnv;
    if size(gmin,1)>=size(mmd,1)
        gmin=gmin(end-size(mmd,1)+1:end,:);
    else
        error('FSDA:FSCorAna:WrongInputOpt','Empirical precalculated envelope cannot be used because it has a smaller size than mmd.');
    end
    nsimul=0;
end


if resc==true
    warmup=500;
    warmup=min([find(mmd(:,1)>round(n/2),1),warmup]);
    % warmup=200;
    % warmup=round(n*0.6-n/4);
    if sum(mmd(1:warmup,2)<gmin(1:warmup,2))>warmup/2   || sum(mmd(1:warmup,2)>gmin(1:warmup,4))>warmup/2
        %     % coeff=mean(gmin(1:warmup,3)-mmd(1:warmup,2));
        coeff=mean(gmin(1:warmup,3))-mean(mmd(1:warmup,2));
        gmin(:,2:end)=gmin(:,2:end)-coeff;
    end

end
thresh=max(gmin(:,end));
% Outlier detection based on Bonferroni threshold
sign=find(mmd(:,2)>thresh,1);
if isempty(sign)
    out.outliers=[];
else
    out.outliers=unique(Un(sign:end,2));
end

% The above was
% fac=sqrt(r)*n;
% mmd(:,2)=(mmd(:,2)./fac).^2;


good=setdiff(1:I,out.outliers);
Ngood=N(good,:);
% ini0 is the sum of elements in the contingency table Nini
finaln=sum(Ngood,'all');
loc=sum(Ngood,1)/finaln;
md=mahalCorAna(ProfileRows,loc);



% Plot minimum Mahalanobis distance with 1%, 50% and conflev envelopes
if plots==1
    figure;

    plot(mmd(:,1),mmd(:,2),'tag','data_mmd');

    % include specified tag in the current plot
    set(gcf,'tag','pl_mmd');
    set(gcf,'Name', 'Monitoring of Minimum (weighted) Mahalanobis distance', 'NumberTitle', 'off');

    lwdenv=2;
    % Superimpose 50% envelope
    line(gmin(:,1),gmin(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');

    show5and95=false;
    if show5and95==true
        % Superimpose 5% and 95% envelope
        line(gmin(:,1),gmin(:,2),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        line(gmin(:,1),gmin(:,4),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        ylabel('')

         if addRowNames == true
               sel=~ismissing(Un(:,end));
             if isempty(label)
                text(mmd(sel,1),mmd(sel,end)*1.05,num2str(Un(sel,end)));
            else
                text(mmd(sel,1),mmd(sel,end)*1.05,label(Un(sel,end)));
            end

        end
    else

        % Superimpose 1% and conflev% envelope
        line(gmin(:,1),gmin(:,5),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        line(gmin(:,1),gmin(:,end),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
        yline(max(gmin(:,end)))
        sel=~ismissing(Un(:,end));

        if addRowNames == true
            if isempty(label)
                text(mmd(sel,1),mmd(sel,end)*1.05,num2str(Un(sel,end)));
            else
                text(mmd(sel,1),mmd(sel,end)*1.05,label(Un(sel,end)));
            end

        end

    end

    xlabel('Subset size m');
    ylabel('Monitoring of minimum (weighted) Mahalanobis distance');


    figure;

    % Plot of total inertia
    plot(ine(:,1),ine(:,2:end),'tag','data_in');

    % include specified tag in the current plot
    set(gcf,'tag','pl_in');
    set(gcf,'Name', 'Monitoring of inertia', 'NumberTitle', 'off');

    lwdenv=2;
    % Superimpose 50% envelope
    line(gine(:,1),gine(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
    % Superimpose alpha1% and alpha2% envelope
    line(gine(:,1),gine(:,2),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    line(gine(:,1),gine(:,4),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');

    xlabel('Subset size m');
    ylabel('Monitoring of inertia');


end

% Store also 1 per cent, 50 per cent  and conflev end per cent envelopes
mmd=[mmd gmin(:,[5 3 end])];

if isempty(mmdEnv)
    % Store also the envelopes for inertia explained
    ine=[ine gine(:,2:4)];
end

out.mmd=mmd;
out.ine=ine;
out.Un=Un;
out.N=N;
out.loc=loc;
out.class='FSCorAna';
out.Y=ProfileRows;
out.md=md;
out.thresh=thresh;
out.conflev=conflev;
out.simulateUnderH0=simulateUnderH0;
out.nsimul=nsimul;
end
%FScategory:MULT-Multivariate
