function out = FSCorAnaeda(N,varargin)
%FSCorAnaeda performs forward search in correspondence analysis with exploratory data analysis purposes
%
%
%<a href="matlab: docsearchFS('FSCorAnaeda')">Link to the help function</a>
%
% Required input arguments:
%
%  N :  contingency table or structure.
%               Array or table of size I-by-J or strucure. If N is a
%               structure it contains the following fields:
%               N.N = contingency table in array format of size I-by-J.
%               N.loc = initial location estimate for the matrix of Profile
%               rows of the contingency table (row vector or length J).
%               Note that input structure N can be conveniently created by
%               function mcdCorAna.
%               If N is not a struct it is possible to specify the unitf
%               forming initial subset with input option bsb.
%                Data Types - single|double
%
% Optional input arguments:
%
% bsb  :       Initial subset. Vector of positive integers containing the
%              indexes of the rows of the contingency table which have to
%              be used to initialize the forward search. If bsb is empty
%              and required input argument is a struct N.loc will be used.
%              If bsb is supplied and N is a struct N.loc is ignored.
%              The default value of bsb is empty, and if N is struct
%              a random subset containing round(n/5) units will be used.
%                 Example - 'bsb',[3 6 8 10 12 14]
%                 Data Types - double
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
%  msg  :       It controls whether to display or not messages
%               about great interchange on the screen. Scalar.
%               If msg==1 (default) messages are displyed on the screen
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
%   out.mmd=    n-init-by-2 matrix which contains the monitoring of minimum
%               MD or (m+1)th ordered MD  at each step of
%               the forward search.
%               1st col = fwd search index (from init to n-1);
%               2nd col = minimum MD;
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
% Copyright 2008-2021.
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

% Input parameters checking
%chkinputM does not do any check if option nocheck=1

if isstruct(N)
    loc=N.loc;
    N=N.N;
    Nisstruct=true;
else
    Nisstruct=false;
end

if istable(N)
    N=table2array(N);
end

% n= total sample size
n=sum(N,'all');

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


init1=floor(n*0.6);
plots=0;
bsb=[];

if nargin > 1
    
    options=struct('init',init1,'plots',plots,'msg',1,'bsb',bsb);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSCorAnaeda:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
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
end

if isempty(bsb)
    ini0=round(n/5);
    if Nisstruct == true
        % In this case initial estimate of location is supplied in input structure RAW.loc;
        % Find vector of means inside subset
        ym=loc;
    else
        % A random subset of ini0 units will be extracted.
        indsamp=randsample(I,I);
        
        cumsumnjdot=cumsum(rtimesn(indsamp));
        indexesCR=find(cumsumnjdot<ini0,1,'last');
        
        
        % Find how many units must be included from
        % row indexesCR+1;
        unitstoADD=ini0-cumsumnjdot(indexesCR);
        bsb=indsamp(1:indexesCR+1);
        
        Niter=N(bsb,:);
        Niter(end,:)=unitstoADD*Niter(end,:)/rtimesn(indsamp(indexesCR+1));
        
        % Compute centroid based on the random extracted ini0 units
        ym  = sum(Niter,1)/ini0;
    end
else
    % Elements of N(bsb,:) are used to find initial estimate of location
    Nini=N(bsb,:);
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
BB=NaN(I,n-init1+1);

% Initialize matrix which will contain the MD monitored in each
% step of the fwd search
MAL=BB;

% Initialize matrix which will contain the means of the variables monitored in
% each step
Loc=cat(2,(init1:n)',NaN(n-init1+1,J));

%  Un is a Matrix whose 2nd column:11th col contains the row(s) just
%  included.
Un = cat(2,(init1+1:n)',NaN(n-init1,10));

%  mmd has three columns
%  1st col = dimension of the subset
%  2nd col min. MD among the units ont belonging to the subset
%  3rd col (m+1) ordered MD
mmd=[(init1:n-1)' zeros(n-init1,1)];


zeron1=false(I,1);
bsbT=zeron1;
bsbT(bsb)=true;

for mm = ini0:n
    
    % Mahalanobis distances from initialloc and Initialshape
    % mahaldist = r.*mahalFS(ProfileRows, initialloc, initialcov);
    MD=mahalCorAna(ProfileRows,ym);
    
    mahaldist=sqrt(MD);
    
    if mm>=init1
        
        BB(bsb,mm-init1+1)=bsb;
        
        % Store the means
        Loc(mm-init1+1,2:end)=ym;
        
        % Store the trace and the determinant
        % detS(mm-init1+1,2:end)=[detcovYb sum(diag(covYb))];
        
        
        % Store MD inside matrix MAL
        MAL(:,mm-init1+1)=mahaldist;
    end
    
    
    if mm<n
        
        % oldbsb=bsb;
        oldbsbT=bsbT;
        % sort MD distances multiplied by row masses
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
            
            % mmd contains minimum of Mahalanobis distances among
            % the units which are not in the subset at step m
            % store minMD and (m+1)th MD
            IndexminMD=find(cumsumnjdot>=mm+1,1,'first');
            
            mmd(mm-init1+1,2)= mahaldistsor(IndexminMD);
            
            
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
    end
end

% Plot minimum Mahalanobis distance with 1%, 50% and 99% envelopes
if plots==1
    figure;
    quant=[0.01;0.5;0.99];
    % Compute theoretical envelops for minimum Mahalanobis distance based on all
    % the observations for the above quantiles.
    disp('Creating empirical confidence band for minimum (weighted) Mahalanobis distance')
    [gmin] = FSCorAnaenvmmd(N,'prob',quant,'init',init1);
    
    plot(mmd(:,1),mmd(:,2),'tag','data_mmd');
    
    % include specified tag in the current plot
    set(gcf,'tag','pl_mmd');
    set(gcf,'Name', 'Monitoring of Minimum (weighted) Mahalanobis distance', 'NumberTitle', 'off');
    
    % Superimpose 1%, 99%, 99.9% envelopes based on all the observations
    lwdenv=2;
    % Superimpose 50% envelope
    line(gmin(:,1),gmin(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
    % Superimpose 1% and 99% envelope
    line(gmin(:,1),gmin(:,2),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    line(gmin(:,1),gmin(:,4),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
    
    xlabel('Subset size m');
    ylabel('Monitoring of minimum (weighted) Mahalanobis distance');
    
end

out.MAL=MAL;
out.BB=BB;
out.mmd=mmd;
% out.detS=detS;
out.Un=Un;
out.N=N;
out.Loc=Loc;
out.class='FSCorAnaeda';
out.Y=ProfileRows;

end
%FScategory:MULT-Multivariate
