function out = FSMeda(Y,bsb,varargin)
%FSMeda performs forward search in multivariate analysis with exploratory data analysis purposes
%
%
%<a href="matlab: docsearchFS('FSMeda')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Variables. Matrix. n x v data matrix; n observations and v variables
%               Rows of Y represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
%               Data Types - single | double
% bsb :         Units forming subset. Vector. List of units forming the initial subset. If bsb=0
%               (default) then the procedure starts with v units randomly
%               chosen else if bsb is not 0 the search will start with
%               m0=length(bsb).
%               Data Types - single | double
%
% Optional input arguments:
%
% init :        It specifies the point where to start monitoring
%               required diagnostics. Scalar. Note that if bsb
%               is supplied init>=length(bsb). If init is not specified it will
%               be set equal to floor(n*0.6).
%                 Example - 'init',50 
%                 Data Types - double
% plots :    It specify whether it is necessary to produce the plots of the
%               monitoring of minMD.
%                 Scalar. If plots=1, a plot of the monitoring of minMD among
%               the units not belonging to the subset is produced on the
%               screen with 1 per cent, 50 per cent and 99 per cent confidence bands
%               else (default), all plots are suppressed.
%               Example - 'plots',0
%               Data Types - double
%  msg  :       It controls whether to display or not messages
%               about great interchange on the screen. Scalar.
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen.
%                 Example - 'msg',0 
%                 Data Types - double
% scaled:     It controls whether to monitor scaled Mahalanobis distances.
%               Scalar. If scaled=1 scaled Mahalanobis distances are
%               monitored during the search.
%                 Example - 'scaled',0 
%                 Data Types - double
% nocheck     : It controls whether to perform checks on matrix Y.Scalar. If nocheck is equal to 1 no check is performed on
%               matrix Y. As default nocheck=0.
%                 Example - 'nocheck',1
%                 Data Types - double
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
%   The output consists of a structure 'out' containing the following fields:
%   out.MAL=        n x (n-init+1) = matrix containing the monitoring of
%               Mahalanobis distances
%               1st row = distance for first unit ......
%               nth row = distance for nth unit.
%    out.BB=        n x (n-init+1) matrix containing the information about the units belonging
%               to the subset at each step of the forward search.
%               1st col = indexes of the units forming subset in the initial step
%               ...
%               last column = units forming subset in the final step (all units)
%   out.mmd=        n-init x 3 matrix which contains the monitoring of minimum
%               MD or (m+1)th ordered MD  at each step of
%               the forward search.
%               1st col = fwd search index (from init to n-1)
%               2nd col = minimum MD
%               3rd col = (m+1)th-ordered MD
%   out.msr=        n-init+1 x 3 = matrix which contains the monitoring of
%               maximum MD or mth ordered MD
%               1st col = fwd search index (from init to n)
%               2nd col = maximum MD
%               3rd col = mth-ordered MD
%    out.gap=       n-init+1 x 3 = matrix which contains the monitoring of
%               the gap (difference between minMD outside subset and max. inside)
%               1st col = fwd search index (from init to n)
%               2nd col = min MD - max MD
%               3rd col = (m+1)th ordered MD - mth ordered distance
%   out.loc=        (n-init+1) x (v+1) matrix containing the monitoring of
%               estimated of the means for each variable in each step of the forward search
%  out.S2cov=       (n-init+1) x (v*(v+1)/2+1) matrix containing the monitoring
%               of the elements of the covariance matrix in each step
%               of the forward search
%               1st col = fwd search index (from init to n)
%               2nd col = monitoring of S(1,1)
%               3rd col = monitoring of S(1,2)
%               ....
%               end col = monitoring of S(v,v)
%  out.detS=        (n-init+1) x (2) matrix containing the monitoring of
%               the determinant of the covariance matrix
%               in each step of the forward search
%    out.Un=        (n-init) x 11 Matrix which contains the unit(s)
%               included in the subset at each step of the fwd search
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Un(1,2) for example contains
%               the unit included in step init+1 Un(end,2) contains the
%               units included in the final step of the search
%     out.Y=        Original data input matrix
%
% See also FSMmmd.m, FSM.m
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
%<a href="matlab: docsearchFS('FSMeda')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:


%{
    % FSMeda with all default options.
    % Run the FS on a simulated dataset by choosing an initial subset
    % formed by the three observations with the smallest Mahalanobis
    % Distance.
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
    [out]=FSMeda(Ycont,bs);
%}

%{
    %% FSMeda with optional arguments.
    % Monitoring the evolution of minimum Mahlanobis distance.
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
    load('swiss_banknotes')
    Y=swiss_banknotes.data;
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
    Y=emilia2001.data;
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
    Y=emilia2001.data;
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

%% Input parameters checking
%chkinputM does not do any check if option nocheck=1
nnargin=nargin;
vvarargin=varargin;
Y = chkinputM(Y,nnargin,vvarargin);

if nargin<2
    error('FSDA:FSMeda:missingInputs','Initial subset is missing')
end

[n,v]=size(Y);
seq=(1:n)';
one=ones(n,1);

hdef=floor(n*0.6);
options=struct('init',hdef,'plots',0,'msg',1,'scaled',0,'nocheck',0);

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
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

scaled=options.scaled;

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
        warning('FSDA:FSMeda:NoFullRank','Unable to randomly sample full rank matrix');
    end
else
end

ini0=length(bsb);

% check init
init1=options.init;
msg=options.msg;

if  init1 <v+1;
    mess=sprintf(['Attention : init1 should be larger than v. \n',...
        'It is set to v+1.']);
    disp(mess);
    init1=v+1;
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


% Initialize matrix which will contain the units forming subset in each
% step of the fwd search
BB=NaN(n,n-init1+1);

% Initialize matrix which will contain the MD monitored in each
% step of the fwd search
MAL=NaN(n,n-init1+1);

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

% Initialize the matrix which contains the values of MD
mala=[seq zeros(n,1)];

if (rank(Y(bsb,:))<v)
    warning('FSDA:FSMeda:NoFullRank','The supplied initial subset is not full rank matrix');
    disp('FSMeda:message: FS loop will not be performed');
    out=struct;
else
    
    for mm = ini0:n
        
        % select subset
        Yb=Y(bsb,:);
        
        
        % Find vector of means inside subset
        % Note that ym is a row vector
        ym=mean(Yb);
        
        
        % Squared Mahalanobis distances computed using QR decomposition
        Ym=Y-one*ym;
        [~,R]=qr(Ym(bsb,:),0);
        
        
        % Remark: u=(Ym/R)' should be much faster than u=inv(R')*Ym';
        u=(Ym/R)';
        % Compute square Mahalanobis distances
        mala(:,2)=sqrt(((mm-1)*sum(u.^2)))';
        
        if scaled==1
            covYb=cov(Yb);
            detcovYb=det(covYb);
            mala(:,2)=mala(:,2)*(detcovYb^(1/(2*v)));
        end
        
        if (mm>=init1);
            BB(bsb,mm-init1+1)=bsb;
            
            % Store the means
            loc(mm-init1+1,2:end)=ym;
            
            % Store the trace and the determinant
            if scaled==1
                detS(mm-init1+1,2:end)=[detcovYb sum(diag(covYb))];
            else
                covYb=cov(Yb);
                detS(mm-init1+1,2:end)=[det(covYb) sum(diag(covYb))];
            end
            
            % Store the elements of the covariance matrix
            aco=triu(covYb);
            aco=aco(abs(aco(:))>1e-15);
            S2cov(mm-init1+1,2:end)=aco';
            
            % Store MD inside matrix MAL
            MAL(:,mm-init1+1)=mala(:,2);
        end
        
        zs=sortrows(mala,2);
        
        if mm>=init1;
            % Store max and mmth ordered MD
            msr(mm-init1+1,2:3)= [max(mala(bsb,2)) zs(mm,2)];
        end
        
        
        if mm<n
            % eval('mm');
            
            if (mm>=init1);
                ncl=setdiff(seq,bsb);
                sncl=sortrows(mala(ncl,2));
                % store minMD and (m+1)th MD
                mmd(mm-init1+1,2:3)=[sncl(1) zs(mm+1,2)];
                % store gap
                gap(mm-init1+1,2:3)=mmd(mm-init1+1,2:3)-msr(mm-init1+1,2:3);
            end
            
            
            
            % store units forming old subset in vector oldbsb
            oldbsb=bsb;
            
            % the dimension of subset increases by one unit.
            % vector bsb contains the indexes corresponding to the units of
            % the new subset
            bsb=zs(1:mm+1,1);
            
            if (mm>=init1);
                unit=setdiff(bsb,oldbsb);
                if (length(unit)<=10)
                    Un(mm-init1+1,2:(length(unit)+1))=unit;
                else
                    if msg==1
                        disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                        disp(['Number of units which entered=' int2str(length(unit))]);
                    end
                    Un(mm-init1+1,2:end)=unit(1:10);
                end
            end
        end
    end % close FS loop
    
    
    if scaled==1
        detcovY=(detcovYb^(1/(2*v)));
        % Rescale MD
        MAL=MAL/detcovY;
        % Rescale minimum MD
        mmd(:,2:3)=mmd(:,2:3)/detcovY;
        % Rescale max MD
        msr(:,2:3)=msr(:,2:3)/detcovY;
        
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
            ylabel('Monitoring of scaled minimum Mahlanobis distance');
        else
            ylabel('Monitoring of minimum Mahlanobis distance');
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
    out.loc=loc;
end

end

