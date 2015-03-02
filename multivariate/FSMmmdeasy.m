function [mmd,Un,varargout] = FSMmmdeasy(Y,bsb,varargin)
%FSMmmdeasy is exactly equal to minMD but much less efficient
%
%<a href="matlab: docsearchFS('FSMmmdeasy')">Link to the help function</a>
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
% mmd:          (n-init) x 2 matrix which contains the monitoring of minimum
%               Mahalanobis distance each step of the forward search.
%               1st col = fwd search index (from init to n-1).
%               2nd col = minimum Mahalanobis distance.
% Un:           (n-init) x 11 Matrix which contains the unit(s) included
%               in the subset at each step of the search.
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one.
%               Un(1,2) for example contains the unit included in step
%               init+1.
%               Un(end,2) contains the units included in the final step
%               of the search.
% varargout :   Matrix containing units belonging to subset in each step of
%               the search (from step init to n)
%
% See also FSMenvmmd.m, FSM.m
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
%<a href="matlab: docsearchFS('FSMmmdeasy')">Link to the help function</a>
% Last modified 06-Feb-2015


% Examples:

%{
    % Run this code to see the output shown in the help file
    n=200;
    v=3;
    m0=4;
    randn('state',123456);
    Y=randn(n,v);
    %Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [fre]=unibiv(Y);
    %create an initial subset with the 4 observations, which fell the smallest
    %number of times outside the robust bivariate ellipses, and with the
    %lowest Mahalanobis Distance.
    fre=sortrows(fre,[3 4]);
    bs=fre(1:m0,1);
    [mmd,Un,BB]=FSMmmd(Ycont,bs,'plots',1);
%}

% rows(Y)
[n,v]=size(Y);
% Initialize matrix which will contain Mahalanobis distances in each step
seq=(1:n)';

%% Input parameters checking

hdef=floor(n*0.6);

options=struct('init',hdef,'plots',0,'msg',1);

if nargin<2
        error('FSDA:FSMmmdeasy:missingInputs','Initial subset is missing')
end

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSMmmdeasy:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
        warning('FSDA:FSMmmdeasy:NoFullRank','Unable to randomly sample full rank matrix');
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


if nargout==3
    % Initialize matrix which will contain the units forming subset in each
    % step of the fwd search
    BB=NaN(n,n-init1+1,'single');
end


%  Un is a Matrix whose 2nd column:11th col contains the unit(s) just
%  included.
Un = cat(2 , (init1+1:n)' , NaN(n-init1,10));

%  mmd has two columns
%  1st col = dimension of the subset
%  2nd col min. Mahalanobis distances among the units
%  which do not belong to the subset
mmd=[(init1:n-1)' zeros(n-init1,1)];

mala=[seq zeros(n,1)];

if (rank(Y(bsb,:))<v)
    warning('FSDA:FSMmmdeasy:NoFullRank','The supplied initial subset is not full rank matrix');
    disp('FS loop will not be performed')
    mmd=NaN;
    Un=NaN;
    varargout={NaN};
    % FS loop will not be performed
else
    
    for mm = ini0:n
        
        % select subset
        Yb=Y(bsb,:);
        
        if (mm>=init1) && nargout==3;
            BB(bsb,mm-init1+1)=bsb;
        end
        
        % Find vector of means inside subset
        % Note that ym is a row vector
        ym=mean(Yb);
        
        % Squared Mahalanobis distances computed using QR decomposition
        % Ym=Y-one*ym;
        Ym = bsxfun(@minus,Y, ym);
        
        [~,R]=qr(Ym(bsb,:),0);
        
        
        % Remark: u=(Ym/R)' should be much faster than u=inv(R')*Ym';
        u=(Ym/R)';
        % Compute square Mahalanobis distances
        mala(:,2)=((mm-1)*sum(u.^2,1))';
        
        zs= sortrows(mala,2);
        
        if mm<n
            % eval('mm');
            
            if (mm>=init1);
                ncl=setdiff(seq,bsb);
                sncl=sortrows(mala(ncl,2));
                % mmd contains minimum of Mahalanobis distances among
                % the units which form the group of potential outliers
                mmd(mm-init1+1,2)=sqrt(sncl(1));
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
                end;
            end
        else
        end;
        
        if nargout==3
            varargout={BB};
        end
        
        % Plot minimum Mahalanobis distance with 1%, 50% and 99% envelopes
        if options.plots==1
            quant=[0.01;0.5;0.99];
            % Compute theoretical envelops for minimum Mahalanobis distance based on all
            % the observations for the above quantiles.
            [gmin] = FSMenvmmd(n,v,'prob',quant,'init',init1);
            plot(mmd(:,1),mmd(:,2));
            
            % Superimpose 1%, 99%, 99.9% envelopes based on all the observations
            lwdenv=2;
            % Superimpose 50% envelope
            line(gmin(:,1),gmin(:,3),'LineWidth',lwdenv,'LineStyle','--','Color','g','tag','env');
            % Superimpose 1% and 99% envelope
            line(gmin(:,1),gmin(:,2),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
            line(gmin(:,1),gmin(:,4),'LineWidth',lwdenv,'LineStyle','--','Color',[0.2 0.8 0.4],'tag','env');
            
            xlabel('Subset size m');
            ylabel('Monitoring of minimum Mahlanobis distance');
        end
        
    end
end

end


