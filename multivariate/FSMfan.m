function out = FSMfan(Y,la0,varargin)
%FSMfan computes confirmatory lrt of a suggested transformation
%It uses the multivariate version of the parametric family of power
%transformations.
%
%<a href="matlab: docsearchFS('FSMfan')">Link to the help function</a>
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
% la0:          Transformation parameters. Vector. Vector of length v=size(Y,2) specifying a reasonable set
%               of transformations for the columns of the multivariate data set.
%               Data Types - single | double
%
% Optional input arguments:
%
%    family :   string which identifies the family of transformations which
%               must be used. Character. Possible values are 'BoxCox' (default) or
%               'YJ'.
%               The Box-Cox family of power transformations equals
%               (y^{\lambda}-1)/\lambda for \lambda not equal to zero, and
%               log(y)
%               if \lambda = 0.
%               The Yeo-Johnson (YJ) transformation is the Box-Cox
%               transformation of y+1 for nonnegative values, and of |y|+1 with
%               parameter 2-\lambda for y negative.
%               The basic power transformation returns y^{\lambda} if \lambda is not
%               zero, and log(\lambda) otherwise.
%               Remark. BoxCox and the basic power family can be used just
%               if input y is positive. YeoJohnson family of
%               transformations does not have this limitation.
%               Example - 'family','YJ' 
%               Data Types - char
%        rf :   confidence level for bivariate ellipses. Scalar. Default is
%               0.9.
%                 Example - 'rf',0.99 
%                 Data Types - double
%       init    : Point where to start monitoring required diagnostics.
%                 Scalar. Note that if bsb is supplied init>=length(bsb). If
%                 init is not specified it will be set equal to
%                 floor(n*0.6).
%                 Example - 'init',50 
%                 Data Types - double
% ColToComp :  variables for which likelihood ratio tests have to be produced. Vector. It is a k x 1 integer vector. For
%               example, if ColToComp = [2 4], the signed likelihood ratio tests
%               are produced for the second and the fourth column of matrix Y. If
%               col.to.compare = '' then all variables (columns of matrix Y) are
%               considered.
%                 Example - 'ColToComp',[1 3] 
%                 Data Types - double
%   laAround:  It specifies for which values
%               of lambda to compute the likelihood ratio test. Scalar. It is a  r x 1 vector. If this argument is
%               omitted, the function produces for each variable specified in
%               ColToComp the likelihood ratio tests associated to the five most
%               common values of lambda [-1, -0.5, 0, 0.5, 1].
%                 Example - 'laAround',[1 0] 
%                 Data Types - double
%   optmin  :   It contains the options dealing with the
%               maximization algorithm. Structure. 
%               Use optimset to set these options.
%               Notice that the maximization algorithm which is used is
%               fminunc if the optimization toolbox is present else is
%               fminsearch.
%                 Example -'optmin.Display','off' 
%                 Data Types - double 
%     speed : It indicates the initial value of
%               the maximization procedure. Scalar. If speed=1 (default) the initial value at step m of
%               the maximization procedure (fminunc or fminsearch) is the
%               final value at step m-1 else it is la0.
%                 Example - 'speed',0
%                 Data Types - double
%   colnames: the names of the variables of the dataset. Cell array of strings. Cell array of strings of length v containing the names of
%               the variables of the dataset. If colnames is empty then the
%               sequence 1:v is created to label the variables.
%                 Example - 'colnames',{'1' '2' '3' '4' '5' '10' '11' '12' '13'}
%                 Data Types - Cell array of strings.
%     signlr:  plots of signed square root
%               likelihood ratios. Scalar. If signlr = 1 (default) plots of signed square root
%               likelihood ratios are produced, else likelihood ratios are
%               produced.
%                 Example - 'signlr',0
%                 Data Types - double
%   plotslrt:   It specifies whether it is necessary to
%               plot the (signed square root) likelihood ratio test. Scalar or structure.
%               If plotslrt is a scalar, the plot of the monitoring of
%               likelihood ratio test is produced on the screen with all
%               default options.
%               If plotslrt is a strucure, it may contain the following fields: 
%                   plotslrt.xlim     = minimum and maximum on the x axis; 
%                   plotslrt.ylim     = minimum and maximum on the y axis; 
%                   plotslrt.LineWidth= Line width of the trajectory of lrt of
%                                       transformation parameters; 
%                   plotslrt.conflev  = vector which defines the confidence
%                                       levels of the horizontal line for
%                                       the likelihood ratio test (default
%                                       is conflev=[0.95 0.99]);
%                   plotslrt.LineWidthEnv= Line width of the horizontal lines; 
%                   plotslrt.Tag      = tag of the plot (default is pl_lrt). 
%                 Example -'plotslrt',1
%                 Data Types - double
%                 Example - 'plotslrt',struct
%                 Data Types - double
%  msg  :       It controls whether to display or not messages
%               about great interchange on the screen. Scalar.
%               If msg==1 (default) messages are displyed on the screen
%               else no message is displayed on the screen.
%                 Example - 'msg',0
%                 Data Types - double
%
%
% Remark:       The user should only give the input arguments that have to
%               change their default value.
%               The name of the input arguments needs to be followed by
%               their value. The order of the input arguments is of no
%               importance.
%
%
% Output:
%
%         out:   structure which contains the following fields
%
%    out.LRT=   Cell of length ColtoComp. Each element of the cell contains the
%               a matrix of size n-init+1 x length(laAround)+1 which
%               contains the monitoring of (signed square root) likelihood
%               ratio for testing H0:\lambda_j=la0_j when all the other
%               variables are transformed as specified in vector la0.
%               More precisely each
%               matrix of size n-init+1 x length(laAround)+1 presents the
%               following structure:
%               1st col = fwd search index (from init to n);
%               2nd col = value of the (signed sqrt) likelihood ratio for
%               testing laj=laAround(1);
%               ...
%               length(laAround)+1 col = value of the (signed sqrt) likelihood ratio for
%               testing laj=laAround(end).
%   out.Exflag= Cell of length ColtoComp. Each element of the cell contains the
%               a matrix of size n-init+1 x length(laAround)+1 which
%               contains the monitoring of the
%               integer identifying the reason why the maximization
%               algorithm terminated. See help page fminunc of the
%               optimization toolbox  for the list of values of exitflag
%               and the corresponding reasons the algorithm terminated.
%               More precisely each
%               matrix of size n-init+1 x length(laAround)+1 presents the
%               following structure: 
%               1st col = fwd search index (from init to n); 
%               2nd col = integer identifying the reason the algorithm terminated
%               when testing laj=laAround(1); 
%               ...
%               length(laAround)+1 col = integer identifying the reason the algorithm terminated
%               when testing laj=laAround(end).
%    out.Un=    Cell of length ColtoComp. Each element of the cell contains the
%               a (sub)cell of size length(laAround). Each element of the (sub)cell
%               contains a  n-init+1 x 11 which informs the order of entry of the units
%               For example Unj=Un{i}{j} refers to ColtoComp(i) and laAround(j)
%               Unj is a (n-init) x 11 matrix which contains the unit(s)
%               included in the subset at each step of the fwd search. 
%               REMARK: in every step the new subset is compared with the
%               old subset. Un contains the unit(s) present in the new
%               subset but not in the old one Unj(1,2) for example contains
%               the unit included in step init+1 Unj(end,2) contains the
%               units included in the final step of the search
%
% See also FSMtra.m, FSM.m
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
%<a href="matlab: docsearchFS('FSMtra')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % FSMfan with all default options.
    % Example with Mussels data.
    load('mussels.mat');
    Y=mussels.data;
    % FS based on with H_0:\lambda=[1 0.5 1 0 1/3]
    [out]=FSMfan(Y,[0.5 0 0.5 0 0]);
%}

%{
    % FSMfan with otpional arguments.
    % Example with Mussels data.
    load('mussels.mat');
    Y=mussels.data;
    % FS based on with H_0:\lambda=[1 0.5 1 0 1/3]
    plotslrt=struct;
    plotslrt.ylim=[-6.2 6.2];
    [out]=FSMfan(Y,[0.5 0 0.5 0 0],'laAround',[-1 -0.5 0 1/3 0.5 1],'init',58,'plotslrt',plotslrt);
    % Compare this plot with Figure 4.24 p. 182 of ARC (2004)
%}


%{
    %% EmiliaRomagna data (demographic variables).
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Extract demographic variables
    Y1=Y(:,[1 2 3 4 5 10 11 12 13]);
    colnames={'1' '2' '3' '4' '5' '10' '11' '12' '13'};
    plotslrt=struct;
    plotslrt.ylim=[-8.2 8.2];
    la0=[0 0.25 0 0.5 0.5 0 0 0.5 0.25];
    ColToComp=[1 3 5 9];
    [out]=FSMfan(Y1,la0,'ColToComp',ColToComp,'plotslrt',plotslrt,'colnames',colnames);
    % Compare the plot Figure 4.35 p. 192 of ARC (2004)
%}

%{
    %% Emilia Romagna data (modified wealth variables), example 1.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify wealth variables
    Y(:,16)=100-Y(:,16);
    Y(:,23)=100-Y(:,23);
    % Extract wealth variables
    Y1=Y(:,[14:23]);
    colnames={'14' '15' '16' '17' '18' '19' '20' '21' '22' '23'};
    plotslrt=struct;
    plotslrt.ylim=[-8.2 8.2];
    la0=[0 1 0.25 1 1 0.5 -0.5 0.25 0.25 -1];
    ColToComp=[1 7];
    [out]=FSMfan(Y1,la0,'ColToComp',ColToComp,'plotslrt',plotslrt,'colnames',colnames);
    % Compare the plot with the two upper panels of Figure 4.38 p. 188 of ARC (2004)
%}

%{
    % Emilia Romagna data (modified wealth variables), example 2.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify wealth variables
    Y(:,16)=100-Y(:,16);
    Y(:,23)=100-Y(:,23);
    % Extract wealth variables
    Y1=Y(:,[14:23]);
    colnames={'14' '15' '16' '17' '18' '19' '20' '21' '22' '23'};
    plotslrt=struct;
    plotslrt.ylim=[-7 7];
    la0=[0 1 0.25 1 1 0.5 -0.5 0.25 0.25 -1];
    ColToComp=[3 9];
    laAround=[0 0.25 1/3 0.5];
    [out]=FSMfan(Y1,la0,'laAround',laAround,'ColToComp',ColToComp,'plotslrt',plotslrt,'colnames',colnames);
    % Compare the plot with the two bottom panels of Figure 4.39 p. 195 of ARC (2004)
%}

%{
    % Emilia Romagna data with Yeo and Johnson parametric family
    load('emilia2001')
    Y=emilia2001.data;

    % Modify wealth variables
    Y(:,16)=100-Y(:,16);
    Y(:,23)=100-Y(:,23);
    % Extract wealth variables
    Y1=Y(:,[14:23]);
    colnames={'14' '15' '16' '17' '18' '19' '20' '21' '22' '23'};
    plotslrt=struct;
    plotslrt.ylim=[-7 7];
    la0=[0 1 0.25 1 1 0.5 -0.5 0.25 0.25 -1];
    ColToComp=[3 9];
    laAround=[0 0.25 1/3 0.5];
    [out]=FSMfan(Y1,la0,'laAround',laAround,'ColToComp',ColToComp,'plotslrt',plotslrt,'colnames',colnames,'family','YJ');
%}

%{
    % Emilia Romagna data (modified work variables), example 1.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21 25 26];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Extract work variables
    Y1=Y(:,[6:9 24:28]);
    colnames={'y6' 'y7' 'y8' 'y9' 'y24' 'y25' 'y26' 'y27' 'y28'};
    la0=[0.25,0,2,-1,0,1.5,0.5,1,1];
    plotslrt=struct;
    plotslrt.ylim=[-8.2 8.2];
    ColToComp=[1:4];
    laAround=[-1 -0.5 0 0.25 0.5 1 2];
    [out]=FSMfan(Y1,la0,'ColToComp',ColToComp,'laAround',laAround,'plotslrt',plotslrt,'colnames',colnames);
    % Compare the plot with Figure 4.43 p. 198 of ARC (2004)
%}


%{
    % Emilia Romagna data (modified work variables), example 2.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify variables 25 and 26
    Y(:,25)=100-Y(:,25);
    Y(:,26)=100-Y(:,26);
    % Extract work variables
    Y1=Y(:,[6:9 24:28]);
    colnames={'y6' 'y7' 'y8' 'y9' 'y24' 'y25' 'y26' 'y27' 'y28'};
    la0=[0.25,0,2,-1,0,0,1.5,1,1];
    plotslrt=struct;
    plotslrt.ylim=[-8.2 8.2];
    ColToComp=[6 7];
    laAround=[-1 -0.5 0 0.25 0.5 1 1.5 2];
    [out]=FSMfan(Y1,la0,'ColToComp',ColToComp,'laAround',laAround,'plotslrt',plotslrt,'colnames',colnames);
    % Compare the plot with Figure 4.44 p. 199 of ARC (2004)
%}


%{
    % Emilia Romagna data (all variables).
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify variables y16 y23 y25 y26
    sel=[16 23 25 26];
    Y(:,sel)=100-Y(:,sel);
                        
    colnames={'y6' 'y7' 'y8' 'y9' 'y24' 'y25' 'y26' 'y27' 'y28'};
    la0demo=[0,0.25,0,0.5,0.5,0,0,0.5,0.25];
    la0weal=[0,1,0.25,1,1,0.5,-0.5,0.25,0.25,-1];
    la0work=[0.25,0,2,-1,0,0,1,1,1];
    la0C1=[la0demo(1:5) la0work(1:4) la0demo(6:9) la0weal la0work(5:9)];
    plotslrt=struct;
    plotslrt.ylim=[-8.2 8.2];
    ColToComp=[8 9 14 25];
    laAround=[-1 -0.5 0 0.25 0.5 1 1.5 2];
    [out]=FSMfan(Y,la0C1,'ColToComp',ColToComp,'laAround',laAround,'plotslrt',plotslrt,'init',100);
    % Compare the plot with Figure 4.45 p. 199 of ARC (2004)
%}

%{
    % Emilia Romagna data (all variables) with Yeo and Johnson parametric
    % family.
    load('emilia2001')
    Y=emilia2001.data;
    % Replace zeros with min values for variables specified in sel
    sel=[6 10 12 13 19 21];
    for i=sel
        Y(Y(:,i)==0,i)=min(Y(Y(:,i)>0,i));
    end
    % Modify variables y16 y23 y25 y26
    sel=[16 23 25 26];
    Y(:,sel)=100-Y(:,sel);
                        
    colnames={'y6' 'y7' 'y8' 'y9' 'y24' 'y25' 'y26' 'y27' 'y28'};
    la0demo=[0,0.25,0,0.5,0.5,0,0,0.5,0.25];
    la0weal=[0,1,0.25,1,1,0.5,-0.5,0.25,0.25,-1];
    la0work=[0.25,0,2,-1,0,0,1,1,1];
    la0C1=[la0demo(1:5) la0work(1:4) la0demo(6:9) la0weal la0work(5:9)];
    plotslrt=struct;
    plotslrt.ylim=[-8.2 8.2];
    ColToComp=[8 9 14 25];
    laAround=[-1 -0.5 0 0.25 0.5 1 1.5 2];
    [out]=FSMfan(Y,la0C1,'ColToComp',ColToComp,'laAround',laAround,'plotslrt',plotslrt,'init',100,'family','YJ');
%}

%% Input parameters checking
% Extract size of the data
[n,v]=size(Y);
% Initialize matrix which will contain Mahalanobis distances in each step
seq=(1:n)';
one=ones(n,1);

if nargin<1
    error('FSDA:FSMfan:missingInputs','Initial data matrix is missing')
end

if nargin<2
    error('FSDA:FSMfan:missingInputs','Vector la0 is missing');
end

if length(la0)~=v
    error('FSDA:FSMfan:WrongLambda',['length of vector la0 (length(la0)=' num2str(length(la0)) ') is not equal' ...
        ' to the number of variables of the dataset (size(Y,2)=' num2str(v) ')']);
end

family='BoxCox';
hdef=floor(n*0.6);
options=struct('rf',0.9,'init',hdef,'ColToComp','','laAround',-1:0.5:1,'onelambda',0,'signlr',1,...
    'speed',1,'optmin',optimset,'msg',1,'colnames','',...
    'plotslrt','','family',family);


UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSMfan:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

%init1=options.init;
if nargin > 1
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    family=options.family;
end

Xnormfamily=strcat('norm',family);
familytra=str2func(Xnormfamily);

if strcmp(family,'BoxCox')
    loglik=@likBoxCox;
elseif strcmp(family,'YJ')
    loglik=@likYJ;
else
    warning('FSDA:FSMfan:WrongFamily','Transformation family which has been chosen is not supported')
    error('FSDA:FSMfan:WrongFamily','Supported values are BoxCox or YaoJohnson')
end


% And check if the optional user parameters are reasonable.
% check init
init1=options.init;
msg=options.msg;

if  init1 <v+1;
    mess=sprintf(['Attention : init1 should be larger than v. \n',...
        'It is set to v+1.']);
    disp(mess);
    init1=v+1;
elseif init1>n;
    mess=sprintf(['Attention : init1 should be smaller than n. \n',...
        'It is set to n.']);
    disp(mess);
    init1=n;
end




% ColToTra = vector which contains the numbers associated to the variables
% which have to expanded using laAround
ColToComp=options.ColToComp;
if isempty(ColToComp)
    ColToComp=(1:v);
elseif max(ColToComp)>v || min(ColToComp)<1
    error('FSDA:FSMfan:WrongColToComp',['The columns to test are not in the range. 1-' num2str(v)]);
end

% Ytr = matrix which contains Box Cox transformed values using la0
Ytr=feval(familytra,Y,1:v,la0);

% laAround= values of transformation parameter(s) to do the expansion
laAround=options.laAround;


% if speed=1 starting value for minimization at step m>init is final value at
% step m-1 else vector la0 is used as starting value in each step
speed=options.speed;

% Options to use in the minimization;
optmin=options.optmin;
optmin.Display='off';
% optmin=optimset('Display','off');
% optmin=optimset('Display','on');
% optmin=optimset(defopt,newopts);
%warning('off','optim:fminunc:SwitchingMethod');

% Check if optimization toolbox is installed in current computer
typemin=exist('fminunc','file');


if ~isempty(options.colnames)
    colnames=options.colnames;
else
    colnames=cellstr(num2str((1:v)'));
end

% Specify whether to store the likelihood ratio test or the its signed
% square root
signlr=options.signlr;

% Confidence level for robust bivariate ellipses
rf=options.rf;

% llaAround = length of vector laAround
llaAround=length(laAround);

% Initialize the matrix which contains the values of MD based on
% transformed obervations Ytr
mala=[seq zeros(n,1)];

% cColToComp = number of columns of the datasets to expand
cColToComp=length(ColToComp);

% LRT is a cell which in position i contains the matrix of loglikelihood
% ratio for variable cColToComp(i) expanded according to the values
% specified in vector laAround
LRT=cell(cColToComp,1);

% Exflag has the same structure of cell LRT but it contains the integers
% identyfing why the maximization terminated
Exflag=LRT;

% Un is a cell of length ColtoComp.
% Each element of the cell contains a set of length(laAround) cells. Each cell contains the matrix
% which informs about
% the order of entry of the units
% For example Unj=Un{2} refers to ColtoComp(2)
% and Unj{3} or Un{2}{3} refers to laAround{3} for ColtoComp(2)
Un=LRT;

% indexjj;
indexjj=1;

% Initialize matrix which stores in each step the integer identifying the
% reason why the algorithm terminated
Exflagjj=cat(2,(init1:n)',NaN(n-init1+1,llaAround));

% Initialize matrix which will contain the likelihood ratio for
% \lambda=\lambda0 monitored in each step
LIKratjj=Exflagjj;


Unjj=cell(llaAround,1);

Unjinit=cat(2,(init1+1:n)',NaN(n-init1,10));

% jj loops over the columns of the datasets specified in ColToComp
for jj=ColToComp;
    
    % Yjj contains in column jj the original untransformd variable and in all the other
    % columns the variables transformed using la0.
    % This matrix will be the input of the maximization in each step
    Yjj=Ytr;
    Yjj(:,jj)=Y(:,jj);
    
    indexlai=1;
    
    for lai=laAround
        
        %  Unj is a Matrix whose 2nd column:11th col contains the unit(s) just
        %  included.
        Unj = Unjinit;
        
        % Ytrlai contains transformed data using lai for column jj and data
        % transformed using la0 for all the other variables
        % This matrix will be used to compute MD at each step and in procedure
        % unibiv to find initial subset
        % Ytrlai=Ytr;
        Ytrlai=feval(familytra,Yjj,jj,lai);
        
        % Find initial subset to initialize the search
        [fre]=unibiv(Ytrlai,'rf',rf);
        fre=sortrows(fre,[3 4 2]);
        bsb=fre(1:v+1,1);
        ini0=v+1;
        
        
        if (rank(Ytrlai(bsb,:))<v)
            warning('FSDA:FSMfan:NoFullRank','The supplied initial subset is not full rank matrix');
            % FS loop will not be performed
            % out=struct;
        else
            
            % lainit starting value for minimization at step m=init
            lainit=la0(jj)+1e-2*randn(1,1);
            laout=lainit;
            
            
            for mm = ini0:n
                
                % Ytrb =subset of transformed values
                Ytrlaib=Ytrlai(bsb,:);
                
                % Yb = subset of matrix Yjj
                % Remember that Yjj contains in column jj the original untransformd variable and in all the other
                % columns the variables transformed using la0.
                % This matrix will be the input of the maximization in each step
                Yb=Yjj(bsb,:);
                
                % Find vector of means of subset of transformed values
                % Note that ym is a row vector
                ym=mean(Ytrlaib);
                
                % Squared Mahalanobis distances computed using QR decomposition
                % based on transformed values
                Ym=Ytrlai-one*ym;
                [~,R]=qr(Ym(bsb,:),0);
                
                
                % Remark: u=(Ym/R)' should be much faster than u=inv(R')*Ym';
                u=(Ym/R)';
                % Compute square Mahalanobis distances on transformed data
                mala(:,2)=sqrt(((mm-1)*sum(u.^2)))';
                
                
                if (mm>=init1);
                    
                    if speed==1;
                        lainit=laout;
                    end
                    
                    
                    try
                        fval='';
                        
                        % In order to find minimum it is possible to use MATLAB
                        % function fminsearch or function fminunc from the optimization
                        % toolbox
                        if typemin==2;
                            [laout,fval,exitflag]  = fminunc(loglik,lainit,optmin);
                        else
                            [laout,fval,exitflag]  = fminsearch(loglik,lainit,optmin);
                        end
                        %laout=5*sin(laout);
                    catch
                        disp(['Warning: non convergence at step mm=' num2str(mm)])
                        disp(['Variable' colnames(jj) ' laAround=' num2str(lai) ])
                        
                        if isempty(fval)
                            exitflag=nan;
                            laout=nan;
                            fval=nan;
                        end
                    end
                    
                    % Store transformed data for variable jj using subset
                    Ytrb0=feval(familytra,Yb,jj,lai);
                    
                    lrat=mm*(log(det(cov(Ytrb0)))-fval);
                    
                    if lrat<0;
                        disp(['Warning: negative lrt, non convergence at step mm=' num2str(mm)])
                        disp(['Variable' colnames(jj) ' laAround=' num2str(lai) ])
                        
                        exitflag=nan;
                        laout=lainit;
                    else
                        
                        
                        % Store reason why minimization stopped
                        Exflagjj(mm-init1+1,indexlai+1)=exitflag;
                        
                        if signlr==1
                            LIKratjj(mm-init1+1,indexlai+1)= sign(laout-lai)*sqrt(lrat);
                        else
                            LIKratjj(mm-init1+1,indexlai+1)=lrat;
                        end
                    end
                end
                
                zs=sortrows(mala,2);
                if mm<n
                    
                    % store units forming old subset in vector oldbsb
                    oldbsb=bsb;
                    
                    % the dimension of subset increases by one unit.
                    % vector bsb contains the indexes corresponding to the units of
                    % the new subset
                    bsb=zs(1:mm+1,1);
                    
                    if (mm>=init1);
                        unit=setdiff(bsb,oldbsb);
                        if (length(unit)<=10)
                            Unj(mm-init1+1,2:(length(unit)+1))=unit;
                        else
                            if msg==1
                                disp(['Warning: interchange greater than 10 when m=' int2str(mm)]);
                                disp(['Number of units which entered=' int2str(length(unit))]);
                            end
                            Unj(mm-init1+1,2:end)=unit(1:10);
                        end
                    end
                end % close if mm<nn
            end % close FS loop
        end % close if relative to the rank
        
        % Store matrix which contains entry order
        Unjj{indexlai}=Unj;
        
        % Store values of lrt associated with lai
        % Unall
        indexlai=indexlai+1;
    end   % Close loop associated with laAround
    
    LRT{indexjj}=LIKratjj;
    Un{indexjj}=Unjj;
    Exflag{indexjj}=Exflagjj;
    indexjj=indexjj+1;
    
end % Close loop associated with columns of the datasets specified in ColToComp

% lik computes the likelihood when different lambdas are possible for
% different variables
    function dZ=likBoxCox(laj)
        Z=Yb;
        
        Gj=exp(mean(log(Yb(:,jj))));
        
        if laj~=0
            Z(:,jj)=(Yb(:,jj).^laj-1)/(laj*(Gj^(laj-1)));
        else
            Z(:,jj)=Gj*log(Yb(:,jj));
        end
        
        dZ=log(det(cov(Z)));
        % disp(dZ);
    end

% lik computes the likelihood when different lambdas are possible for
% different variables
    function dZ=likYJ(laj)
        
        Z=Yb;
        
        nonnegs = Yb(:,jj) > 0;
        negs = ~nonnegs;
        % YJ transformation is the Box-Cox transformation of
        % y+1 for nonnegative values of y
        if laj ~=0
            Z(nonnegs,jj)= ((Yb(nonnegs,jj)+1).^laj-1)/laj;
        else
            Z(nonnegs,jj)= log(Yb(nonnegs,jj)+1);
        end
        
        % YJ transformation is the Box-Cox transformation of
        %  |y|+1 with parameter 2-lambda for y negative.
        if 2-laj~=0
            Z(negs,jj) = - ((-Yb(negs,jj)+1).^(2-laj)-1)/(2-laj);
        else
            Z(negs,jj) = -log(-Yb(negs,jj)+1);
        end
        
        % transformation is normalized so that its
        % Jacobian will be 1
        Z(:,jj)=Z(:,jj) * (exp(mean(log(   (1 + abs(Yb(:,jj))).^(2 * nonnegs - 1)) )))^(1 - laj);
        
        dZ=log(det(cov(Z)));
        % disp(dZ);
    end



%% Plot of (signed sqrt root) LRT
plotslrt=options.plotslrt;

if ~isempty(plotslrt)
    
    if  cColToComp==1
        nr=1;
        nc=1;
    elseif cColToComp==2
        nr=2;
        nc=1;
    elseif cColToComp<=4;
        nr=2;
        nc=2;
    elseif cColToComp<=6
        nr=3;
        nc=2;
    elseif cColToComp<=9
        nr=3;
        nc=3;
    elseif cColToComp<=16
        nr=4;
        nc=4;
    elseif cColToComp<=25
        nr=5;
        nc=5;
    else
        error('FSDA:FSMfan:TooManyVars','plot of lrt cannot be displayed for more than 25 variables')
    end
    
    if isstruct(plotslrt)
        
        fplotslrt=fieldnames(plotslrt);
        
        d=find(strcmp('xlim',fplotslrt));
        if d>0
            xlimx=plotslrt.xlim;
        else
            xlimx='';
        end
        
        d=find(strcmp('ylim',fplotslrt));
        if d>0
            ylimy=plotslrt.ylim;
        else
            if signlr==1;
                ylimy=[-5 5];
            else
                ylimy=[0 25];
            end
        end
        
        d=find(strcmp('LineWidth',fplotslrt));
        if d>0
            LineWidth=plotslrt.LineWidth;
        else
            LineWidth=2;
        end
        
        % LineWidthEnv = line width of the
        % horizontal lines associated with the asymptotic confidence levels
        % for the likelihood ratio test
        
        d=find(strcmp('LineWidthEnv',fplotslrt));
        if d>0
            LineWidthEnv=plotslrt.LineWidthEnv;
        else
            LineWidthEnv=1;
        end
        
        d=find(strcmp('conflev',fplotslrt));
        if d>0
            conflev=plotslrt.conflev;
        else
            conflev=0.99;
        end
        
        d=find(strcmp('Tag',fplotslrt));
        if d>0
            tag=plotslrt.Tag;
        else
            tag='pl_lrt';
        end
        
        % Specify the line type for the trajectories of MLE of
        % transformation parameters
        d=find(strcmp('LineStyle',fplotslrt));
        if d>0
            LineStyle=plotslrt.LineStyle;
        else
            slin=repmat({'-';'--';':';'-.'},ceil(llaAround/4),1);
            LineStyle=slin(1:llaAround);
        end
        
        d=find(strcmp('FontSize',fplotslrt));
        if d>0
            FontSize=plotslrt.FontSize;
        else
            FontSize=12;
        end
        
        
    else
        
        xlimx='';
        if signlr==1;
            ylimy=[-5 5];
        else
            ylimy=[0 25];
        end
        
        LineWidth=2;
        LineWidthEnv=1;
        tag='pl_fanm';
        conflev=0.99;
        FontSize=12;
        slin=repmat({'-';'--';':';'-.'},ceil(llaAround/4),1);
        LineStyle=slin(1:llaAround);
        
    end
    
    % Plot of the likelihood ratio
    % Specify where to send the output of the monitoring of lrt
    hlrt=findobj('-depth',1,'tag',tag);
    if (~isempty(hlrt))
        clf(hlrt);
        figure(hlrt)
        axes;
    else
        figure;
        set(gcf,'Name','Likelihood ratio');
    end
    
    ij=1;
    
    for jj=ColToComp;
        
        LIKrat=LRT{ij};
        subplot(nr,nc,ij)
        
        plot1=plot(LIKrat(:,1),LIKrat(:,2:end),'LineWidth',LineWidth);
        
        if ~isempty(xlimx)
            xlim(xlimx);
        end
        
        if ~isempty(ylimy)
            ylim(ylimy);
        end
        
        % Add labels at the end of the search only if they are inside the
        % limits
        for i=1:llaAround;
            if LIKrat(end,i+1)>=ylimy(1) && LIKrat(end,i+1)<=ylimy(2)
                text(n,LIKrat(end,i+1),num2str(laAround(i)),'FontSize',FontSize,'HorizontalAlignment','Left');
            end
        end
        % Specify the line type for the trajectories of lRT of
        % transformation parameters
        set(plot1,{'LineStyle'},LineStyle);
        
        
        if ij==cColToComp-nc+1 || ij==cColToComp
            xlabel('Subset size m');
        end
        
        ij=ij+1;
        
        ylabel(['Var: ' colnames{jj} ]);
        
        if signlr==1
            qua=sqrt(chi2inv(conflev,1));
            quant=[-qua qua];
        else
            quant=chi2inv(conflev,1);
        end
        v=axis;
        for i=1:length(quant)
            line(v(1:2),[quant(i),quant(i)] ,'color','r','LineWidth',LineWidthEnv,'Tag','env');
        end
        
        % Apply tag to current figure
        set(gcf,'Tag',tag)
        
    end
end
%% Store quantities in structure out
% Initialize out
out=struct;

% MLEtra=MLE of tramsformation parameters
out.LRT=LRT;
% Exflag = reason why maximization algorithm stopped
out.Exflag=Exflag;
% Un = Units entering the subset
out.Un=Un;
end

