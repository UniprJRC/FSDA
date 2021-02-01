function [RAW,REW, varargout] = mcdCorAna(N,varargin)
%mcdCorAna computes Minimum Covariance Determinant in correspondence analysis
%
%<a href="matlab: docsearchFS('mcdCorAna')">Link to the help function</a>
%
%
%  Required input arguments:
%
%       N    :    Contingency table (default) or n-by-2 input dataset.
%                 2D Array or Table.
%                 2D array or table which contains the input contingency
%                 table (say of size I-by-J) or the original data matrix X.
%                 In this last case N=crosstab(X(:,1),X(:,2)) or
%                 N=crosstab(X(:,1),X(:,2)) if X is in table format. As
%                 default procedure assumes that the input is a contingency
%                 table.
%               Data Types - table, or array
%
%
%  Optional input arguments:
%
%
%      bdp    : Breakdown point. Scalar. (Number between 0
%               and 0.5) or if it an integer greater than 1 bdp is the
%               number of data points which have to determine the fit
%                The default value is 0.5.
%               Example - 'bdp',1/4
%               Data Types - double
%
%      nsamp  : Number of subsamples. Scalar. Number of subsamples of size
%               J which have to be extracted (if not given, default =
%               1000).
%               Example - 'nsamp',10000
%               Data Types - double
%
%    refsteps : Number of refining iterations. Scalar. Number of refining
%               iterations in each subsample (default = 3).
%               refsteps = 0 means "raw-subsampling" without iterations.
%               Example - 'refsteps',10
%               Data Types - double
%
%     reftol  : Refining steps tolerance. Scalar. Tolerance for the refining steps.
%               The default value is 1e-6;
%               Example - 'reftol',1e-8
%               Data Types - double
%
%refstepsbestr: Number of refining iterations. Scalar. Number of refining iterations
%               for each best subset (default = 50).
%               Example - 'refstepsbestr',10
%               Data Types - double
%
% reftolbestr : Tolerance for refining steps. Scalar. Value of tolerance for the
%               refining steps for each of the best subsets.
%               The default value is 1e-8;
%               Example - 'reftolbestr',1e-8
%               Data Types - double
%
%      bestr  : Number of best solutions to store. Scalar. Number of "best locations"
%               to remember from the subsamples. These will be later iterated until
%               convergence (default=5)
%               Example - 'bestr',10
%               Data Types - double
%
%     conflev : Confidence level. Scalar. Number between 0 and 1 containing
%               confidence level which is used to declare units as outliers.
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/I, 1-0.025/I, 1-0.01/I (simultaneous alpha).
%               Default value is 0.99 per cent simultaneous
%               Example - 'conflev',0.99
%               Data Types - double
%
%
%       plots : Plot on the screen. Scalar or structure.
%               If plots is a structure or scalar equal to 1, generates:
%               (1) a plot of Mahalanobis distances against index number. The
%               confidence level used to draw the confidence bands for
%               the MD is given by the input option conflev. If conflev is
%               not specified a nominal 0.975 confidence interval will be
%               used.
%               (2) a scatter plot matrix with the outliers highlighted.
%               If plots is a structure it may contain the following fields
%                   plots.labeladd = if this option is '1', the outliers in the
%                       spm are labelled with their unit row index. The
%                       default value is labeladd='', i.e. no label is
%                       added.
%                   plots.nameY = cell array of strings containing the labels of
%                       the variables. As default value, the labels which
%                       are added are Y1, ...Yv.
%               Example - 'plots',1
%               Data Types - double or structure
%
%      Lr : row labels. Cell.
%               Cell of length n containing the labels of the rows.
%                   Example - 'label',{'UK' ...  'IT'}
%                   Data Types - cell
%
%      Lc : column labels. Cell.
%               Cell of length n containing the labels of the rows.
%                   Example - 'label',{'UK' ...  'IT'}
%                   Data Types - cell
%        msg  : Display or not messages on the screen.
%               Scalar. If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the final
%               estimator else no message is displayed on the screen.
%               Example - 'msg',1
%               Data Types - double
%
%     tolMCD  : Tolerance to declare a subset as singular. Scalar. The
%               default value of tolMCD is exp(-50*v).
%               Example - 'tolMCD',1e-20
%               Data Types - double
%
% findEmpiricalEnvelope : Empirical Confidence level. Boolean.
%               If findEmpiricalEnvelope is true (default is false) the
%               empirical envelope for each Mahalanobis distance of
%               each Profile row of the contingency table is computed, else
%               the empirical envelopes are found just if input option
%               plots=1. In case findEmpiricalEnvelope is false the
%               theoretical envelope is based on the quantiles of the
%               following scaled gamma distribution
%               chi2inv(conflev,(J-1)*(I-1)/I)/n,I,1)./r;
%               Example - 'findEmpiricalEnvelope',true
%               Data Types - Boolean
%
%
%
%  Output:
%
%  The output consists of two structures RAW and REW. RAW refers to raw
%  MCD, on the other hand, REW refers to reweighted MCD
%
%
%         RAW:   structure which contains the following fields
%
%         RAW.h    = scalar. The number of observations that have
%                    determined the MCD estimator
%         RAW.loc  = 1 x J  vector containing raw MCD location of the data
%         RAW.cov  = robust MCD estimate of
%                    covariance matrix. It is the raw MCD covariance matrix
%                    (multiplied by a finite sample correction factor and
%                    an asymptotic consistency factor).
%           RAW.obj= The determinant of the raw MCD covariance matrix.
%           RAW.bsb= k x 1 vector containing the rows of matrix N which
%                    contributed to the computation of the MCS estimate of location
%           RAW.md = I x 1 vector containing the estimates of the robust
%                    Mahalanobis distances (in squared units). This vector
%                    contains the distances of each observation from the
%                    raw MCD location of the data, relative to the raw MCD
%                    scatter matrix diag(raw MCD location)
%     RAW.outliers = A vector containing the list of the rows declared as
%                    outliers using confidence level specified in input
%                    scalar conflev
%      RAW.conflev = Confidence level that was used to declare outliers and
%                    to do reweighting.
%      RAW.singsub = Number of subsets without full rank. Notice that
%                    out.singsub > 0.1*(number of subsamples) produces a
%                    warning
%      RAW.weights = I x 1 vector containing the estimates of the weights.
%                    Weights assume values in the interval [0 1].
%                    Weight is 1 if the
%                    associated row fully contributes to compute
%                    centroid and covariance matrix. If for a particular row
%                    weight is 0.7 it means that the associated row contributes
%                    with 70 per cent of its row mass. 0 weight for a
%                    particular row it means that the associated row does
%                    not participate at all.
%                    Note that sum(N,2)'*RAW.weights=h
%            RAW.N = Original contingency table in array format.
%            RAW.Ntable = Original contingency table in table format.
%            RAW.Y = array I-by-J containing matrix of Profile Rows.
%            RAW.EmpEnv=array of size I-by-1 containing empirical envelopes for
%                   each Mahalanobis distance if input option
%                   findEmpiricalEnvelope is true or scalar containing
%                   quantile which ahse been used to declare the outliers.
%        RAW.class = 'mcdCorAna'
%
%         REW:   structure which contains the following fields
%           REW.md = I x 1 vector containing the estimates of the robust
%                    Mahalanobis distances (in squared units). This vector
%                    contains the distances of each observation from the
%                    rewighted MCD location of the data, relative to the reweighted MCD
%                    scatter matrix diag(reweighted MCD location)
%      REW.weights = I x 1 vector containing the estimates of the weights.
%                    Weights assume values 0 or 1. Weight is 0 if the
%                    associated row has been declared outlier after reweighting.
%     REW.outliers = A vector containing the list of the rows declared as
%                    outliers using confidence level specified in input
%                    scalar conflev
%            REW.Y = array I-by-J containing matrix of Profile Rows.
%            REW.EmpEnv=array of size I-by-1 containing empirical envelopes for
%                   each Mahalanobis distance if input option
%                   findEmpiricalEnvelope is true or scalar containing
%                   quantile which ahse been used to declare the outliers.
%        REW.class = 'mcdCorAna'
%
%
%  Optional Output:
%
%            C     : matrix of size nsamp-by-J which contains the indices
%                    of the subsamples extracted for
%                    computing the estimate.
%
% More About:
%
% MCDcorAna computes the MCD estimator for a contingency table.  This
% estimator is given by the subset of s Profile rows with smallest
% covariance determinant.  The MCD location estimate is then the mean of
% those h Profile points.
% The default value of h is roughly 0.5n (where n is the total number of
% observations), but the user may choose each value between n/2 and n.
%
%
% See also: mcd.m, CorAna.m
%
% References:
%
% Rousseeuw, P.J. and Van Driessen, K. 1999. A fast algorithm for the
% minimum covariance determinant estimator. Technometrics, 41:212?223.
% Greenacre, M.J. (1993), "Correspondence Analysis in Practice", London,
% Academic Press.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mcdCorAna')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % mcdCorAna with option plots=1.
  N=[134    76    43    50    49
    173    62    20    23    16
    67    76    48    36    23
    11    21    31    36    52
    25    32    57    60    58
    32    42    40    67    67
    20    35    31    41    41
    10    16    23    23    24
    54    28    29    30    23
    12    19    14    15    20
     9    10    14    20    23
    52    43    38    47    54
    21    36    33    30    36
    85    74    55    31    22
     3     8    12    12    25
    28    33    40    31    45
     9    17    23    19    34
    18    36    44    35    40
    12    24    22    25    37
    16    32    35    39    38
    28    39    36    41    54
     3    15    22    25    24
    30    40    28    20    26
     8    10    12    13    17
     2     1     2     3     3
    29    10    16     8     9
    47    51    29    19    12
     7    19    20    26     9];
    rowlab={'GB' 'SK' 'BG' 'IE' 'BE' 'ES' 'PL' 'FI' 'GR' 'HU' 'SI' 'NL' 'IT' 'RO'...
     'AT' 'FR' 'HR' 'SE' 'CZ' 'DK' 'DE' 'LT' 'PT' 'EE' 'LU' 'MT' 'LV' 'CY'};
    collab={'x1' 'x2' 'x3' 'x4' 'x5'};
    Ntable=array2table(N,'RowNames',rowlab,'VariableNames',collab);
    RAW=mcdCorAna(Ntable,'plots',1);
%}


%{
    %% mcdCorAna with bdp=0.
  N=[  69    46    41    13    22    18
    29    52    45     3     5     3
    19    55    47     2     3     1
    50    22    19     8    10     7
    25    38    33     2     4     3
    30     2     1    45     8     2
    35     6     5    32     5     2
    28    12     7     7     5     4
    26    12    11    11     4     3
    21     6     4     3     3     2];
    rowlab={'Teens' 'PicksYouUp' 'Energy' 'EnjoyLife' ...
        'WhenTired' 'Kids' 'Fun' 'Refreshes' ...
        'CheersYouUp' 'Relax'};
    collab={'Coke' 'V' 'RedBull' 'Fanta' 'Pepsi' 'DietCoke'};
    Ntable=array2table(N,'RowNames',rowlab,'VariableNames',collab);
    RAW=mcdCorAna(Ntable,'bdp',0);
    % Note that in this case RAW.md is equal to
    % out.OverviewRows.Inertia./out.OverviewRows.Mass
    % out = output from traditional correspondence analysis.
    out=CorAna(Ntable,'dispresults',false,'plots',0);
    d2=out.OverviewRows.Inertia./out.OverviewRows.Mass;
    disp('Square distance of each row profile from the centroid')
    disp([RAW.md d2])
%}

%{
    % mcdCorAna with option plots.
  N=[134    76    43    50    49
    173    62    20    23    16
    67    76    48    36    23
    11    21    31    36    52
    25    32    57    60    58
    32    42    40    67    67
    20    35    31    41    41
    10    16    23    23    24
    54    28    29    30    23
    12    19    14    15    20
     9    10    14    20    23
    52    43    38    47    54
    21    36    33    30    36
    85    74    55    31    22
     3     8    12    12    25
    28    33    40    31    45
     9    17    23    19    34
    18    36    44    35    40
    12    24    22    25    37
    16    32    35    39    38
    28    39    36    41    54
     3    15    22    25    24
    30    40    28    20    26
     8    10    12    13    17
     2     1     2     3     3
    29    10    16     8     9
    47    51    29    19    12
     7    19    20    26     9];
    rowlab={'GB' 'SK' 'BG' 'IE' 'BE' 'ES' 'PL' 'FI' 'GR' 'HU' 'SI' 'NL' 'IT' 'RO'...
     'AT' 'FR' 'HR' 'SE' 'CZ' 'DK' 'DE' 'LT' 'PT' 'EE' 'LU' 'MT' 'LV' 'CY'};
    collab={'x1' 'x2' 'x3' 'x4' 'x5'};
    Ntable=array2table(N,'RowNames',rowlab,'VariableNames',collab);
    RAW=mcdCorAna(Ntable,'plots',1);
%}

%{
    %% Raw and reweighted MCD.
    % mcdCorAna with option findEmpirical.
  N=[134    76    43    50    49
    173    62    20    23    16
    67    76    48    36    23
    11    21    31    36    52
    25    32    57    60    58
    32    42    40    67    67
    20    35    31    41    41
    10    16    23    23    24
    54    28    29    30    23
    12    19    14    15    20
     9    10    14    20    23
    52    43    38    47    54
    21    36    33    30    36
    85    74    55    31    22
     3     8    12    12    25
    28    33    40    31    45
     9    17    23    19    34
    18    36    44    35    40
    12    24    22    25    37
    16    32    35    39    38
    28    39    36    41    54
     3    15    22    25    24
    30    40    28    20    26
     8    10    12    13    17
     2     1     2     3     3
    29    10    16     8     9
    47    51    29    19    12
     7    19    20    26     9];
    rowlab={'GB' 'SK' 'BG' 'IE' 'BE' 'ES' 'PL' 'FI' 'GR' 'HU' 'SI' 'NL' 'IT' 'RO'...
     'AT' 'FR' 'HR' 'SE' 'CZ' 'DK' 'DE' 'LT' 'PT' 'EE' 'LU' 'MT' 'LV' 'CY'};
    collab={'x1' 'x2' 'x3' 'x4' 'x5'};
    Ntable=array2table(N,'RowNames',rowlab,'VariableNames',collab);
    [RAW,REW]=mcdCorAna(Ntable,'plots',1);
%}

%% Beginning of code

if ~isempty(varargin)
    UserOptions=varargin(1:2:length(varargin));
    checkdatamatrix = strcmp(UserOptions,'datamatrix');
    if sum(checkdatamatrix)
        datamatrix = varargin{2*find(checkdatamatrix)};
    else
        datamatrix=false;
    end
else
    datamatrix=false;
end

% If input is a datamatrix it is necessary to construct the contingency
% table
if datamatrix == true
    if istable(N)
        [N,~,~,labelsTab] =crosstab(N{:,1},N{:,2});
    else
        [N,~,~,labelsTab] =crosstab(N(:,1),N(:,2));
    end
    [I,J]=size(N);
    % default labels for rows of contingency table
    Lr=labelsTab(1:I,1);
    % default labels for columns of contingency table
    Lc=labelsTab(1:J,2);
    if ~verMatlab
        % Make valid names
        Lr=matlab.lang.makeValidName(Lr);
        Lc=matlab.lang.makeValidName(Lc);
    end
    Ntable=array2table(N,'RowNames',Lr,'VariableNames',Lc);
else
    [I,J]=size(N);
    if istable(N)
        Ntable=N;
        Lr = Ntable.Properties.RowNames;
        Lc = Ntable.Properties.VariableNames;
        N=N{:,:};
    else
        % default labels for rows of contingency table
        Lr=cellstr(strcat('r',num2str((1:I)')));
        % default labels for columns of contingency table
        Lc=cellstr(strcat('c',num2str((1:J)')));
        Ntable=array2table(N,'RowNames',Lr,'VariableNames',Lc);
        
    end
    
end


%grand total
n=sum(sum(N));


% P = correspondence matrix  containing relative frequencies
P = (1/n) * N;

% r= vector which contains row masses = centroids of the column profiles
r  = sum(P,2) ;

rtimesn=r*n;

% Row profiles (equation 4.14)
% Dr=diag(r);
% ProfilesRows = Dr^(-1) * P;
% Every row of matrix P is divided by the row total
ProfilesRows=P./r;

% default value of break down point
bdpdef=0.5;

% If the number of all possible subsets is <10000 the default is to extract
% all subsets otherwise just 1000.
% Notice that we use bc, a fast version of nchoosek. One may also use the
% approximation floor(exp(gammaln(n+1)-gammaln(n-p+1)-gammaln(p+1))+0.5)
ncomb=bc(I,J);
nsampdef=min(1000,ncomb);

% default value of number of refining iterations (C steps) for each
% extracted subset
refstepsdef=3;
% default value of tolerance for the refining steps convergence for  each
% extracted subset
reftoldef=1e-6;
% default value of number of best locs to remember
bestrdef=5;
% default value of number of refining iterations (C steps) for best subsets
refstepsbestrdef=50;
% default value of tolerance for the refining steps convergence for best
% subsets
reftolbestrdef=1e-8;

% Tolerance to declare a subset as singular
% It was set to
% tolMCDdef=exp(-50*v);
% but this is useless, as the roundoff level is eps = 2^(-52)
tolMCDdef=eps('double');

% Boolean which specified whether it is necessary to find the empirical
% confidence level for each row of the contingency table.
findEmpiricalEnvelope=false;
conflev=1-0.01/I;

% store default values in the structure options
options=struct('nsamp',nsampdef,'refsteps',refstepsdef,'bestr',bestrdef,...
    'reftol',reftoldef,...
    'refstepsbestr',refstepsbestrdef,'reftolbestr',reftolbestrdef,...
    'bdp',bdpdef,'plots',0,'conflev',conflev,...
    'msg',1,'tolMCD',tolMCDdef,'findEmpiricalEnvelope',findEmpiricalEnvelope,'Lr',{Lr},'Lc',{Lc});

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:mcdCorAna:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

bdp = options.bdp;              % break down point
refsteps = options.refsteps;    % refining steps
bestr = options.bestr;          % best locs for refining steps till convergence
nsamp = options.nsamp;          % subsamples to extract
reftol = options.reftol;        % tolerance for refining steps

refstepsbestr=options.refstepsbestr;  % refining steps for the best subsets
reftolbestr=options.reftolbestr;      % tolerance for refining steps for the best subsets

% tolMCD threshold under which the determinant of the covariance matrix is
% thought to be singular
tolMCD=options.tolMCD;

% If msg =1 the total estimated time to compute MCD is printed on the screen
msg=options.msg;

% Use or not empirical envelope to detect the outliers.
findEmpiricalEnvelope=options.findEmpiricalEnvelope;

%Labels of row profiles.
Lr = options.Lr;
Lc = options.Lc;

% Initialize the matrices which contain the best "bestr" estimates of
% location, indexes of subsets, cov matrices and objective function
bestlocs = zeros(bestr, J);

bestobjs = Inf * ones(bestr,1);

% singsub = scalar which will contain the number of singular subsets which
% are extracted (that is the subsets of size p which are not full rank)
singsub=0;


% hmin is the minimum number of observations whose covariance determinant
% will be minimized.
hmin=floor(2*floor((n+J+1)/2)-n+2*(n-floor((n+J+1)/2))*(0.5));

if bdp>1
    h=bdp;
else
    h=floor(2*floor((n+J+1)/2)-n+2*(n-floor((n+J+1)/2))*(1-bdp));
end


if h < hmin
    error('FSDA:mcdCorAna:Wrongh',['The MCD must cover at least ' int2str(hmin) ' observations.'])
elseif h > n
    error('FSDA:mcdCorAna:Wrongh','h is greater than the number of non-missings and non-infinites.')
end

% bestsubset is the matrix which will contain the indexes of the bestr
% subsets. Each row refers to a subset.
% Remark: note that each subset should have v+1 elements. However, due to
% the fact that if the subset is singulat we continue adding randomly
% elements to it up to when it becomes non singular, it is possible that
% certain subsets have more than v+1 elements
bestsubset = zeros(bestr, h,'int8');

% write in structure RAW the value of h
RAW=struct;
RAW.h=h;

seq=1:n;
conflev = options.conflev;

%% h==n is the case in which there is no trimming (classical case)

if h==n
    if msg
        disp(['The MCD estimates are equal to the classical estimates h=n=',num2str(h)]);
    end
    %  REW.method=char(REW.method,msgn);
    
    % cprime = row vector of column masses = centroid of row profiles
    cprime  = sum(N,1)/n;
    
    RAW.loc=cprime;
    RAW.cov=diag(cprime);
    RAW.obj=prod(cprime);
    md=mahalCorAna(ProfilesRows,cprime);
    
    RAW.md=md;
    RAW.class='mcdCorAna';
    bsbh=1:I;
else
    
    %% Extract in the rows of matrix C the indexes of all required subsets
    [C,nselected] = subsets(nsamp,I,J,ncomb,1);
    % Store the indices in varargout
    if nargout==3
        varargout={C};
    end
    % initialise and start timer.
    tsampling = ceil(min(nselected/100 , 1000));
    time=zeros(tsampling,1);
    
    % ij is a scalar used to ensure that the best first bestr non singular
    % subsets are stored
    ij=1;
    for i = 1:nselected
        if i <= tsampling, tic; end
        
        % extract a subset of size v+1 .
        index = C(i,:);
        
        % size Yj is J-by-J
        Nj = N(index,:);
        
        %grand total
        nred=sum(Nj,'all');
        
        % Column masses = centroids of the row profiles using subset
        cj  = sum(Nj,1)/nred;
        
        % Find weighted estimate of the mean and cov
        locj = cj;        % centroid of subset
        Sj = diag(cj);    % covariance of subset
        
        % Check if the subset is in general position (rank<v)
        if det(Sj)< tolMCD
            singsub = singsub + 1;
        end
        
        
        % Function IRWLSmcd performs refsteps concentration steps of IRLS on elemental
        % start. Input:
        % - N = contingency table of dimension I-by-J
        % - ProfilesRows = matrix of profile rows
        % - rtimes = sum(N,1) = vector of row masses multiplied by n
        % - locj = row vector containing (robust) centroid
        % - refsteps = number of refining iterations
        % - reftol = tolerance for convergence of refining iterations
        % outIRWLS = IRWLSmcd(N,ProfilesRows, rtimesn, locj, Sj, h, refsteps, reftol);
        outIRWLS = IRWLSmcd(N,ProfilesRows, rtimesn, locj, h, refsteps, reftol);
        
        % If the value of the objective function is smaller than tolMCD
        % we have a perfect fit situation, that is there are h observations
        % that lie on the hyperplane.
        if outIRWLS.obj < tolMCD
            return
        end
        
        % The output of IRWLSmult is a structure containing centroid, cov
        % matrix and estimate value of the objective function (which has
        % been minimized, that is |cov|
        locrw = outIRWLS.loc;
        objrw = outIRWLS.obj;
        
        
        % to find s, save first the best bestr scales and shape matrices
        % (deriving from non singular subsets) and, from iteration bestr+1
        % (associated to another non singular subset), replace the worst scale
        % with a better one as follows
        if ij > bestr
            % from the second step check whether new loc and new shape belong
            % to the top best loc; if so keep loc and shape with
            % corresponding scale.
            
            
            if  objrw < max(bestobjs)
                % Find position of the maximum value of bestscale
                [~,ind] = max(bestobjs);
                bestobjs(ind) = objrw;
                bestlocs(ind,:) = locrw;
                % best subset associated with minimum value
                % of the objective function
                bestsubset(ind,1:length(index))=index;
            else
                
            end
        else
            bestobjs(ij) = objrw;
            bestlocs(ij,:) = locrw;
            bestsubset(ij,1:length(index)) = index;
            ij=ij+1;
        end
        
        
        % Write total estimation time to compute final estimate
        if i <= tsampling
            
            % sampling time until step tsampling
            time(i)=toc;
        elseif i==tsampling+1
            % stop sampling and print the estimated time
            if msg==1
                fprintf('Total estimated time to complete MCD: %5.2f seconds \n', nselected*median(time));
            end
        end
        
        
    end
    if singsub==nselected
        error('FSDA:mcdCorAna:NoFullRank','No subset had full rank. Please increase the number of subsets or check your design matrix X')
    end
    
    if singsub/nselected>0.1
        disp('------------------------------')
        disp(['Warning: Number of subsets without full rank equal to ' num2str(100*singsub/nsamp) '%'])
    end
    
    % perform C-steps on best 'bestr' solutions, till convergence or for a
    % maximum of refstepsbestr steps using a convergence tolerance as specified
    % by scalar reftolbestr
    
    % this is to ensure that the condition tmp.scale < superbestscale in the
    % next if statement is satisfied at least once
    superbestobj = Inf;
    for i=1:bestr
        tmp = IRWLSmcd(N,ProfilesRows, rtimesn, bestlocs(i,:), h, refstepsbestr, reftolbestr);
        
        if tmp.obj < superbestobj
            superbestobj    = tmp.obj;
            superbestloc    = tmp.loc;
            bsbh=tmp.bsb;
            % weights = tmp.weights;
        end
    end
    
    RAW.class   = 'mcdCorAna';
    RAW.obj   = superbestobj;       % value of the objective function
    
    % RAW.bsb is the list of the rows of contingency table used to compute best centroid
    RAW.bsb=bsbh;
    
    % RAW.bsbh is the list of the rows of contingency table which
    % determined the final fit to obtain the best h units.
    % RAW.bsbh=bsbh;
    
    RAW.loc= superbestloc;
    RAW.cov = diag(superbestloc);
    RAW.obj = superbestobj;
    
    % Mahalanobis distances (in squared units) on Profile Rows matrix
    % md=mahalFS(ProfilesRows,superbestloc,diag(superbestloc));
    md=mahalCorAna(ProfilesRows,superbestloc);
    
    % Store vector of Mahalanobis distances (in squared units)
    RAW.md = md;
    
    
end


%  Store confidence level
RAW.conflev=conflev;

% Store total number of singular subsets
RAW.singsub=singsub;
RAW.N=N;
RAW.Ntable=Ntable;

% weights of each row of the contingency table for the calculation of the
% centroid.
weights=zeros(I,1);
weights(bsbh(1:end-1))=1;
lastunit=bsbh(end);
weights(lastunit)=(h-sum(N(bsbh(1:end-1),:),'all'))/sum(N(lastunit,:));
RAW.weights=weights;

if findEmpiricalEnvelope == true
    % nrowt = column vector containing row marginal totals
    nrowt=sum(N,2);
    % ncolt = row vector containing column marginal totals
    ncolt=sum(N,1);
    disp('Finding empirical bands')
    if conflev<=0.99
        nsimul=200;
    elseif conflev <=0.999
        nsimul=1000;
    else
        disp('Empirical band is very extreme: running 5000 replicates')
        nsimul=5000;
    end
    
    mmdStore=zeros(I,nsimul);
    if nargout>1
        mmdStoreR=zeros(I,nsimul);
    end
    if nargout>1
        alsorew=true;
    else
        alsorew=false;
    end
    parfor j=1:nsimul
        % Generate random contingency table using row and column marginals
        % of the current table
        out1=rcontFS(I,J,nrowt,ncolt,'nocheck',true,'algorithm','144');
        Nsim=out1.m144;
        % Nsim=out1.m159;
        if alsorew==false
            TMP=mcdCorAna(Nsim,'plots',0,'msg',0,'bdp',h);
            mmdStore(:,j)=TMP.md;
        else
            [TMP,TMPr]=mcdCorAna(Nsim,'plots',0,'msg',0,'bdp',h);
            mmdStore(:,j)=TMP.md;
            mmdStoreR(:,j)=TMPr.md;
        end
    end
    
    % Sort rows of matrix mmdStore (MD in squared units)
    mmdStore=sort(mmdStore,2);
    EmpEnv=mmdStore(:,ceil(nsimul*conflev));
    if alsorew==true
        mmdStoreR=sort(mmdStoreR,2);
        EmpEnvR=mmdStoreR(:,ceil(nsimul*conflev));
    end
    
    
else
    EmpEnv=repmat(chi2inv(conflev,(J-1)*(I-1)/I)/n,I,1)./r;
    EmpEnvR=EmpEnv;
end

weightsboo=md <= EmpEnv;
RAW.outliers=seq(~weightsboo);
% Matrix Profile Rows is stored in stucutre RAW
RAW.Y=ProfilesRows;
RAW.EmpEnv=EmpEnv;

plo=options.plots;

% Plot Mahalanobis distances with outliers highlighted
if isstruct(plo) || (~isstruct(plo) && plo~=0)
    if h<n
        laby='Raw MCD squared Mahalanobis distances';
    else
        laby='Squared Mahalanobis distances';
    end
    malindexplot(RAW.md,EmpEnv,'conflev',conflev,'laby',laby,...
        'numlab',RAW.outliers,'tag','rawmcd','label',Lr);
    
    figure('Tag','pl_spm_outliers');
    group=ones(I,1);
    if ~isempty(RAW.outliers)
        group(RAW.outliers)=2;
    end
    spmplot(ProfilesRows,group,plo);
    set(gcf,'Name',' Raw MCD: scatter plot matrix of row profiles with outliers highlighted');
end


% Reweighting part

if nargout>1
    % size Yj is J-by-J
    Nrew = N(weightsboo,:);
    
    %grand total
    nred=sum(Nrew,'all');
    
    % Column masses = centroids of the row profiles using subset
    % crew= final centroid after reweighting
    crew  = sum(Nrew,1)/nred;
    
    % Compute Mahalanobis distances using locrew
    md=mahalCorAna(ProfilesRows,crew);
    
    REW=struct;
    % Store vector of Mahalanobis distances (in squared units)
    REW.md = md;
    REW.weights=md <= EmpEnv;
    
    REW.outliers=seq(~REW.weights);
    % Matrix Profile Rows is stored in stucutre RAW
    REW.Y=ProfilesRows;
    REW.EmpEnv=EmpEnvR;
    REW.class   = 'mcdCorAna';
    
end

% Plot Mahalanobis distances with outliers highlighted
if (isstruct(plo) || (~isstruct(plo) && plo~=0)) && nargout>1
    
    laby='Reweighted squared Mahalanobis distances';
    malindexplot(REW.md,EmpEnvR,'conflev',conflev,'laby',laby,...
        'numlab',REW.outliers,'tag','rewmcd','label',Lr);
    
    figure('Tag','pl_spm_outliers(REW)');
    group=ones(I,1);
    if ~isempty(REW.outliers)
        group(REW.outliers)=2;
    end
    spmplot(ProfilesRows,group,plo);
    set(gcf,'Name',' Rew MCD: scatter plot matrix of profile rows with outliers highlighted');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Subfunctions called by main function mcd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% subfunction IRWLSmcd
    function outIRWLS = IRWLSmcd(N, ProfileRows, r, initialloc, h, refsteps, reftol)
        %IRWLSmcd (iterative reweighted least squares) does refsteps refining steps from initialloc
        % for refsteps times or till convergence.
        %
        %  Required input arguments:
        %
        % - N = contingency table of dimension I-by-J
        % - ProfilesRows = matrix of profile rows
        % - r = sum(N,1) = vector of row masses multiplied by n
        % - initialloc = row vector containing (robust) centroid
        % - h = number of observations to use
        % - refsteps  = scalar, number of refining (IRLS) steps
        % - reftol   = relative convergence tolerance for the fully iterated
        %               best candidates. Deafult value is 1e-7
        %
        %  Optional input arguments:
        %
        %
        %  Output:
        %
        %  The output consists of a structure 'outIRWLS' containing:
        %      outIRWLS.loc     : v x 1 vector. Estimate of location after refsteps
        %                         refining steps.
        %      outIRWLS.obj     : scalar. Value of the objective function after refsteps refining
        %                         steps.
        %      outIRWLS.bsb     : vector. Indexes of the rows of
        %                         contingency table which contributed to the
        %                         calculation of loc. If the length of
        %                         bsb is q it means that the rows of
        %                         the contingency table bsb(1:m-1)
        %                         contributed to loc with weight 1, and row bsb(m)
        %                         contributed with a fraction
        %                         (h-sum(N(bsb(1:q-1),:),'all')/sum(N(bsb(q),:))
        %                         of its mass
        %
        
        loc = initialloc;
        % Mahalanobis distances (in squared units) from initialloc and Initialshape
        % (multiplied by masses)
        mahaldist2=r.*mahalCorAna(ProfileRows,initialloc);
        
        iter = 0;
        locdiff = 9999;
        
        while ( (locdiff > reftol) && (iter < refsteps) )
            iter = iter + 1;
            
            [~,sortdist]=sort(mahaldist2);
            
            % The rows sortdist(1:indexesCR) will be completely
            % represented;
            cumsumnjdot=cumsum(r(sortdist));
            indexesCR=find(cumsumnjdot<h,1,'last');
            
            if isempty(indexesCR)
                indexesCR=0;
                unitstoADD=h;
            else
                % Find how many units must be included from
                % row indexesCR+1 of the initial contingency table;
                unitstoADD=h-cumsumnjdot(indexesCR);
            end
            
            bsb=sortdist(1:indexesCR+1);
            Niter=N(bsb,:);
            Niter(end,:)=unitstoADD*Niter(end,:)/r(sortdist(indexesCR+1));
            newloc      = sum(Niter,1)/h;
            
            % Compute MD
            mahaldist2=r.*mahalCorAna(ProfileRows,newloc);
            % mahaldist2=mahalCorAna(ProfileRows,newloc);
            
            
            % locdiff is linked to the tolerance
            locdiff = norm(newloc-loc,1)/norm(loc,1);
            loc = newloc;
            
        end
        obj         = prod(newloc);
        outIRWLS = struct('loc',newloc,'obj',obj,'bsb',bsb);
    end

end
%FScategory:MULT-Multivariate
