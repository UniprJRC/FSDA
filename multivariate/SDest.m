function [out, varargout] = SDest(Y,varargin)
%SDest computes Stahel-Donoho robust estimator of dispersion/location
%
%<a href="matlab: docsearchFS('sdest')">Link to the help function</a>
%
% Required input arguments:
%
%    Y: Input data. Matrix. Data matrix containing n observations on v variables
%       $Y=(y_1^T, ..., y_i^T, ..., y_n^T)^T$. 
%       $y_i^T$ of size 1-by-v is the ith row of matrix Y.
%       Rows of Y represent observations, and columns
%       represent variables.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       automatically be excluded from the computations.
%
% Optional input arguments:
%
%      nsamp   : Number of random directions. Scalar. 
%                Scalar defining number of random directions which have
%                to be extracted (if not given, default = 1000)
%                 Example - 'nsamp',1000 
%                 Data Types - single | double
%      jpcorr  : Subsamples additional size. Scalar. 
%                Integer value greater or equal 0. 
%                If jpcorr=0 subsamples of size v are extracted. Each
%                subsample gives rise to only 1 direction (called dir). 
%                If jpcorr=1 subsamples of size v+1 are extracted and
%                they give rise to (v+1) directions. 
%                If jpcorr = 2, 3, 4, ... subsamples of size v+jpcorr are
%                extracted. We remove the jpcorr-1 units from the subset
%                which have the largest jpcorr-1 Mahalanobis distances so
%                we obtain a subsample of v+1 units which, as before, gives
%                rise to (v+1) directions. 
%                In the particular case when
%                jpcorr=2 we end up with Juan and Prieto suggestion, that
%                is the unit which has the largest MD is removed from each
%                subsample. The default value of jpcorr is 0.
%                 Example - 'jpcorr',1 
%                 Data Types - single | double
%     conflev : Confidence level which is
%               used to declare units as outliers. Scalar. 
%               Scalar between 0 and 1. 
%               Usually conflev=0.95, 0.975 0.99 (individual alpha)
%               or 1-0.05/n, 1-0.025/n, 1-0.01/n (simultaneous alpha).
%               Default value is 0.975. 
%                 Example - 'conflev',0.95 
%                 Data Types - single | double
%      margin : Marginal projections. Scalar. 
%               Scalar which specifies if it is necessary to consider
%               marginal projections. Scalar margin specifies up to which
%               order marginal projections must be considered. For example
%               margin = 2 implies that all possible univariate and
%               bivariate marginal projections are considered. Default
%               value is 0 (marginal projections are not considered).
%                 Example - 'margin',2 
%                 Data Types - single | double
%               Remark: note that if margin>0, data are
%               preliminary standardized using medians and MADs.
%      weight : Value to use in the weight function. String. 
%               Value to use in the weight function to transform
%               the outlyingness measure of each observation into a weight.
%               Possible values are 'huber', 'tukey', 'zch', and 'mcd'.
%               Default is 'mcd'.
%               If weight = 'mcd', the [n/2] observations with the smallest
%               distance get weight 1. An asymptotic consistency factor is
%               applied to the estimated covariance matrix.
%               If weight = huber, the weights are determined according to
%               the following formula
%               $w(r) = \min(1, (r/c)^{-q})$ 
%               with: 
%               $c=\min(\sqrt{\chi^2_{v,0.5}},4)$ 
%               in Van Aelst, Vandervieren and
%               Willems, "Stahel-Donoho Estimators with Cellwise Weights",
%               (J STAT COMPUT SIM, 2011) (option c='hdim', see below); 
%               $c=\sqrt{\chi^2_{v,0.95}}$ in Todorov and Filzmoser, "An Object
%               Oriented Framework for Robust Multivariate Analysis", (J OF
%               STAT SOFTWARE, 2009) (option c='sdim', see below);
%               q = scalar see below.
%               If weight = 'tukey' the Tukey Biweight function is applied,
%               where weights are given by:
%               \[
%                   w(r)= 
%                   \left\{
%                   \begin{array}{c}
%                    [1-(r/c)^2]^2 \qquad if \qquad |r| \leq c \\
%                     0 \qquad \mbox{otherwise} 
%                   \end{array}
%                  \right.
%                 \]
%               with constant $c$ computed to obtain a prefixed nominal
%               breakdown point (bdp). 
%               If weight = 'zch', weights are computed according to Zuo,
%               Cui and He's formula (Zuo, Cui and He, "ON THE STAHEL-DONOHO ESTIMATOR AND DEPTH-WEIGHTED
%               MEANS OF MULTIVARIATE DATA", Annals of Statistics (2004)):
%               
%                 \[
%                   w(PD)= 
%                   \left\{
%                   \begin{array}{c}
%                     (\exp\{-K(1-PD/C)^2\} - \exp\{-K\})/(1-\exp\{-K\}) \qquad\qquad if \qquad PD < c \\
%                     1 \qquad \mbox{otherwise} 
%                   \end{array}
%                  \right.
%                 \]
%
%               where:
%               - $PD$ is the Projection Depth: $PD = 1/(1+r)$
%               ($r$=outlyingness measure);
%               - $C=Median(PD)$;
%               - $K$ is a positive tuning parameter.
%                 Example - 'weight','zch' 
%                 Data Types - char
%            q: Constant to be used in the Huber weight function. Scalar. 
%               The default value of q is 2 (see Maronna and Yohai, 1995).
%                 Example - 'q',2
%                 Data Types - single | double
%            c: Scale parameter. String. If c='hdim' (high dimensions) the scale parameter c in the
%               Huber weight function is given by:
%               $c= min(\sqrt{\chi^2_{v,0.5}},4)$. If c='sdim' (small
%               dimensions), parameter c is given by:
%               $c=\sqrt{\chi^2_{v,0.95}}$. The default is 'hdim'.
%                 Example - 'c','hdim' 
%                 Data Types - char
%          nbp: Nominal breakdown point. Scalar. Nominal breakdown point to be fixed in the Tukey
%               biweight function to obtain the thresold value c (0<nbp<1). The default
%               value of npb is 0.5.
%                 Example - 'nbp',0.6
%                 Data Types - single | double
%            K: Constant to be used in Zuo, Cui and He's family of weights. Scalar. 
%               The default of K is 3.
%                 Example - 'K',3
%                 Data Types - single | double
%      nocheck: Check input arguments. Scalar. If nocheck is equal to 1 no check is performed on
%               matrix Y. As default nocheck=0.
%                 Example - 'nocheck',1
%                 Data Types - single | double
%       plots : Plot on the screen. Scalar or structure.
%               If plots is a structure or scalar equal to 1, generates: 
%               (1) a plot of robust Mahalanobis distances against index number. The 
%               confidence level used to draw the confidence bands for
%               the MD is given by the input option conflev. If conflev is
%               not specified a nominal 0.975 confidence interval will be
%               used. 
%               (2) a scatter plot matrix with the outliers highlighted. 
%               (3) a scatter plot of robust Mahalanobis distances against observation weights (i.e. the
%               outlyingness measure transformed according to the weight
%               function). 
%               If plots is a structure it may contain the following
%               fields: 
%                   plots.labeladd = if this option is '1', the outliers in the
%                       spm are labelled with their unit row index. The
%                       default value is labeladd='', i.e. no label is
%                       added.
%                   plots.nameY = cell array of strings containing the labels of
%                       the variables. As default value, the labels which are
%                       added are Y1, ...Yv.
%                 Example - 'plots',1
%                 Data Types - single | double
%        msg  : Level of output to display. Scalar. scalar which controls whether to display or not messages
%               on the screen If msg==1 (default) messages are displayed
%               on the screen about estimated time to compute the final estimator
%               else no message is displayed on the screen
%                 Example - 'msg',1
%                 Data Types - single | double
%       ysave : save input matrix Y. Scalar. Scalar that is set to 1 to request that the data matrix Y
%               is saved into the output structure out. This feature is
%               meant at simplifying the use of function malindexplot.
%               Default is 0, i.e. no saving is done.
%                 Example - 'ysave',1
%                 Data Types - single | double
%      dirsave: save directions. Scalar. scalar that is set to 1 to request that the all directions
%               for all extracted subsets are saved. If dirsave=1 out
%               structure will contain a field named Dir
%                 Example - 'dirsave',1
%                 Data Types - single | double
%  rstprojsave: save robust standardized projection scores. Scalar. 
%               Scalar that is set to 1 to request that the robust
%               standardized projection scores associated to each direction
%               are saved for each subset. If projsave=1 out structure will
%               contain a field named RstProj
%                 Example - 'rstprojsave',1
%                 Data Types - single | double
%      projloc: Type of location. String. String with possible values 'median' (default) and 'mean'
%               This option controls the type of location  (robust
%               estimator of scale) to use for the projections for
%               each subset. The projections are defined as $d^T \times y_i$
%               where d is a v-by-1 vector containing a particular direction
%              ($d^T$ is its transpose) (to make estimator location invariant).
%                 Example - 'projloc',1
%                 Data Types - single | double
%    projscale: Type of standardization. String. string with possible values
%               'mad' (default), 'sn', 'qn' and 'std'.
%               This option controls the type of standardization  (robust
%               estimator of scale) to use for the centered projections for
%               each subset (to make estimator scale invariant).
%               'mad' uses median absolute deviations from the medians (see
%               file mad.m of statistics toolbox for further details)
%               'sn' uses a robust version of Gini's average difference
%               (see file Sn.m of FSDA toolbox for further details)
%               'qn' uses first quartile of interpoint distances
%               $|x_i-x_j|$
%               (see file Qn.m of FSDA toolbox for further details)
%               'std' uses non robust standard deviations
%               (see file std.m of statistics toolbox for further details)
%               The two estimators 'sn' and 'qn' have been introduced by
%               Rousseeuw and Croux (1993), JASA, as alternatives to MAD.
%                 Example - 'projscale',1
%                 Data Types - single | double
%
% Output:
%
%  out :     A structure containing the following fields
%
%     out.maxdir   = n x v matrix which contains the direction maximizing
%                    the robust standardized projection scores for each
%                    unit of the sample (that is the direction which
%                    produces the outlyingness measure for each
%                    observation)
%         out.loc  = 1 x v  vector containing SD estimate of location
%         out.cov  = v x v matrix containing robust estimate of the
%                    covariance matrix
%          out.md  = n x 1 vector containing the estimates of the robust
%                    Mahalanobis distances (in squared units)
%     out.outliers = A vector containing the list of the units declared as
%                    outliers using confidence level specified in input
%                    scalar conflev
%      out.conflev = scalar, confidence level that was used to declare outliers
%      out.weights = n x 1 vector containing the estimates of the weights
%                    Weights assume values between 0 and 1.
%                    If input string weight is 'mcd' the weights are 0 or
%                    1. More precisely the [n/2] associated with the
%                    smallest [n/2] measure of outlyingness get weight 1
%            out.Y = Data matrix Y. The field is present if option ysave
%                    is set to 1.
%          out.Dir = nsamp-by-v matrix or nsamp-by-v-by-(v+1) array which
%                    contains for each subset the direction vectors. More in
%                    details, if jpcorr=0 Dirsave is a nsamp-by-v matrix else
%                    Dirsave is a 3D array of dimension nsamp-by-v-by-(v+1).
%                    Remember that in this last case v+1 directions are
%                    considered for each subset. The field is present only
%                    if option dirsave is set to 1.
%      out.RstProj = n-by-nsamp matrix or n-by-nsamp-by-(v+1) 3D array
%                    which contains the robust standardized projection
%                    scores for each unit of the sample for the nsamp
%                    projections (or for the nsamp*(v+1) projections if
%                    jpcorr>0). The field is present only if option
%                    rstprojsave is set to 1.
%        out.class = 'SD'
%
%  Optional Output:
%
%      C     : nsamp-by-v+jpcorr matrix containing the indices of the
%             subsamples extracted for computing the SD estimate.
%             jpcorr=0,2, 3, ... (see input option jpcorr)
%
%
% More About:
%               A "robust standardized" projection
%               score along direction vector d is defined as follows: 
%
%               $Rstproj_i = \frac{d^T \times y_i -med_j(d^T \times y_j)}{MAD_j(d^T \times y_j)} \;\;\; i=1,\cdots,n$; 
%
%              where $med_j(d^T \times y_j)$ and $MAD_j(d^T \times y_j)$ are respectively
%              the median and the modified MAD  $j=1,2,\cdots,n$.
%              With our two input options  projloc and projscale it is
%              possible to use alternative estimators of location and scale
%              to standardize $d^T \times y_i$: 
%
%               $Rstproj_i = \frac{d^T \times y_i -projloc(d^T \times y_j)}{projscale(d^T \times y_j)} \;\;\; i=1,\cdots,n$; 
%
%              The outlying measure for unit i ($outl_i$) is defined as: 
%
%              $outl_i = sup_{d \in R^v} (Rstproj_i) \;\;\; i=1, \cdots, n$; 
%
%              This outlyingness measure is based on the idea that for any
%              multivariate outlier, one can always find a one-dimensional
%              projection for which the observation is a univariate
%              outlier. 
%              Remark: note that outl_i(d) = outl_i(c*d) where c
%              is a positive scalar therefore it is not necessary to
%              rescale d to unit norm.
%
% See also: MCD, Smult, MMmult, FSM
%
%
% References:
%
% Stahel, W.A. (1981), "Breakdown of covariance estimators", Research
%        Report 31, Fachgruppe für Statistik, E.T.H. Zürich, Switzerland.
% Donoho, D.L. (1982), "Breakdown Properties of Multivariate Location
%        Estimators", Ph.D. dissertation, Harvard University.
% Maronna, R.A. and Yohai, V.J. (1995), "The behavior of the Stahel-Donoho
%         robust multivariate estimator", Journal of the American
%         Statistical Association, 90, 329-341.
% Juan J. and Prieto F.J. (1995) Journal of Computational and Graphical
%         Statistics, 4, 319-334.
% Maronna, R.A., Martin D. and Yohai V.J. (2006),Robust Statistics, Theory
%         and Methods, Wiley,New York.
%
%
% Acknowledgements: 
%
% This function follows the lines of MATLAB code developed during the
% years by many authors.
%
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('sdest')">Link to the help function</a>
% Last modified 06-Feb-2015
%

% Examples

%{
    % SDest with all default options.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out]=SDest(Ycont);
%}

%{
    %% SDest with optional arguments.
    % SDest with v+1 directions for each subsample (jpcorr=1).
    % Produce plot of robust Mahalanobis distances
    n=50;
    v=5;
    randn('state', 1256);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+4;
    [out]=SDest(Ycont,'jpcorr',1,'plots',1);
%}

%{
    % SDest with exctracted subsamples.
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+3;
    [out,C]=SDest(Ycont);
%}

%{
    %% SDest jpcorr equal to 1.
    % v+1 directions for each subsample. Produce plot of robust Mahalanobis distances.
    n=50;
    v=5;
    randn('state', 1256);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(2:20,3)=100;
    [out]=SDest(Ycont,'jpcorr',1,'margin',3,'plots',1);
%}

%{
    % SDest with with Juan and Prieto adjustment.
    % jpcorr>1. Subsamples of size equal to v+5 are initially drawn (jpcorr=5).
    % In each of them, the 4 units (=jpcorr-1) having the four
    % largest MDs are then discarded from the subsample.
    % For each subsample v+1 directions are computed.
    % Produce plot of robust Mahalanobis distances.
    n=50;
    v=5;
    randn('state', 1256);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(2:20,3)=5;
    [out]=SDest(Ycont,'jpcorr',5,'plots',1,'nsamp',100000);

    %  SDest with directions and robust standardized projection scores saved
    [out]=SDest(Ycont,'jpcorr',0,'plots',1,'nsamp',1000,'dirsave',1,'rstprojsave',1);

    % Compare the output of SD with the one produced by the fwd search
    outFS=FSM(Ycont,'init',20);

    % SDest with different choices for the weight functions.
    n=100;
    v=5;
    randn('state', 1256);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    rndind = randsample(n,10);
    Ycont(rndind,1)=20;
    [outhub]=SDest(Ycont,'jpcorr',2, 'weight','huber', 'plots',1);
    [outtuk]=SDest(Ycont,'jpcorr',2, 'weight','tukey', 'plots',1);
    [outtuk025]=SDest(Ycont,'jpcorr',2, 'weight','tukey', 'nbp',0.25, 'plots',1);
    [outzch]=SDest(Ycont,'jpcorr',2, 'weight','zch', 'plots',1);
    [outzchK50]=SDest(Ycont,'jpcorr',2, 'weight','zch', 'K',50, 'plots',1);

    % SDest with alternatives to the MAD.
    n=50;
    v=5;
    randn('state', 1256);
    Y=randn(n,v);
    [outMAD]=SDest(Y,'jpcorr',5,'plots',1,'nsamp',10000);
    [outSN]=SDest(Y,'jpcorr',5,'plots',1,'nsamp',10000, 'projscale', 'sn');
    [outQN]=SDest(Y,'jpcorr',5,'plots',1,'nsamp',10000, 'projscale', 'qn');
%}


%% Beginning of code
[n,v]=size(Y);

% default values of subsamples to extract
ncomb=bc(n,v);
nsampdef=min(10000,ncomb);

% store default values in the structure options
options=struct('nsamp',nsampdef,'weight','mcd',...
    'jpcorr',0,'q',2,'c','hdim','nbp',0.5,'K',3,'margin',0,'plots',0,...
    'conflev',0.975,'nocheck',0,'msg',1,'ysave',0,'dirsave',0,'rstprojsave',0,...
    'projloc','median','projscale','mad');

% check user options and update structure options
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:SDest:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

nsamp = options.nsamp;          % subsamples to extract
msg=options.msg;                %

% If margin >0 data are preliminary standardized using medians and MADs
margin=options.margin;

if margin>0
    med=median(Y);
    Ytilde = bsxfun(@minus,Y, med);
    mads=median(abs(Ytilde));
    Y = bsxfun(@rdivide, Ytilde, mads);
end

%% Extract in the rows of matrix C the indexes of all required subsets
% jpcorr Juan and Prieto correction
jpcorr=options.jpcorr;
if jpcorr-floor(jpcorr)>0
    error('FSDA:SDest:WrongJpcorr','jpcorr must be an integer value [0, 1 2 ....]')
end

if jpcorr==0
    [C,nselected] = subsets(nsamp,n,v,ncomb,msg);
    % elseif jpcorr==1
    %     ncomb=bc(n,v+jpcorr);
    %     [C,nselected] = subsets(nsamp,n,v+1,ncomb,msg);
    %     eyev=eye(v);
else % jpcorr==1 jpcorr =2 or jpcorr =3 or jpcorr =4 or jpcorr =5 ...
    ncomb=bc(n,v+jpcorr);
    [C,nselected] = subsets(nsamp,n,v+jpcorr,ncomb,msg);
    eyev=eye(v);
end

% dirsave = option which specifies if it is necessary to store all
% directions
dirsave=options.dirsave;
if dirsave==1
    if jpcorr==0
        % Dir=matrix(#subsamples-by-#vars) - only one direction
        Dir=zeros(nsamp,v);
    else
        % Dir=3D-array(#subsamples-by-#vars)-by-#directions
        Dir=zeros(nsamp,v,v+1);
    end
end

% rstprojsave = option which specifies if it is necessary to store all
% robust standardized projection measures
rstprojsave=options.rstprojsave;
if rstprojsave ==1
    % RstProj: n-by-nsamp matrix or n-by-nsamp-by-(v+1)
    if jpcorr==0
        RstProj=zeros(n,nsamp);
    else
        RstProj=zeros(n,nsamp,v+1);
    end
end

% Store the indices of the subsets in varargout
if nargout==2
    varargout={C};
end

% initialise and start timer.
tsampling = ceil(min(nselected/100 , 1000));
time=zeros(tsampling,1);

% Initialize the matrices which contain estimates of
% location, index of subsets, shape matrices and scales

projloc=options.projloc;
projscale=options.projscale;
if strcmp(projscale,'mad')
    % n1 and n2 ordered positions which will be necessary to compute modified
    % MAD
    n1 = ceil((n+v-1)/2);
    n2 = floor((n+v-1)/2)+1;
    nb = (n+v-1)/2;
    % beta = constant necessary to rescale the modified MAD
    beta = norminv(0.5*(nb/n+1),0,1);
end
% outlvec = vector which will contain the maximum of the robust standardized
% projection score for each observation, that is the outlying measure for
% each observations
outlvec = zeros(n,1);

% maxdir = n x v matrix which will contain the direction maximizing the
% robust standardized projection score (that is the outlyingness measure)
% for each observation
maxdir = zeros(n,v);

% onev1 will be used to find the orthogonal direction to the space spanned
% by subset
onev1=ones(v,1);
% half will be used to find the median
half = floor(n/2);

for i = 1:nselected
    
    if i <= tsampling, tic; end
    
    % find a subset of size v
    index = C(i,:);
    Yj = Y(index,:);
    
    if jpcorr==0;
        % The most efficient way to find the vector associated to the direction
        % which is orthogonal to the vectors which form matrix subs
        % where subs is
        % 0 ... 0
        % (y_2-y_1)'     y_2' = second row of matrix Yj
        % ...
        % (y_v-y_1)'     y_v' = vth row of matrix Yj
        % is the following
        % dir=Yj\ones(v,1);
        dir= Yj\onev1;   % dir = inv(Yj)*ones(v,1)
        
        % Remark: dir satisfies the following system of equations
        %      Yj*dir =k *ones(v,1) or without loss of generality
        %      Yj*dir = ones(v,1)
        %      or
        %      subs*dir = zeros(v,1)
        %      or
        %     (y_1-y_z)' * dir = r
        %     (y_2-y_z)' * dir = r
        %     ...
        %     (y_v-y_z)' * dir = r
        %     where y_z is any vector in R_v
        %     and r is a scalar
        %
        
        if dirsave == 1
            % Remark: it is not necessary to normalize dir
            dir=dir./sqrt(dir'*dir);
            % Store the direction associated to the subset
            Dir(i,:)=dir;
        end
        
        % An alternative (but very inefficient) way to find dir is as follows)
        %{
                % subs=
                % \left
                % (y_2-y_1)'
                % ...
                % (y_v-y_1)'
                % \right
                subs = bsxfun(@minus,Yj(2:v,:), Yj(1,:));
        
                % base = v-by-(v-1) matrix which defines a orthonormal basis
                base = orth(subs');
                % I -base*base'=mproj = idempotent projection matrix
                %
                % Each column of matrix mproj
                % produces the vector associated to the
                % direction orthogonal to the space spanned by base
                % Remark: the columns of mproj (once normalized to unit norm)
                % are all equal
        
                eyev=eye(v);
                mproj = eyev-base*base';
                % Select the column of matrix mproj which has the largest
                % norm
                % Remark: one could take any column of matrix mproj,
                % but in order to avoid numerical erros it is better to
                % take the one associated with % the largest norm
                hvec = sum(mproj.^2);
        
                % Remark: the columns of mproj (once normalized to unit norm)
                % are all equal  (apart from the sign)
                [~,ind]=max(hvec);
        
                dir = mproj(:,ind);
                norm = sqrt(dir'*dir);
                dir = dir/norm;
                % It is possible to verify this using the code below
                % all columns of matrix mproj1 are equal  (apart from the sign)
                mproj1=mproj;
                stand1=sqrt(diag(mproj'*mproj));
                for ii=1:v
                    mproj1(:,ii)=mproj(:,ii)./stand1(ii);
                end
        %}
        
        % Note that subs*dir  =zeros(v-1,1) where
        % subs = bsxfun(@minus,Yj(2:v,:), Yj(1,:));
        % and Yj*dir = ones(v,1)
        
        % vector containing the scores (using dir projection)
        Q = Y*dir;
        
        % Every observation gets a measure of outlyingness r
        % r = sup_a |Y*dir -loc(Y*dir)|/scale(Y*dir)
        % r = sup_a |Q -loc(Q) |/ scale(Q)
        % The default for "loc" and "scale" are median and modified mad respectively
        % loc(Y*dir) = loc(Q) = median of projected points along direction dir
        % scale(Y*dir) = modified MAD of projected points along direction
        % dir
        
        
        if strcmp(projloc,'median') % Center with median
            
            % Find the median of Q=Y*dir (it is much better to compute the median
            % directly rather than computing function median of stat toolbox)
            Qsor = sort(Q);
            % half = floor(n/2);
            
            me = Qsor(half+1);
            if 2*half == n       % Average if even number of elements
                me =(Qsor(half)+me)/2;
            end
        else % Center with mean
            me=sum(Q)/n;
        end
        
        % Qtilde = centered projected points
        % Qtilde = projected points - estimate of location of projected points
        Qtilde = Q-me;
        
        % Qtildeabs = absolute values of centered projected points
        Qtildeabs = abs(Qtilde);
        
        if strcmp(projscale,'mad')
            % Modified MAD
            ordprojs = sort(Qtildeabs);
            MADmod = (ordprojs(n1)+ordprojs(n2))/(2*beta);
            newoutlvec = Qtildeabs/MADmod;
        elseif strcmp(projscale,'sn')
            SnEst = Sn(Q);
            newoutlvec = Qtildeabs/SnEst;
        elseif strcmp(projscale,'qn')
            QnEst = Qn(Q);
            newoutlvec = Qtildeabs/QnEst;
        elseif strcmp(projscale,'std')
            if strcmp(projloc,'median')
                stdEst = std(Q);
            else
                % if measure of location was the mean then instead of
                % using inefficient function std compute standard
                % deviations using signres (deviations from the mean)
                stdEst = sqrt(sum(Qtilde.^2, 1)/(n-1));
            end
            newoutlvec = Qtildeabs/stdEst;
        else
            warning('FSDA:SDest:WrongScale','Supplied scale measure to standardize scaled projections is not in the list')
            error('FSDA:SDest:WrongScale','You must supply as scale measure one of the following strings ''mad'' ''qn'' ''sn'' or ''std''')
        end
        
        if rstprojsave==1
            % Store in matrix RstProj the n robust standardized projection
            % scores with the sign
            % RstProj(:,i)=newoutlvec;
            if strcmp(projscale,'mad')
                RstProj(:,i)=Qtilde/MADmod;
            elseif strcmp(projscale,'sn')
                RstProj(:,i)=Qtilde/SnEst;
            elseif strcmp(projscale,'qn')
                RstProj(:,i)=Qtilde/QnEst;
            elseif strcmp(projscale,'std')
                RstProj(:,i)=Qtilde/stdEst;
            end
        end
        
        
        % Now check for every unit if the current projection dir produces a
        % greater measure of robust standardized projection score for each unit
        
        % BELOW IS OLD INEFFICIENT CODE
        %     for ii = 1:n
        %         if newoutlvec(ii)>=outlvec(ii)
        %             maxdir(ii,:) = dir';
        %             outlvec(ii) = newoutlvec(ii);
        %         end
        %     end
        %
        % repl = Boolean vector which contains the indexes which have to be
        % replaced
        repl=newoutlvec>=outlvec;
        
        if any(repl)>0
            maxdir(repl,:)=0;
            maxnew=bsxfun(@plus,maxdir(repl,:), dir');
            maxdir(repl,:)=maxnew;
            outlvec(repl)=newoutlvec(repl);
        end
        
    else
        % Juan and Prieto adjustment
        if jpcorr>1
            
            % Yj is a matrix of v+jpcorr-times-v observations
            % Yj = subsample of v+jpcorr units
            
            % Compute the MD using mean and cov based on Yj
            % Inefficiente
            % mah=mahalFS(Yj,mean(Yj),cov(Yj));
            
            meaYj=sum(Yj)/(v+jpcorr);
            % Covariance matrix of Yj
            % Note that cov would recompute the sample means; code below is more
            % efficient
            Ymeac = bsxfun(@minus,Yj,meaYj);
            % It is not necessary to divide by (v+jpcorr-1) because we simply
            % need to find the unit with the largest MD
            covYj = (Ymeac' * Ymeac)/(v+jpcorr-1);
            mah=mahalFS(Yj,meaYj,covYj);
            
            
            [~,indmah]=sort(mah);
            
            % Remove the row of Yj associated with the largest (jpcorr-1) md
            % Yj becomes of size (v+1)-by-v
            Yj(indmah(end-jpcorr+2:end),:)=[];
        end
        
        WDiffk = bsxfun(@minus,Yj(2:v+1,:), Yj(1,:));
        
        % D1 contains (in the columns) v directions
        % Inefficient code to find the inverse
        % D1=pinv(WDiffk);
        % Efficient code to find the inverse
        D1=WDiffk\eyev;
        % where eyev=eye(v);
        
        % q sums every row of matrix D1 (extra direction)
        q=sum(D1,2);     % direction v+1
        
        % Matrix Dk is v-by v+1 matrix and contains (in the columns) v+1 directions
        Dk=[D1 q];
        
        
        % Each column of matrix Dk contains a direction (there are (v+1)
        % directions, v+1 = (v+1 choose v) )
        % The first v columns of Dk (matrix D1) are the solutions of the system of
        % equations WDiffk*D1=I_v     D1=inv(WDiffk)*I_v=inv(WDiffk)
        % The last column of Dk (say q) is the unique solution of the system of equations
        % WDiffk*q=ones(v,1) because q is orthogonal to the rows of WDiffk, and is
        % orthogonal to the the last v rows of matrix Yj, therefore
        % q=inv(WDiffk)*ones(v,1)=D1*ones(v,1)=sum(D1,2) = sum of the rows of
        % matrix D1
        % 1st column of Dk = direction orthogonal to the space spanned by
        % rows 1, 3, 4, ..., v+1 of matrix Yj
        %      that is direction orthogonal to y_1, y_3, y_4, ..., y_{v+1}
        % 2nd column of Dk = direction orthogonal to the space spanned by rows 1, 2, 4, ..., v+1 of matrix Yj
        %      that is direction orthogonal to y_1, y_2, y_4, ..., y_{v+1}
        % ...
        % vth column of Dk = direction orthogonal to the space spanned by rows 1, 2, 3, ..., v of matrix Yj
        %      that is direction orthogonal to y_1, y_2, y_4, ..., y_{v}
        % (v+1)th column of Dk = direction orthogonal to the space spanned by rows 2, 3, 4, ..., v+1 of matrix Yj
        %      that is direction orthogonal to y_2, y_3, y_4, ..., y_{v+1}
        % Note that WDiffk*Dk = [I_v ones(v,1)]
        % Note also that Yj*Dk (if for example v=5) is a matrix
        % (v+1)-by-(v+1) (6 x 6) as follows
        %	k1	k2	k3	k4	k5	f
        %	a	k2	k3	k4	k5	k6
        %	k1	b	k3	k4	k5	k6
        %	k1	k2	c	k4	k5	k6
        %	k1	k2	k3	d	k5	k6
        %	k1	k2	k3	k4	e	k6
        
        % where a, b, c, d, e, f, k1, ..., k6 are scalars in R
        
        % Q is of dimension n x (v+1)
        % Projection of the n points into the v+1 directions (find scores)
        Q=Y*Dk;
        
        if dirsave==1
            % Standardize each direction to unit norm in such a way that
            % diag(Dk'*Dk)=[1 ... 1] of length v+1
            for j=1:v+1
                Dk(:,j)=Dk(:,j)/sqrt(sum(Dk(:,j).^2));
            end
            
            % Store the v+1 directions for subset i inside matrix Dir
            Dir(i,:,:)=Dk;
        end
        
        % Center matrix Q with the medians (default) or with means
        if strcmp(projloc,'median') % Center with medians
            % me=median(Q);
            Qsor = sort(Q,1);
            % Use vectorized method with column indexing.
            me = Qsor(half+1,:);
            if 2*half == n
                me = 0.5*(Qsor(half,:)+me);
            end
        else % Center with means
            me=sum(Q,1)/n;
        end
        
        % Qtilde = centered projected points
        % Qtilde = projected points - estimate of location of projected points
        Qtilde = bsxfun(@minus,Q, me);
        % Qtildeabs = absolute values of centered projected points
        Qtildeabs = abs(Qtilde);
        
        if strcmp(projscale,'mad')
            ress = sort(Qtildeabs);
            % half = floor(n/2);
            mad = ress(half+1,:);
            if 2*half == n       % Average if even number of elements
                mad =(ress(half,:)+mad)/2;
            end
            % Divide by beta (asymptotic consistency factor for MAD)
            mad=mad/beta;
            % Standardize with MADs
            % newoutlmat(:,j)=|(Q(:,j)-median(Q(:,j)))/MAD(Q(:,j))|
            newoutlmat = bsxfun(@rdivide, Qtildeabs, mad);
            % bsxfun instruction is the efficient way of doing
            %    newoutlmat=diag(1../mad)*res';
            %    newoutlmat=newoutlmat';
        elseif strcmp(projscale,'sn')
            SnEst = Sn(Q);
            newoutlmat = bsxfun(@rdivide, Qtildeabs, SnEst);
        elseif strcmp(projscale,'qn')
            QnEst = Qn(Q);
            newoutlmat = bsxfun(@rdivide, Qtildeabs, QnEst);
        elseif strcmp(projscale,'std')
            if strcmp(projloc,'median')
                stdEst = std(Q);
            else
                % if measure of location was the mean then instead of
                % using inefficient function std compute standard
                % deviations using signres (deviations from the mean)
                stdEst = sqrt(sum(Qtilde.^2, 1)/(n-1));
            end
            newoutlmat = bsxfun(@rdivide, Qtildeabs, stdEst);
        else
            disp('Supplied scale measure to standardize scaled projections is not in the list')
            error('FSDA:SDest:WrongScale','You must supply as scale measure one of the following strings ''mad'' ''qn'' ''sn'' or ''std''')
        end
        
        
        if rstprojsave==1
            if strcmp(projscale,'mad')
                %RstProj(:,i,:)=newoutlmat;
                RstProj(:,i,:)=bsxfun(@rdivide, Qtilde, mad);
            elseif strcmp(projscale,'sn')
                RstProj(:,i,:)=bsxfun(@rdivide, Qtilde, SnEst);
            elseif strcmp(projscale,'qn')
                RstProj(:,i,:)=bsxfun(@rdivide, Qtilde, QnEst);
            elseif strcmp(projscale,'std')
                RstProj(:,i,:)=bsxfun(@rdivide, Qtilde, stdEst);
            end
        end
        
        % Find the largest element of each row and store the associated
        % direction
        [newoutlvec,maxinds]=max(newoutlmat,[],2);
        
        repl=newoutlvec>=outlvec;
        if any(repl)>0
            % Store direction of maximum outlyingness
            maxdir(repl,:)=Dk(:,maxinds(repl))';
            
            % Store outlyingness measure
            outlvec(repl)=newoutlvec(repl);
        end
        
        
    end
    
    % Write total estimation time to compute final estimate
    if i <= tsampling
        
        % sampling time until step tsampling
        time(i)=toc;
    elseif i==tsampling+1
        % stop sampling and print the estimated time
        if msg==1
            fprintf('Total estimated time to complete Stahel-Donoho estimator: %5.2f seconds \n', nselected*median(time));
        end
    end
    
end


% Standardize direction of maximum outlyingness
% for Juan and Prieto method
if jpcorr~=0;
    nor=sqrt(sum(maxdir.^2,2));
    maxdir = bsxfun(@rdivide, maxdir, nor);
end

% If necessary do the marginal projections
% For example if margin =3 consider
% each marginal univariate projection
% 1 0 ..... 0
% 0 1 ..... 0
% .........
% 0 0 ..... 1
% each bivariate marginal projection
% 1 1 0 .....0 0
% 0 1 1 .....0 0
% 0 0 1 .....0 0
% ......
% 0 0 0 .....1 1
% each trivariate marginal projection
% 1 1 1 .....0 0 0
% 0 1 1 .....0 0 0
% 0 0 1 .....0 0 0
% ......
% 0 0 0 .....1 1 1

if margin >0
    
    seqv=1:v;
    for j=1:margin
        
        marg=combsFS(seqv,j);
        normj=sqrt(j);
        
        for im=1:size(marg,1);
            % Ysel = extract just the columns of Y specified by marg(im,:)
            Yselori=Y(:,marg(im,:));
            
            % Project using original sign for all the columns of Yselori
            Q=sum(Yselori,2)/normj;
            
            % OLD VERSION TO DELETE
            %             Qtildeabs = abs(Q-median(Q));
            %             ordprojs = sort(Qtildeabs);
            %             % Modified MAD
            %             MADmod = (ordprojs(n1)+ordprojs(n2))/(2*beta);
            %             newoutlvec = Qtildeabs/MADmod;
            
            if strcmp(projloc,'median') % Center with median
                
                % Find the median of Q=Y*dir (it is much better to compute the median
                % directly rather than computing function median of stat toolbox)
                Qsor = sort(Q);
                % half = floor(n/2);
                
                me = Qsor(half+1);
                if 2*half == n       % Average if even number of elements
                    me =(Qsor(half)+me)/2;
                end
            else % Center with mean
                me=sum(Q)/n;
            end
            
            % Qtilde = centered projected points
            % Qtilde = projected points - estimate of location of projected points
            Qtilde = Q-me;
            
            % Qtildeabs = absolute values of centered projected points
            Qtildeabs = abs(Qtilde);
            
            if strcmp(projscale,'mad')
                % Modified MAD
                ordprojs = sort(Qtildeabs);
                MADmod = (ordprojs(n1)+ordprojs(n2))/(2*beta);
                newoutlvec = Qtildeabs/MADmod;
            elseif strcmp(projscale,'sn')
                SnEst = Sn(Q);
                newoutlvec = Qtildeabs/SnEst;
            elseif strcmp(projscale,'qn')
                QnEst = Qn(Q);
                newoutlvec = Qtildeabs/QnEst;
            elseif strcmp(projscale,'std')
                if strcmp(projloc,'median')
                    stdEst = std(Q);
                else
                    % if measure of location was the mean then instead of
                    % using inefficient function std compute standard
                    % deviations using signres (deviations from the mean)
                    stdEst = sqrt(sum(Qtilde.^2, 1)/(n-1));
                end
                newoutlvec = Qtildeabs/stdEst;
            else
                warning('FSDA:SDest:WrongScale','Supplied scale measure to standardize scaled projections is not in the list')
                error('FSDA:SDest:WrongScale','You must supply as scale measure one of the following strings ''mad'' ''qn'' ''sn'' or ''std''')
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            for ii = 1:n
                if newoutlvec(ii)>=outlvec(ii)
                    maxdir(ii,:) = zeros(1,v);
                    maxdir(ii,marg(im,:)) = 1;
                    outlvec(ii) = newoutlvec(ii);
                end
            end
            
            
            
            % Add the projections with the sign
            % For every row of marg(im,:) it is necessary an additional
            % loop to cope not only with +1...+1 but also with the -1
            % For example
            % if marg(im,:)= 1 2 3 The first three column can be
            % projected not only with +1 +1+1, but also with
            % -1 1 1
            % 1 -1 1
            % 1 1 -1
            % if marg(im,:)= 1 2 3 4 The first four columns can be
            % projected not only with +1 +1+1 +1, but also with
            % -1  1  1  1
            %  1 -1  1  1
            %  1  1 -1  1
            %
            % -1 -1  1  1
            % -1  1 -1  1
            % -1  1  1 -1
            %  1 -1 -1  1
            %  1 -1  1 -1
            %  1  1 -1 -1
            
            % The columns of Yselori can enter the projection with sign +1
            % or -1
            
            if j>1
                sizYselori=size(Yselori,2);
                
                for sign=1:fix(sizYselori/2);
                    
                    
                    onesign=ones(sizYselori,1);
                    
                    % if j==2 and sign=1 the only one possible case is
                    % +1 -1
                    if j==2 && sign==1
                        signc=1;
                    else
                        signc=combsFS(1:sizYselori,sign);
                    end
                    
                    for iii=1:size(signc,1)
                        onesignc=onesign;
                        onesignc(signc(iii,:))=-1;
                        
                        
                        % Project the selected column(s) taking care of the sign
                        Ysel=Yselori*onesignc;
                        Q=Ysel/normj;
                        
                        % OLD CODE TO DELETE
                        %                         Qtildeabs = abs(Q-median(Q));
                        %                         ordprojs = sort(Qtildeabs);
                        %                         % Modified MAD
                        %                         MADmod = (ordprojs(n1)+ordprojs(n2))/(2*beta);
                        %                         newoutlvec = Qtildeabs/MADmod;
                        %%%%%%
                        
                        if strcmp(projloc,'median') % Center with median
                            
                            % Find the median of Q=Y*dir (it is much better to compute the median
                            % directly rather than computing function median of stat toolbox)
                            Qsor = sort(Q);
                            % half = floor(n/2);
                            
                            me = Qsor(half+1);
                            if 2*half == n       % Average if even number of elements
                                me =(Qsor(half)+me)/2;
                            end
                        else % Center with mean
                            me=sum(Q)/n;
                        end
                        
                        % Qtilde = centered projected points
                        % Qtilde = projected points - estimate of location of projected points
                        Qtilde = Q-me;
                        
                        % Qtildeabs = absolute values of centered projected points
                        Qtildeabs = abs(Qtilde);
                        
                        if strcmp(projscale,'mad')
                            % Modified MAD
                            ordprojs = sort(Qtildeabs);
                            MADmod = (ordprojs(n1)+ordprojs(n2))/(2*beta);
                            newoutlvec = Qtildeabs/MADmod;
                        elseif strcmp(projscale,'sn')
                            SnEst = Sn(Q);
                            newoutlvec = Qtildeabs/SnEst;
                        elseif strcmp(projscale,'qn')
                            QnEst = Qn(Q);
                            newoutlvec = Qtildeabs/QnEst;
                        elseif strcmp(projscale,'std')
                            if strcmp(projloc,'median')
                                stdEst = std(Q);
                            else
                                % if measure of location was the mean then instead of
                                % using inefficient function std compute standard
                                % deviations using signres (deviations from the mean)
                                stdEst = sqrt(sum(Qtilde.^2, 1)/(n-1));
                            end
                            newoutlvec = Qtildeabs/stdEst;
                        else
                            disp('Supplied scale measure to standardize scaled projections is not in the list')
                            error('FSDA:SDest:WrongScale','You must supply as scale measure one of the following strings ''mad'' ''qn'' ''sn'' or ''std''')
                        end
                        
                        %%%%%%%%%%%
                        for ii = 1:n
                            if newoutlvec(ii)>=outlvec(ii)
                                maxdir(ii,:) = zeros(1,v);
                                maxdir(ii,marg(im,:)) = onesignc';
                                outlvec(ii) = newoutlvec(ii);
                            end
                        end
                        
                    end   % loop associeated with iii
                    
                end % loop associated with sign
            end
        end % loop associated with im=1:size(marg,1)
    end % loop associated with j=1:margin
end % loop associated with if margin>0

% constant q of huber weight function
q = options.q;
% string c of huber weight function
c = options.c;
% constant nbp of Tukey biweight function
nbp = options.nbp;
% constant K of Zuo, Cui and He's weight function
K = options.K;

weight=options.weight;
% Huber weight function
if strcmp(weight,'huber')
    if strcmp(c, 'hdim')
        c1 = min(sqrt(chi2inv(0.50,v)),4);
    elseif strcmp(c, 'sdim')
        c1 = sqrt(chi2inv(0.95,v));
    end
    weights = huber(outlvec, c1, q);
    factor=1;
    %  Tukey biweight function
elseif strcmp(weight, 'tukey')
    c2 = TBbdp(nbp,v);
    weights = TBwei(outlvec, c2);
    factor=1;
    % Zuo, Cui and He's family of weights
elseif strcmp(weight, 'zch')
    weights = zch(outlvec, K);
    factor=1;
    % mcd weights
else
    [~,iwei]=sort(outlvec);
    weights=zeros(n,1);
    weights(iwei(1:half))=1;
    factor=consistencyfactor(half,n,v);
    
end



% Find robust estimate of location (using weights previously found)
% loc = \sum_{i=1}^n y_i'*w_i / \sum_{i=1}^n w(d_i)
sumw=sum(weights);
loc = sum(bsxfun(@times,Y,weights),1)/sumw;
% Inefficient code is loc = sum(repmat(weights',1,v).*Y) / sum(weights);

% Res = n x v matrix which contains deviations from the robust estimate
% of location
Res = bsxfun(@minus,Y, loc);
cov1= (Res')*bsxfun(@times,Res,weights)/sumw;
% Inefficient code is as follows
%    Y1 = repmat(sqrt(weights'),1,v) .* (Y - repmat(loc,n,1));
%   cov = (Y1'*Y1) / sum(weights);

%Output:

%=============================================================================

out.class   = 'SD';
out.loc     =  loc; %robust estimate of location


out.cov     =  factor*cov1; %robust estimate of covariance matrix
out.maxdir  = maxdir; % direction which produces the largest projection
out.weights = weights;
out.md = mahalFS(Y,out.loc,out.cov);

% Store in output structure the outliers found with confidence level conflev
conflev = options.conflev;
seq = 1:n;
out.outliers = seq(out.md > chi2inv(conflev,v) );
out.conflev = conflev;

if dirsave==1
    out.Dir=Dir;
end

if rstprojsave==1
    out.RstProj=RstProj;
end

plo=options.plots;

% Plot Mahalanobis distances with outliers highlighted
if isstruct(plo) || (~isstruct(plo) && plo~=0)
    
    laby='SD Mahalanobis distances';
    malindexplot(out.md,v,'conflev',conflev,'laby',laby,'numlab',out.outliers);
    
    figure('Tag','pl_spm_outliers');
    group=ones(n,1);
    if ~isempty(out.outliers)
        group(out.outliers)=2;
    end
    spmplot(Y,group,plo);
    set(gcf,'Name',' SD estimator: scatter plot matrix with outliers highlighted');
    
    figure('Tag','weights_outliers');
    col = 'b';
    symb = 'x';
    if ~isempty(out.outliers)
        group(out.outliers)=2;
        col = ['b', 'r'];
        symb = ['x', 'o'];
    end
    %     gscatter(weights, out.md, group,'br','xo');
    gscatter(weights, out.md, group, col, symb);
    legend('Normal units');
    if ~isempty(out.outliers)
        legend('Normal units', 'Outliers');
    end
    xlabel('observation weights');
    ylabel('SD Mahalanobis distances');
    title('SDest: plot of robust MDs against weights');
    set(gcf,'Name',' SD estimator: weights-by-Mahalanobis distances plot');
    
end

if options.ysave
    out.Y = Y;
end

end



function w = huber(x, c, q)

u = x./c;
w = (u<=1) + (u + 1.e-9).^(-q).*(u>1);

end

function rawconsfac=consistencyfactor(h,n,v)
a=chi2inv(h/n,v);
rawconsfac=(h/n)/(chi2cdf(a,v+2));
end


function wz = zch(x, K)
n = length(x);
half = floor(n/2);
xadd = x+1;
% PD = Projection Depth
PD = 1./xadd;
PDsor = sort(PD);
%  Median of PD
C = PDsor(half+1);
if 2*half == n
    C = (PDsor(half)+C)/2;
end

u = PD./C;
expK = exp(-K);
b = 1-expK;
% expdiff = (expK).^((1-u).*(1-u)) - expK;
% expdiff = expK*( expK.^(u.*(u-2)) - 1);
expdiff = exp(-K*(1-u).^2) - expK;
wzarg = expdiff./b;
wz = (u>=1) + wzarg.*(u<1);
end

%FScategory:MULT-Multivariate