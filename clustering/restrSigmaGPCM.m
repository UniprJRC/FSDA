function [Sigma, lmd, OMG, GAM]  = restrSigmaGPCM(SigmaB, niini, pa, nocheck, lmd, OMG)
%restrSigmaGPCM computes constrained covariance matrices for the 14 GPCM specifications
%
%
%<a href="matlab: docsearchFS('restrSigmaGPCM')">Link to the help function</a>
%
%
%  This routine applies the constraints to the covariance matrices using the
%  specifications contained in input structure pa.
%
%
% Required input arguments:
%
%   SigmaB : initial unconstrained covariance matrices. p-by-p-by-k array.
%            p-by-p-by-k array containing the k covariance matrices for the
%            k groups.
%               Data Types - single|double
%
%   niini  : sizes of the groups. Vector. Vector of length k containing
%           the size of the groups. ninini can be either a row or a columns
%           vector.
%               Data Types - single|double
%
%      pa  : Constraints to apply and model specification. Structure.
%            Structure containing the following fields:
%             pa.pars= type of Gaussian Parsimonious Clustering Model. Character.
%               A 3 letter word in the set:
%               'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',
%               'EVI','VEI','EEI','VII','EII'.
%               The field pa.pars is compulsory. All the other fields are
%               non necessary. If they are not present they are set to
%               their default values.
%             pa.cdet = scalar in the interval [1 Inf) which specifies the
%               the restriction which has to be applied to the determinants.
%               If pa.cdet=1 all determinants are forced to be equal.
%               See section More About for additional details.
%             pa.shw = scalar in the interval [1 Inf) which specifies the
%               the restriction which has to be applied to the elements of
%               the shape matrices inside each group. If pa.shw=1 all diagonal
%               elements of the shape matrix of cluster j (with j=1, ...,
%               k) will be equal.
%             pa.shb = scalar in the interval [1 Inf) which specifies the
%               the restriction which has to be applied to the elements of
%               the shape matrices across each group.
%             pa.maxiterS = positive integer which specifies the maximum
%               number of iterations to obtain the restricted shape matrix.
%               This parameter is used by routine restrshapepars. The
%               default value of pa.maxiterS is 5.
%             pa.maxiterR = positive integer which specifies the maximum
%               number of iterations to obtain the common rotation matrix
%               in presence of varying shape.
%               This parameter is used by routine cpcV. The
%               default value of pa.maxiterR is 20.
%          pa.maxiterDSR = positive integer which specifies the maximum
%               number of iterations to obtain the requested restricted
%               determinants, shape matrices and rotation. For all
%               parametrizations  pa.maxiterDSR is set to 1 apart from for
%               the specifications 'VVE', 'EVE' and 'VEE'. The default
%               value of pa.maxiterDSR is 20.
%           pa.tolS=tolerance to use to exit the iterative procedure for
%               estimating the shape. Scalar. The
%               iterative procedures stops when the relative difference of
%               a certain output matrix is smaller than itertol in two consecutive
%               iterations. The default value of pa.tol is 1e-12.
%      pa.zerotol = tolerance value to declare all input values equal to 0
%               in the eigenvalues restriction routine (file restreigen.m)
%               or in the final reconstruction of covariance matrices.
%               The default value of zerotol is 1e-10.
%          pa.msg = boolean which if set equal to true enables to monitor
%               the relative change of the estimates of lambda Gamma and
%               Omega in each iteration. The default value of pa.msg is
%               false, that is nothing is displayed in each iteration.
%           pa.k  = the number of groups.
%           pa.v  = the number of variables.
%   pa.userepmat  = scalar, which specifies whether to use implicit
%                   expansion or bsxfun.  pa.userepmat =2 implies implicit
%                   expansion, pa.userepmat=1 implies use of bsxfun. The
%                   default is to use implicit expansion (faster)
%                   if verLessThanFS(9.1) is false and bsxfun if MATLAB is
%                   older than 2016b.
%               Data Types - struct
%
%    nocheck : check in input option pa. Boolean. Specify whether it is
%               necessary to check the input fields in
%               previous input option pa. If nocheck is
%               false (default is true) no check is performed on input
%               structure pa.
%               Data Types - Boolean
%
%  Optional input arguments:
%
%     lmd  : determinants. Vector.
%             Initial estimates of (constrained) determinants
%                 Example - [ 2 4 6]
%                 Data Types - double
%      OMG : rotation matrices. p-by-p-by-k array.
%             p-by-p-by-k array containing the preliminary estimates of the
%             rotation matrices for the k groups. If common rotation is
%             imposed (third letter is equal to E),
%             OMG(:,:,1)=...=OMG(:,:,k).
%                 Example - .5*hadamard(4)
%                 Data Types - double
%
% Output:
%
%
%    Sigma  : constrained covariance matrices. p-by-p-by-k array.
%             p-by-p-by-k array containing the k covariance matrices for
%             the k groups. See section 'More About' for the notation for
%             the eigen-decomposition of the component covariance matrices.
%     lmd  : restricted determinants. Vector.
%             Row vector of length $k$ containing restricted determinants.
%             More precisely, the $j$-th element of lmd contains
%             $\lambda_j^{1/p}$. The elements of lmd satisfy the
%             constraint pa.cdet in the sense that $\max(lmd) / \min(lmd)
%             \leq pa.cdet^{(1/p)}$. In other words, the ratio between the
%             largest and the smallest determinant is not greater than
%             pa.cdet. All the elements of vector lmd are equal if
%             modeltype is E** or if pa.cdet=1.
%      OMG : constrained rotation matrices. p-by-p-by-k array.
%             p-by-p-by-k array containing the k rotation matrices for
%             the k groups. If common rotation is imposed (third letter is
%             equal to E), OMG(:,:,1)=...=OMG(:,:,k).
%     GAM : constrained shape matrix. 2D array.
%           Matrix of size p-by-k containing in
%           column j the elements on the main diagonal of shape matrix
%           $\Gamma_j$. The elements of GAM satisfy the following
%           constraints:
%           The product of the elements of each column is equal to 1.
%           The ratio among the largest elements of each column is
%           not greater than pa.shb.
%           The ratio among the second largest elements of each column is
%           not greater than pa.shb.
%           ....
%           The ratio among the smallest elements of each column is
%           not greater than pa.shb.
%           The ratio of the elements of each column is not greater than
%           pa.shw.
%
%
% More About:
% The notation for the eigen-decomposition of the
% component covariance matrices is as follows
%
% \[
% \Sigma_j= \lambda_j^{1/p} \Omega_j \Gamma_j \Omega_j'  \qquad j=1, 2, \ldots, k
% \]
% The dimension of matrices $\Omega_j$ (rotation) and $\Gamma_j$ (shape) is $p\times p$.
%
% $c_{det}=$ scalar, constraint associated with the determinants.
%
% $c_{shw}=$ scalar, constraint inside each group of the shape matrix.
%
% $c_{shb}=$ scalar, constraint among groups of the shape matrix.
%
% Note that if you impose equal volumes $c_{det}=1$. Similarly, if you
% impose a spherical shape $c_{shw}=1$.
%
% We also denote with
%
%   [1] $\Sigma_B$ the 3D array of size $p\times p \times k$ containing the
%   empirical covariance matrices of the $k$ groups, before applying the
%   constraints coming from the 14 parametrizations. In the code $\Sigma_B$
%   is called $SigmaB$. The $j$-th slice of this 3D array of size $p\times
%   p$ is denoted with symbol $\hat \Sigma_j$.
%   [2] $\Omega$ the 3D array of size $p\times p \times k$ containing the
%   rotation matrices of the $k$ groups. In the code $\Omega$  is called
%   $OMG$. The $j$-th slice of this 3D array of size $p\times p$ is called
%   $\hat \Omega_j$.
%   [3] $\Gamma$ the $p\times k$ matrix containing in column $j$, with
%   $j=1, 2, \ldots, k$, the diagonal elements of matrix $\Gamma_j$ (shape
%   matrix of group j). In the code matrix $\Gamma$ is called GAM.
%   After the application of this routine, the product of the elements of
%   each column of matrix GAM is equal to 1.
%   The ratio among the largest (second largest, ...smallest) elements of
%   each column is not greater than $c_{shb}$ (pa.shb).
%   The ratio of the elements of each column is not greater than
%   $c_{shw}$ (pa.shb). All the columns of matrix GAM are equal if the second
%   letter of modeltype is E. All the columns of matrix GAM are
%   equal to 1 if the second letter of modeltype is I.
%   [4] niini the vector of length $k$ containing the number of units
%   (weights) associated to each group.
%   [5]  $\lambda$ = the vector of length $p$ containing in the $j$-th
%   position $\lambda_j^{1/p}=|\Sigma_j|^{1/p}$. In the code  vector
%   $\lambda$ is called $lmd$.
%   The elements of lmd satisfy the constraint pa.cdet in the sense that
%   $\max(lmd) / \min(lmd) \leq pa.cdet^{(1/p)}$. In other words, the
%   ratio between the largest and the smallest determinant is not
%   greater than pa.cdet. All the elements of vector lmd are equal
%   if modeltype is E** or if $c_{det}=1$ (pa.cdet=1).
%
%
%
% See also restrshapeGPCM, restrdeterGPCM, restreigen, tclust
%
%
% References:
%
%   Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2019),
%   Robust parsimonious clustering models. Submitted.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('restrSigmaGPCM')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%{
    % Use of restrSigmaGPCM with just one output argument.
    % Generate 3 cov matrices of size 2-by-2
    rng('default')
    rng(1500)
    Sig = [1 .5; .5 2];
    df = 10;
    k=2;
    v=2;
    Sigma=zeros(v,v,k);
    for j=1:k
        Sigma(:,:,j) = wishrnd(Sig,df)/df;
    end
    niini=2;
    % Apply model EEE.
    pa=struct;
    pa.pars='EEE';
    SigmaNEW  = restrSigmaGPCM(Sigma, niini, pa);
%}

%{
    %% Use of restrSigmaGPCM with all default options.
    % Generate 3 cov matrices of size 2-by-2
    rng('default')
    rng(1500)
    Sig = [1 .5; .5 2];
    df = 10;
    k=3;
    v=2;
    Sigma=zeros(v,v,k);
    for j=1:k
        Sigma(:,:,j) = wishrnd(Sig,df)/df;
    end
    niini=[20 150 200];
    % Apply model VVE.
    pa=struct;
    pa.pars='VVE';
    [SigmaNEW, lmd, OMG, GAM]  = restrSigmaGPCM(Sigma, niini, pa);
%}

%{
    %% Use function genSigmaGPCM to generate the covariance matrices.
    v=3;
    k=5;
    pa=struct;
    pa.pars='EVI';
    S=genSigmaGPCM(v, k, pa);
    niini=100*ones(1,k);
    pa.pars='VVE';
    [SigmaNEW, lmd, OMG, GAM]  = restrSigmaGPCM(S, niini, pa);
%}

%{
    %% Generate ad hoc cov matrices.
    k=7; v=20; n=100;
    rng('default')
    seed=1141;
    add=ones(v,v)+diag(1:v);
    niini= round(100*mtR(k,0,seed));
    Sigma=zeros(v,v,k);
    for j=1:k
        Sigma(:,:,j)=cov(reshape(mtR(n*v,1,-1),n,v)).*add;
    end
    sph=struct;
    sph.pars='EVV';
    niini=100*ones(k,1);
    [SigmaNEW, lmd, OMG, GAM]  = restrSigmaGPCM(Sigma, niini', sph);
%}

%% Beginning of code

% Set default values for tolerances and maximum number of iterations.
maxiterDSRdef=20;
tolDSRdef=1e-5;
maxiterSdef=20;
tolSdef=1e-5;
maxiterRdef = 20;
tolRdef = 1e-5;
shwdef=100;
shbdef=100;
cdetdef=100;
zerotoldef=1e-10;

if nargin<4
    nocheck=false;
end


% pa = structure containing modeltype, number of iterations .....
fpa=fieldnames(pa);

% if nocheck is false check that the fieldnames of input structure pa are
% those specified in input structure options.
if nocheck == false
    options=struct('maxiterDSR',maxiterDSRdef,'tolDSR',tolDSRdef,'maxiterS',maxiterSdef,'tolS',tolSdef, ...
        'maxiterR',maxiterRdef,'tolR',tolRdef,'shw',shwdef,'shb',shbdef,...
        'cdet',cdetdef,'zerotol',zerotoldef,'pars','','userepmat','');
    chkoptions(options,fpa)
end

% SigmaB = p-times-p-times-k = empirical covariance matrix
Sigma=SigmaB;
k=length(niini);
v=size(Sigma,1);
pa.v=v;
pa.k=k;

if nargin<5
    % OMG = initialize 3D array containing rotation matrices
    OMG=zeros(size(Sigma));
end

% Tolerance associated to the maximum of the elements which have to be
% constrained. If the maximum  is smaller than zerotol restreigen
% procedures returns in output what has been given in input. For example,
% if the elements which have to be constrained are the eigenvalues of the
% covariance matrices and the max of the eigenvalues is smaller than
% zerotol it means that all n points are concentrated in k points and there
% is a perfect fit therefore no further changes on the eigenvalues is
% required.
d=max(strcmp('zerotol',fpa));
if d==0
    pa.zerotol=zerotoldef;
end
zerotol=pa.zerotol;

% pa.pars = character vector with three letters specifying the type of the
% 14 constraints (i.e. EEE, CVVV, EVE, ...)
pars=pa.pars;

% Select cases in which maxiterDSR>1 and specify relative tolerance
% EVE VEE VVE VVV VEV VVI and VEI require iterations
% All the other specification do not
d=max(strcmp('maxiterDSR',fpa));
if d==0
    pa.maxiterDSR=maxiterDSRdef;
end

% Model which require more than one iteration in the main loop
modDSR={'EVE' 'VEE' 'VVE' 'VVV' 'VEV' 'VVI' 'VEI'};
if max(strcmp(pars,modDSR))>0
    maxiterDSR=pa.maxiterDSR;
else
    maxiterDSR = 1;
end

d=max(strcmp('tolDSR',fpa));
if d==0
    tolDSR=tolDSRdef;
else
    tolDSR=pa.tolDSR;
end

d=max(strcmp('maxiterS',fpa));
if d==0
    pa.maxiterS=maxiterSdef;
end

% Cases with different shape (they require iteration for the shape)
diffSHP={'VVE','EVE','VVV','EVV','VVI','EVI'};
if max(strcmp(pars,diffSHP))==0
    pa.maxiterS=1;
end

d=max(strcmp('tolS',fpa));
if d==0
    pa.tolS=tolSdef;
end

d=max(strcmp('tolDSR',fpa));
if d==0
    pa.maxiterR=maxiterRdef;
end

eqROTdiffSHP={'EVE','VVE'};
if max(strcmp(pars,eqROTdiffSHP))>0
    % maxiterR=pa.maxiterR;
    cpc=true;
else
    cpc=false;
    pa.maxiterR=1;
end

d=max(strcmp('tolR',fpa));
if d==0
    pa.tolR=tolRdef;
end

d=max(strcmp('shw',fpa));
if d==0
    pa.shw=shwdef;
end

d=max(strcmp('shb',fpa));
if d==0
    pa.shb=shbdef;
end

d=max(strcmp('cdet',fpa));
if d==0
    pa.cdet=cdetdef;
end

d=max(strcmp('msg',fpa));
if d==1
    msg=pa.msg;
else
    msg=false;
end


d=max(strcmp('userepmat',fpa));
if d==0
    verLess2016b=verLessThanFS(9.1);
    if verLess2016b == true
        pa.userepmat=1;
    else
        pa.userepmat=2;
    end
end

%% Parameters not set by the user

if strcmp(pars(3),'E') || strcmp(pars(3),'I')
    pa.sortsh=1;
    GAMfc=zeros(v,k);
else
    pa.sortsh=0;
end

if strcmp(pars(1),'E')
    pa.cdet=1;
end

% If OMG is identity, shape restriction parameter within groups is set to 1
if strcmp(pars(2),'I')
    pa.shw = 1;
end

% if Equal shape is imposed shape restriction parameter between groups is
% set to 1
if strcmp(pars(2),'E')
    pa.shb = 1;
end

%% Initialization part
if strcmp(pars(3),'E')
    % In the common principal components case it is necessary to find
    % initial values for OMG (rotation), and lmd (unconstrained
    % determinants)
    if nargin<5
        % fifth argin is lm and 6th argin is OMG
        [lmd, OMG]  = initR(SigmaB, niini, pa);
    else
        %         % Initialize lmd
        %         lmd = NaN(1,k);
        %
        %         if strcmp(pars(1),'V')
        %             for j=1:k
        %                 % lmd(j)=exp(log(det(SigmaB(:,:,j)))/v);
        %                 lmd(j) = (det(SigmaB(:,:,j))) ^ (1 / v);
        %             end
        %         else
        %             lmd = ones(1,k);
        %         end
    end
    
    % In presence of variable shape
    % compute Wk and wk once and for all.
    % Wk(:,:,j) contains (n_j/n) \Sigma_j
    % wk(j) contains largest eigenvalue of Wk(:,:,j)
    % These two matrices will be used inside routine cpcV
    if (strcmp(pars,'VVE') || strcmp(pars,'EVE'))
        Wk=zeros(v,v,k);
        wk=zeros(k,1);
        sumnini=sum(niini);
        for j=1:k
            Wk(:,:,j)  = (niini(j) /sumnini)  * SigmaB(:,:,j);
            wk(j) = max(eig(Wk(:,:,j)));
        end
    end
elseif  strcmp(pars(3),'V')
    if nargin<5
        % Initialize lmd
        lmd = ones(1,k);
        
        % Find initial (and final value for OMG)
        for j=1:k
            [V,eigunsorted]= eig(SigmaB(:,:,j));
            diageigunsorted=diag(abs(eigunsorted));
            % Sort eigenvalues from largest to smallest and reorder the columns
            % of the matrix of eigenvectors accordingly
            [~,ordeig]=sort(diageigunsorted,'descend');
            V=V(:,ordeig);
            OMG(:,:,j)=V;
            
            if strcmp(pars(1),'V')
                % lmd(j) = (det(SigmaB(:,:,j))) ^ (1 / v);
                lmd(j) = (prod(diageigunsorted)) ^ (1 / v);
            end
        end
    end
else % The remaining case is when **I
    if nargin<5
        % Find initial (and final value for OMG).
        % in this case OMG is a 3D arry contaning identity matrices
        eyep=eye(pa.v);
        for j=1:k
            OMG(:,:,j)=eyep;
        end
        
        % Initialize lmd
        lmd = ones(1,k);
        if strcmp(pars(1),'V')
            for j=1:k
                lmd(j) = (det(SigmaB(:,:,j))) ^ (1 / v);
            end
        end
    end
end

GAM=ones(v,k);
% Immediately apply the restriction on vector lmd
% if ~isequal(lmd,ones(1,k))
%
%     [lmd]=restrdeterGPCM(GAM, OMG, SigmaB, niini, pa);
% end

OMGold=OMG(:,:,1);
GAMold=9999;
lmdold=GAMold;

%% Beginning of iterative process

% Apply the iterative procedure to find Det, Shape and Rot matrices
iter=0;
diffglob = Inf;

if msg==true
    disp(['   Iter' '    diff_lmd' '    diff_GAM' '   diff_OMG'])
end

while ( (diffglob > tolDSR) && (iter < maxiterDSR) )
    iter=iter+1;
    
    % In the **E case (except for the case EEE) it is necessary to update in each step of the
    % iterative procedure OMG
    if iter>1 &&  strcmp(pars(3),'E')
        
        if cpc == true
            % Variable shape: update OMG (rotation)
            % parameter pa.maxiterR is used here
            [OMG]  = cpcV(lmd, GAM, OMG(:,:,1), Wk, wk, pa);
        else
            % Equal shape: update OMG (rotation)
            [OMG]  = cpcE(lmd, SigmaB, niini, pa);
            % In all the other cases OMG is not updated
        end
        
        % diffOMG is the relative sum of squares of the differences between
        % the element of matrix Omega2D in two consecutive iterations
        OMGnew=OMG(:,:,1);
        diffOMG=abs((v-sum(sum( (OMGnew'*OMGold).^2 )) )/v);
        OMGold=OMGnew;
    else
        diffOMG=0;
    end
    
    
    % Update GAM
    [GAM] =restrshapeGPCM(lmd, OMG, SigmaB, niini, pa);
    % GAMf=GAM;
    if pa.sortsh==1
        for j=1:k
            GAMfc(:,j)=sort(GAM(:,j),'descend');
        end
    else
        GAMfc=GAM;
    end
    
    % GAMnew = new values of matrix GAM in vectorized form
    GAMnew=GAMfc(:);
    % diff = (new values of GAM - old values of GAM)
    diff=GAMnew-GAMold;
    % relative sum of squares of the differences
    diffGAM=diff'*diff/(GAMold'*GAMold);
    GAMold=GAMnew;
    
    % Update determinants in case of varying determinants (apart from VII)
    % Update lmd
    [lmd]=restrdeterGPCM(GAM, OMG, SigmaB, niini, pa);
    
    % lmdnew = new values of vector lmd
    lmdnew=lmd(:);
    % diff = (new values of lmd) - (old values of lmd)
    diff=lmdnew-lmdold;
    % relative sum of squares of the differences
    difflmd=diff'*diff/(lmdold'*lmdold);
    lmdold=lmdnew;
    
    diffglob=max([difflmd, diffGAM,diffOMG,]);
    if msg==true
        disp([iter difflmd, diffGAM, diffOMG])
    end
end

% Check if all is well
codenonZero = ~any(isnan(GAM(:))) && max(max(GAM)) > zerotol;

if  codenonZero
    % Reconstruct the cov matrices using final values of lmd, OMG and GAM
    for j=1:k
        Sigma(:,:,j) = lmd(j) * OMG(:,:,j)* diag(GAM(:,j))* (OMG(:,:,j)');
    end
else
    % In this case final Sigma is the identity matrix replicated k times
    eyep=eye(v);
    for j=1:k
        Sigma(:,:,j) = eyep;
    end
end

end
%FScategory:CLUS-RobClaMULT
