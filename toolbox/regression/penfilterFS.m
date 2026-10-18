function [out,Aready] = penfilterFS(y,K,N,varargin)
%penfilterFS fits a generic quadratic penalized trend on a selected subset
%
%<a href="matlab: docsearchFS('penfilterFS')">Link to the help function</a>
%
% penfilterFS is a low-level fitting routine for quadratic smoothers in the
% observation space. Given a symmetric positive semidefinite penalty matrix
% K, it solves
%
%   min_{m,beta} ||S(y-m-X*beta)||^2 + lambda*m'*K*m
%                            + beta'*betaPenalty*beta,
%
% where S selects the observations in bsb. If X is empty, the criterion
% reduces to
%
%   min_m ||S(y-m)||^2 + lambda*m'*K*m.
%
% Let N be a matrix whose columns span null(K), and let B0 span the null
% space of betaPenalty. Conditional identification on bsb requires
%
%   rank([N(bsb,:) X(bsb,:)*B0]) = size(N,2)+size(B0,2).
%
% The matrix N is therefore a first-class input of the routine. This avoids
% recomputing a numerical null space in repeated robust-fitting calls and
% makes the elemental-subset dimension explicit. If K is positive definite,
% use N=zeros(T,0) or N=[].
%
% penfilterFS intentionally performs fitting only. Variance estimation and
% prediction intervals depend on the stochastic interpretation attached to
% a particular penalty and are left to wrappers such as hpfilterFS.
%
%  Required input arguments:
%
%    y :        Time series to analyze. Vector or timetable.
%               Vector with T elements containing the response. Missing
%               values are allowed outside bsb. If y is a timetable, the
%               first variable is used and row times are used in plots.
%                 Data Types - double | timetable
%
%    K :        Quadratic penalty matrix. T-by-T symmetric positive
%               semidefinite matrix. The penalty contribution is
%               lambda*m'*K*m. K can be dense or sparse.
%               For the Hodrick--Prescott penalty K=D2'*D2, where D2 is
%               the second-difference matrix.
%                 Data Types - double | sparse double
%
%    N :        Basis of null(K). T-by-r matrix, where r is the nullity of
%               K. The columns of N must span the complete null space of K.
%               The basis need not be orthonormal; penfilterFS internally
%               orthonormalizes it. If K is positive definite, use [] or a
%               T-by-0 matrix.
%                 Data Types - double
%
%  Optional input arguments:
%
%    bsb :      Indices of observations used in the fit. Vector.
%               If bsb is empty (default), all finite observations in y are
%               used. All y(bsb) values must be finite. The minimum number
%               of selected observations is size(N,2) plus the nullity of
%               betaPenalty, except
%               when this number is zero, in which case at least one
%               observation is required.
%                 Example - 'bsb',(1:80)'
%                 Data Types - double
%
%    X :        Additional linear components. Matrix or empty.
%               T-by-q matrix entering the fitted signal as m+X*beta.
%               Directions penalized by betaPenalty need not be identified
%               by the data. The unpenalized directions must be identified
%               jointly with null(K). The default is empty.
%                 Example - 'X',[sin(2*pi*t/12) cos(2*pi*t/12)]
%                 Data Types - double | empty
%
% betaPenalty : q-by-q positive-semidefinite matrix P_beta applied to
%               beta.
%               Its scale is included in the matrix. The contribution to
%               the objective is beta'*P_beta*beta. If empty (default),
%               beta is unpenalized.
%                 Example - 'betaPenalty',Pseason
%                 Data Types - double or sparse double
%
% betaPenaltyNullBasis : Optional basis for the full null space of
%               betaPenalty. These coefficient directions must be
%               identified by the data. If omitted, penfilterFS computes
%               the basis. Supplying a precomputed basis avoids repeated
%               null-space decompositions in repeated fitting calls.
%                 Example - 'betaPenaltyNullBasis',Nseason
%                 Data Types - double
%
% lambda :      Smoothing parameter. Positive scalar. The criterion is
%               RSS + lambda*m'*K*m + beta'*betaPenalty*beta. The default
%               value is 1. No automatic frequency scaling is used because
%               the scale of a generic K is user-defined.
%                 Example - 'lambda',1600
%                 Data Types - double
%
% nocheck :     Skip structural checks on K and N. Boolean.
%               If nocheck=false (default), penfilterFS checks symmetry and
%               approximate positive semidefiniteness of K, full column
%               rank of N, and verifies that K*N is numerically zero.
%               Set nocheck=true only after K and N have been validated,
%               for example inside repeated LTSts concentration steps.
%               Identification on the selected subset is always checked.
%                 Example - 'nocheck',true
%                 Data Types - logical
%
% plots :       Plot observed values and fitted components. Boolean.
%               If plots=true, the observed series and penalized trend are
%               shown. When X is nonempty, the complete fitted signal is
%               also shown. The default is false.
%                 Example - 'plots',true
%                 Data Types - logical
%
%  Output:
%
%    out :      Structure containing the following fields:
%
%          out.mhat      = T-by-1 estimated penalized trend m.
%          out.beta      = q-by-1 coefficient vector associated with X.
%                          Empty if X is empty.
%          out.yhat      = T-by-1 complete fitted signal m+X*beta.
%          out.residuals = T-by-1 residual vector y-out.yhat.
%          out.bsb       = indices used in the fit.
%          out.lambda    = smoothing parameter used in the fit.
%          out.rss       = residual sum of squares computed on bsb.
%          out.roughness = value m'*K*m.
%          out.penalty   = value lambda*m'*K*m.
%          out.betaPenaltyValue = value beta'*betaPenalty*beta.
%          out.objective = out.rss+out.penalty+out.betaPenaltyValue.
%          out.nullity   = number of columns of N.
%          out.rankAug   = rank of [N(bsb,:) X(bsb,:)*B0].
%          out.class     = 'penfilterFS'.
%
%    Aready :   Optional second output containing the factorization of
%               W+lambda*K. It is intended for model-specific wrappers
%               that need additional calculations based on the same fit.
%
%  More About:
%
% The matrix
%
%   A = W + lambda*K,       W=S'*S,
%
% is common to both the trend-only and joint trend-regression problems.
% If X is present, m is profiled out exactly. With Xb=X(bsb,:),
%
%   P = W-W*A^{-1}*W,
%
% and beta solves
%
%   (X'*P*X+betaPenalty)*beta = X'*P*y.
%
% The implementation never forms P explicitly. This is the same profiling
% algebra used in hpfilterFS, but with the HP matrix D2'*D2 replaced by a
% generic fixed positive semidefinite K.
%
% The supplied N must span the complete null space of K. penfilterFS checks
% that the supplied columns lie in null(K), but it does not compute the
% nullity of a large sparse K in order to verify completeness. Supplying an
% incomplete N can therefore lead to a singular A for some subsets; this is
% detected when the penalized system is factorized.
%
% See also: hpfilterFS, LTSts
%
% References:
%
% Green, P. J. and Silverman, B. W. (1994), "Nonparametric Regression and
% Generalized Linear Models: A Roughness Penalty Approach", Chapman and
% Hall, London.
%
% Wahba, G. (1990), "Spline Models for Observational Data", SIAM,
% Philadelphia.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('penfilterFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%{
   %% HP penalty: check equivalence with hpfilterFS in fitting-only mode.
   rng(1000)
   T=120;
   t=(1:T)';
   e=ones(T,1);
   D2=spdiags([e -2*e e],0:2,T-2,T);
   K=D2'*D2;
   N=[ones(T,1) t];
   X=[sin(2*pi*t/12) cos(2*pi*t/12)];
   y=0.02*t+X*[3;-2]+randn(T,1);
   bsb=(1:100)';

   outPen=penfilterFS(y,K,N,'bsb',bsb,'X',X,'lambda',1600);
   outHP=hpfilterFS(y,'bsb',bsb,'X',X,'lambda',1600,...
       'fitonly',true);

   assert(max(abs(outPen.mhat-outHP.mhat))<1e-9)
   assert(max(abs(outPen.yhat-outHP.yhat))<1e-9)
   assert(abs(outPen.objective-outHP.objective)<1e-7)
%}

%{
   %% Third-difference penalty.
   % The null space is formed by polynomials of degree at most two.
   rng(12345)
   T=100;
   t=(1:T)';
   e=ones(T,1);
   D3=spdiags([e -3*e 3*e -e],0:3,T-3,T);
   K3=D3'*D3;
   x=(t-mean(t))/T;
   N3=[ones(T,1) x x.^2];
   y=2+0.03*t+4*sin(2*pi*t/T)+randn(T,1);
   out=penfilterFS(y,K3,N3,'lambda',1e5,'plots',true);
%}

%{
   %% Positive definite penalty.
   rng(123)
   T=80;
   K=speye(T);
   N=zeros(T,0);
   y=randn(T,1);
   out=penfilterFS(y,K,N,'lambda',2);
%}

%% Beginning of code

if nargin<3
    error('FSDA:penfilterFS:missingInputs', ...
        'Required input arguments y, K and N are missing.')
end

% Default options.
bsb=[];
X=[];
lambda=1;
betaPenalty=[];
betaPenaltyNullBasis=[];
nocheck=false;
plots=false;

if ~isempty(varargin)
    if rem(length(varargin),2)~=0
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'Optional arguments must be supplied as name/value pairs.')
    end

    options=struct('bsb',bsb,'X',X,'lambda',lambda, ...
        'betaPenalty',betaPenalty, ...
        'betaPenaltyNullBasis',betaPenaltyNullBasis, ...
        'nocheck',nocheck,'plots',plots);
    names=fieldnames(options);

    for i=1:2:length(varargin)
        name=varargin{i};
        if ~(ischar(name) || (isstring(name) && isscalar(name)))
            error('FSDA:penfilterFS:WrongInputOpt', ...
                'Optional argument names must be strings or character vectors.')
        end
        name=char(name);
        if ~any(strcmp(name,names))
            error('FSDA:penfilterFS:WrongInputOpt', ...
                'Unknown option ''%s''.',name)
        end
        options.(name)=varargin{i+1};
    end

    bsb=options.bsb;
    X=options.X;
    lambda=options.lambda;
    betaPenalty=options.betaPenalty;
    betaPenaltyNullBasis=options.betaPenaltyNullBasis;
    nocheck=options.nocheck;
    plots=options.plots;
end

% Convert timetable input. The generic penalty has no periodicity to infer.
if istimetable(y)
    if width(y)<1
        error('FSDA:penfilterFS:WrongInput', ...
            'Input timetable y must contain at least one variable.')
    end
    rowTimes=y.Properties.RowTimes;
    y=y{:,1};
    isTT=true;
else
    y=y(:);
    rowTimes=[];
    isTT=false;
end

T=length(y);
if T<1
    error('FSDA:penfilterFS:WrongInput', ...
        'Input y must contain at least one observation.')
end
seq=(1:T)';

% Check scalar logical options.
if ~(islogical(nocheck) && isscalar(nocheck))
    if isnumeric(nocheck) && isscalar(nocheck) && ismember(nocheck,[0 1])
        nocheck=logical(nocheck);
    else
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'nocheck must be a logical scalar.')
    end
end
if ~(islogical(plots) && isscalar(plots))
    if isnumeric(plots) && isscalar(plots) && ismember(plots,[0 1])
        plots=logical(plots);
    else
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'plots must be a logical scalar.')
    end
end

% Check smoothing parameter.
if ~(isnumeric(lambda) && isreal(lambda) && isscalar(lambda) && ...
        isfinite(lambda) && lambda>0)
    error('FSDA:penfilterFS:WrongInputOpt', ...
        'lambda must be a positive finite scalar.')
end

% Validate K.
if ~(isnumeric(K) && isreal(K) && ismatrix(K) && ...
        size(K,1)==T && size(K,2)==T)
    error('FSDA:penfilterFS:WrongInput', ...
        'K must be a real numeric T-by-T matrix.')
end
if issparse(K)
    kvals=nonzeros(K);
else
    kvals=K(:);
end
if any(~isfinite(kvals))
    error('FSDA:penfilterFS:WrongInput', ...
        'All stored entries of K must be finite.')
end

% Validate N and convert an empty input to a T-by-0 matrix.
if isempty(N)
    N=zeros(T,0);
else
    if ~(isnumeric(N) && isreal(N) && ismatrix(N) && size(N,1)==T)
        error('FSDA:penfilterFS:WrongInput', ...
            'N must be a real numeric matrix with T rows.')
    end
    if any(~isfinite(N(:)))
        error('FSDA:penfilterFS:WrongInput', ...
            'All entries of N must be finite.')
    end
end
r=size(N,2);

% Structural checks on K and N. These can be skipped in repeated calls once
% a fixed penalty has been validated.
if ~nocheck
    nK=max(1,norm(K,'fro'));
    symerr=norm(K-K','fro')/nK;
    if symerr>1e-10
        error('FSDA:penfilterFS:NonSymmetricPenalty', ...
            'K must be symmetric. Relative symmetry error is %g.',symerr)
    end
end

% Use the symmetric part exactly. For a valid input this changes K only at
% round-off level and ensures a symmetric penalized system.
K=(K+K')/2;

if ~nocheck
    % Approximate PSD check. A tiny diagonal perturbation permits exact
    % semidefiniteness but still detects materially negative eigenvalues.
    scaleK=max(1,norm(K,1));
    tolPSD=1e-10*scaleK;
    [~,pPSD]=chol(K+tolPSD*speye(T));
    if pPSD~=0
        error('FSDA:penfilterFS:NonPSDPenalty', ...
            'K must be positive semidefinite.')
    end

    if r>0
        if rank(full(N))<r
            error('FSDA:penfilterFS:WrongNullSpace', ...
                'The columns of N must be linearly independent.')
        end
        % Orthonormalize before checking the null-space residual so that
        % the tolerance is invariant to the scale of the supplied basis.
        [N,~]=qr(N,0);
        nullerr=norm(K*N,'fro')/max(1,norm(K,'fro')*sqrt(r));
        if nullerr>1e-8
            error('FSDA:penfilterFS:WrongNullSpace', ...
                ['The columns of N must lie in null(K). Relative ' ...
                'null-space residual is %g.'],nullerr)
        end
    end
else
    % Even when checks are skipped, orthonormalization is inexpensive
    % because r is normally small and improves subset rank diagnostics.
    if r>0
        [N,~]=qr(N,0);
    end
end

% Validate bsb. If omitted, use all finite observations.
if isempty(bsb)
    bsb=find(isfinite(y));
else
    if ~isnumeric(bsb) || ~isvector(bsb) || any(~isfinite(bsb)) || ...
            any(bsb~=fix(bsb))
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'bsb must contain finite integer indices.')
    end
    bsb=unique(bsb(:));
    if any(bsb<1) || any(bsb>T)
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'bsb contains indices outside 1..T.')
    end
    if any(~isfinite(y(bsb)))
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'All values y(bsb) must be finite.')
    end
end
nbsb=length(bsb);
if nbsb==0
    error('FSDA:penfilterFS:WrongInputOpt', ...
        'At least one finite observation must belong to bsb.')
end

% Validate optional linear design matrix.
if isempty(X)
    X=zeros(T,0);
else
    if ~(isnumeric(X) && isreal(X) && ismatrix(X) && size(X,1)==T)
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'X must be a real numeric matrix with T rows.')
    end
    if any(~isfinite(X(:)))
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'All entries of X must be finite.')
    end
end
q=size(X,2);
Xb=X(bsb,:);

% Validate the coefficient penalty and its unpenalized directions.
[betaPenalty,betaPenaltyNullBasis]=validateBetaPenalty( ...
    betaPenalty,betaPenaltyNullBasis,q);
qnull=size(betaPenaltyNullBasis,2);

% Conditional identification on the selected subset.
p0=r+qnull;
if nbsb<p0
    error('FSDA:penfilterFS:RankDeficient', ...
        ['The number of observations in bsb must be at least ' ...
        'size(N,2)+the nullity of betaPenalty.'])
end
if p0>0
    Aid=[N(bsb,:) Xb*betaPenaltyNullBasis];
    rankAug=rank(full(Aid));
    if rankAug<p0
        error('FSDA:penfilterFS:RankDeficient', ...
            ['The null directions of K and betaPenalty are not ' ...
            'identifiable on bsb. The augmented matrix ' ...
            '[N(bsb,:) X(bsb,:)*betaPenaltyNullBasis] must have ' ...
            'full column rank.'])
    end
else
    rankAug=0;
end

% Selection via diagonal weights W=S'*S.
w=zeros(T,1);
w(bsb)=1;
W=spdiags(w,0,T,T);

% Penalized system. Under the stated identification condition and when N
% spans null(K), A is positive definite.
A=W+lambda*K;
try
    Aready=decomposition(A,'chol');
catch ME
    error('FSDA:penfilterFS:SingularPenalizedSystem', ...
        ['Unable to factorize W+lambda*K as positive definite. Check ' ...
            'that K is positive semidefinite, that N spans the complete ' ...
            'null space of K, and that the unpenalized coefficient ' ...
            'directions are identifiable on bsb. ' ...
        'Original message: %s'],ME.message)
end

% Never form W*y directly because y can contain NaNs outside bsb.
Wy=zeros(T,1);
Wy(bsb)=y(bsb);
AinvWy=Aready\Wy;

if q==0
    beta=zeros(0,1);
    mhat=AinvWy;
else
    % Profile m out of the joint penalized least-squares criterion.
    WX=zeros(T,q);
    WX(bsb,:)=Xb;
    AinvWX=Aready\WX;

    Ry=y(bsb)-AinvWy(bsb);
    RX=Xb-AinvWX(bsb,:);

    G=Xb'*RX+betaPenalty;
    G=(G+G')/2;
    g=Xb'*Ry;

    % The augmented-rank check above addresses structural
    % non-identifiability. Backslash is retained so MATLAB can diagnose
    % severe numerical ill-conditioning in the usual way.
    beta=G\g;
    mhat=AinvWy-AinvWX*beta;
end

% Complete fitted signal and residuals.
yhat=mhat+X*beta;
residuals=y-yhat;

% Objective components.
r_obs=y(bsb)-yhat(bsb);
rss=full(r_obs'*r_obs);
Km=K*mhat;
roughness=full(mhat'*Km);

% Tiny negative values can arise solely from round-off when K is PSD.
if roughness<0 && abs(roughness)<=1e-10*max(1,norm(mhat)^2*norm(K,1))
    roughness=0;
end
if roughness<0
    error('FSDA:penfilterFS:NegativePenalty', ...
        ['The computed quadratic form m''*K*m is negative. Check that K ' ...
        'is positive semidefinite.'])
end
penalty=lambda*roughness;
if q>0
    betaPenaltyValue=full(beta'*betaPenalty*beta);
else
    betaPenaltyValue=0;
end
if betaPenaltyValue<0 && abs(betaPenaltyValue)<= ...
        1e-10*max(1,norm(beta)^2*norm(betaPenalty,1))
    betaPenaltyValue=0;
end
if betaPenaltyValue<0
    error('FSDA:penfilterFS:NegativeCoefficientPenalty', ...
        'The computed coefficient penalty is negative.')
end
objective=rss+penalty+betaPenaltyValue;

% Store output.
out=struct;
out.mhat=mhat;
out.beta=beta;
out.yhat=yhat;
out.residuals=residuals;
out.bsb=bsb;
out.lambda=lambda;
out.rss=rss;
out.roughness=roughness;
out.penalty=penalty;
out.betaPenaltyValue=betaPenaltyValue;
out.objective=objective;
out.nullity=r;
out.rankAug=rankAug;
out.class='penfilterFS';

% Optional plot.
if plots
    if isTT
        xplot=rowTimes;
        xlab='Time';
    else
        xplot=seq;
        xlab='Observation';
    end

    figure('Tag','penfilterFS');
    plot(xplot,y,'.','DisplayName','Observed');
    hold on
    plot(xplot,mhat,'LineWidth',1.5,'DisplayName','Penalized trend');
    if q>0
        plot(xplot,yhat,'--','LineWidth',1.2,'DisplayName','Fitted signal');
    end
    hold off
    xlabel(xlab)
    ylabel('y')
    legend('Location','best')
    box on
    grid on
end

end

function [Pbeta,Nbeta] = validateBetaPenalty(Pbeta,Nbeta,q)
% Validate the coefficient-block penalty and its complete null basis.

if ~(isnumeric(Pbeta) && isreal(Pbeta) && ismatrix(Pbeta))
    error('FSDA:penfilterFS:WrongInputOpt', ...
        'betaPenalty must be a real numeric matrix or empty.')
end
if isempty(Pbeta)
    Pbeta=sparse(q,q);
else
    if ~isequal(size(Pbeta),[q q]) || ...
            any(~isfinite(nonzeros(Pbeta)))
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'betaPenalty must be a finite real q-by-q numeric matrix.')
    end
    Pbeta=sparse(double(Pbeta));
    penaltyScale=norm(Pbeta,'fro');
    symmetryTolerance=1e-10*max(penaltyScale,realmin);
    if norm(Pbeta-Pbeta','fro')>symmetryTolerance
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'betaPenalty must be symmetric.')
    end
    Pbeta=(Pbeta+Pbeta')/2;

    % Sparse Cholesky with a scale-relative shift checks PSD without
    % forming a dense eigendecomposition.
    penaltyScale=norm(Pbeta,1);
    if q>0 && penaltyScale>0
        psdTolerance=q*eps(penaltyScale);
        [~,cholFlag]=chol(Pbeta+psdTolerance*speye(q));
        if cholFlag~=0
            error('FSDA:penfilterFS:WrongInputOpt', ...
                'betaPenalty must be positive semidefinite.')
        end
    end
end

if ~(isnumeric(Nbeta) && isreal(Nbeta) && ismatrix(Nbeta))
    error('FSDA:penfilterFS:WrongInputOpt', ...
        'betaPenaltyNullBasis must be a real numeric matrix or empty.')
end
nullBasisSupplied=~isempty(Nbeta);
if isempty(Nbeta)
    if q==0
        Nbeta=zeros(0,0);
    elseif nnz(Pbeta)==0
        Nbeta=speye(q);
    else
        Nbeta=null(full(Pbeta));
    end
elseif ~(size(Nbeta,1)==q && size(Nbeta,2)<=q && ...
        all(isfinite(Nbeta(:))))
    error('FSDA:penfilterFS:WrongInputOpt', ...
        'betaPenaltyNullBasis must be a finite real matrix with q rows.')
end

qnull=size(Nbeta,2);
if nullBasisSupplied && qnull>0
    Ncheck=double(full(Nbeta));
    if rank(Ncheck)<qnull
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'The columns of betaPenaltyNullBasis must be linearly independent.')
    end
    basisColumnNorms=sqrt(sum(Ncheck.^2,1));
    if any(~isfinite(basisColumnNorms)) || any(basisColumnNorms==0)
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'The columns of betaPenaltyNullBasis must be nonzero and finite.')
    end
    Ncheck=bsxfun(@rdivide,Ncheck,basisColumnNorms);
    [Ncheck,~]=qr(Ncheck,0);

    penaltyScale=norm(Pbeta,'fro');
    if penaltyScale>0
        nullerr=norm((Pbeta/penaltyScale)*Ncheck,'fro');
        if nullerr>1e-8*sqrt(qnull)
            error('FSDA:penfilterFS:WrongInputOpt', ...
                'betaPenaltyNullBasis must lie in the null space of betaPenalty.')
        end

        completionMatrix=full(Pbeta/penaltyScale)+Ncheck*Ncheck';
        completionMatrix=(completionMatrix+completionMatrix')/2;
        [~,completionFlag]=chol(completionMatrix);
        if completionFlag~=0 && rank(completionMatrix)<q
            error('FSDA:penfilterFS:WrongInputOpt', ...
                'betaPenaltyNullBasis must span the full null space of betaPenalty.')
        end
    elseif qnull~=q
        error('FSDA:penfilterFS:WrongInputOpt', ...
            'betaPenaltyNullBasis must span the full null space of betaPenalty.')
    end

    % Use an orthonormal version downstream so subset rank checks do not
    % depend on the arbitrary scale or conditioning of the supplied basis.
    Nbeta=Ncheck;
end
end
