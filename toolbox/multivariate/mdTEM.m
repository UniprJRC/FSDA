function out = mdTEM(Y, varargin)
%mdTEM EM algorithm with trimming (TEM) for data with missing values.
%
%<a href="matlab: docsearchFS('mdTEM')">Link to the help function</a>
%
% The algorithm:
%  - At each iteration compute adjusted partial Mahalanobis distances;
%  - Rank them and set weights w_i = 1 for the lowest n*(1-alpha) rows,
%    else 0;
%  - Run E-step and M-step using these weights;
%  - Apply a consistency correction for the truncation induced by
%    trimming (see input option consistencyfactor);
%  - Repeat until convergence or maxiter.
%
% Required input arguments:
%
% Y :           Input data. Matrix. n x p data matrix; n observations and
%               p variables possibly with missing values (NaN's). Rows of
%               Y represent observations, and columns represent variables.
%               Data Types - single | double
%
% Optional input arguments:
%
%       alpha   : proportion to trim. Real number in the interval [0 0.5]
%                 or empty value. At each iteration the lowest
%                 n*(1-alpha) adjusted distances are retained. If alpha is
%                 empty the default value is 0.5.
%                 Example - 'alpha',0.1
%                 Data Types - single | double
%
%       mus     : initial mean. p x 1 vector or empty double. If empty
%                 (default), column means ignoring NaN values are used.
%                 Example - 'mus',[]
%                 Data Types - single | double
%
%       sigs    : initial covariance/scatter matrix. p x p matrix or empty
%                 double. If empty, a covariance matrix computed after
%                 conditional mean replacement by the initial column means
%                 is used.
%                 Example - 'sigs',eye(p)
%                 Data Types - single | double
%
%       maxiter : maximum number of iterations. Positive integer. The
%                 default value is 100.
%                 Example - 'maxiter',50
%                 Data Types - single | double
%
%       tol     : tolerance for convergence. Positive real number. The
%                 default value is 1e-5.
%                 Example - 'tol',1e-10
%                 Data Types - single | double
%
%   tol_sigma   : use tolerance for both location and scatter. Boolean. If
%                 true, both location and scatter changes enter the
%                 convergence criterion (default true).
%                 Example - 'tol_sigma',false
%                 Data Types - logical
%
%   method      : method used to rescale the partial distances. String
%                 scalar or character vector. Possible values are:
%
%                 'pri'      = principled EM rescaling (default),
%                              d2_partial + (p-pobs).
%                 'expScale' = expectation scaling,
%                              d2_partial*(p/pobs).
%                 'zMap'     = standardization mapping,
%                              p + sqrt(2*p)*(d2_partial-pobs)/sqrt(2*pobs).
%                 'detMap'   = determinant-based rescaling.
%                 'chiMap'   = chi-square quantile mapping.
%                 'betaMap'  = Beta quantile mapping.
%                 'impMD'    = MD on EM-imputed data.
%                 Example - 'method','chiMap'
%                 Data Types - string scalar | char vector
%
% consistencyfactor : treatment of the truncation bias of the scatter
%                 estimate. Character vector or string scalar. Possible
%                 values are:
%
%                 'global'   = (default) one Gaussian Tallis factor at the
%                              full dimension p and global retained
%                              fraction h/n, applied to the whole scatter.
%
%                 'pattern'  = Gaussian pattern-wise Tallis correction.
%                              Within pattern g, trimming is converted to a
%                              cutoff a_g on the partial squared-distance
%                              scale and
%                                 k_g=F_{p_g+2}(a_g)/F_{p_g}(a_g).
%                              Only the data-driven conditional second-
%                              moment part is divided by k_g; the
%                              conditional-variance term is not corrected.
%
%                 'adaptive' = distribution-adaptive pattern-wise radial
%                              correction. For each observed pattern g the
%                              factor is estimated from complete reference
%                              observations as
%                                 k_g = mean(Q_g | Q_g<=a_g)/p_g,
%                              where Q_g is computed with muref and
%                              sigmaref projected on the variables observed
%                              in pattern g. Under elliptical MCAR, this
%                              makes TEM target the same affine-equivariant
%                              scatter functional as the reference fit. The
%                              reference scatter need not equal the
%                              covariance matrix. The current implementation
%                              supports 'pri', 'expScale', 'zMap', 'chiMap',
%                              'betaMap' and 'impMD' with this option.
%
%                 'weighted' = one scalar Gaussian correction equal to the
%                              information-weighted average of the exact
%                              pattern-wise Tallis factors.
%
%                 'none'     = no correction.
%
%                 Example - 'consistencyfactor','adaptive'
%                 Data Types - char | string
%
%   Yref         : complete reference data used only when
%                  consistencyfactor='adaptive'. nref x p finite matrix.
%                  If empty, the complete rows of Y are used.
%                  Example - 'Yref',Ycc
%                  Data Types - single | double
%
%   muref        : reference location used by the adaptive factors. p x 1
%                  vector or empty. If empty, the mean of Yref is used.
%                  When mdTEM is initialized by a robust complete-case fit,
%                  muref should normally be set to that fit's location.
%                  Example - 'muref',RAW.loc(:)
%                  Data Types - single | double
%
%   sigmaref     : reference scatter used by the adaptive factors. p x p
%                  positive definite matrix or empty. If empty, cov(Yref,1)
%                  is used. When mdTEM is initialized by a robust complete-
%                  case fit, sigmaref should normally be set to that fit's
%                  scatter matrix.
%                  Example - 'sigmaref',RAW.cov
%                  Data Types - single | double
%
% adaptivepool   : pool adaptive projected radii across patterns having the
%                  same observed dimension. Boolean. Under ellipticity this
%                  is asymptotically valid and can stabilize small samples.
%                  Default is false.
%                  Example - 'adaptivepool',true
%                  Data Types - logical
%
% adaptiveminref : minimum number of reference radii inside a pattern
%                  cutoff required by the adaptive correction. Positive
%                  integer. Default is 20. If a retained pattern does not
%                  meet this requirement, mdTEM stops with an informative
%                  error rather than silently replacing the factor.
%                  Example - 'adaptiveminref',10
%                  Data Types - single | double
%
%   condmeanimp  : also return conditional-mean imputed values. Boolean.
%                  Default is false.
%                  Example - 'condmeanimp',true
%                  Data Types - logical
%
%   stochimp     : also return stochastic imputed values. Boolean. Default
%                  is false.
%                  Example - 'stochimp',true
%                  Data Types - logical
%
%   storeobj     : store the trimmed sum of the smallest adjusted distances
%                  in each iteration. Boolean. Default is true.
%                  Example - 'storeobj',false
%                  Data Types - logical
%
% Output:
%
% out : structure which contains the following fields
%       out.loc       = final location estimate.
%       out.cov       = final scatter estimate.
%       out.iter      = number of iterations to convergence.
%       out.weights   = n x 1 vector of final 0/1 trimming weights.
%       out.Yimp      = empty or conditional-mean imputed data.
%       out.stochYimp = empty or stochastically imputed data.
%       out.obj       = empty or objective values by iteration.
%       out.kfactor   = scalar factor actually applied to the whole scatter
%                       for 'global', 'weighted' and 'none', or the
%                       information-weighted average of the pattern factors
%                       for 'pattern' and 'adaptive'.
%       out.kinfo     = table of final pattern-wise quantities. For
%                       'pattern' it contains pobs, nkept, athr, gammag and
%                       kg. For 'adaptive' it additionally contains nref,
%                       nrefkept and factorType. Empty for 'global'/'none'.
%       out.consistencyfactor = consistency-factor option used.
%
% More About:
%
% The Gaussian pattern correction is exact under a Gaussian MCAR model at
% the population target. The adaptive correction replaces the Gaussian
% radial truncated moment by its empirical complete-case analogue. Under
% an elliptical MCAR model and a reference fit converging to S0=tau*Sigma,
% the adaptive factor converges to the radial factor appropriate to S0, so
% the TEM and reference scatters have the same population target.
%
% The adaptive reference location and scatter are frozen during the TEM
% iteration. This is essential: they define the complete-case scatter
% functional that TEM is asked to reproduce. The TEM iterates themselves
% continue to determine distances, retained rows and the common threshold.
%
% Patterns with very few retained TEM observations can make the Gaussian
% pattern factor unstable; the existing Gaussian stabilization is retained.
% Adaptive factors are not silently shrunk. Use adaptivepool=true or reduce
% adaptiveminref if the complete reference sample is too small.
%
% See also: mdEM, mdImputeCondMean, mdPartialMD, mdPartialMD2full, mcd
%
% References:
%
% Little, R. J. A., & Rubin, D. B. (2019). Statistical Analysis with
% Missing Data (3rd ed.). Hoboken, NJ: John Wiley & Sons.
% Tallis, G. M. (1963). Elliptical and radial truncation in normal samples.
% Annals of Mathematical Statistics, 34, 940-944.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdTEM')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Adaptive radial correction with a robust complete-case reference.
    rng(7)
    p=5; n=1000; rho=0.5; alpha=0.25;
    Sigma=(1-rho)*eye(p)+rho*ones(p);
    Yfull=randn(n,p)*chol(Sigma);
    miss=rand(n,p)<0.12;
    Y=Yfull;
    Y(miss)=NaN;
    cc=all(~isnan(Y),2);
    [RAW,~]=mcd(Y(cc,:),'bdp',alpha,'msg',0,'plots',0);
    out=mdTEM(Y,'alpha',alpha,'mus',RAW.loc(:),'sigs',RAW.cov, ...
        'consistencyfactor','adaptive','Yref',Y(cc,:), ...
        'muref',RAW.loc(:),'sigmaref',RAW.cov);
    disp(out.kinfo)
%}

%% Beginning of code
alpha=0.5;
mus=[];
sigs=[];
maxiter=100;
tol=1e-5;
tol_sigma=true;
condmeanimp=false;
stochimp=false;
storeobj=true;
method='pri';
consistencyfactor='global';
Yref=[];
muref=[];
sigmaref=[];
adaptivepool=false;
adaptiveminref=20;

if nargin>1
    options=struct('storeobj',storeobj,'alpha',alpha,'mus',mus,'sigs',sigs, ...
        'maxiter',maxiter,'tol',tol,'tol_sigma',tol_sigma, ...
        'condmeanimp',condmeanimp,'stochimp',stochimp,'method',method, ...
        'consistencyfactor',consistencyfactor,'Yref',Yref,'muref',muref, ...
        'sigmaref',sigmaref,'adaptivepool',adaptivepool, ...
        'adaptiveminref',adaptiveminref);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdTEM:WrongInputOpt', ...
                ['Number of supplied options is invalid. Probably values ' ...
                'for some parameters are missing.']);
        end
        aux.chkoptions(options,UserOptions)
    end

    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    alpha=options.alpha;
    mus=options.mus;
    sigs=options.sigs;
    maxiter=options.maxiter;
    tol=options.tol;
    tol_sigma=options.tol_sigma;
    condmeanimp=options.condmeanimp;
    stochimp=options.stochimp;
    method=string(options.method);
    storeobj=options.storeobj;
    consistencyfactor=char(string(options.consistencyfactor));
    Yref=options.Yref;
    muref=options.muref;
    sigmaref=options.sigmaref;
    adaptivepool=logical(options.adaptivepool);
    adaptiveminref=options.adaptiveminref;
end

validCF={'global','pattern','adaptive','weighted','none'};
if ~any(strcmpi(consistencyfactor,validCF))
    error('FSDA:mdTEM:WrongInputOpt', ...
        ['Option ''consistencyfactor'' must be one of ''global'', ', ...
         '''pattern'', ''adaptive'', ''weighted'' or ''none''.']);
end
consistencyfactor=lower(consistencyfactor);

if ~isscalar(alpha) || ~isfinite(alpha) || alpha < 0 || alpha > 0.5
    error('FSDA:mdTEM:WrongAlpha','alpha must be a scalar in [0,0.5].');
end
if ~isscalar(maxiter) || maxiter < 1 || maxiter ~= floor(maxiter)
    error('FSDA:mdTEM:WrongMaxiter','maxiter must be a positive integer.');
end
if ~isscalar(tol) || ~isfinite(tol) || tol <= 0
    error('FSDA:mdTEM:WrongTol','tol must be a positive finite scalar.');
end
if ~isscalar(adaptiveminref) || adaptiveminref < 1 || ...
        adaptiveminref ~= floor(adaptiveminref)
    error('FSDA:mdTEM:WrongAdaptiveMinref', ...
        'adaptiveminref must be a positive integer.');
end

if strcmp(consistencyfactor,'adaptive') && method=="detMap"
    error('FSDA:mdTEM:AdaptiveDetMapNotSupported', ...
        ['The adaptive correction is not yet implemented for ''detMap'' ' ...
        'because that adjustment depends on the iterated scatter matrix.']);
end

if storeobj==true
    obj=zeros(maxiter,1);
else
    obj=[];
end
[n,p]=size(Y);

% initialize mus and sigs if not provided
if isempty(mus)
    mus=mean(Y,1,"omitmissing")';
else
    mus=mus(:);
end
if isempty(sigs)
    X0=Y;
    for j=1:p
        miss=isnan(X0(:,j));
        X0(miss,j)=mus(j);
    end
    sigs=cov(X0,1);
end
sigs=(sigs+sigs')/2;

if numel(mus)~=p || ~isequal(size(sigs),[p p])
    error('FSDA:mdTEM:WrongInitialFit', ...
        'mus and sigs have incompatible dimensions.');
end

% missingness mask and patterns are intrinsic to the data
nanY=isnan(Y);

% Set and freeze the reference distribution used by adaptive factors.
if strcmp(consistencyfactor,'adaptive')
    if isempty(Yref)
        complete=all(~nanY,2);
        Yref=Y(complete,:);
    end
    if isempty(Yref) || size(Yref,2)~=p || any(~isfinite(Yref(:)))
        error('FSDA:mdTEM:WrongAdaptiveYref', ...
            ['Yref must be a nonempty finite complete nref-by-p matrix ' ...
            'when consistencyfactor is ''adaptive''.']);
    end
    Yref=double(Yref);

    if isempty(muref)
        muref=mean(Yref,1)';
    else
        muref=muref(:);
    end
    if numel(muref)~=p || any(~isfinite(muref))
        error('FSDA:mdTEM:WrongAdaptiveMuref', ...
            'muref must be a finite p-by-1 vector.');
    end

    if isempty(sigmaref)
        sigmaref=cov(Yref,1);
    end
    if ~isequal(size(sigmaref),[p p]) || any(~isfinite(sigmaref(:)))
        error('FSDA:mdTEM:WrongAdaptiveSigmaref', ...
            'sigmaref must be a finite p-by-p matrix.');
    end
    sigmaref=(sigmaref+sigmaref')/2;
    [~,flagref]=chol(sigmaref);
    if flagref~=0
        error('FSDA:mdTEM:NonPDAdaptiveSigmaref', ...
            'sigmaref must be positive definite.');
    end
end

dif=Inf;
iter=0;
keep_count=max(0,floor(n*(1-alpha)));
w=zeros(n,1);
kinfo=[];
kfactor=1;
cthr=NaN;
d2_adj=NaN(n,1);

while (dif>tol) && (iter<maxiter)
    iter=iter+1;
    mus_old=mus;
    sigs_old=sigs;

    if method=="impMD"
        Yimp=mdImputeCondMean(Y,mus,sigs);
        d2_adj=mahalFS(Yimp,mus',sigs);
        poss=sum(~nanY,2);
    elseif method=="detMap"
        [d2,poss]=mdPartialMD(Y,mus,sigs);
        d2_adj=mdPartialMD2full(d2,p,poss,'method',method,'Y',Y,'Sigma',sigs);
    else
        [d2,poss]=mdPartialMD(Y,mus,sigs);
        d2_adj=mdPartialMD2full(d2,p,poss,'method',method);
    end

    nan_mask=isnan(d2_adj);
    [~,idx_sorted]=sort(d2_adj,'ascend','MissingPlacement','last');
    keep_idx=idx_sorted(1:min(keep_count,sum(~nan_mask)));
    if isempty(keep_idx)
        error('FSDA:mdTEM:NoRetainedRows', ...
            'No row has a finite adjusted distance at the current iteration.');
    end

    w=zeros(n,1);
    w(keep_idx)=1;
    mm=sum(w);
    cthr=max(d2_adj(keep_idx));

    switch consistencyfactor
        case 'pattern'
            kinfo=localPatternFactors(nanY,w,poss,cthr,p,n,sigs,method);
            [mus,sigs,kfactor]=localCorrectedStep(Y,nanY,w,mus,sigs,kinfo);

        case 'adaptive'
            kinfo=localAdaptiveFactors(nanY,w,cthr,p,n,sigs,method, ...
                Yref,muref,sigmaref,adaptivepool,adaptiveminref);
            bad=~isfinite(kinfo.kg) | kinfo.kg<=0;
            if any(bad)
                ibad=find(bad,1);
                error('FSDA:mdTEM:FewAdaptiveReferencePoints', ...
                    ['Adaptive factor unavailable for a retained pattern ' ...
                    '(pobs=%d, nkept=%d, nrefkept=%d). Increase the ' ...
                    'complete reference sample, set adaptivepool=true, or ' ...
                    'reduce adaptiveminref.'], ...
                    kinfo.pobs(ibad),kinfo.nkept(ibad),kinfo.nrefkept(ibad));
            end
            [mus,sigs,kfactor]=localCorrectedStep(Y,nanY,w,mus,sigs,kinfo);

        otherwise
            [T1,T2]=aux.NAcompute_expected_stats(Y,mus,sigs,w);
            [mus,sigs]=aux.NAmaximization_step(T1,T2,w);

            switch consistencyfactor
                case 'global'
                    a=chi2inv(mm/n,p);
                    kfactor=(n./mm).*chi2cdf(a,p+2);
                    kinfo=[];
                case 'weighted'
                    kinfo=localPatternFactors(nanY,w,poss,cthr,p,n,sigs,method);
                    kfactor=localWeightedFactor(kinfo);
                case 'none'
                    kfactor=1;
                    kinfo=[];
            end
            if isfinite(kfactor) && kfactor>0
                sigs=sigs/kfactor;
            end
    end

    if storeobj==true
        obj(iter)=sum(d2_adj(keep_idx))/kfactor;
    end

    mu_diff=max(abs(mus(:)-mus_old(:)));
    sigma_diff=max(abs(sigs(:)-sigs_old(:)));
    if tol_sigma
        dif=max(mu_diff,sigma_diff);
    else
        dif=mu_diff;
    end
end

%% EM imputation of missing values
if condmeanimp==true
    Yimp=mdImputeCondMean(Y,mus,sigs);
else
    Yimp=[];
end

if stochimp==true
    stochYimp=mdImputeStochastic(Y,mus,sigs);
else
    stochYimp=[];
end

if storeobj==true
    obj=obj(1:iter);
end

out.loc=mus;
out.cov=sigs;
out.iter=iter;
out.weights=w;
out.Yimp=Yimp;
out.stochYimp=stochYimp;
out.obj=obj;
out.kfactor=kfactor;
out.kinfo=kinfo;
out.consistencyfactor=consistencyfactor;
out.cthr=cthr;
out.adjustedD2=d2_adj;

end

%% ------------------------------------------------------------------
function kinfo=localPatternFactors(nanY,w,poss,cthr,p,n,Sigma,method)
% Exact Gaussian Tallis factors, one per retained missingness pattern.
keep=w>0;
patt=unique(nanY(keep,:),'rows');
G=size(patt,1);

pobs=zeros(G,1);
nkept=zeros(G,1);
athr=zeros(G,1);
gammag=zeros(G,1);
kg=nan(G,1);

for g=1:G
    pg=patt(g,:);
    rows=keep & all(nanY==pg,2);
    obs=~pg;
    pobs(g)=sum(obs);
    nkept(g)=sum(rows);

    if pobs(g)==0
        athr(g)=0;
        gammag(g)=0;
        kg(g)=NaN;
        continue
    end

    athr(g)=localInvAdjust(cthr,pobs(g),p,n,Sigma,obs,method);

    if athr(g)<=0 || ~isfinite(athr(g))
        gammag(g)=0;
        kg(g)=NaN;
    else
        gammag(g)=chi2cdf(athr(g),pobs(g));
        if gammag(g)<=0
            kg(g)=NaN;
        else
            kg(g)=chi2cdf(athr(g),pobs(g)+2)/gammag(g);
        end
    end
end

kinfo=table(pobs,nkept,athr,gammag,kg);

% Existing Gaussian stabilization for small retained patterns.
kbar=localWeightedFactor(kinfo);
bad=~isfinite(kinfo.kg) | kinfo.kg<=0 | kinfo.kg>1 | ...
    kinfo.nkept<max(2,kinfo.pobs);
kinfo.kg(bad)=kbar;
end

%% ------------------------------------------------------------------
function kinfo=localAdaptiveFactors(nanY,w,cthr,p,n,Sigma,method, ...
    Yref,muref,sigmaref,poolByDimension,minRef)
% Distribution-adaptive radial factors estimated from complete references.
%
% The factor for pattern g is mean(Q_g | Q_g<=a_g)/p_g. Reference radii
% are always evaluated with the frozen complete-case reference fit.
keep=w>0;
patt=unique(nanY(keep,:),'rows');
G=size(patt,1);

pobs=zeros(G,1);
nkept=zeros(G,1);
athr=zeros(G,1);
gammag=nan(G,1);
kg=nan(G,1);
nref=zeros(G,1);
nrefkept=zeros(G,1);

for g=1:G
    pg=patt(g,:);
    rows=keep & all(nanY==pg,2);
    obs=~pg;
    pobs(g)=sum(obs);
    nkept(g)=sum(rows);

    if pobs(g)==0
        athr(g)=0;
        gammag(g)=0;
        continue
    end

    athr(g)=localInvAdjust(cthr,pobs(g),p,n,Sigma,obs,method);
    if athr(g)<=0 || ~isfinite(athr(g))
        gammag(g)=0;
        continue
    end

    if poolByDimension
        sameDim=find(sum(~patt,2)==pobs(g));
        qcell=cell(numel(sameDim),1);
        for jj=1:numel(sameDim)
            obsjj=~patt(sameDim(jj),:);
            qcell{jj}=localReferenceRadii(Yref,muref,sigmaref,obsjj);
        end
        q=vertcat(qcell{:});
    else
        q=localReferenceRadii(Yref,muref,sigmaref,obs);
    end

    nref(g)=numel(q);
    inside=q<=athr(g);
    nrefkept(g)=sum(inside);
    gammag(g)=mean(inside);
    if nrefkept(g)>=max(minRef,pobs(g)+1)
        kg(g)=mean(q(inside))/pobs(g);
    end
end

factorType=repmat("adaptive",G,1);
kinfo=table(pobs,nkept,athr,gammag,kg,nref,nrefkept,factorType);
end

%% ------------------------------------------------------------------
function q=localReferenceRadii(Yref,muref,sigmaref,obs)
% Compute squared projected Mahalanobis radii under the frozen reference.
X=Yref(:,obs)-muref(obs)';
S=sigmaref(obs,obs);
S=(S+S')/2;
[R,flag]=chol(S);
if flag~=0
    error('FSDA:mdTEM:NonPDAdaptiveProjection', ...
        'A projected adaptive reference scatter matrix is not positive definite.');
end
Z=X/R;
q=sum(Z.^2,2);
end

%% ------------------------------------------------------------------
function kbar=localWeightedFactor(kinfo)
% Information-weighted average of positive finite pattern-wise factors.
ok=isfinite(kinfo.kg) & kinfo.kg>0 & kinfo.nkept>0;
if ~any(ok)
    kbar=1;
    return
end
wgt=kinfo.nkept(ok).*kinfo.pobs(ok);
kbar=sum(wgt.*kinfo.kg(ok))/sum(wgt);
if ~isfinite(kbar) || kbar<=0
    kbar=1;
end
end

%% ------------------------------------------------------------------
function a=localInvAdjust(c,pg,p,n,Sigma,obs,method)
% Inverse of the adjusted-distance mapping for a fixed missingness pattern.
method=string(method);

switch method
    case "pri"
        a=c-(p-pg);

    case "expScale"
        a=c*pg/p;

    case "zMap"
        a=pg+sqrt(pg/p)*(c-p);

    case "detMap"
        gfull=exp(localLogdetSPD(Sigma)/p);
        gobs=exp(localLogdetSPD(Sigma(obs,obs))/pg);
        a=c*(pg/p)*(gobs/gfull);

    case "chiMap"
        u=chi2cdf(c,p);
        u=min(max(u,eps),1-eps);
        a=chi2inv(u,pg);

    case "betaMap"
        cn=(n-1)^2/n;
        if n<=p+1 || n<=pg+1
            a=c;
            return
        end
        u=min(max(c/cn,0),1-eps);
        al=betacdf(u,p/2,(n-p-1)/2);
        al=min(max(al,eps),1-eps);
        a=cn*betainv(al,pg/2,(n-pg-1)/2);

    case "impMD"
        a=c;

    otherwise
        a=c;
end
end

%% ------------------------------------------------------------------
function [musNew,sigsNew,kbar]=localCorrectedStep(Y,nanY,w,mus,sigs,kinfo)
% Corrected E- and M-step, grouped by retained missingness pattern.
%
% Only the data-driven conditional second-moment term is divided by k_g.
% The conditional covariance contribution is left uncorrected. Location is
% not corrected. kinfo and patt are constructed from the same retained rows
% and therefore have the same row ordering.
p=size(Y,2);
mus=mus(:);
keep=w>0;

S=zeros(p,p);
m1=zeros(p,1);
h=0;
patt=unique(nanY(keep,:),'rows');

if height(kinfo)~=size(patt,1)
    error('FSDA:mdTEM:InternalPatternMismatch', ...
        'Internal pattern-factor table does not match retained patterns.');
end

for g=1:size(patt,1)
    pg=patt(g,:);
    rows=keep & all(nanY==pg,2);
    hg=sum(rows);
    if hg==0
        continue
    end

    o=find(~pg);
    m=find(pg);
    if isempty(o)
        continue
    end

    kg=kinfo.kg(g);
    if ~isfinite(kg) || kg<=0
        error('FSDA:mdTEM:InvalidPatternFactor', ...
            'A retained pattern has an invalid consistency factor.');
    end

    Z=Y(rows,o)-mus(o)';
    Szz=(Z'*Z)/kg;
    sZ=sum(Z,1)';

    S(o,o)=S(o,o)+Szz;
    m1(o)=m1(o)+sZ;

    if ~isempty(m)
        Soo=sigs(o,o);
        Ag=sigs(m,o)/Soo;
        Cg=sigs(m,m)-Ag*sigs(o,m);
        Cg=(Cg+Cg')/2;

        SzzAg=Szz*Ag';
        S(o,m)=S(o,m)+SzzAg;
        S(m,o)=S(m,o)+SzzAg';
        S(m,m)=S(m,m)+Ag*SzzAg+hg*Cg;
        m1(m)=m1(m)+Ag*sZ;
    end

    h=h+hg;
end

if h==0
    musNew=mus;
    sigsNew=sigs;
    kbar=1;
    return
end

m1=m1/h;
S=S/h;
musNew=mus+m1;
sigsNew=S-(m1*m1');
sigsNew=(sigsNew+sigsNew')/2;
kbar=localWeightedFactor(kinfo);
end

%% ------------------------------------------------------------------
function val=localLogdetSPD(S)
S=(S+S')/2;
[R,flag]=chol(S);
if flag==0
    val=2*sum(log(diag(R)));
else
    val=log(max(det(S),realmin));
end
end

%FScategory:MULT-MissingData
