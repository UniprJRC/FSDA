function out = hpfilterFS(y, varargin)
%hpfilterFS HP filter with missing/excluded observations via selection matrix.
%
%<a href="matlab: docsearchFS('hpfilterFS')">Link to the help function</a>
%
% Model interpretation (Gaussian).
%
%   $y=(y_1, y_2, \ldots, y_T)'$ is the observed time series. $y_bsb$ is a
%   subset of units from $y$ of length $n_{bsb}$. 
%
%   $Sy=y_{bsb} = Sm + S\epsilon$,        $\epsilon \sim N(0, \sigma_\epsilon^2 I_{nbsb})$;
%   $D m = u$,                 $u   \sim N(0, \sigma^2_\epsilon (1/\lambda) I_T)$.
%
%   HP ratio: $\lambda=\sigma^2_\epsilon/\sigma^2_u$,
%             the greater $\lambda$, the smoother is the trend.
%
%   $D$ is the Second-difference matrix of size (T-2 x T): 
%       each row has [1 -2 1].
%   $S$ is the selection matrix of size nbsb x T. Note that $SS'=I_{nbsb}$.  
%
%   Conditioning on observed subset y_bsb = S y,
%   using W = S'S (diag 0/1).
%
% Posterior mean:
%   mhat = $\hat m= E(m|y_{bsb})=argmin_m ||S(y-m)||^2 + \lambda ||D m||^2$
%         = $(W + \lambda D'D) \ (W y)$.
%
%   Posterior cov:
%  $Cov(m|y_{bsb}) = \sigma^2_\epsilon  (W + \lambda D'D)^{-1}$
%
%  Required input arguments:
%
%    y:         Time series to analyze. Vector or timetable. A row or a column vector
%               with T elements, which contains the time series. Note that
%               y may contain missing values. If y is a timetable then
%               times of the timetable are shown in the plot of y and
%               fitted values.
%                 Data Types - double or timetable
%
%
%  Optional input arguments:
%
%
%   bsb :   Indices of observed points to condition on. Vector of
%           length>=3. Indices of observed points used in the conditioning.
%           If bsb is empty the the first 80 per cent units are used
%               Example - 'bsb',1:round(T/0.7): T=length(y);
%               Data Types - double
%
%      conflev : confidence level for the confidence bands. Scalar.
%                A number between 0 and 1 which defines the confidence
%                level which is used to produce the bands. The default
%                value of conflev is 0.99.
%               Example - 'conflev',0.999
%               Data Types - double
%
%   lambda : HP smoothing parameter. Numeric scalar or empty value.
%           The greater is the value of lambda the greater is the degree of
%           smoothness which is used to estimate the trend. The default
%           value of lambda which is used depends on the periodicity of the
%           series (see optional input parameter s).
%           If lambda is empty (default) the value of lambda
%           which is used depends on the Ravn–Uhlig scaling rule (adjust by
%           the 4th power of the observation frequency ratio):
%           \[
%               \lambda(s)=1600(s/4​)^4
%           \]
%           where $s$ is the number of observations per year
%           Therefore, lambda=1600 for quarterly data, lambda=129600
%           for monthly data...
%               Example - 'lambda',1000;
%               Data Types - double or empty.
%
%  predint : prediciton interval. String or empty value.
%            String which specifies for which units to compute the variance
%            of the prediction interval, that is the variance of the
%            estimated trend plus the variance of the estimated noise
%            component. If predint="nobsb" variance of the prediction is
%            computed just for the units not beloing to bsb. If predint is
%            "all" variance of the prediciton interval is computed for all
%            the units. If predint="" variance of prediction interval is
%            not computed.
%               Example - 'predint',"all";
%               Data Types - string or empty.
%
%     s    : length of seasonal period. Numeric scalar or empty value.
%           For monthly data s=12 (default), for quarterly data s=4,
%               Example - 's',52;
%               Data Types - double or empty.
%
% sigma2_eps: method to use to estimate residual variance. String of
%             scalar numeric value. If sigma2_eps is a string possible
%             values are: "MLaug" "MAPjef" "MAPig" "dfREML".
%             Given:
%             $RSS=||y_{bsb} - S \hat m||^2$.
%             and Qhat:
%             $\hat Q =RSS + \lambda ||D \hat m||^2$;
%             $K = n_{bsb} + (n-2)$;
%             $MLaug= ML (augmented) = \hat Q /K$.
%             MAPjef= MAP (maximum a posteriori estimate) with Jeffreys prior.
%               $MAPjef=\hat Q / (K+2)$.
%             MAPig = MAP (maximum a posteriori
%               estimate) with inverse gamma prior with parameters a0 and b0.
%               $MAPig= (b_0 + 0.5\hat Q)/(a_0 + 1 + 0.5K)$.
%            dfREML = df-REML-like (marginal smoother likelihood).
%               $dfREML=RSS / (nbsb - df(lambda))$.
%               $df(lambda)=trace(H)
%               = trace(S A^{-1} S') = trace(A^{-1} S' S) = trace(A^{-1} W)$.
%               The estimate of the $trace(A^{-1} W)$ is via Hutchinson:
%               $tr(M) ≈ (1/niter) \sum_{i=1}^{niter} z'_i M z_i$.
%               $z_i$ are Rademacher random. $niter$ is fixed to 50.
%               variables, equal to +1 or −1 with equal probability.
%            On the other hand, if sigma2_eps is a numeric scalar it is
%            possible to supply the prior value.
%               Example - 'sigma2_eps',20;
%               Data Types - string or scalar double.
%
%      ig_a0 :  prior value of a in inverse gamma prior. Numeric scalar.
%               Prior value of parameter a of the inverse gamma to be used
%               in case sigma2_eps="MAPig".
%               Example - 'ig_a0',10;
%               Data Types - scalar double.
%
%
%      ig_b0 :  prior value of b in IG. Numeric scalar.
%               Prior value of parameter b of the inverse gamma to be used
%               in case sigma2_eps="MAPig".
%               Example - 'ig_a0',10;
%               Data Types - scalar double.
%
%      plots  : plot on the screeen of HP trend estimate. Boolean.
%               If plots = true a plot with the real time series with fitted
%               values and trend estimate will appear on
%               the screen. This plot is tagged forecastTS.
%               The confidence bands which are shown depend on the input option
%               predint.
%               The default value of plot is false, that is no plot is shown on
%               the screen.
%               Example - 'plots',true;
%               Data Types - logical.
%
%
% Output:
%
%
%  out :     A structure containing the following fields
%
%   out.mhat        = n x 1 posterior mean of y. Estimated HP trend
%   out.bsb         = indexes of the units used in the fit.
%   out.sigma2_ML   = estimate of residual variance using augmented
%                     likelihood approach.
%   out.sigma2_Jef  = MAP estimate of residual variance using Jeffreys
%                     prior.
%   out.sigma2_IG   = MAP estimate of residual variance based on
%                     inverse gamma prior.
%   out.sigma2_dfREML= MAP estimate of residual variance based on
%                     marginal smoother likelihood.
%   out.predVar    = n x 1 predictive variance or missing. This is the
%                     estimated trend variance + observation noise
%                     variance. The observation noise variance depends on
%                     the value of input option sigma2_eps. out.predVar is
%                     a scalar missing if input option predint is empty.
%                     out.PredVar is populated just for the units not
%                     belonging to bsb (input option predint is "nobsb"),
%                     or for all the units (input option predint is "all")
%   out.PI_low    = n x 1 predictive variance or missing. This is the
%                    lower band of the confidence interval.  out.PI_low is
%                     a scalar missing if input option predint is empty.
%                     out.PI_low is populated just for the units not
%                     belonging to bsb (input option predint is "nobsb"),
%                     or for all the units (input option predint is "all")
%   out.PI_high    = n x 1 predictive variance or missing. This is the
%                    upper band of the confidence interval.  out.PI_high is
%                     a scalar missing if input option predint is empty.
%                     out.PI_high is populated just for the units not
%                     belonging to bsb (input option predint is "nobsb"),
%                     or for all the units (input option predint is "all")
%
%
% See also hpfilter, LTSts, supsmu
%
% References:
%
% Hutchinson, M. F. (1989), A stochastic estimator of the trace of the influence
% matrix for laplacian smoothing splines, "Communications in Statistics
% Simulation and Computation", Vol. 18, pp. 1059–1076.
%
% Hodrick, R. J., & Prescott, E. C. (1997), Postwar U.S. Business Cycles:
% An Empirical Investigation, "Journal of Money, Credit and Banking",
% Vol. 29, pp. 1–16.
%
% Ravn, M. O., & Uhlig, H. (2002), On Adjusting the Hodrick–Prescott Filter
% for the Frequency of Observations, "Review of Economics and Statistics",
% Vol. 84, pp. 371–376.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('hpfilterFS')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Call to hpfilterFS with all default arguments.
    Mdl = arima('Constant',0,'D',1,'MA',{0.5},'Variance',100);
    n = 60;
    y = simulate(Mdl,n);      
    % In this case just the conditional mena 
    out = hpfilterFS(y); 
%}

%{
    % Call to hpfilterFS with optional argument bsb.
    rng(1000)
    Mdl = arima('Constant',0,'D',1,'MA',{0.5},'Variance',100);
    n = 150;
    y = simulate(Mdl,n); 
    bsb = (1:round(n*0.9))';      
    
    % In this case just the conditional mena 
    out = hpfilterFS(y,'bsb',bsb); 
%}

%{
    %% Call to hpfilterFS with optional arguments bsb and predint.
    rng(1000)
    Mdl = arima('Constant',0,'D',1,'MA',{0.5},'Variance',100);
    n = 150;
    y = simulate(Mdl,n); 
    bsb = (1:round(n*0.9))';      
    % Prediction interval on all the observations (included and excluded)
    out = hpfilterFS(y,'bsb',bsb,'predint','all','plots',true); 
%}

%{
    % Check equality with output of function hpfilter.
    load Data_SchwertStock
    TTM = rmmissing(DataTimeTableMth);
    % Aggregate the monthly data in the timetable to quarterly measurements.
    TTQ = convert2quarterly(TTM);
    % Apply the Hodrick-Prescott filter to all variables in the quarterly timetable. The default smoothing parameter value is 1600. Display the last few observed components.
    TQTT = hpfilter(TTQ);
    TQTTchk=TQTT;
    for j=1:size(TQTT,2)
        outj=hpfilterFS(TTQ{:,j},'plots',0);
        TQTTchk{:,j}=outj.mhat;
    end
    maxdiff=max(abs(TQTT{:,:}-TQTTchk{:,:}),[],"all");
    assert(maxdiff<1e-10,"output of hpfilterFS different from hpfilter")
%}

%% Beginning of code

if nargin<1
    error('FSDA:hpfilterFS:missingInputs','Required input argument is missing.')
end

conflev = 0.99;
s=4;
bsb=[];
plots=false;
lambda=1600;
ig_a0=[];
ig_b0=[];

% Method to use to estimate the residual variance
sigma2_eps="MLaug";
predint=[];

if nargin > 1

    options=struct('bsb',bsb,'plots',plots,'lambda',lambda,'conflev',conflev, ...
        'sigma2_eps',sigma2_eps,'s',s,'ig_a0',ig_a0,'ig_b0',ig_b0,'predint',predint);

    for i=1:2:(length(varargin)-1)
        options.(varargin{i})=varargin{i+1};
    end


    s=options.s;
    plots=options.plots;
    bsb=options.bsb;
    lambda=options.lambda;
    conflev =options.conflev;
    sigma2_eps=options.sigma2_eps;
    ig_a0=options.ig_a0;
    ig_b0=options.ig_b0;
    predint=options.predint;
end

if istimetable(y)
    isTT=true;
    rowTimes=y.Properties.RowTimes;
    % Given rowTimes find time series periodicity
    stent=findTimeSeriesPeriodicity(rowTimes);
    % Overwrite the value of s with found periodicity
    if ~isempty(stent)
        s=stent;
    end
    y=y{:,1};
else
    y = y(:);
    isTT=false;
end

% default value of lambda depending on the sampling frequency. We use
% Ravn–Uhlig scaling rule (adjust by the 4th power of the observation
% frequency ratio)
if isempty(lambda)
    if s==4
        lambda=1600;
    else
        lambda = 1600 * (s/4)^4;
    end
end



T = length(y);
seq=(1:T)';
if T < 3
    error('FSDA:hpfilterFS:WrongInputOpt','Need n >= 3 for HP second differences.');
end

if isempty(bsb)
    bsb=seq;
else
    % Validate bsb
    bsb = bsb(:);
    bsb = unique(bsb);
    if any(bsb < 1) || any(bsb > T)
        error('FSDA:hpfilterFS:WrongInputOpt','bsb contains indices outside 1..n.');
    end
end
nbsb=length(bsb);

% Selection via diagonal weights: W = S'S (n x n), with 1 on observed
% entries
w = zeros(T,1);
w(bsb) = 1;
W = spdiags(w, 0, T, T);   % sparse diagonal

% Second-difference matrix D (T-2 x T): each row has [1 -2 1]
e = ones(T,1);
D = spdiags([e -2*e e], 0:2, T-2, T);

% % Build P = D'*D as 5-diagonal sparse matrix (T x T)
% d0 = [1; 5; repmat(6,T-4,1); 5; 1];
% d1b = [-2; repmat(-4,T-3,1); -2; 0];
% d1a = [-2; -2; repmat(-4,T-3,1);-2];
% d2 = e;
% DtD = spdiags([d2 d1b d0 d1a d2], [-2 -1 0 1 2], T, T);

% Solve for posterior mean of trend:
A = W + lambda*(D'*D);
b = W*y;
Aready = decomposition(A, 'chol');

% Predictive mean for y is E[y|data] = E[m|data] since E[eps]=0
m_hat = Aready \ b;

% Compute RSS ||y_bsb - S m_hat||^2
r_obs = y(bsb) - m_hat(bsb);
RSS = full(r_obs' * r_obs);
% Compute Qhat =  RSS + lambda ||D m_hat||^2
Qhat = RSS + lambda * full((D*m_hat)'*(D*m_hat));

% Estimate of sigma^2 (fixed lambda)
%   1) ML (augmented)         : Qhat / K, K = nObs + (n-2)
K=nbsb + (T-2);
sigma2_ML   = Qhat / K;
%   2) MAP (Jeffreys)         : Qhat / (K+2)
sigma2_Jef  = Qhat / (K + 2);
%   3) MAP (Inv-Gamma a0,b0)  : (b0 + 0.5*Qhat)/(a0 + 1 + 0.5*K)
if isempty(ig_a0)
    ig_a0=1e-03;
end
if isempty(ig_b0)
    ig_b0=1e-03*var(y(bsb));
end
sigma2_IG   = (ig_b0 + 0.5*Qhat) / (ig_a0 + 1 + 0.5*K);
%   4) df-REML-like (smoother): RSS / (nObs - df(lambda)),
% df(lambda)=trace(H)
% trace(H) = trace(S A^{-1} S') = trace(A^{-1} S' S) = trace(A^{-1} W)
% Estimate trace(A^{-1} W) via Hutchinson: tr(M) ≈ (1/R) Σ z' M z, z_i=±1
% This methods correspond to a marginal smoother likelihood
nTrace=50;
df = hutch_trace_AinvW(Aready, W, nTrace);

if (nbsb - df) <= 0
    warning('nbsb - df(lambda) <= 0. df-based estimator not defined; returning NaN.');
    sigma2_dfREML = NaN;
else
    sigma2_dfREML = RSS / (nbsb - df);
end

if sigma2_eps=="MLaug"
    sigma2_eps = sigma2_ML;
elseif sigma2_eps=="MAPjef"
    sigma2_eps = sigma2_Jef;
elseif sigma2_eps=="MAPig"
    sigma2_eps = sigma2_IG;  % Use Inverse-Gamma estimate for sigma^2
elseif sigma2_eps=="dfREML"
    sigma2_eps=sigma2_dfREML;
else
    % Check that the value of sigma2 if it is a string is one
    % of the above values
    if isstring(sigma2_eps) & ~ismember(sigma2_eps,["MLaug" "MAPjef" "MAPig" "dfREML"])
        error('FSDA:hpfilterFS:InvalidInput', 'sigma2 must be a numeric scalar or a valid string option.');
    end
    if isnumeric(sigma2_eps) && isscalar(sigma2_eps)
        % use prior value
    else
        error('FSDA:hpfilterFS:InvalidInput', 'sigma2 must be a numeric scalar or a valid string option.');
    end
end

% Note that using var(y - y_hat) is not coherent with the model, because
% those residuals mix observation noise with trend-estimation error in a
% way that depends on the smoother.

% --- Prediction intervals need Var(y | data)
% Posterior covariance of m is sigma_eps^2 * A^{-1}
% We only need diagonal elements of A^{-1} for missing indices.
%
% Efficient trick: compute columns of A^{-1} corresponding to missing indices:
% Solve A * X = E, where E has columns e_{idxMiss(j)}.
% Then var_m(idxMiss(j)) = X(idxMiss(j), j)

if ~isempty(predint)
    switch predint
        case "all"
            idxMiss=seq;
        case "nobsb"
            % Forecasts for excluded indices:
            idxMiss=setdiff(seq,bsb) ;
        otherwise
            % predint must be equal to all or to nobsb
            % otherwise produce an error
            error('FSDA:hpfilterFS:InvalidInput', 'predint must be "all" or "nobsb" or []');
    end


    blockSize=1000;
    var_m = diag_inv_sparse_block(Aready, blockSize,idxMiss);

    % kMiss=length(idxMiss);
    %
    % E = sparse(idxMiss, 1:kMiss, 1, n, kMiss);   % n x kMiss selection columns
    % % One factorization, many RHS:
    % Aready = decomposition(A, 'chol');                % uses sparse Cholesky if available
    % X = Aready \ E;                                   % X = A^{-1} E
    %
    % % Extract posterior variance of m at missing indices
    % var_m = full( diag( X(idxMiss, :) ) );        % (A^{-1})_{ii} for i in idxMiss

    % niter=10000;
    % [dhat, info] = diag_inv_hutchinson(A, niter, 1e-8, 300, struct('type','ict','droptol',1e-3));
    % var_mCHK =  dhat;    % posterior Var(m_t|data) = sigma^2 * (A^{-1})_tt


    % Predictive variance adds observation noise variance
    predVar=NaN(T,1);
    predVar(idxMiss) = var_m+ sigma2_eps;

    % z-quantile
    z = norminv((1+conflev)/2);

    PI_low=NaN(T,1);
    PI_high=PI_low;
    PI_low(idxMiss)  = full(m_hat(idxMiss)) - z*sqrt(predVar(idxMiss));
    PI_high(idxMiss) = full(m_hat(idxMiss)) + z*sqrt(predVar(idxMiss));
else
    predVar = [];
    PI_low = [];
    PI_high = [];
end

if plots==true
    figure;
    if isTT ==true
        seq=rowTimes;
    end
    plot(seq, y, 'k-');
    hold on;
    if T<100
        plot(seq(bsb), y(bsb), 'o');
    end
    plot(seq, m_hat, 'b-', 'LineWidth', 1.5);
    if ~isempty(predint)
        plot(seq, PI_low, 'r--');
        plot(seq, PI_high, 'r--');
        if T<100
            legend('y','values of bsb','HP mean (all t)','PI low','PI high');
        else
            legend('y','HP mean (all t)','PI low','PI high');
        end
        title('y, HP-based predictive mean and prediction intervals');
    else
        if T<100
            legend('y','values of bsb','HP mean (all t)');
        else
            legend('y','HP mean (all t)');
        end
        title('HP-based predictive mean');
    end
    legend('AutoUpdate','off','Location','best')
    grid on;
    if max(bsb)<T
        xline(seq(max(bsb))+0.5)
    end

end


out = struct();
out.mhat       = full(m_hat);
out.bsb         = bsb;

% Store estimates of residual variance
out.sigma2_ML  =sigma2_ML;
out.sigma2_Jef =sigma2_Jef;
out.sigma2_IG  =sigma2_IG;
out.sigma2_dfREML=sigma2_dfREML;

% Store prediciton intervals
out.predVar    = predVar; % predictive variance
out.PI_low     = PI_low;
out.PI_high    = PI_high;


end

function est = hutch_trace_AinvW(A, W, niter)
% Estimate trace(A^{-1} W) using Hutchinson with niter Rademacher probes.
n = size(W,1);
acc = 0;
for i = 1:niter
    z = sign(randn(n,1));       % Rademacher +/-1
    v = W*z;
    x = A \ v;                  % x = A^{-1} W z
    acc = acc + (z' * x);
end
est = acc / niter;
end


% function [dhat, info] = diag_inv_hutchinson(A, niter, tol, maxit, ichol_opts)
% %DIAG_INV_HUTCHINSON  Approximate diag(inv(A)) for large sparse SPD A.
% %
% % dhat = (1/niter) * sum_r (x_r .* z_r), where A x_r = z_r and z_r are Rademacher.
% %
% % Inputs:
% %   A          : sparse SPD matrix (n x n)
% %   niter      : number of probes (e.g. 50-300)
% %   tol        : PCG tolerance (default 1e-8)
% %   maxit      : PCG max iterations (default 200)
% %   ichol_opts : options struct for ichol (default: type='ict', droptol=1e-3)
% %
% % Outputs:
% %   dhat : approx diagonal of inv(A)
% %   info : diagnostics (mean iters, flags, etc.)
%
% if nargin < 3 || isempty(tol), tol = 1e-8; end
% if nargin < 4 || isempty(maxit), maxit = 200; end
% if nargin < 5 || isempty(ichol_opts)
%     ichol_opts = struct('type','ict','droptol',1e-3);
% end
%
% n = size(A,1);
% dhat = zeros(n,1);
%
% % Preconditioner (very important for speed)
% L = ichol(A, ichol_opts);
%
% flags = zeros(niter,1);
% relres = zeros(niter,1);
% iters = zeros(niter,1);
%
% for r = 1:niter
%     z = sign(randn(n,1));  % Rademacher (+/-1)
%     [x,flag,rr,it] = pcg(A, z, tol, maxit, L, L');
%     flags(r) = flag;
%     relres(r) = rr;
%     iters(r) = it;
%
%     dhat = dhat + (x .* z);
% end
%
% dhat = dhat / niter;
%
% info = struct();
% info.flags = flags;
% info.relres = relres;
% info.iters = iters;
% info.mean_iters = mean(iters);
% info.fail_rate = mean(flags~=0);
% end

function d = diag_inv_sparse_block(A, blockSize, idx)
%DIAG_INV_SPARSE_BLOCK  Exact diagonal of inv(A) via block solves.
%
% d = diag_inv_sparse_block(A)
% d = diag_inv_sparse_block(A, blockSize)
% d = diag_inv_sparse_block(A, blockSize, idx)
%
% Inputs:
%   A         : sparse SPD matrix
%   blockSize : number of RHS columns per block (default 200)
%   idx       : optional subset of indices for which diagonal is needed
%
% Output:
%   d : diagonal entries (length n if idx not given, else length(idx))
%
% Method:
%   For each block of indices I:
%       Solve A X = E_I
%       Extract diag entries from X(I,:)

if nargin < 2 || isempty(blockSize)
    blockSize = 200;
end

n = A.MatrixSize(1);

if nargin < 3 || isempty(idx)
    idx = (1:n)';
else
    idx = idx(:);
end

k = length(idx);
d = zeros(k,1);

% Process in blocks
for startIdx = 1:blockSize:k
    stopIdx = min(startIdx + blockSize - 1, k);

    I = idx(startIdx:stopIdx);
    nb = length(I);

    % Build sparse identity block E_I
    E = sparse(I, 1:nb, 1, n, nb);

    % Solve A X = E_I
    X = A \ E;

    % Extract diagonal elements (A^{-1})_{ii}
    % For each column j, diagonal entry is X(I(j), j)
    d(startIdx:stopIdx) = full(diag(X(I,:)));
end
end

function s=findTimeSeriesPeriodicity(rowTimes)

% Compute day differences (numeric)
dayDiffs = days(diff(rowTimes));
% Use the modal (most common) rounded day step to be robust
dayStep = mode(round(dayDiffs));

% If data are sampled daily/weekly, dayStep will be small (<= 28)
if dayStep <= 28
    % Classify daily vs weekly (allow small jitter)
    if abs(dayStep - 1) <= 1        % 0..2 -> daily
        s = 365;
    elseif abs(dayStep - 7) <= 2    % 5..9 -> weekly
        s = 52;
    else
        s='';
    end
else
    % Compute month differences between consecutive times:
    yr = year(rowTimes);
    mo = month(rowTimes);
    monthsDiff = (yr(2:end) - yr(1:end-1)) * 12 + (mo(2:end) - mo(1:end-1));

    % Use the modal (most common) month step to be robust to occasional gaps
    monthStep = mode(monthsDiff);

    % Map month step to s:
    % 12 months -> yearly       -> s = 1
    % 6  months -> semiannual   -> s = 2
    % 4  months -> triannual    -> s = 3
    % 3  months -> quarterly    -> s = 4
    % 1  month  -> monthly      -> s = 12
    switch monthStep
        case 12
            s = 1;
        case 6
            s = 2;
        case 4
            s = 3;
        case 3
            s = 4;
        case 1
            s = 12;
        otherwise
            s = ''; % unknown / not one of the supported frequencies
    end
end
end


%FScategory:REG-Regression

