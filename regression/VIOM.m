function [out] = VIOM(y,X,dw,varargin)
%VIOM computes weights estimates under a Variance-Inflation Outlier Model
%using MLE or Restricted MLE (REMLE)
%
%<a href="matlab: docsearchFS('VIOM')">Link to the help function</a>
%
%
%  Required input arguments:
%
%    y:         Response variable. Vector. A vector with n elements that
%               contains the response
%               variable.  It can be either a row or a column vector.
%    X :        Predictor variables. Matrix. Data matrix of explanatory
%               variables (also called 'regressors')
%               of dimension (n x p-1). Rows of X represent observations, and
%               columns represent variables.
%
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%
%    dw:        Units flagged as possible VIOM outliers to be downweighted.
%
% Optional input arguments:
%
%   intercept:  Indicator for constant term. Scalar.
%               If 1, a model with constant term will be fitted (default),
%               else no constant term will be included.
%               Example - 'intercept',1
%               Data Types - double
%   mult:       Indicator for joint weights estimate. Scalar.
%               If mult==1 the weights are jointly estimated by iterative REML. 
%               Default is mult==0 and singularly optimal weights are 
%               estimated using REML closed form solution.
%               Example - 'mult',1
%               Data Types - double
%   trim:       Units flagged as possible MSOM outliers, they are forced to 
%               have 0 weights. By default no units are trimmed, i.e. trim==[]. 
%               Vector of labels. 
%               Example - 'trim',[1,2,3] 
%               Data Types - double
%   trsh:       Threshold on residuals. Scalar. 
%               If thrsh>0 all the (standard) residuals greater than trsh 
%               are set to 0. [[TBA:modify to studentized or scaled residuals]]
%               If trsh<1 all the estimated weights smaller than trsh are
%               forced to be 0. 
%               If trsh==0 (default option) no weights are forced to be 0.
%               (Note: It might be useful to reduce the computational burden).
%               Example - 'trsh',5 
%               Data Types - double
%   cook:       Use Cook et al. (1982) formula to estimate single weights 
%               using MLE. Scalar. Default cook==0 and Thompson (1985) formula 
%               based on REMLE is used.
%               Example - 'cook',1 
%               Data Types - double
%
%  Output:
%
%  out :     A structure containing the following fields
% 
%           out.w    =  n x 1 vector of weights.
%           out.beta =  p x 1 vector of estimated coefficients based on WLS.
%
%
% See also: FSR, FSRcore
%
% References:
%
% Cook, R. Dennis, Norton Holschuh, and Sanford Weisberg. "A note on an 
%       alternative outlier model." Journal of the Royal Statistical Society: 
%       Series B (Methodological) 44.3 (1982): 370-376.
% Thompson, Robin. "A note on restricted maximum likelihood estimation with 
%       an alternative outlier model." Journal of the Royal Statistical 
%       Society: Series B (Methodological) 47.1 (1985): 53-55.
% Gumedze, Freedom N. "Use of likelihood ratio tests to detect outliers 
%       under the variance shift outlier model." Journal of Applied 
%       Statistics 46.4 (2019): 598-620.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('VIOM')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % VIOM with default input.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)*2;
    [out]=VIOM(y,X, 1:5);
    out.w(1:5);
    out.beta;
%}

%{
    % VIOM with optional arguments.
    % Use MLE for single weights.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)*2;
    [out]=VIOM(y,X,1:5,'cook',1);
    out.w(1:5);
    out.beta;
%}

%{
    %% VIOM with optional arguments.
    % Use MLE for single weights and with pre-specified trimming.
    n=200;
    p=3;
    randn('state', 123456);
    X=randn(n,p);
    y=randn(n,1);
    y(1:5)=y(1:5)*2;
    y(6:15)=y(6:15)+10;
    [out]=VIOM(y,X,1:5,'mult',1,'trim',6:15);
    out.w(1:15);
    out.beta;
%}


%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% User options

options=struct('intercept', 1, 'mult', 0, ...
    'trsh', 0, 'trim', [], 'cook', 0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('wrong input for VIOM.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 3
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

% intercept = options.intercept;
mult = options.mult;
trsh = options.trsh;
trim = options.trim;
cook = options.cook;

% parameters check [TBA more]
[m,q]=size(dw);
if min(m,q)>1
    error('FSDA:chkinputR:Wrongy','dw is not one-dimensional.');
elseif q~=1
    % If dw is a row vector it is transformed in a column vector
    dw=dw';
end
[m,q]=size(trim);
if min(m,q)>1
    error('FSDA:chkinputR:Wrongy','trim is not one-dimensional.');
elseif m>0 || q>0
    if m==1
        % If trim is a row vector it is transformed in a column vector
        trim=trim';
    end
end

% parameters check [TBA more checks]

%%

% clean units
cll = setdiff(1:n, [dw; trim]);
% weights vector
w = ones(n, 1);
if ~isnan(trim)
    w(trim) = 0;
elseif isnan(trim) 
    trim = [];
end

% [TBA: optional input beta vector]
% [TBA: possible double trimming: on res and w]
% remove residuals>trsh in order treat them as MSOM since the begin
% (it allows us to save some computing time)    
if trsh >= 1
    e = y - X * regress(y(cll), X(cll,:));
    ind_excl = dw(abs(e(dw)) > trsh);
    w(ind_excl) = 0;
else
    ind_excl = [];
end
% remove excluded units
dww = setdiff(dw, ind_excl);   
n_dww = length(dww);
        
% some outliers have to be downweighted        
if any(~isnan(dww) & ~isnan(dw))

    % joint estimate all weights
    if mult == 1

            % downweight ALL outliers 
            % NOTE: we use REMLE in a mixed model as Thompson (1985)
            dum = dummyvar(1:n);
            dum = dum(:, dww);
            dum = mat2cell(dum, n, ones(1,n_dww));
            % [TBA: we might use fminunc as optimizer]
            lme = fitlmematrix(X, y, dum, [], 'CovariancePattern', ...
                cellstr(string(repmat('Diagonal', n_dww, 1))), ...
                'FitMethod', 'REML', 'Exclude', [ind_excl, trim]); 
                    % 'Optimizer', 'fminunc'); [TBA]
                    % [~,~,STATS] = randomEffects(lme); 
            [V,S,~] = covarianceParameters(lme);
            w(dww) = 1 ./ (1 + cell2mat(V) ./ S);
            if trsh <1 && trsh>0
                % set to 0 weights smaller than a trsh 
                w(dw(w(dw)<trsh)) = 0; 
                W = w .* eye(n);
                beta = (X'*W*X)\X'*W*y;
            else
                beta = fixedEffects(lme);
                % check that it is equal to:
%                 W = w .* eye(n);
%                 beta = (X'*W*X)\X'*W*y;
            end
    
    % singularly estimate each weight
    % i.e. assuming that one outlier per time is included in the robust fit
    else

        % downweight EACH outlier independently
        % NOTE: we use the closed form solution of Thompson, 1985
        % we simulate the inclusion of each outlier independently
        nn = n - length(trim) - length(dw) + 1;

        % initialize result 
        % leverage values
        h = nan(length(dw), 1);
        % std. deviation
        sigma = nan(length(dw), 1);
        % 
        r_incl_i = nan(length(dw), 1);

        % needed quantities to simulate the inclusion of each outlier
        % ([TBA] we might simply use MSOM updating formulas)
        % ([TBA] we might extract these quantities from FSR)
        for out_i = 1:length(dw)
            % i-th outlier
            out_i_id = dw(out_i);
            % its inclusion in the fit
            ind_incl_i = sort([cll, out_i_id]);
            ind_incl_i = ind_incl_i(~isnan(ind_incl_i));
            % i-th beta estimates
            beta_i = regress(y(ind_incl_i), X(ind_incl_i,:));
            % all residuals
            r_tot = y - X * beta_i;
            % i-th res
            r_incl_i(out_i) = r_tot(out_i_id); 
            % i-th leverage value
            h(out_i) = X(out_i_id,:) * ...
                ((X(ind_incl_i,:)' * X(ind_incl_i,:)) ...
                \ X(out_i_id,:)');          
            % i-th std dev est
            sigma(out_i) = sqrt((r_tot(ind_incl_i)'*r_tot(ind_incl_i)) / (nn-p));
        end

        % stud residual
        t = r_incl_i ./ (sigma .* sqrt(1 - h));

        % enforce possible solution
        if sum(t.^2 <= 1) > 0
            warning('In estimating single weights some units got a weight=1');
            warning('t^2<1 in VIOM for unit: %s \n, ', string(dw(t.^2 <= 100)));
            t = t(t.^2 > 1);
            h = h(t.^2 > 1);
        end
        
        % Cook et al. (1982) formula
        if cook==1
            t = sqrt(((n - p) / n)) .* r_incl_i ./ (sigma .* sqrt(1 - h));
            b = t.^2 .* (n + 2 .* h - 1) -2 * h * (n - p);
            del = sqrt(t.^2 .* (4 * n * h .* (t.^2 - n + p) + (n - 1)^2 .* (t.^2)));
            a = (1 - h) .* (n - p - t.^2);
            out1 = (b + del) ./ (2 * a);
            out2 = (b - del) ./ (2 * a);
            w(dw) = max(out1, out2);
        
        % using REML
        else
            % thompson (1985) formula
            % w(dw(t.^2 > 1)) = ((nn - p) .* (t.^2 - h) - t.^2 .* (1 - h)) ./ ...
            %     ((nn - p) .* (1 - h) - t.^2 .* (1 - h));
            % check that it is equivalent to Gumzede (2019)
            w(dw(t.^2 > 1)) = 1 + ((nn - p) .* (t.^2 - 1)) ./ ...
              ((nn - p - t.^2) .* (1 - h));
        
        end
      
        % enforce real part        
        if ~isreal(w)
            warning('w is no real in VIOM')
        end
        % est weights
        w(dw(t.^2 > 1)) = real(1./w(dw(t.^2 > 1))); 
        if trsh <1 && trsh >0
            % set to 0 weights smaller than a trsh 
            w(dw(w(dw)<trsh)) = 0; 
        end
        W = w .* eye(n);
        % [TBA: WLS update using general Sherman formula]
        beta = (X'*W*X)\X'*W*y;          
                    
    end

else
    W = w .* eye(n);
    beta = (X'*W*X)\X'*W*y;
end

out.w = w;
out.beta = beta;

end
