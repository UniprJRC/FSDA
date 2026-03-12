function d2_adj = NApartialMD2full(d2, p, pobs, varargin)
% NApartialMD2full  rescale partial squared distances to full-dimension.
%
%<a href="matlab: docsearchFS('NApartialMD2full')">Link to the help function</a>
%
% Required input arguments:
%
%   d2    : squared partial distances. Vector. 
%           Vector of length n, squared distances computed on observed dims
%           (from MD_partial)
%   p     : full dimension. Scalar. p is the number of variables in the
%           original data matrix.
%   pobs  : number of observed variables. Vector.
%           vector of length n containing the number of observed variables
%           (per row)
%
%  Optional input arguments:
%
%  typeresc : method to rescale. Positive integer.
%             Method which must be used to rescale the variables.
%             typeresc=1 use of asymptotic Chi2 distribution.
%             typeresc=2 use of exact Beta distribution.
%             typeresc=3 use of expectation correction
%             typeresc=4 use of standardization correction.
%                 Example - 'typeresc',1
%                 Data Types - positive integer
%
%  Output:
%
%   d2_adj : column vector of length n, where each element is the adjusted
%            squared distance:
%            For example if typeresc=1,
%            d2_adj(i) = chi2inv( chi2cdf(d2(i), poss(i)), p );
%            for i=1, 2, ..., n.
%           Notet that if poss(i) is 0 or d2(i) is NaN, result is NaN.
%
%
% Copyright 2008-2026.
% Written by FSDA team
%
%
% See also: NApartialMD.m, NAem.m, NAimputeConditionalmean
%
%
% References:
%
% Little, R. J. A., & Rubin, D. B. (2019). Statistical Analysis with
% Missing Data (3rd ed.). Hoboken, NJ: John Wiley & Sons.
% van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.).
% Boca Raton, FL: Chapman & Hall/CRC (Taylor & Francis Group).
% Templ, M. (2025). Visualization and Imputation of Missing Values: With
% Applications in R. Cham, Switzerland: Springer Nature.
%
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('NApartialMD2full')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

%{
    % Using beta rescaling.
    p = 5;                % number of variables
    n = 100;             % number of observations
    rho = 0.9;            % target pairwise correlation (0<rho<1)
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);      % upper-triangular such that Sigma = R'*R
    missRate = 0.01;     % MCAR missing probability per entry
    
    % Generate samples ~ N(0,Sigma)
    Yfull = randn(n,p) * R;   % Strong positive correlation between the vars
    missMask = rand(n,p) < missRate;
    mu=mean(Yfull)';
    S=cov(Yfull);
    Y=Yfull;
    Y(missMask) = NaN;
    [d2partial,pobs]=NApartialMD(Y,mu,Sigma);
    d2_adj=NApartialMD2full(d2partial,p,pobs);
%}

%% Beginning of code

d2 = d2(:);
pobs = pobs(:);
typeAdj=2;

if nargin>3
    options=struct('typeAdj',typeAdj);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:NApartialMD2full:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        aux.chkoptions(options,UserOptions)
    end

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    typeAdj=options.typeAdj;
end

if numel(d2) ~= numel(pobs)
    error('FSDA:partialMW2full:WrongInputOpt','d2 and pobs must have same length');
end

n = numel(d2);
d2_adj = nan(n,1);

if typeAdj==1 % rescale using Chi2

    for i = 1:n
        p1 = pobs(i);
        if isnan(d2(i)) || p1 <= 0
            d2_adj(i) = NaN;
            continue
        end

        c = chi2cdf(d2(i), p1);
        c = min(max(c, eps), 1 - eps); % avoid 0 or 1
        d2_adj(i) = chi2inv(c, p);

    end
elseif typeAdj==2 % rescale using Beta
    c1 = (n-1)^2 / n;
    epsv = eps;

    % parameters full
    a2 = p/2;
    b2 = (n-p-1)/2;


    for i=1:n

        p1 = pobs(i);

        if isnan(d2(i)) || p1<=0 || n<=p1+1
            d2_adj(i)=NaN;
            continue
        end

        % parameters subset
        a1 = p1/2;
        b1 = (n-p1-1)/2;

        % scaled variable
        u = d2(i)/c1;
        u = min(max(u,0),1-epsv);

        % probability under subset exact distribution
        alpha = betacdf(u,a1,b1);
        alpha = min(max(alpha,epsv),1-epsv);

        % invert under full exact distribution
        v = betainv(alpha,a2,b2);

        % rescale back
        d2_adj(i) = c1*v;
    end

elseif typeAdj==3 % rescale using Expectation
    idx = pobs>0 & ~isnan(d2);
    d2_adj(idx) = d2(idx) .* (p ./ pobs(idx));

elseif typeAdj ==4 % rescale using Standardization
    idx = pobs>0 & ~isnan(d2);
    z = (d2(idx) - pobs(idx)) ./ sqrt(2*pobs(idx));
    d2_adj(idx) = p + sqrt(2*p).*z;

else
    error('FSDA:partialMW2full:WrongInputOpt','typeresc must 1, 2, 3 or 4')
end

end

% Below is an attempt to write a precision chi2cdf using vpa
% digits(70);
% chi2cdf_hp = @(xval,kval) vpa(1 - igamma( sym(kval)/2, sym(xval)/2 ) / gamma( sym(kval)/2 ));
% % chi2inv_hp = @(pval,kval) local_chi2inv_vpa(pval, kval);

% compute cdf at df = p1, then inverse at df = p
% guard numerical edge cases: ensure CDF is strictly in (0,1)
% c = chi2cdf_hp(d2(i), p1);


% if c>0.9999999999999
%     c=0.9999999999999;
% end
% c = min(max(c, eps), 1 - eps); % avoid 0 or 1
% try
% d2_adj(i) = chi2inv_hp(c, p);
% catch
%     dd=1;
% end