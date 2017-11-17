function  [obj, varargout1, varargout2]=estepFS(log_lh)
%estepFS performs e-step for Gaussian mixture distribution
%
%<a href="matlab: docsearchFS('estepFS')">Link to the help function</a>
%
%   obj = estepFS(log_lh) returns value of the loglikelihood of mixture model.
%
%   obj is equal to
%   \begin{equation}\label{mixlik}
%   obj = \log   \left( \prod_{i=1}^n  \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)    \right).
%   \end{equation}
%
%   or
%
%   \begin{equation}\label{mixlik}
%   obj =  \sum_{i=1}^n  \log   \left( \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)    \right).
%   \end{equation}
%
%
%   k = number of components of the mixture
%   \pi_j = component probabilitites
%   \theta_j = parameters of the j-th component of the mixture
%
%  Required input arguments:
%
%       log_lh: n-by-p matrix containing the log of component conditional
%               density weighted by the component probability.
%               log_lh = log( \pi_j \phi (y_i; \; \theta_j))
%
%
%  Output:
%
%         obj  : scalar. Value of the log likelihood (see above) of
%                mixture models
%
%  Optional Output:
%
%   varargout1 : n-by-k matrix containing posterior probabilities
%                varargout1 i=1, ..., n and j=1, ..., k is the posterior
%                probability of observation i belonging to component j of
%                the mixture.
%                Posterior probabilities are computed as
%                varargout1(i,j)= exp( log_lh(i,j))/ (\sum_{j=1^k} exp(log_lh(i,j)))
%
%   varargout2 : n-by-1 vector which contains the contributions
%                to the loglikelihood of the mixture model of each observation
%                More precisely,
%                lodpdf = \log   \left( \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)    \right)
%                       = \log   \left( \sum_{j=1}^k  \exp( \log_lh ) \right)
%
%
%
% DETAILS. Formally a mixture model corresponds to the mixture distribution that
% represents the probability distribution of observations in the overall
% population. Mixture models are used
% to make statistical inferences about the properties of the
% sub-populations given only observations on the pooled population, without
% sub-population-identity information.
% Mixture modeling approaches assume that data at hand $y_1, ..., y_n$ in
% $R^p$ come from a probability distribution with density given by the sum of k components
% $\sum_{j=1}^k  \pi_j \phi( \cdot, \theta_j)$
% with $\phi( \cdot, \theta_j)$ being the
% $p$-variate  (generally multivariate normal) densities with parameters
% $\theta_j$, $j=1, \ldots, k$. Generally $\theta_j= (\mu_j, \Sigma_j)$
% where $\mu_j$ is the population mean  and   $\Sigma_j$ is the covariance
% matrix for component $j$.
% \pi_j is the (prior) probability of component j
%
% References:
%
%   McLachlan, G.J.; Peel, D. (2000). Finite Mixture Models. Wiley. ISBN 0-471-00626-2
%
% Copyright 2008-2017.
%
%<a href="matlab: docsearchFS('estepFS')">Link to the help function</a>
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%
%{
%      Generate two Gaussian normal distributions
%      and do not produce plots:
        MU1 = [1 2];
        SIGMA1 = [2 0; 0 .5];
        MU2 = [-3 -5];
        SIGMA2 = [1 0; 0 1];
        n1=100;
        n2=200;
        Y = [mvnrnd(MU1,SIGMA1,n1);mvnrnd(MU2,SIGMA2,n2)];
        k=2;
        PI=[1/3;2/3];
        MU=[MU1;MU2];
        SIGMA=zeros(2,2,2);
        SIGMA(:,:,1)=SIGMA1;
        SIGMA(:,:,2)=SIGMA2;
        ll=zeros(n1+n2,1);
            for j=1:k
                ll(:,j)= log(PI(j)) +  logmvnpdfFS(Y,MU(j,:),SIGMA(:,:,j));
            end
        % Compute posterior probabilities
        [~,postprob]=estepFS(ll);

%}

%% Beginning of code

maxll = max (log_lh,[],2);
% minus maxll to avoid underflow
post = exp(bsxfun(@minus, log_lh, maxll));

% density is a size(log_lh,1)-by-1 vector which contains the contribution
% of each observation to the likelihood  (each element is divided by the
% maximum of each row of input matrix log_lh to avoid overflow)
% density(i) is \sum_{j=1}^k \pi_j \phi(x_i| \theta_j)     / exp(maxll(i))
% i=1, 2, ..., size(log_lh,1)
density = sum(post,2);

if nargout>1
    %normalize posteriors
    % varargout1 is a size(log_lh,1)-by-p matrix which contains posterior
    % probabilities
    varargout1 = bsxfun(@rdivide, post, density);
end

% logpdf is a size(log_lh,1)-by-1 vector which contains the contributions
% to the loglikelihood of the mixture model of each observation
% Note that give that maxll was substracted now it must be added
% More precisely,
% lodpdf = \log   \left( \sum_{j=1}^k \pi_j \phi (y_i; \; \theta_j)    \right)
logpdf = log(density) + maxll;

if nargout>2
    varargout2=logpdf;
end

% obj is the value of the loglikelihood of the mixture model
% obj = \sum_{i=1}^n logpdf_i
obj = sum(logpdf) ;
end