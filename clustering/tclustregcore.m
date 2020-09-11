function [bopt,sigma2opt,nopt,postprobopt,muXopt,sigmaXopt,vopt,subsetopt,idxopt,webeta,cstepopt,webetaopt,Beta_all, obj_all]...%,retained_idopt
    =tclustregcore(y,X,RandNumbForNini,reftol,refsteps,mixt,equalweights,h,nselected,k,restrfact,restrfactX,alphaLik,alphaX,...
    seqk,NoPriorNini,msg,C,intercept,cwm,wtype_beta,we,wtype_obj,zigzag,wei,cup,pstar)

% This function is called by tclustregeda and it is not intended to be called directly

% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

% Groups with less than skipthin_th units are not considered for thinning
skipthin_th = 5;

%monitors evolution of parameters in refining steps (just in case of
%intercept through the origin and one regression coefficient)
monitor = 0;

% Prior estep before entering concentration steps loop
warmingup=true;

[n,p]=size(X);
sigma2ini = ones(1,k);

% bopt     = beta parameters of each group obtained in the best subset
bopt      = zeros(p,k);
% sigma2opt= sigma parameters of each group obtained in the best subset
% sigma2opt = NaN(1,k);

%initializations for clusterwise regression
if cwm==1
    % muX = matrix of size (p-intercept)-by-k which contains centroids of
    % each extracted subset
    % sigmaX = matrix of size (p-intercept)-by-(p-intercept)-k which
    % contains covariance matrices of each extracted subset
    muX       = zeros(k,p-intercept);
    sigmaX    = zeros(p-intercept,p-intercept,k);
    
    % U = 3D array which contains the eigenvectors of covariance matrix of
    % the explanatory variables for each group
    U         = sigmaX;
    
    % Lambda_pk = matrix which contains the eigenvalues of sigmaX referred
    % to the gruops (first column is associated with group 1....)
    Lambda_pk = muX';
end

%wbetaopt = vector of {0,1} weights with 1 to identify units that are not
%thinned.  The vector is the one giving rise to the optimal solution.
%The intermediate vector, wbeta, is used in the first trimming step to
%compute the likelihood contribution on the units that are not thinned.
webetaopt = ones(n,1);

%wobjopt = vector of {0,1} weights with 1 to identify units that are not
%thinned.  The vector is the one giving rise to the optimal solution.
%The intermediate vector, wbeta, is used in the first trimming step to
%compute the likelihood contribution on the units that are not thinned.
%weobjopt = ones(n,1);

%Remark. The assignement of the thinned units can be found with:
%A = [find(wbetaopt==0)' , idxopt(find(wbetaopt==0))']
%where the first column of A is the row index of the thinned unit and the
%second one is the cluster assignement.

%idxopt = (nx1) vector of values in {1, ... , k, -1 , -2} with cluster
%assignements (1,...,k) and trimmed units (-1 for first level trimming and
%-2 for second level trimming).
idxopt    = zeros(n,1);


% nopt = (kx1) vector containing the number of observations in each of
% the k groups in the optimal subset. The trimmed units are not counted.
nopt       = zeros(1,k);

% postprob = (nxk) matrix of posterior probabilities of the obtimal subset
postprobopt = zeros(n,k);

if monitor == 1
    Beta_all   = NaN(k,refsteps,nselected);
    obj_all    = NaN(nselected, refsteps);
else
    Beta_all=[];
    obj_all =[];
end

% penalized objective function
penal_obj = 0;

% current and best objective function values
vopt = -1e+20;
muXopt    = NaN(k,p-intercept);
%sigmaXopt =  X sigma.
sigmaXopt = NaN(p-intercept,p-intercept,k);

% index of the best concentration step
cstepopt = 0;

% verLess2016b is true if current version is smaller than 2016b
verLess2016b=verLessThanFS(9.1);

if verLess2016b == true
    userepmat=1;
else
    userepmat=2;
end

% tolrestreigen = tolerance to use in function restreigen
tolrestreigen=1e-08;

% noconv = scalar linked to the number of times in which there was no
% convergence
% noconv=0;


% if mixt>=1
%     % log_lh = h-by-k matrix containing the log of component conditional
%     %          density weighted by the component probability.
%     % log_lh = log( \pi_j \phi (y_i; \; \theta_j))
%     log_lh=zeros(h,k);
% end

% timer to monitor performances.
tsampling = ceil(min(nselected/5 , 50));
time      = zeros(tsampling,1);

% ll      = loglikelihood for each observations in each group
ll        = zeros(n,k);

%%  RANDOM STARTS

for i =1:nselected
    
%     switch msg
%         case 1
%             % monitor time execution
%             if msg==1
%                 if i <= tsampling
%                     tstart = tic;
%                 end
%             end
%         case 2
%             % monitor iteration step
%             if msg
%                 disp(['Iteration ' num2str(i)]);
%             end
%     end
    
    
    % ltkg = becomes 1 if a particular subset leads to less than k groups
    ltkg=0;
    
    % Beta = matrix of beta coefficients (j-col refers to j-th group)
    Beta = zeros(p,k);
    
    % to replicate results, fix the seed by uncommenting the line below
    % rng(1234);
    
    %%% -- Initialization of mixing proportions
    if NoPriorNini==1
        % if initial mixing proportions are not supplied by the user, they
        % are randomly assigned, making sure that:
        % a) minimum group size is positive
        % b) sum of group sizes is equal to h
        niin=0;
        itermax=0;
        while min(niin)==0 && itermax<1000
            randk=rand(k,1);
            niin=floor(h*randk/sum(randk));
            diffh=sum(niin)-h;
            [~,imin]=min(niin);
            niin(imin)=niin(imin)-diffh;
            itermax=itermax+1;
        end
        if itermax ==1000
            error('FSDA:tclustregcore:WrongInput','Initialization of the group proportions failed')
        end
        niini=niin;
    else
        %check the initial mixing proportions supplied by the user
        randk=RandNumbForNini(:,i);
        itermax=0;
        
        % Make sure that the sum of niin is h
        niin=floor(h*randk/sum(randk));
        if sum(niin)<h
            niin(1)=niin(1)+h-sum(niin);
        end
        % Make sure that minimum group size is strictly positive
        while min(niin)==0 && itermax<1000
            ntoreplace=seqk(niin==0);
            niin(ntoreplace(1))=1;
            [~,posmax]=max(niin);
            niin(posmax)=niin(posmax)-1;
            itermax=itermax+1;
        end
        niini=niin;
    end
    
    %%% -- Initial regression parameters estimation
    
    % Estimate the regression parameters Beta and the Var-Cov matrix for
    % each of the k *initial* groups provided by the user or randomly
    % generated.
    
    %index: the i-th line of C, containing the i-th subset units
    index  = C(i,:);
    for j = 1:k
        %selj: the units of the current subset i, belonging to the group j
        ilow   = (j-1)*p+1;
        iup    = j*p;
        selj   = index(ilow:iup);
        
        %Xb and yb: X and y of the current subset i belonging to the group j
        Xb     = X(selj,:);
        yb     = y(selj,:);
        
        %If the model is without intercept and Xb is zero, the regression
        %cannot be computed. Some jittering is applied to Xb
        if intercept == 0 && isscalar(Xb) && Xb == 0
            Xb = Xb + 0.0001*abs(randn(1,1));
        end
        
        %Beta and sigma2ini: estimation of regression parameters
        Beta(:,j) = Xb\yb;
        if length(yb)==1
            sigma2ini(j)= var(y);
        else
            sigma2ini(j) =var(yb);
        end
        
        % clusterwise regression: initialize X-distribution parameters
        if cwm==1
            
            % muX = matrix of size (p-intercept)-by-k which contains the k
            % centroids of the current subset i
            % sigmaX = matrix of size (p-intercept)-by-(p-intercept)-k which
            % contains covariance matrices of the current subset i
            muX(j,:)=mean(Xb(:,intercept+1:end),1);
            if length(yb)==1
                sigmaX(:,:,j)=var(X);
            else
                sigmaX(:,:,j)=cov(Xb(:,intercept+1:end));
            end
            
            % Lambdaj and Uj: eigenvalue eigenvector decomposition in the
            % current subset i, for group j
            [Uj,Lambdaj] = eig(sigmaX(:,:,j));
            % Store eigenvectors and eigenvalues of group j
            U(:,:,j)=Uj;
            Lambda_pk(:,j)=diag(Lambdaj);
        end
        
    end
    
    %%% -- Application of the eigenvector-eigenvalue restriction
    % $$ \frac{ max_{g=1,\ldots,G} \sigma_g^2 \pi_g}{min_{g=1,\ldots,G} \sigma_g^2 \pi_g} < restrfact $$
    
    %restrict the eigenvalues according to the constraint specified in restrfact:
    if equalweights==1
        sigma2ini= restreigen(sigma2ini,ones(k,1),restrfact,tolrestreigen,userepmat);
        if cwm==1
            autovalues= restreigen(Lambda_pk,ones(k,1),restrfactX,tolrestreigen,userepmat);
        end
    else
        sigma2ini= restreigen(sigma2ini,niini,restrfact,tolrestreigen,userepmat);
        if cwm==1
            autovalues= restreigen(Lambda_pk,niini,restrfactX,tolrestreigen,userepmat);
        end
    end
    
    % CWM: re-construct the covariance matrices keeping into account the
    % constraints on the eigenvalues
    if cwm==1
        for j=1:k
            sigmaX(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
            % Alternative code: more elegant but slower because diag is a
            % built in function
            % sigmaX(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
        end
    end
    
    if warmingup==true
        %%% -- Classification *E-Step* (before concentration steps):
        %%% -- -- Log-likelihood of all observations based on the estimated regression parameters
        
        % equalweight: each group has the same weight, $1/k$.
        if equalweights == 1
            for jj = 1:k
                ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                if cwm==1
                    ll(:,jj)=ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj));
                end
            end
            % the group weights (niini) are estimated
        else
            for jj = 1:k
                ll(:,jj) = log((niini(jj)/sum(niini))) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                if cwm==1
                    ll(:,jj)=ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj));
                end
            end
        end
        
        %%% -- -- Posterior probabilities of all observations
        
        % Mixture likelihood model
        
        if mixt > 0
            %E-step is run to compute the posterior probabilities of all
            %observations
            [~,postprob] = estepFS(ll,verLess2016b);
            % idx: (nx1) vector indicating the group for which each observation
            % has the largest posterior probability. It takes values in  {1,
            % ... , k}; at the end of the algorithm it will take values in {1,
            % ... , k,0, -1 , -2}, respectively for group assignement, thinned
            % units, first and second trimmed units.
            [~,idx]= max(postprob,[],2);
            
            %classification likelihood model
        else %  mixt == 0
            zeronk=zeros(n,k);
            % idx: (nx1) vector indicating the group for which each observation
            % has the largest likelihood. It takes values in  {1, ... , k};
            % At the end it will take values in {1, ... , k,0, -1 ,
            % -2}, respectively for group assignement, thinned units, first and
            % second trimmed units.
            [~,idx] = max(ll,[],2);
            postprob = zeronk;
            for j=1:k
                postprob(idx==j,j)=1;
            end
        end
    end
    
    
    %% -- CONCENTRATION STEPS
    
    indold = zeros(n,1)-1;
    postprobold = zeros(n,k);
    mudiff=1e+15;
    cstep=0;
    while mudiff > reftol && cstep <= refsteps-1
        cstep = cstep + 1;
        
        
        if warmingup==false
            %%% -- Classification *E-Step* (inside concentration steps):
            %%% -- -- Log-likelihood of all observations based on the estimated regression parameters
            
            % equalweight: each group has the same weight, $1/k$.
            if equalweights == 1
                for jj = 1:k
                    ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                    if cwm==1
                        ll(:,jj)=ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj));
                    end
                end
                % the group weights (niini) are estimated
            else
                for jj = 1:k
                    ll(:,jj) = log((niini(jj)/sum(niini))) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                    if cwm==1
                        ll(:,jj)=ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj));
                    end
                end
            end
            
            %%% -- -- Posterior probabilities of all observations
            
            % Mixture likelihood model
            
            if mixt > 0
                %E-step is run to compute the posterior probabilities of all
                %observations
                [~,postprob,disc] = estepFS(ll,verLess2016b);
                % idx: (nx1) vector indicating the group for which each observation
                % has the largest posterior probability. It takes values in  {1,
                % ... , k}; at the end of the algorithm it will take values in {1,
                % ... , k,0, -1 , -2}, respectively for group assignement, thinned
                % units, first and second trimmed units.
                [~,idx]= max(postprob,[],2);
                
                %classification likelihood model
            else %  mixt == 0
                zeronk=zeros(n,k);
                % idx: (nx1) vector indicating the group for which each observation
                % has the largest likelihood. It takes values in  {1, ... , k};
                % At the end it will take values in {1, ... , k,0, -1 ,
                % -2}, respectively for group assignement, thinned units, first and
                % second trimmed units.
                [disc,idx] = max(ll,[],2);
                postprob = zeronk;
                for j=1:k
                    postprob(idx==j,j)=1;
                end
            end
        end
        
        %%% -- -- Observation weighting, according to wtrim option (*componentwise thinning* is applied here)
        
        % wbeta = vector of weights in {0,1}; wbeta = 1 identifies units
        % that are not thinned and contribute to the likelihood in the
        % first trimming step.
        
        switch wtype_beta
            
            case 0
                %%% -- -- * case 0: no observation weighting: wbeta is constant
                webeta = ones(n,1);
                weobj = ones(n,1);
                
            case 1
                %%% -- -- * case 1: wbeta contains user-specific weights
                webeta = we;
                if strcmp(wtype_obj,'user')
                    weobj =  webeta;
                else
                    weobj = ones(n,1);
                end
                
            case 2
                %%% -- -- * case 2: Componentwise Thinning.
                % wbeta contains the retention probabilities on yhat,
                % estimated with FSDA function wthin.
                
                % kx1 vector for groups that are too small to be thinned
                small_group     = zeros(k,1);
                
                webeta = ones(n,1);
                weobj = ones(n,1);
                for jj=1:k
                    % Boolean index of units in group j
                    groupj=idx==jj;
                    if  sum(groupj) > skipthin_th
                        % update wbeta if the group has more than
                        % skipthin_th units.
                        Xj   = X(groupj,:);
                        yhat = Xj*Beta(:,jj);
                        % thinning
                        [Zt , pretain]  = wthin(yhat,'cup',cup,'pstar',pstar);
                        webeta(groupj) = pretain;
                        if strcmp(wtype_obj,'0')
                            weobj(groupj) = ones(sum(groupj),1);
                        elseif strcmp(wtype_obj,'Z')
                            weobj(groupj) = Zt;
                        elseif strcmp(wtype_obj,'w')
                            weobj(groupj) = pretain;
                        elseif strcmp(wtype_obj,'wZ')
                            weobj(groupj) = pretain .* Zt;
                        else
                            error('FSDA:tclustregcore:WrongInput','wtype_obj option not correct')
                        end
                        % pretain: the retention probabilities are based on
                        % the predicted values (yhat) estimated at the
                        % previous concentration step. REMARK: trimmed and
                        % non-trimmed units are both considered.
                    else
                        % The group is too small: skip thinning
                        if sum(groupj) > 0
                            small_group(jj) = 1;
                        end
                    end
                    
                end
                % For not-empty groups that are too-small to run thinning,
                % wbeta is the median of the weights of the other groups.
                if sum(small_group) > 0
                    % REMARK: was nanmedian
                    medianweights  = median(webeta(ismember(idx,find(small_group==0))));
                    webeta(groupj) = medianweights;
                    if strcmp(wtype_obj,'0')
                        weobj(groupj) = ones(sum(groupj),1);
                    elseif strcmp(wtype_obj,'Z')
                        weobj(groupj) = ones(sum(groupj),1);
                    elseif strcmp(wtype_obj,'w')
                        weobj(groupj) = medianweights;
                    elseif strcmp(wtype_obj,'wZ')
                        weobj(groupj) = medianweights;
                    else
                        error('FSDA:tclustregcore:WrongInput','wtype_obj option not correct')
                    end
                    
                end
                
            case 3
                %%% -- -- * case 3: Componentwise Thinning.
                
                % weights are the posterior probabilities multiplied by
                % the bernoulli weights
                
                webeta = ones(n,1);
                weobj = ones(n,1);
                % find indices of units in group jj
                for jj=1:k
                    %The first concentration step uses idx which takes
                    %values in {1, ... , k}.
                    %The following concentration steps, where idx takes
                    %values in {1, ... , k, 0, -1 , -2}, use idx_ne0 which
                    %takes values  {1, ... , k, -1 , -2}
                    
                    if cstep == 1
                        ijj = idx==jj;
                    else
                        ijj = idx_ne0==jj;
                    end
                    
                    % update wbeta if the group has more than
                    % skipthin_th units.
                    if  sum(ijj)> skipthin_th
                        % Bernoulli weights based on density estimated
                        % on the component predicted values in the
                        % previous step.
                        Xj          = X(ijj,:);
                        yj          = y(ijj,:);
                        yhat        = Xj*Beta(:,jj);

                        [Zt , ~]    = wthin(yhat,'cup',cup,'pstar',pstar);%,'bandwidth',0.9);
                        tmp(jj,1)=sum(ijj);
                        tmp(jj,2)=sum(Zt==0);
                        webeta(ijj) = Zt;
                        if strcmp(wtype_obj,'Z')
                            weobj(ijj) = Zt;
                        elseif strcmp(wtype_obj,'0')
                            weobj(ijj) = ones(length(sum(ijj)),1);
                            
                        else
                            error('FSDA:tclustregcore:WrongInput','wtype_obj option not correct')
                        end
                        % count the thinned observations
                        % nthinned = nthinned + sum(Wt == 0);
                        % REMARK: group weights do not consider the thinned
                        % units. To be discussed.
                        niini(jj) = sum(Zt > 0);
                        
                    end
                    
                end
                
            case 4
                %%% -- -- * case 4: Bivariate tandem thinning
                %wbeta = Wt4;
                webeta = ones(n,1);
                weobj = ones(n,1);
                
            case 5 % TO BE IMPLEMENTED
                
            case 6 % TO BE IMPLEMENTED
                
        end
        if monitor == 1
            Beta_all(1:k,cstep) = Beta';
        end
        
        % Mean of the weights must be 1.
        mean_wbeta = mean(webeta);
        webeta      = webeta/mean_wbeta;
        
        mean_wobj = mean(weobj);
        weobj      = weobj/mean_wobj;
        
        if sum(isnan(Beta(:)))>0
            %%% -- -- _Stop if one of the current beta parameters is undefined, and go to another subset_
            break
        end
        
        
        if warmingup==true
            %% -- -- Classification *E-Step* (inside concentration steps)
            
            %%% -- -- -- Log-likelihood of all observations based on the estimated regression parameters
            
            % equalweight: each group has the same weight, 1/k.
            if equalweights == 1
                for jj = 1:k
                    ll(:,jj) = log((1/k)) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                    if cwm==1
                        ll(:,jj)=  ll(:,jj)+ logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj));
                    end
                end
                % the group weights (niini) are estimated
            else
                for jj = 1:k
                    ll(:,jj) = log((niini(jj)/sum(niini))) + logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj));
                    if cwm==1
                        ll(:,jj)=  ll(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj));
                    end
                end
            end
            
            %%% -- -- -- Posterior probabilities of all observations
            
            % Mixture likelihood model
            if mixt == 2
                % Store previously computed posterior probabilities
                postprobold=postprob;
                
                %E-step is run to compute the posterior probabilities of all
                %observations
                [~,postprob,disc] = estepFS(ll,verLess2016b);
                % idx: (nx1) vector indicating the group for which each observation
                % has the largest posterior probability. It takes values in  {1,
                % ... , k}; at the end of the algorithm it will take values in {1,
                % ... , k,0, -1 , -2}, respectively for group assignement, thinned
                % units, first and second trimmed units.
                [~,idx]= max(postprob,[],2);
                
                %classification likelihood model
            else
                
                % idx: (nx1) vector indicating the group for which each
                % observation has the largest log-likelihood as in point
                % 2.1 of appendix of Garcia Escudero et al. It takes values
                % in  {1, ... , k}; at the end of the algorithm it will
                % take values in {1, ... , k,0, -1 , -2}, respectively for
                % group assignement, thinned units, first and second
                % trimmed units.
                [disc,idx] = max(ll,[],2);
                
                postprob=zeronk;
                for j=1:k
                    postprob(idx==j,j)=1;
                end
            end
        end
        
        % Assign infinite weights to thinned units so that they can
        % never be trimmed
        disc(webeta==0)=1e+20;
        
        % Sort the n likelihood contributions and save them in qq
        [~,qq] = sort(disc,'ascend');
        
        % remember to which group the thinned units belong.
        if wtype_beta == 3
            % idx_ne0, computed only for wtrim = 3, is a (nx1) vector of
            % cluster assignements, taking values in  {1, ... , k}; at the
            % end it will take values in {1, ... , k,  -1 , -2}.
            
            % idx is a (nx1) vector of cluster assignements, taking values
            % in  {1, ... , k, 0}; at the end it will take values in {1,
            % ... , k, 0, -1 , -2}.
            % idx could be always obtained as idx=idx_ne0.*wbeta
            
            idx_ne0=idx;
            idx(webeta==0)=0;
        end
        
        
        %% -- -- First level trimming
        
        % Order the weights according to qq
        wbetaordered=webeta(qq);
        
        % Find cumlative sum of weights
        cumsumww = cumsum(wbetaordered);
        
        % qqunassigned_small = a n-by-1 Boolean vector
        % containing true for the units which
        % have to be trimmed)
        qqunassigned_small = cumsumww <= n*alphaLik;
        
        % qqunassigned = indexes of units subject to first
        % level trimming
        qqunassigned = qq(qqunassigned_small);
        
        %indmax_before_tr should be saved in order to be able
        %to understand which groups the two trimmings affect more
        %indmax_before_tr = idxtt; TO DELETE
        
        
        % idx takes values in {1, ... , k, 0,  -1}. At the end it will
        % takes values in {1, ... , k, 0,  -1, -2}.
        idx(qqunassigned)=-1;
        if wtype_beta==3
            % idx_ne0 takes values in {1, ... , k,  -1}. At the end it
            % will takes values in {1, ... , k,  -1, -2}.
            idx_ne0(qqunassigned)=-1;
        end
        % postprob = vector nxk taking values in {(0 1], 0}:, (0 1]
        % are the posterior probabilities referring to not trimmed
        % units and 0 to trimmed units.
        postprob(qqunassigned,:)=0;
        
        % Z = vector nxk taking values in {[0 1], 0}: (0 1] are the
        % posterior probabilities referring to not trimmed and not
        % thinned units, 0 to trimmed or thinned units.
        Z=bsxfun(@times,postprob,webeta);
        
        %% -- -- Second level of trimming or CWM
        
        % FS or MCD are used to find units to trim
        
        %alphaX < 1: is the the fixed (<0.5) or adaptive (>0.5)
        %proportion of units to trim
        if alphaX<1 && alphaX>0
            for jj=1:k
                % groupjind = indices of units belonging to group j
                groupjind = find(idx==jj);
                %Xjnointercept = X of units belonging to group j
                Xjnointercept  = X(groupjind,intercept+1:end);
                %Xjnointercept = X and possible intercept of units belonging to group j
                Xj=X(groupjind,:);
                %number of units in group j
                njj=size(Xj,1);
                
                % hj = number of units to retain for group j after
                % second trimming level.
                if alphaX>0.5
                    %For adaptive trimming hj is the number of units in group j
                    hj=njj;
                else
                    %For (traditional) fixed trimming hj is the number of
                    %units in group j - trimmed units
                    hj = floor(njj*(1-alphaX));
                end
                
                %check number of units in absolute terms: (1) the size
                %of group after second level trimming is not smaller
                %than p and (2) the size of a group before second level
                %trimming is at least greater than p+2 (necessary to
                %run the FS)
                if hj >=p && njj>p+2
                    
                    % check number of units in relation to the number
                    % of variables
                    if njj/(p-intercept)>10
                        
                        %The MCD is applied only when p=1, because in
                        %this case it is faster than the FS.
                        if p-intercept==1
                            
                            %R2016a has introduced robustcov, which could be used here as below.
                            %Remember however that mcd returns the squared distances, i.e. RAW.md = mah.^2.
                            %{
                                [~,~,mah,~,~] = robustcov(Xjnointercept,'Method','fmcd','NumTrials',nsampmcd,'OutlierFraction',alpha2b,'BiasCorrection',1); %
                                plot(1:ni(jk),mah.^2,'-r',1:ni(jk),RAW.md,'-b');
                                close;
                            %}
                            
                            if alphaX>0.5
                                %In adaptive second trimming step, units
                                %outside the confidence bands are trimmed.
                                [~,REW]      = mcd(Xjnointercept,'msg',0,'conflev',1-(1-alphaX)/njj,'betathresh',1);
                                if isfield(REW,'outliers')
                                    trimj=REW.outliers;
                                else
                                    trimj = [];
                                end
                                
                            else
                                %In traditional fixed second trimming step,
                                %units with the largest njj-hj residuals
                                %are trimmed.
                                [~,REW]      = mcd(Xjnointercept,'msg',0,'betathresh',1);
                                
                                if isfield(REW,'md')
                                    [~,indmdsor] = sort(REW.md);
                                    trimj=indmdsor(hj+1:end);
                                else
                                    trimj = [];
                                end
                            end
                        else
                            % With multiple explanatory variables, the
                            % Forward Search is faster than MCD.
                            
                            if alphaX>0.5
                                %In adaptive second trimming step, units
                                %outside the confidence bands are trimmed.
                                outj=FSM(Xjnointercept,'nocheck',1,'init',round(njj*0.9),'msg',0,'bonflev',alphaX,'plots',0);
                                
                                trimj=outj.outliers;
                                if isnan(trimj)
                                    trimj=[];
                                end
                            else
                                %In traditional fixed second trimming step,
                                %units with the largest njj-hj residuals
                                %are trimmed.
                                [~,BBsel]=FSMbsb(Xjnointercept,0,'bsbsteps',hj,'init',hj,'nocheck',1,'msg',0);
                                seqj=1:njj;
                                % BBsel contains a NAN for the units
                                % not belonging to subset in step hj
                                trimj=seqj(isnan(BBsel));
                            end
                        end
                    else
                        trimj=[];
                    end
                    
                    % idx takes values in {1, ... , k, 0,  -1,-2}.
                    idx(groupjind(trimj))=-2;
                    if wtype_beta==3
                        % idx_ne0 takes values in {1, ... , k, -1,-2}.
                        idx_ne0(groupjind(trimj))=-2;
                    end
                else
                    %%% -- -- _Stop if we end up with less than k groups, and go to another subset_
                    ltkg = 1;
                    break
                end
            end
            
            % Z = vector nxk taking values in {[0 1], 0}: (0 1] are the
            % posterior probabilities referring to not trimmed and not
            % thinned units, 0 to trimmed or thinned units.
            Z(idx==-2,:)=0;
            
            % postprob = vector nxk taking values in {(0 1], 0}:, (0 1]
            % are the posterior probabilities referring to not trimmed
            % units and 0 to trimmed units.
            postprob(idx==-2,:)=0;
        end
        
        %%% -- -- _Stop if we end up with less than k groups, and go to another subset_
        if ltkg==1
            Beta=NaN; %#ok<NASGU>
            break
        end
        
        %% -- -- *M-Step*: find beta coefficients and sigma2 using weighted regression
        for jj=1:k
            
            %sqweights = weights (for beta estimation) of observations
            %belonging to group jj. Note that the weights of the
            %trimmed and thinned units are 0, while the weights of the
            %other observations are the sqrt of the posterior
            %probabilities.
            sqweights = sqrt(Z(:,jj));
            % nj = sum of the weights (posterior probabilities of not
            % trimmed and not thinned observations) for group jj.
            nj=sum(Z(:,jj));
            % ninini = kx1 vector containing the sum of the weights for each group
            niini(jj)=nj;
            
            %check if current groups is large enought to do estimation
            if nj>p+1
                %multiple X and y for the observations weights (posterior
                %probabilities of not trimmed and not thinned observations)
                % sqweights=sqweights.*(X(:,end).^0.585);
                Xw = bsxfun(@times, X, sqweights);
                yw = y .* sqweights;
                
                % breg = estimate of beta of group jj from (re)weighted regression
                % (RWLS)
                breg = Xw\yw;
                % outp=regressH(yw,Xw,X(:,end),'maxiter',10);
                % breg=outp.Beta(1:end-1,1);
                % Beta = estimate of beta of all group from (re)weighted regression
                % (RWLS).
                Beta(:,jj)=breg;
                
                % sigma2 = estimate of sigma2 of group jj after
                % weighted regression.
                res2=(yw-Xw*breg).^2;
                sigma2=sum(res2)/nj;
                
                % sigma2ini = estimate of sigma2 of all group from
                % (re)weighted regression (RWLS).
                sigma2ini(jj) = sigma2;
                
                %CWM
                if cwm ==1
                    %muX = estimate of centroids of X
                    muX(jj,:)=sum(bsxfun(@times, X(:,intercept+1:end), Z(:,jj)),1)/sum(Z(:,jj))';
                    %sigmaX = estimate of sigma of X
                    sigmaX(:,:,jj)= (Xw(:,intercept+1:end)'*Xw(:,intercept+1:end))/sum(Z(:,jj))-muX(jj,:)'*(muX(jj,:));
                end
            else
                %%% -- -- _Stop if we end up with less than k groups_
                ltkg=1;
                break
            end
        end % loop on groups
        
        %%% -- -- _Stop if we end up with less than k groups_
        if ltkg==1
            Beta=NaN; %#ok<NASGU>
            break
        end
        
        %%% -- -- Application of the eigenvector-eigenvalue restriction
        
        % CWM
        if cwm==1
            for j=1:k
                % Eigenvalue eigenvector decomposition for group j
                [Uj,Lambdaj] = eig(sigmaX(:,:,j));
                % Store eigenvectors and eigenvalues of group j
                U(:,:,j)=Uj;
                Lambda_pk(:,j)=diag(Lambdaj);
            end
        end
        
        % equalweights =1: equal proportions are supplied for group sizes
        if equalweights==1
            sigma2ini= restreigen(sigma2ini,ones(k,1),restrfact,tolrestreigen,userepmat);
            
            if cwm==1
                autovalues= restreigen(Lambda_pk,ones(k,1),restrfactX,tolrestreigen,userepmat);
            end
            % equalweights = 0: proportions are set equal to the group sizes
        else
            sigma2ini= restreigen(sigma2ini,niini,restrfact,tolrestreigen,userepmat);
            if cwm==1
                autovalues= restreigen(Lambda_pk,niini,restrfactX,tolrestreigen,userepmat);
            end
            
        end
        
        % CWM
        if cwm==1
            % Covariance matrices are reconstructed keeping into account the
            % constraints on the eigenvalues
            for j=1:k
                sigmaX(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
                
                % Alternative code: in principle more efficient but slower
                % because diag is a built in function
                % sigmaX(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
            end
        end
        
        
        %% -- -- Computation of the value of the target function (Decide if thinned units have to enter or not in the obj)
        obj = 0;
        
        %identify possible empty groups
        not_empty_g = seqk(~( niini <= p + 1 ));
        
        %%% -- -- -- Classification likelihood
        if mixt == 0
            
            %loop on not-empty groups
            for jj = not_empty_g
                
                % resjj = residuals (weighted) of all units from group jj
                resjj=(y-X*Beta(:,jj))./wei;
                
                % equalweights =1: equal proportions are supplied for group sizes
                if equalweights ==1
                    %the next if-then-else statement is an experiment to
                    %decide if the thinned units must be considered or
                    %not in the obj function. The same if-then-else is
                    %not present for equalweights == 0
                    
                    %In case of bernoulli thinning, if the flag
                    %obj_with_thinned == 1, the objective function is
                    %computed excluding only trimmed units, by using
                    %matrix postprob. REMARK: when classification
                    %likelihood is used (mixt = 2), postprob is a
                    %boolean matrix.
                    
                    %postprob: row i-column j:
                    % ---1 if observation i belongs to group j,
                    % ---0 if observation i belongs to another group than j or if it is trimmed
                    % Note that thinned observations are assigned to
                    % the group for which the log-likelihood is highest
                    
                    %weobj depends from the choice of the input
                    %parameter wtype_obj. The defaults is to assign to
                    %each observation weight 1.
                    
                    obj = obj + log(1/k) +...
                        sum(logmvnpdfFS(resjj,0,sigma2ini(jj)).*postprob(:,jj).*weobj(:)) ;
                    %sum_dens(jj,cstep) = sum(logmvnpdfFS(y-X*Beta(:,jj),0,sigma2ini(jj)));
                    
                    
                    % equalweights = 0: proportions are set equal to the group sizes.
                else
                    %the objective function is computed excluding both
                    %thinned and trimmed (matrix Z)
                    obj = obj + niini(jj)*log(niini(jj)/sum(niini)) +...
                        sum(logmvnpdfFS(resjj,0,sigma2ini(jj)).*postprob(:,jj).*weobj(:));
                end
                
                %CWM
                if cwm==1
                    %the objective function is computed excluding both
                    %thinned and trimmed (matrix Z)
                    obj=obj+sum(logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj)).*Z(:,jj));
                end
                
            end
            
            %%% -- -- -- Mixture likelihood (check penalization)
        else
            
            %initialization of log-likelihood
            log_lh=NaN(n,size(not_empty_g,2));
            
            for jj = 1:k
                
                % resjj = residuals (weighted) of all units from group jj
                resjj=(y-X*Beta(:,jj))./wei;
                
                % equalweights =1: proportions are set equal to the
                % group sizes.
                if equalweights ==1
                    log_lh(:,jj) = ...
                        log(1/k) + (logmvnpdfFS(resjj,0,sigma2ini(jj) ) );
                    
                    % equalweights = 0: equal proportions are supplied for
                    % group sizes.
                else
                    log_lh(:,jj) = ...
                        log(niini(jj)/sum(niini)) + (logmvnpdfFS(...
                        resjj,0,sigma2ini(jj) ) );
                end
                %CWM
                if cwm==1
                    log_lh(:,jj)=log_lh(:,jj)+logmvnpdfFS(X(:,(intercept+1):end),muX(jj,:),sigmaX(:,:,jj));
                end
            end
            
            %the objective function is computed excluding both thinned
            %and trimmed.
            log_lh(idx<=0,:)=[];
            obj = estepFS(log_lh,verLess2016b);
        end
        
        if penal_obj == 1
            % penalized term (to be corrected)
            ptermAIC = 2*nParam/mean_wbeta;
            %Penalization of objective function
            obj = 2*obj - ptermAIC;
        end
        if monitor ==1
            %obj_all = (nselected x refsteps) vector containing the obj values
            %of all subsets and all concentration steps.
            obj_all(i,cstep) = obj;
            
            Beta_all(1:k,cstep) = Beta;
        end
        
        
        %%% -- -- Update the 'optimal' target value and parameters
        
        if mixt>0
            % if mixt >0 stopping criterion is referred to postprob
            mudiff=sum(sum(abs(postprob-postprobold)))/n;
            % disp(mudiff)
        else
            % if mixt=0 stopping criterion is referred to no modiification in the classification
            mudiff=sum(abs(indold-idx)>0)/n;
            % disp(mudiff)
        end
        
        
        % This is done every step if the function is non-monotonic (i.e.
        % when second-level triming and/or thinning are applied);
        % otherwise, this is done once at the final step.
        
        if zigzag == 1 || mudiff<reftol || (zigzag == 0  && cstep == refsteps) || vopt==-1e+20
            
            % obj >= vopt: to check an increase in the target value.
            if obj >= vopt && sum(isnan(Beta(:))) ==0

                %cstepopt = cstep with the largest obj
                cstepopt  = cstep;
                %subsetopt = subset with the largest obj
                subsetopt  = i;
                %vopt     = value of obj in the optimal cstep
                vopt      = obj;
                %bopt     = value of regression parameters in the optimal cstep
                bopt      = Beta;
                %nopt     = value of groups sizes in the optimal cstep
                %(trimming units exluded)
                nopt      = niini;
                %sigma2opt= value of groups variances in the optimal cstep
                sigma2opt  = sigma2ini;
                %wbetaopt= value of weights (thinning probabilities) in
                %the optimal cstep
                webetaopt  = webeta;
                if wtype_beta==5
                    %w4trimopt_obj_5  = w4trim_obj_5;
                    %TO BE IMPLEMENTED
                end
                %idxopt = (nx1) vector of {-1,-2, 1, ..., k}
                %decomment the next lines if in idx you want the id of thinned
                %units
                %                 if wtype_beta==3
                %                     idxopt = idx_ne0;
                %                 else
                %                     idxopt = idx;
                %                 end
                idxopt = idx;
                %retained_idopt = idx;
                % postprob = vector nxk taking values in {(0 1], 0}:, (0 1]
                % are the posterior probabilities referring to not trimmed
                % units and 0 to trimmed units.
                if mixt == 2
                    postprobopt = postprob;
                else
                    [~,postprobopt,~] = estepFS(ll,verLess2016b);
                    postprobopt(idx<=0,:)=0;
                end
                
                %CWM
                if cwm==1
                    %muXopt = X centroids.
                    muXopt    = muX;
                    %sigmaXopt =  X sigma.
                    sigmaXopt = sigmaX;
                end
            end
            
        end % End of zigzag if
        
        
        
    end % Loop on the concentration steps
    
    %% -- END OF CONCENTRATION STEPS
    
    % monitor time execution
    if msg==1
        if i <= tsampling
            % sampling time until step tsampling
            time(i)=toc(tstart);
        elseif i==tsampling+1
            % stop sampling and print the estimated time
            fprintf('Total estimated time to complete tclustreg: %5.2f seconds \n', nselected*median(time));
        end
    end
end % end of loop over the nsamp subsets

%%  END OF RANDOM STARTS

chkExistence=exist('sigma2opt','var');
if chkExistence==0
    sigma2opt=NaN;
    nopt=NaN;
    postprobopt=NaN;
    muXopt=NaN;
    sigmaXopt=NaN;
    % vopt
    subsetopt=NaN;
    idxopt=NaN;
    %retained_idopt=NaN;
    webeta=NaN;
    webetaopt=NaN;
    cstepopt=NaN;
end

end