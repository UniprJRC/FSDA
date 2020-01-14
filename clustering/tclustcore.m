function out  = tclustcore(Y,Cini,Sigmaini,Niini,reftol,refsteps,mixt,equalweights,h,nselected,k,restrnum,restrfactor,userepmat,nParam)
% This function is called by tclusteda and it is not intended to be called directly

% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

[n,v]=size(Y);

% callmex is a Boolean which is equal to true if the mex file exists
callmex=existFS('DfM');

% tolrestreigen = tolerance to use in function restreigen
tolrestreigen=1e-08;

% noconv = scalar linked to the number of times in which there was no
% convergence
noconv=0;

%Initialize the objective function (trimmed variance) by a
%large  value
vopt=-1e+30;
% fullsol = vector which stores value of the objective function in each
% iteration
fullsol=zeros(nselected,1);

persistent ll eyev Y0tmp ey onev1 Lambda_vk U

if isempty(ll)
    ll=zeros(n,k);
end

if isempty(Y0tmp)
    Y0tmp=zeros(n,v);
end

% Create an identity matrix which will be used in fucntion logmvnpdfFS
if isempty(eyev)
    eyev=eye(v);
end

% The covariances are given initially by k identity matrices
if isempty(ey)
    ey=eye(v,v);
end

if isempty(onev1)
    onev1=ones(v,1);
end

% Lambda_vk = matrix which will contain in column j the v (unrestricted)
% eigevalues of covariance matrix of group j (j=1, ..., k)
if isempty(Lambda_vk)
    Lambda_vk=ones(v,k);
end

if isempty(U)
    U=zeros(v,v,k);
end

% sigmaopt = 3 dimensional array which will contain the estimates of the
% covariance matrices for the best solution
sigmaopt=zeros(v,v,k);

if mixt>=1
    % log_lh = h-by-k matrix containing the log of component conditional
    %          density weighted by the component probability.
    % log_lh = log( \pi_j \phi (y_i; \; \theta_j))
    log_lh=zeros(h,k);
end

for i=1:nselected
    cini=Cini{i};
    sigmaini=Sigmaini{i};
    niini=Niini{i};
    iter=0;
    mudiff=1e+15;
    
    postprob=0;
    ind=0;
    
    % refsteps "concentration" steps will be carried out
    while ( (mudiff > reftol) && (iter < refsteps) )
        iter = iter + 1;
        if equalweights
            % In this case we are (ideally) assuming equally sized groups
            for j=1:k
                ll(:,j)= logmvnpdfFS(Y,cini(j,:),sigmaini(:,:,j),Y0tmp,eyev,n,v,0);
            end
        else
            
            % In this case we allow for different group weights or we are
            % assuming a mixture model
            for j=1:k
                ll(:,j)= log(niini(j)/h) +  logmvnpdfFS(Y,cini(j,:),sigmaini(:,:,j),Y0tmp,eyev,n,v,0);
                % Line above is faster but equivalent to
                % ll(:,j)= (niini(j)/h)*mvnpdf(Y,cini(j,:),sigmaini(:,:,j));
            end
            
        end
        
        if mixt==2
            
            postprobold=postprob;
            
            [~,postprob,disc]=estepFS(ll);
            
            % Sort the n likelihood contributions
            % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
            [~,qq]=sort(disc,'descend');
            
            % qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
            % (weighted) likelihood contributions
            qqunassigned=qq((h+1):n);
            qq=qq(1:h);
            
            % Ytri = n(1-alpha)-by-v matrix associated with the units
            % which have the largest n(1-alpha) likelihood contributions
            Ytri=Y(qq,:);
            
            postprob(qqunassigned,:)=0;
            
            % M-step update of niini
            % niini = numerator of component probabilities
            niini=(sum(postprob))';
            
        else
            indold=ind;
            
            % In this part we select the untrimmed units.
            % They are those which have the n(1-alpha) largest values among the
            % maxima of each row of matrix ll.
            % vector disc of length(n) contains the (weighted) contribution of
            % each unit to the log likelihood.
            [disc,ind]= max(ll,[],2);
            
            % Sort the n likelihood contributions
            % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
            [~,qq]=sort(disc,'descend');
            
            
            % qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
            % (weighted) likelihood contributions
            qqunassigned=qq(h+1:end);
            qq=qq(1:h);
            
            % try
            % Ytri = n(1-alpha)-by-v matrix associated with the units
            % which have the largest n(1-alpha) likelihood contributions
            Ytri=Y(qq,:);
            %catch
            %    jjj=1;
            %end
            
            % Ytriind = grouping indicator vector (of size n(1-alpha))
            % associated to Ytri
            groupind=ind(qq);
            
            % ind is the identifier vector
            % trimmed units have a value of ind=0
            ind(qqunassigned)=0;
        end
        
        
        if mixt == 1
            %  expll=exp(ll(qq,:));
            %  sumll=sum(expll,2);
            %  postprob=bsxfun(@rdivide,expll,sumll);
            
            % E-step: computation of posterior probabilities for untrimmed
            % units. In the context of mixture models posterior
            % probabilities will be used to estimate new component
            % probabilities of the mixtures, new centroids and new
            % covariance matrices
            
            postprobold=postprob;
            [~,postprob]=estepFS(ll);
            
            postprob(qqunassigned,:)=0;
            
            % M-step update of niini
            % niini = numerator of component probabilities
            niini=(sum(postprob))';
            
        end
        
        % M-step: parameters are updated
        % Matrix cini contains estimates of the new k centroids
        % Array sigmaini contains estimates of the new covariance matrices
        
        for j=1:k
            
            if mixt>=1
                % Matrix cini is updated using weighted means. The weights
                % are given by the posterior probabilities.
                % Note that Y is used instead of Ytri because posterior
                % probabilities for unassigned units are 0.
                cini(j,:)= sum(bsxfun(@times, Y, postprob(:,j)),1)/niini(j);
                
                if niini(j)>0
                    Ytric = bsxfun(@minus,Y,cini(j,:));
                    
                    sqweights = postprob(:,j).^(1/2);
                    
                    % Ytric = [X(:,1).*sqweights   X(:,2).*sqweights ...   X(:,end).*sqweights]
                    Ytric = bsxfun(@times, Ytric, sqweights);
                    
                    sigmaini(:,:,j) = (Ytric' * Ytric) / niini(j);
                    
                    % Eigenvalue eigenvector decomposition for group j
                    [Uj,Lambdaj] = eig(sigmaini(:,:,j));
                    
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                else
                    sigmaini(:,:,j)=ey;
                    U(:,:,j)=ey;
                    Lambda_vk(:,j)=onev1;
                end
                
            else  % This is the "crisp assignment" setting
                
                % Boolean index of units forming group j
                groupj=groupind==j;
                
                % Size of group j
                niini(j)=sum(groupj);
                
                % Group j values
                Ytrij=Ytri(groupj,:);
                
                % Means of group j
                cini(j,:)=sum(Ytrij)/niini(j);
                
                % niini=sum(Ytri(:,v+1)==j);
                if niini(j)>0
                    % Covariance of group j:
                    % sigmaini(:,:,j)=cov(Ytrij);
                    % cov would recompute the sample means; code below is
                    % more efficient
                    
                    % Important remark: DfM is a mex file with the compiled
                    % code of an efficient method to compute the following
                    % element by element operation:
                    % Ytrij = bsxfun(@minus,Ytrij,cini(j,:));
                    % The mex has been compiled for the following
                    % platforms: 32 and 64 bit MS Windows, 64 bit Linux
                    % and 64 bit MacOS. However, if you experience an error
                    % in correspondence to the DfM execution, you should
                    % comment the DfM line below and uncomment the bsxfun
                    % instruction above. In contexts where this is called
                    % many times, this solution is much more performant.
                    if callmex==true
                        DfM(Ytrij,cini(j,:),Ytrij,niini(j),v);
                    else
                        Ytrij = bsxfun(@minus,Ytrij,cini(j,:));
                    end
                    sigmaini(:,:,j) = (Ytrij' * Ytrij) / niini(j);
                    
                    % Eigenvalue eigenvector decomposition for group j
                    [Uj,Lambdaj] = eig(sigmaini(:,:,j));
                    % Store eigenvectors and eigenvalues of group j
                    U(:,:,j)=Uj;
                    Lambda_vk(:,j)=diag(Lambdaj);
                    
                else
                    sigmaini(:,:,j)=ey;
                    U(:,:,j)=ey;
                    if restrnum==1
                        Lambda_vk(:,j)=onev1;
                    else
                        Lambda_vk(j)=1;
                    end
                end
                
            end
            
        end
        
        % Lambda_vk is a v-by-k  matrix whose jth column contains the
        % unrestricted eigenvalues of cov. matrix of group j   j=1, ..., k
        % The row below is just to avoid numerical problems
        Lambda_vk(Lambda_vk<0)=0;
        if restrnum==1
            autovalues=restreigen(Lambda_vk,niini,restrfactor,tolrestreigen,userepmat);
            
            
        elseif restrnum==2
            % Restriction on the determinants
            autovalues=restrdeter(Lambda_vk,niini,restrfactor,tolrestreigen,userepmat);
        end
        
        % Covariance matrices are reconstructed keeping into account the
        % constraints of the eigenvalues
        for j=1:k
            sigmaini(:,:,j) = U(:,:,j)*diag(autovalues(:,j))* (U(:,:,j)');
            
            % Alternative code: in principle more efficient but slower
            % because diag is a built in function
            % sigmaini(:,:,j) = bsxfun(@times,U(:,:,j),autovalues(:,j)') * (U(:,:,j)');
        end
        
        
        % Calculus of the objective function (E-step)
        % oldobj=obj;
        obj = 0;
        
        if mixt>=1
            
            % Likelihood for mixture modelling
            
            %   log_lh is a h-by-k matrix where k is the number of Gaussian components of the mixture
            %   log_lh is the log of component conditional density weighted by the component
            %   probability.  The probability of j-th component is niini(j)/h
            %   log_lh(i,j) is log (Pr(point i|component j) * Prob( component j))
            
            for j=1:k
                log_lh(:,j)=  log(niini(j)/h)+logmvnpdfFS(Ytri,cini(j,:),sigmaini(:,:,j),Y0tmp(1:h,:),eyev,h,v,0);
            end
            
            % obj contains the value of the log likelihood for mixture models
            obj=estepFS(log_lh);
            
        else
            
            % Likelihood for  crisp clustering
            for j=1:k
                % disp(ni)
                if niini(j)>0
                    if equalweights
                        % we simply sum the log of the densities for the untrimmed
                        % units
                        obj=obj+ sum(logmvnpdfFS(Ytri(groupind==j,:),cini(j,:),sigmaini(:,:,j)));
                    else
                        % niini(j)*log(niini(j)/h) is the so called entropy
                        % term which allows for different group weights
                        
                        niinij=niini(j);
                        obj=obj+ niini(j)*log(niinij/h)+sum(logmvnpdfFS(Ytri(groupind==j,:),cini(j,:),sigmaini(:,:,j),Y0tmp(1:niinij,:),eyev,niinij,v,0));
                    end
                end
            end
            
        end
        
        if mixt>0
            % if mixt >0 stopping criterion is referred to postprob
            mudiff=sum(sum(abs(postprob-postprobold)))/n;
            % disp(mudiff)
        else
            % if mixt=0 stopping criterion is referred to no modiification in the classification
            mudiff=sum(abs(indold-ind)>0)/n;
            % disp(mudiff)
        end
        
        %disp(num2str(iter))
        %disp(mudiff)
        %
        % Alternative stopping criterion was based  on the relative
        % modification of the objective function.
        %                  mudiff =abs(oldobj-obj)/abs(obj);
        %                  disp(['Iteration ' num2str(iter)])
        %                  disp([oldobj-obj]/abs(obj))
        %                  disp('monit')
        
        if iter==refsteps
            noconv=noconv+1;
        end
        
    end
    
    % Store value of the objective function for iteration i
    fullsol(i)=obj;
    
    % Store the centroids and the value of the objective function
    if obj>=vopt
        % vopt = value of the objective function in correspondence of the
        % best centroids
        vopt=obj;
        % muopt = matrix containing best centroids
        muopt=cini;
        % nopt = vector containing sizes of the groups
        nopt=niini;
        % format long;
        %disp(index)
        %disp(obj)
        
        % sigmaopt
        sigmaopt=sigmaini;
        
    end
end
notconver=noconv/nselected;

%% Store quantities in out structure

% Procedure to order the non-empty components
if any(any(isnan(muopt)))
    
    % restore apropriate order of the components
    NanGroups = isnan(muopt(:,1)); % missing components
    
    % order of the components in nopt, muopt and sigmaopt
    nopt = [nopt(~NanGroups); nopt(NanGroups)];
    muopt = [muopt(~NanGroups,:); muopt(NanGroups,:)];
    sigmaopt(:,:,NanGroups) = NaN; % assign NaN on the empty clusters
    sigmaopt = cat(3, sigmaopt(:,:,~NanGroups), sigmaopt(:,:,NanGroups));
end

% With the best obtained values for the parameters, we compute the final
% assignments and parameters

% construct the  log of component conditional density weighted by the
% component probability.
% ll = log( \pi_j \phi (y_i; \; \theta_j))
% Get the likelihood for each point with each component
% ll is a n by k matrix,
% if equalweights is false
% ll(i,j) is log( (n_j/h) * f(x_i|\theta_j))
% else if  equalweights is true
% ll(i,j) is log( f(x_i|\theta_j))
% f(x_i|\theta_j) is multivariate normal with theta_j =(mu_j, \Sigma_j)
if equalweights
    for j=1:k
        if any(~isnan(muopt(j,:)))
            ll(:,j) = logmvnpdfFS(Y,muopt(j,:),sigmaopt(:,:,j),Y0tmp,eyev,n,v,0);
        else
            % avoid the computation for empty components and assign NaN
            ll(:,j) = NaN;
        end
    end
else
    for j=1:k
        if any(~isnan(muopt(j,:)))
            ll(:,j) = log(nopt(j)/h) + logmvnpdfFS(Y,muopt(j,:),sigmaopt(:,:,j),Y0tmp,eyev,n,v,0);
        else
            % avoid the computation for empty components and assign NaN
            ll(:,j) = NaN;
        end
    end
end

% matrix ll forms the input to compute both the MIXTURE and the
% CLASSIFICATION LIKELIHOOD

% postprob n x k containing posterior probabilities
% logpdf n x 1 vector containg the n contributions to the log
% likelihood of mixture models
[~,postprob,logpdf]=estepFS(ll);

% %%%%%
% In this part we select the untrimmed units
% They are those which have the n(1-alpha) largest values among the
% maxima of each row of matrix ll
% vector disc of length(n) contains the (weighted) contribution of
% each unit to the log likelihood
% idx = n x 1 vector containing the final assignments
% disc = n x 1 vector which contains the likelihood of each unit to
% the closest cluster

[disc,idx]= max(ll,[],2);

% Sort the n likelihood contributions
% qq contains the orderd (weighted) likelihood contributions
[~,qq]=sort(disc,'descend');

% %%%%%

% Find final trimmed and untrimmed units for final classification
if mixt>=1
    
    % Sort the n likelihood contributions
    % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
    [~,qqmixt]=sort(logpdf,'descend');
    
    unassignedmixt=qqmixt(h+1:n);
    assignedmixt=qqmixt(1:h);
    
    % Store in vector idx the cluster associated to the highest posterior
    % probability
    [~,idxmixt]=max(postprob,[],2);
    idxmixt(unassignedmixt)=0;
    
    postprob(unassignedmixt,:)=0;
    % Remark:
    % If there was full convergence sum(logpdf(assigned)) = vopt
else
    
    unassigned=qq(h+1:n);
    % Assign observations to clusters and assign a 0 value to trimmed ones
    idx(unassigned)=0;
    postprob(unassigned,:)=0;
end


% Compute AIC and BIC

if mixt>=1
    % Compute value of the maximized MiXTURE log likelihood
    [NlogLmixt]=estepFS(ll(assignedmixt,:));
    
    % NlogLmixt is the negative of the maximized MIXTURE LOG-LIKELIHOOD
    % Note that if there was convergence NlogL should be exactly equal to
    % -vopt
    NlogLmixt = -NlogLmixt;
end

% Note that disc(qq(1:h)) is the contribution to the CLASSIFICATION
% loglikelihood of the untrimmed units
loglik=disc(qq(1:h));

% NlogL is the negative of the CLASSIFICATION LOG-LIKELIHOOD  of the
% untrimmed units
% NlogL=-sum(max(ll(untrimmed units,[],2));
% Note that if there was convergence NlogL should be exactly equal to
% -vopt
NlogL =-sum(loglik);

% Compute INFORMATION CRITERIA
logh=log(h);

if mixt>0
    % MIXMIX = BIC which uses parameters estimated using the mixture loglikelihood
    % and the maximized mixture likelihood as goodness of fit measure (New BIC)
    MIXMIX  = 2*NlogLmixt +nParam*logh;
    
    % MIXCLA = BIC which uses the classification likelihood based on
    % parameters estimated using the mixture likelihood (New ICL)
    MIXCLA  = 2*NlogL +nParam*logh;
    
else
    % CLACLA = BIC which uses parameters estimated using the classification
    % likelihood and the maximized classification likelihood as goodness of fit
    % measure (New New)
    CLACLA  = 2*NlogL +nParam*logh;
    
end


% Store the assignments in matrix out. Unassigned units have an assignment
% equal to 0
if mixt>=1
    out.idx=idxmixt;
else
    out.idx=idx;
end

% Store value of the objective function (maximized trimmed log likelihood)
out.obj=vopt;

if out.obj==-1e+25
    warning('FSDA:tclust:NoConvergence','The result is artificially constrained due to restr.fact = 1')
end


out.fullsol=fullsol;
if mixt>0
    out.MIXMIX=MIXMIX;
    out.MIXCLA=MIXCLA;
else
    out.CLACLA=CLACLA;
end

% Store n x k matrix containing posterior probability
% of each row from each component (cluster)
out.postprob=postprob;

% Store the fraction of subsamples without convergence.
out.notconver=notconver;
% Store robust estimate of final centroids of the groups
out.muopt=muopt;

% Store robust estimate of final covariance matrix of the groups
out.sigmaopt=sigmaopt;

end
