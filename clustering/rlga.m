function out = rlga(X,k,alpha,varargin)
%rlga performs robust linear grouping analysis
%
%<a href="matlab: docsearchFS('rlga')">Link to the help function</a>
%
%  Required input arguments:
%
%       X   : input data matrix. Matrix. Input data as matrix of size
%             n-by-p
%       k   : number of clusters. Scalar. Scalar which specifies the number
%             of clusters.
%    alpha  : a numeric value between 0 and 0.5. Scalar. For the robust estimate
%             of RLGA, specifying the percentage of points to be trimmed.
%             alpha must be a number in the interval [0 0.5].
%
%  Optional input arguments:
%
%
%     biter : number of different starting hyperplanes. Scalar. An integer
%             for the number of different starting hyperplanes to try.
%               Example - 'biter',1
%               Data Types - double
%
%     niter : Number of iterations. Scalar. An integer for the number of
%             iterations to attempt for convergence.
%               Example - 'niter',1
%               Data Types - double
%
%   showall :  Level of display. Logical. If true then display all the
%               outcomes, not just the best one.
%               Example - 'showall','true'
%               Data Types - char
%
%    stand  : Data standardization. Boolean. If true the X
%               matrix is standardized using the standard deviation before
%               fitting. Logical.
%               Example - 'stand','true'
%               Data Types - char
%
%    silent : Text output. Logical. If true, produces no text output during
%               processing. The default value is false.
%               Example - 'silent','true'
%               Data Types - char
%
%    plots  : plot on the screen. Scalar. If plots=1 a plot is showed on the
%               screen with the final allocation (and if size(X,2)==2 with the
%               lines associated to the groups).
%               Example - 'plots',1
%               Data Types - double
%
%  Output:
%
%         out:   structure which contains the following fields
%
%         out.cluster   = vector containing the cluster memberships.
%         out.ROSS      = the Residual Orthogonal Sum of Squares for the solution.
%         out.converged = logical. True if at least one solution has converged.
%         out.nconverg  = the number of converged solutions (out of biter starts).
%         out.hpcoeff   = best hyerplane
%         out.X	        = the (scaled if selected) dataset.
%         out.scaled    = logical. Is the data set scaled?
%         out.k         = the number of clusters to be found.
%         out.alpha     = the trimming percentage.
%         out.biter     = the biter setting used.
%         out.niter	    = the niter setting used.
%         out.class     = 'rlga'.
%
% See also: lga.m
%
% References:
%
% Garcia-Escudero, L. A., Gordaliza, A., San Martin, R., Van Aelst, S. and Zamar, R. (2009), 
% Robust linear clustering, "Journal of the Royal Statistical Society:
% Series B", Vol. 71, pp. 301-318.
% [doi: 10.1111/j.1467-9868.2008.00682.x]
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('rlga')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit


%
% Examples:
%
%{
    % rlga with all default options.
    rng(123); % this leads to ROSS = 5.1298;
    X=load('X.txt');
    out=rlga(X,3,0.05);
%}
%{
    % rlga with niter = 500 and biter = 1000.
    X=load('X.txt');
    out=rlga(X,4,0.05,'niter',500,'biter',1000);
%}
%{
    %% Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids = 0.01. Use all default options.
    rng(372,'twister');
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=500;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);

    % run rlga
    out=rlga([y,X(:,2:end)],k,0.01);
%}
%

%% Beginning of code

[n,d]=size(X);

n1= ceil(n/k);
p=0.95;
biterdef=ceil(log(1-p)/log(1-nchoosekFS(n1,d)^k/nchoosekFS(n1*k,k*d)));

niterdef=10;

options=struct('biter',biterdef,'niter',niterdef,'showall',false,...
    'stand',true,'silent',false,'plots',1);

% Write in structure 'options' the options chosen by the user
for i=1:2:length(varargin)
    options.(varargin{i})=varargin{i+1};
end

% Default value of p

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
       
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:rlga:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:rlga:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

if alpha>0.5
    error('FSDA:rlga:WrongAlpha','Error:: alpha must be in [0, 0.5]');
end

biter=options.biter;
niter=options.niter;
silent=options.silent;

if ~silent
    disp(['RLGA Algorithm k =', num2str(k), ' biter =', num2str(biter), ...
        ' niter =', num2str(niter), ''])
end

stand=options.stand;

if stand
    sigma = std(X);
    X = bsxfun(@rdivide, X, sigma);
end

% hpcoef is a a 3D array of size
% number of cluster-times-number of variables+1-times-biter
% each row is a hyperplane
hpcoef=zeros(k,d+1,biter);

for j=1:biter
    % Choose starting clusters
    clindex=reshape(randsample(1:n,k*d),k,d);
    for i=1:k
        hpcoef(i,:,j)=lgaorthreg(X(clindex(i,:),:));
    end
end

% nnode = an integer of many CPUS to use for parallel processing. Defaults
% to NULL i.e. no parallel processing.
% nnode=options.nnode;
% if isempty(nnode)
% not parallelrows 1:n contain the information about group assignment
% row n+1 contains information about convergence (1) or non convergence (0)
% row n+2 contains information about ROSS
outputs=zeros(n+2,biter);
for j=1:biter
    outputs(:,j) = rlgaiterate(hpcoef(:,:,j),X, k, d, n, niter,alpha);
end
% else
%     % parallel
%
% end

% Find the number of converged results
nconverg = sum(outputs(n+1,:));
if nconverg == 0
    warning('FSDA:rlga:NoConvergence','rlga failed to converge for any iteration')
end

showall=options.showall;

if ~showall
    % remove any columns of matrix outputs which contains NAs
    outputs = outputs(:, sum(isnan(outputs))==0);
    
    
    if (nconverg ~= 0)
        outputs = outputs(:,outputs(n+1,:)==1);
    end
    
    % using the instruction below we take the column of outputs which is
    % associated to the first minimum
    [~,minoutindcol] = min(outputs(n+2,:));
    outputs = outputs(:,minoutindcol);
    
    % Instructions below are to check whether there is more than one minimum
    %
    %{
        minout=min(outputs(n+2,:));
        minoutindcol=outputs(n+2,:)==minout;
        outputs = outputs(:,minoutindcol);

        if size(outputs,2) > 1
            outputs = lgaCheckUnique(outputs);
        end
    %}
    
end

if (~silent)
    disp('Finished.')
end

if showall
    disp('Returning all outputs')
    ROSS = outputs(n+2,:);
    hp =nan;
else
    % Fit the best hyerplane(s) with ROSS
    hp =nan(k,d+1);
    for i=1:k
        hp(i,:)=lgaorthreg(X(outputs(1:n)== i,:));
    end
    
    ROSS=rlgacalculateROSS(hp, X, d, outputs(1:n));
end

out.cluster=outputs(1:n,:);
out.ROSS=ROSS;
out.converged=outputs(n+1,:);
out.nconverg=nconverg;
out.X=X;
out.hpcoeff=hp;
out.biter=biter;
out.niter=niter;
out.scaled=stand;
out.k=k;
out.alpha=alpha;
out.class='rlga';

if options.plots
    % the plot requires the data X, the hyperplan information hp, and the
    % output structure of lga or rlga
    linear_grouping_plot(X,hp,out);
end

%% inner funtions

%% lgaorthreg performs orthogonal regression
    function yorthreg = lgaorthreg(X)
        
        mu = mean(X);
        y = bsxfun(@minus,X, mu);
        
        [~,~,V] = svd(y);
        emat=V(:,size(y,2))';
        
        yorthreg=[emat, sum(emat.*mu)];
        
    end

%% Do "concentration steps"

    function outputsj=rlgaiterate(hpcoef, xsc, k, d, n, niter,alpha)
        
        % give the function the inital set of hyperplanes (in hpcoef)
        groups = rlgadodist(xsc, hpcoef, d, n,alpha);
        iter = 0;
        tabgrs=tabulate(groups);
        if min(tabgrs(:,2)< d) || length(unique(groups)) < k
            iter = niter+10;
        end
        
        converged = false;
        while (~converged && iter < niter)
            iter = iter+1;
            oldgroups = groups;
            for ii=1:k
                hpcoef(ii,:) = lgaorthreg(xsc(groups==ii,:));
            end
            groups = rlgadodist(xsc, hpcoef, d, n, alpha);
            
            tabgrs=tabulate(groups);
            if min(tabgrs(:,2))< d || length(unique(groups)) < k
                iter = (niter+10);
                % if there aren't enough obs in a group
            end
            
            if (isequal(oldgroups,groups))
                converged = true;
            end
        end
        iROSS     = rlgacalculateROSS(hpcoef, xsc, d, groups);
        outputsj = [groups;converged;iROSS];
        
    end

%% rlgadodist calculates the (orthogonal) Residuals for different hyerplanes

    function indmin=rlgadodist(y, coeff, d, n,alpha)
        % Find which hyperplane each observation is closest to. Then takes the
        % smallest alpha of them (setting the rest to zero)
        dist = (y * (coeff(:,1:d)')-ones(n,1)*(coeff(:,d+1)') ).^2;
        [distmin,indmin]=min(dist,[],2);
        indmin( distmin > quantile(distmin, 1 - alpha) ) = 0;
    end

%% rlgacalculateROSS calculates the total Residual Orthogonal Sum of Squares for a given grouping

    function outROSS = rlgacalculateROSS(hpcoef, xsc, d, groups)
        
        z    = bsxfun(@minus,xsc *(hpcoef(:,1:d)'), hpcoef(:,d+1)');
        dist = z.^2;
        
        seldis = groups>0;
        dist   = dist(seldis,:);
        groups = groups(seldis);
        
        % extract elements of dist
        seq1 = (1:length(groups))';
        
        outROSS=sum(dist((groups-1) * length(groups) + seq1));
    end

%% lgaCheckUnique check if there is more than one minimum

% used for debugging

    function xbest = lgaCheckUnique(x) %#ok<DEFNU>
        function zfin=CheckUniqueRand (z)
            zfin=sum(sum(z.^2))-0.5*(sum((sum(z,2)').^2)+ sum((sum(z,1)').^2));
        end
        
        dd =size(x,2);
        
        index = true(dd,1);
        for ii=1:(dd-1)
            for jj =(ii+1):dd
                y = crosstab(x(:,ii), x(:,jj));
                z = CheckUniqueRand(y);
                if (z==0)
                    index(jj) =  false;
                end
            end
        end
        if (sum(index) > 1)
            % In the incredibly unlikely situation....
            warning('FSDA:rlga:DuplicateSol','More than one unique solutions with identical ROSS -  returning first solution only')
            [~,index]=max(index);
            
        end
        
        xbest=x(:,index);
    end

end
%FScategory:CLUS-RobClaREG
