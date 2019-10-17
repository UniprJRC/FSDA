function out = lga(X,k,varargin)
%lga performs linear grouping analysis
%
%
%<a href="matlab: docsearchFS('lga')">Link to the help function</a>
%
%  Required input arguments:
%
%       X   : input data matrix. Matrix. Input data as matrix of size
%            n-by-p
%       k   : number of clusters. Scalar. Scalar which specifies the number
%             of clusters.
%
%  Optional input arguments:
%
%
%     biter : Hyperplane number. Integer. number of different starting hyperplanes to try.
%               Example - 'biter',1
%               Data Types - double
%     niter : Number of iterations. Positive integer. Number of iterations to attempt for convergence.
%               Example - 'niter',1
%               Data Types - double
%   showall : Type of display. Logical. If true then display all the outcomes, not just the best one.
%               Example - 'showall','true'
%               Data Types - char
%    stand  : Data standardization. Logical. If true standardize the X matrix with the standard
%             deviation before fitting.
%               Example - 'stand','true'
%               Data Types - char
%    silent : Text ouptut. Logical. If true, produces no text output during processing. The
%               default value is false.
%               Example - 'silent','true'
%               Data Types - char
%    plots  : plot on the screen. Scalar. If plots=1 a plot is showed on the screen with the
%             final allocation (and if size(X,2)==2 with the lines
%             associated to the groups).
%               Example - 'plots',1
%               Data Types - double
%
%  Output:
%
%         out:   structure which contains the following fields
%
%             out.cluster = vector containing the cluster memberships.
%                out.ROSS = the Residual Orthogonal Sum of Squares for the solution.
%           out.converged = logical. True if at least one solution has converged.
%            out.nconverg = the number of converged solutions (out of biter starts).
%                   out.X = the (scaled if selected) dataset.
%                   out.hpcoeff =  best hyerplane
%              out.scaled = logical. Is the data scaled?
%                   out.k = the number of clusters to be found.
%               out.biter = the biter setting used.
%               out.niter = the niter setting used.
%               out.class = 'lga'.
%
% See also: rlga.m
%
% References:
%
% Van Aelst, S. and Wang, X. and Zamar, R. and Zhu, R. (2006), Linear
% Grouping Using Orthogonal Regression, "Computational Statistics and Data
% Analysis", Vol. 50, pp. 1287-1312.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('lga')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
%
% Examples:
%
%{
    %% lga with all default options.
    X=load('X.txt');
    out=lga(X,3);
%}
%{
    % lga with niter = 1000 and biter = 3000.
    X=load('X.txt');
    out=lga(X,4,'niter',1000,'biter',3000);
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

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:lga:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:lga:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

% Write in structure 'options' the options chosen by the user
for i=1:2:length(varargin)
    options.(varargin{i})=varargin{i+1};
end

biter=options.biter;
niter=options.niter;
silent=options.silent;

if ~silent
    disp(['LGA Algorithm k =', num2str(k), ' biter =', num2str(biter), ...
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
    %     clindex(1,:)=[191 185];
    %     clindex(2,:)=[51, 27];
    %     clindex(3,:)=[8, 35];
    
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
    outputs(:,j)=lgaiterate(hpcoef(:,:,j),X, k, d, n, niter);
end
% else
%
%     % parallel code to be introduced
%
% end


% Find the number of converged results
nconverg = sum(outputs(n+1,:));
if nconverg == 0
    warning('FSDA:lga:NoConvergence','lga failed to converge for any iteration');
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
    [~,minoutindcol]=min(outputs(n+2,:));
    outputs=outputs(:,minoutindcol);
    
    % Below we check whether there is more than one minimum
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
    
    ROSS=lgacalculateROSS(hp, X, n, d, outputs(1:n));
    
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
out.class='lga';

plots=options.plots;

if plots
    if d==2
        % spmplot(X,cluster)
        gscatter(X(:,1),X(:,2),out.cluster);
        v=axis';
        hold('on');
        for i=1:k
            a= hp(i, 3)/hp(i, 2);
            b= -hp(i, 1)/hp(i, 2);
            plot(v(1:2),a+b*v(1:2));
        end
    else
        spmplot(X,out.cluster);
    end
    
end

    function yorthreg=lgaorthreg(X)
        % Perform orthogonal regression.
        mu = mean(X);
        y = bsxfun(@minus,X, mu);
        
        %y <- scale(x, scale=FALSE)
        [~,~,V] = svd(y);
        emat=V(:,size(y,2))';
        
        yorthreg=[emat, sum(emat.*mu)];
        % emat <- svd(y)$v[,dim(y)[2]]
        % return(c(emat, emat %*% attr(y, 'scaled:center')))
    end

    function outputsj=lgaiterate(hpcoef, xsc, k, d, n, niter)
        
        % give the function the inital set of hyperplanes (in hpcoef)
        groups = lgadodist(xsc, hpcoef, d, n);
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
            groups = lgadodist(xsc, hpcoef, d, n);
            
            tabgrs=tabulate(groups);
            if min(tabgrs(:,2))< d || length(unique(groups)) < k
                iter = (niter+10);
                % if there aren't enough obs in a group
            end
            
            if (isequal(oldgroups,groups))
                converged = true;
            end
        end
        ROSS=lgacalculateROSS(hpcoef, xsc, n, d, groups);
        outputsj=[groups;converged;ROSS];
        
    end

    function indmax=lgadodist(y, coeff, d, n)
        % This function calculates the (orthogonal) Residuals for different hyerplanes,
        % and returns the closest for each observation
        
        dist = (y * (coeff(:,1:d)')-ones(n,1)*(coeff(:,d+1)') ).^2;
        [~,indmax]=max(-dist,[],2);
    end

    function ROSS=lgacalculateROSS(hpcoef, xsc, n, d, groups)
        % This function calculates the total
        % Residual Orthogonal Sum of Squares for a given grouping
        z = bsxfun(@minus,xsc *(hpcoef(:,1:d)'), hpcoef(:,d+1)');
        dist = z.^2;
        
        seq=(1:n)';
        ROSS=sum(dist((groups-1) * n + seq));
        
        % Next lines are alternative ways to extract the elements of dist
        % dist(sub2ind(size(dist), indici(:,1), indici(:,2)))
        % indici=[(1:n)' groups];
        % dist((indici(:,2)-1) * size(dist,1) + indici(:,1))
        % dist((groups-1) * n + seq)
    end


    function xbest = lgaCheckUnique(x)
        % function used above to check whether there is more than one minimum
        function zfin=CheckUniqueRand (z)
            zfin=sum(sum(z.^2))-0.5*(sum((sum(z,2)').^2)+ sum((sum(z,1)').^2));
        end
        
        d =size(x,2);
        
        index = true(d,1);
        for ii=1:(d-1)
            for jj =(ii+1):d
                y = crosstab(x(:,ii), x(:,jj));
                z = CheckUniqueRand(y);
                if (z==0)
                    index(jj) =  false;
                end
            end
        end
        if (sum(index) > 1)
            % In the incredibly unlikely situation....
            warning('FSDA:lga:MultipleSolutions','More than one unique solutions with identical ROSS -  returning first solution only')
            [~,index]=max(index);
        end
        
        xbest=x(:,index);
    end

end
%FScategory:CLUS-RobClaREG
