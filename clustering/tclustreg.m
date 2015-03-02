function out = tclustreg(X,k,factor,alpha1,alpha2,varargin)
%tclustreg performs robust linear grouping analysis
%
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help function</a>
%
%  Required input arguments:
%
%     SO FAR THERE IS JUST MATRIX X
%
%    y: A vector with n elements that contains the response variable. y can
%       be both a row of column vector.
%    X: Data matrix of explanatory variables (also called 'regressors') of
%       dimension (n x p-1). Rows of X represent observations, and columns
%       represent variables.
%       k   : scalar number of clusters
%     factor: scalar. This is the constant c controlling the scatter constrain
%
%  Optional input arguments:
%
%
%   intercept : If 1, a model with constant term will be fitted (default),
%               if 0, no constant term will be included.
%      niter  : an integer for the number of iterations to attempt for convergence
%      Kiter  : integer. Number of concentrarion steps
%      plots  : scalar. If plots=1 a plot is showed on the screen with the
%               final allocation (and if size(X,2)==2 with the lines
%               associated to the groups)
%
%  Output:
%
%  The output consists of a structure 'out' containing the following fields:
%
%   bopt =  regression parameters
%   sigmaopt are the estimated group variances
%   numopt are the number of observations in each cluster after the second trimming
%   vopt is the value of the target function
%   asig1 is the cluster assigments after first trimming ('0' means a trimmed observation...)
%   asig2 is the (-final-) cluster assigments after second trimming ('0' means a trimmed observation...)
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('rlga')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
%
%
% Examples:
%
%{
    X=load('X.txt');
    out=lga(X,3);
    out=tclustreg(X,3,5,0.1,0.1);
%}
%{
    load fishery;
    X=fishery.data;
    % some jittering is necessary because duplicated units are not treated
    % in tclustreg: this needs to be addressed
    X=X+0.000001*randn(677,2);
    out=lga(X,3);
    out=tclustreg(X,3,5,0.01,0.01,'intercept',0);
%}
%

%% Beginning of code

% Check if optimization toolbox is installed in current computer
typemin=exist('fminunc','file');

if typemin ~=2
    error('FSDA:tclustreg:MissingOptToolbox','This function requires the optimization toolbox')
end

[n,p]=size(X);

niterdef=20;

options=struct('intercept',1,'niter',niterdef,'Ksteps',10,...
    'plots',1,'output',false);

% Write in structure 'options' the options chosen by the user
for i=1:2:length(varargin);
    options.(varargin{i})=varargin{i+1};
end

% Default value of p

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:tclustreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    % Remark: the nocheck option has already been dealt by routine
    % chkinputR
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:tclustreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
end

intercept=options.intercept;

if intercept==1
    X=[ones(n,1) X];
else
    p=p-1;
end

y=X(:,end);
X=X(:,1:end-1);


Ksteps=options.Ksteps;
niter=options.niter;

notrim = floor(n*(1-alpha1));
trimm = n-notrim;

% Total trimming after second trimming
trimm2 = floor(n*(1-alpha1)*(1-alpha2));

% Initialize vector and matrices
ll = zeros(n,k);

ni =ones(1,k);
sigmaopt =ni;
bopt = zeros(k,p);
numopt = 1:k;

%Initialize the objective function through a very small value
vopt = -1e+20;


%  Randon restarts
for iter=1:niter
    
    if options.output==true
        disp(['Iteration' num2str(iter)])
    end
    
    % Initial sigmas
    sigmaini =ones(1,k);
    
    % Initial n's
    niini = floor(n*(1-alpha1)*sigmaini/k);
    
    % Initialize betas
    % bini is the matrix which will contains the estimate of \hat beta for
    % each group. Column 1 refers to first group, ... , Column k refers to
    % kth group
    bini =zeros(p,k);
    
    % Search for k*p random points avoiding group degeneracies for the x's
    
    subs=randsampleFS(n,k*p);
    % subs=[138, 17,  68,   3,  41,  32];
    
    Xb = X(subs,:);
    yb = y(subs);
    
    % Check that all submatrices of the groups are full rank
    degen = 100;
    while degen == 100
        degen=1;
        for j=1:k
            if  abs(det(  Xb((1+(j-1)*p):(j*p),1:p)  )) < 1e-50
                degen = 100;
            end
        end
        if degen==1
            break
        end
    end
    
    % Initial betas are obtained solving a linear system
    for j=1:k
        Xbj=Xb((1+(j-1)*p):(j*p),:);
        ybj=yb((1+(j-1)*p):(j*p));
        bini(:,j) =Xbj\ybj;
    end
    
    % Concentration steps
    
    indold = zeros(n,1)-1;
    for t=1:Ksteps
        
        % Discriminant functions for the assignments
        for jk=1:k
            ll(:,jk) = (niini(jk)/n)*normpdf(y-X*bini(:,jk),0,sqrt(sigmaini(jk)));
        end
        
        
        % In this part we select the untrimmed units
        % They are those which have the n(1-alpha) largest values among the
        % maxima of each row of matrix ll
        % vector disc of length(n) contains the (weighted) contribution of
        % each unit to the log likelihood
        [disc,indll]= max(ll,[],2);
        
        % Sort the n likelihood contributions
        % qq contains the largest n*(1-alpha) (weighted) likelihood contributions
        [~,qq]=sort(disc,'descend');
        
        
        % qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
        % (weighted) likelihood contributions
        qq=qq(1:n-trimm);
        
        % Ytri = n(1-alpha)-by-v matrix associated with the units
        % which have the largest n(1-alpha) likelihood contributions
        Xtri=X(qq,:);
        ytri=y(qq,:);
        indtri=indll(qq);
        % xmod = matrix with length(qq) rows which contains
        % 1st-pth column  explanatory vairiables
        % (p+1)-th column response
        % final column is
        xmod=[Xtri ytri indtri];
        
        %         xmod = as.matrix(cbind(X[qq,],ind[qq]))
        %
        % If any cluster is void or with fewer than p+1 elements we stop the concentration steps (control != 0)...
        for jj=1:k
            ni(jj) = sum(indtri==jj);
        end
        
        
        control = sum( ni <= p+2 );  % Domenico : it was p+1
        
        if control==0
            %Calculus of the new k regression parameters and the new sigmas
            xmodtemp =zeros(n,p+2);
            indxmodtemp=0;
            
            for jk=1:k
                % xmodj Data points in each group
                xmodj = xmod(xmod(:,end)==jk,:);
                
                % Perform the second trimming
                % disp(num2str(length(xmodj)))
                
                % qqs contains contains the indexes of untrimmed units for
                % group j
                if alpha2==0
                    qqs = 1:ni(jk);
                else
                    % Find the units with the smallest h distances
                    % Apply mcd on the x space (without the intercept if
                    % present)
                    if intercept
                        RAW = mcd(xmodj(:,2:p),'bdp',alpha2,'msg',0);
                    else
                        RAW = mcd(xmodj(:,1:p),'bdp',alpha2,'msg',0);
                    end
                    [~,indmdsor]=sort(RAW.md);
                    qqs=indmdsor(1:floor(ni(jk)*(1-alpha2)));
                end
                
                % Update betas through ordinary least squares regression
                xxx = xmodj(qqs,1:p);
                yyy = xmodj(qqs,p+1);
                ni(jk) = length(yyy);
                breg = xxx\yyy;
                bini(:,jk) = breg;
                
                % now find residuals
                residuals=yyy-xxx*breg;
                % Update sigmas through the mean square residuals
                sigmaini(jk) =sum(residuals.^2)/ni(jk);
                
                % Update weights
                niini(jk) = ni(jk);
                xmodtemp((indxmodtemp+1):(indxmodtemp+ni(jk)),:)=xmodj(qqs,:);
                indxmodtemp=indxmodtemp+ni(jk);
            end
            
            % New xmod
            xmod = xmodtemp(1:indxmodtemp,:);
            
            % If the scatters do not satisfy the restriction then a
            % quadratic programming problem is solved
            sigmaini = (quadi(sigmaini.^(-1), factor)).^(-1);
        end
        
        
        
        % This criterium serves to stop if two consecutive concentration steps has the same result
        if indll == indold
            break
        else
            indold = indll;
        end
    end
    
    
    % After the concentration steps, we compute the value of the target function
    obj =0;
    for jk=1:k
        if control==0,
            yj=xmod(xmod(:,end)==jk,end-1);
            Xj=xmod(xmod(:,end)==jk,1:end-2);
            obj=obj+ niini(jk)*log(niini(jk)/trimm2)+sum(log(normpdf(yj-Xj*bini(:,jk),0,sqrt(sigmaini(jk)))));
        else
            obj=obj-10^10;
        end
    end
    
    % Change the 'optimal' target value and 'optimal' parameters if a increase in the target value is achieved
    if (obj >= vopt)
        vopt = obj;
        bopt = bini;
        numopt = niini;
        sigmaopt = sigmaini;
    end
    disp(['iter ' num2str(iter)]);
end

%Last part of the program prepares some graphs and numerical outputs for the program
% Assignment vectors: asig.1 are the clusters after the first trimming and asig.2 after the second
asig1 =zeros(n,1);
asig2 =asig1;


for jk=1:k
    ll(:,jk) = (numopt(jk)/n)*normpdf(y-X*bopt(:,jk),0,sqrt(sigmaopt(jk)));
end

[dist,indll]= max(ll,[],2);

% Sort the n likelihood contributions
% qq contains the largest n*(1-alpha) (weighted) likelihood contributions
[val,qq]=sort(dist,'descend');
val=val(n-trimm);



for jk=1:k
    asig1((indll==jk) & (dist>=val)) = jk;
end

plots=options.plots;

% A graph summarizing the results is given if dimension is equal p = 2
if (p<=2 && plots)
    plot(X(dist<val,end),y(dist<val),'o','color','r')
    hold('on')
end

% qq = vector of size h which contains the indexes associated with the largest n(1-alpha)
% (weighted) likelihood contributions
qq=qq(1:n-trimm);
xmod=[X(qq,:) y(qq) indll(qq)];

for jk=1:k
    
    booljk=xmod(:,end)==jk;
    qqk = qq(booljk);
    ni(jk) = sum(booljk);
    xmodjk = xmod(booljk,:);
    
    
    if alpha2==0
        qqs = 1:ni(jk);
    else
        
        if intercept
            RAW = mcd(xmodjk(:,2:p),'bdp',alpha2,'msg',0);
        else
            RAW = mcd(xmodjk(:,1:p),'bdp',alpha2,'msg',0);
        end
        [~,indmdsor]=sort(RAW.md);
        qqs=indmdsor(1:floor(ni(jk)*(1-alpha2)));
    end
    xxx = xmodjk(qqs,1:end-2);
    yyy = xmodjk(qqs,end-1);
    
    qqsn=setdiff(1:ni(jk),qqs);
    xxx0 = xmodjk(qqsn,1:end-2);
    yyy0 = xmodjk(qqsn,end-1);
    
    
    if (p<=2 && plots)
        plot(xxx(:,end),yyy,'.w');
        % points of each component (pch=k+2)
        text(xxx(:,end),yyy,num2str(jk*ones(length(yyy),1)) , ...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle');
        % second level trimming points
        plot(xxx0(:,end),yyy0,'*','color','c')
    end
    
    qqf = qqk(qqs);
    asig2(qqf) = jk;
    
    if (p<=2 && plots)
        reg=xxx\yyy;
        %         v=axis';
        %         plot(v(1:2),reg(1)+reg(2)*v(1:2))
        v = [min(X(:,end)) max(X(:,end))];
        if intercept==1
            plot(v,reg(1)+reg(2)*v)
        elseif intercept==0
            plot(v,reg*v)
        end
    end
end


%   Return the values for the function
%   bopt are the regression parameters
%   sigmaopt are the estimated group variances
%   numopt are the number of observations in each cluster after the second trimming
%   vopt is the value of the target function
%   asig1 is the cluster assigments after first trimming ('0' means a trimmed observation...)
%   asig2 is the (-final-) cluster assigments after second trimming ('0' means a trimmed observation...)
out=struct;
out.bopt=bopt;
out.sigmaopt=sigmaopt;
out.numopt=numopt;
out.vopt=vopt;
out.asig1=asig1;
out.asig2=asig2;



% quadi is the subfunction which prepares the quantities to call the matlab quadratic
% programming routine quadprog,
    function   gnew=quadi(gg,factor)
        if size(gg,1)>1
            gg=gg';
        end
        
        % gnew will the new scatters
        gnew = gg;
        
        if (length(gg)>1)
            %Sort scatters
            [ggsor,ggsorind] = sort(gg);
            
            % g(1) = smallest sigma
            % ...
            % g(end) = largest sigma
            g = ggsor;
            
            maximun = 10^5;
            
            % Constant "c" defining the scatter constraint
            factor = factor+0.0001;
            
            % nscat is the number of scatter parameters
            nscat = length(g);
            
            Amat =zeros(nscat,nscat);
            % rr = 1:nscat;
            
            for ii =1:(nscat-1)
                Amat(ii,ii) = -1;
                Amat(ii,ii+1) =1;
            end
            
            % Definition of the quadratic problem
            Amat(nscat,1) = factor;
            Amat(nscat,nscat) = -1;
            Vmat = diag([ones(nscat,1);zeros(nscat,1)]);
            dvec = - [g,zeros(1,nscat)];
            bvec = zeros(1,nscat);
            uvecmax = maximun+zeros(1,2*nscat);
            uvecmin = zeros(1,2*nscat);
            
            Amat = [Amat,-1*eye(nscat)];
            
            % Solve this quadratic problem
            % a = quadprog(Vmat,dvec,[],[],Amat,bvec',uvecmin,uvecmax,g,'Algorithm','interior-point-convex');
            
            option = optimoptions('quadprog','algorithm','interior-point-convex','Display','off');
            a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],option);
            
            %a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],'algorithm','interior-point-convex','Display','iter');
            %a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],'algorithm','active-set');
            
            gnew =a(1:nscat);
            
            %Original order
            gnew(ggsorind) = gnew;
            
        end
    end
end