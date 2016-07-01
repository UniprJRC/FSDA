function [out] = tclustreg(y,X,k,restrfact,alpha1,alpha2,varargin)
%tclustreg performs robust clustering in regression 
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help function</a>
%
%  Required input arguments:
%
%         y : Response variable. Vector. 
%             A vector with n elements that contains the response variable.
%             y can be either a row or a column vector.
%             Data Types - single|double
%
%         X : Explanatory variables (also called 'regressors'). Matrix.
%             Data matrix of dimension $(n \times p-1)$. Rows of X represent
%             observations, and columns represent variables. Missing values
%             (NaN's) and infinite values (Inf's) are allowed, since
%             observations (rows) with missing or infinite values will
%             automatically be excluded from the computations.
%             Data Types - single|double
%
%         k : Number of clusters. Scalar. 
%             This is a guess on the number of data groups.
%             Data Types - single|double
%
% restrfact : Scatter constraint. Scalar. 
%            This is a constant c controlling the differences among
%            group scatters. The value 1 is the strongest restriction.
%            Data Types - single|double
%
%   alpha1 : Trimming level. Scalar. 
%            alpha1 is a value between 0 and 0.5 or an  integer specifying
%            the number of observations which have to be trimmed. If
%            alpha=0 there is no trimming. More in detail, if 0<alpha1<1
%            clustering is based on h=fix(n*(1-alpha1)) observations.
%            Else if alpha1 is an integer greater than 1 clustering is
%            based on h=n-floor(alpha1).
%            Data Types - single|double
%
%   alpha2 : Second-level trimming. Scalar. 
%            alpha2 is a value between 0 and 0.5, usually smaller than
%            alpha1. If alpha2=0 there is no second-level trimming.
%            Data Types - single|double
%
%
%  Optional input arguments:
%
%intercept : Indicator for constant term. Scalar. If 1, a model with
%            constant term will be fitted (default), if 0, no constant
%            term will be included.
%            Example - 'intercept',1 
%            Data Types - double
%
%    niter : Number of random starts. Scalar. An integer for the number 
%            of iterations to attempt for convergence.
%            Example - niter = 20 
%            Data Types - double
%
%    Kiter : Number of concentrarion steps. Scalar. An integer for the 
%            number of concentration steps. 
%            Example - niter = 10 
%            Data Types - double
%
%    plots : Plot on the screen. Scalar. A flag to control the
%            generation of the plots.
%            If plots=1 a plot is showed on the screen with the
%            final allocation (and if size(X,2)==2 with the lines
%            associated to the groups)
%            Example - 'plots',1 
%            Data Types - double
%
%  Output:
%
%  out :  structure containing the following fields
%
%   out.bopt        = $p-1 \times k$ matrix containing the regression
%                     parameters.
%   out.sigmaopt    = $k$ row vector containing the estimated group
%                     variances.
%   out.numopt      = $k$ column vector containing the number of
%                     observations in each cluster after the second
%                     trimming.
%   out.vopt        = Scalar. The value of the target function.
%   out.asig1       = $n$ vector containing the cluster assigments after 
%                     first trimming ('0' means a trimmed observation).
%   out.asig2       = $n$ vector containing the final cluster assigments 
%                     after second trimming ('0' means a trimmed
%                     observation).
%
%
% See also: tclust, tkmeans, estepFS
%
% References:
%
% Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar, A.
% (2008), "A General Trimming Approach to Robust Cluster Analysis". Annals
% of Statistics, Vol.36, 1324-1345. Technical Report available at
% www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustreg')">Link to the help page for this function</a>
% Last modified 11-06-2016
%
%
%
% Examples:
%
%{
    X=load('X.txt');
    out=lga(X,3);

    y1=X(:,end);
    X1=X(:,1:end-1);

    out=tclustreg(y1,X1,3,5,0.1,0.1);
%}
%{
    load fishery;
    X=fishery.data;
    % some jittering is necessary because duplicated units are not treated
    % in tclustreg: this needs to be addressed
    X=X+0.000001*randn(677,2);
    out=lga(X,3);
    out=rlga(X,3,0.5);

    y1=X(:,end);
    X1=X(:,1:end-1);

    out=tclustreg(y1,X1,3,5,0.01,0.01,'intercept',0);
%}

%{
    % Generate mixture of regression using MixSimReg, with an average
    % overlapping at centroids =0,001. Use all default options.
    p=3;
    k=2;
    Q=MixSimreg(k,p,'BarOmega',0.001);
    n=400;
    [y,X,id]=simdatasetreg(n,Q.Pi,Q.Beta,Q.S,Q.Xdistrib);
    spmplot([y X(:,2:end)],id);
    out=tclustreg(y,X,2,50,0.01,0.01,'intercept',1);

%}


%% Input parameters checking

nnargin=nargin;
vvarargin=varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

%% Few internal parameters and basic checks

% default number of random starts
niterdef = 20;

% default number of concentration starts
Ksteps  = 10;

%Initialize the objective function (optimized during the random starts)
%through a very small value
vopt = -1e+20;

% this is just for rotating colors in the plots
clrdef='bkmgyrcbkmgyrcbkmgyrcbkmgyrcbkmgyrcbkmgyrcbkmgyrc';

% Check if optimization toolbox is installed in current computer
typemin=exist('fminunc','file');

if typemin ~=2
    error('FSDA:tclustreg:MissingOptToolbox','This function requires the optimization toolbox')
end

%% User options

options=struct('intercept',1,'niter',niterdef,'Ksteps',Ksteps,...
    'plots',1,'output',false);

if nargin > 6
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)  
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustreg:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if all the specified optional arguments were present in
        % structure options. Remark: the nocheck option has already been dealt
        % by routine chkinputR.
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:tclustreg:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end

end
% Graph summarizing the results, yes/no; For p=2 the plot is specific for
% tclustreg. Otherwise spmplot is used
plots       = options.plots;

% Intercept, yes/no
intercept   = options.intercept;

% Number of variables without considering the constant term. It is used for
% deciding the type of plot.
if intercept == 1 
    v = p-1;
else
    v = p;
end

% lines below are obsolete: replaced by chkinputR
% if intercept==1
%     X=[ones(n,1) X];
% else
%     p=p-1;
% end
% 
% y=X(:,end);
% X=X(:,1:end-1);

% Concentration steps
Ksteps      = options.Ksteps;

% Random starts
niter       = options.niter;

% First level trimming
notrim      = floor(n*(1-alpha1));
trimm       = n-notrim;

% Total trimming after second trimming
trimm2 = floor(n*(1-alpha1)*(1-alpha2));

%% Initialize vector and matrices

ll          = zeros(n,k);
ni          = ones(1,k);
sigmaopt    = ni;
bopt        = zeros(k,p);
numopt      = 1:k;

%%  Random starts

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
    nameYY =zeros(p,k);
    
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
            if  abs(det(  Xb((1+(j-1)*p):(j*p),1:p) )) < eps
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
        nameYY(:,j) =Xbj\ybj;
    end
    
    %% Concentration steps
    
    indold = zeros(n,1)-1;
    for t=1:Ksteps
        
        % Discriminant functions for the assignments
        for jk=1:k
            ll(:,jk) = (niini(jk)/n)*normpdf(y-X*nameYY(:,jk),0,sqrt(sigmaini(jk)));
        end
        
        
        % In this part we select the untrimmed units. They are those which
        % have the n(1-alpha) largest values among the maxima of each row
        % of matrix ll.
        % vector disc of length(n) contains the (weighted) contribution of
        % each unit to the log likelihood.
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
        % 1st-pth column  explanatory variables
        % (p+1)-th column response
        xmod=[Xtri , ytri , indtri];
        
        % If any cluster is void or with fewer than p+1 elements we stop the concentration steps (control != 0)...
        for jj=1:k
            ni(jj) = sum(indtri==jj);
        end
        
        control = sum( ni <= p+2 );  % Domenico : was p+1
        
        if control==0
            
            % new k regression parameters and new sigmas
            xmodtemp    = zeros(n,p+2);
            indxmodtemp = 0;
            
            for jk=1:k
                % xmodj Data points in each group
                xmodj = xmod(xmod(:,end)==jk,:);
                
                % Perform the second trimming
                
                % qqs contains contains the indexes of untrimmed units for
                % group j
                if alpha2==0
                    qqs = 1:ni(jk);
                else
                    % Find the units with the smallest h distances.
                    % Apply mcd on the x space (without the intercept if
                    % present).
                    % REMARK: This is by far the computationally most 
                    % expensive instruction of tclustreg. More precisely,
                    % the dominant expensive function inside mcd is IRWLSmcd.
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
                nameYY(:,jk) = breg;
                
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
            sigmaini = (quadi(sigmaini.^(-1), restrfact)).^(-1);
        end
           
        % Stop if two consecutive concentration steps have the same result
        if indll == indold
            break
        else
            indold = indll;
        end
    end
    
    %% Concentration steps concluded
    
    % Now compute the value of the target function
    obj =0;
    for jk=1:k
        if control==0,
            yj=xmod(xmod(:,end)==jk,end-1);
            Xj=xmod(xmod(:,end)==jk,1:end-2);
            obj = obj + niini(jk)*log(niini(jk)/trimm2) +...
                  sum(log(normpdf(yj-Xj*nameYY(:,jk),0,sqrt(sigmaini(jk)))));
        else
            obj=obj-10^10;
        end
    end
    
    % Change the 'optimal' target value and 'optimal' parameters if an
    % increase in the target value is achieved
    if (obj >= vopt)
        vopt = obj;
        bopt = nameYY;
        numopt = niini;
        sigmaopt = sigmaini;
    end
    disp(['iter ' num2str(iter)]);
end

%% Prepares the output structure and some variables for the plots

% Assignment vectors:
% - asig.1 will contain the clusters after the first trimming
% - asig.2 will contain the clusters after after the second trimming
asig1 = zeros(n,1);
asig2 = asig1;

% log-likelihoods for each unit and group
for jk=1:k
    ll(:,jk) = (numopt(jk)/n)*normpdf(y-X*bopt(:,jk),0,sqrt(sigmaopt(jk)));
end

% indll: for each unit, group with best log-likelihood
[dist,indll] = max(ll,[],2);

% Sort the n likelihood contributions;
[val,qq] = sort(dist,'descend');

% qq is updated to be a vector of size h which contains the indexes
% associated with the largest n(1-alpha) (weighted) likelihood
% contributions
qq  = qq(1:n-trimm);

% boolean vectors indicating the good and outlying units
val = val(n-trimm);
b_good = (dist>=val);
b_outl = (dist <val);

% asig1: grouping variable for good units
for jk=1:k
    asig1((indll==jk) & b_good) = jk;
end

xmod = [X(qq,:) y(qq) indll(qq)];

xxx0_all = [];
yyy0_all = [];

% go over the groups
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
        [~,indmdsor] = sort(RAW.md);
        qqs = indmdsor(1:floor(ni(jk)*(1-alpha2)));
    end
    
    % good units of the current group
    xxx = xmodjk(qqs,1:end-2);
    yyy = xmodjk(qqs,end-1);
    
    % second level trimming units of the current group
    qqsn=setdiff(1:ni(jk),qqs);
    xxx0 = xmodjk(qqsn,1:end-2);
    yyy0 = xmodjk(qqsn,end-1);
    
    % collect all second level trimming units in a same group, for plotting
    xxx0_all = [xxx0_all ; xxx0(:,end)]; %#ok<AGROW>
    yyy0_all = [yyy0_all ; yyy0];        %#ok<AGROW>
    
    % plot good units allocated to the current group
    if (plots && v < 2) 
        
        % initialize figure
        if jk == 1
            fh = figure('Name','TclustReg plot','NumberTitle','off','Visible','on');
            ah=gca(fh);
            hold on;
            xlabel('X');
            ylabel('y');
            title('TclustReg clustering','Fontsize',14);
            %print(fh,[graphs 'PCvariance.png'],'-dpng');
            %close(fh);
        end
        
        group_label = ['Group ' num2str(jk)];
        plot(xxx(:,end),yyy,'.w','DisplayName',group_label);
        % units of each component (pch=k+2)
        text(xxx(:,end),yyy,num2str(jk*ones(length(yyy),1)),...
            'DisplayName',group_label , ...
            'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'Color',clrdef(jk));
        % % second level trimming points
        % % non worth having them by group
        % plot(xxx0(:,end),yyy0,'*','color','c',...
        % 'DisplayName','Level-2 trim');
    end
    
    qqf = qqk(qqs);
    asig2(qqf) = jk;
    
    % plot regression lines
    if (plots && v < 2)
        reg=xxx\yyy;
        vv = [min(X(:,end)) max(X(:,end))];
        if intercept==1
            plot(vv,reg(1)+reg(2)*vv,...
                'DisplayName',['fit of group '  num2str(jk)],...
                'Color',clrdef(jk));
        elseif intercept==0
            plot(vv,reg*vv,...
                'DisplayName',['fit of group '  num2str(jk)],...
                'Color',clrdef(jk));
        end
    end
    
end

if plots
    
    if v < 2
    % Plot the group of outliers
    plot(X(b_outl,end),y(b_outl),'o','color','r',...
        'DisplayName','Trimmed units');
    
    % second level trimming points
    plot(xxx0_all,yyy0_all,'*','color','c',...
        'DisplayName','L2 trimmed units');
    
    % position the legends and make them clickable
    lh=legend('show');
    %set(lh,'FontSize',14);
    axis('manual');
    legstr = get(lh,'String');
    clickableMultiLegend(legstr,'FontSize',14,'Location','northwest');
    %[hleg, hobj, hout, mout] = clickableMultiLegend(legstr,'FontSize',14,'Location','northwest');
    %Unfortunately custom markers for line objects are not possible in MATLAB
    %set(hobj(10),'Marker','1','Color','k');
    end

else
    % in this case p > 2 and a standard spmplot is used
    
    if intercept
        YY = [X(:,2:end),y];
    else
        YY = [X,y];
    end
    
    % axis labels
    nameYY = cellstr([repmat('X',size(YY,2)-1,1) , num2str((1:size(YY,2)-1)')]);
    nameYY = [nameYY ; 'y'];
    nameYY = nameYY';
    plo=struct;
    plo.nameY=nameYY;
    
    % group names in the legend
    group = cell(size(asig2,1),1);
    group(asig2==0) = {'Trimmed units'};
    for iii = 1:k
        group(asig2==iii) = {['Group ' num2str(iii)]};
    end
    
    % scatterplot
    hout = spmplot(YY,group,plo,'hist');
    
    %group_l = cellstr([repmat('Group',k,1) , num2str((1:k)')]);
    %group_l = ['Trimmed units' ; group];
    %[hleg, hobj, hout, mout] = legend((hout(1,end,:)));
end

%%  Set the output structure

out             = struct;
out.bopt        = bopt;
out.sigmaopt    = sigmaopt;
out.numopt      = numopt;
out.vopt        = vopt;
out.asig1       = asig1;
out.asig2       = asig2;

%   bopt        are the regression parameters
%   sigmaopt    are the estimated group variances
%   numopt      are the number of observations in each cluster after the
%               second trimming
%   vopt        is the value of the target function
%   asig1       is the cluster assigments after first trimming ('0' means a
%               trimmed observation)
%   asig2       is the (-final-) cluster assigments after second trimming
%               ('0' means a trimmed observation)

%% subfunction quadi
%  prepares the quantities to call the matlab quadratic programming routine
%  quadprog
    function gnew = quadi(gg,factor)
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
            
            % FSDATOAPP:tclustreg:DF
            % Remark: for compatibilty with old version of MATLAB we use
            % intruction optimset. However recent versions of Matlab accept
            % function optimoptions as follows
            % option = optimoptions('quadprog','algorithm','interior-point-convex','Display','off');
            option = optimset('OutputFcn','quadprog','algorithm','interior-point-convex','Display','off');
            
            a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],option);
            %a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],'algorithm','interior-point-convex','Display','iter');
            %a = quadprog(Vmat,dvec,[],[],Amat,bvec,uvecmin,uvecmax,[],'algorithm','active-set');
            
            gnew =a(1:nscat);
            
            %Original order
            gnew(ggsorind) = gnew;
            
        end
    end
end
%FScategory:CLUS-RobClaREG
