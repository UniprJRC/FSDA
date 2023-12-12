function [GAMc]  = restrshapeExact(k, v, shw, shb, verLess2016b)
%restrshapeExact computes constrained Gamma (shape) matrix with exact constraints.
%
%
%<a href="matlab: docsearchFS('restrshapeExact')">Link to the help function</a>
%
%
% The purpose is to generate a constrained shape matrix in such a way that
% the maximum ratio between the largest/smallest element of each column is
% equal to shw and the maximum ratio between the largest/smallest element
% of each row is equal to shb. Finally the product of the element of each
% column is equal to 1.
%
% Required input arguments:
%
%            k: number of groups (components). Scalar.
%               Desired number of groups.
%               Data Types - double
%            v: number of dimensions (variables). Scalar.
%               Desired number of variables.
%               Data Types - double
%     shw  : within groups shape constraint. Scalar greater or equal 1.
%           Constraint to impose inside each group. For example, if shw is
%           3 the maximum ratio between the largest/smallest of each column
%           of output matrix GAMc will be equal to 3.
%               Data Types - double
%     shb  : across groups shape constraint. Scalar greater or equal 1.
%           Constraint to impose among groups. For example, if shb is 5 the
%           maximum ratio between the largest/smallest of each row of
%           matrix GAMc will be equal to 5. Note that shb must be smaller
%           or equal than shw.
%               Data Types - double
%
%
%  Optional input arguments:
%
%
% verLess2016b : is a boolean which is true if current version is less than
%                   2016b. If this argument is omitted routine
%                   automatically checks the MATLAB version.
%               Example - true
%               Data Types - double
%
% Output:
%
%     GAMc : constrained shape matrix. Matrix of size v-by-k containing in
%           column j the elements on the main diagonal of shape matrix
%           $\Gamma_j$, $j=1, 2, \ldots, k$. The elements of GAMc satisfy
%           the following constraints:
%           1) the product of the elements of each column is equal to 1;.
%           2) The maximum ratio among the largest element divided by the smallest
%           3) element of each column is equal to shw;
%           4) The maximum ratio among the largest element divided by the smallest
%           element of each row is equal to shb. The tolerance of the
%           procedure is 0.01.
%           5) All the elements of GAMc are strictly positive.
%
%
% See also tclust, restrshapeGPCM, tclustreg
%
% References:
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('restrshapeExact')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples:
%
%{
    %% restrshapeExact with all the default options.
    % two groups and 8 variables
    k=2;
    v=8;
    shb=20;
    shw=100;
    Sh=restrshapeExact(k, v, shw, shb, 2);
    assert(abs(max(max(Sh,[],1)./min(Sh,[],1))-shw)<0.01)
    empshb=max(max(Sh,[],2)./min(Sh,[],2));
    assert(abs(empshb-shb)<0.1)
%}

%{
    % Generate 100 replicates of a shape matrix.
    shb=1.5;
    shw=4;
    v=5;
    k=20;
    for j=1:100
        disp(j)
        Sh=restrshapeExact(v,k, shw, shb, 2);
        assert(abs(max(max(Sh,[],1)./min(Sh,[],1))-shw)<0.01)
            empshb=max(max(Sh,[],2)./min(Sh,[],2));
            assert(abs(empshb-shb)<0.1)
    end
%}


%% Beginning of code
if shw <1 || shb <1
    error('FSDA:restrshapeExact:Wwrongconstraints','shb and shw cannot be smaller than 1');
end

if nargin<5
    % verLess2016b is true if current version is smaller than 2016b
    verLess2016b=verLessThanFS('9.1');
end

if verLess2016b ==true
    userepmat=1;
else
    userepmat=2;
end

maxiterS=5;
zerotol=1e-12;
tolS=0.01;

% ij=0;
% iterall=zeros(10000,1);

GAMini=zeros(v,k);
diffGAM = 9999;
stoploop=false;
iterouterloop=0;
maxiter=50000;
onedivv=1/v;
while stoploop==false && iterouterloop<maxiter
    iterouterloop=iterouterloop+1;
    for j=1:k
        gam=sort(rand(v,1));
        gam=gam./(prod(gam)^onedivv);
        
        if  j==1
            condshw=false;
            while condshw==false
                if max(gam)/min(gam)>=shw
                    condshw=true;
                else
                    gam=rand(v,1);
                    gam=gam./(prod(gam)^onedivv);
                end
            end
            gam=sort(gam);
        end
        
        
        if  j==2
            condshb=false;
            itercondshb=0;
            maxitercondshb=100000;
            while condshb==false && itercondshb<maxitercondshb
                itercondshb=itercondshb+1;
                if max(gam./GAMini(:,1))>=shb
                    condshb=true;
                else
                    
                    % Generate beta random numbers with parameters 0.5 and
                    % 0.5 in order to increase the probability of large
                    % values for the ratio max/min.
                    % Generate gamma random values and take ratio of the
                    % first to the sum.
                    g1 = randg(0.5,v,1);
                    g2 = randg(0.5,v,1);
                    gam = g1 ./ (g1 + g2);
                    % The 3 above instructions are much faster than
                    % gam=betarnd(0.5,0.5,v,1);
                    gam=sort(gam);
                    
                    % gam=sort(rand(v,1));
                    gam=gam./(prod(gam)^onedivv);
                end
            end
            if itercondshb==maxitercondshb
                error('FSDA:restrshapeExact:Wwrongconstraints','shb is too large in relation to shw please decrease it');
            end
        end
        GAMini(:,j)=gam;
    end
    
    lamGAMc = GAMini;
    
    % Apply the iterative procedure to find constrained \Gamma matrix
    iter=0;
    % diffGAM value of the relative sum of squares of the difference between
    % the element of matrix \Gamma in two consecutive iterations
    diffshw=0;
    while ( (diffGAM > tolS) && (iter < maxiterS) )
        iter = iter + 1;
        
        % In this stage GAM(:,j) is diag(\lambda_j^(1/p) \times \Gamma_j )
        GAM =  lamGAMc;
        
        if diffshw<0
            GAM(end,:)=GAM(end,:)*(1+rand(1,1));
        end
        
        ratempshw=max(GAM,[],1)./min(GAM,[],1);
        booshw=ratempshw>shw;
        
        % Apply eigenvalue restriction inside each group using constraining
        % parameter shw. The ratio of each column of matrix GAM is not
        % greater than shw. Note that restreigen is called just if it is
        % needed.
        for j=1:k
            if booshw(j)==true
                GAM(:,j) = restreigen(GAM(:,j),1,shw,zerotol,userepmat);
            end
        end
        
        lmd=(prod(GAM,1)).^onedivv;
        GAM=GAM./repmat(lmd,v,1);
        
        if k>1
            % Apply restriction between groups
            % The elements of each column of GAM are sorted from largest to smallest
            % The ranks of the orginal ordering of each column is store in matrix Ord
            % The ratio of each row of matrix GAMc is not greater than shb
            niini=100*rand(k,1);
            boo=max(GAM,[],2)./min(GAM,[],2)>shb;
            for i=1:v
                if boo(i)==true
                    GAM(i,:) = restreigen(GAM(i,:), niini, shb, zerotol,userepmat);
                end
            end
            
            lmd=(prod(GAM,1)).^onedivv;
            GAM=GAM./repmat(lmd,v,1);
            
            empshb=max(max(GAM,[],2)./min(GAM,[],2));
            diffshb=empshb-shb;
        else
            diffshb=0;
        end
        ratempshw=max(GAM,[],1)./min(GAM,[],1);
        empshw=max(ratempshw);
        diffshw=empshw-shw;
        diffGAM=max([abs(diffshb) abs(diffshw)]);
    end
    
    % ij=ij+1;
    % iterall(ij)=iter;
    
    % disp(iter)
    if diffGAM<tolS
        stoploop=true;
    end
    
end
if iterouterloop==maxiter
    if diffshw>-0.1
    error('FSDA:restrshapeExact:Wwrongconstraints','shb is too large in relation to shw please decrease shb or increase shw');
    else
    error('FSDA:restrshapeExact:Wwrongconstraints','shb is too small in relation to shw please increase shb or decrease shw');
    end
else
    %  disp(iterouterloop)
end

% disp(iterouterloop)
GAMc=GAM;
end
%FScategory:CLUS-RobClaMULT