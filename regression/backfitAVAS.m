function [newtX, rsq]=backfitAVAS(ty,tX,X,w,M,l,rsq,maxit,sw,p,delrsq,bsb,outliers,PredictorOrderR2)
%backfitAVAS contains the backfitting algorithm (inner loop) inside avas function
%
% This funtion is not intended to be called directly
%
% Copyright 2008-2021.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit
n=length(ty);
if nargin<13
    bsb=1:n;
    outliers=[];
    PredictorOrderR2=false;
elseif nargin<14
    PredictorOrderR2=false;
end

% ($nit$ = counter for number of iterations in the inner loop).
nit=0;
lfinishInnerLoop=1;

% Iterate until e^2(\theta, \phi_1, ..., \phi_p) fails to decrease
while lfinishInnerLoop==1 % Beginning of Inner Loop
    rsqi=rsq;
    nit=nit+1;
    
    % Find $e = ty- tX_1 \cdots - tX_p$.
    %  $e$ is the $n \times 1$ vector which contains the residuals
    %  from transformed $y$ and transformed $X$. In what follows $Z_1,
    %  \ldots, Z_5$ are scratch variables (each of dimension $n \times
    %  1$).
    e=ty-sum(tX,2);
    
    % Loop through all the explanatory variables
    % and transform them
    %   do 420 i=1,p (in Fortran code)
    
    if PredictorOrderR2 ==true && p>1
        varin=false(p,1);
        varout=true(p,1);
        ChooseBestR2=true;
        wbsb=w(bsb);
        sqweightsbsb = wbsb.^(1/2);
        tybsb=ty(bsb);
        tywbsb = tybsb .* sqweightsbsb;
        tywav=(tybsb.*wbsb)/sum(wbsb);
        tyminusav=(tybsb-tywav).*sqweightsbsb;
        DEVtotbsb=tyminusav'*tyminusav;
        interceptbsb=ones(length(bsb),1);
        seqp=1:p;
    else
        ChooseBestR2=false;
    end
    % avas is dependent on the order of the predictors in X.
    
    for  ij=1:p
        
        if ChooseBestR2 == true
            tXbsbVarout=tX(bsb,varout);
            % Indexes (interger numbers in 1, ...p) containing the
            % variables not included yet in the transformation 
            tXcandIndexes=seqp(varout);
            R2cand=zeros(length(tXcandIndexes),1);
            % Matrix tX referred to the explanatory variables already
            % included
            tXbsbVarin=tX(bsb,varin);
            for jj=1:size(tXbsbVarout,2)
                tXjjbsb=[interceptbsb tXbsbVarin tXbsbVarout(:,jj)];
                tXjjw=tXjjbsb.*sqweightsbsb;
                % estimate of beta from weighted regression (RWLS)
                bjj=tXjjw\tywbsb;
                DEVresbsb=sum(((tybsb-tXjjbsb*bjj).^2).*wbsb);
                % Note that ty is standardized
                R2jj=1-DEVresbsb/DEVtotbsb;
                R2cand(jj)=R2jj;
            end
            [~,indmax]=max(R2cand);
            j=tXcandIndexes(indmax);
            % Variable j is selected as the one producing the best R2 given
            % those already in
            varin(j)=true;
            varout(j)=false;
        else
            j=ij;
        end
        
        % sorted indexes of the j-th column of X
        ordXj=M(:,j);
        if ~isempty(outliers)
            ordXj=intersect(ordXj,bsb,'stable');
        end
        
        %   Find eplustXjOrdXj (eplustXjOrdXj is Z_1 is the
        %   original fortran code)
        %   (note that the observations of all the matrices involved are
        %   ordered using the j-th column of matrix $M$).
        %   Therefore to be precise   $eplustXjOrdXj= e(M(:,j)) + tX_i(M(:,j))$.
        %   $eplustXjOrdXj$ is the $n \times 1$ vector
        %   which contains the residuals from transformed $y$ and transformed $X$, $tX$ excluding
        %   variable $X_j$ with values ordered according to $X_j$.
        %   \[
        %   eplustXjOrdXj=  ty(M(:,j))- \sum_{i \ne j } tX_j(M(:,j)), \qquad i=1, 2, \ldots, p
        %   \]
        %    The first element of $eplustXjOrdXj$ is associated with the
        %    lowest value of $X_j$, the second element is associated with
        %    the second lowet value of $X_j$ ...,
        eplustXjOrdXj=e(ordXj)+tX(ordXj,j);
        
        
        % xjord = original expl. variable X(:,j) sorted
        Xjord=X(ordXj,j);
        
        % wOrdXj  corresponding weights (z4 in the fortran code)
        wOrdXj=w(ordXj);
        
        %  smo=smothr(abs(l(i)),z(:,2),z(:,1),z(:,4));
        % Find smoothed values of X(:,i)
        smoOrdXj=smothr(abs(l(j)),Xjord,eplustXjOrdXj,wOrdXj);
        
        % Weighted average of
        sm=sum(smoOrdXj.*wOrdXj)/sw;
        % z3 = z3- mean(z3)
        tXOrdXj=smoOrdXj-sm;
        % sv=sum(((z(:,1)-z(:,3)).^2).*w(m(:,i)));
        
        % Compute residual sum of squares
        sv=sum(((eplustXjOrdXj-tXOrdXj).^2).*wOrdXj);
        
        % Convert residual sum of squares into R2
        % Remember that y has been standardized
        sv=1-sv/sw;
        
        %%%%%%%% backfitAVAS does not have this if
        %         if sv <= rsq
        %         else
        % rsq = the multiple R-squared value for the transformed values
        rsq=sv;
        
        tX(ordXj,j)=tXOrdXj;
        e(ordXj)=eplustXjOrdXj-tXOrdXj;
        %         end
        
        % Now compute fitted values for transformed variable Xj in correspondence of the outlierss
        [XjordUni,XjOrdInd]=unique(Xjord);
        tXjoutliers=interp1(XjordUni,tXOrdXj(XjOrdInd),X(outliers,j),'linear','extrap');
        tX(outliers,j)=tXjoutliers;
        
    end       % End of the loop for the explanatory variables
    
    % Condition to exit the Inner Loop
    % There is just one explnatory variable ||
    % The change in R2 is smaller than delrsq
    % The maximum number of iterations has been achieved
    if (p == 1 || abs(rsq-rsqi) <= delrsq || nit == maxit)
        lfinishInnerLoop=0;
    end
end % End of the Inner Loop

% Store new tranformed values for the explanatory variables (also included
% are the estimated transformed values for the rows declared as outliers).
newtX=tX;
end
