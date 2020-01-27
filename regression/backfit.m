function newtX=backfit(ty,tX,X,w,M,l,rsq,maxit,sw,p,delrsq)
% Backfitting algorithm (inner loop) inside ace function

% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit


% ($nit$ = counter for number of iterations in the inner loop).
nit=0;
lfinishInnerLoop=1;

% Iterate until e^2(\theta, \phi_1, ..., \phi_p) fails to decrease
while lfinishInnerLoop==1 % Beginning of Inner Loop
    rsqi=rsq;
    nit=nit+1;
    
    % Find $Z_5 = ty- tX_1 \cdots - tX_p$.
    %  $Z_5$ is the $n \times 1$ vector which contains the residuals
    %  from transformed $y$ and transformed $X$. In what follows $Z_1,
    %  \ldots, Z_5$ are scratch variables (each of dimension $n \times
    %  1$).
    z5=ty-sum(tX,2);
    
    % Loop through all the explanatory variables
    % and transform them
    %   do 420 i=1,p (in Fortran code)
    for  j=1:p
        
        % sorted indexes of the j-th column of X
        ordXj=M(:,j);
        
        %   Find $Z_1= Z_5 + tX_i$ (note that the observations of all the matrices involved are ordered using the ith column of matrix $M$).
        %   Therefore to be precise   $Z_1= Z_5(M(:,i)) + tX_i(M(:,i))$.
        %   $Z_1$ is the $n \times 1$ vector
        %   which contains the residuals from transformed $y$ and transformed $X$ excluding
        %   variable $X_i$ with values ordered according to $X_i$.
        %   \[
        %   Z_1=  ty(M(:,i))- \sum_{j \ne i } tX_j(M(:,i)), \qquad j=1, 2, \ldots, p
        %   \]
        %    The first element of $Z_1$ is associated with the lowest value of $X_i$, the second element is associated with the second lowet value of $X_i$ ...,
        z1=z5(ordXj)+tX(ordXj,j);
        
        % z2 = original expl. variable X(:,j) sorted
        z2=X(ordXj,j);
        
        % z4 = corresponding weights
        z4=w(ordXj);
        
        %  smo=smothr(abs(l(i)),z(:,2),z(:,1),z(:,4));
        % Find smoothed values of X(:,i)
        smo=smothr(abs(l(j)),z2,z1,z4);
        
        % Weighted average of
        sm=sum(smo.*z4)/sw;
        % z3 = z3- mean(z3)
        z3=smo-sm;
        % sv=sum(((z(:,1)-z(:,3)).^2).*w(m(:,i)));
        
        % Compute residual sum of squares
        sv=sum(((z1-z3).^2).*z4);
        
        % Convert residual sum of squares into R2
        % Remember that y has been standardized
        sv=1-sv/sw;
        
        if sv <= rsq
            % lfinishInnerLoop=0;
        else
            % rsq = the multiple R-squared value for the transformed values
            rsq=sv;
            
            tX(ordXj,j)=z3;
            z5(ordXj)=z1-z3;
        end
        
    end       % End of the loop for the explanatory variables
    
    % Condition to exit the Inner Loop
    % There is just one explnatory variable ||
    % The change in R2 is smaller than delrsq
    % The maximum number of iterations has been achieved
    if (p == 1 || rsq-rsqi <= delrsq || nit == maxit)
        lfinishInnerLoop=0;
    end
end % End of the Inner Loop

% Store new tranformed values for the exaplanatory variables
newtX=tX;
end
