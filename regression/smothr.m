function ysmo=smothr(l,x,y,w)
%smothr produces smoothed values with constraints
%
%
%<a href="matlab: docsearchFS('smothr')">Link to the help page for this function</a>
%
%   This function is used in each step of the iterative procedure for ACE
%   but can be called directly when it is necessary to smooth a set of
%   values. Note that the x values must be non decreasing.
%
%
% Required input arguments:
%
%       l :  type of transformation. Scalar. Scalar which
%           specifies how the type of transformation.
%           l=1 => transformation can also be non monotone. In this case
%                  the supersmoother is initially applied. In presence of
%                  equal values of x the unweighted arithmetic mean of the
%                  smoothed values is returned.
%           l=2 => j-th variable assumes circular (periodic) values
%                 in the range (0.0,1.0) with period 1.0.
%           l=3 => transformation is to be monotone. In this case the
%                  supersmoother is initially applied. Monotonic
%                  transformation is forced applying isotonic regression to
%                  (1) the output of the supersmoother and to the flipped
%                  upside down (2) output of the supersmoother. The choice
%                  between solution (1) and (2) is made taked the output
%                  which is closest to the output of the supersmoother.
%                  Closeness is measured in terms of sum
%                  of squares of residuals. Equal consecutive values
%                  smoothed values are replaced by linearly iterpolated values.
%                  In presence of equal values of x, the unweighted
%                  arithmetic mean of the final smoothed values is
%                  returned.
%           l=4 => transformation is to be linear. In this case the smoothed
%                  values are simply the fitted values from least squares
%                  fit.
%           l=5 => the predictor variable is categorical. In this case the smoothed
%                  values are simply the (weighted) values of y in
%                  correspondence of each value of x.
%    x   :      Predictor variable sorted. Vector.  Ordered abscissa values.
%               Note that the x values are assumed non decreasing.
%    y  :       Response variable. Vector. Response variable which has to
%               be smoothed, specified as
%               a vector of length n, where n is the number of
%               observations.
%
%
%  Optional input arguments:
%
%       w  : weights for the observations. Vector. Row or column vector of
%           length n containing the weights associated to each
%           observations. If w is not specified we assum $w=1$ for $i=1,
%           2, \ldots, n$.
%           Example - 'w',1:n
%           Data Types - double
%
%
%  Output:
%
%         ysmo:  smoothed values. Vector. A vector with the same dimension
%               of y containing smoothed values, that is the y values on
%               the fitted curve. The smoothed values come from linear
%               regression if input value l=4. The smoothed values are
%               monotonic if input value l=3. The smoothed values can also
%               be non monotonic if input value l=1;
%
%
% See also: ace.m, supsmu.m, avas.m, rlsmo.m
%
% References:
%
% Breiman, L. and Friedman, J.H. (1985), Estimating optimal transformations
% for multiple regression and correlation, "Journal of the American
% Statistical Association", Vol. 80, pp. 580-597.
% Wang D.  and Murphy M. (2005), Identifying nonlinear relationships
% regression using the ACE algorithm, "Journal of Applied Statistics", Vol.
% 32, pp. 243-258.
% Friedman, J.H. (1984), A variable span scatterplot smoother. Laboratory
% for Computational Statistics, Stanford University, Technical Report No. 5.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('smothr')">Link to the help function</a>
%
%$LastChangedDate:: 2018-06-08 01:27:40 #$: Date of the last commit

% Examples:

%{
    %% Compare 4 different smoothers.
    % The data give the speed of cars and the distances taken to stop. Note
    % that the data were recorded in the 1920s.
    % The first column of X is speed while the second is the time to stop.
    X=[ 4    2
        4   10
        7    4
        7   22
        8   16
        9   10
        10   18
        10   26
        10   34
        11   17
        11   28
        12   14
        12   20
        12   24
        12   28
        13   26
        13   34
        13   34
        13   46
        14   26
        14   36
        14   60
        14   80
        15   20
        15   26
        15   54
        16   32
        16   40
        17   32
        17   40
        17   50
        18   42
        18   56
        18   76
        18   84
        19   36
        19   46
        19   68
        20   32
        20   48
        20   52
        20   56
        20   64
        22   66
        23   54
        24   70
        24   92
        24   93
        24  120
        25   85];
    x=X(:,1);
    y=X(:,2);
    % Compare the output
    subplot(2,2,1)
    % Non monotonic output l=1
    plot(x,y,'o')
    hold('on')
    l=1;
    ysmo=smothr(l,x,y);
    plot(x,ysmo,'-','LineWidth',2)
    title('l=1: non monotonic transformation')
    ylabel('Distance')

    subplot(2,2,2)
    % Impose monotonic output
    % Input option l=3
    plot(x,y,'o')
    hold('on')
    l=3;
    ysmo=smothr(l,x,y);
    plot(x,ysmo,'-','LineWidth',2)
    title('l=1: monotonic transformation')

    subplot(2,2,3)
    % Impose monotonic output
    % Input option l=4
    plot(x,y,'o')
    hold('on')
    l=4; % Impose a linear smoother
    ysmo=smothr(l,x,y);
    plot(x,ysmo,'-','LineWidth',2)
    title('l=4: linear transformation')
    xlabel('Speed')
    ylabel('Distance')

    subplot(2,2,4)
    % Impose monotonic output
    plot(x,y,'o')
    hold('on')
    ysmo=supsmu(x,y);
    plot(x,ysmo,'-','LineWidth',2)
    title('Supersmoother with all the default options')
    xlabel('Speed')
    % ylabel('Distance')
%}

%{
    % An example where the predicted variable is categorical.
    seed=20;
    n=5;
    y1=exp(-0.5+0.5*mtR(n,1,seed));
    y2=exp(0.5+0.5*mtR(n,1,-seed));
    y=[y1;y2];
    X=[-0.5*ones(n,1); 0.5*ones(n,1)];
    X(9:10)=1;
    ysmo=smothr(5,X,y);
    plot(X,y,'o')
    hold('on')
    plot(X,ysmo)
    xlabel('Variable with 3 levels')
    ylabel('Original and smoothed y values')
%}

%% Beginning of code
n=length(y);

if nargin<4
    w=ones(n,1);
end

if l==4 % Transformation is forced to be linear
    %     % Approach to find Xsmo based on loops
    %     % Old way of programming
    %        n=length(y);
    %         Xsmo=zeros(n,1);
    %         sm=0.0;
    %         sw=sm;
    %         b=sw;
    %         d=b;
    %         for j=1:n
    %             sm=sm+w(j)*x(j)*y(j);
    %             sw=sw+w(j)*x(j)^2;
    %             b=b+w(j)*x(j);
    %             d=d+w(j);
    %         end
    %         a=sm/(sw-(b^2)/d);
    %         b=b/d;
    %         for j=1:n
    %             Xsmo(j)=a*(x(j)-b);
    %         end
    
    % Compute LS line without loops
    %     yori=y;
    %     xori=x;
    %     fitlm(xori,yori,'Weights',1:n )
    
    % The formulae below assume that x and y are mean centered.
    sumw=sum(w);
    sumxw=sum(x.*w);
    % meanx = mean of x
    meanx=sumxw/sumw;
    % The formulae below assume that y is mean centered.
    meany=sum(y.*w)/sumw;
    y=y-meany;
    sumxyw=sum(w.*x.*y);
    sumx2w=sum((x.^2).*w);
    a=sumxyw/(sumx2w-(sumxw^2)/sumw);
    % The slope of the LS line is a
    % The intercept of the LS line is -a*meanx+meany
    ysmo=a*(x-meanx)+meany;
    
    % The approach below to find ysmo is slower
    %         X=[ones(length(x),1) x];
    %         wsqrt=sqrt(w);
    %         Xw=X.*wsqrt;
    %         yw=yori.*wsqrt;
    %         bb=Xw\yw;
    %         ysmo=bb(1)+x*bb(2);
elseif l==5 % variable is categorical 
    ysmo=zeros(n,1);
   % In this case the smoothed values are equal to the local weighted
   % arithmetic mean of y in correspondence of the equal consecutive
   % values of x. Note that the x values are preliminarly ordered therefere
   % the weighted arithmetic mean of equal consecutive values considers all the
   % values of x which have a particular value.
    j0=1;
    salta=0;
    for j=1:n-1
        if salta==0
            sm=y(j)*w(j);
            sw=w(j);
        end
        
        if x(j+1) <=x(j)
            sm=sm+y(j+1)*w(j+1);
            sw=sw+w(j+1);
            salta=1;
            if j==n-1
                sm=sm/sw;
                ysmo(j0:j+1)=sm;
            end
        else
            salta=0;
            sm=sm/sw;
            ysmo(j0:j)=sm;
            j0=j+1;
        end
    end
    
else % Transformation must not necessarily be linear
    
    alpha=5;
    % Apply the super smoother and find the smoothed values of y
    if l~=2
        Period=Inf;
    else
        % Note that in this case the x values are assumed to be in [0, 1] and of period 1.
        % Therefore x(end+1)=x(1)+period x(end+2) = x(2)+period ..... x(end+k-1)=x(k-1)
        % y(end+1)=y(1)  y(end+2) = y(2) ..... y(end+k-1)=y(k-1), that is
        % x = [x;x(1:k-1)+period];
        % y = [y;y(1:k-1)];
        % The smoothed values in the first k-1 positions are added in the
        % last k-1 positions
        Period=1;
    end
    [ysmo]=supsmu(x,y,'Alpha',alpha,'Weights',w,'Period',Period);
    
    
    if l==3 % if l=3 force the transformation to be monotonic
        
        scr1=ysmo;
        scr2=flip(ysmo,1);
        
        % Smoothed y values forced to be monotonic
        % Force scr1 and scr2 to be monotonically non decreasing
        % Replace consecutive non decreasing values with the corresponding
        % arithmetic mean (isotonic regression)
        scr1=montne(scr1);
        scr2=montne(scr2);
        
        % Choose between scr1 and scr2: the one which is closer in terms of
        % residual sum of squares
        sm=sum((ysmo-scr1).^2);
        scr2flipped=flip(scr2,1);
        sw=sum((ysmo-scr2flipped).^2);
        
        if sm<sw
            ysmo=scr1;
        else
            ysmo=scr2flipped;
        end
        
        % Replace equal consecutive values of ysmo with linearly interpolated values
        for j=1:n-1
            j0=j;
            jp=j;
            if ysmo(j+1)==ysmo(j)
                
                while  ysmo(jp+1)==ysmo(jp)
                    jp=jp+1;
                    if jp==n
                        break
                    end
                end
                
                a=0.0;
                if j0 >1
                    a=0.5*(ysmo(j0)-ysmo(j0-1));
                end
                b=0.0;
                if jp<n
                    b=0.5*(ysmo(jp+1)-ysmo(jp));
                end
                d=(a+b)/(jp-j0);
                if a== 0.0 || b == 0.0
                    d=2.0*d;
                end
                if a ==0.0
                    a=b;
                end
                ii=(j0:jp)';
                ysmo(ii)=ysmo(ii)-a+d*(ii-j0);
            end
        end % End of loop to replace equal consecutive values of ysmo
    end % End of l=3 (monotonic transformation)
    
    % Loop which deals with equal consecutive values of x
    % row 397 of Fortran program. The values of ysmo in correspondence of
    % the equal values of x are equal to the simple average of the values
    % of y.
    j0=1;
    salta=0;
    for j=1:n-1
        if salta==0
            sm=ysmo(j);
        end
        
        if x(j+1) <=x(j)
            sm=sm+ysmo(j+1);
            salta=1;
            if j==n-1
                sm=sm/(j-j0+2);
                ysmo(j0:j+1)=sm;
            end
        else
            salta=0;
            sm=sm/(j-j0+1);
            ysmo(j0:j)=sm;
            j0=j+1;
        end
    end
end

end
%FScategory:REG-Transformations