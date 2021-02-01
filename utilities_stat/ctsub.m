function a=ctsub(x,y,z)
%ctsub computes numerical integration from x(1) to z(i) of y=f(x).
%
%<a href="matlab: docsearchFS('ctsub')">Link to the help page for this function</a>
%
% For a function defined $y=f(x)$ with $n$ pairs $(x_1,y_1)$ $\ldots$
% $(x_n,y_n)$, with $x_1 \leq x_2 \leq, \ldots, \leq x_n$ this routine
% computes the (approximate) integral using the trapezoidal rule
% \[
% a_i = \int_{x_1}^{z_i} f(x) dx
% \]   
% For further details see "more about".  
%
% Required input arguments:
%
%    x   :      Predictor variable sorted. Vector. Vector of length n
%               containing ordered abscissa values.
%               Note that the x values are assumed non decreasing.
%    y  :       Response variable. Vector. Vector of length n
%               containing ordinate values.
%    z  :       Upper limits of integration. Vector. Vector of length n
%               containing the upper integration limits.
%
%  Optional input arguments:
%
% Output:
%
%    a   :      Result of numerical integration. Vector.  Vector of length
%               n containing the results of the n numerical integrations. 
%               The $i$-th element of $a_i$ with $i=1, 2, \ldots, n$ is 
%               equal to:
%               \[
%               a_i =
%               \int_{x_1}^{z_i} f(x) dx
%               \]   
%
% More About:
% 
% This function estimates the integral of Y via the trapezoidal method.
% For a function defined $y=f(x)$ with $n$ pairs  $(x_1,y_1)$ $\ldots$
% $(x_n,y_n)$, with $x_1 \leq x_2 \leq, \ldots, \leq x_n$,
% if $x_i<z_i \leq x_{i+1}$, $i=1, \ldots, n-1$, this routine computes the
% (approximate) integral  using the trapezoidal rule:
% \[
% a_i =
% \int_{x_1}^{z_i} f(x) dx
% \]
% More precisely
% \begin{equation}\label{ai}
% a_i= \sum_{j=2}^{i} 0.5 (x_j-x_{j-1})(y_j+y_{j-1})
% +0.5 (z_i-x_i) \left\{ 2y_i+(z_i-x_i) \frac{y_{i+1}-y_i}{x_{i+1}-x_i} \right\}
%  \end{equation}
% The last term of the equation is the area of the trapezoid with
% coordinates $(x_i, y_i)$,$(z_i, f(z_i))$  and 
% \[
% f(z_i)=y_i+(z_i-x_i) \frac{y_{i+1}-y_i}{x_{i+1}-x_i} 
% \]
% is found by linear interpolation.
% 
% If $z_i>x_n$ the function for $x > x_n$ is assumed constant and equal to
% $y_n$ therefore to the expression of $a_i$ computed as described above we
% must add $(z_i-x_n) y_n$.
% On the other hand, if $z_i<x_1$ the function for $x<x_1$ is assumed
% constant and equal to $y_1$ therefore $a_i$ is computed as:
% \[
% a_i = -(x_1-z_i) y_1
% \]
% Note that $a_i$ in this last case (if $y_1$ is positive) becomes negative.
%
% This routine in called in every step of the outer loop by function avas.m
% in order to compute a new set of transformed values for the response
% which have approximately constant variance.
% Inside avas, the $x$ coordinates are the fitted values ordered, the $y$
% coordinates are the reciprocal of the smoothed absolute values of
% residuals sorted using the ordering of fitted values, while the upper
% values of the range of integration are given by the elements of vector
% $ty$ sorted using the ordering of fitted values. The output is the new
% vector $ty$ with the elements ordered using $ordyhat$. 
%
% See also: avas, rlsmo
%
% References:
% Tibshirani R. (1987), Estimating optimal transformations for regression,
% "Journal of the American Statistical Association", Vol. 83, 394-405.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('ctsub')">Link to the help page for this function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%{
    % Transform a linear relation.
    n=1000;
    x=sort(randn(n,1))+10;
    y=1+2*x+ randn(n,1);
    % y(50:70)=-5;
    z=randn(n,1);
    a=ctsub(x,y,z);
    subplot(2,1,1)
    plot(x,y)
    subplot(2,1,2)
    plot(x,a)
%}

%{
    % Comparison with MATLAB function cumtrapz.
    n=1000;
    x=sort(randn(n,1));
    y=randn(n,1);
    % If the third argument of ctsub is equal to the first
    % argument then the output of cumtraps and ctsub is exactly the same.
    res=cumtrapz(x,y);
    res1=ctsub(x,y,x);
    disp(max(abs(res-res1)))
%}

%{
    % Apply ctsub to heteroskedastic data.
    n=1000;
    x=10*sort(randn(n,1))+10;
    % The variance of y depends on x.
    y=1+1*x+ 10*x.*randn(n,1);
    % Upper limits of integration.
    z=rand(n,1)*100-50;
    a=ctsub(x,y,z);
    subplot(2,1,1)
    plot(x,y)
    title('Original data')
    subplot(2,1,2)
    plot(x,a)
    title('Transformed data after the variance stabilizing transformation')
%}
    
% Note that trapz([1/3 pi],[sin(1/3) sin(pi)]) =
% ctsub([1/3 pi],[sin(1/3) sin(pi)],[[1/3 pi]])


%% Beginning of code
n =length(z);
a =zeros(n,1);
for i=1:n
    
    if(z(i)<= x(1))
        a(i) =(z(i)-x(1))*y(1);
    else
        j =1;
        a(i) =0;
        while (j<=n && z(i)>x(j))
            if j > 1
                a(i) =a(i)+(x(j)-x(j-1))*(y(j)+y(j-1))/2;
            end
            j =j+1;
        end
        
        if z(i) <= x(n)
            a(i) = a(i)+.5*(z(i)-x(j-1))*(2*y(j-1)+(z(i)-x(j-1))*(y(j)-y(j-1))/(x(j)-x(j-1)));
        else
            a(i) =a(i)+(z(i)-x(n))*y(n);
        end
    end
end
%FScategory:UTISTAT

