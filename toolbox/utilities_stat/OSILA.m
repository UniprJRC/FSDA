function out=OSILA(X,k,alpha,n)
%OSILA finds the k-th order statistic
%
% <a href="matlab: docsearchFS('OSILA')">Link to the help function</a>
%
% OSILA is an algorithm for computing order statistics through an approach 
% based on random sampling. For arrays of dimension bigger than 10^4,
% it allows to achieve the exact solution remarkably faster than the naive
% approach (i.e. sorting all the elements and take the k-th one).
%
% Required input arguments:
%
%   X:  a set of N numbers. Vector. Vector containing a set of N numbers.
%                 Data Type - double
%   k:  order statistic index. Scalar. An integer between 1 and N indicating 
%       the desired order statistic.
%                 Data Type - double
%
%
%  Optional input arguments:
%
%  alpha: confidence level, it provides an upper value for the probability 
%        to have more than one iteration in the algorithm. If not provided 
%        (or empty), it is set to 0.99.  
%                   Example - 0.99
%                   Data Types - double
%  n: dimension of the random sample. If not provided (or n=0 or empty), it
%       is calculated through the subfunction 'FindOptimalDimension' (see
%       Section 3 of the reference).  
%                   Example - 1000
%                   Data Types - double
%
% Output:
%
% out: k-th order statistic. Scalar.
%                 Data Type - double
%
%
% See also:  quickselectFS
%
% References:
%
% Cerasa, A. (2023). "Order statistics in large arrays (OSILA): a simple 
% randomised algorithm for a fast and efficient attainment of the order 
% statistics in very large arrays." Computational Statistics, p. 1-26.
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('OSILA')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:
%
%{
    %% OSILA: standard application - alpha=0.99 (default value), n calculated.
    N=10^7;
    Y=randn(N,1);
    k=34580;
    out=OSILA(Y,k);
    % Check the result
    sorY=sort(Y);
    disp([out,sorY(k)])
    disp(['Distance = ' num2str(abs(out-sorY(k)))])
%}

%{
    %% OSILA with alpha supplied and n calculated.
    N=10^7;
    Y=randn(N,1);
    k=34580;
    n=50;
    s1=tic;
    out1=OSILA(Y,k,0.01);
    t1=toc(s1);
    s2=tic;
    out2=OSILA(Y,k);
    t2=toc(s2);
    disp(['Execution time with n = 50: ' num2str(t1)])
    disp(['Execution time with n optimal: ' num2str(t2)])
    disp(['Execution time difference: ' num2str(t1-t2)])
%}

%{
    %% OSILA with n (sample dimension) supplied and alpha empty (default value).
    N=10^7;
    Y=randn(N,1);
    k=34580;
    n=50;
    s1=tic;
    out1=OSILA(Y,k,[],n);
    t1=toc(s1);
    s2=tic;
    out2=OSILA(Y,k);
    t2=toc(s2);
    disp(['Execution time with n = 50: ' num2str(t1)])
    disp(['Execution time with n optimal: ' num2str(t2)])
    disp(['Execution time difference: ' num2str(t1-t2)])
%}

%{
    %% OSILA: Execution time difference with respect to the naive procedure.
    N=10^7;
    Y=randn(N,1);
    k=34580;
    s1=tic;
    out=OSILA(Y,k);
    t1=toc(s1);
    s2=tic;
    sorY=sort(Y);
    out=sorY(k);
    t2=toc(s2);
    disp(['Execution time OSILA: ' num2str(t1)])
    disp(['Execution time naive: ' num2str(t2)])
    disp(['Execution time difference: ' num2str(t1-t2)])
%}



%% Beginning of code
if nargin<3 || isempty(alpha)
    alpha=0.99;
end
N=numel(X);
zAlpha=norminv(alpha);
if nargin<4 || isempty(n) || n==0
    n=FindOptimalDimension(N,k,zAlpha);
end
j0=min(round((n+1)*k/(N+1)),n); % Between 1 and n
j0=max(j0,1);
out=Inf;
intervalX=[NaN NaN];
intervalPos=[-Inf Inf];
coin=randperm(N,n);
xCoin=X(coin);
[x0,xCoinOrd]=OrderStatistics(xCoin,j0);
posLessX0=(X<=x0);
k0=sum(posLessX0);
if k0==k
    out=x0;
elseif k0>k
    intervalX(2)=x0;
    intervalPos(2)=k0;
    while isnan(intervalX(1))
        j0=findJ0down(k0,j0,k,zAlpha);
        x0=xCoinOrd(j0);
        posLessX0=(X<=x0);
        k0=sum(posLessX0);
        if k0<=k
            intervalX(1)=x0;
            intervalPos(1)=k0;
        else
            intervalX(2)=x0;
            intervalPos(2)=k0;
            if j0==1
                intervalX(1)=-Inf;
                intervalPos(1)=1;
            end
        end
    end
else
    intervalX(1)=x0;
    intervalPos(1)=k0;
    while isnan(intervalX(2))
        j0=j0+findJ0up(N,n,k0,j0,k,zAlpha);
        x0=xCoinOrd(j0);
        posLessX0=(X<=x0);
        k0=sum(posLessX0);
        if k0>=k
            intervalX(2)=x0;
            intervalPos(2)=k0;
        else
            intervalX(1)=x0;
            intervalPos(1)=k0;
            if j0==n
                intervalX(2)=Inf;
                intervalPos(2)=N;
            end  
        end
    end    
end 
if out==Inf
    newX=X(X>=intervalX(1) & X<=intervalX(2));
    newX=sort(newX);
    newPos=k-intervalPos(1)+1;
    out=newX(newPos);
end
end


function [ordStat,xSort]=OrderStatistics(x,pos)
xSort=sort(x);
ordStat=xSort(pos);
end


function j0=findJ0down(k0,j0,k,zAlpha)
A=k0/j0;
B2=(zAlpha^2)*(k0*(k0-j0))/((j0^2)*(j0+1));
a=-A^2-B2;
b=B2*j0+2*A*k;
c=-k^2;
r1=(-b+sqrt(b^2-4*a*c))/(2*a);
j0=max(floor(r1),1);
end



function j0=findJ0up(N,n,k0,j0,k,zAlpha)
A=(N-k0+1)/(n-j0+1);
B2=(zAlpha^2)*(N-k0+1)*(N-k0-n+j0)/(((n-j0+1)^2)*(n-j0+2));
a=A^2+B2;
b=2*A*(k0-k)-n*B2+j0*B2-B2;
c=(k0-k)^2;
r1=(-b+sqrt(b^2-4*a*c))/(2*a);
j0=min(ceil(r1),n-j0);
end


function n=FindOptimalDimension(N,k,zAlpha)
if k==N
    k=N-1;
end
A=sqrt((2*k/pi)*(N-k+1)/(N+1));
B=zAlpha/4;
a=2*A*B*(k^(-0.5)+(N-k)^(-0.5));
b=A+B*k^0.5-B*k^(-0.5)+B*(N-k)^0.5;
c=0;
d=0;
e=-2*N;
p=(8*a*c-3*b^2)/(8*a^2);
q=12*a*e-3*b*d+c^2;
s=27*a*d^2-72*a*c*e+27*(b^2)*e-9*b*c*d+2*c^3;
S=(8*(a^2)*d-4*a*b*c+b^3)/(8*a^3);
delta0=(s+sqrt(s^2-4*q^3))/2;
delta0=sign(delta0)*abs(delta0)^(1/3);
Q=0.5*sqrt(-(2/3)*p+(1/(3*a))*(delta0+q/delta0));
b4a=-b/(4*a);
QpSplus=0.5*sqrt(-4*Q^2-2*p+S/Q);
r1=b4a-Q+QpSplus;
n=round(N/r1^2);
end

%FScategory:UTISTAT




