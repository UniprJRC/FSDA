function y = corrpdf(r, rho, n)
%corrpdf correlation coefficient probability density function
%
% Density function of the correlation coefficient under the null hypothesis
% that the n bivariate observations come from a bivariate normal
% distribution with correlation parameter rho
%
%<a href="matlab: docsearchFS('corrpdf')">Link to the help function</a>
%
% Required input arguments:
%
%    r:         Value at which the pdf must be evaluated.
%               Scalar, vector or matrix or 3D array.
%               See "More About:" for details about the distribution
%               of the correlation coefficient.
%               Data Types - single | double
%    rho :      value of the correlation coefficient in the population.
%               Scalar, vector or matrix or 3D array. If rho is not a
%               scalar all the 3 input arguments (r,rho and n) must have
%               the same size or just the numel of one of the 3 input
%               arguments must be greater than 1
%               Data Types - single | double
%    n :        sample size.
%               Scalar, vector or matrix or 3D array. If n is not a scalar
%               all the 3 input arguments (r,rho and n) must have
%               the same size or just the numel of one of the 3 input
%               arguments must be greater than 1
%               Data Types - single | double
%
% Optional input arguments:
%
%
%  Output:
%
% y:    PDF value. Scalar, vector or matrix or 3D array of the same size
%               of input arguments x, rho and n. This is the probability
%               density function for r, given n bivariate data coming from
%               bivariate normal distribution with parameter rho.
%
%
%
% See also: corrcdf
%
% References:
%
% Das Gupta, S. (1980). Distribution of the Correlation Coefficient,
% in: Fienberg, S.E., Hinkley, D.V. (eds) R.A. Fisher: An Appreciation, 
% Lecture Notes in Statistics, vol 1. Springer, New York, NY. 
% https://doi.org/10.1007/978-1-4612-6079-0_3
%
% Acknowledgements:
%
% For additional information see
% https://mathworld.wolfram.com/CorrelationCoefficientBivariateNormalDistribution.html
% This function follows the lines of MATLAB code developed by Xu Cui,
% https://www.alivelearn.net/?p=709 Stanford University and the file
% exchange submission Joshua Carmichael (2022). sample correlation
% distribution function
% https://www.mathworks.com/matlabcentral/fileexchange/45785-sample-correlation-distribution-function/
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('corrpdf')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%

%{
    % An example where x, rho and n are all scalars.
    % Find f(r=x|rho=0.1|n=12)
    x=0;
    rho=0.1;
    n=12;
    xs=num2str(x);
    rhos=num2str(rho);
    ns=num2str(n);
    dens=corrpdf(x,rho,n);
    disp(['f(r=' xs '|rho=' rhos, ', n=' ns ')'])
    disp(dens);
%}

%{
    %An example where x is not scalar.
    x=-1:0.01:1;
    rho=0.8;
    n=12;
    rhos=num2str(rho);
    ns=num2str(n);
    dens=corrpdf(x,rho,n);
    % disp(['f(r=x|rho=' rhos, ', n=' ns ')'])
    plot(x,dens)
    xlabel("x")
    ylabel(['f(r=x|rho=' rhos, ', n=' ns ')'])
%}

%{
    %% Show density of r given different values for rho.
    x=-1:0.01:1;
    rho=[0:0.2:0.8 0.9]';
    n=12;
    rhos=num2str(rho);
    ns=num2str(n);
    close all
    for i=1:length(rho)
        nexttile
        dens=corrpdf(x,rho(i),n);
        plot(x,dens)
        title(['f(r|n=' ns ', rho=' rhos(i,:) ')'])
        xlabel("r")
    end
    
    sgtitle('Density of the correlation coefficient for different values of rho')
%}

%{
    %% An example where rho is not scalar.
    x=0.3;
    rho=(0:0.1:0.8)';
    n=12;
    xs=string(x);
    rhos=string(rho);
    ns=string(n);
    Dens=corrpdf(x,rho,n);
    nameRows="f(r="+xs+"|rho="+ rhos+ ", n="+ ns+ ")=";
    nameRowsT=array2table(Dens,"RowNames",nameRows);
    disp(nameRowsT)
%}

%{
    % An example where n is not scalar.
    x=0.3;
    rho=0';
    n=(5:5:50)';
    xs=string(x);
    rhos=string(rho);
    ns=string(n);
    Dens=corrpdf(x,rho,n);
    nameRows="f(r="+xs+"|rho="+ rhos+ ", n="+ ns+ ")=";
    nameRowsT=array2table(Dens,"RowNames",nameRows);
    disp(nameRowsT)
%}

%{
    % An example where n, r and rho have the same dimension.
    x=(0.2:0.1:0.6)';
    rho=(0.1:0.1:0.5)';
    n= (10:10:50)';
    xs=string(x);
    rhos=string(rho);
    ns=string(n);
    Dens=corrpdf(x,rho,n);
    nameRows="f(r="+xs+"|rho="+ rhos+ ", n="+ ns+ ")=";
    nameRowsT=array2table(Dens,"RowNames",nameRows);
    disp(nameRowsT)
%}

%% Beginning of code

if numel(r)>1 && numel(n)==1 &&  numel(rho)==1
    n=repmat(n,size(r));
    rho=repmat(rho,size(r));
    y=zeros(size(rho));
elseif numel(r)==1 && numel(n)>1 &&  numel(rho)==1
    r=repmat(r,size(n));
    rho=repmat(rho,size(n));
    y=zeros(size(n));
elseif numel(r)==1 && numel(n)==1 &&  numel(rho)>1
    r=repmat(r,size(rho));
    n=repmat(n,size(rho));
    y=zeros(size(rho));
else
    % In this case rho, x and n all have the same size
    y=zeros(size(rho));
end

boo=n<120;

% For the values of n smaller than 120
rb=r(boo);
nb=n(boo);
rhob=rho(boo);

y(boo) = (nb-2).*gamma(nb-1) .* ((1-rhob.^2).^((nb-1)/2)).* (1-rb.^2).^((nb-4)/2);
y(boo) =y(boo)./ (sqrt(2*pi) .* gamma(nb-1/2) .* (1-rhob.*rb).^(nb-3/2)); % .* ...
y(boo) =y(boo).*(1+ 1/4*(rhob.*rb+1)./(2*nb-1) + 9/16*(rhob.*rb+1).^2 ./ (2*nb-1)./(2*nb+1));

% For the values of n greater or equal than 120
rnb=r(~boo);
nnb=n(~boo);
rhonb=rho(~boo);
y(~boo) =  (1-rhonb.^2).^((nnb-1)/2) .* (1-rnb.^2).^((nnb-4)/2).* (nnb-2);
y(~boo) =y(~boo)./  (sqrt(2*pi) .* (1-rhonb.*rnb).^(nnb-3/2)) .* nnb.^(-1/2);
y(~boo) =y(~boo).*(1+ 1/4*(rhonb.*rnb+1)./(2*nnb-1) + 9/16*(rhonb.*rnb+1).^2 ./ (2*nnb-1)./(2*nnb+1));

% Set to 0 the density for the values of r outside the admissible range [-1
% 1]
y(r>1)              = 0;
y(r<-1)             = 0;
y(~isfinite(y))     = 0;
end
%FScategory:UTISTAT