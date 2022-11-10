function bw = bwe(X, bwopt)
%bwe estimates the bandwidth smoothing parameter for kernel density estimation.
%
%<a href="matlab: docsearchFS('bwe')">Link to the help page for this function</a>
%
%  Required input arguments:
%
%   X :          Input data. Vector or matrix. The data to be smoothed by
%                kernel density estimation.
%
%  Optional input arguments:
%
%   bwopt :      Estimation method. String. Default is Scott's rule.
%                Other options are:
%                - 'normal', the normal reference rule, applied only for
%                  d=1. It is valid if the underlying density being
%                  estimated is Gaussian.
%                - 'robust', is the normal reference rule, applicable in
%                  presence of outliers, again for d=1.
%                Data Types - char
%                Example - 'method','robust'
%
%
%  Output:
%
%   bw :        bandwidth estimate. Vector or Scalar. It is a scalar if the
%               data is uni-dimensional, otherwise is a vector with a
%               bandwidth value for each dimension.
%
%
%
%
%  See also: ksdensity
%
%
% References:
%
% Bowman, A.W. and Azzalini, A. (1997), "Applied Smoothing
% Techniques for Data Analysis", Oxford University Press.
% Silverman, B.W. (1998), "Density Estimation for Statistics and Data
% Analysis", Chapman & Hall/CRC, London. [pp. 48]
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('bwe')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
%
% Examples:
%
%{
    % Bandwidth and kernel density estimates for a univariate normal sample.
    
    % the normal probability density function
    npdf = @(x) (exp(-0.5*x.^2)/sqrt(2*pi));
    % normal kernel density
    nkde = @(x,unidata,h) mean(npdf((x-unidata)/h)/h); % kernel density

    % a univariate normal sample
    unidata = randn(200,1);

    % bandwidth estimation
    h = bwe(unidata);

    % plot of kernel density with estimated bandwidth
    warning('off');
    fplot(@(x) nkde(x,unidata,h),[-10,10],'r')

    % plot of the true density
    hold on
    fplot(@(x) (npdf(x)) ,[-10,10],'k')

    % plot of the data
    plot(unidata,npdf(unidata),'xb')
    warning('on');

    axis manual;
    title(['estimated bandwidth: ' num2str(h) ]);
    legend('estimated density','true density','data');
%}

%{
    % Bandwidth and kernel density estimates for a univariate mixture of two normals.
    % The smoothing is shown for various bandwidth values.

    % the normal probability density function
    npdf = @(x) (exp(-0.5*x.^2)/sqrt(2*pi));
    % normal kernel density
    nkde = @(x,unidata,h) mean(npdf((x-unidata)/h)/h); % kernel density

    % mixture of two univariate normal samples
    unidata = [randn(100,1)-5 ; randn(100,1)+5];

    i=0;
    for smfact = 1:3:7
        i=i+1;
        % bandwidth estimation
        h = bwe(unidata) / smfact;
        subplot(3,1,i);
        % plot of kernel density with estimated bandwidth
        warning('off');
        fplot(@(x) nkde(x,unidata,h),[-10,10],'r')
        % plot of the true density
        hold on;
        fplot(@(x) (npdf(x-5) + npdf(x+5)),[-10,10],'k')
        % plot of the data
        plot(unidata,(npdf(unidata-5) + npdf(unidata+5)),'xb')
        warning('on');
        if i == 1
            xlabel(['bw0 = ' num2str(h) ' (estimated from the data)' ]);
        else
            xlabel(['bw0 / ' num2str(i) ' = ' num2str(h) ]);
        end
    end

%}

%{
    % Bandwidth and kernel density estimates for a bivariate dataset.

    load fishery;
    X = fishery{:,:};
    X = X+10^(-8)*abs(randn(677,2)); % some jittering to avoid dplicate points
    h = bwe(X)
    h = bwe(X,'scott')
    h = bwe(X,'normal')
    h = bwe(X,'robust')
%}

%% Beginning of code

% bandwidth selection
%  Remark: ksdensity uses by default Scott's rule

[n,d] = size(X);
% units must be along the rows.
if n<d
    X = X';
    [n,d] = size(X);
end

maxX = max(X);
minX = min(X);

% Scott's rule is called by default. The 'if' statement is introduced for
% computational efficiency reasons in case bwe has to called many times.
% This is to avoid the string comparison in the 'switch' statement.
if nargin<2 %|| d > 1
    
    %  Scott's rule (optimal for normal distribution, with robust sigma estimate)
    %  Scott p. 164
    sig = mad(X,1) / 0.6745;
    if sig <= 0, sig = maxX - minX; end
    if sig > 0
        bw  = sig * (4/((d+2)*n))^(1/(d+4));
    else
        bw = 1;
    end
    
    if d==1 , bw=1.06 * bw ; end % Univ. Normal rule - Scott p. 144 + p. 154
    
else
    
    switch bwopt
        case 'scott'
            %Scott's rule is the default (Scott p. 164): just call again
            %the function with no bwopt option.
            bw = bwe(X);
            
        case 'normal'
            % Silverman's (1986) rule of thumb, applied for d = 1. It is
            % valid if the underlying density being estimated is Gaussian.
            bw = 1.06  * std(X) * n^(-1/5);
            
        case 'robust'
            % As for Silverman's rule, but in presence of outliers.
            bw = 0.786 * iqr(X) * n^(-1/5);
            
        otherwise
            % It the user's bandwidth option does not exist, bwe is called
            % recursevely without options, so that the default is applied.
            bw = bwe(X);
            
            % Remark: other methods to be considered in the future include
            % Diggle's method and Sheater and Jones(1991) and Wand and
            % Jones (1995) pp.74-75
            
    end
    
end

end

%FScategory:UTISTAT
