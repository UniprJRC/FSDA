function F = kdebiv(X,contourtype,cmap)
%kdebiv computes (and optionally plots) a kernel smoothing estimate for bivariate data.
%
%<a href="matlab: docsearchFS('kdebiv')">Link to the help function</a>
%
% This function is introduced in FSDA for use with MATLAB releases older
% than R2016a, when function ksdensity was only supporting one-dimensional
% data. For R2016a and subsequent releases, kdebiv uses ksdensity.
%
%  Required input arguments:
%
%   X: two-column matrix with the bi-variate data sample on which a
%      probability density estimate is computed. Matrix. The density is
%      estimated on a grid of points covering the range of the data,
%      created using MATLAB function meshgrid.
%      Data Types - single | double.
%
%
%   contourtype: Plot on the screen. String. Takes one of these string:
%               - contourtype = 'contour' superimposes to the bivariate
%                 scatterplot a contour plot.
%               - contourtype = 'contourf' superimposes to the bivariate
%                 scatterplot a filled contour plot. The colormap of the
%                 filled contour is based on grey levels.
%                 Data Types - char.
%
%   cmap:       - Three-column matrix of values in the range [0,1]
%                 representing a colormap. Matrix. A personalized colormap
%                 is used to plot the contour. Each row of 'plots' is
%                 an RGB triplet that defines one color.
%                 Data Types - double.
%
% Optional input arguments:
%
%
%  Output:
%
%   F =         F is the vector of density values. Matrix. The estimate is 
%               based on the normal kernel function, using the window
%               parameter (bandwidth) that is a function of the number of
%               points and dimension in X.
%
% See also: ksdensity
%
%
% References:
%
%    A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
%    Techniques for Data Analysis," Oxford University Press.
%
%   
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('kdebiv')">Link to the help function</a>

% Examples:

%{
      %% Density plots for a mixture of two normal distributions.
      X1 = [0+.5*randn(100,1)   5+2.5*randn(100,1)];
      X2 = [1.75+.25*randn(40,1) 8.75+1.25*randn(40,1)];
      X = [X1 ; X2];
    
      % A filled contour plot.
      contourtype = 'contourf';
      F1 = kdebiv(X,contourtype);

      % The default is a standard (not filled) contour plot.
      figure;
      F2 = kdebiv(X,[],'hot');

      % A filled contour plot with personalized colormap: note the last
      % line of cmap (1 1 1), which is added to obtain a white background
      % in the low densit areas.
      figure;
      cmap =   [0, 0, 0.3 ; 0, 0, 0.4 ;  0, 0, 0.5 ; 0, 0, 0.6 ;  0, 0, 0.8 ; 0, 0, 1.0 ; 1, 1, 1 ; 1, 1, 1 ];

      F3 = kdebiv(X,contourtype,cmap);

      % Superimpose data points to the last density plot.
      hold on;
      plot(X1(:,1),X1(:,2),'xr' , X2(:,1),X2(:,2),'oc');

      % just to position the figures in "cascade".
      cascade;
%}

%% check options
switch nargin
    case 1
        plot_contour = 0;
    case 2
        plot_contour = 1;
        cmap         = 'gray';
        if ~isempty(contourtype) && ~(ischar(contourtype) && max(strcmp(contourtype,{'contourf' , 'contour' , 'surf' , 'mesh'})))
            contourtype = 'contour';
        end
    case 3
        plot_contour = 1;
        if isempty(contourtype) || ~(ischar(contourtype) && max(strcmp(contourtype,{'contourf' , 'contour' , 'surf' , 'mesh'})))
            contourtype = 'contour';
        end
        if isempty(cmap) || (~ischar(cmap) && ~((size(cmap,2) == 3 && (min(min(cmap))>=0 && max(max(cmap))<=1))))
            cmap        = 'gray';
        end
end

% get data
xx = X(:,1);
yy = X(:,2);

if verLessThan('matlab','9.0')
    % The kernel smoothing function is estimated here from scratch
    
    % Estimate the bandwidth using Scott's rule (optimal for normal
    % distribution). It is a rule of thumb  suggested by Bowman and
    % Azzalini (1997), p.31.
    xy       = [xx,yy];
    [nn , d] = size(xy);
    sig      = mad(xy,1,1) / 0.6745;
    bw       = sig * (4/((d+2)*nn))^(1/(d+4));
    %bw = median(abs(xy-repmat(median(xy),nn,1)))/0.6745*(1/nn)^(1/6);
    
    % rule of thumb for the number of bins
    nbins = round(8 * log(nn))/2; %nbins = 50;
    
    %             % bivariate Gaussian kernel function
    %             K = @(sigma,x,y) exp(-(x.^2+y.^2)/2/sigma.^2 );
    %             [dx,dy]=meshgrid(xx,yy);
    %             weight = K(bw,dx,dy)/sum(sum(K(bw,dx,dy)));
    %             Ysmooth = conv2(xy,weight,'same');
    
    
    % Compute a two-dimensional histogram
    xy_max   = max(xy);
    xy_min   = min(xy);
    xy_lim   = [-inf -inf inf inf];
    xy_max   = min([xy_max+3*bw ; xy_lim(3:4)]);
    xy_min   = max([xy_min-3*bw ; xy_lim(1:2)]);
    ed1     = linspace(xy_min(1),xy_max(1),nbins+1);
    ed2     = linspace(xy_min(2),xy_max(2),nbins+1);
    
    nbins   = [nbins , nbins];
    
    %%
    % The histograms have 200 bins max in both directions.
    %             minx  = min(xx,[],1);  maxx  = max(xx,[],1);
    %             miny  = min(yy,[],1);  maxy  = max(yy,[],1);
    %             nbins = [min(numel(unique(xx)),200) , min(numel(unique(yy)),200)];
    %
    %             ed1  = linspace(minx, maxx, nbins(1)+1);
    %             ed2  = linspace(miny, maxy, nbins(2)+1);
    %%
    
    xi1  = ed1(1:end-1) + .5*diff(ed1);
    ed1  = [-Inf ed1(2:end-1) Inf];
    
    xi2  = ed2(1:end-1) + .5*diff(ed2);
    ed2  = [-Inf ed2(2:end-1) Inf];
    
    % The xy space is cut into rectangles and the number of
    % observations in each rectangle is counted.
    [nn,~] = size(xx);
    bin    = zeros(nn,2);
    [~,bin(:,2)] = histc(xx,ed1);
    [~,bin(:,1)] = histc(yy,ed2);
    H  = accumarray(bin,1,nbins([2 1])) ./ nn;
    
    % subfunction smooth1D smoothes the two-dimensional histogram.
    % lambda is a smoothing parameter. Smaller lambda values provide
    % smoother results.
    lambda = 10;
    G  = expsm(H ,nbins(2)/lambda);
    F = expsm(G',nbins(1)/lambda)';
    
else
    % The MATLAB ksdensity follows. It is based on:
    %   A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
    %   Techniques for Data Analysis," Oxford University Press.
    [F,xi] = ksdensity([xx yy],'Support','unbounded');
    xi1 = xi(:,1);
    xi2 = xi(:,2);
end

% Call subfunction computeGrid to interpolate the estimated
% density on the grid xi1 and xi1. Uses meshgrid and griddata.
[xq,yq,F] = computeGrid(xi1 , xi2 , F);

if plot_contour
    
    % For plotting reasons, we do not want zero values
    FF = F;
    FF(FF==0) = min(FF(FF~=0));
    
    mymap = colormap(cmap);
    switch(contourtype)
        case 'surf'
            % undocumented: produces a surface plot
            surf(xq,yq,1-FF,'EdgeAlpha',0);
        case 'mesh'
            % undocumented: produces a mesh plot
            mesh(xq,yq,1-FF,'EdgeAlpha',1,'FaceAlpha',0);
            mymap(1,:) = [1 1 1];
        case 'contour'
            contour(xq,yq,1-FF,'Clipping','off');
        case 'contourf'
            contourf(xq,yq,1-FF,'Clipping','off');
    end
    
    colormap(mymap);
    
end

    function [xq,yq,z] = computeGrid(x1,x2,fout)
        % computeGrid is a subfunction used to interpolate function
        % fout on the points in the grid given by x1 and x2
        x = linspace(min(x1),max(x1));
        y = linspace(min(x2),max(x2));
        % define a data grid
        [xq,yq] = meshgrid(x,y);
        orig_state = warning;
        warning('off','all');
        % Interpolate the scattered data on the grid
        z = griddata(x1,x2,fout,xq,yq);
        warning(orig_state);
    end

    function Z  = expsm(GG,lambda)
        % A rough but fast smoothing of the two-dimensional histogram
        [m,~]   = size(GG);
        E       = eye(m);
        D1      = diff(E,1);
        D2      = diff(D1,1);
        P       = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
        Z       = (E + P) \ GG;
    end

end