function F = kdebiv(X,varargin)
%kdebiv computes (and optionally plots) a kernel smoothing estimate for bivariate data.
%
%<a href="matlab: docsearchFS('kdebiv')">Link to the help function</a>
%
% This function is introduced in FSDA to support MATLAB releases older
% than R2016a, when function ksdensity was only addressing one-dimensional
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
% Optional input arguments:
%
%   contourtype: Plot on the screen. String. Takes one of these strings:
%               - contourtype = 'contour' generates a contour plot.
%               - contourtype = 'contourf' generates a filled contour plot. 
%               - contourtype = 'surf' generates a surf plot.
%               - contourtype = 'mesh' generates a mesh plot.
%               Unless specified otherwise, the colormap of the plots is
%               based on grey levels.
%               Data Types - char
%               Example - 'contourtype','contourf'
%
%   cmap:       Three-column matrix with colormap values in the range
%               [0,1]. Matrix. A personalized colormap is used to plot 
%               the contour.  Each row of 'plots' is an RGB triplet that
%               defines one color.
%                 Data Types - char | double
%                 Example - 'cmap','gray'
%                 Example - 'cmap',[0, 0, 0.3 ; 0, 0, 0.4 ;  0, 0, 0.5 ]
%
%
%  Output:
%
%   F :         F is the vector of density values. Matrix. The estimate is
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
% Copyright 2008-2016.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('kdebiv')">Link to the help function</a>

% Examples:

%{
      %% Density plots for a mixture of two normal distributions.
      X1 = [0+.5*randn(150,1)   5+2.5*randn(150,1)];
      X2 = [1.75+.25*randn(60,1) 8.75+1.25*randn(60,1)];
      X = [X1 ; X2];

      % A filled contour plot obtained using colormap 'cmap' = 'summer'.
      F1 = kdebiv(X,'contourtype','contourf','cmap','summer');
      title('A filled contour plot obtained using colormap ''summer''');

%}

%{
      %% A standard (not filled) contour plot obtained using colormap 'cmap' = 'hot'.
      figure;
      F2 = kdebiv(X,'cmap','hot');
      title('A standard (not filled) contour plot obtained using colormap ''hot''');

      % A filled contour plot with personalized colormap: note the last
      % line of cmap (1 1 1), which is added to obtain a white background
      % in the low densit areas.
      figure;
      % Data points, with associated clickable legends.
      plot(X1(:,1),X1(:,2),'xr' , X2(:,1),X2(:,2),'oc');
      clickableMultiLegend('group 1','group 2');
      % superimpose the contour plot
      hold on;
      cmap =   [0, 0, 0.3 ; 0, 0, 0.4 ;  0, 0, 0.5 ; 0, 0, 0.6 ;  0, 0, 0.8 ; 0, 0, 1.0 ; 1, 1, 1 ];
      F3 = kdebiv(X,'cmap',cmap);
      title('A filled contour plot with personalized colormap and data point superimposed');

      % just to position the figures in "cascade".
      cascade;

%}

%{
      % Just to test surf and mesh plots
      figure;
      F4 = kdebiv(X,'contourtype','surf');
      figure;
      F5 = kdebiv(X,'cmap',summer,'contourtype','surf');
      figure;
      F6 = kdebiv(X,'contourtype','mesh');
      figure;
      F7 = kdebiv(X,'cmap',summer,'contourtype','mesh');
%}

%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[X, nn , d] = chkinputM(X,nnargin,vvarargin);

options     = struct('contourtype','contour','cmap','gray');
UserOptions = varargin(1:2:length(varargin));
if ~isempty(UserOptions) && (length(varargin) ~= 2*length(UserOptions))
    error('FSDA:kdebiv:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
end

if nargin>1
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

contourtype = options.contourtype;
cmap        = options.cmap;

if ~isempty(UserOptions)
    cmaptypes    = {'white','flag','prism','colorcube','lines','pink','copper','bone','gray','winter','autumn','summer','spring','parula' , 'jet' , 'hsv' , 'hot' , 'cool'};
    contourtypes = {'contourf' , 'contour' , 'surf' , 'mesh'};
    
    if  isempty(contourtype) || ~(ischar(contourtype) && max(strcmp(contourtype,contourtypes)))
        contourtype = 'contour';
    end
    if isempty(cmap) || ...
            (~ischar(cmap) && ~((size(cmap,2) == 3 && (min(min(cmap))>=0 && max(max(cmap))<=1)))) ...
             || (~ischar(cmap) && max(strcmp(cmap,cmaptypes)))
        cmap        = 'gray';
    end
    
    plot_contour = 1;

else
    plot_contour = 0;
end

%% kernel smoothing estimation

% rule of thumb for the number of bins
nbins = round(8 * log(nn)/2); 
nbins = [nbins , nbins];

if verLessThan('matlab','9.0')
    % The kernel smoothing function is estimated from scratch if MATLAB
    % release is before R2016a
    
    %{
    % Estimate the bandwidth using Scott's rule (optimal for normal
    % distribution):
    % bw = median(abs(xy-repmat(median(xy),nn,1)))/0.6745*(1/nn)^(1/6);
    % It is a rule of thumb  suggested by Bowman and Azzalini (1997), p.31.
    sig      = mad(X,1,1) / 0.6745;
    bw       = sig * (4/((d+2)*nn))^(1/(d+4));
    
    % bivariate c function
    K  = @(sigma,x,y) exp(-(x^2+y^2)/2/norm(sigma) );
    xx = linspace(min(X(:,1)),max(X(:,1)));
    yy = linspace(min(X(:,2)),max(X(:,2)));
    [dx,dy] = meshgrid(xx,yy);
    weight  = K(bw,dx,dy)./sum(sum(K(bw,dx,dy)));
    Ysmooth = conv2(X,weight,'same');
    
    %}
        
    % Compute a two-dimensional histogram without using histc
    %{
    xy_max   = max(X);
    xy_min   = min(X);
    xy_lim   = [-inf -inf inf inf];
    xy_max   = min([xy_max+3*bw ; xy_lim(3:4)]);
    xy_min   = max([xy_min-3*bw ; xy_lim(1:2)]);
    ed1      = linspace(xy_min(1),xy_max(1),nbins(1)+1);
    ed2      = linspace(xy_min(2),xy_max(2),nbins(1)+1);
                
    xi1  = ed1(1:end-1) + .5*diff(ed1);
    ed1  = [-Inf ed1(2:end-1) Inf];
    
    xi2  = ed2(1:end-1) + .5*diff(ed2);
    ed2  = [-Inf ed2(2:end-1) Inf];
    
    % Cut the xy space into cells and count the number of units in each cell.
    bin    = zeros(nn,2);
    [~,bin(:,2)] = histc(X(:,1),ed1);
    [~,bin(:,1)] = histc(X(:,2),ed2);
    
    H  = accumarray(bin,1,nbins([2 1])) ./ nn;
    %}
    
    %Compute a two-dimensional histogram using hist3
    [H,C] = hist3(X,nbins) ;%./ nn
    
    % position of the bin centers
    xi1 = C{1,1};
    xi2 = C{1,2};
                
    % subfunction smooth1D smoothes the two-dimensional histogram.
    % lambda is a smoothing parameter. Smaller lambda values provide
    % smoother results.
    lambda = 10;
    G  = expsm(H ,nbins(2)/lambda);
    F  = expsm(G',nbins(1)/lambda)';

else
    % The MATLAB ksdensity follows. It is based on:
    %   A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
    %   Techniques for Data Analysis," Oxford University Press.
    [F,xi] = ksdensity(X,'Support','unbounded');
    xi1 = xi(:,1);
    xi2 = xi(:,2);
end

% interpolate the estimated density on the grid xi1 and xi1, using meshgrid
% and griddata.
%[xq,yq,F] = computeGrid(xi1 , xi2 , F);

xx = linspace(min(xi1),max(xi1),nbins(1));
yy = linspace(min(xi2),max(xi2),nbins(2));
% define a data grid
[xq,yq] = meshgrid(xx,yy);
% Interpolate the scattered data on the grid
F = griddata(xi1,xi2,F,xq,yq);
        

%% Now plot the countour

if plot_contour
    
    % For plotting reasons, we do not want zero values
    FF = F;
    FF(FF==0) = min(FF(FF~=0));
    
    mymap = colormap(cmap);
    switch(contourtype)
        case 'surf'
            % undocumented: produces a surface plot
            surf(xq,yq,FF,'EdgeAlpha',0);
        case 'mesh'
            % undocumented: produces a mesh plot
            mesh(xq,yq,FF,'EdgeAlpha',1,'FaceAlpha',0);
            mymap(1,:) = [1 1 1];
        case 'contour'
            contour(xq,yq,1-FF,'Clipping','off');
        case 'contourf'
            contourf(xq,yq,1-FF,'Clipping','off');
    end
    
    colormap(mymap);
    
end

%% subfunctions needed for the density estimate if before R2016a

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
