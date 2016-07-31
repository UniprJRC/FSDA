function [F,Xi,bw] = kdebiv(X,varargin)
%kdebiv computes (and optionally plots) a kernel smoothing estimate for bivariate data.
%
%<a href="matlab: docsearchFS('kdebiv')">Link to the help function</a>
%
% This function is introduced in FSDA to support MATLAB releases older than
% R2016a, when function ksdensity was only addressing one-dimensional data.
% For R2016a and subsequent releases, kdebiv uses ksdensity. Otherwise, the
% function computes a nonparametric estimate of the probability density
% function based on a normal kernel and using a bandwidth estimated as a
% function of the number of points in X.
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
%               Data Types - char | double
%               Example - 'cmap','gray'
%               Example - 'cmap',[0, 0, 0.3 ; 0, 0, 0.4 ;  0, 0, 0.5 ]
%
%   pdfmethod:  Density estimation method. Supported options are 'matlab' 
%               and 'fsda'. 
%               - 'matlab' (default) uses the default approach implemented
%                  in the MATLAB ksdensity function, using a normal kernel. 
%               - 'fsda' computes a nonparametric estimate of the 
%                 probability density function based on a normal kernel 
%                 and using a bandwidth estimated as a function of the 
%                 number of points in X.
%               Independently from the choice of the user, the function
%               switches automatically to 'fsda' in case the user is using
%               releases older than R2016a, when function ksdensity was
%               only addressing one-dimensional data.
%               Data Types - char
%               Example - 'pdfmethod','fsda'
%
%  Output:
%
%   F :         Vector of density values. Matrix. The estimate od F is
%               based on the normal kernel function, using the window
%               parameter (bandwidth) that is a function of the number of
%               points and dimension in X.
%
%   Xi :        Grid of evaluation points. Matrix. 2d matrix of equally-spaced
%               points where the normal kernel function has been evaluated.
%
%   bw :        Bandwidth value. Vector. The bandwidth used for the density
%               estimation.
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
      [F1,Xi,bw] = kdebiv(X,'contourtype','contourf','cmap','summer');
      title('A filled contour plot obtained using colormap ''summer''');
      hold on
      plot(X(:,1),X(:,2),'rx')

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

%{
      %% Test option 'method'.
      [F,Xi,bw]      = kdebiv(X,'pdfmethod','fsda');
      disp(['fsda: estimated bandwidth over x axis is ' num2str(bw(1))]);
      disp(['fsda: estimated bandwidth over y axis is ' num2str(bw(2))]);
      if ~verLessThan('matlab','9.0')
        [F2,Xi2,bw2] = kdebiv(X,'pdfmethod','matlab');
        disp(['matlab (ksdensity): estimated bandwidth over x axis is ' num2str(bw2(1))]);
        disp(['matlab (ksdensity): estimated bandwidth over y axis is ' num2str(bw2(2))]);
      end
%}

%{
    % test option 'method', with plots.
    close all;

    figure
    [F2,Xi2,bw2] = kdebiv(X,'cmap','gray','pdfmethod','fsda');
    hold on
    plot(X(:,1),X(:,2),'rx');
    title('pdfmethod = fsda');

    figure
    [F3,Xi3,bw3] = kdebiv(X,'cmap','gray','pdfmethod','matlab');
    hold on
    plot(X(:,1),X(:,2),'rx')
    title('pdfmethod = matalb');

    figure
    [F4,Xi4,bw4] = kdebiv(X,'cmap','gray','pdfmethod','quick_and_dirty');
    hold on
    plot(X(:,1),X(:,2),'rx')
    title('pdfmethod = quick and dirty (remark: to be fixed)');

    figure
    [F5,Xi5,bw5] = kdebiv(X,'cmap','gray','pdfmethod','independence');
    hold on
    plot(X(:,1),X(:,2),'rx')
    title('pdfmethod = independence');

    cascade;

%}

%{
     % Just a speed test 
     if ~verLessThan('matlab','9.0')
        tt = 0; tt2=0;
        for i = 1 : 20
              X1 = [0+.5*randn(150,1)   5+2.5*randn(150,1)];
              X2 = [1.75+.25*randn(60,1) 8.75+1.25*randn(60,1)];
              X = [X1 ; X2];

              t = tic;
              [F,Xi,bw] = kdebiv(X,'pdfmethod','fsda');
              tt = tt+toc(t);

              t2 = tic;
              [F,Xi,bw] = ksdensity(X);
              tt2 = tt2+toc(t2);
        end
        disp(['kdebiv    time = ' num2str(tt)] );
        disp(['ksdensity time = ' num2str(tt2)] );
    end
%}



%% Input parameters checking
nnargin=nargin;
vvarargin=varargin;
[X, nn , d] = chkinputM(X,nnargin,vvarargin);

if d ~= 2
    error('FSDA:kdebiv:WrongInput','This function applies to bivariate data only!');
end

options     = struct('contourtype','contour','cmap','gray','pdfmethod','matlab');
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
pdfmethod   = options.pdfmethod;

if ~isempty(UserOptions)
    
    % do plot only if the user has selected an optional plotting option
    if  max(strcmp('contourtype',UserOptions)) || max(strcmp('cmap',UserOptions))
        plot_contour = 1;
    else
        plot_contour = 0;
    end
    
    cmaptypes    = {'white','flag','prism','colorcube','lines','pink','copper','bone','gray','winter','autumn','summer','spring','parula' , 'jet' , 'hsv' , 'hot' , 'cool'};
    contourtypes = {'contourf' , 'contour' , 'surf' , 'mesh'};
    pdfmethodtypes  = {'matlab','fsda','quick_and_dirty','independence'};
    
    if  isempty(contourtype) || ~(ischar(contourtype) && max(strcmp(contourtype,contourtypes)))
        contourtype = 'contour';
    end
    if isempty(cmap) || ...
            (~ischar(cmap) && ~((size(cmap,2) == 3 && (min(min(cmap))>=0 && max(max(cmap))<=1)))) ...
            || (~ischar(cmap) && max(strcmp(cmap,cmaptypes)))
        cmap        = 'gray';
    end
    if  isempty(pdfmethod) || ~(ischar(pdfmethod) && max(strcmp(pdfmethod,pdfmethodtypes)))
        pdfmethod = 'matlab';
    end
    
else
    % do not make plots if the user do not selected an optional plotting option
     plot_contour = 0;
end


%% kernel smoothing parameters

% estimate the number of bins (default is Freedman-Diaconis like rule):
nbins = binsnum(X,'sqrt');
nbins = min(nbins , nn/2);
nbins = [nbins , nbins];

%% kernel smoothing estimation.

if verLessThan('matlab','9.0')
    method = 'fsda';
else
    method = pdfmethod;
end

switch method
    
    case 'matlab'
        
        % Recent releases of MATLAB ksdensity (or equivalently mvksdensity)
        % support bivariate density and we therefore use it. The estimate is
        % based on a normal kernel function evaluated at equally-spaced points,
        % xi, that cover the range of the data in X. This is done along the
        % lines of:
        %   A.W. Bowman and A. Azzalini (1997), "Applied Smoothing
        %   Techniques for Data Analysis," Oxford University Press.
        
        %[F,xi,bw] = ksdensity(X,'Support','unbounded');
        [F,xi,bw] = mvksdensity(X);
        xx1 = xi(:,1);
        xx2 = xi(:,2);
        
    case 'fsda'
        
        % nonparametric estimate of the probability density function
        % based on a normal kernel and using a bandwidth estimated as a
        % function of the number of points in X.
        
        % data points
        x1 = X(:,1);
        x2 = X(:,2);
           
        % generate a vector of 100 evenly spaced points between x1 and x2
        xx1 = linspace(min(x1),max(x1));
        xx2 = linspace(min(x2),max(x2));

        m1 = length(xx1);
        m2 = length(xx2);

        % A bandwidth estimate over original data X(:,1) (x direction)
        bw(1) = bwe(xx1);
        % A bandwidth estimate over original data X(:,2) (y direction)
        bw(2) = bwe(xx2);

        % prepare a rectangular grid over xi1 and xi2
        %[ggridx2,ggridx1] = meshgrid(xx2,xx1);
        [ggridx1,ggridx2] = meshgrid(xx1,xx2);
        ggridx1 = repmat(ggridx1, [1,1,nn]);
        ggridx2 = repmat(ggridx2, [1,1,nn]);
        
        % mean estimates over xi1 and xi2
        mu1(1,1,:) = x1;
        mu1 = repmat(mu1,[m1,m2,1]);
        mu2(1,1,:) = x2;
        mu2 = repmat(mu2,[m1,m2,1]);

        % Normal density estimate over the grid
        F = sum(normpdf(ggridx1,mu1,bw(1)) .* normpdf(ggridx2,mu2,bw(2)), 3) / nn;
        
    case 'independence'
        % Calculate combined x-y pdf under assumption of independence
        
        % data points
        x1 = X(:,1);
        x2 = X(:,2);
        
        % indepennt estimates of the density
        [pdfx1 , xx1]= ksdensity(x1);
        [pdfx2 , xx2]= ksdensity(x2);
        
        % Create 2-d grid of function values
        [pdfxx1,pdfxx2] = meshgrid(pdfx1,pdfx2);
        
        % Calculate combined pdf
        F = pdfxx1.*pdfxx2;
        
        bw = [];
        
    case 'quick_and_dirty'
        % The kernel smoothing function is estimated from scratch if MATLAB
        % release is before R2016a
        
        %Compute a two-dimensional histogram using hist3
        [H,C] = hist3(X,nbins) ;%./ nn
        
        % position of the bin centers
        xx1 = C{1,1};
        xx2 = C{1,2};
        
        % subfunction smooth1D smoothes the two-dimensional histogram.
        % lambda is a smoothing parameter. Smaller lambda values provide
        % smoother results.
        lambda = 10;
        G  = expsm(H ,nbins(2)/lambda);
        F  = expsm(G',nbins(1)/lambda)';
        
        bw = [];
end

% equally-spaced points where the normal kernel function has been evaluated
Xi = [xx1(:) , xx2(:)];

%% Now plot the countour

if plot_contour
    
    % interpolate the estimated density on the grid xi1 and xi1, using
    % meshgrid and griddata.

    % control of the axis limits
    xmin = min(X(:,1)); xmax = max(X(:,1));
    ymin = min(X(:,2)); ymax = max(X(:,2));
    deltax = (xmax - xmin) / 10;     
    deltay = (ymax - ymin) / 10;

    % generate a vector of 100 evenly spaced points between x1 and x2
    %xx = linspace(min(xx1),max(xx1));
    %yy = linspace(min(xx2),max(xx2));
    
    % generate a vector of 100 evenly spaced points between the data limits
    xx = linspace(xmin-deltax,xmax+deltax);
    yy = linspace(ymin-deltay,ymax+deltay);
    
    % define a data grid on the evenly spaced points
    [xq,yq] = meshgrid(xx,yy);
    % Interpolate the scattered data on the grid
    F = griddata(xx1,xx2,F,xq,yq);

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
    
    % colormap
    colormap(mymap);
    
end

%% subfunctions needed for the density estimate if before R2016a

    function [nbin] = binsnum(A,rule)
        % binsnum estimates the "optimal" number of bins for a two
        % dimensional histogram or for kernel density estimation.
        % A is the data matrix and 'rule', if given, is the rule to be
        % applied (default is Freedman?Diaconis).
        %{
            rng(123,'twister');
            N = 250;
            A = random('normal',0,1,[N,2]);  % Generate a N-by-2 matrix with N(0,1)
            A(:,2) = A(:,2) * 5;             % Make the second dimension more variable

            [nbins_sqrtN] = binsnum(A,'sqrt');
            [nbins_str]   = binsnum(A,'sturges');
            [nbins_fd_1]  = binsnum(A,'fd');
            [nbins_fd_2]  = binsnum(A,'fdn');
        
            % Plot the results / Make bivariate histograms
            figure(1);
            subplot(2,2,1);
            hist3(A,[ nbins_sqrtN nbins_sqrtN] );
            title('Square Root rule');

            subplot(2,2,2);
            hist3(A,[ nbins_str nbins_str] );
            title('Sturges formula rule');

            subplot(2,2,3);
            hist3(A,[ nbins_fd_1 nbins_fd_1]);
            title('Freedman?Diaconis-like rule');

            subplot(2,2,4);
            hist3(A,[ nbins_fd_2 nbins_fd_2]);
            title('Freedman?Diaconis rule on the norms');
        %}
        
        
        if nargin < 2 || (nargin == 2 && (ischar(rule) && sum(strcmp(rule,{'sqrt','sturges','fd','fdn'}))))
            rule = 'sqrt';
        end
        
        N = size(A,1);
        
        switch rule
            case 'sqrt'
                % The sqrt(N) rule:
                nbin = floor(sqrt(N));
                
            case 'sturges'
                % The Sturges formula:
                nbin = ceil(log2(N) +1);
                
            case 'fd'
                % The Freedman?Diaconis-like choice:
                IQRs = iqr(A);              % Get the IQ ranges across each dimension
                Hs = 2* IQRs* N^(-1/3);     % Get the bandwidths across each dimension
                Ranges = range(A);          % Get the range of values across each dimension
                % Get the suggested number of bins along each dimension
                nbins_dim1 = ceil(Ranges(1)/Hs(1)); % 12 here
                nbins_dim2 = ceil(Ranges(2)/Hs(2)); % 15 here
                % Get the maximum of the two
                nbin = max( [nbins_dim1, nbins_dim2]);
                
            case 'fdn'
                % The Freedman-Diaconis choice on the norms
                Norms   = sqrt(sum(A.^2,2));            % Get the norm of each point in th 2-D sample
                H_norms = 2* iqr(Norms)* N^(-1/3);      % Get the "norm" bandwidth
                nbin    = ceil(range(Norms)/ H_norms);  % Get number of bins
                
            case 'tobeverified'
                % rule of thumb for the number of bins
                nbin = round(8 * log(N)/2);
                
            otherwise
                % If nothing is provided, the sqrt(N) rule is used:
                nbin = floor(sqrt(N));
                
        end
        
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
