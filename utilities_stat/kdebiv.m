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
%  XI:          Evaluation points of the estimated density.
%               Matrix. In this case the density is estimated using X and evaluated on XI.
%               Data Types - single | double.
%               Example - 'XI',X
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
%   F :         Density values. Vector. The estimate of F is based on the
%               normal kernel function, using the window parameter
%               (bandwidth) that is a function of the number of points and
%               dimension in X.
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
      % A standard (not filled) contour plot obtained using colormap 'cmap' = 'hot'.
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
      close all;

      % Just to test surf and mesh plots.
      figure;
      F4 = kdebiv(X,'contourtype','surf');
      figure;
      F5 = kdebiv(X,'cmap',summer,'contourtype','surf');
      figure;
      F6 = kdebiv(X,'contourtype','mesh');
      figure;
      F7 = kdebiv(X,'cmap',summer,'contourtype','mesh');

      cascade;
%}

%{
      % Test option 'method'.
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

    if ~verLessThan('matlab','9.0')

        figure
        [F3,Xi3,bw3] = kdebiv(X,'cmap','gray','pdfmethod','matlab');
        hold on
        plot(X(:,1),X(:,2),'rx')
        title('pdfmethod = matalb');

        figure
        [F4,Xi4,bw4] = kdebiv(X,'cmap','gray','pdfmethod','histsmooth');
        hold on
        plot(X(:,1),X(:,2),'rx')
        title('pdfmethod = histogram smoothing (remark: to be fixed)');

        figure
        [F5,Xi5,bw5] = kdebiv(X,'cmap','gray','pdfmethod','independence');
        hold on
        plot(X(:,1),X(:,2),'rx')
        title('pdfmethod = independence');

    else

        disp('For this MATLAB release, only ''fsda'' option can be used' );

    end

    cascade;

%}

%{
     % Just a speed test
     if ~verLessThan('matlab','9.0')
        tt = 0; tt2=0;
        for i = 1 : 100
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

options     = struct('XI',[],'contourtype','contour','cmap','gray','pdfmethod','matlab');
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
XI          = options.XI;

if ~isempty(UserOptions)
    
    % do plot only if the user has selected an optional plotting option
    if  max(strcmp('contourtype',UserOptions)) || max(strcmp('cmap',UserOptions))
        plot_contour = 1;
    else
        plot_contour = 0;
    end
    
    cmaptypes    = {'white','flag','prism','colorcube','lines','pink','copper','bone','gray','winter','autumn','summer','spring','parula' , 'jet' , 'hsv' , 'hot' , 'cool'};
    contourtypes = {'contourf' , 'contour' , 'surf' , 'mesh'};
    pdfmethodtypes  = {'matlab','fsda','histsmooth','independence'};
    
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

%% kernel smoothing estimation.

if verLessThan('matlab','9.0') %&& strcmp(pdfmethod,'matlab')
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
        
        if isempty(XI)
            [F,xi,bw] = mvksdensity(X);
        else
            [F,xi,bw] = mvksdensity(X,XI);
        end
        xi1 = xi(:,1);
        xi2 = xi(:,2);
        
    case 'fsda'
        
        % nonparametric estimate of the probability density function
        % based on a normal kernel and using a bandwidth estimated as a
        % function of the number of points in X.
        
        % data points
        X1 = X(:,1);
        X2 = X(:,2);
        
        % Generate a vector of m evenly spaced points between x1 and x2.
        % To be coherent with MATLAB (ksdensity/mvksdensity) we set m=30.
        m  = 30;
        x1 = linspace(min(X1),max(X1),m);
        x2 = linspace(min(X2),max(X2),m);
        
        % Prepare a rectangular grid over x1 and x2. Note that gx1 and gx2
        % are the equally-spaced points where the normal kernel function
        % has been evaluated
        [gx1 , gx2] = meshgrid(x1,x2);
        ggx1 = repmat(gx1, [1,1,nn]);
        ggx2 = repmat(gx2, [1,1,nn]);
        
        % Normal density means are the original points themselves. To
        % prepare the case (in future releases) of different grid sizes in
        % the two dimensions, we use m1 and m2 instead of just m.
        m1 = m; m2 = m;
        mu1(1,1,:) = X1;
        mu1 = repmat(mu1,[m1,m2,1]);
        mu2(1,1,:) = X2;
        mu2 = repmat(mu2,[m1,m2,1]);
        
        % Normal density standard deviations are given by the bandwidths,
        % which are estimated over the two directions of the original data,
        % X(:,1) and X(:,2)
        bw(1) = bwe(X1);
        bw(2) = bwe(X2);
        
        % Normal density estimate over the grid
        F = sum(normpdf(ggx1,mu1,bw(1)) .* normpdf(ggx2,mu2,bw(2)), 3) / nn;
        
        % Values to be returned by kdebiv
        xi1 = gx1(:);
        xi2 = gx2(:);
        F = F(:);
        
        if ~isempty(XI)
            if verLessThan('matlab', '8.1')
                Fpdfe = TriScatteredInterp(xi1,xi2,F); %#ok<DTRIINT>
            else
                Fpdfe = scatteredInterpolant(xi1,xi2,F);
            end
            F  = Fpdfe(XI(:,1),XI(:,2)); 
        end
        
        
    case 'independence'
        % Calculate combined x-y pdf under assumption of independence
        
        % data points
        X1 = X(:,1);
        X2 = X(:,2);
        
        % indepennt estimates of the density
        [pdfx1 , xxi1]= ksdensity(X1);
        [pdfx2 , xxi2]= ksdensity(X2);
        
        % Create 2-d grid of function values
        [pdfxx1,pdfxx2] = meshgrid(pdfx1,pdfx2);
        
        % Calculate combined pdf
        F = pdfxx1.*pdfxx2;
        
        % Values to be returned by kdebiv
        F = F(:);
        [xi1 , xi2] = meshgrid(xxi1,xxi2);
        xi1 = xi1(:);
        xi2 = xi2(:);
        
        bw = [];
        
    case 'histsmooth'
        % A histogram smoothing method which does not make use of a model density estimate
        
        % Estimate the number of bins (default is Freedman-Diaconis like rule):
        nbins = binsnum(X,'sqrt');
        nbins = min(nbins , nn/2);
        nbins = [nbins , nbins];
        
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
        
        bw = [];
end

% equally-spaced points where the normal kernel function has been evaluated
Xi = [xi1(:) , xi2(:)];

%% Now plot the countour

if plot_contour
    
    % interpolate the estimated density on the grid xi1 and xi1, using
    % meshgrid and griddata.
    
    % control of the axis limits
    xmin = min(X(:,1)); xmax = max(X(:,1));
    ymin = min(X(:,2)); ymax = max(X(:,2));
    deltax = 0;%(xmax - xmin) / 10;
    deltay = 0;%(ymax - ymin) / 10;
    
    % generate a vector of 100 evenly spaced points between the data limits
    xx = linspace(xmin-deltax,xmax+deltax);
    yy = linspace(ymin-deltay,ymax+deltay);
    
    % define a data grid on the evenly spaced points
    [xq,yq] = meshgrid(xx,yy);
    % Interpolate the scattered data on the grid
    if isempty(XI)
        FF = griddata(xi1,xi2,F,xq,yq);
    else
        FF = griddata(XI(:,1),XI(:,2),F,xq,yq);
    end
    
    % For plotting reasons, we do not want zero values
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

%% Subfunctions

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
                % The Freedman-Diaconis-like choice:
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
        E       = eye(size(GG));
        D1      = diff(E,1);
        D2      = diff(D1,1);
        P       = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
        Z       = (E + P) \ GG;
    end

    % case 'fsdaori'
    %
    %     % nonparametric estimate of the probability density function
    %     % based on a normal kernel and using a bandwidth estimated as a
    %     % function of the number of points in X.
    %
    %     % data points
    %     X1 = X(:,1);
    %     X2 = X(:,2);
    %
    %     % generate a vector of 100 evenly spaced points between x1 and x2
    %     m   = 100;
    %     x1 = linspace(min(X1),max(X1),m);
    %     x2 = linspace(min(X2),max(X2),m);
    %
    %     m1 = length(x1);
    %     m2 = length(x2);
    %
    %     % A bandwidth estimate over original data X(:,1) (x direction)
    %     bw(1) = bwe(x1);
    %     % A bandwidth estimate over original data X(:,2) (y direction)
    %     bw(2) = bwe(x2);
    %
    %     % prepare a rectangular grid over xi1 and xi2
    %     %[ggridx2,ggridx1] = meshgrid(xx2,xx1);
    %     [gx1,gx2] = meshgrid(x1,x2);
    %     gx1 = repmat(gx1, [1,1,nn]);
    %     gx2 = repmat(gx2, [1,1,nn]);
    %
    %     % mean estimates over xi1 and xi2
    %     mu1(1,1,:) = X1;
    %     mu1 = repmat(mu1,[m1,m2,1]);
    %     mu2(1,1,:) = X2;
    %     mu2 = repmat(mu2,[m1,m2,1]);
    %
    %     % Normal density estimate over the grid
    %     F = sum(normpdf(gx1,mu1,bw(1)) .* normpdf(gx2,mu2,bw(2)), 3) / nn;

end

%FScategory:UTISTAT
