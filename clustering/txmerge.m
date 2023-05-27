function out = txmerge(Y, k, g, varargin)
% txmerge performs a (hierarchical) merging of the inflated number of 
% components found by tkmeans or TCLUST
%
% <a href="matlab: docsearchFS('txmerge')">Link to the help function</a>
%
% The function txmerge performs either a hierarchical merging of the k 
% components found by tkmeans/TCLUST into g groups, or if g is a decimal 
% number between 0 and 1 and DEMP is used as a distance it performs 
% the merging phase according to such threshold.
%
% Required input arguments:
%
% Y :           Input data. Matrix. n x v data matrix. n observations and v
%               variables. Rows of Y represent observations, and columns represent
%               variables. Missing values (NaN's) and infinite values
%               (Inf's) are allowed, since observations (rows) with missing
%               or infinite values will automatically be excluded from the
%               computations.
%               Data Types - single | double
%
% k:            Number of components searched by tkmeans/TCLUST algorithms.
%               Integer scalar.
%               Data Types - single | double
%
% g:            Merging rule. Scalar. Number of groups retained by the
%               hierarchical agglomeration phase, or threshold of the 
%               pairwise overlap values (i.e. omegaStar) if 0<g<1 and dist=1.
%               Data Types - single | double
%
%
% Optional input arguments:
%
% tkm:          Using tkmerge or tcmerge. Scalar. Scalar. tkm=1 (default) 
%               relies on tkmeans to find an inflated number of clusters,
%               tkm=0 is used to rely on TCLUST instead.
%               Example - 'tkm',0
%               Data Types - double
%
% dist:         Distance between clusters. Scalar, char. Its value 
%               indicates the merging rule for the initial number of 
%               (inflated) clusters. If dist=0 the distance between centroids 
%               is Euclidean (default), if dist=1 directly estimated
%               misclassification probabilities  (DEMP) are used, else use
%               a character according to MATLAB pdist function.
%               Example - 'dist', 'squaredeuclidean'
%               Data Types - single | double | string
%
% alpha:        Global trimming level. Scalar. alpha is a scalar between 0
%               and 0.5. If alpha=0 (default) tkmeans reduces to kmeans and
%               TCLUST reduces to MCLUST.
%               Example - 'alpha', 0.05
%               Data Types - single | double
%
% linkagearg:   Linkage used for hierarchical agglomeration. Single linkage 
%               is the default, see the MATLAB linkage function for other options.
%               Example - 'linkagearg', 'weights'
%
% auto:         Automatic trimming level detection. Scalar. It is set to 1 to 
%               overwrite the prespecified alpha parameter, or it is equal to 0
%               to use alpha as trimming level (default). 
%               Example - 'auto',1
%               Data Types - double
%
% plots    :    Plot on the screen. Scalar, char, or struct.
%               - If plots=0 (default) no plot is produced.
%               - If plots=1, the components merged are shown using the
%               spmplot function. In particular:
%                   * for v=1, an histogram of the univariate data.
%                   * for v=2, a bivariate scatterplot.
%                   * for v>2, a scatterplot matrix.
%               When v>=2 plots offers the following additional features
%               (for v=1 the behaviour is forced to be as for plots=1):
%               - plots='contourf' adds in the background of the bivariate
%                 scatterplots a filled contour plot. The colormap of the
%                 filled contour is based on grey levels as default.
%                 This argument may also be inserted in a field named 'type'
%                 of a structure. In the latter case it is possible to
%                 specify the additional field 'cmap', which changes the
%                 default colors of the color map used. The field 'cmap'
%                 may be a three-column matrix of values in the range [0,1]
%                 where each row is an RGB triplet that defines one color.
%                 Check the colormap function for additional informations.
%               - plots='contour' adds in the background of the bivariate
%                 scatterplots a contour plot. The colormap of the contour
%                 is based on grey levels as default. This argument may
%                 also be inserted in a field named 'type' of a structure.
%                 In the latter case it is possible to specify the additional
%                 field 'cmap', which changes the default colors of the
%                 color map used. The field 'cmap' may be a three-column
%                 matrix of values in the range [0,1] where each row is an
%                 RGB triplet that defines one color.
%                 Check the colormap function for additional informations.
%               - plots='ellipse' superimposes confidence ellipses to
%                 each group in the bivariate scatterplots. The size of the
%                 ellipse is chi2inv(0.95,2), i.e. the confidence level used
%                 by default is 95%. This argument may also be inserted in
%                 a field named 'type' of a structure. In the latter case it
%                 is possible to specify the additional field 'conflev',
%                 which specifies the confidence level to use and it is a
%                 value between 0 and 1.
%               - plots='boxplotb' superimposes on the bivariate scatterplots
%                 the bivariate boxplots for each group, using the boxplotb
%                 function. This argument may also be inserted in a field
%                 named 'type' of a structure.
%               REMARK - The labels<=0 are automatically excluded from the
%                        overlaying phase, considering them as outliers.
%                   Example - 'plots', 1
%                   Data Types - single | double | string
%
%   txOpt:      tkmeans/TCLUST optional arguments. Structure. Empty structure 
%               (default) or structure containing optional input arguments 
%               for tkmeans or TCLUST. See tkmeans and tclust functions.
%               Example - 'txOpt.reftol', 0.0001
%               Data Types - struct
%
%   txOut:      Saving tkmeans/TCLUST output structure. Scalar. It is set 
%               to 1 to save the output structure of tkmeans/TCLUUST into 
%               the output structure of txmerge. Default is 0, i.e. no saving.
%               Example - 'txOut', 1
%               Data Types - single | double
%
% Ysave:        Saving Y. Scalar. Scalar that is set to 1 to request that
%               the input matrix Y is saved into the output structure out.
%               Default is 0, i.e. no saving is done.
%               Example - 'Ysave',1
%               Data Types - double
%
% Output:
%
%    out:   structure which contains the following fields
%
%           out.PairOver 	= Distance matrix among the k components
%                             found by tkmeans/TCLUST.
%
%           out.mergID      = Label for each unit. It is a vector with n
%                             elements which assigns each unit to one of
%                             the groups obtained according to the merging
%                             algorithm applied.
%                             REMARK - out.mergID=0 denotes trimmed units.
%
%           out.txOut       = Output from tkmeans function. The structure
%                             is present if option txOut is set to 1.
%
%           out.Y           = Original data matrix Y. This field is present
%                             only if option Ysave is set to 1.
%
%
%
% See also: dempk, tkmeans, clusterdata, tclusteda, overlapmap
%
%
% References:
%
%   Insolia, L., Perrotta, D. (2023). Tk-Merge: Computationally Efficient 
%   Robust Clustering Under General Assumptions. Advances in Intelligent 
%   Systems and Computing, vol 1433. Springer, Cham. 
%
% Copyright 2008-2021.
% Written by FSDA team
%
% <a href="matlab: docsearchFS('txmerge')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example using txmerge with euclidean distances.
    close all
    % Specify k cluster in v dimensions with n obs
    k = 10;
    v = 2;
    n = 5000;
    % Generate homogeneous and spherical clusters
    rng(100, 'twister');
    out = MixSim(k, v, 'sph', true, 'hom', true, 'int', [0 10], 'Display', 'off', 'BarOmega', 0.05, 'Display','off');
    % Simulating data
    [X, id] = simdataset(n, out.Pi, out.Mu, out.S);
    % Plotting data
    gscatter(X(:,1), X(:,2), id);
    str = sprintf('Simulated data with %d groups in %d dimensions and %d units', k, v, n);
    title(str,'Interpreter','Latex');
    % merging algorithm based on hierarchical clustering
    g = 3;
    txsol = txmerge(X, k*5, g, 'dist', 1, 'plots', 'contourf');
%}

%{
    %% Example using txmerge with additional arguments in the call to tkmeans.
    close all
    % Specify k cluster in v dimensions with n obs
    g = 3;
    v = 2;
    n = 5000;
    % null trimming and noise level
    alpha0 = 0;
    % restriction factor
    restr = 30;
    % Maximum overlap
    maxOm = 0.005;
    % Generate heterogeneous and elliptical clusters
    rng(500, 'twister');
    out = MixSim(g, v, 'sph', false, 'restrfactor', restr, 'int', [0 10], ...
        'Display', 'off', 'MaxOmega', maxOm, 'Display','off');
    % Simulating data
    [X, id] = simdataset(n, out.Pi, out.Mu, out.S);
    % Plotting data
    gg = gscatter(X(:,1), X(:,2), id);
    str = sprintf('Simulated data with %d groups in %d dimensions and %d \nunits, with restriction factor %d and maximum overlap %.2f', ...
        g, v, n, restr, maxOm);
    title(str,'Interpreter','Latex', 'fontsize', 12);
    set(findobj(gg), 'MarkerSize',10);
    legend1 = legend(gca,'show');
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',14, 'Location', 'northwest')
    % number of components searched by tkmeans
    k = g * 6;
    % additional input for tkmeans
    txOpt = struct;
    txOpt.reftol = 0.0001;
    txOpt.msg = 1;
    tkmplots = struct;
    tkmplots.type = 'contourf';
    tkmplots.cmap = [0.3 0.2 0.4; 0.4 0.5 0.5; 0.1 0.7 0.9; 0.5 0.3 0.8; 1 1 1];
    txOpt.plots = tkmplots;
    txOpt.nomes = 0;
    % saving tkmeans output
    txOut = 1;
    txsol = txmerge(X, k, g, 'txOpt', txOpt, 'plots', 'ellipse');
    cascade;
%}

%{
    %% Example using txmerge based on TCLUST and 'weights' linkage
    close all
    % Specify k cluster in v dimensions with n obs
    g = 3;
    v = 2;
    n = 5000;
    % null trimming and noise level
    alpha0 = 0;
    % restriction factor
    restr = 30;
    % Maximum overlap
    maxOm = 0.005;
    % Generate heterogeneous and elliptical clusters
    rng(500, 'twister');
    out = MixSim(g, v, 'sph', false, 'restrfactor', restr, 'int', [0 10], ...
        'Display', 'off', 'MaxOmega', maxOm, 'Display','off');
    % Simulating data
    [X, id] = simdataset(n, out.Pi, out.Mu, out.S);
    % Plotting data
    gg = gscatter(X(:,1), X(:,2), id);
    str = sprintf('Simulated data with %d groups in %d dimensions and %d \nunits, with restriction factor %d and maximum overlap %.2f', ...
        g, v, n, restr, maxOm);
    title(str,'Interpreter','Latex', 'fontsize', 12);
    set(findobj(gg), 'MarkerSize',10);
    legend1 = legend(gca,'Group 1','Group 2','Group 3');
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',12, 'Location', 'northwest')
    % number of components searched by tkmeans
    k = g * 3;
    % additional input for clusterdata (i.e. hierOpt)
    linkagearg = 'weights';
    txsol = txmerge(X, k, g, 'tkm', 1,'linkagearg', linkagearg, 'plots', 'ellipse');
    cascade;
%}

%{
    %% Example using txmerge Euclidean distances or DEMP in the presence of contamination.
    close all
    % Specify k cluster in v dimensions with n obs
    g = 3;
    v = 2;
    n = 5000;
    % 10 percent trimming and uniform noise
    alpha = 0.1;
    noise = alpha*n;
    % restriction factor
    restr = 30;
    % Maximum overlap
    maxOm = 0.005;
    % Generate heterogeneous and elliptical clusters
    rng(500, 'twister');
    out = MixSim(g, v, 'sph', false, 'restrfactor', restr, 'int', [0 10], ...
        'Display', 'off', 'MaxOmega', maxOm, 'Display','off');
    % Simulating data
    [X,id] = simdataset(n, out.Pi, out.Mu, out.S, 'noiseunits', noise);
    % Plotting data
    gg = gscatter(X(:,1), X(:,2), id);
    str = sprintf('Simulating %d groups in %d dimensions and %d units with %d%s \nuniform noise, setting a restriction factor %d and maximum overlap %.2f', ...
        g, v, n, alpha*100, '\%', restr, maxOm);
    title(str,'Interpreter','Latex', 'fontsize', 10);
    set(findobj(gg), 'MarkerSize',10);
    legend1 = legend(gca,'Outliers','Group 1','Group 2','Group 3');
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',12, 'Location', 'northwest')
    % fixing the number of components searched by tkmeans
    k = g * 6;
    % txmerge with hierarchical merging and trimming equal to the level of noise
    txsol1 = txmerge(X, k, g, 'alpha', alpha, 'plots', 'contourf');
    % txmerge using a cutoff g to detect the clusters based on DEMP
    g = 0.05;
    txsol2 = txmerge(X, k, g, 'alpha', alpha, 'dist', 1', 'plots', 'contourf');
    cascade;
%}

%{
    %% Example using txmerge on the M5 dataset using different strategies.
    close all
    Y = load('M5data.txt');
    id = Y(:,3);
    Y = Y(:, 1:2);
    g = max(id);
    n = length(Y);
    noise = length(Y(id==0, 1));
    v = 2; % dimensions
    id(id==0) = -1; % changing noise label
    gg = gscatter(Y(:,1), Y(:,2), id);
    str = sprintf('M5 data set with %d groups in %d dimensions and \n%d units where %d%s of them are noise', g, v, n, noise/n*100, '\%');
    title(str,'Interpreter','Latex', 'fontsize', 12);
    set(findobj(gg), 'MarkerSize',12);
    legend1 = legend(gca,'Outliers','Group 1','Group 2','Group 3');
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',12, 'Location', 'northwest')
    % number of components to search
    k = g*5;
    % null trimming and noise level
    alpha0 = 0;
    % mimimum overlap cut-off value between pair of merged components
    txsol1= txmerge(Y, k, g, 'alpha', alpha0, 'txOut', 1, 'plots', 1);
    % setting alpha equal to noise level (usually not effective here)
    alpha = noise/n;
    txsol2= txmerge(Y, k, g, 'alpha', alpha, 'txOut', 1, 'plots', 1);
    % setting alpha greater than the noise level 
    txsol3 = txmerge(Y, k, g, 'alpha', alpha+0.04, 'txOut', 1, 'plots', 1);
    % using DEMP instead (usually effective)
    txsol3 = txmerge(Y, k, g, 'alpha', alpha+0.04, 'txOut', 1, 'dist', 1, 'plots', 1);
    cascade;
%}

%% Beginning of code 

% Input parameters checking
nnargin = nargin;
vvarargin = varargin;
Y = chkinputM(Y,nnargin,vvarargin);
[~, v] = size(Y);

%% Deleting possible NaN and Inf values
% noNaN = ~sum(isnan(Y), 2);
% Y = Y(noNaN, :);
%
% noInf = ~sum(isinf(Y), 2);
% Y = Y(noInf, :);

%% Optional parameters

% txmerge optional input values
options = struct('tkm', 1, 'dist', 0, 'alpha', 0, 'linkagearg', 'single', ...
    'auto', 0, 'plots', 0, 'txOut', 0, 'txOpt', struct(), 'Ysave', 0);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    
    % Check if all the specified optional arguments were present
    % in structure options
    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
    end
    
end

% use optional arguments specified by the user
if nargin > 3
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    % And check if the optional user parameters are reasonable.
    if options.alpha<0 || options.alpha>0.5
        error('FSDA:WrongAlpha','alpha must be a scalar in the interval [0 0.5]')
    end
    %     if ~(options.plots == 0 || options.plots == 1)
    %         error('FSDA:PlotsOptions','Plot options must be a scalar equal to 0 or 1')
    %     end
    if ~(options.Ysave == 0 || options.Ysave == 1)
        error('FSDA:Ysave','Ysave must be a scalar equal to 0 or 1')
    end
    if ~(options.txOut == 0 || options.txOut == 1)
        error('FSDA:txOut','txOut must be a scalar equal to 0 or 1')
    end
    if ~isstruct(options.txOpt)
        error('FSDA:tkmeansCall','txOpt must be a structure')
    end
end

% Default values for the optional parameters are set inside structure 'options'
tkm = options.tkm;
dist = options.dist;
alpha = options.alpha;
linkagearg = options.linkagearg;
auto = options.auto;
plots = options.plots;
txOut = options.txOut;
txOpt = options.txOpt;
Ysave = options.Ysave;


%% Calling tkmeans (with default or optional arguments)
if tkm == 1
    
    if auto == 1
        % restrfactor.pars = 'EEE';
        restrfactor = 1;
        outEDA = tclusteda(Y,k,0.4:-0.05:0,restrfactor,'msg',0,'plots',1);
        figure
        plot(outEDA.Amon(:,1), outEDA.Amon(:,2))
        set(gca,'XDir','reverse');
        xlabel('\alpha')
        ylabel('ARI index')
        set(gca,'FontSize',16)
        [~, mloc] =  max(outEDA.Amon(:,2));
        alpha = outEDA.Amon(mloc, 1);
    end
    
    if isempty(fieldnames(txOpt))
        % default
        clu = tkmeans(Y, k, alpha);
    else
        % add optional arguments to tkmeans
        txOptNames = fieldnames(txOpt);
        txOptCell = cell(1,2*length(txOptNames));
        for i = 1:numel(txOptNames)
            first=(i-1)*2+1;
            txOptCell{first}=txOptNames{i};
            txOptCell{first+1} = txOpt.(txOptNames{i});
        end
        clu = tkmeans(Y, k, alpha, txOptCell{:});
    end

    if isstruct(clu.emp)
        % change the output used. Use the empirical values when the clustering
        % algorithm did not converge
        clu.muopt = clu.emp.muemp;
        clu.sigmaopt =  clu.emp.sigmaemp;
        clu.siz = clu.emp.sizemp;
        clu.idx = clu.emp.idxemp;
    end

    % index non-empty clusters 
    % eliminate any empty component or outliers
    clu.siz = clu.siz(clu.siz(:,1)>0,:);   % exclude outliers and NaN in clu.siz
    % Eliminate any empty component (NaN)
    if sum(isnan(clu.muopt(:,1)))>0
        indd = ~isnan(clu.muopt(:,1));
        clu.muopt = clu.muopt(indd,:);
        clu.sigmaopt = clu.sigmaopt(:,:,indd);
    end
    
%% Calling TCLUST (with default or optional arguments)
elseif tkm == 0
    restrfactor = 64;
    if auto == 1
        outEDA = tclusteda(Y,k,0.4:-0.05:0,restrfactor,'msg',0,'plots',1);
        figure
        plot(outEDA.Amon(:,1), outEDA.Amon(:,2))
        set(gca,'XDir','reverse');
        xlabel('\alpha')
        ylabel('ARI index')
        set(gca,'FontSize',16)
        [~, mloc] =  max(outEDA.Amon(:,2));
        alpha = outEDA.Amon(mloc, 1);
    end
    
    if isempty(fieldnames(txOpt))
        % default
        clu = tclust(Y, k, alpha, restrfactor);
    else
        % add optional arguments to TCLUST
        txOptNames = fieldnames(txOpt);
        txOptCell = cell(1,2*length(txOptNames));
        for i = 1:numel(txOptNames)
            first=(i-1)*2+1;
            txOptCell{first}=txOptNames{i};
            txOptCell{first+1} = txOpt.(txOptNames{i});
        end
        clu = tclust(Y, k, alpha, restrfactor, txOptCell{:});
    end
else
    error('FSDA:tkmCall','tkm must be either 1 (tkmerge) of 0 (TCLUST)')
end
        
%% DEMP
if dist == 1
    % Empirical values obtained by tkmeans
    [clu.OmegaMap, ~, ~, ~, ~] = overlap(k, v, clu.siz(:,3), clu.muopt, clu.sigmaopt);

    % Sum of all pairs of misclassification probaility to obtain an overlap matrix
    overM = triu(clu.OmegaMap,1)+(tril(clu.OmegaMap,-1))';
    overMsave = overM; % to save it

    % DEMP-K algorithm with hierarchical merging (if g is an integer scalar)
    if g >= 1 && mod(g, 1)==0 && g<=k
        % transform the triangular overlaps matrix in a similarities vector
        candInt = triu(overM)'+ overM; % obtain a symmetrical overlap matrix
        candVec = squareform(candInt); % vector form

        % trasform the similarity measure (i.e. overlap) in a dissimilarity one
        % (i.e. similar to a distance)
        candVec = 1-candVec; % distance vector
    end

%% centroids' distances based on other measures
elseif dist == 0
    candVec = pdist(clu.muopt, 'euclidean');
else
    candVec = pdist(clu.muopt, dist);
end
                                    
    if g >= 1 && mod(g, 1)==0 && g<=k
        % hieararchical clustering
        if ~any(strcmp(linkagearg, {'single', 'nearest'; ...
                'complete', 'farthest'; ...
                'average',  'upgma'; ...
                'weighted', 'wpgma'; ...
                'centroid', 'upgmc'; ...
                'median',   'wpgmc'; ...
                'ward''s',  'incremental'}))
            linkagearg = 'single';
        end
        
        % Construct clusters from the overlap "distance" vector. It follows the
        % MATLAB function clusterdata (see the MATLAB help page for additional
        % informations).
        %Hier = clusterdataFS(candVec, g, linkagearg);
        Z    = linkage_txmerge(candVec,linkagearg);
        Hier = cluster(Z,'maxclust',g);
        
        % initialize mergID with components found by tkmeans
        mergID = clu.idx;
        % update mergID with the new merged labels
        for i = 1:max(Hier)
            ids_i = find(Hier==i);
            merges_i = ismember(mergID, ids_i);
            mergID(merges_i) = i+10000; % trivial strategy to get unique ID [To Enhance]
        end
        
        % restore proper ID values
        mergID(mergID~=0) = mergID(mergID~=0)-10000;
    end
    
%% merging components using a cut-off if 0<g<1 (i.e. omegastar)
if 0<g && g<1
    
    omegaStar = g;
    
    % Merging phase
    % set to nan elements below the main diagonal of overM
    d1 = size(overM, 1);
    In = logical(tril(ones(d1)));
    overM(In) = nan;
    
    % find max pairwise overlap in overM
    mas = max(overM(:));
    [r,c] = find(overM==mas); % problem: if some entries have exactly the same max values
    
    % initialize mergMat (merging matrix)
    % if overM has some 0 these rows are not considered in mergMat
    mergMat = zeros(k*(k-1)/2 - sum(overM(:)==0), 3);
    
    % save the max pairwise overlap value and its indexes,
    mergMat(1,1) = overM(r, c);
    mergMat(1,2) = r;
    mergMat(1,3) = c;
    % exclude it setting it equal to NaN
    overM(r,c) = nan;
    
    % initialize cand0: it will be filled with groups closer to some component
    % already included in mergMat(:, 2:3)
    cand0 = nan(size(overM));
    cand0(r,:) = overM(r,:);
    cand0(:,c) = overM(:,c);
    cand0(:,r) = overM(:,r);
    cand0(c,:) = overM(c,:);
    
    % starting index value which goes row-wise in mergMat after the max is
    % found and placed in the 1st row
    i = 2;
    
    % initialize totRC (total rows and columns), used to index greatest pairwise overlap
    totRC = 1:length(overM);
    
    % overlap values sorted and saved according to nearest groups already
    % found (i.e. with greater overlap)
    while ~all(isnan(overM(:)) | overM(:) == 0) 
        % stopping rule: cand matrix has been totally evaluated or the
        % overlap is 0
        
        if any(cand0(:))~=0 % it could happen that cand0 remains only with zeros
            % find max values and indexes according to cand0
            mas = max(cand0(:));
            [rr,cc] = find(cand0==mas);
            
            % save the results and exclude the corresponding values from overM
            mergMat(i,1) = mas;
            mergMat(i,2) = rr;
            mergMat(i,3) = cc;
            overM(rr,cc) = nan;
            
            % update indexes
            r = ind2sub(totRC, unique([r rr]));
            c = ind2sub(totRC, unique([c cc]));
            
            % update cand0
            cand0(r,:) = overM(r,:);
            cand0(:,r) = overM(:,r);
            cand0(:,c) = overM(:,c);
            cand0(c,:) = overM(c,:);
            
            % next step
            i=i+1;
            
        else % if there are no components with overlap>0 in selected rows and columns of
            % cand0 the search will restart for the remaining elements in overM matrix
            
            mas = max(overM(:));
            [rr,cc] = find(overM==mas);
            
            % add an empty row in matrix to distinguish these disjoint groups
            % from the previous, useful later to obtain the variable groups
            % unefficient method, to enhance
            mergMat(i,:) = 0;
            i = i+1;
            
            % save the results and delete the values from cand
            mergMat(i,1) = mas;
            mergMat(i,2) = rr;
            mergMat(i,3) = cc;
            overM(rr,cc) = nan;
            
            % update indexes
            r = ind2sub(totRC, unique([r rr]));
            c = ind2sub(totRC, unique([c cc]));
            
            % update cand0
            cand0(r,:) = overM(r,:);
            cand0(:,r) = overM(:,r);
            cand0(:,c) = overM(:,c);
            cand0(c,:) = overM(c,:);
            
            % next step
            i=i+1;
        end
    end
    
    % Index of the pairs of components above the threshold
    ind = mergMat(:,1) >= omegaStar;
    % REMARK - Now all the 1 in a row have to be merged together, a 0
    % after them indicate a disjoint cluster
    
    % Initialize variables to obtain a vector labelling each cluster
    % according to the rule used
    ng = 0; % initial value to count merged clusters
    ind_Prev = 0; % label obtained from the previous pair evaluated in the loop
    groups = zeros(length(ind), 1); % objective vector
    
    % Assign merged groups with progressive enumeration to gorups (zeros
    % point out separate clusters from the previous ones)
    for iii = 1:length(ind)
        ind_i = ind(iii);
        if ind_i == 0
            % Components not to merge
            groups(iii) = 0;
        elseif  ind_i == ind_Prev
            % Components to merge
            groups(iii) = groups(iii-1);
        else
            % Subsequent cluster
            ng = ng + 1;
            groups(iii) = ng;
        end
        % Update previous index
        ind_Prev = ind_i;
    end
    
    % Get unique values of the components to merge and store them in label cell array
    label = cell(max(groups),1);
    for iii = 1:max(groups)
        labelSol = unique(mergMat(groups==iii, 2:3));
        if size(labelSol, 1) == 1
            % If is a column vector (i.e. just 2 components to merge) transpose it
            labelSol = labelSol';
        end
        label{iii} = {labelSol};
    end
    
    % Find possible non-merged clusters
    % Specify all possible single clusters to find
    singleOnes = 1:k;
    % Loop through all the groups excluding the merged ones
    for iii = 1:max(groups)
        eachMergClu = cell2mat(label{iii});
        % Find non-merged components
        singleOnes = setdiff(singleOnes, eachMergClu, 'sorted');
    end
    
    % assign singleOnes at the end of labels array
    if ~isempty(singleOnes)
        label{max(groups)+length(singleOnes)} = [];  % initialize [[Needed?]]
        label(max(groups)+1:max(groups)+length(singleOnes)) = num2cell(singleOnes);
    end
    
    % get lebels for the merged groups referred to original units
    mergID = clu.idx; %initialization
    % looping through each component of merged and non merged components
    for j = 1:length(label)
        if iscell(label{j})
            % for merged components
            mergID(ismember(clu.idx, cell2mat(label{j}))) = j;
        else
            % for single clusters
            mergID(ismember(clu.idx, label{j})) = j;
        end
    end
    
    % Check for errors in g
elseif mod(g, 1)~=0 && dist~=1
    error('FSDA:txmerge:wrongInputs','Fractional g is only allowed for DEMP.')
elseif g >= 1 && g<=k && mod(g, 1)~=0 && dist~=1
    error('FSDA:txmerge:wrongInputs','For hierarchical clustering the argument ''g'' has to be an integer.')
elseif g>k
    error('FSDA:txmerge:wrongInputs','The argument g has to be smaller than the number of components k.')
end


%% Plotting phase

if ~(isscalar(plots) && plots==0)
    
    figure();
    % set(gcf,'numbertitle','off'. 'name','Merging of the components found by kmeans');
    
    if v == 1
        % Univariate case: plot the histogram
        histFS(Y(mergID>0), length(Y(mergID>0)), mergID(mergID>0));
        str = sprintf('Histograms of the merged components');
        title(str,'Interpreter','Latex');
        legend(num2str(unique(mergID(mergID>0))))
        % clickableMultiLegend(num2str(unique(mergID(mergID>0)))) % [To Adjust]
        
    elseif v >= 2
        
        % Bivariate plot, optionally with confidence ellipses, density
        % countours or bivariate boxplot
        
        % define what to superimpose on the plot
        if ischar(plots)
            overlay.type = plots;
        elseif isstruct(plots)
            overlay = plots;
        elseif plots==1
            % if plots=1 do not add anything to the scatter plot
            overlay ='';
        end
        
        % exclude outliers if present (when plots is char or struct)
        if any(mergID<=0) && ~isempty(overlay)
            overlay.include = true(length(unique(mergID)), 1);
            overlay.include(unique(mergID)<=0) = false;
        end
        
        % Add default axes label [Add the possibility to be specified by the user?]
        plo.labeladd=1;
        
        % start with a black color when there are some outliers
        if any(mergID<=0)
            plo.clr = 'kbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcykbrmgcy';
            plo.clr = plo.clr(1:length(unique(mergID)));
        end
        
        % differentiate for bivariate and multivariate data
        if v==2
            undock = [2 1];
        else
            undock = '';
        end
        
        % String for the title
        if g>1
            str = sprintf('Scatter plot of the merged components using hierarchical clustering');
        else
            str = sprintf('Scatter plot of the merged components using the threshold %s=%.4f', '$\omega^*$', g);
        end
        
        % bivariate scatter or scatter matrix
        spmplot(Y, 'group', mergID, 'plo', plo, 'undock', undock, 'overlay', overlay);
        title(str,'Interpreter','Latex');
    end
end


%% Getting all ouputs and store them in out structure

out=struct;

if dist == 1
    % triangular overlap matrix
    out.PairOver = overMsave;
end

% labels for units according to merging
out.mergID = mergID;

if txOut
    % when the output for tkmeans is requested
    out.txOut = clu;  % store it
end

if Ysave
    % store original data matrix
    out.Y = Ysave;
end


end

%%

function Z = linkage_txmerge(Y, method)
% linkage_txmerge create hierarchical cluster tree. It is a modified version
% of the MATLAB function linkage (see the MATLAB help page for additional
% informations).

% this is to speedup the string comparison later, in the switch statement
method = method(1:2);

n = size(Y,2);
m = ceil(sqrt(2*n)); % (1+sqrt(1+8*n))/2, but works for large n
if isa(Y,'single')
    Z = zeros(m-1,3,'single'); % allocate the output matrix.
else
    Z = zeros(m-1,3); % allocate the output matrix.
end

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  N denotes how many points are contained in each cluster.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;

% Square the distances so updates are easier.  The cluster heights will be
% square-rooted back to the original scale after everything is done.
if any(strcmp(method,{'ce' 'me' 'wa'}))
    Y = Y .* Y;
end

for s = 1:(n-1)
    if strcmp(method,'av')
        p = (m-1):-1:2;
        I = zeros(m*(m-1)/2,1);
        I(cumsum([1 p])) = 1;
        I = cumsum(I);
        J = ones(m*(m-1)/2,1);
        J(cumsum(p)+1) = 2-p;
        J(1)=2;
        J = cumsum(J);
        W = N(R(I)).*N(R(J));
        [v, k] = min(Y./W);
    else
        [v, k] = min(Y);
    end
    
    i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
    j = k - (i-1)*(m-i/2)+i;
    
    Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A
    
    % Update Y. In order to vectorize the computation, we need to compute
    % all the indices corresponding to cluster i and j in Y, denoted by I
    % and J.
    I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables
    U = [I1 I2 I3];
    I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
    J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];
    
    switch method
        case 'si' % single linkage
            Y(I) = min(Y(I),Y(J));
        case 'co' % complete linkage
            Y(I) = max(Y(I),Y(J));
        case 'av' % average linkage
            Y(I) = Y(I) + Y(J);
        case 'we' % weighted average linkage
            Y(I) = (Y(I) + Y(J))/2;
        case 'ce' % centroid linkage
            K = N(R(i))+N(R(j));
            Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v)./K)./K;
        case 'me' % median linkage
            Y(I) = (Y(I) + Y(J))/2 - v /4;
        case 'wa' % Ward's linkage
            Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
                N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
    end
    J = [J i*(m-(i+1)/2)-m+j]; %#ok<AGROW>
    Y(J) = []; % no need for the cluster information about j.
    
    % update m, N, R
    m = m-1;
    N(n+s) = N(R(i)) + N(R(j));
    R(i) = n+s;
    R(j:(n-1))=R((j+1):n);
end

if any(strcmp(method,{'ce' 'me' 'wa'}))
    Z(:,3) = sqrt(Z(:,3));
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);

end

%FScategory:CLUS-RobClaMULT
% © 2021 GitHub, Inc.
% Terms
% Privacy
% Security
% Status
% Docs
% Contact GitHub
% Pricing
% API
% Training
% Blog
% About
% Loading complete