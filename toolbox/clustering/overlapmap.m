function [out] = overlapmap(D, varargin)
%overlapmap produces an interactive overlap map
%
% <a href="matlab: docsearchFS('overlapmap')">Link to the help function</a>
%
% The function overlapmap plots the ordered pairwise overlap values between
% components. These components are ordered according to a specific rule:
% first the closest pair is plotted in the lowest left corner, then the
% components closer to the ones already included are plotted (when all of
% them have a zero overlap value with the ones already included, the
% closest pair between all the remaining ones is inserted). The overlap map
% can either shows with different colors the closeness between components
% (i.e. in a descriptive manner), or it becomes an interactive plot with a
% left click on the color bar, which find and visualize the closest
% components according to a specific threshold value $ \omega^* $ (i.e.
% omegaStar), which specifies the minimum paiwise overlap threshold value
% used to merge the components. The interactive process ends with a right
% click on the white grid in the upper left corner of the plot, it also
% updates the results creating in the workspace a new variable
% 'userOverlap'. See the More About section for further informations.
%
% Required input arguments:
%
%     D :       Informations to compute the overlap matrix. Structure.
%               D is a structure which can have the following fields (not
%               all of them are strictly required).
%
%   Admissable fields for the structure D:
%
%       D.idx = Label of the units. Vector. It is a vector with n
%               elements which assigns each unit to one of the k groups.
%               REMARK - labels<=0 denotes trimmed units.
%
%       D.Y  =  Input data. Matrix. Data matrix containining n
%               observations on v variables. Rows of Y represent
%               observations, and columns represent variables. Missing
%               values (NaN's) and infinite values (Inf's) are allowed,
%               since observations (rows) with missing or infinite values
%               will automatically be excluded from the computations.
%               When this field is specified the algorithm evaluate the
%               statistics of interest to obtain the overlap matrix, it
%               also allows the user to obtain additional plots when the
%               interaction is closed (using spmplot). When this field is
%               not specified the fields D.sigmaopt, D.muopt and D.siz are
%               required.
%
% D.sigmaopt  = v-by-v-by-k covariance matrices of the groups.
%
%
%     D.muopt = k-by-v matrix containing cluster centroid locations.
%
%       D.siz = Matrix or vector. If it is a matrix of size
%               k-by-3, where:
%               1st col = labels of the k components.
%               2nd col = number of observations in each component.
%               3rd col = percentage of observations in each
%               component.
%       REMARK: in case there is a field structure named emp containing the
%               same informations, these ones will be used
%                   Data Types - struct
%
% Optional input arguments:
%
% omegaStar:    Pairwise overlap threshold. Scalar. It is the value between
%               pairs of components considered disjunct if their overlap is
%               below omegaStar. If specified, these components would be
%               highlighted in the overlap map with an 'X' mark.
%               The default value is 0 (i.e. all components should be merged).
%                   Example - 'omegaStar', 0.01
%                   Data Types - single | double
%
% plots    :    Additional plot on the screen. Scalar, char or struct.
%               This arguments requires the presence of the field D.Y.
%               - If plots=0 (default) no additional plot is produced.
%               - If plots=1, at the end of the interaction with the overlap
%               map (i.e. right click on the white grid), the components
%               merged are shown using the spmplot function. In particular:
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
%                If plots is a struct it may contain the following fields:
%                 plots.type = a char specifying the type of superimposition
%                    Choices are 'contourf', 'contour', 'ellipse' or 'boxplotb'.
%               REMARK - The labels<=0 are automatically excluded from the
%                          overlaying phase, considering them as outliers.
%                   Example - 'plots', 1
%                   Data Types - single | double | string
%
% userColors:     Color used for the color map. Matrix or string. Check the
%                 colormap function for more informations.
%                   Example - 'userColors', winter
%                   Data Types - single | double | string
%
% Output:
%
%    out:   A structure containing the following fields
%
%           out.Ghat        = Estimated number of clusters in the data.
%
%           out.PairOver 	= Pairwise overlap triangular matrix (sum of
%                             misclassification probabilities) among
%                             components found by tkmeans.
%
%           out.mergID      = Label for each unit. It is a vector with n
%                             elements which assigns each unit to one of
%                             the groups obtained.
%                             REMARK - out.mergID<=0 denotes trimmed units.
%
%           out.merged      = Cell array containing the labels of the
%                             components merged together.
%
%           out.single      = Vector containing the labels of single clusters
%                             found, i.e. not merged with any other component.
%
%Optional Output:
%
%   userOverlap     : Updating of the results. Structure. userOverlap
%                     is obtained when the interaction with the overlap
%                     map is closed and is added in the Workspace.
%                     It contains the following fields, which represent
%                     an update of their corresponding variable in
%                     the structure out:
%                     - userOverlap.omegaStar = update of out.omegaStar.
%                     - userOverlap.Ghat = update of out.Ghat.
%                     - userOverlap.merged = update of out.merged.
%                     - userOverlap.single = update of out.single.
%
%
% More About:
%
% In the code 'overM' represents a triangular matrix, denoted as $ \Omega $,
% which contains the pairwise overlap values. The merging phase starts
% searching the maximum pairwise overlap value in $ \Omega $, i.e. $ \max
% (\Omega_{k k'}) $, and then deletes this value (e.g. setting it to NaN).
% This new matrix obtained is denoted as $ \Omega' $. The respective rows and
% columns corresponding to the element deleted in $ \Omega' $ are placed in
% a new matrix $ \Omega'' $. The algorithm progressively continue the same
% process, searching the highest pairwise overlap value in the components
% closest to the ones previously found, i.e. in the respective rows or
% columns of the components $ k $ and $ k' $. When the latter are all
% zeros, the process starts again considering the remaining values in
% $ \Omega $.
%
% The values $ \max(\Omega'_{k k'}) $ and the respective $ k $ and $ k' $
% labels are sequentially saved in a $ k(k-1)/2 \times 3 $ matrix MergMat.
%
% See also: dempk, tkmeans, tclust
%
%
% References:
%
%   Melnykov, V., Michael, S. (2020), Clustering Large Datasets by Merging
%   K-Means Solutions, Journal of Classification, Vol. 37, pp. 97–123,
%   https://doi.org/10.1007/s00357-019-09314-8
%
%   Melnykov, V. (2016), Merging Mixture Components for Clustering Through
%   Pairwise Overlap, "Journal of Computational and Graphical Statistics",
%   Vol. 25, pp. 66-90.
%
% Acknowledgements:
%
%   ...
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('overlapmap')">Link to the help page for this function</a>

%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example using tkmeans on geyser data.
    
    close all
    Y = load('geyser2.txt');
    k = 3;

    % using tkmeans
    out = tkmeans(Y, k*2, 0.05, 'plots', 1);
    overl_1 = overlapmap(out);

    % using tkmeans for a higher number of components
    out2 = tkmeans(Y, k*4, 0.05, 'plots', 1);
    overl_2 = overlapmap(out2);

    cascade;
%}


%{
    %% Example using M5data with tclust and tkmeans, specifying an
    % initial threshold omegaStar, a colormap, and allowing for additional
    % interactive plots.

    close all
    rng('default')
    rng(2)
    Y=load('M5data.txt');
    gscatter(Y(:,1),Y(:,2), Y(:,3))
    k = 3;

    out = tkmeans(Y(:,1:2), k*5, 0.2, 'plots', 'ellipse', 'Ysave', true);
    overl = overlapmap(out, 'omegaStar', 0.025, 'plots', 'contour', 'userColors', winter);
    
    rng('default')
    if verLessThan('matlab', '8.5')
        rng(5)
    else
        rng(1)
    end
    out_2 = tclust(Y(:,1:2), k*2, 0.2, 1, 'plots', 'contourf', 'Ysave', true);
    overl_2 = overlapmap(out_2, 'omegaStar', 0.0025, 'plots', 'contourf', 'userColors', summer);

    cascade;
%}

%{
    % Example using simdataset to create homogeneous and spherical clusters.
    % This output is used as input for the overlap map and then also
    % tkmeans and tclust solutions, for a higher number of components.
    
    close all

    % Specify k cluster in v dimensions with n obs
    k = 8;
    v = 4;
    n = 5000;
    % Generate 8 homogeneous spherical clusters
    rng('default')
    rng(10);
    out = MixSim(k, v, 'sph', true, 'hom', true, 'int', [0 10], 'Display', ...
        'off', 'MaxOmega', 0.005, 'Display','off');
    % 5 percent noise
    alpha0 = 0.05*n;
    % Simulating data
    [X, id] = simdataset(n, out.Pi, out.Mu, out.S, 'noiseunits', alpha0);
    % Plotting data
    figure;
    spmplot(X, 'group', id);
    str = sprintf('Simulated data with %d groups in %d dimensions and %d units', k, v, n);
    title(str,'Interpreter','Latex');

    % overlap map on simdataset output
    Inputs.Y = X;
    Inputs.idx = id;
    overlapmap(Inputs, 'plots', 'contourf');

    % overlap map on tkmeans solution for simdataset output
    out = tkmeans(X, k*4, 0.05, 'plots', 'contourf', 'Ysave', true);
    overlapmap(out, 'plots', 'contourf');

    out = tclust(X, 10, 0.05, 100, 'plots', 'contour', 'Ysave', true);
    overlapmap(out, 'plots', 'contourf');

    cascade;
%}

%{
    %% Example using simdataset to create heterogeneous and elliptical
    % clusters and using tkmeans output as input for the overlap map.
    
    close all
    % Specify k cluster in v dimensions with n obs
    k = 3;
    v = 2;
    n = 50000;
    % restriction factor
    restr = 30;
    % Maximum overlap
    maxOm = 0.005;
    % Generate heterogeneous and elliptical clusters
    rng('default')
    rng(500, 'twister');
    out = MixSim(k, v, 'sph', false, 'restrfactor', restr, 'int', [0 10], ...
        'Display', 'off', 'MaxOmega', maxOm, 'Display','off');
    % null noise
    alpha0 = 0;
    % Simulating data
    [X, id] = simdataset(n, out.Pi, out.Mu, out.S, 'noiseunits', alpha0);
    % Plotting data
    gg = gscatter(X(:,1), X(:,2), id);
    str = sprintf('Simulated data with %d groups in %d dimensions and %d units, \n with restriction factor %d and maximum overlap %.2f', ...
        k, v, n, restr, maxOm);
    title(str,'Interpreter','Latex');
    
    % use tkmeans for a larger number of cluster and without trimming
    tkm = tkmeans(X, k*3, 0,'plots', 2,'Ysave',true, 'plots', 'ellipse');
    
    % overlap map with interctive mode
    overl = overlapmap(tkm, 'omegaStar', 0.01, 'plots', 'contourf');

    cascade;
%}

%{
    %% Example using simdataset to create homogeneous and spherical clusters
    % and using tkmeans.

    clear variables; close all
    % Specify k cluster in v dimensions with n obs
    k = 10;
    v = 2;
    n = 5000;
    % Generate homogeneous and spherical clusters
    rng('default')
    rng(100, 'twister');
    out = MixSim(k, v, 'sph', true, 'hom', true, 'int', [0 10], 'Display', 'off', 'BarOmega', 0.05, 'Display','off');
    % Simulating data
    [X, id] = simdataset(n, out.Pi, out.Mu, out.S);
    % Plotting data
    gscatter(X(:,1), X(:,2), id);
    str = sprintf('Simulated data with %d groups in %d dimensions and %d units', k, v, n);
    title(str,'Interpreter','Latex');
    clickableMultiLegend(num2str((1:k)'));

    % use tkmeans for a larger number of cluster and without trimming
    tkm = tkmeans(X, k*3, 0,'plots', 2,'Ysave',true, 'plots', 'ellipse');
 
    % overlap map with interctive mode
    out = overlapmap(tkm, 'omegaStar', 0.01, 'plots', 'contourf');

    cascade;
%}
%% Beginning of code

% ... Checking required inputs ...

%% Optional parameters

% creating option structure
options = struct('omegaStar', 0, 'plots', 0, 'userColors', autumn);

[varargin{:}] = convertStringsToChars(varargin{:});
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

if nargin > 1
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    % Check if the specified optional arguments are reasonable
    if ~isscalar(options.omegaStar) || options.omegaStar<0 || options.omegaStar>1
        error('FSDA:WrongOmegaStar','The argument OmegaStar must be a scalar in the interval [0 1].')
    end
    if ~((isscalar(options.plots) && (options.plots == 0 || options.plots == 1)) ||  ischar(options.plots) ||  isstruct(options.plots))
        error('FSDA:Wrongplots','The argument plots is wrongly specified.')
    end
    if size(options.userColors, 2)~=3
        error('FSDA:WronguserColors','The argument userColors must be a 3 columns matrix.')
    end
    
end

% Assign the values for the optional arguments
omegaStar = options.omegaStar;
plots = options.plots;
userColors = options.userColors;

%% Overlap Matrix

% Unique clusters
if isfield(D,'emp') && isstruct(D.emp)
    % use the empirical values computed by tkmeans or tclust
    id = D.emp.idxemp;
else
    % use the field idx assigned in D
    id = D.idx;
end
UniqID = unique(id);
% Length of the dataset
n = length(id);
% Actual components (excluding trimmed units)
k = sum(UniqID>0);

if isfield(D, 'muopt') && isfield(D, 'sigmaopt') && isfield(D, 'siz')
    % If D is already a structure containing the needed informations
    
    if ~isstruct(D.emp)
        % use the optimal arguments
        Mu = D.muopt;                   % Centroids
        S = D.sigmaopt;                 % Covariance matrices
        Pi = D.siz(D.siz(:, 1)>0, 3);   % Mixing proportions
    else
        % use the empirical quantities
        Mu = D.emp.muemp;                   % Centroids
        S = D.emp.sigmaemp;                 % Covariance matrices
        Pi = D.emp.sizemp(D.emp.sizemp(:, 1)>0, 3);   % Mixing proportions
    end
    
    % Dimension of the dataset
    v = size(Mu, 2);
    
    % Eliminate any empty component
    if sum(isnan(Mu(:,1)))>0
        indd = ~isnan(Mu(:,1));
        Mu = Mu(indd,:);
        S = S(:,:,indd);
    end
    
else
    % Y are the data used to compute the required statistics
    Y = D.Y;
    
    % Dimension of the dataset
    v = size(Y, 2); % should be checked before that the columns are the variables?
    
    % Excluding labels <=0, i.e. outliers
    if any(UniqID<=0)
        idIn = id(id>0);
        Yin = Y(id>0,:);
    else
        idIn=id;
        Yin = Y;
    end
    
    % Evaluate Mu and S for each component
    Mu = zeros(k,v);
    S = zeros(v,v,k);
    for i=1:k
        Mu(i,:) = mean(Yin(idIn==i,:));
        S(:,:,i) = cov(Yin(idIn==i,:));
    end
    
    % Evaluate Pi
    Pi = tabulate(idIn);
    Pi = Pi(:, 3);
    
end

% check consistency
if size(Mu, 1)~=k || size(S, 3)~=k
    error('FSDA:overlapmap:wrongInputs','The number of rows in k, Mu and S do not correspond');
end

% Give an error if any covariance matrix of the components is close to zero
% (because the eigenvalues wouldn't be computable)
if any(sum(sum(S))<=eps)
    errS = UniqID(sum(sum(S))<=eps);
    error('FSDA:overlapmap:notfullrank','The component number %d has a null coariance matrix. Change the components.', errS);
end

% Evaluation of misclassification probabilities
try
    [OmegaMap, ~, ~, ~, ~] = overlap(k, v, Pi, Mu, S);
catch
    error('FSDA:overlapmap:wrongInputs','Numerical issues in computing the overlap or empy groups: try with different parameters.');
end

% Sum of misclassification probailities in order to obtain an overlap matrix
overM = triu(OmegaMap,1)+(tril(OmegaMap,-1))';
overMsave = overM; % to store it


%% Order of the components

% Set to nan elements below the main diagonal of overM
d1 = size(overM, 1);
In = logical(tril(ones(d1)));
overM(In) = nan;

% Initialize mergMat (merging matrix).
% If overM has some zeros these rows are not considered in mergMat (i.e. not overlapping at all)
mergMat = zeros(k*(k-1)/2 - sum(overM(:)==0), 3);
% Otherwise, and should be commented the % 2nd statement of the while loop later
% mergMat = zeros(k*(k-1)/2, 3);

% Initialize cand0: it will be filled with the overlap values of the
% components closer to some component already included in mergMat(:,2:3)
% (i.e. the corresponding rows/columns in overM)
cand0 = nan(size(overM));

% Starting index value which goes row-wise in mergMat
i = 1;

% Initialize indexes that store the location of the maximum overlap value
r = [];
c = [];

% totRC (total rows and columns) is used to store the indexes of the
% components already included in mergMat and hence added or to add in cand0
totRC = 1:length(overM);

% Components sorted according to their pairwise overlap values
while ~all(isnan(overM(:)) | overM(:) == 0)  % 2nd statement previously described to comment if needed
    
    % Stopping rule: the overlap matrix overM has been totally evaluated
    %(i.e. only NaN) or the remained overlaps are 0
    
    if any(cand0(:))~=0
        
        % When cand0 does not contain only zeros (the opposite happens when
        % there are no close components to the ones already found or at the
        % 1st iteration of the algorithm), find the maximum overlap value
        % and its indexes in cand0
        mas = max(cand0(:));
        [rr,cc] = find(cand0==mas, 1, 'first'); % choose the first when some entries are the same
        
    else
        
        % If there are no components with overlap>0 in cand0 on a specific
        % step (e.g. the 1st) the maximum and its indexes are researched
        % among the elements of overM.
        % (REMIND 1 - This way to approach the probelem is not specified by
        % the authors of the paper for the step after the 1st, but is the
        % easiest strategy found)
        % (REMIND 2 - Possible problem: if some entries have exactly the
        % same max values how to choose? We are arbitrarily choosing the
        % first value found)
        
        mas = max(overM(:));
        [rr,cc] = find(overM==mas, 1, 'first'); % choose the first when some entries are the same
        
        if i>1
            
            % When it happens ater the 1st iteration, we trivially
            % add an empty row in mergMat to keep track that the skipping
            % process happened (i.e. to distinguish these disjoint groups
            % from the previous found). It will be useful later to run the
            % algorithm to merge the groups.
            % (REMIND - not such an efficient method probably, to enhance)
            mergMat(i,:) = 0;
            i = i+1;                            % skip line in mergMat
            % mergMat = [mergMat; zeros(1, 3)];    % augment mergMat rows [USELESS?]
            
        end
    end
    
    % Save the max pairwise overlap value found and its indexes in
    % the respective columns of mergMat
    mergMat(i,1) = mas;
    mergMat(i,2) = rr;
    mergMat(i,3) = cc;
    % Exclude the max value in overM (i.e. set to NaN)
    overM(rr,cc) = nan;
    
    % Update indexes
    r = ind2sub(totRC, unique([r rr]));
    c = ind2sub(totRC, unique([c cc]));
    
    % Update cand0 (i.e. add the new rows/columns found)
    cand0(r,:) = overM(r,:);
    cand0(:,r) = overM(:,r);
    cand0(:,c) = overM(:,c);
    cand0(c,:) = overM(c,:);
    
    % Next iteration
    i=i+1;
    
end


%% Preliminaries to obtain the overlap map

% ord is a vector containing the correct descending order of the pairs
% of components found according to the specific rule used
resh = reshape(mergMat(:,2:3)', size(mergMat, 1)*2, 1);
ord = unique(resh, 'stable');
% exclude possible zeros assigned trivially in mergMat(:, 2:3)
ord(ord==0) = [];

% if some components are missing in ord (e.g thay had 0 overlap with all
% the others components found) add them at the end of ord (apart from outliers)
if sum(~ismember(UniqID(UniqID>0), ord))
    misComp = find(~ismember(UniqID(UniqID>0), ord));
    ord(end+1:end+length(misComp)) = misComp;
end

% Re-load the overlap triangular matrix
overM = overMsave;
% Obtain a symmetrical overlap matrix
candInt = triu(overM)'+ overM;
% Reorder the overlap matrix according to the specific rule used
overMorder = zeros(size(overM)); % initialize
for i = 1:length(overM)
    for j = 1:length(overM)
        overMorder(i, j) = candInt(ord(i), ord(j));
    end
end
% Set to NaN the elements below the main diagonal
overMorder(In) = nan;

% Estimate the number of clusters
if omegaStar>0
    % If omegaStar is specified call the sub-function mergComp
    [Ghat, label, singleOnes, mergID] =  mergComp(omegaStar, mergMat, k, id);
else
    % Without use of omegaStar (i.e. all components should ideally be merged)
    Ghat = 1;
end


%% Save the outputs in out structure

% Estimate of the number of cluster
out.Ghat = Ghat;

% Triangular overlap matrix
out.PairOver = overMsave;

% store merged groups according to a cut-off omegaStar
if omegaStar > 0
    out.mergID = mergID;        % vector of length n containing the new ID found
    out.merged = label;         % cell containing merged components
    out.single = singleOnes;    % vector containing id of the non merged components
end


%% Overlap map

if isreal(overMorder)
    figure();
    set(gcf,'name','Overlap map');
    
    % Subplot 1: main plot of overlap values
    s1 = subplot(2,1,1);
    % Produce the plot
    surface(zeros(size(overMorder)),overMorder(:, 2:end),'EdgeColor','none');
    axis tight;
    colormap(flipud(userColors));
    ax = gca;
    if verLessThan('matlab', '8.5')
        set(ax,'XAxisLocation','top',...
            'YTick',1:length(ord(2:end)),'YTickLabel',[],...
            'XTick',1:length(ord(2:end)),'XTickLabel',[] );
    else
        ax.XAxisLocation    = 'top';
        ax.YTick            = 1:length(ord(2:end));
        ax.YTickLabel       = [];
        ax.XTick            = 1:length(ord(2:end));
        ax.XTickLabel       = [];
    end
    ylabel('Ordered components','FontSize',14,'Interpreter','Latex');
    % Add the text about the corresponding ordered clusters (resized according
    % to the dimension of the map) [To Enhance]
    axes(ax);
    if k<=4
        text(ones(1, k)+0:k, 1:k, num2str(ord(1:end)), ...
            'FontWeight','bold','FontSize',18,'HorizontalAlignment','right','VerticalAlignment','baseline');
    elseif k<=7
        text(ones(1, k)/1.1+0:k, 0.1+(1:k), num2str(ord(1:end)), ...
            'FontWeight','bold','FontSize',16,'HorizontalAlignment','right','VerticalAlignment','baseline');
    elseif k<=15
        text(ones(1, k)/1.2+0:k, 0.3+(1:k), num2str(ord(1:end)), ...
            'FontWeight','bold','FontSize',14,'HorizontalAlignment','right','VerticalAlignment','baseline');
    elseif k<=35
        text(ones(1, k)/2+0:k, 0.3+(1:k), num2str(ord(1:end)), ...
            'FontWeight','bold','FontSize',12,'HorizontalAlignment','right','VerticalAlignment','baseline');
    else
        text(ones(1, k)/2.5+0:k, 0.3+(1:k), num2str(ord(1:end)), ...
            'FontWeight','bold','FontSize',10,'HorizontalAlignment','right','VerticalAlignment','baseline');
    end
    % Add box and grid
    box 'on';
    grid on;
    
    % Subplot 2: extract maximum value in each column of subplot 1
    s2 = subplot(2,1,2);
    surface(zeros(2, size(overMorder,1)),max(overMorder(:, 2:end)),'EdgeColor','black');
    axis tight;
    ax2 = gca;
    set(ax2,'YTick', [], 'XTick', []);
    % ~verLessThan('matlab', '8.5')
    %     ax2.YTick = [];
    %     ax2.XTick = [];
    
    % Moving and resizing subplots
    if verLessThan('matlab', '8.5')
        s1p = get(s1,'Position');
        s1p(1) = s1p(1)-0.08;
        s1p(4) = s1p(4)*2.3;
        s1p(2) = s1p(2)*0.25;
        set(s1,'Position',s1p);
        
        s2p = get(s2,'Position');
        s2p(4) = s2p(4)*0.2;
        s2p(2) = s2p(2)*0.4;
        s2p(3) = s1p(3);
        s2p(1) = s2p(1)-0.035;
        set(s2,'Position',s2p);
    else
        s1.Position(1) = s1.Position(1)-0.08; % moving left
        s1.Position(4) = s1.Position(4)*2.3; % increasing size 230%
        s1.Position(2) = s1.Position(2)*0.25; % moving down
        s2.Position(4) = s2.Position(4)*0.2; % reducing size to 20%
        s2.Position(2) = s2.Position(2)*0.4; % moving down
        s2.Position(3) = s1.Position(3); % equalize horizontally
        s2.Position(1) = s2.Position(1)-0.08; % moving left
    end
    % Colormap position in the plot (in normalized units)
    colorPos = [0.84875 0.0455512465373961 0.045875 0.885503231763619];
    % Add a common colorbar for both subplots
    if max(overMorder(:))>0.1
        % Set a specific scale when the max overlap is not so small [Useful?]
        if verLessThan('matlab', '8.5')
            rmo = (max(overMorder(:))-0.1)/5;
            rmo = round(rmo*100)/100;
            co = colorbar('Position', colorPos);
            set(co,'YTick',[0.01, 0.04:0.02:0.1, 0.12:rmo:max(overMorder(:))],...
                'YTickLabel',{num2str([0.01, 0.04:0.02:0.1, 0.12:rmo:max(overMorder(:))]')}, ...
                'Position', colorPos, 'YLim', [0, max(overMorder(:))+0.001]);
        else
            rmo = round((max(overMorder(:))-0.1)/5,2);
            co = colorbar('Ticks',[0.01, 0.04:0.02:0.1, 0.12:rmo:max(overMorder(:))],...
                'TickLabels',{num2str([0.01, 0.04:0.02:0.1, 0.12:rmo:max(overMorder(:))]')}, ...
                'Position', colorPos, 'Limits', [0, max(overMorder(:))+0.001]);
        end
        
    elseif verLessThan('matlab', '8.5')
        % Use default scale (to avoid errors)
        co = colorbar('Position', colorPos);
        
    else
        % Use default scale (to avoid errors)
        co = colorbar('Position', colorPos, 'Limits', [0, max(overMorder(:))+0.001]);
    end
    caxis([0 max(max(overMorder))]); % set equal for both subplots
    % Add text
    if verLessThan('matlab', '8.4')
        axes(co)
        ylabel('Pairwise overlap values','FontSize',14,'Interpreter','Latex');
    else
        co.Label.String = 'Pairwise overlap values';
        co.Label.FontSize = 14;
        co.Label.Interpreter = 'Latex';
    end
    
    % Add figure title
    str = sprintf('%d groups in %d dimensions and %d units', k, v, n); % Overlap map with ..
    title(ax, str,'Interpreter','Latex','fontsize',14);
    
    % Set callback to start the interaction with the plot, activated by a
    % left click on the colorbar
    set(co,'ButtonDownFcn', {@colormapClickSx, k, id, omegaStar, Ghat, overMorder, mergMat, ax, ax2, co, D, plots});
    
else
    msgbox('Numerical issues in computing the overlap (complex values are generated): try with different parameters.');
end

end

%% Left click on the Colormap
% This function runs as there is a left click on the colormap and activate
% the interactive use of the overlap map.

function colormapClickSx(~, ~, k, id, omegaStar, Ghat, overMorder, mergMat, ax, ax2, co, D, plots)

if strcmp(get(gcf,'SelectionType'),'normal') % left click
    
    % Show a message on the console
    disp(' ');
    disp('-----------------------------------------');
    disp('Left Click: Interactive Overlap Map ''on''');
    disp('-----------------------------------------');
    disp(' ');
    
    % To avoid problems caused by multiple clicks delete these objects
    delete(findobj('Tag','GlobalTextHandle'));
    delete(findobj('Tag','GlobalLineHandle'));
    delete(findobj('Tag','GlobalXHandle'));
    
    % Plot omegaStar and Ghat values in the upper left corner
    % Position on the Y axes adapted to the number of groups [To Enhance]
    if k<=6
        TextPos = k-0.5;
    elseif k<10
        TextPos = k-1.5;
    elseif k<40
        TextPos = k-3;
    else
        TextPos = k-10;
    end
    str = sprintf('%s = %.3f\n%s = %d', '$\omega^{*}$', omegaStar, '$\hat{G}$', Ghat);
    axes(ax);
    text(1.1, TextPos, str, 'FontWeight','bold','FontSize',14,'Interpreter', 'Latex', 'Tag', 'GlobalTextHandle');
    
    % Add line to the colorbar according to the threshold omegaStar chosen
    if verLessThan('matlab', '8.5')
        % Line to be added
        %         line_axes = axes('position', get(co, 'Position'), 'color', 'none', 'visible','off');
        %         line(get(line_axes, 'XLim'), omegaStar*[1 1], 'color', 'white', 'parent', line_axes, ...
        %             'color','k','LineWidth',3, 'Tag','GlobalLineHandle');
    else
        line_axes = axes('position', co.Position, 'ylim', co.Limits, 'color', 'none', 'visible','off');
        line(line_axes.XLim, omegaStar*[1 1], 'color', 'white', 'parent', line_axes, ...
            'color','k','LineWidth',3, 'Tag','GlobalLineHandle');
    end
    
    % Highlights the merged components below the threshold with a 'X' mark
    below = max(overMorder(:, 2:end)) < omegaStar;  % index of the max for each columns
    leng = 0.5+1:k;                                 % position on the x axes
    highl = leng(below);                            % index of the marks to add
    axes(ax2);  % deletes the line if ~verLessThan('matlab', '8.5')
    text(highl-0.1, 0.5+ones(1, length(highl)) , 'x', ...
        'FontWeight','bold','FontSize',20, 'Tag','GlobalXHandle',...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    
    % Add a slider and its callback
    uicontrol(gcf, 'Style','slider', 'Units','normalized','Position',...
        get(co, 'Position')-[-0.018 0.0175 0.035 -0.032], 'value', omegaStar, ...
        'min',0, 'max',max(overMorder(:)+0.001),'SliderStep',[0.001 0.010], 'Tag','slider1', ...
        'Callback', {@MoveSlider,  mergMat, k, id, overMorder, ax, ax2, co, TextPos, D, plots});
    
    % Set callback to close the interaction with the plot, disabled by a
    % right click on the white grid
    % Adjust interaction when verLessThan('matlab', '8.5')
    %     if verLessThan('matlab', '8.5')
    %         ax = axes('position', get(ax, 'Position'), 'color', 'none', 'visible','off');
    %         set(ax,'ButtonDownFcn', {@colormapClickDx, omegaStar, mergMat, k, id, D, plots});
    %     else
    set(ax,'ButtonDownFcn', {@colormapClickDx, omegaStar, mergMat, k, id, D, plots});
    %     end
end
end


%% Slider function
% This function runs as the slider is moved and updates the statistics of
% interest.

function [omegaStar, Ghat] = MoveSlider(slider, ~, mergMat, k, id, overMorder, ax, ax2, co, TextPos, D, plots)

% Update omegaStar
omegaStar = get(slider, 'Value');

% Update results calling mergComp
[Ghat, ~, ~, ~] = mergComp(omegaStar, mergMat, k, id);

% Update the visualization of omegaStar and Ghat values on the plot
delete(findobj('Tag','GlobalTextHandle'));
str = sprintf('%s = %.3f\n%s = %d', '$\omega^{*}$', omegaStar, '$\hat{G}$', Ghat);
axes(ax);
text(1.1, TextPos, str, 'FontWeight','bold','FontSize',14,'Interpreter','Latex', 'Tag', 'GlobalTextHandle');

% Update the line to colorbar according to the threshold omegaStar
delete(findobj('Tag','GlobalLineHandle'));
if verLessThan('matlab', '8.5')
    % Line to be added
    %     h_axes = axes('position', get(co, 'Position'), 'ylim', caxis(co), 'color', 'none', 'visible','off');
    %     line(get(h_axes, 'XLim'), omegaStar*[1 1], 'color', 'white', 'parent', h_axes,'color','k','LineWidth', 3, 'Tag', 'GlobalLineHandle');
else
    h_axes = axes('position', co.Position, 'ylim', co.Limits, 'color', 'none', 'visible','off');
    line(h_axes.XLim, omegaStar*[1 1], 'color', 'white', 'parent', h_axes,'color','k','LineWidth', 3, 'Tag', 'GlobalLineHandle');
end
% Highlights the merged components below the threshold with an 'X' mark
delete(findobj('Tag', 'GlobalXHandle'));
below = max(overMorder(:, 2:end)) < omegaStar;
leng = 0.5 + 1:length(overMorder(2:end,:)-1);
highl = leng(below);

axes(ax2); % deletes the line if ~verLessThan('matlab', '8.5')
text(highl-0.1, 0.5+ones(1, length(highl)) , 'x', 'FontWeight','bold', 'FontSize', 20,...
    'Tag', 'GlobalXHandle', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% Preserve the callbacks
set(ax,'ButtonDownFcn', {@colormapClickDx, omegaStar, mergMat, k, id, D, plots});
set(co,'ButtonDownFcn', {@colormapClickSx, k, id, omegaStar, Ghat, overMorder, mergMat, ax, ax2, co, D, plots})

end


%% Right click on the main plot [[[[ To enhance ]]]
% This function runs as there is a right click on the grid of the plot and
% close the interactive process showing and updating the results obtained.

function [userOverlap] = colormapClickDx(~,~,omegaStar, mergMat, k, id, D, plots)

if strcmp(get(gcf,'SelectionType'), 'alt')   % Right click
    
    % Evaluate the solution [To Enhance, it may be included in MoveSlider function]
    [Ghat, label, singleOnes, mergID] =  mergComp(omegaStar, mergMat, k, id);
    
    % Show a message on the console regarding the end of interaction
    disp(' ');
    disp('--------------------------------------------');
    disp('Right Click: Interactive Overlap Map ''off''');
    disp('--------------------------------------------');
    disp(' ');
    
    % Show messages on the console regarding interactive results obtained
    % Merged clusters
    disp(' ');
    disp('--------------------------------------------');
    if ~isempty(label)
        disp(' ');
        for jj =1:length(label)
            disp(['Components merged in cluster ', num2str(jj), ':']);
            disp(cell2mat(label{jj})');
        end
    else
        disp('No components merged.');
    end
    
    % Non-merged clusters
    disp(' ');
    if ~isempty(singleOnes)
        disp('Single clusters found:');
        disp(singleOnes);
    else
        disp('No single clusters found.');
    end
    disp('--------------------------------------------');
    disp(' ');
    
    % Delete Slider
    delete(findobj('Tag', 'slider1'));
    
    % Create a new variable 'userOverlap' in the workspace which updates
    % and saves the results of interest
    userOverlap = struct();
    userOverlap.Ghat = Ghat;
    userOverlap.omegaStar = omegaStar;
    userOverlap.merged = label;
    userOverlap.single = singleOnes;
    userOverlap.mergID = mergID;
    assignin('base', 'userOverlap', userOverlap);
    % The same output is saved on the UserData property of the figure
    % containing the overlap map
    set(gcf,'UserData',userOverlap);
    
    % Produce an additional plot of the components merged using spmplot
    if ~(isscalar(plots) && plots==0) && isfield(D, 'Y')
        % When plots is specified and the data Y are given
        
        figure;
        
        % To avoid errors calling spmplot for the first time (i.e. on call back figure)
        hp1 = pan(gcbf);
        set(hp1, 'ActionPreCallback', 'none');
        
        % When Y data are given produce a scatter plot
        Y = D.Y;
        [~, v] = size(Y);
        
        if v == 1
            % Univariate case: plot the histogram
            histFS(Y(mergID>0), length(Y(mergID>0)), mergID(mergID>0));
            str = sprintf('Histograms of the merged components');
            title(str,'Interpreter','Latex', 'fontsize', 12);
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
            if v>1
                str = sprintf('Scatter plot of the merged components');
            else
                str = sprintf('Scatter plot matrix of the merged components');
            end
            
            % bivariate scatter or scatter matrix
            spmplot(Y, 'group', mergID, 'plo', plo, 'undock', undock, 'overlay', overlay);
            title(str,'Interpreter','Latex'); % , 'fontsize', 12
        end
    end
end

end


%% This is the function which actually compute the merging of the components.

function [Ghat,label, singleOnes, mergID] = mergComp(omegaStar, mergMat, k, id)

% number of units
n = length(id);

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

% Estimated number of groups in the data set
Ghat = max(groups) + length(singleOnes);

% Create the new id vector denoted mergID
mergID = zeros(n, 1);
% Assign new labels to the merged clusters
if max(groups) > 0
    for iii = 1:max(groups)
        posit = ismember(id, cell2mat(label{iii}));
        mergID(posit) = iii;
    end
else
    % To be able to evaluate the single clusters correctly later
    iii = 0;
end

% Assign new labels to the non-merged clusters
if ~isempty(singleOnes)
    for jjj = iii+1:iii+length(singleOnes)
        posit = id==singleOnes(jjj-iii);
        mergID(posit) = jjj;
    end
end

end

%FScategory:CLUS-RobClaMULT
