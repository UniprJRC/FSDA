function [hax, out] = sphereprojplot(Y, varargin)
%sphereprojplot plots high-dimensional directional data projected on a 3-dimensional space
%
%<a href="matlab: docsearchFS('sphereprojplot')">Link to the help page for this function</a>
%
%   sphereprojplot is the counterpart of sphereplot for unit vectors in
%   dimension p>3 (for example document embeddings, where p is in the
%   hundreds). The data are projected orthogonally on a 3-dimensional
%   linear subspace of R^p, chosen automatically or by the user, and the
%   projections are drawn with sphereplot inside the unit sphere of that
%   subspace.
%
%   This is a separate function because, differently from sphereplot,
%   the plot is NOT an exact representation of the data, and it must be
%   read in a different way:
%   1) an orthogonal projection shortens the vectors, so that the points
%   lie inside the sphere and not on its surface. The squared distance of
%   a point from the origin is the fraction of that unit which is "seen"
%   in the chosen subspace; the average of this quantity is reported in
%   the title of the plot and in the output (out.retained);
%   2) the within-group variability is spread over p-1 directions, of
%   which only 2 are visible. The groups therefore look much tighter
%   than they are;
%   3) the confidence regions are drawn as the projection of the
%   p-dimensional caps (see "More About"): a circle on the sphere which
%   delimits a flat disc. Units inside the cap are projected on the side
%   of the disc which is opposite to the origin.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-p matrix of unit vectors (one
%               per row), with p>3.
%
%  Optional input arguments:
%
%        proj : Projection. Character vector, vector of 3 indices or
%               p-by-3 matrix. It defines the 3-dimensional subspace.
%               'means' = the subspace which contains the mean directions
%               of the groups (the centres of the regions in cap if they
%               are given, otherwise the directions of the resultant
%               vectors of the groups in group). This is the subspace
%               where the groups are best separated. With 3 groups it is
%               determined by the 3 mean directions; with fewer groups it
%               is completed with the principal directions of the part of
%               the data which is orthogonal to the means; with more
%               than 3 groups the 3 principal directions of the set of
%               mean directions are used.
%               'pca' = the first 3 principal directions of the data,
%               i.e. the eigenvectors of Y'*Y associated with the 3
%               largest eigenvalues. Y is not centred: the subspace must
%               contain the origin, which is the centre of the sphere.
%               This is the subspace which retains the largest part of
%               the squared length of the units.
%               'auto' (default) = 'means' if at least two groups or one
%               region are given, 'pca' otherwise.
%               A vector of 3 distinct integers between 1 and p selects 3
%               of the original coordinates.
%               A p-by-3 matrix defines the subspace spanned by its
%               columns; the columns are orthonormalised if necessary.
%                 Example - 'proj','pca'
%                 Data Types - char | double
%       group : Grouping variable. Vector of length n of numeric labels,
%               as in sphereplot (0 = trimmed units).
%                 Example - 'group',out.idx
%                 Data Types - double
%        cap : Confidence regions to superimpose. Structure or cell array.
%               Default is empty (no region). Two alternatives.
%               1) von Mises-Fisher caps: structure with fields:
%                  cap.mu = k-by-p matrix of mean directions.
%                  cap.kappa = vector of length k of concentrations.
%                  The caps are computed at level conflev.
%               2) Cell array of length k, as out.ellipse of tclustvmf:
%                  cap{j} = matrix with p columns of points on the
%                  boundary of the p-dimensional cap of group j.
%                  In this case supply also option center.
%                 Example - 'cap',struct('mu',out.muopt,'kappa',out.kappaopt)
%                 Data Types - struct | cell
%      center : Centres of the regions. k-by-p matrix of mean directions,
%               for example out.muopt of tclustvmf. Used only when cap is
%               a cell array. If it is empty the centre of a region is
%               the direction of the mean of its boundary points.
%                 Example - 'center',out.muopt
%                 Data Types - double
%     conflev : Confidence level(s). Scalar or vector with elements in
%               (0,1). Default 0.95. Used only when cap is a structure.
%                 Example - 'conflev',[0.5 0.95]
%                 Data Types - double
%        npts : Number of points of each contour. Scalar. Default 200.
%                 Example - 'npts',400
%                 Data Types - double
%  noisegroup : Label of the noise units. Scalar. Label of the units of the
%               noise component, as in sphereplot.
%                 Example - 'noisegroup',3
%                 Data Types - double
%   plottitle : Add a title with the type of projection and the
%               percentage of squared length retained. Logical. Default
%               true.
%                 Example - 'plottitle',false
%                 Data Types - logical
%
%   All the other options of sphereplot (capfill, capalpha, capflat,
%   plotcenter, showpoints, labels, clr, markersize, spherealpha, grid,
%   rotate, legend, h) can be used and are passed to sphereplot.
%
%  Output:
%
%         hax : Axes handle. Handle. Handle of the axes which contain the plot.
%         out : Structure with the following fields.
%               out.V = p-by-3 matrix with orthonormal columns which
%                    defines the projection: the points which are
%                    plotted are Y*out.V. It can be passed as option
%                    proj in other calls, to show different sets of data
%                    in the same subspace.
%               out.Yproj = n-by-3 matrix of projected data Y*out.V.
%               out.retained = scalar in [0,1]. Average over the units of
%                    the squared norm of the projection, i.e. fraction of
%                    the squared length of the units retained by the
%                    projection.
%               out.retainedunit = n-by-1 vector with the squared norm of
%                    the projection of each unit.
%               out.type = character vector with the type of projection
%                    ('means', 'pca', 'coordinates' or 'user').
%
% More About:
%
% PROJECTION OF A CAP
%
% The confidence region of a vMF group in R^p is the cap
% {y: ||y||=1, mu'*y>=t}. Let V be the p-by-3 matrix of the projection,
% m=V'*mu the projection of the mean direction, and x=V'*y the
% projection of y. The points y of the cap which have projection x are
% those for which mu'*y = m'*x + (component outside the subspace) >= t.
% When mu lies in the subspace (||m||=1, which is the case for
% proj='means' with up to 3 groups) the projection of the cap is exactly
% the set {x: ||x||<=1, m'*x>=t}: the part of the unit ball beyond the
% plane m'*x=t. Its boundary is made of a cap of the unit sphere and of a
% flat disc, which have in common the circle {||x||=1, m'*x=t}.
% sphereprojplot draws this circle and (with capfill) shades the disc.
% When mu does not lie in the subspace, the circle which is drawn is
% {||x||=1, m'*x=t}, that is the intersection between the projected
% region and the surface of the sphere. If t>=||m|| this intersection is
% empty and the region is not drawn.
%
% THRESHOLD OF THE CAP
%
% t is such that P(mu'*y>=t)=conflev when y is vMF(mu,kappa) in dimension
% p, and it is computed by vMFcap (second output argument).
%
% References:
%
%   Mardia, K.V. and Jupp, P.E. (2000), Directional Statistics, Wiley,
%   Chichester.
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
% See also: sphereplot, tclustvmf, vMFrnd, vMFcap
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('sphereprojplot')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Three vMF groups in dimension p=100.
    % The data are shown in the subspace of the three mean directions,
    % with the projection of the 95 per cent caps.
    rng(1)
    p = 100; n = 200;
    mu = zeros(3,p); mu(1,1) = 1; mu(2,2) = 1; mu(3,3) = 1;
    kappa = [400; 200; 300];
    Y = [vMFrnd(n, mu(1,:), kappa(1)); vMFrnd(n, mu(2,:), kappa(2)); vMFrnd(n, mu(3,:), kappa(3))];
    group = [ones(n,1); 2*ones(n,1); 3*ones(n,1)];
    sphereprojplot(Y, 'group', group, 'cap', struct('mu', mu, 'kappa', kappa));
%}
%
%{
    %% Output of tclustvmf in high dimension.
    % Two groups and uniform contamination in dimension p=50. The plot
    % uses the estimated mean directions and concentrations.
    rng(2)
    p = 50; n = 200;
    mu1 = [1 zeros(1,p-1)]; mu2 = [0 1 zeros(1,p-2)];
    Y = [vMFrnd(n, mu1, 150); vMFrnd(n, mu2, 80); vMFrnd(30, mu1, 0)];
    out = tclustvmf(Y, 2, 0.08, 50, 'msg', 0);
    [~, outproj] = sphereprojplot(Y, 'group', out.idx, ...
        'cap', struct('mu', out.muopt, 'kappa', out.kappaopt));
    fprintf('Squared length retained by the projection: %.1f%%\n', 100*outproj.retained);
%}
%
%{
    %% Comparison of the projections.
    % The same data in the subspace of the means, in the subspace of the
    % first 3 principal directions and in the subspace of 3 coordinates
    % chosen by the user.
    rng(3)
    p = 30; n = 150;
    mu = zeros(3,p); mu(1,1:2) = [1 1]; mu(2,2:3) = [1 -1]; mu(3,4) = 1;
    Y = [vMFrnd(n, mu(1,:), 100); vMFrnd(n, mu(2,:), 100); vMFrnd(n, mu(3,:), 100)];
    group = [ones(n,1); 2*ones(n,1); 3*ones(n,1)];
    figure
    sphereprojplot(Y, 'group', group, 'proj', 'means', 'h', subplot(1,3,1), 'legend', false);
    sphereprojplot(Y, 'group', group, 'proj', 'pca', 'h', subplot(1,3,2), 'legend', false);
    sphereprojplot(Y, 'group', group, 'proj', [1 2 4], 'h', subplot(1,3,3), 'legend', false);
%}
%
%{
    % New data in the subspace of a previous plot.
    % The matrix out.V of the first call is used as projection in the
    % second call.
    rng(4)
    p = 100; n = 200;
    mu = zeros(2,p); mu(1,1) = 1; mu(2,2) = 1;
    Y = [vMFrnd(n, mu(1,:), 300); vMFrnd(n, mu(2,:), 300)];
    group = [ones(n,1); 2*ones(n,1)];
    [~, out] = sphereprojplot(Y, 'group', group);
    Ynew = [vMFrnd(50, mu(1,:), 300); vMFrnd(50, mu(2,:), 300)];
    sphereprojplot(Ynew, 'group', [ones(50,1); 2*ones(50,1)], 'proj', out.V);
%}
%
%{
    % Data without groups: principal directions.
    rng(5)
    p = 20;
    Y = vMFrnd(500, ones(1,p), 60);
    [~, out] = sphereprojplot(Y);
    disp(out.type)
%}

%% Beginning of code

if nargin < 1 || isempty(Y)
    error('FSDA:sphereprojplot:missingInputs', 'Input data matrix is missing')
end
[n, p] = size(Y);
if p <= 3
    error('FSDA:sphereprojplot:WrongDim', ...
        'sphereprojplot is for data with more than 3 columns: use sphereplot for p=2 or p=3')
end

% Options of this function; all the others are passed to sphereplot.
% group and noisegroup are used here and also passed to sphereplot.
options = struct('proj', 'auto', 'cap', [], 'center', [], 'conflev', 0.95, ...
    'npts', 200, 'group', [], 'noisegroup', [], 'plottitle', true);
fwd = {};
if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
    if mod(length(varargin), 2) ~= 0
        error('FSDA:sphereprojplot:WrongInputOpt', ...
            'Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    for i = 1:2:length(varargin)
        name = varargin{i};
        if isfield(options, name)
            options.(name) = varargin{i+1};
            if any(strcmp(name, {'group','noisegroup'}))
                fwd = [fwd, varargin(i:i+1)]; %#ok<AGROW>
            end
        else
            fwd = [fwd, varargin(i:i+1)]; %#ok<AGROW>
        end
    end
end
proj       = options.proj;
cap        = options.cap;
center     = options.center;
conflev    = options.conflev(:)';
npts       = options.npts;
group      = options.group;
noisegroup = options.noisegroup;
plottitle  = options.plottitle;

if isempty(conflev) || any(conflev <= 0 | conflev >= 1)
    error('FSDA:sphereprojplot:WrongConflev', 'conflev must have elements in the interval (0,1)')
end

%% Regions in R^p: centres Cp (k-by-p) and thresholds thr (k-by-nlev)
Cp = zeros(0, p); thr = zeros(0, 1); okreg = false(0, 1);
if ~isempty(cap)
    if isstruct(cap) && isfield(cap, 'mu')
        Cp = cap.mu;
        if isvector(Cp)
            Cp = Cp(:)';
        end
        kap = cap.kappa(:);
        if size(Cp, 2) ~= p || numel(kap) ~= size(Cp, 1)
            error('FSDA:sphereprojplot:WrongCap', ...
                'cap.mu must be a k-by-p matrix and cap.kappa a vector of length k')
        end
        Cp = Cp./sqrt(sum(Cp.^2, 2));
        okreg = ~isnan(kap) & ~any(isnan(Cp), 2);
        thr = NaN(size(Cp, 1), numel(conflev));
        for j = 1:size(Cp, 1)
            if okreg(j)
                for l = 1:numel(conflev)
                    [~, thr(j,l)] = vMFcap(Cp(j,:), kap(j), conflev(l), 0);
                end
            end
        end
    elseif isstruct(cap) && isfield(cap, 'gamma3')
        error('FSDA:sphereprojplot:WrongCap', ...
            'Kent regions are defined only for p=3: use sphereplot')
    elseif iscell(cap)
        kc = numel(cap);
        okreg = false(kc, 1);
        for j = 1:kc
            okreg(j) = ~isempty(cap{j}) && size(cap{j}, 2) == p && ~any(isnan(cap{j}(:)));
        end
        if isempty(center)
            Cp = zeros(kc, p);
            for j = 1:kc
                if okreg(j)
                    m = mean(cap{j}, 1);
                    okreg(j) = norm(m) > 1e-8;
                    Cp(j,:) = m/max(norm(m), 1e-8);
                end
            end
        else
            Cp = center;
            if isvector(Cp)
                Cp = Cp(:)';
            end
            if size(Cp, 2) ~= p || size(Cp, 1) ~= kc
                error('FSDA:sphereprojplot:WrongCenter', 'center must be a k-by-p matrix')
            end
            okreg = okreg & ~any(isnan(Cp), 2);
            Cp = Cp./sqrt(sum(Cp.^2, 2));
        end
        % all the points of the boundary have the same inner product with
        % the centre: the threshold of the cap
        thr = NaN(kc, 1);
        for j = 1:kc
            if okreg(j)
                thr(j) = mean(cap{j}*Cp(j,:)');
            end
        end
    else
        error('FSDA:sphereprojplot:WrongCap', ...
            'cap must be a structure with fields mu and kappa or a cell array of boundaries')
    end
end

%% Mean directions of the groups (used by the projection 'means')
if ~isempty(group)
    group = double(group(:));
    if numel(group) ~= n
        error('FSDA:sphereprojplot:WrongGroup', 'group must have as many elements as the rows of Y')
    end
    ids = unique(group(~isnan(group)));
    ids = ids(ids ~= 0);
    if ~isempty(noisegroup)
        ids = ids(~ismember(ids, noisegroup));
    end
else
    group = ones(n, 1);
    ids = 1;
end
M = zeros(0, p);
for i = 1:numel(ids)
    R = sum(Y(group == ids(i), :), 1);
    if norm(R) > 1e-10
        M(end+1, :) = R/norm(R); %#ok<AGROW>
    end
end
if any(okreg)
    % the centres of the regions have priority
    Mdir = Cp(okreg, :);
else
    Mdir = M;
end

%% Projection matrix V (p-by-3, orthonormal columns)
if ischar(proj)
    type = lower(proj);
    if strcmp(type, 'auto')
        if size(Mdir, 1) >= 2 || any(okreg)
            type = 'means';
        else
            type = 'pca';
        end
    end
    switch type
        case 'means'
            if isempty(Mdir)
                error('FSDA:sphereprojplot:NoMeans', ...
                    'Projection ''means'' requires option group or option cap')
            end
            % principal directions of the set of mean directions
            [~, S, W] = svd(Mdir, 'econ');
            s = diag(S);
            r = min(sum(s > 1e-8*s(1)), 3);
            V = W(:, 1:r);
            if r < 3
                % completed with the principal directions of the part of
                % the data which is orthogonal to the means
                Yres = Y - (Y*V)*V';
                V = [V, principalDirections(Yres, 3-r)];
            end
            axlab = {'proj 1','proj 2','proj 3'};
        case 'pca'
            V = principalDirections(Y, 3);
            axlab = {'PC 1','PC 2','PC 3'};
        otherwise
            error('FSDA:sphereprojplot:WrongProj', ...
                'proj must be ''auto'', ''means'', ''pca'', a vector of 3 indices or a p-by-3 matrix')
    end
elseif isnumeric(proj) && isvector(proj) && numel(proj) == 3
    if any(proj < 1 | proj > p | proj ~= fix(proj)) || numel(unique(proj)) ~= 3
        error('FSDA:sphereprojplot:WrongProj', ...
            'As a vector, proj must contain 3 distinct integers between 1 and p')
    end
    V = zeros(p, 3);
    V(sub2ind([p 3], proj(:)', 1:3)) = 1;
    type = 'coordinates';
    axlab = {sprintf('dim %d',proj(1)), sprintf('dim %d',proj(2)), sprintf('dim %d',proj(3))};
elseif isnumeric(proj) && isequal(size(proj), [p 3])
    V = proj;
    if norm(V'*V - eye(3)) > 1e-8
        [V, ~] = qr(V, 0);
    end
    type = 'user';
    axlab = {'proj 1','proj 2','proj 3'};
else
    error('FSDA:sphereprojplot:WrongProj', ...
        'proj must be ''auto'', ''means'', ''pca'', a vector of 3 indices or a p-by-3 matrix')
end

%% Projected data
Yproj = Y*V;
retainedunit = sum(Yproj.^2, 2);
retained = mean(retainedunit(~isnan(retainedunit)));

%% Projected regions: circles {||x||=1, m'*x=t} (see More About)
caps = {}; cen = zeros(0, 3); cg = zeros(0, 1);
phi = linspace(0, 2*pi, npts)';
for j = 1:size(Cp, 1)
    for l = 1:size(thr, 2)
        Bjl = [];
        chat = [1 0 0];
        if okreg(j)
            m = Cp(j,:)*V;
            normm = norm(m);
            if normm > 1e-10 && abs(thr(j,l)) < normm
                chat = m/normm;
                ct = thr(j,l)/normm;
                [~, imin] = min(abs(chat));
                e = zeros(1, 3); e(imin) = 1;
                e1 = e - (e*chat')*chat;
                e1 = e1/norm(e1);
                e2 = cross(chat, e1);
                Bjl = ct*chat + sqrt(1 - ct^2)*(cos(phi)*e1 + sin(phi)*e2);
            end
        end
        caps{end+1, 1} = Bjl; %#ok<AGROW>
        cen(end+1, :) = chat; %#ok<AGROW>
        cg(end+1, 1) = j; %#ok<AGROW>
    end
end

%% Plot
if isempty(caps)
    hax = sphereplot(Yproj, fwd{:});
else
    hax = sphereplot(Yproj, 'cap', caps, 'center', cen, 'capgroup', cg, ...
        'capflat', true, fwd{:});
end
xlabel(hax, axlab{1}); ylabel(hax, axlab{2}); zlabel(hax, axlab{3});
if plottitle
    title(hax, {sprintf('Projection of p=%d dimensions on 3 (%s)', p, type), ...
        sprintf('%.1f%% of the squared length of the units is retained', 100*retained)});
end

out = struct;
out.V = V;
out.Yproj = Yproj;
out.retained = retained;
out.retainedunit = retainedunit;
out.type = type;

end

%% Local functions

function V = principalDirections(X, q)
% Eigenvectors of X'*X associated with the q largest eigenvalues
X = X(~any(isnan(X), 2), :);
G = X'*X;
G = (G + G')/2;
[E, D] = eig(G);
[~, ord] = sort(diag(D), 'descend');
V = E(:, ord(1:q));
end

%FScategory:VIS-Clu
