function hax = sphereplot(Y, varargin)
%sphereplot plots directional data on the unit circle or on the unit sphere
%
%<a href="matlab: docsearchFS('sphereplot')">Link to the help page for this function</a>
%
%   sphereplot displays the rows of the n-by-p matrix Y, which are unit
%   vectors, on the unit circle (p=2) or on the unit sphere (p=3). The
%   units can be coloured according to a grouping variable (for example
%   the classification produced by tclustvmf or tclustkent, where 0
%   denotes the trimmed units), and the confidence regions of the groups
%   can be superimposed: circular caps for von Mises-Fisher groups,
%   elliptical caps for Kent groups. The regions can be drawn as simple
%   contours or as shaded areas lying on the sphere. When p=3 the sphere
%   can be rotated with the mouse.
%
%   sphereplot is the plotting function of all the FSDA functions for
%   directional data: the random number generators vMFrnd and Kentrnd,
%   the functions vMFcap and Kentcap which compute the confidence
%   regions, and the clustering functions tclustvmf and tclustkent. For
%   data in more than 3 dimensions (e.g. embeddings) see sphereprojplot,
%   which projects the data on a 3-dimensional space and then calls
%   sphereplot.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-2 or n-by-3 matrix of unit
%               vectors (one per row). Y can be empty if only the
%               confidence regions have to be drawn.
%
%  Optional input arguments:
%
%       group : Grouping variable. Vector of length n of numeric labels.
%               Units with the same label have the same colour. Label 0
%               identifies trimmed units/outliers, which are shown as red
%               crosses. If the labels are positive integers, group j
%               has the same colour as the j-th confidence region. If
%               group is empty (default) all units belong to one group.
%                 Example - 'group',out.idx
%                 Data Types - double
%         cap : Confidence regions to superimpose. Structure or cell
%               array. Default is empty (no region). Three alternatives.
%               1) von Mises-Fisher caps: structure with fields:
%                  cap.mu = k-by-p matrix of mean directions.
%                  cap.kappa = vector of length k of concentrations.
%                  The caps are computed at level conflev.
%               2) Kent elliptical caps (p=3 only): structure with fields:
%                  cap.gamma3 = k-by-3 matrix of poles.
%                  cap.gamma1 = k-by-3 matrix of major axes.
%                  cap.kappa = vector of length k of concentrations.
%                  cap.beta = vector of length k of ovalness parameters.
%                  The regions are computed with Kentcap at level conflev.
%               3) Cell array of length k of boundaries already computed,
%                  for example out.ellipse of tclustvmf and tclustkent or
%                  the output of vMFcap and Kentcap: cap{j} is a matrix with
%                  p columns containing points of the boundary of region j
%                  (ordered along the contour if p=3; the two end points
%                  of the arc if p=2). See also options center and capgroup.

%                 Example - 'cap',struct('mu',[1 0 0],'kappa',20)
%                 Data Types - struct | cell
%      center : Centres of the regions. k-by-p matrix. Used only when cap
%               is a cell array: row j is the mean direction (pole) of
%               region j, for example out.muopt of tclustvmf or
%               out.gamma3 of tclustkent. If it is empty (default) the
%               centre of a region is taken as the direction of the mean
%               of its boundary points, which is correct for regions
%               smaller than a hemisphere (half circle if p=2).
%                 Example - 'center',out.muopt
%                 Data Types - double
%     conflev : Confidence level(s) of the regions. Scalar or vector with
%               elements in (0,1). Default 0.95. Used only when cap is a
%               structure. If it is a vector, the nested regions of all
%               the levels are drawn for each group.
%                 Example - 'conflev',[0.5 0.9 0.99]
%                 Data Types - double
%    capgroup : Group to which each region refers. Vector of positive
%               integers with one element for each region. Region j is
%               drawn with the colour of group capgroup(j). Used only
%               when cap is a cell array; the default is 1:k (region j
%               refers to group j). It is useful to draw several regions
%               for the same group.
%                 Example - 'capgroup',[1 1 2]
%                 Data Types - double
%     capfill : Shaded regions. Logical. If true (default) the regions
%               are shaded with the colour of the group (a patch lying on
%               the sphere if p=3, a circular sector if p=2); if false
%               only the contours are drawn.
%                 Example - 'capfill',false
%                 Data Types - logical
%    capalpha : Transparency of the shaded regions. Scalar in [0,1].
%               Default 0.25.
%                 Example - 'capalpha',0.5
%                 Data Types - double
%     capflat : Flat shading (p=3). Logical. If false (default) the
%               shaded area lies on the surface of the sphere. If true
%               the shaded area is the flat surface inside the sphere
%               delimited by the contour. This is the appropriate choice
%               when the points are projections which lie inside the
%               sphere, and it is the default of sphereprojplot.
%                 Example - 'capflat',true
%                 Data Types - logical
%  plotcenter : Show the centres of the regions. Logical. Default true:
%               a black cross on the sphere (p=3) or a radius of the
%               circle (p=2).
%                 Example - 'plotcenter',false
%                 Data Types - logical
%        npts : Number of points of the contours computed when cap is a
%               structure. Scalar. Default 200.
%                 Example - 'npts',400
%                 Data Types - double
%  showpoints : Show the units. Logical. Default true. If it is false
%               only the confidence regions are drawn (Y is still used to
%               choose the point of view) and the legend refers to the
%               regions.
%                 Example - 'showpoints',false
%                 Data Types - logical
%      labels : Names of the groups in the legend. Cell array of
%               character vectors or string array. labels{j} is the name
%               of group j (of region j if the units are not shown).
%               Default is 'Group j' ('Region j').
%                 Example - 'labels',{'Sport','Politics'}
%                 Data Types - cell | string
%         clr : Colours of the groups. Matrix with 3 columns of RGB
%               triplets; row j is the colour of group j and of its
%               regions (colours are recycled if necessary). Default is a
%               palette of 14 colours which does not contain red,
%               reserved for trimmed units.
%                 Example - 'clr',[0 0 1; 0 0.6 0]
%                 Data Types - double
%  markersize : Size of the markers. Scalar. Default 12.
%                 Example - 'markersize',20
%                 Data Types - double
%  noisegroup : Label of the units assigned to the uniform noise
%               component (k+1 in the output of tclustvmf and tclustkent
%               with option noise). These units are shown in gray. Default
%               is empty (no noise group).
%                 Example - 'noisegroup',3
%                 Data Types - double
% spherealpha : Transparency of the sphere (p=3). Scalar in [0,1].
%               Default 0.25. With 1 the sphere is opaque and the units on
%               the far side are hidden.
%                 Example - 'spherealpha',1
%                 Data Types - double
%        grid : Reference grid. Logical. Default false. If true,
%               meridians and parallels every 30 degrees are drawn on the
%               sphere (p=3), or the two axes inside the circle (p=2).
%                 Example - 'grid',true
%                 Data Types - logical
%      rotate : Interactive rotation (p=3). Logical. If true (default)
%               the 3D rotation mode of the figure is switched on, so
%               that the sphere can be rotated by dragging it with the
%               mouse.
%                 Example - 'rotate',false
%                 Data Types - logical
%      jitter : Radial displacement of the points when p=2. Scalar.
%               Default 0.05. Points on a circle overlap heavily: unit i
%               is therefore drawn at a distance from the origin between
%               1 and 1+jitter (the displacement is a deterministic
%               function of i; no random number is generated). Use 0 to
%               draw all the points exactly on the circle.
%                 Example - 'jitter',0
%                 Data Types - double
%      legend : Add the legend. Logical. Default is true if there is
%               more than one set of units (groups, trimmed units, noise)
%               or of regions, and false otherwise.
%                 Example - 'legend',false
%                 Data Types - logical
%           h : Axes where to draw. Axes handle. Default is empty: a new
%               figure is created.
%                 Example - 'h',subplot(1,2,1)
%                 Data Types - matlab.graphics.axis.Axes
%
%  Output:
%
%         hax : Axes handle. Handle. Handle of the axes which contain the plot.
%
% More About:
%
% SHADED REGIONS ON THE SPHERE
%
% A region is shaded with a surface which lies on the sphere: each point
% b of its boundary is joined to the centre c of the region with an arc
% of great circle, c*cos(s*theta)+u*sin(s*theta) for s in [0,1], where
% theta is the angle between c and b and u is the unit vector orthogonal
% to c in the plane of c and b. This requires the region to be
% star-shaped with respect to its centre, which is always true for vMF
% caps and for the regions of unimodal Kent distributions. With option
% capflat the boundary points are instead joined with straight segments
% to their mean, which gives the flat surface delimited by the contour.
%
% DIMENSION
%
% Only p=2 and p=3 are supported, because only in these cases the plot
% is an exact representation of the data. For p>3 use sphereprojplot.
%
% References:
%
%   Mardia, K.V. and Jupp, P.E. (2000), Directional Statistics, Wiley,
%   Chichester.
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
% See also: sphereprojplot, vMFrnd, Kentrnd, vMFcap, Kentcap, tclustvmf, tclustkent
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('sphereplot')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Data on the sphere with all default options.
    rng(1)
    Y = vMFrnd(300, [1 0 0], 20);
    sphereplot(Y);
%}
%
%{
    %% A circular (vMF) and an elliptical (Kent) group with their 95 per cent regions.
    % The boundaries are computed with vMFcap and Kentcap and passed as a
    % cell array.
    rng(2)
    mu = [1 0 0]; kappavmf = 40;
    g3 = [0 1 0]; g1 = [0 0 1]; kappa = 80; beta = 32;
    Y = [vMFrnd(300, mu, kappavmf); Kentrnd(300, kappa, beta, g3, g1)];
    group = [ones(300,1); 2*ones(300,1)];
    caps = {vMFcap(mu, kappavmf), Kentcap(g1, [], g3, kappa, beta)};
    sphereplot(Y, 'group', group, 'cap', caps, 'center', [mu; g3], ...
        'labels', {'von Mises-Fisher','Kent'});
%}
%
%{
    % Shaded regions versus contours only.
    rng(3)
    g3 = [1 0 0; 0 1 0]; g1 = [0 0 1; 1 0 0];
    kappa = [100; 50]; beta = [40; 15];
    Y = [Kentrnd(200, kappa(1), beta(1), g3(1,:), g1(1,:));
         Kentrnd(200, kappa(2), beta(2), g3(2,:), g1(2,:))];
    group = [ones(200,1); 2*ones(200,1)];
    cap = struct('gamma3', g3, 'gamma1', g1, 'kappa', kappa, 'beta', beta);
    figure
    sphereplot(Y, 'group', group, 'cap', cap, 'h', subplot(1,2,1));
    title('capfill=true')
    sphereplot(Y, 'group', group, 'cap', cap, 'capfill', false, 'h', subplot(1,2,2));
    title('capfill=false')
%}
%
%{
    %% Regions only: a concentrated and a dispersed vMF distribution.
    % No data: only the 95 per cent caps are drawn, and the legend
    % refers to the regions.
    cap = struct('mu', [1 0 0; 0 1 0], 'kappa', [40; 8]);
    sphereplot([], 'cap', cap, 'labels', {'kappa=40','kappa=8'});
%}
%
%{
    % Nested regions at several confidence levels.
    rng(4)
    g3 = [1 0 0]; g1 = [0 0 1]; kappa = 60; beta = 25;
    Y = Kentrnd(500, kappa, beta, g3, g1);
    cap = struct('gamma3', g3, 'gamma1', g1, 'kappa', kappa, 'beta', beta);
    sphereplot(Y, 'cap', cap, 'conflev', [0.5 0.9 0.99], 'capalpha', 0.15, 'markersize', 8);
    title('Kent regions at levels 0.5, 0.9 and 0.99')
%}
%
%{
    % From a circular to an elongated cap: effect of beta.
    % kappa is fixed; the cap is circular for beta=0 and becomes more
    % elongated as beta approaches kappa/2.
    rng(5)
    kappa = 60; betas = [0 15 25];
    g3 = [1 0 0]; g1 = [0 0 1];
    figure
    for i = 1:3
        Y = Kentrnd(400, kappa, betas(i), g3, g1);
        cap = struct('gamma3', g3, 'gamma1', g1, 'kappa', kappa, 'beta', betas(i));
        sphereplot(Y, 'cap', cap, 'markersize', 6, 'h', subplot(1,3,i));
        title(sprintf('\\kappa=%g, \\beta=%g', kappa, betas(i)))
    end
%}
%
%{
    %% Units inside and outside a fitted region.
    % A Kent distribution is fitted to the data; the units outside the
    % fitted 95 per cent region are highlighted.
    rng(6)
    n = 500;
    Y = Kentrnd(n, 80, 30, [1 0 0], [0 1 1]);
    [g1,g2,g3,kap,bet] = KentMstepFS(Y, ones(n,1));
    [P, thresh] = Kentcap(g1, g2, g3, kap, bet, 0.95);
    inside = logKentpdfFS(Y, g1, g2, g3, kap, bet) >= thresh;
    sphereplot(Y, 'group', 2 - inside, 'cap', {P}, 'center', g3, ...
        'labels', {'inside','outside'});
    title(sprintf('%.1f%% of the units inside the fitted 95%% region', 100*mean(inside)))
%}
%
%{
    % Output of tclustkent: groups, trimmed units and noise component.
    % Two elliptical groups and uniform contamination. Units assigned to
    % the noise component (label 3) are in gray, trimmed units (label 0)
    % are red crosses.
    rng(7)
    Y = [Kentrnd(200, 100, 40, [1 0 0], [0 0 1]);
         Kentrnd(200,  50, 15, [0 1 0], [1 0 0]); vMFrnd(60, [0 0 1], 0)];
    out = tclustkent(Y, 2, 0.02, 20, 'noise', true, 'msg', 0);
    sphereplot(Y, 'group', out.idx, 'noisegroup', 3, ...
        'cap', out.ellipse, 'center', out.gamma3);
    title('tclustkent with noise component')
%}
%
%{
    % Output of tclustvmf.
    % The caps are computed from the estimated parameters at level 0.9.
    rng(8)
    Y = [vMFrnd(200, [1 0 0], 40); vMFrnd(200, [0 1 1], 15); vMFrnd(30, [0 0 1], 0)];
    out = tclustvmf(Y, 2, 0.08, 50, 'msg', 0);
    sphereplot(Y, 'group', out.idx, 'conflev', 0.9, ...
        'cap', struct('mu', out.muopt, 'kappa', out.kappaopt));
%}
%
%{
    % Data on the circle (p=2).
    % Two von Mises groups with their 95 per cent arcs.
    rng(9)
    mu = [1 0; -1 1]; kappa = [10; 4];
    Y = [vMFrnd(200, mu(1,:), kappa(1)); vMFrnd(200, mu(2,:), kappa(2))];
    group = [ones(200,1); 2*ones(200,1)];
    sphereplot(Y, 'group', group, 'cap', struct('mu', mu, 'kappa', kappa));
%}
%
%{
    % Opaque sphere with meridians and parallels, personalised colours.
    rng(10)
    g3 = [1 1 0]; g1 = [0 0 1]; kappa = 60; beta = 25;
    Y = Kentrnd(400, kappa, beta, g3, g1);
    cap = struct('gamma3', g3, 'gamma1', g1, 'kappa', kappa, 'beta', beta);
    sphereplot(Y, 'cap', cap, 'clr', [0 0 0.6], 'spherealpha', 1, ...
        'grid', true, 'capalpha', 0.15, 'markersize', 8);
%}

%% Beginning of code

if nargin < 1
    error('FSDA:sphereplot:missingInputs', 'Input data matrix is missing')
end

% Default palette: MATLAB default axes colours (dark red replaced by dark
% gray) plus some additional well separated colours. Red is reserved for
% trimmed units.
clrdef = [0.0000 0.4470 0.7410
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250
    0.4940 0.1840 0.5560
    0.4660 0.6740 0.1880
    0.3010 0.7450 0.9330
    0.2500 0.2500 0.2500
    0.0000 0.5000 0.5000
    0.6000 0.4000 0.2000
    0.5000 0.0000 0.5000
    0.5000 0.5000 0.0000
    0.0000 0.0000 0.5000
    1.0000 0.6000 0.7840
    0.0000 0.5000 0.0000];

options = struct('group', [], 'cap', [], 'center', [], 'conflev', 0.95, ...
    'capgroup', [], 'capfill', true, 'capalpha', 0.25, 'capflat', false, ...
    'plotcenter', true, 'npts', 200, 'showpoints', true, 'labels', [], ...
    'clr', clrdef, 'markersize', 12, 'noisegroup', [], 'spherealpha', 0.25, ...
    'grid', false, 'rotate', true, 'jitter', 0.05, 'legend', [], 'h', []);

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:sphereplot:WrongInputOpt', ...
            'Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    inpchk = isfield(options, UserOptions);
    WrongOptions = UserOptions(inpchk == 0);
    if ~isempty(WrongOptions)
        disp(strcat('Non existent user option found->', char(WrongOptions{:})))
        error('FSDA:sphereplot:NonExistInputOpt', ...
            'In total %d non-existent user options found.', length(WrongOptions));
    end
    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
end

group       = options.group;
cap         = options.cap;
center      = options.center;
conflev     = options.conflev(:)';
capgroup    = options.capgroup;
capfill     = options.capfill;
capalpha    = options.capalpha;
capflat     = options.capflat;
plotcenter  = options.plotcenter;
npts        = options.npts;
showpoints  = options.showpoints;
labels      = options.labels;
clr         = options.clr;
markersize  = options.markersize;
noisegroup  = options.noisegroup;
spherealpha = options.spherealpha;
gridOpt     = options.grid;
rotateOpt   = options.rotate;
jitter      = options.jitter;
legendOpt   = options.legend;
hax         = options.h;

if ischar(labels)
    labels = {labels};
end

%% Dimension
[n, p] = size(Y);
if n == 0
    % no data: the dimension is taken from the regions
    p = dimFromCap(cap, center);
end
if p ~= 2 && p ~= 3
    error('FSDA:sphereplot:WrongDim', ...
        'Y must have 2 or 3 columns. For data in higher dimension use sphereprojplot.')
end
if size(clr, 2) ~= 3
    error('FSDA:sphereplot:WrongClr', 'clr must be a matrix with 3 columns (RGB triplets)')
end
if isempty(conflev) || any(conflev <= 0 | conflev >= 1)
    error('FSDA:sphereplot:WrongConflev', 'conflev must have elements in the interval (0,1)')
end
nlev = numel(conflev);
colfun = @(j) clr(mod(j-1, size(clr,1)) + 1, :);
outlierColor = [1 0 0];
noiseColor   = [0.6 0.6 0.6];

%% Groups
if isempty(group)
    group = ones(n, 1);
    nogroup = true;
else
    group = double(group(:));
    nogroup = false;
    if numel(group) ~= n
        error('FSDA:sphereplot:WrongGroup', 'group must have as many elements as the rows of Y')
    end
end
ids = unique(group(~isnan(group)));
isspecial = ids == 0;
if ~isempty(noisegroup)
    isspecial = isspecial | ismember(ids, noisegroup);
end
idsreg = ids(~isspecial);
% index of the colour of each regular group
if all(idsreg > 0 & idsreg == fix(idsreg))
    cidx = idsreg;
else
    cidx = (1:numel(idsreg))';
end

%% Confidence regions
% For region j: B{j} = boundary (p=3); ang(j) = half width of the arc
% (p=2); C(j,:) = centre; cgrp(j) = group (colour) to which it refers;
% valid(j) = false for regions which must be skipped
B = cell(0, 1); C = zeros(0, p); ang = zeros(0, 1);
valid = false(0, 1); cgrp = zeros(0, 1);
if ~isempty(cap)
    if isstruct(cap) && (isfield(cap, 'mu') || isfield(cap, 'gamma3'))
        iskent = ~isfield(cap, 'mu');
        if iskent
            % Kent elliptical caps
            if p ~= 3
                error('FSDA:sphereplot:WrongCap', 'Kent regions are defined only for p=3')
            end
            if ~all(isfield(cap, {'gamma1','kappa','beta'}))
                error('FSDA:sphereplot:WrongCap', ...
                    'For Kent regions cap must have fields gamma3, gamma1, kappa and beta')
            end
            C0 = rowsOrVector(cap.gamma3, p);
            G1 = rowsOrVector(cap.gamma1, p);
            bet = cap.beta(:);
        else
            % von Mises-Fisher caps
            C0 = rowsOrVector(cap.mu, p);
            G1 = zeros(size(C0));
            bet = zeros(size(C0, 1), 1);
        end
        k0 = size(C0, 1);
        kap = cap.kappa(:);
        if size(G1,1) ~= k0 || numel(kap) ~= k0 || numel(bet) ~= k0
            error('FSDA:sphereplot:WrongCap', ...
                'The fields of cap must refer to the same number of regions')
        end
        C0 = C0./sqrt(sum(C0.^2, 2));
        ok0 = ~isnan(kap) & ~isnan(bet) & ~any(isnan(C0), 2) & ~any(isnan(G1), 2);
        for j = 1:k0
            for l = 1:nlev
                Bjl = []; angjl = 0;
                if ok0(j)
                    if iskent
                        Bjl = Kentcap(G1(j,:), [], C0(j,:), kap(j), bet(j), conflev(l), npts);
                    elseif p == 3
                        Bjl = vMFcap(C0(j,:), kap(j), conflev(l), npts);
                    else
                        % half width of the arc
                        [~, tjl] = vMFcap(C0(j,:), kap(j), conflev(l), 0);
                        angjl = acos(tjl);
                    end
                end
                B{end+1, 1} = Bjl; %#ok<AGROW>
                C(end+1, :) = C0(j,:); %#ok<AGROW>
                ang(end+1, 1) = angjl; %#ok<AGROW>
                valid(end+1, 1) = ok0(j); %#ok<AGROW>
                cgrp(end+1, 1) = j; %#ok<AGROW>
            end
        end
    elseif iscell(cap)
        % boundaries supplied by the user
        B = cap(:);
        kc = numel(B);
        valid = false(kc, 1);
        for j = 1:kc
            valid(j) = ~isempty(B{j}) && size(B{j}, 2) == p && ~any(isnan(B{j}(:)));
        end
        if isempty(center)
            C = zeros(kc, p);
            for j = 1:kc
                if valid(j)
                    m = mean(unique(round(B{j}, 10), 'rows'), 1);
                    if norm(m) < 1e-8
                        error('FSDA:sphereplot:CenterRequired', ...
                            'The centre of region %d cannot be found from its boundary: use option ''center''', j)
                    end
                    C(j,:) = m/norm(m);
                end
            end
        else
            C = rowsOrVector(center, p);
            if size(C, 1) ~= kc
                error('FSDA:sphereplot:WrongCenter', 'center must have one row for each region')
            end
            valid = valid & ~any(isnan(C), 2);
            C = C./sqrt(sum(C.^2, 2));
        end
        ang = zeros(kc, 1);
        if p == 2
            for j = 1:kc
                if valid(j)
                    ang(j) = max(acos(min(max(B{j}*C(j,:)', -1), 1)));
                end
            end
        end
        if isempty(capgroup)
            cgrp = (1:kc)';
        else
            cgrp = capgroup(:);
            if numel(cgrp) ~= kc || any(cgrp < 1 | cgrp ~= fix(cgrp))
                error('FSDA:sphereplot:WrongCapgroup', ...
                    'capgroup must contain a positive integer for each region')
            end
        end
    else
        error('FSDA:sphereplot:WrongCap', ...
            ['cap must be a structure with fields mu and kappa (vMF), a structure with fields ' ...
            'gamma3, gamma1, kappa and beta (Kent) or a cell array of boundaries'])
    end
end
kc = numel(valid);

% The regions are described in the legend when the units are not shown
showpoints = showpoints && n > 0;
reglegend = ~showpoints;
reglabelled = false(max([cgrp; 0]), 1);

%% Axes
if isempty(hax)
    figure;
    hax = axes;
end
washold = ishold(hax);
nsets = 0;

if p == 3
    %% Sphere
    [sx, sy, sz] = sphere(60);
    surf(hax, sx, sy, sz, 'FaceColor', [0.9 0.9 1], 'EdgeColor', 'none', ...
        'FaceAlpha', spherealpha, 'HandleVisibility', 'off');
    hold(hax, 'on')

    if gridOpt
        t = linspace(0, 2*pi, 181)';
        gridColor = [0.65 0.65 0.75];
        for lon = 0:30:150
            plot3(hax, cos(t)*cosd(lon), cos(t)*sind(lon), sin(t), '-', ...
                'Color', gridColor, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
        for lat = -60:30:60
            plot3(hax, cosd(lat)*cos(t), cosd(lat)*sin(t), sind(lat)*ones(size(t)), '-', ...
                'Color', gridColor, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end

    % Confidence regions
    for j = 1:kc
        if ~valid(j)
            continue
        end
        col = colfun(cgrp(j));
        Bj = B{j};
        cj = C(j,:);
        % closed contour
        if norm(Bj(1,:) - Bj(end,:)) > 1e-8
            Bj = [Bj; Bj(1,:)]; %#ok<AGROW>
        end
        if capfill
            s = linspace(0, 1, 25);
            if capflat
                % straight segments from the mean of the contour
                c0 = mean(Bj(1:end-1,:), 1);
                Xc = c0(1) + (Bj(:,1) - c0(1))*s;
                Yc = c0(2) + (Bj(:,2) - c0(2))*s;
                Zc = c0(3) + (Bj(:,3) - c0(3))*s;
            else
                % arcs of great circle from the centre to the contour
                ct = min(max(Bj*cj', -1), 1);
                th = acos(ct);
                U = Bj - ct*cj;
                normU = sqrt(sum(U.^2, 2));
                normU(normU < 1e-12) = 1;
                U = U./normU;
                % slightly outside the sphere to avoid overlapping surfaces
                rho = 1.002;
                Xc = rho*(cos(th*s)*cj(1) + sin(th*s).*U(:,1));
                Yc = rho*(cos(th*s)*cj(2) + sin(th*s).*U(:,2));
                Zc = rho*(cos(th*s)*cj(3) + sin(th*s).*U(:,3));
            end
            surf(hax, Xc, Yc, Zc, 'FaceColor', col, 'FaceAlpha', capalpha, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hline = plot3(hax, Bj(:,1), Bj(:,2), Bj(:,3), '-', 'Color', col, ...
            'LineWidth', 2, 'HandleVisibility', 'off');
        if reglegend && ~reglabelled(cgrp(j))
            set(hline, 'HandleVisibility', 'on', 'DisplayName', ...
                setLabel(labels, cgrp(j), sprintf('Region %d', cgrp(j))));
            reglabelled(cgrp(j)) = true;
            nsets = nsets + 1;
        end
        if plotcenter
            plot3(hax, 1.004*cj(1), 1.004*cj(2), 1.004*cj(3), '+', 'Color', 'k', ...
                'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end

    % Units
    if showpoints
        for i = 1:numel(idsreg)
            sel = group == idsreg(i);
            if nogroup
                dname = setLabel(labels, 1, 'Data');
            else
                dname = setLabel(labels, cidx(i), sprintf('Group %g', idsreg(i)));
            end
            scatter3(hax, Y(sel,1), Y(sel,2), Y(sel,3), markersize, 'o', ...
                'MarkerFaceColor', colfun(cidx(i)), 'MarkerEdgeColor', 'none', 'DisplayName', dname);
            nsets = nsets + 1;
        end
        if ~isempty(noisegroup)
            sel = ismember(group, noisegroup);
            if any(sel)
                scatter3(hax, Y(sel,1), Y(sel,2), Y(sel,3), markersize, 'o', ...
                    'MarkerFaceColor', noiseColor, 'MarkerEdgeColor', 'none', 'DisplayName', 'Noise');
                nsets = nsets + 1;
            end
        end
        sel = group == 0;
        if any(sel)
            scatter3(hax, Y(sel,1), Y(sel,2), Y(sel,3), 2*markersize, 'x', ...
                'MarkerEdgeColor', outlierColor, 'LineWidth', 1, 'DisplayName', 'Trimmed');
            nsets = nsets + 1;
        end
    end

    axis(hax, 'equal')
    axis(hax, [-1.05 1.05 -1.05 1.05 -1.05 1.05])
    xlabel(hax, 'X'); ylabel(hax, 'Y'); zlabel(hax, 'Z');
    % point of view: along the mean direction of what is plotted
    if n > 0
        R = sum(Y(~any(isnan(Y), 2), :), 1);
    elseif any(valid)
        R = sum(C(valid,:), 1);
    else
        R = zeros(1, 3);
    end
    if norm(R) > 1e-6*max(n, 1)
        view(hax, R/norm(R));
    else
        view(hax, 3);
    end
    axis(hax, 'vis3d')
    if rotateOpt
        % the sphere can be rotated with the mouse
        try
            rotate3d(hax, 'on');
        catch
            % rotation mode not available (e.g. no display)
        end
    end

else
    %% Circle
    t = linspace(0, 2*pi, 361);
    plot(hax, cos(t), sin(t), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
        'HandleVisibility', 'off');
    hold(hax, 'on')
    if gridOpt
        plot(hax, [-1 1], [0 0], ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
        plot(hax, [0 0], [-1 1], ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    end

    % Confidence arcs
    for j = 1:kc
        if ~valid(j)
            continue
        end
        col = colfun(cgrp(j));
        a0 = atan2(C(j,2), C(j,1));
        tt = linspace(a0 - ang(j), a0 + ang(j), 200);
        if capfill
            patch(hax, [0 cos(tt)], [0 sin(tt)], col, 'FaceAlpha', capalpha, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hline = plot(hax, cos(tt), sin(tt), '-', 'Color', col, 'LineWidth', 3, ...
            'HandleVisibility', 'off');
        if reglegend && ~reglabelled(cgrp(j))
            set(hline, 'HandleVisibility', 'on', 'DisplayName', ...
                setLabel(labels, cgrp(j), sprintf('Region %d', cgrp(j))));
            reglabelled(cgrp(j)) = true;
            nsets = nsets + 1;
        end
        if plotcenter
            plot(hax, [0 C(j,1)], [0 C(j,2)], '-', 'Color', col, 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
    end

    % Units: deterministic radial displacement in [1, 1+jitter]
    if showpoints
        rad = 1 + jitter*mod((1:n)'*0.6180339887, 1);
        for i = 1:numel(idsreg)
            sel = group == idsreg(i);
            if nogroup
                dname = setLabel(labels, 1, 'Data');
            else
                dname = setLabel(labels, cidx(i), sprintf('Group %g', idsreg(i)));
            end
            scatter(hax, rad(sel).*Y(sel,1), rad(sel).*Y(sel,2), markersize, 'o', ...
                'MarkerFaceColor', colfun(cidx(i)), 'MarkerEdgeColor', 'none', 'DisplayName', dname);
            nsets = nsets + 1;
        end
        if ~isempty(noisegroup)
            sel = ismember(group, noisegroup);
            if any(sel)
                scatter(hax, rad(sel).*Y(sel,1), rad(sel).*Y(sel,2), markersize, 'o', ...
                    'MarkerFaceColor', noiseColor, 'MarkerEdgeColor', 'none', 'DisplayName', 'Noise');
                nsets = nsets + 1;
            end
        end
        sel = group == 0;
        if any(sel)
            scatter(hax, rad(sel).*Y(sel,1), rad(sel).*Y(sel,2), 2*markersize, 'x', ...
                'MarkerEdgeColor', outlierColor, 'LineWidth', 1, 'DisplayName', 'Trimmed');
            nsets = nsets + 1;
        end
    end

    axis(hax, 'equal')
    lim = 1.1 + jitter;
    axis(hax, [-lim lim -lim lim])
    xlabel(hax, 'X'); ylabel(hax, 'Y');
end

%% Legend
if isempty(legendOpt)
    legendOpt = nsets > 1;
end
if legendOpt && nsets > 0
    legend(hax, 'show', 'Location', 'best');
end

if ~washold
    hold(hax, 'off')
end

end

%% Local functions

function M = rowsOrVector(M, p)
% A vector of length p is a single direction (one row)
if isvector(M) && numel(M) == p
    M = M(:)';
end
if size(M, 2) ~= p
    error('FSDA:sphereplot:WrongDim', ...
        'Directions of the confidence regions must have %d columns, as Y', p)
end
end

function p = dimFromCap(cap, center)
% Dimension of the space when Y is empty
p = 0;
if isstruct(cap) && isfield(cap, 'mu')
    if isvector(cap.mu)
        p = numel(cap.mu);
    else
        p = size(cap.mu, 2);
    end
elseif isstruct(cap) && isfield(cap, 'gamma3')
    p = 3;
elseif iscell(cap) && ~isempty(cap)
    for j = 1:numel(cap)
        if ~isempty(cap{j})
            p = size(cap{j}, 2);
            break
        end
    end
elseif ~isempty(center)
    p = size(center, 2);
end
end

function s = setLabel(labels, idx, default)
% Name of group/region idx: labels{idx} if supplied, default otherwise
s = default;
if iscell(labels) && idx >= 1 && idx <= numel(labels) && idx == fix(idx) ...
        && ~isempty(labels{idx})
    s = char(labels{idx});
end
end

%FScategory:VIS-Clu
