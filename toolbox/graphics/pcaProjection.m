function pcaProjection(Y,varargin)
%pcaProjection shows projections of 3D points in PCA
%
%<a href="matlab: docsearchFS('pcaProjection')">Link to the help function</a>
%
%   This function is intended for didactic purposes. It illustrates the
%   details of the orthogonal projection of trivariate points onto the space of the
%   first 3 principal components and displays the principal ellipsoid in both
%   the original and projected spaces.
%
%  Required input arguments:
%
% Y :           Input data. 2D array or table.
%               n x p data matrix; n observations and p variables. Rows of
%               Y represent observations, and columns represent variables.
%                Data Types - single|double
%
%  Optional input arguments:
%
%
%   DataVars :  Table variables for which to compute PCA and projection of points.
%               Numeric vector or string array of cell array of character
%               vectors of length 3. For examples if 'DataVars' is [3 5 6],
%               PCA is done for variables 3, 5 and 6. If 'DataVars' is
%               ["Name1" "Name4" "Name8"] variable with these names inside the
%               input table are used. Note that if DataVars is not
%               specified the first 3 variables are used.
%               Example - 'DataVars',[2 4 5]
%               Data Types -  character vector | string array | cell array of character vectors | vector of positive integers | logical vector
%
%
%
%     conflev   :  confidence level of the ellipsoid. Scalar. A number in the interval (0 1)
%                  which specifies the confidence level of the ellipsoid which encloses
%                  the points. The confidence level is computed under normality assumption
%                   and it is based on the chi2 with 3 degrees of freedom.
%                   Example - 'conflev',0.99
%                   Data Types - single, double
%
%      AddAxes   : show the axes. Boolean vector of length 3. Boolean vector which
%                  specifies whether to show or hide the PC axis in the original
%                  or transformed space. For example if AddAxes [true false false]
%                  just the first principal axis is shown. The default of
%                  addAxes is true(3,1) that is all the 3 axes are shown.
%                   Example - 'AddAxes',[false false true]
%                   Data Types - logical
%
% AddConstantPlane: show the constant plane spanned by the first two PCs.
%                   Boolean scalar. If this option is true (default) the constant plane
%                   spanned by the first two PCs is shown on the screen in the plot of
%                   projected points in the space of the first two PCs.
%                   Example - 'AddConstantPlane',false
%                   Data Types - logical
%
% LineWidthAxes   : line width of the principal axes. Scalar.
%                   Width of the lines associated with the PC axis (when
%                   they are visible). The default value of LineWidthAxes
%                   is 3.
%                   Example - 'LineWidthAxes',1
%                   Data Types - single, double
%
%    standardize : standardize data. boolean. Boolean which specifies
%               whether to standardize the variables, that is we operate on
%               the correlation matrix (default) or simply remove column
%               means (in this last case we operate on the covariance
%               matrix).
%                   Example - 'standardize',false
%                   Data Types - boolean
%
%TextDensityPercentage:  Percentage of text data to show. 
%                        Scalar from 0through 100. Percentage of text data
%                        to show, specified as a scalar from 0 through 100.
%                        To show all text, set TextDensityPercentage to
%                        100. To show no text, set TextDensityPercentage to
%                        0. The default value of TextDensityPercentage is
%                        60. Note that Text Analytic toolbox is required
%                   Example - 'TextDensityPercentage',80
%                   Data Types - single or double
%
%
%
% Output:
%
%
%
%
% See also: pcaFS, biplotFS
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('pcaProjection')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    %% Call to pcaProjection with all default options.
    load citiesItaly.mat
    pcaProjection(citiesItaly)
%}

%{
    % Call to pcaProjection with option DataVars.
    load citiesItaly.mat
    % Select variables in position 2 6 and 7.
    pcaProjection(citiesItaly,'DataVars',[2 6 7])
%}

%{
    % Call to pcaProjection with option TextDensityPercentage.
    % If TextDensityPercentage=100 then all labels are shown
    load citiesItaly.mat
    pcaProjection(citiesItaly,'DataVars',[2 6 7],'TextDensityPercentage',0)
%}


%{
    % Call to pcaProjection with option conflev.
    % conflev controls the confidence level for the ellipsoid
    load citiesItaly.mat
    pcaProjection(citiesItaly,'DataVars',[1 2 7],'conflev',0.5)
%}


%{
    % Call to pcaProjection with option AddAxes.
    % conflev controls the confidence level for the ellipsoid
    close all
    load citiesItaly.mat
    % Remove all axes from plots
    AddAxes=[false false false];
    pcaProjection(citiesItaly(:,[ 3 5 6]),'AddAxes',AddAxes)
%}


%{
    % Call to pcaProjection with option addPCaxes and LineWidthAxes.
    close all
    load citiesItaly.mat
    % Just show first principal axe
    AddAxes=[true false false];
    % Line width of the axes
    l=4;
    pcaProjection(citiesItaly(:,[ 3 5 6]),'AddAxes',AddAxes,'LineWidthAxes',l)
%}

%{
    %% Call to pcaProjection with option standardize.
    load citiesItaly.mat
    % This is the effect of non standardized the data when
    % standardization is needed!
    pcaProjection(citiesItaly(:,[ 1 2 5]),'standardize',false)
%}

%% Beginning of code
[n,v]=size(Y);
standardize=true;
conflev=0.95;
DataVars=1:3;
AddAxes=true(3,1);
AddConstantPlane=true;
LineWidthAxes=3;
TextDensityPercentage =60;

if nargin>1
    options=struct('standardize',standardize,...
        'conflev',conflev,'AddAxes',AddAxes,'LineWidthAxes',LineWidthAxes,...
        'AddConstantPlane', AddConstantPlane, ...
        'DataVars',DataVars,'TextDensityPercentage',TextDensityPercentage);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)


        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:pcaProjection:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end

        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:pcaProjection:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end


    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    standardize=options.standardize;
    conflev=options.conflev;
    DataVars=options.DataVars;
    TextDensityPercentage=options.TextDensityPercentage;
    AddConstantPlane=options.AddConstantPlane;
    AddAxes=options.AddAxes;
    LineWidthAxes=options.LineWidthAxes;

    if length(AddAxes)<3
        AddAxes=[AddAxes(:); false; false];
    end
end

if istable(Y)
    Y=Y(:,DataVars);
    varnames=Y.Properties.VariableNames;
    rownames=Y.Properties.RowNames;
    if isempty(rownames)
        rownames=cellstr(num2str((1:n)','%d'));
    end
    Y=table2array(Y);
else
    if isstring(DataVars)
        error('FSDA:pcaProjection:WngInput','Input is not table and VarNames have been supplied')
    end
    Y=Y(:,DataVars);
    varnames=cellstr(num2str((1:v)','Y%d'));
    rownames=cellstr(num2str((1:n)','%d'));
end


% X = array on which to perform projection
X=Y;

% meaX = row vector of arithmetic means
meaX = mean(X);
% Xtilde = matrix of deviations from the mean
Xtilde = X - meaX;
if standardize==true
    Xtilde = Xtilde ./ std(Xtilde);
end

% S = covariance (correlation) matrix calculated in matrix form
S = Xtilde'*Xtilde/(n-1);

% Eigenvalues and eigenvectors of S
[Vini,Lambdaini] = eig(S);

[~,ord] = sort(diag(Lambdaini),'descend');
% La is the  vector with the sorted eigenvalues
La = diag(Lambdaini(ord,ord));
Gam=sqrt(La);
% V contains corresponding eigenvectors
V = Vini(:,ord);

%% Scatter 3D with  principal line
% Principal line = line associated with the direction of maximum variability

% 3D scatter plot of the data
scatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3))
hold('on')


hold on
for j=1:3
    if AddAxes(j) == true
        addLinePCj(V,Xtilde,j,LineWidthAxes);
    end
end


xlabel(varnames(1))
ylabel(varnames(2))
zlabel(varnames(3))
axis equal
title("3D scatter plot with the line associated with first PC")

textscatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3),rownames,'TextDensityPercentage',TextDensityPercentage)
xlabel(varnames(1))


%% New figure with projections into the line of the first PC
figure
v1=V(:,1);
v2=V(:,2);
v3=V(:,3);
Xhat=Xtilde*(v1*v1') ;

scatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3))
hold('on')
plot3([Xhat(:,1) Xtilde(:,1)]', [Xhat(:,2) Xtilde(:,2)]', [Xhat(:,3) Xtilde(:,3)]', 'r')

for j=1:3
    if AddAxes(j) == true
        addLinePCj(V,Xtilde,1,LineWidthAxes)
    end
end

title("3D scatter plot with main line and orthogonal projections")

textscatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3),rownames,'TextDensityPercentage',TextDensityPercentage)
xlabel(varnames(1))
ylabel(varnames(2))
zlabel(varnames(3))
axis equal


%% 3D scatter with projections onto first two PCs plane
figure

Xhat2=Xtilde*((v1*v1')+(v2*v2')) ;

scatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3))
hold('on')
for j=1:3
    if AddAxes(j) == true
        addLinePCj(V,Xtilde,j,LineWidthAxes);
    end
end

plot3([Xhat2(:,1) Xtilde(:,1)]', [Xhat2(:,2) Xtilde(:,2)]', [Xhat2(:,3) Xtilde(:,3)]', 'r')
axis equal
title("3D scatter plot: orthogonal projections onto first 2PC plane of first 2 PCs")

if AddConstantPlane == true
    constantplane(V(:,3),0)
end

textscatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3),rownames, ...
    'TextDensityPercentage',TextDensityPercentage);
xlabel(varnames(1))
ylabel(varnames(2))
zlabel(varnames(3))

%% Scatter 3d in the space of the first 3 PCs

figure
scatter3(Xtilde*v1,Xtilde*v2,Xtilde*v3)
hold('on')
axis equal
constantplane('z',0);

title('3D scatter plot of projections onto the first three principal components');

textscatter3(Xtilde*v1,Xtilde*v2,Xtilde*v3,rownames,'TextDensityPercentage',TextDensityPercentage);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

color='k';

for j=1:3
    if AddAxes(j) == true
        add3Daxes(j,color,LineWidthAxes);
    end
end

%% Ellpsoid in the space of the PCs
figure
scatter3(Xtilde*v1,Xtilde*v2,Xtilde*v3)


hold('on')
axis equal

facProp=sqrt(chi2inv(conflev,3));
[xr,yr,zr]=ellipsoid(0,0,0,facProp*Gam(1),facProp*Gam(2),facProp*Gam(3));

surf(xr,yr,zr,'FaceColor','none')
textscatter3(Xtilde*v1,Xtilde*v2,Xtilde*v3,rownames,'TextDensityPercentage',TextDensityPercentage);


for j=1:3
    if AddAxes(j) == true
        add3Daxes(j,color,LineWidthAxes);
    end
end

xlabel('PC1'); ylabel('PC2'); zlabel('PC3');


%% Ellipsoid in the original space
figure
scatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3))
hold('on')

xr1=xr(:);
yr1=yr(:);
zr1=zr(:);
XX=[xr1 yr1 zr1];
XX1=XX*V';
nfac=sqrt(size(XX1,1));
surf(reshape(XX1(:,1),nfac,nfac),reshape(XX1(:,2),nfac,nfac),reshape(XX1(:,3),nfac,nfac),'FaceAlpha',0.1,'FaceColor','none');
textscatter3(Xtilde(:,1),Xtilde(:,2),Xtilde(:,3),rownames,'TextDensityPercentage',TextDensityPercentage);

for j=1:3
    if AddAxes(j) == true
        addLinePCj(V,XX1,j,LineWidthAxes);
    end
end

axis equal
xlabel(varnames(1))
ylabel(varnames(2))
zlabel(varnames(3))

end

function addLinePCj(V,Xtilde,j,lwd)
vj=V(:,j);
Xhatj = Xtilde*(vj*vj');
% Take two points to draw this line
[~,indminXj] = min(Xhatj(:,2));
[~,indmaxXj] = max(Xhatj(:,2));

% Add the line associated with jth PC
line([Xhatj(indminXj,1); Xhatj(indmaxXj,1)], [Xhatj(indminXj,2); Xhatj(indmaxXj,2)], ...
    [Xhatj(indminXj,3); Xhatj(indmaxXj,3)],'LineWidth',lwd);
textscatter3(Xhatj(indmaxXj,1), Xhatj(indmaxXj,2), ...
    Xhatj(indmaxXj,3),"PC"+j);
end

function add3Daxes(j,color,lwd)
h = gca;
if j==1
    plot3(h.XLim, [0 0], [0, 0], color,'LineWidth',lwd);
elseif j==2
    plot3([0, 0], h.YLim, [0 0], color,'LineWidth',lwd);
elseif j==3
    plot3([0, 0], [0 0], h.ZLim, color,'LineWidth',lwd);
else
end
end