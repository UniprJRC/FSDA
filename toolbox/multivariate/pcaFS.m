function out=pcaFS(Y,varargin)
%pcaFS performs Principal Component Analysis (PCA) on raw data.
%
%<a href="matlab: docsearchFS('pcaFS')">Link to the help function</a>
%
%   The main differences with respect to MATLAB function pca are:
%   1) accepts an input X also as table;
%   2) produces in table format the percentage of the variance explained
%      single and cumulative of the various components and the associated
%      scree plot, in order to decide about the number of components to
%      retain.
%   3) returns the loadings in table format and shows them graphically.
%   4) provides guidelines about the automatic choice of the number of
%       components;
%   5) returns the communalities for each variable with respect to the
%       first k principal components in table format;
%   6) retuns the orthogonal distance ($OD_i$) of each observation to the
%      PCA subspace. For example, if the subspace is defined by the first two
%      principal components, $OD_i$ is computed as:
%      \[
%        OD_i=|| z_i- V_{(2)} V_{(2)}' z_i ||
%      \]
%      where z_i is the i-th row of the original centered data matrix $Z$ of
%      dimension $n \times v$ and $V_{(2)}=(v_1 v_2)$ is the matrix of size
%      $p \times 2$ containing the first two eigenvectors of $Z'Z/(n-1)$. The
%      observations with large $OD_i$ are not well represented in the space of
%      the principal components.
%   7)  returns the score distance $SD_i$ of each observation. For example,
%      if the subspace is defined by the first two principal components,
%      $SD_i$ is computed as:
%      \[
%        SD_i=\sqrt{(z_i'v_1)^2/l_1+ (z_i'v_2)^2/l_2 }
%      \]
%      where $l_1$ and $l_2$ are the first two eigenvalues of $Z'Z/(n-1)$.
%   8) calls app biplotFS which enables to obtain an interactive biplot in
%      which points, rowslabels or arrows can be shown or hidden. This app
%      also gives the possibility of controlling the length of the arrows
%      and the position of the row points through two interactive slider
%      bars. In the app it is also possible to color row points depending
%      on the orthogonal distance ($OD_i$) of each observation to the PCA
%      subspace. If optional input argument bsb or bdp is specified, it is
%      possible to have in the app two tabs which enable the user to select
%      the breakdown point of the analysis or the subset size to use in the
%      svd. The units which are declared as outliers or the units outside
%      the subset are shown in the biplot with filled circles.
%
%
%  Required input arguments:
%
% Y :           Input data. 2D array or table.
%               n x p data matrix; n observations and p variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
%  Optional input arguments:
%
%
%     biplot     : launch app biplotFS. Scalar. If biplot is 1
%                   (default) app biplotFS is automatically launched. With
%                   this app it is possible to show in a dynamic way the
%                   rows points (PC coordinates), the arrows, the row
%                   labels and control with a scrolling bar the length of
%                   the arrows and the spread of row points.
%                   Example - 'biplot',0
%                   Data Types - double
%
%  dispresults   : show the results in the command window. If dispresults
%                   is true, the percentage of variance explained together
%                   with the loadings, the criteria for deciding the number
%                   of components to retain and the 5 units with the
%                   largest score and orthogonal distance (combined) are
%                   shown in the command window.
%                   Example - 'dispresults',false
%                    Data Types - char
%
%    Latitude    : Latitude coordinates in degrees of the rows. nx1 vector containint the
%                   numerical values of latitudes for each row.
%                   If vectors Latitude and Longitude are present once
%                   brushing of the units in biplotAPP is done, the
%                   geobubble of the brushed units automatically
%                   appears. The size of the bubble is proportional to the
%                   value of the first principal component while the color
%                   of the bubbles depends on the value of the second
%                   principal component.
%                   Example - 'Latitude',[54 43.2 ... 47]
%                    Data Types - vector with elements in the range [–90, 90]
%
%    Longitude    : Longitude coordinates in degrees of the rows. nx1 vector containint the
%                   numerical values of latitudes for each row.
%                   If vectors Latitude and Longitude are present once
%                   brushing of the units in biplotAPP is done, the
%                   geobubble of the brushed units automatically
%                   appears. The size of the bubble is proportional to the
%                   value of the first principal component while the color
%                   of the bubbles depends on the value of the second
%                   principal component.
%                   Example - 'Longitude',[54 43.2 ... 47]
%                    Data Types - vector with elements in the range [–90, 90]
%
%  NumComponents : the number of components desired. Specified as a
%                  scalar integer $k$ satisfying $0 < k \leq p$. When
%                  specified, pcaFS returns the first $k$ columns of
%                  out.coeff and out.score. If NumComponents is not
%                  specified pcaFS returns the minimum number of components
%                  which cumulatively enable to explain a percentage of
%                  variance which is equal at least to $0.95^p$. If this
%                  threshold is exceeded already by the first PC, pcaFS
%                  still returns the first two PCs.
%                   Example - 'NumComponents',2
%                    Data Types - char
%
%      plots     : plots on the screen. Scalar. If plots is 1 (default) it is
%                   possible to show on the screen the scree plot of the
%                   variance explained, the plot of the loadings for the
%                   first two PCs.
%                   Example - 'plots',0
%                   Data Types - double
%
%       robust   : robust principal components. boolean or struct. If
%               robust is a scalar boolean equal to true (default), FS is
%               applied and using a Bonferronized confidence level units
%               are declared as outliers. If robust is a struct it is
%               possible to choose between FS or MCD and the structure may
%               contain the following fields:
%           robust.bsb = units forming subset on which to
%                  compute the svd. bsb can be either a
%                  numeric vector of length m (m<=n) containing the list of the
%                  units (e.g. 1:50) or a logical vector of length n
%                  containing the true for the units which have to be used
%                  in the calculation of svd. For example bsb=true(n,1),
%                  bsb(13)=false; excludes from the svd unit number 13.
%                  The other units are projected into this robust space.
%                  Note that is robust.bsb is supplied all the other
%                  options are ignored.
%           robust.conflev = scalar in the interval (0 1) which contains
%                   the confidence level to declare the units as outliers.
%                   if robust.conflev is not specified or if it is missing
%                   a 0.99 Bonferronized confidence level is used.
%           robust.class = class can be either 'FS' (default) or 'MCD' or 'MM'
%                   in order to choose between Forward search, or minimum
%                   covariance determinant or MM estimators.
%           robust.bdp= break down point to use in the 'FS' or 'MCD' estimator or in the
%                   preliminary S estimator which is used to compute the MM.
%                   The default is bdp=0.4; It measures the fraction of outliers the algorithm should
%               resist. In this case any value greater than 0 but smaller
%               or equal than 0.5 will do fine.
%           robust.eff = nominal efficiency to use for the MM estimator
%                       (note that  this field is used just if robust.class
%                       is 'MM')
%                 Example - 'bsb',[2 10:90 93]
%                 Data Types - double or logical
%
%  ShapeFile  : name of ShapeFile or geospatial table containg shapes.
%               Character or string or geotable.
%               Name of the ShapeFile containing the containing the
%               geometric details of the rows. The ShapeFile, which is
%               loaded using function readgeotable, must have n rows and
%               the n rows must have the same order of the n rows of Y. The
%               default value of shapefile is empty that is we assume that
%               no shapefile is given. If ShapeFile is given an additional
%               GUI containing the areas colored using
%               the first PC is shown.
%                   Example - 'ShapeFile','shapefileName'
%                    Data Types - char or string or geotable.
%               Remark: note that this option can be used just is the
%               "Mapping toolbox" is installed.
%
%    standardize : standardize data. boolean. Boolean which specifies
%               whether to standardize the variables, that is we operate on
%               the correlation matrix (default) or simply remove column
%               means (in this last case we operate on the covariance
%               matrix).
%                   Example - 'standardize',false
%                   Data Types - boolean
%
% Output:
%
%
%         out:   structure which contains the following fields
%
%out.Rtable = p-by-p correlation matrix in table format.
%
%     out.class = character vector containing pcaFS.
%
% out.explained = p \times 3 matrix containing respectively
%                1st col = eigenvalues;
%                2nd col = Explained Variance (in percentage)
%                3rd col = Cumulative Explained Variance (in percentage)
%
%out.explainedT = the same as out.explained but in table format.
%
%out.coeff=  p-by-NumComponents matrix containing the ordered eigenvectors
%           of the correlation (covariance matrix) in table format.
%            First column is referred to first eigenvector ...
%            Note that out.coeff'*out.coeff= I_NumComponents.
%
%out.coeffT = the same as out.coeff but in table format.
%
%out.loadings=p-by-NumComponents matrix containing the correlation
%             coefficients between the original variables and the first
%             NumComponents principal components.
%
%out.loadingsT = the same as out.loadings but in table format.
%
% out.score= the principal component scores. The rows of out.score
%            correspond to observations, columns to components. The
%            covariance matrix of out.score is $\Lambda$ (the diagonal
%            matrix containing the eigenvalues of the correlation
%            (covariance matrix).
%
% out.scoreT = the same as outscore but in table format.
%
% out.communalities = matrix with p-by-2*NumComponents-1 columns.
%               The first NumComponents columns contain the communalities
%               (variance extracted) by the the first NumComponents
%               principal components. Column NumComponents+1 contains the
%               communalities extracted by the first two principal
%               components. Column NumComponents+2 contains the
%               communalities extracted by the first three principal
%               components...
%
%  out.communalitiesT= the same as out.communalities but in table format.
%
%  out.orthDist = orthogonal distance from PCA subspace.
%                 Column vector of length n containing the orthogonal
%                 distance of each observation from the PCA subspace.
%
%  out.scoreDist = score distance from centroid.
%                 Column vector of length n containing the score
%                 distance of each observation from the PCA subspace.
%                 The analysis of out.orthDist and out.scoreDist reveals
%                 the good leverage points, the orthogonal outliers and the
%                 bad leverage points.
%                 Good leverage points: points which lie close to the PCA
%                 space but far from the regular observations. Good
%                 leverage points have a large score distance and low
%                 orthogonal distance. Orthogonal outliers are points which
%                 have a large orthogonal distance to the PCA subspace but
%                 cannot be seen when we look only at their projection on
%                 the PCA subspace. Bad leverage points are points which
%                 have a large orthogonal distance and whose projection on
%                 the PCA subspace is remote from the typical projections.
%                 These points have a large score distance and a large
%                 orthogonal distance.
%
%
% See also: pca, biplotFS
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('pcaFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    % Use of pcaFS with citiesItaly dataset.
    load citiesItaly;
    % Use all default options
    out=pcaFS(citiesItaly);
%}

%{
    %% Use of pcaFS on the ingredients dataset.
    load hald
    % Operate on the covariance matrix.
    out=pcaFS(ingredients,'standardize',false,'biplot',0);
%}

%{
    %% Use of pcaFS with option robust true.
    load citiesItaly;
    % Use all default options
    out=pcaFS(citiesItaly,'robust',true);
%}

%{
    %% Use of pcaFS with options Latitude and Longitude.
    load citiesItaly2024.mat
    X=citiesItaly2024;
    % Retrieve Latitude and Longitude of each province
    LatLong=X.Properties.UserData{2};
    Latitude=LatLong(:,1);
    Longitude=LatLong(:,2);
    out=pcaFS(X,'Latitude',Latitude,'Longitude',Longitude);
%}

%{
    %% Use of pcaFS with option ShapeFile.
    % Note that this option requires the Mapping Toolbox to be installed.
    a=struct2table(ver);
    MappingInstalled=any(string(a{:,1})=="Mapping Toolbox");
    if MappingInstalled ==true
        load citiesItaly2024.mat
        X=citiesItaly2024;
        ShapeFile=X.Properties.UserData{1};
        out=pcaFS(X,"ShapeFile",ShapeFile,'biplot',0);
    else
           disp('This option requires that the "mapping toolbox" is installed') 
    end
%}


%% Beginning of code
[n,v]=size(Y);
plots=1;
standardize=true;
biplot=1;
dispresults=true;
NumComponents=[];
robust=false;
Latitude=[];
Longitude=[];
ShapeFile='';

if nargin>1
    options=struct('plots',plots, ...
        'standardize',standardize,'biplot', biplot,...
        'dispresults',dispresults,'NumComponents',NumComponents,...
        'robust',robust,'Latitude',Latitude,'Longitude',Longitude, ...
        'ShapeFile',ShapeFile);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)


        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:pcaFS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end

        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:pcaFS:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end


    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    plots=options.plots;
    standardize=options.standardize;
    biplot=options.biplot;
    dispresults=options.dispresults;
    NumComponents=options.NumComponents;
    robust=options.robust;
    Longitude=options.Longitude;
    Latitude=options.Latitude;
    ShapeFile=options.ShapeFile;
end

if istable(Y)
    varnames=Y.Properties.VariableNames;
    rownames=Y.Properties.RowNames;
    if isempty(rownames)
        rownames=cellstr(num2str((1:n)','%d'));
    end
    Y=table2array(Y);
else
    varnames=cellstr(num2str((1:v)','Y%d'));
    rownames=cellstr(num2str((1:n)','%d'));
end

bsb=[];

if isstruct(robust) || islogical(robust)

    if isstruct(robust)
        if isfield(robust,'bsb') % User has chosen a prespecified subset of units to use

            bsb=false(n,1);
            bsb(robust.bsb)=true;

        elseif  isfield(robust,'class')
            bsb=true(n,1);

            class=robust.class;
            if strcmp(class,'MCD')
                % User has chosen MCD
                if isfield(robust,'conflev')
                    conflev=robust.conflev;
                else
                    conflev=1-0.01/n;
                end
                if isfield(robust,'bdp')
                    bdp=robust.bdp;
                else
                    bdp=0.5;
                end
                outMCD=mcd(Y,'bdp',bdp,'conflev',conflev,'plots',0);
                bsb(outMCD.outliers)=false;
                if sum(bsb)==n
                    nooutlier();
                end


            elseif strcmp(class,'MM') % user has chosen MM estimators
                if isfield(robust,'eff') % User has specified eff of MM estimators
                    eff=robust.eff;
                else
                    eff=0.95;
                end

                if isfield(robust,'conflev')
                    conflev=robust.conflev;
                else
                    conflev=1-0.01/n;
                end

                outMM=MMmult(Y,'eff',eff,'conflev',conflev);
                bsb(outMM.outliers)=false;

                if sum(bsb)==n
                    nooutlier();
                end

            elseif strcmp(class,'FS') % user has chosen FS
                if isfield(robust,'conflev')
                    conflev=robust.conflev;
                else
                    conflev=0.99;
                end
                outFS=FSM(Y,'bonflev',conflev,'plots',0,'msg',0);
                bsb(outFS.outliers)=false;
                if sum(bsb)==n
                    nooutlier();
                end

            else
                error('FSDA:pcaFS:wrongCLASS','Admissible classes are "FS", "MCD" or "MM"')
            end
        else
            error('FSDA:pcaFS:missingCLASS','Specify a class among "FS", "MCD" or "MM"')
        end

    else % in this case robust is a scalar boolean

        if robust == true
            bsb=true(n,1);
            outFS=FSM(Y,'bonflev',0.99,'msg',0,'plots',0);
            bsb(outFS.outliers)=false;
            if sum(bsb)==n
                nooutlier();
            end
        end
    end
end

%{
    Ybsb=Y(bsb,:);
    nbsb=size(Ybsb,1);

    center=mean(Ybsb);
    if standardize==true
        dispersion=std(Ybsb);
        % Create matrix of standardized data
        Z=(Y-center)./dispersion;
    else
        % Create matrix of deviations from the means
        Z=Y-center;
    end

    % [~,S,loadings]=svd(Z./sqrt(n-1),0);
    % Z=(Y-mean(Y))*loadings;

    Ztable=array2table(Z,'RowNames',rownames,'VariableNames',varnames);

    % Correlation (Covariance) matrix in table format
    Zbsb=Z(bsb,:);
    R=corr(Zbsb);
    Rtable=array2table(R,'VariableNames',varnames,'RowNames',varnames);

    sigmas=sqrt(diag(R));

    % svd on matrix Z.
    [~,Gamma,V]=svd(Zbsb,'econ');
    Gamma=Gamma/sqrt(nbsb-1);

    % \Gamma*\Gamma = matrice degli autovalori della matrice di correlazione
    La=Gamma.^2;
    la=diag(La);

    %% Explained variance
    sumla=sum(la);
    explained=[la 100*(la)/sumla 100*cumsum(la)/sumla];
    namerows=cellstr([repmat('PC',v,1) num2str((1:v)')]);
    namecols={'Eigenvalues' 'Explained_Variance' 'Explained_Variance_cum'};
    explainedT=array2table(explained,'RowNames',namerows,'VariableNames',namecols);
    if isempty(NumComponents)
        NumComponents=find(explained(:,3)>100*0.95^v,1);
        if NumComponents==1 && v>1
            disp('The first PC already explains more than 0.95^v variability')
            disp('In what follows we still extract the first 2 PCs')
            NumComponents=2;
        end
    end

    % labels of the PCs
    pcnames=cellstr(num2str((1:NumComponents)','PC%d'));

    V=V(:,1:NumComponents);
    La=La(1:NumComponents,1:NumComponents);
    VT=array2table(V,'RowNames',varnames','VariableNames',pcnames);


    %% Loadings
    loadings=V*sqrt(La)./sigmas;
    loadingsT=array2table(loadings,'RowNames',varnames','VariableNames',pcnames);


    %% Principal component scores
    score=Z*V;
    scoreT=array2table(score,'RowNames',rownames,'VariableNames',pcnames);


    %% Communalities
    commun=loadings.^2;
    labelscum=cellstr([repmat([pcnames{1} '-'],NumComponents-1,1) char(pcnames{2:end})]);
    communcum=cumsum(loadings.^2,2);
    communwithcum=[commun communcum(:,2:end)];
    if isempty(labelscum{1})
        varNames=pcnames;
    else
        varNames=[pcnames; labelscum];
    end

    if verLessThanFS('9.7')
        varNames=matlab.lang.makeValidName(varNames);
    end
    communwithcumT=array2table(communwithcum,'RowNames',varnames,...
        'VariableNames',varNames);

    %% Orthogonal distance to PCA subspace based on k PC
    Res=Z-score*V';
    orthDist=sqrt(sum(Res.^2,2));


    %% Score distance in PCA subspace of dimension k
    larow=diag(La)';
    scoreDist=sqrt(sum(score.^2./larow,2));

    % Find the 5 units with the largest value of the combination between
    % orthogonal and score distance
    DD=[orthDist,scoreDist];
    mDD=zeros(1,2);
    distM=mahalFS(DD,mDD,cov(DD));
    [~,indsor]=sort(distM,1,"descend");
    selu=indsor(1:5);
%}



%{

    if dispresults == true
        format bank
        if standardize == true
            disp('Initial correlation matrix')
        else
            disp('Initial covariance matrix')
        end
        disp(Rtable)

        disp('Explained variance by PCs')
        disp(explainedT)

        disp('Loadings = correlations between variables and PCs')
        disp(loadingsT)

        disp('Communalities')
        disp(communwithcumT)
        format short

        disp('Units with the 5 largest values of (combined) score and orthogonal distance')
        disp(selu')
    end

    if plots==1

        %% Explained variance through Pareto plot
        figure('Name','Explained variance')
        [h,axesPareto]=pareto(explained(:,1),namerows);
        % h(1) refers to the bars h(2) to the line
        h(1).FaceColor='g';
        linelabels = string(round(100*h(2).YData/sumla,2));
        text(axesPareto(2),h(2).XData,h(2).YData,linelabels,...
            'Interpreter','none');
        xlabel('Principal components')
        ylabel('Explained variance (%)')

        %% Plot loadings
        xlabels=categorical(varnames,varnames);
        figure('Name','Loadings')

        for i=1:NumComponents
            subplot(NumComponents,1,i)
            b=bar(xlabels, loadings(:,i),'g');
            title(['Correlations with PC' num2str(i)])
            xtips=b(1).XData;
            ytips=b(1).YData;
            % The alternative instructions below only work from MATLAB
            % 2019b
            %   xtips = b.XEndPoints;
            %   ytips = b.YEndPoints;
            barlabels = string(round(loadings(:,i),2));
            text(xtips,ytips,barlabels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
            title(['Correlations  with PC' num2str(i)])
        end

        %% Plot of orthogonal distance (Y) versus score distance (X)
        figure('Name','OutlierMap')
        group1=repelem("Normal units",n,1);
        group1(bsb==0)="Outliers";
        scatterboxplot(scoreDist,orthDist,'group',group1);
        xlabel('Score distance')
        ylabel('Orth. dist. from PCA subspace')
        text(1.01,0,['Good' newline 'leverage' newline 'points'],'Units','normalized')
        text(-0.05,0,['Normal' newline 'units'],'Units','normalized','HorizontalAlignment','right')
        text(0.05,-0.05,['Normal' newline 'units'],'Units','normalized','HorizontalAlignment','left')
        text(-0.05,1.05,'Orthogonal outliers','Units','normalized','HorizontalAlignment','left')
        text(0.95,1.05,'Bad leverage points','Units','normalized','HorizontalAlignment','right')
        text(1.01,0.95,['Bad' newline 'leverage' newline 'points'],'Units','normalized','HorizontalAlignment','left')
        text(scoreDist(selu),orthDist(selu),rownames(selu),'HorizontalAlignment','left','VerticalAlignment','bottom');
        % Good leverage points: points which lie close to the PCA space but far
        % from the regular observations.
        % Orthogonal outliers points: points which have a large orthogonal distance
        % to the PCA space but cannot be seen when we look only at their
        % projection on the PCA subspace.
        % Bad leverage points: points which have a large orthogonal distance
        % and whose projection on the PCA subspace is remote from the typical
        % projections.

    end
%}


% The part above which has been commented is done inside the function below
% In this way the function which does PCA computations can be called
% without pcaFS
[Ztable,Rtable,explained,explainedT,V,VT,loadings,loadingsT,communwithcum,communwithcumT, ...
    score,scoreT,orthDist,scoreDist]=aux.computePCA(Y,bsb,rownames,varnames,standardize,NumComponents,dispresults,plots,Latitude,Longitude,ShapeFile);

out=struct;
out.Rtable=Rtable;
out.explained=explained;
out.explainedT=explainedT;
out.coeff=V;
out.coeffT=VT;
out.loadings=loadings;
out.loadingsT=loadingsT;
out.communalities=communwithcum;
out.communalitiesT=communwithcumT;
out.score=score;
out.scoreT=scoreT;
out.orthDist=orthDist;
out.scoreDist=scoreDist;
out.class='pcaFS';


if biplot==1
    biplotAPP(Ztable,'standardize',standardize,'bsb',bsb,'Latitude',Latitude,'Longitude',Longitude,'ShapeFile',ShapeFile)
end

end

% Message to show when robust option has been used but no outlier has been
% found
function nooutlier()
disp('--------------------------------------------')
disp('No outlier has been found by robust analysis')
disp('                              ')
disp('To explore a different confidence level or a different robust method')
disp('Click on the Check box "Robust analysis" in the FSDA dynamic bixplot')
disp('--------------------------------------------')
end

%FScategory:MULT-Multivariate

