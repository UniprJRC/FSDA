function  [Ztable,Rtable,explained,explainedT,V,VT,loadings,loadingsT,communwithcum,communwithcumT,score,scoreT,orthDist,scoreDist] = ...
    computePCA(Y,bsb,rownames,varnames,standardize,NumComponents,dispresults,plots,Latitude,Longitude,ShapeFile)
% Compute all PCA quantities (this function is not intended to be called
% directly)
[n,v]=size(Y);
if isempty(bsb)
    bsb=true(n,1);
end
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
    % Delete figure if it already exists

    delete(findobj(0, 'type', 'figure','tag','pl_eigen'));
    figure('Name','Explained variance','Tag','pl_eigen')
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
    delete(findobj(0, 'type', 'figure','tag','pl_loadings'));
    figure('Name','Loadings','Tag','pl_loadings')

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
    delete(findobj(0, 'type', 'figure','tag','pl_OutlierMap'));
    figure('Name','OutlierMap','tag','pl_OutlierMap')
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

    %% Show the plot of latitude and longitude

    if ~isempty(Latitude) &&  ~isempty(Longitude) 
        delete(findobj(0, 'type', 'figure','tag','pl_latlong'));
        figure('Name','geobubble','tag','pl_latlong')

        % Divide the second PC into 4 classes
        cate=discretize(score(:,2),linspace(min(score(:,2)),max(score(:,2)),4),'categorical');

        % Add title and legend titles
        gb=geobubble(Latitude,Longitude,score(:,1),cate,'Basemap','topographic','MapLayout','maximized');
        gb.Title = "Geobubble of all units";
        gb.SizeLegendTitle = "First PC (size)";
        gb.ColorLegendTitle = "Second PC (color)";
    elseif ~isempty(Latitude) ||  ~isempty(Longitude) 
            error('FSDA:pcaFS:WrongInputOpt','Both Latitude and Longitude must be given.');
    else
    end

    if ~isempty(ShapeFile)
        LatLongAPP(Ztable,score,ShapeFile)
    end


    % % ShapeFile
    % if ~isempty(ShapeFile)
    %
    %
    %     delete(findobj(0, 'type', 'figure','tag','pl_geoplot'));
    %     figure('Name','geoplot','tag','pl_geoplot')
    %     % Specify the variable for which the colors have to be set
    %     colorIntensity = score(:,1);
    %     Xgeotable = readgeotable(ShapeFile);
    %
    %     % Define a colorbar from green to red
    %     numColors=50;
    %     customGradient = interp1([0;1],[0 1 0; 1 0 0],linspace(0,1,numColors));
    %
    %
    %     % Check if Mapping toolbox is installed
    %     instTOO=struct2table(ver);
    %     if any(string(instTOO{:,1})=="Mapping Toolbox")
    %
    %         geoaxes
    %         % Add to the geotable the Variable which contains the colors
    %         Xgeotable.ColorVar=colorIntensity;
    %         geoplot(Xgeotable,'ColorVariable','ColorVar');
    %         colormap(customGradient);
    %         h=colorbar;
    %         h.Label.String='Color based on first PC';
    %         geobasemap('topographic')
    %
    %         % Add labels in the centroids of the shapes
    %         % estract XY pairs coords of enclosing poligons for each row of input table
    %         Xtab= geotable2table(Xgeotable, ["X" "Y"]);
    %         % Obtain the averate latitute and longitude for each row
    %         info = shapeinfo(ShapeFile);
    %         proj = info.CoordinateReferenceSystem;
    %         arrlonglat=zeros(n,2);
    %         for i = 1:n
    %             [lat, lon] = projinv(proj, Xtab.X{i},Xtab.Y{i});  % Convert to lat/lon
    %             arrlonglat(i,1:2)=[mean(lat,'omitmissing') mean(lon,'omitmissing')];
    %         end
    %
    %         % Show the labels
    %         text(arrlonglat(:,1), arrlonglat(:,2),rownames)
    %
    %     else % mapping toolbox is not installed
    %
    %         fc=customGradient(round(rescale(colorIntensity,1,numColors)),:);
    %         polygon=mapshow(Xgeotable, "FaceAlpha",0.9);
    %         for i=1:n
    %             set(polygon.Children(i), 'FaceColor', fc(i,:),"FaceAlpha",0.9);
    %         end
    %         colormap(customGradient);
    %         rang=[min(colorIntensity) max(colorIntensity)];
    %         h=colorbar;
    %         h.Label.String='Color based on first PC';
    %         soglie=linspace(min(rang),max(rang),6);
    %         valori=(soglie-min(soglie))/(max(soglie)-min(soglie));
    %         h.Ticks=valori;
    %         h.TickLabels=string(soglie);
    %
    %         % Add labels in the centroids of the shapes
    %         arrlonglat=zeros(n,2);
    %         % estract XY pairs coords of enclosing poligons for each row of input table
    %         Xtab= geotable2table(Xgeotable, ["X" "Y"]);
    %
    %         for i=1:n
    %             % prima creiamo degli oggetti poligonali
    %             polObj= polyshape(Xtab.X{i,1}, ...
    %                 Xtab.Y{i,1});
    %             [centX,centY]=centroid(polObj);
    %             % popoliamo i vettori
    %             arrlonglat(i,1:2)=[centX centY];
    %         end
    %         text(arrlonglat(:,1), arrlonglat(:,2),rownames)
    %
    %
    %     end
    % end
end