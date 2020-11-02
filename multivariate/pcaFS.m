function out=pcaFS(Y,varargin)
%pcaFS performs Principal Component Analysis (PCA) on raw data.
%
%<a href="matlab: docsearchFS('pcaFS')">Link to the help function</a>
%
%   The main differences with respect to MATLAB function pca are:
%   1) accepts an input X also as table;
%   2) produces in table format the percentage of the variance explained
%      single and cumulative of the various components and the associated
%      scree plot in order to decide about the number of components to
%      retain.
%   3) returns the loadings in table format and shows them graphically.
%   4) provides guidelines about the automatic choice of the number of
%       components;
%   5) returns the communalities for each variable with respect to the
%       first k principal components in table format; 
%   5) calls app biplotFS which enables to obtain an interactive biplot in
%      which points, rowslabels or arrows can be shown or hidden. This app
%      also gives the possibility of controlling the length of the arrows
%      and the position of the row points through two interactive slider
%      bars.
%
%  Required input arguments:
%
% Y :           Input data. 2D array or table.
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
%  Optional input arguments:
%
%    standardize : standardize data. boolean. Boolean which specifies
%               whether to standardize the variables, that is we operate on
%               the correlation matrix (default) or simply remove column
%               means (in this last case we operate on the covariance
%               matrix).
%                   Example - 'standardize',false
%                   Data Types - boolean
%
%      plots     : plots on the screen. boolean. If plots is true (default) it is
%                   possible to show on the screen the scree plot of the
%                   variance explained, the plot of the loadings for the
%                   first two PCs.
%                   Example - 'plots',false
%                   Data Types - boolean
%
%     biplot     : launch app biplotFS. boolean. If biplot is true
%                   (default) app biplotFS is automatically launched. With
%                   this app it is possible to show in a dynamic way the
%                   rows points (PC coordinates), the arrows, the row
%                   labels and control with a scrolling bar the length of
%                   the arrows and the spread of row points.
%                   Example - 'biplot',false
%                   Data Types - boolean
%
%
%  dispresults   : show the results in the command window. If dispresults
%                   is true, the percentage of variance explained together
%                   with the loadings and the criteria for deciding the
%                   number of components to retain is shown in the command
%                   window.
%                   Example - 'dispresults',false
%                    Data Types - char
%
%  NumComponents : the number of components desired. Specified as a
%                  scalar integer $k$ satisfying $0 < k \leq v$. When
%                  specified, pcaFS returns the first $k$ columns of
%                  out.coeff and out.score. If NumComponents is not
%                  specified pcaFS returns the minimum number of components
%                  which cumulatively enable to explain a percentage of
%                  variance which is equal at least to $0.95^v$. If this
%                  threshold is exceeded already by the first PC, pcaFS
%                  still returns the first two PCs.
%                   Example - 'NumComponents',2
%                    Data Types - char
%
%
% Output:
%
%
%         out:   structure which contains the following fields
%
%out.Rtable = v-by-v correlation matrix in table format.
%
% out.explained = v \times 3 matrix containing respectively
%                1st col = eigenvalues;
%                2nd col = Explained Variance (in percentage)
%                3rd col = Cumulative Explained Variance (in percentage)
%
%out.explainedT = the same as out.explained but in table format.
%
%out.coeff=  v-by-NumComponents matrix containing the ordered eigenvectors
%           of the correlation (covariance matrix) in table format.
%            First column is referred to first eigenvector ...
%            Note that out.coeff'*out.coeff= I_NumComponents.
%
%out.coeffT = the same as out.coeff but in table format.
%
%out.loadings=v-by-NumComponents matrix containing the correlation
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
% out.communalities = matrix with v-by-2*NumComponents-1 columns.
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
% See also: pca, biplotFS
%
% References:
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('pcaFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    % Use of pcaFS with creditrating dataset.
    creditrating = readtable('CreditRating_Historical.dat','ReadRowNames',true);
    % Use all default options
    out=pcaFS(creditrating(1:100,1:6))
%}

%{
    %% use of pcaFS on the ingredients dataset.
    load hald
    % Operate on the covariance matrix.
    out=pcaFS(ingredients,'standardize',false,'biplot',false);
%}


%% Beginning of code
[n,v]=size(Y);
plots=true;
standardize=true;
biplot=true;
dispresults=true;
NumComponents=[];

if nargin>1
    options=struct('plots',plots, ...
        'standardize',standardize,'biplot', biplot,...
        'dispresults',dispresults,'NumComponents',NumComponents);
    
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
end

if istable(Y)
    varnames=Y.Properties.VariableNames;
    rownames=Y.Properties.RowNames;
    Y=table2array(Y);
else
    varnames=cellstr(num2str((1:v)','Y%d'));
    rownames=cellstr(num2str((1:n)','%d'));
end

if standardize==true
    % Create matrix of standardized data
    Z=zscore(Y);
else
    % Create matrix of deviations from the means
    Z=Y-mean(Y);
end
Ztable=array2table(Z,'RowNames',rownames,'VariableNames',varnames);


% Correlation (Covariance) matrix in table format
R=cov(Z);
Rtable=array2table(R,'VariableNames',varnames,'RowNames',varnames);

sigmas=sqrt(diag(R));

% svd on matrix Z.
[~,Gamma,V]=svd(Z,'econ');
Gamma=Gamma/sqrt(n-1);

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
    if NumComponents==1
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
communwithcumT=array2table(communwithcum,'RowNames',varnames,...
    'VariableNames',[pcnames; labelscum]);

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
end

if plots==true
    
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
        
        xtips = b.XEndPoints;
        ytips = b.YEndPoints;
        barlabels = string(round(loadings(:,i),2));
        text(xtips,ytips,barlabels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
        title(['Correlations  with PC' num2str(i)])
    end
end
if biplot==true
    biplotFS(Ztable,'standardize',standardize)
end

end

%FScategory:MULT-Multivariate

