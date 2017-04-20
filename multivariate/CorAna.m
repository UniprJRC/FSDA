function out=CorAna(N, varargin)
%CorAna performs correpondence analysis
%
%
%<a href="matlab: docsearchFS('CorAna')">Link to the help function</a>
%
%  Required input arguments:
%
%
%       N    :    Contingency table (default) or n-by-2 input datasets. Matrix or Table.
%                 Matrix or table which contains the input contingency table (say of
%                 size I-by-J) or the original data matrix.
%                 In this last case N=crosstab(N(:,1),N(:,2)).
%                 As default procedure assumes that the input is a
%                 contingency table.
%
%
%  Optional input arguments:
%
%       k    :    Number of dimensions to retain. Scalar.
%                 Scalar which contains the number of dimensions to retain.
%               Example - 'k',2
%               Data Types - double
%       Lr   :  Vector of row labels. Cell.
%               Cell containing the labels of the rows of the input
%               contingency matrix N. This option is unnecessary if N is a table.
%               because in this case  Lr=N.Properties.RowNames;
%               Example - 'Lr',{'a' 'b' 'c'}
%               Data Types - cell array of strings
%       Lc   :  Vector of column labels. Cell.
%               Cell containing the labels of the columns of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table because in this case Lc=N.Properties.VariableNames;
%               Example - 'Lc',{'c1' c2' 'c3' 'c4'}
%               Data Types - cell array of strings
%       Sup :  Structure containing indexes of supplementary rows or columns. Structure.
%              Structure with the followin fields.
%               Sup.r = indexes or cell array of strings or table containing
%                       supplementary rows. If indexes or cell array of
%                       strings are supplied, we assume that supplementary
%                       rows belong to contingency table N.
%                       For example, if Sup.r=[2 5] we use rows 2 and 5 of
%                       the input contingency table as supplementary rows.
%                       Of course the length of Sup.r must be smaller than
%                       the number of rows of the contigencey matrix
%                       divided by 2. If, on the other hand, Sup.r is a
%                       table of size size(N) supplementary rows do not belong
%                       to N.
%               Sup.c = indexes or cell array of string containing
%                       supplementary columns, or table. If indexes or cell array of
%                       strings are supplied, we assume that supplementary
%                       columns belong to contingency table N.
%                       For example, if Sup.c={'Smokers' 'NonSmokers'} we
%                       use columns of the contingencey table labelled
%                       'Smokers' and 'NonSmokers' of the input contingency
%                       table N as supplementary columns.
%                       Of course the length of Sup.c must be smaller than
%                       the number of columns of the contigencey matrix
%                       divided by 2. If, on the other hand, Sup.c is a
%                       table of size size(N) supplementary columns  do not belong
%                       to N.
%      datamatrix   :  data matrix or contingency table. Scalar or empty value. If
%                       datamatrix is 1 the first input argument N is
%                       forced to be interpreted as a data matrix, else
%                       if the input argument different from 1, N
%                       is treated as a contingency table. The default value
%                       of datamatrix is 0, that is the procedure
%                       automatically considers N as a contingency table
%               Example - 'datamatrix',1
%               Data Types - single | double
%       plots : Plot on the screen. Scalar or structure.
%               If plots = 1, a plot the plot which shows  the Principal
%               coordinates of rows and columns is shown on the screen. If
%               plots is a structure it may contain the following fields:
%               plots.alpha = type of plot. scalar in the interval [0 1] or
%               a string identifyind the type of coordinates to use in the
%               plot.
%               If plots.alpha='rowprincipal' The row points are in
%                   principal coordinates and the column coordinates are
%                   standard coordinates. Distances between row points are
%                   (approximated) chi-squared distances
%                   (row-metric-preserving). The position of the row points
%                   are at the weighted average of the column points.
%                   'rowprincipal' can also be specified setting
%                   plots.alpha=1.
%               If plots.alpha='colprincipal' the column
%                   coordinates are referred to as principal coordinates
%                   and the row coordinates as standard coordinates.
%                   Distances between column points are (approximated)
%                   chi-squared distances (row-metric-preserving). The
%                   position of the column points are at the weighted
%                   average of the row points. 'rowprincipal' can also be
%                   specified setting plots.alpha=0.
%               If plots.alpha='symbiplot'  The row and column coordinates
%                   are scaled similarly. The sum of weighted squared
%                   coordinates for each dimension is equal to the
%                   corresponding singular values.  These coordinates are often
%                   called symmetrical coordinates. This representation is
%                   particularly useful if one is primarily interested in
%                   the relationships between categories of row and column
%                   variables rather than in the distances among rows or
%                   among columns. 'symbiplot' can also be specified
%                   setting plots.alpha=0.5;
%               If plots.alpha='bothprincipal' both the rows and columns
%                   are depicted in principal coordinates. Such a plot is often
%                   referred to as a symmetrical plot or French symemtrical model. Note that such a
%                   symmetrical plot does not provide a feasible solution
%                   in the sense that it does not approximately approximate
%                   matrix $D_r^{-0.5}(P-rc')D_c^{-0.5}$.
%               If plots.alpha='rowgab'  rows are in principal coordinates
%                   and columns are in standard coordinates multiplied by
%                   the mass. This biplot has been suggested by Gabriel and
%                   Odoroff (1990).
%               If plots.alpha='colgab'  columns are in principal coordinates
%                   and rows are in standard coordinates multiplied by the
%                   mass. This biplot has been suggested by Gabriel and
%                   Odoroff (1990).
%               If plots.alpha='rowgreen'  rows are in principal
%                   coordinates and columns are in standard coordinates
%                   multiplied by square root of the mass.
%                  This biplot has been suggested by Greenacre and
%                  incorporates the contribution of points. In this
%                  display, points that contribute very little to the
%                  solution, are close to the center of the biplot and are
%                  relatively unimportant to the interpretation. This
%                  biplot is often referred as contribution biplot because
%                  it shows visually the most contributing points
%                  (Greenacre 2006b).
%               If plots.alpha='colgreen' columns in principal coordinates
%                   and rows in standard coordinates multiplied by the square
%                   root of the mass.
%                  This biplot has been suggested by Greenacre and
%                  incorporates the contribution of points. In this
%                  display, points that contribute very little to the
%                  solution, are close to the center of the biplot and are
%                  relatively unimportant to the interpretation. This
%                  biplot is often referred as contribution biplot because
%                  it shows visually the most contributing points
%                  (Greenacre 2006b).
%               If plots.alpha=scalar in the interval [0 1], row
%                   coordinates are given by $D_r^{-1/2} U \Gamma^\alpha$
%                   and column coordinates are given by $D_c^{-1/2} V
%                   \Gamma^{1-\alpha}$. Note that for any choice of $alpha$
%                   the matrix product $ D_r^{-1/2} U \Gamma^\alpha (D_c^{-1/2} V
%                   \Gamma^{1-\alpha})^T$ optimally approximates matrix
%                   $D_r^{-0.5}(P-rc')D_c^{-0.5}$, in the sense that the
%                   sum of squared differences between $D_r^{1/2}
%                   D_r^{-1/2} U \Gamma^\alpha (D_c^{-1/2} V
%                   \Gamma^{1-\alpha})^T D_c^{1/2}$ and
%                   $D_r^{-0.5}(P-rc')D_c^{-0.5}$ is as small as possible.
%              plots.dim = vector with two elements which specifies which
%               dimensions to show in the factor map. The default is to
%               show to first two dimensione therefore plots.dim=[1 2]
%                 Example - 'plots',1
%                 Data Types - scalar double or struct
%
%  Output:
%
%  out :     A structure containing the following fields
%
%
% 		out.N         =   Matrix containing Contigency table
% 		out.SupRow    =   Contigency table related to supplementary rows
% 		out.SupCol    =   Contigency table related to supplementary columns
% 		out.Lr        =   cell containing row labels.
% 		out.Lc        =   cell containing column labels
%       out.labels.sr =   Row supplementary labels
%       out.labels.sc =   Column supplementary labels
% 		out.I         =   Number of rows
% 		out.J         =   Number of columns
% 		out.n         =   Grand total
% 		out.P         =   Correspondence matrix (proportions)
% 		out.r         =   Row masses (vector)
% 		out.Dr        =   Row masses (diagonal matrix)
% 		out.c         =   Column masses (vector); c is the centroid of row profiles
% 		out.Dc        =   Column masses (diagonal matrix)
% 		out.RP        =   Rows profiles
% 		out.CP        =   Column profiles
% 		out.K         =   Maximun number of dimensions
% 		out.k         =   Number of retained dimensions
% 		out.A         =   Standarized residuals
% 		out.chis      =   Chi-square statistic for the contingency table
% 		out.ti        =   Total inertia
% 		out.CV        =   Cramer's V
% 		out.eig       =    Principal inertias [Raw Percentage Acumulate_percentage]
%       out.aev       =    Average explained variance (dimensions explaining less variance should be excluded from the map)
% 		out.F         =    Principal coordinates of rows
% 		out.G         =    Principal coordinates of columns
% 		out.X         =    Standard coordinates of rows
% 		out.Y         =    Standard coordinates of columns
% 		out.H         =    Symmetrical coordinates of rows
% 		out.Z         =    Symmetrical coordinates of columns
%
% 		out.ARC       =    Absolute row contributions: Contributions of points to axes.
% 		out.ACC       =    Absolute column contributions: Contributions of points to axes.
% 		out.INr       =    Inertia of each row [raw percentage]
% 		out.INc       =    Inertia of each column [raw percentage]
% 		out.RRC       =    Relative row contributions: Contributions of axes to points; squared correlations of points with axes
% 		out.RCC       =    Relative column contributions: Contributions of axes to points; squared correlations of points with axes
% 		out.QLTr      =    Measure of quality of display for each row in the reduced k-dimensional map. The qualities are equivalent to the communalities in PCA
% 		out.QLTc      =    Measure of quality of display for each column in the reduced k-dimensional map. The qualities are equivalent to the communalities in PCA
% 		out.Br        =    Correlation of rows with axes. Similar to component loadings in PCA
% 		out.Bc        =    Correlation of columns with axes. Similar to component loadings in PCA
%
% 		out.Fsup      =    Principal coordinates of supplementary columns
% 		out.Xsup      =    Standard coordinates of supplementary columns
% 		out.Hsup      =    Symmetrical coordinates of supplementary columns
%
% 		out.Gsup     =     Principal coordinates of supplementary columns
% 		out.Ysup      =    Standard coordinates of supplementary columns
% 		out.Zsup     =     Symmetrical coordinates of supplementary columns
%
%                   For a given axis, the standard and principle co-ordinates are related as follows:
%                         P = sqrt(eigenvalue) X S
%                         S: the standard coordinate
%                         P: the principal coordinate of a row (or a column) on the axis
%                eigenvalue: the eigenvalue of the axis
%
% See also crosstab
%
% References:
%
%   Greenacre.....
%   Wiley.
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('CorAna')">Link to the help function</a>
% Last modified 14-06-2016
%
%

%% Beginning of code

% The input is certainly a
UserOptions=varargin(1:2:length(varargin));

if ~isempty(varargin)
    checkdatamatrix = strcmp(UserOptions,'datamatrix');
    if sum(checkdatamatrix)
        datamatrix = varargin{2*find(checkdatamatrix)};
    else
        datamatrix='';
    end
else
    datamatrix='';
end

% If input is a datamatrix it is necessary to construct the contingency
% table
if datamatrix==1
    [N,chi2,p,labels] =crosstab(N(:,1),N(:,2));
    [I,J]=size(N);
    % default labels for rows of contingency table
    Lr=labels(1:I,1);
    % default labels for columns of contingency table
    Lc=labels(1:J,2);
    
else
    [I,J]=size(N);
    
    %Size of N
    % default labels for rows of contingency table
    Lr=cellstr(strcat('r',num2str((1:I)')));
    % default labels for columns of contingency table
    Lc=cellstr(strcat('c',num2str((1:J)')));
end


out.N = N;

% default value for supplementary units
Sup='';
k=2;
plots=1;

% Dimensions to show in the plot. The default is to show the first two
% dimensions

options=struct('k',k,'Sup',Sup,'plots',plots,'datamatrix','');
options.Lr=Lr;
options.Lc=Lc;

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:CorAna:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    k=options.k;
    Lr=options.Lr;
    Lc=options.Lc;
    Sup=options.Sup;
    plots=options.plots;
end

% Extract labels for rows and columns
if istable(N)
    Lc=N.Properties.VariableNames;
    Lr=N.Properties.RowNames;
    N=table2array(N);
else
    
    if isempty(Lr)
        Lr=cellstr(num2str((1:I)'));
    else
        % Check that the length of Lr is equal to I
        if length(Lr)~=I
            error('Wrong length of row labels')
        end
    end
    
    if isempty(Lc)
        Lc=cellstr(num2str((1:J)'));
    else
        % Check that the length of Lc is equal to J
        if length(Lc)~=J
            error('Wrong length of column labels')
        end
    end
end

% Nred will contain the contingency table after removing supplementary rows and
% columns (if supplementary rows and columns belong to the table)
Nred=N;

labels=struct;

if ~isempty(Sup)
    % if Sup.r (Sup.c) is a cell or is character or it is a numeric
    % vector, then the supplementary rows (columns) belong to the the
    % actual contingency table, else if Sup.r (Sup.c) is a Matlab table
    % then the supplementary units do not belong to the actual
    % contingency table N
    
    if isfield(Sup,'r')
        out.SupRow = Sup.r;
        if iscell(Sup.r) || ischar(Sup.r)
            % find the indexes of the rows of matrix N to delete (rows to use as supplementary rows)
            Indexesr=zeros(length(Sup.r),1);
            for i=1:length(Sup.r)
                if iscell(Sup.r)
                    Indexesr(i)=find(strcmp(Sup.r{i},Lr),1);
                else
                    Indexesr(i)=find(strcmp(Sup.r(i),Lr),1);
                end
            end
            
            labels.sr=Lr;
            labels.sr = labels.sr(Indexesr);
            % Delete the labels of contingency table associated to supplementary rows
            Lr(Indexesr)=[];
            
            
        elseif istable(Sup.r)
            Indexesr='';
            Nsupr=table2array(Sup.r);
            labels.sr=Sup.r.Properties.RowNames;
        else
            
            
            Indexesr=Sup.r;
            if min(Indexesr)<1 || max(Indexesr)> size(N,1)
                error('FSDA:CorAna:wrongInput',['Numeric indexes of supplementary rows must be integers between 1 and ' num2str(size(N,1))])
            end
            labels.sr=Lr;
            labels.sr = labels.sr(Indexesr);
            
            % Delete the labels of contingency table associated to supplementary rows
            Lr(Indexesr)=[];
            
        end
        
        
        out.labels.sr = labels.sr;
    else
        Indexesr='';
        
    end
    
    if isfield(Sup,'c')
        out.SupCol = Sup.c;
        if iscell(Sup.c) || ischar(Sup.c)
            % find the indexes of the rows to delete (rows to use as supplementary rows)
            Indexesc=zeros(length(Sup.c),1);
            for i=1:length(Sup.c)
                if iscell(Sup.c)
                    Indexesc(i)=find(strcmp(Sup.c{i},Lc),1);
                else
                    Indexesc(i)=find(strcmp(Sup.c(i),Lc),1);
                end
            end
            labels.sc=Lc;
            labels.sc = labels.sc(Indexesc);
            
            % Delete the labels of contingency table associated to supplementary rows
            Lc(Indexesc)=[];
            
        elseif istable(Sup.c)
            Indexesc='';
            Nsupc=table2array(Sup.c);
            labels.sc=Sup.c.Properties.VariableNames;
        else
            
            Indexesc=Sup.c;
            if min(Indexesc)<1 || max(Indexesc)> size(N,2)
                error('FSDA:CorAna:wrongInput',['Numeric indexes of supplementary columns must be integers between 1 and ' num2str(size(N,2))])
            end
            
            
            labels.sc=Lc;
            labels.sc = labels.sc(Indexesc);
            % Delete the labels of contingency table associated to supplementary rows
            Lc(Indexesc)=[];
            
        end
        
        
        out.labels.sc = labels.sc;
    else
        Indexesc='';
    end
    
    
    % if ~isempty(Indexesr) this means that supplementary rows belong to
    % matrix N
    if ~isempty(Indexesr)
        % Contingency table referred to supplementary rows.
        if ~isempty(Indexesc)
            Nsupr=N(Indexesr,Indexesc);
        else
            Nsupr=N(Indexesr,:);
        end
    end
    
    % if ~isempty(Indexesc) this means that supplementary columns belong to
    % matrix N
    if ~isempty(Indexesc)
        % Contingency table referred to supplementary columns.
        if ~isempty(Indexesr)
            Nsupc=N(Indexesr,Indexesc);
        else
            Nsupc=N(:,Indexesc);
        end
    end
    
    % Delete the rows and columns of contingency table associated to
    % supplementary rows and/or associated to supplementary columns
    Nred(Indexesr,:)=[];
    Nred(:,Indexesc)=[];
    
end


[I,J]=size(Nred);
out.I=I;            %number of active rows (excluding supplementary rows)
out.J=J;            %number of active columns (excluding supplementary columns)

%vectors of ones
onesI = ones(I,1);
onesJ = ones(J,1);

%grand total
n=sum(sum(Nred));
out.n=n;            %grand total

%correspondence matrix P
P = (1/n) * Nred;
out.P=P;            %Correspondence matrix (proportions)
out.Lr=Lr;
out.Lc=Lc;
out.labels=labels;

% r= vector which contains row masses
r  = P * onesJ ;
Dr = diag(r);
out.r = r;          %row masses (vector)
out.Dr = Dr;        %row masses (diagonal matrix)

%column masses
% Centroid of the row profiles   r' * Dr^(-1) * P = 1' * P = c'

c  = (onesI' * P)';
Dc = diag(c);
out.c = c;          %column masses (vector); c is the centroid of row profiles
out.Dc = Dc;        %column masses (diagonal matrix)


%Rows profiles (equation 4.14)
RP = Dr^(-1) * P;
out.RP = RP;        %rows profiles

%Columns profiles  (equation 4.14)
CP = P * Dc^(-1);
out.CP = CP;        %column profiles

%Number of dimensions

K = min(I-1,J-1);
out.K = K;          %maximun number of dimensions

% k= number of retained dimensions
out.k = k;

% Standarized residuals
%
% A = Dr^(1/2) * (Dr^(-1) * P - Ir * c') * Dc^(-1/2) = Dr^(-1/2) * (P - r * c') * Dc^(-1/2);  aij = sqrt( (pij - ricj)^2 /ricj ) =(pij - ricj)/sqrt(ricj)
% Chi-square distances
%
A =  Dr^(-1/2) * (P - r * c') * Dc^(-1/2);
out.A= A;           % Standarized residuals

%
% SVD of A = U*Gam*V'

[U,Gam,V] = svd(A);
Gam=Gam(1:K,1:K);
U=U(:,1:K);
V=V(:,1:K);

% Chi-square statistic for the contingence table
%
%  chis = n * sum_i sum_j (pij - ricj)^2 / ricj
%
chis = onesI'* (A.*A) * onesJ * n;
out.chis = chis;    % Chi-square statistic for the contingence table

% Total inertia
%
% ti = sum_i sum_j (pij - ricj)^2 / ricj = chis/n
%
ti = chis/n;
out.ti = ti;        %total inertia

out.CV = sqrt(chis/(n*(min(I,J)-1)));   %Cramer's V

% Principal inertias
%
% eigenvalues of A'A
% percentages of the total inertia

% api = cumulative proportion of explained inertia
api=cumsum(diag(Gam.^2))/ti;

out.eig = [diag(Gam.*Gam) diag(Gam.*Gam / ti) api] ;      % Principal inertias [Raw Percentage Acumulate_percentage]

out.aev = 100/K;                        %Average explained variance (dimensions explaining less variance should be excluded from the map)

% Principal coordinates of rows  (alpha=1 for the rows)
RowsPri = Dr^(-1/2) * U*Gam;
out.F = RowsPri;                               % Principal coordinates of rows

% Principal coordinates of columns (Dc^(-1/2)*Gamma*v)
ColsPri = Dc^(-1/2) * V*Gam;
out.G = ColsPri;                               % Principal coordinates of columns

% Standard coordinates of rows
% X = F * Gam^(-1) = Dr^(-1/2) * U * Gam * Gam^(-1) = Dr^(-1/2) * U
RowsSta = Dr^(-1/2) * U ;
out.RowsSta = RowsSta;                               % Standard coordinates of rows

% Standard coordinates of columns
% Y = G * Gam^(-1) = Dc^(-1/2) * V * Gam * Gam^(-1) = Dc^(-1/2) * V
ColsSta = Dc^(-1/2) * V ;
out.ColsSta=ColsSta;                                 % Standard coordinates of columns

% Symmetrical coordinates of rows
RowsSym = Dr^(-1/2) * U*Gam^(1/2);
out.RowsSym = RowsSym;                               % Symmetrical coordinates of rows

% Symmetrical coordinates of columns
ColsSym = Dc^(-1/2) * V*Gam^(1/2);
out.ColsSym=ColsSym;                                 % Symmetrical coordinates of columns


% Contributions of rows to inertia
ARC = Dr*RowsSta.*RowsSta;
out.ARC = ARC;                          % Absolute row contributions: Contributions of points to axes.

% Contributions of columns to inertia
ACC = Dc*ColsSta.*ColsSta;
out.ACC = ACC;                          % Absolute column contributions: Contributions of points to axes.

%Inertia of each row
INr = Dr*(RowsPri.*RowsPri)*ones(K,1);
%Inertia relative to their total (used in the output)
INr_p=INr / trace(Gam.*Gam);
out.INr = [INr INr_p];                   %Inertia of each row [raw percentage]

%Inertia of each column
INc = Dc*(ColsPri.*ColsPri)*ones(K,1);
%Inertia relative to their total (used in the output)
INc_p=INc / trace(Gam.*Gam);
out.INc = [INc INc_p];                   %Inertia of each column [raw percentage]

%Squared correlations of rows with axes
QCOR_r = diag(INr)^(-1)*Dr*(RowsPri.*RowsPri);
out.RRC=QCOR_r;                          %Relative row contributions: Contributions of axes to points; squared correlations of points with axes

%Squared correlations of columns with axes
QCOR_c = diag(INc)^(-1)*Dc*(ColsPri.*ColsPri);
out.RCC=QCOR_c;                          %Relative columns contributions: Contributions of axes to points; squared correlations of points with axes

%Measure of quality of display for each row in the reduced k-dimensional map. The qualities are equivalent to the communalities in PCA
QLTr = QCOR_r(:,1:k)*ones(k,1);
out.QLTr = QLTr;                         %Measure of quality of display for each row in the reduced k-dimensional map. The qualities are equivalent to the communalities in PCA

%Measure of quality of display for each column in the reduced k-dimensional map. The qualities are equivalent to the communalities in PCA
QLTc = QCOR_c(:,1:k)*ones(k,1);
out.QLTc = QLTc;                         %Measure of quality of display for each column in the reduced k-dimensional map. The qualities are equivalent to the communalities in PCA

%Correlation of rows with axes. Similar to component loadings in PCA
Br = sign(RowsPri).*sqrt(QCOR_r);
out.Br = Br;                             %Correlation of rows with axes. Similar to component loadings in PCA

%Correlation of columns with axes. Similar to component loadings in PCA
Bc = sign(ColsPri).*sqrt(QCOR_c);
out.Bc = Bc;                             %Correlation of columns with axes. Similar to component loadings in PCA



if exist('Sup','var')
    %Supplementary rows
    if isfield(Sup,'r')
        % The sum of each row of h must be equal to 1
        % h=Nsup(Indexesr,:)/(diag(sum(Nsup(Indexesr,:))));
        h=bsxfun(@rdivide,Nsupr,sum(Nsupr,2));
        
        RowsPrisup=h*ColsSta;                                %Principal coordinates of supplementary rows
        RowsStasup=h*ColsSta*Gam^(-1);                       %Standard coordinates of supplementary rows
        RowsSymsup=h*ColsSta*Gam^(-1/2);                     %Symmetrical coordinates of supplementary rows
        out.Fsup=RowsPrisup;
        out.Xsup=RowsStasup;
        out.Hsup=RowsSymsup;
        %       rrc=(Gsup(:,1:k).*Gsup(:,1:k))/trace(Gsup(:,1:k)'*Gsup(:,1:k));
    end
    
    %Supplementary columns
    if isfield(Sup,'c')
        % The sum of each column of h must be equal to 1
        h=Nsupc/(diag(sum(Nsupc)));
        ColsPrisup=h'*RowsSta;                              %Principal coordinates of supplementary columns
        ColsStasup=h'*RowsSta*Gam^(-1);                     %Standard coordinates of supplementary columns
        ColsSymsup=h'*RowsSta*Gam^(-1/2);                   %Symmetrical coordinates of supplementary columns
        out.Gsup=ColsPrisup;
        out.Ysup=ColsStasup;
        out.Zsup=ColsSymsup;
        %     rrc=(Gsup(:,1:k).*Gsup(:,1:k))/trace(Gsup(:,1:k)'*Gsup(:,1:k));
    end
end

if isstruct(plots) || plots==1
    FontName='Times';
    FontSize=12;
    
    d1=1;
    d2=2;
    if isstruct(plots)
        % This anonymous function anables to extract the variable name to a
        % string
        ExtractVariableName=@(x) inputname(1);
        
        if isfield(plots,'alpha')
            if strcmp(plots.alpha,'rowprincipal')
                typeR='RowsPri'; % rows are in principal coordinates
                typeC='ColsSta';        % columns are in standard coordinates
                titl='Rows principal coordinates, and column standard coordinates  $\alpha=1$,$X=D_r^{-1/2}U\Gamma $ and $Y= D_c^{-1/2} V$';
                
            elseif strcmp(plots.alpha,'colprincipal')
                typeR='RowsSta'; % rows are in standard coordinates
                typeC='ColsPri';        % columns are in principal coordinates
                titl='Rows standard coordinates, and column principal coordinates $\alpha=0$, $X=D_r^{-1/2}U $ and $G= D_c^{-1/2} V \Gamma$';
                
            elseif strcmp(plots.alpha,'symbiplot')
                % equivalent to alpha=0.5
                typeR='RowsSym';        % rows are in symmetrical coordinates
                typeC='RowsSym';        % columns are in symmetrical coordinates
                titl='Biplot symmetrical model $\alpha=0.5$ $X=D_r^{-1/2}U\Gamma^{1/2} $ and $Y= D_c^{-1/2} \Gamma V^{1/2}$';
                
            elseif strcmp(plots.alpha,'bothprincipal')
                typeR='RowsPri';        % rows are in principal coordinates
                typeC='ColsPri';        % columns are in principal coordinates
                titl='French symmetrical model: rows and cols in principal coordinates. Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$';
                
            elseif strcmp(plots.alpha,'rowgab')
                %  If plots.alpha='rowgab'  rows are in principal coordinates
                %  and columns are in standard coordinates multiplied by the
                %  mass.
                typeR='RowsPri';        % rows are in principal coordinates
                ColsStaDc=Dc*ColsSta;
                typeC=ExtractVariableName(ColsStaDc);
                titl='Rows principal coordinates, and column standard coordinates times masses  $\alpha=1$,$X=D_r^{-1/2}U\Gamma $ and $Y= D_c^{1/2} V$';
                
                
            elseif strcmp(plots.alpha,'colgab')
                % If plots.alpha='colgab'  columns are in principal coordinates
                % and rows are in standard coordinates multiplied by the
                % masses.
                RowsStaDr=Dr*RowsSta;
                typeR=ExtractVariableName(RowsStaDr);
                typeC='ColsPri';        % columns are in principal coordinates
                titl='Rows standard coordinates multiplied by masses, and column principal coordinates , $X=D_r^-1/2}U $ and $Y= D_c^{-1/2} V \Gamma$';
                
                
            elseif strcmp(plots.alpha,'rowgreen')
                %  If plots.alpha='rowgreen'  rows are in principal
                %  coordinates and columns are in standard coordinates
                %  multiplied by square root of the mass.
                typeR='RowsPri';        % rows are in principal coordinates
                ColsStaDcSqrt=(Dc^(1/2))*ColsSta;
                typeC= ExtractVariableName(ColsStaDcSqrt);
                titl='Rows principal coordinates, and column standard coordinates times sqrt of masses $X=D_r^{-1/2}U\Gamma $ and $Y= V$';
                
                
            elseif strcmp(plots.alpha,'colgreen')
                %  If plots.alpha='colgreen' columns in principal coordinates
                %  and rows in standard coordinates multiplied by the square
                %  root of the mass.
                RowsStaDrSqrt=(sqrt(Dr))*RowsSta;
                typeR=ExtractVariableName(RowsStaDrSqrt);
                typeC='ColsPri';        % columns are in principal coordinates
                titl='Rows standard coordinates times sqrt of masses, and column principal coordinates, $X=U $ and $G= D_c^{-1/2} V \Gamma$';
                
            else
                if isnumeric(plots.alpha)
                    if plots.alpha>=0 && plots.alpha<=1
                        RowsAlpha= Dr^(-1/2) * U*Gam^plots.alpha;
                        ColsAlpha= Dc^(-1/2) * V*Gam^(1-plots.alpha);
                        typeR=ExtractVariableName(RowsAlpha);
                        typeC=ExtractVariableName(ColsAlpha);
                        titl=['$\alpha='  num2str(plots.alpha) '\qquad   X=D_r^{-1/2} U \Gamma^{' num2str(plots.alpha) '}$ and $Y= D_c^{-1/2} V \Gamma^{1-'  num2str(plots.alpha) '}$'];
                    else
                        error('Value of plots.alpha must lie in the interval [0 1]')
                    end
                else
                    listStrings={'rowprincipal'; 'colprincipal'; 'symbiplot'; 'bothprincipal'; 'rowgab'; 'colgab'; 'rowgreen'; 'colgreen'};
                    warning(['Input string ''' plots.alpha ''' is  not found'])
                    disp('Possible strings are')
                    disp(listStrings)
                    error('Please use one of the above strings')
                end
            end
        end
    else
        typeR='RowsPri';        % rows are in principal coordinates
        typeC='ColsPri';        % columns are in principal coordinates
        titl='';
    end
    symbolrows='o';
    symbolcols='^';
    symbolsuprows='o';
    symbolsupcols='^';
    colorrows='b';
    colorcols='r';
    colorsuprows='b';
    colorsupcols='r';
    
    d1str=num2str(d1);
    d2str=num2str(d2);
    
    % case 1
    figure
    hold('on')
    % Plot row points
    propR=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolrows ,''',''Color'',''', colorrows , '''');
    eval(['plot(' typeR '(:,' d1str '),' typeR '(:,' d2str '),' propR ')'])
    % plot(H(:,d1),H(:,d2),'LineStyle','none','Marker',symbolrows,'Color',colorrows)
    
    % Plot column points
    propC=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolcols ,''',''Color'',''', colorcols, '''');
    eval(['plot(' typeC '(:,' d1str '),' typeC '(:,' d2str '),' propC ')'])
    % plot(Z(:,d1),Z(:,d2),'LineStyle','none','Marker',symbolcols,'Color',colorcols)
    
    % Add labels for row points and column points
    %     text(H(:,d1),H(:,d2),Lr)
    %     text(Z(:,d1),Z(:,d2),Lc)
    eval(['text(' typeR '(:,' d1str '),' typeR '(:,' d2str '),Lr,''Interpreter'',''None'')'])
    eval(['text(' typeC '(:,' d1str '),' typeC '(:,' d2str '),Lc,''Interpreter'',''None'')'])
    
    title(titl,'Interpreter','Latex');
    
    % Labels for axes
    xlabel(['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',out.eig(d1,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    ylabel(['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',out.eig(d2,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    
    % Add points associated to supplementary rows
    if isstruct(Sup) && isfield(Sup,'r')
        propsupR=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolsuprows ,''',''Color'',''', colorsuprows , ''',''MarkerFaceColor'',''', colorsuprows ,'''');
        eval(['plot(' typeR 'sup(:,d1),' typeR 'sup(:,d2),' propsupR ')'])
        eval(['text(' typeR 'sup(:,d1),' typeR 'sup(:,d2),labels.sr,''Interpreter'',''None'')'])
        
        %         plot(Hsup(:,d1),Hsup(:,d2),'Marker',symbolsuprows,'Color',colorsuprows,'MarkerFaceColor',colorsuprows)
        %         text(Hsup(:,d1),Hsup(:,d2),Ls.r)
    end
    
    % Add points associated to supplementary columns
    if isstruct(Sup) && isfield(Sup,'c')
        propsupC=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolsupcols ,''',''Color'',''', colorsupcols , ''',''MarkerFaceColor'',''', colorsuprows ,'''');
        eval(['plot(' typeC 'sup(:,d1),' typeC 'sup(:,d2),' propsupC ')'])
        eval(['text(' typeC 'sup(:,d1),' typeC 'sup(:,d2),labels.sc,''Interpreter'',''None'')'])
        %         plot(Zsup(:,d1),Zsup(:,d2),'Marker',symbolsuprows,'Color',colorsuprows,'MarkerFaceColor',colorsuprows)
        %         text(Zsup(:,d1),Zsup(:,d2),Ls.r)
    end
    
    % Make axis equal and add cartesian axes
    axis(gca,'equal')
    axis(gca,'equal')
    vv=axis;
    line([vv(1);vv(2)],[0;0])
    line([0;0],[vv(3);vv(4)])
    
    %{
    % case 2
    figure
    hold('on')
    plot(F(:,d1),F(:,d2),'Marker',symbolrows,'Color',colorrows)
    plot(Y(:,d1),Y(:,d2),'Marker',symbolcols,'Color',colorcols)
    axis(gca,'equal')
    vv=axis;
    line([vv(1);vv(2)],[0;0])
    line([0;0],[vv(3);vv(4)])
    text(F(:,d1),F(:,d2),Lr)
    text(Y(:,d1),Y(:,d2),Lc)
    title('Row principal coordinates, and column standard coordinates  $\alpha=1$,$X=D_r^{-1/2}U\Gamma $ and $Y= D_c^{-1/2} V$','Interpreter','Latex');
    % Labels for axes
    xlabel(['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',output.eig(d1,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    ylabel(['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',output.eig(d2,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    
    % Add points associated to supplementary rows
    if isstruct(Sup) && isfield(Sup,'r')
        plot(Fsup(:,d1),Fsup(:,d2),'Marker',symbolsuprows,'MarkerFaceColor','auto')
        text(Fsup(:,d1),Fsup(:,d2),Ls.r)
    end
    
    
    % case 3
    figure
    hold('on')
    plot(X(:,d1),X(:,d2),'Marker',symbolrows,'Color',colorrows)
    plot(G(:,d1),G(:,d2),'Marker',symbolcols,'Color',colorcols)
    axis(gca,'equal')
    vv=axis;
    line([vv(1);vv(2)],[0;0])
    line([0;0],[vv(3);vv(4)])
    text(X(:,d1),X(:,d2),Lr)
    text(G(:,d1),G(:,d2),Lc)
    title('Row standard coordinates, and column principal coordinates $\alpha=0$, $X=D_r^{-1/2}U $ and $G= D_c^{-1/2} V \Gamma$','Interpreter','Latex');
    % Labels for axes
    xlabel(['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',output.eig(d1,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    ylabel(['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',output.eig(d2,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    % Add points associated to supplementary rows
    if isstruct(Sup) && isfield(Sup,'r')
        plot(Xsup(:,d1),Xsup(:,d2),'Marker',symbolsuprows,'Color',colorsuprows,'MarkerFaceColor','auto')
        text(Xsup(:,d1),Xsup(:,d2),Ls.r)
    end
    

    
    % case 4
    figure
    hold('on')
    plot(F(:,d1),F(:,d2),'Marker',symbolrows,'Color',colorrows)
    plot(G(:,d1),G(:,d2),'Marker',symbolcols,'Color',colorcols)
    axis(gca,'equal')
    vv=axis;
    line([vv(1);vv(2)],[0;0])
    line([0;0],[vv(3);vv(4)])
    text(F(:,d1),F(:,d2),Lr)
    text(G(:,d1),G(:,d2),Lc)
    title('French symmetrical model: plot of $F=D_r^{-1/2}U \Gamma$ and $G= D_r^{-1/2} V \Gamma$','Interpreter','Latex')
    % Labels for axes
    xlabel(['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',output.eig(d1,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    ylabel(['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',output.eig(d2,2)*100),'%)'],'FontName', FontName, 'FontSize', FontSize);
    
    % Add points associated to supplementary rows
    if isstruct(Sup) && isfield(Sup,'r')
        plot(Fsup(:,d1),Fsup(:,d2),'Marker',symbolsuprows,'Color',colorsuprows,'MarkerFaceColor','auto')
        text(Fsup(:,d1),Fsup(:,d2),Ls.r)
    end
    %}
    %    Symetric plot represents the row and column profiles simultaneously
    %    in a common space (Bendixen, 2003). In this case, only the distance
    %    between row points or the distance between column points can be
    %    really interpreted.
    %   The distance between any row and column items is not meaningful! You can
    %   only make a general statements about the observed pattern. In order to
    %   interpret the distance between column and row points, the column profiles
    %   must be presented in row space or vice-versa. This type of map is called
    %   asymmetric biplot
end

dispresults=1;
if dispresults==1
    disp('ROW POINTS')
    disp(['Results for dimension: ' d1str])
    Tabresults=array2table(eval(strcat('[',typeR,'(:,', d1str ,') ARC(:,', d1str ,')    QCOR_r(:,', d1str ,')         ]')));
    Tabresults.Properties.RowNames=Lr;
    Tabresults.Properties.VariableNames={'Scores', 'ctr' 'cos2'};
    disp(Tabresults)
    
    disp(['Results for dimension: ' d2str])
    Tabresults=array2table(eval(strcat('[',typeR,'(:,', d2str ,') ARC(:,', d2str ,')    QCOR_r(:,', d2str ,')         ]')));
    Tabresults.Properties.RowNames=Lr;
    Tabresults.Properties.VariableNames={'Scores', 'ctr' 'cos2'};
    disp(Tabresults)
    
    disp('COLUMN POINTS')
    disp(['Results for dimension: ' d1str])
    Tabresults=array2table(eval(strcat('[',typeC,'(:,', d1str ,') ACC(:,', d1str ,')    QCOR_c(:,', d1str ,')         ]')));
    Tabresults.Properties.RowNames=Lc;
    Tabresults.Properties.VariableNames={'Scores', 'ctr' 'cos2'};
    disp(Tabresults)
    
    disp(['Results for dimension: ' d2str])
    Tabresults=array2table(eval(strcat('[',typeC,'(:,', d2str ,') ACC(:,', d2str ,')    QCOR_c(:,', d2str ,')         ]')));
    Tabresults.Properties.RowNames=Lc;
    Tabresults.Properties.VariableNames={'Scores', 'ctr' 'cos2'};
    disp(Tabresults)
    
    %TODO
    %     disp('SUPPLEMENTARY ROW POINTS')
    %     disp(['Results for dimension: ' d1str])
    %    if isstruct(Sup) && isfield(Sup,'r')
    %         propsupR=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolsuprows ,''',''Color'',''', colorsuprows , ''',''MarkerFaceColor'',''', colorsuprows ,'''');
    %         eval(['plot(' typeR 'sup(:,d1),' typeR 'sup(:,d2),' propsupR ')'])
    %         eval(['text(' typeR 'sup(:,d1),' typeR 'sup(:,d2),labels.sr,''Interpreter'',''None'')'])
    %
    %
    %        disp(['Results for dimension: ' d1str])
    %     Tabresults=array2table(eval(strcat('[',typeR,'(:,', d1str ,') ARC(:,', d1str ,')    QCOR_r(:,', d1str ,')         ]')));
    %     Tabresults.Properties.RowNames=Lr;
    %     Tabresults.Properties.VariableNames={'Scores', 'ctr' 'cos2'};
    %     disp(Tabresults)
    %
    %     disp(['Results for dimension: ' d2str])
    %     Tabresults=array2table(eval(strcat('[',typeR,'(:,', d2str ,') ARC(:,', d2str ,')    QCOR_r(:,', d2str ,')         ]')));
    %     Tabresults.Properties.RowNames=Lr;
    %     Tabresults.Properties.VariableNames={'Scores', 'ctr' 'cos2'};
    %     disp(Tabresults)
    %
    %    end
    
    % ctr= relative contribution or each row to the inertia of the
    % dimension.
    % The rows  with the larger value of ctr are those which contribute the
    % most to the definition of the dimensions. If the row contributions
    % were uniform, the expected value would be 1/nrow(contingeny table)
    % For a given dimension, any row with a
    % contribution larger than this threshold could be considered as
    % important in contributing to that dimension.
    
    % cos2= contribution of dimension to the inertia of point
    % (where inertia of point is is the squared distance of point d_i^2 to the centroid)
    % The quality of representation of the rows on the factor map is called
    % the squared cosine (cos2) or the squared correlations.
    % The cos2 measures the degree of association between rows/columns and
    % a particular axis.
    % The values of the cos2 are comprised between 0 and 1
    % If a row item is well represented by two dimensions, the sum of the
    % cos2 is closed to one.
end

end
