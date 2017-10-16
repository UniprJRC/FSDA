function out=CorAna(N, varargin)
%CorAna performs correspondence analysis
%
%<a href="matlab: docsearchFS('CorAna')">Link to the help function</a>
%
%  Required input arguments:
%
%       N    :    Contingency table (default) or n-by-2 input datasets.
%                 Matrix or Table.
%                 Matrix or table which contains the input contingency
%                 table (say of size I-by-J) or the original data matrix.
%                 In this last case N=crosstab(N(:,1),N(:,2)). As default
%                 procedure assumes that the input is a contingency table.
%
%  Optional input arguments:
%
%       k    :  Number of dimensions to retain. Scalar.
%               Scalar which contains the number of dimensions to retain.
%               The default value of k is 2.
%               Example - 'k',2
%               Data Types - double
%       Lr   :  Vector of row labels. Cell.
%               Cell containing the labels of the rows of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table. because in this case  Lr=N.Properties.RowNames;
%               Example - 'Lr',{'a' 'b' 'c'}
%               Data Types - cell array of strings
%       Lc   :  Vector of column labels. Cell.
%               Cell containing the labels of the columns of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table because in this case Lc=N.Properties.VariableNames;
%               Example - 'Lc',{'c1' c2' 'c3' 'c4'}
%               Data Types - cell array of strings
%       Sup :  Structure containing indexes or names of supplementary rows
%              or columns. Structure.  Structure with the followin fields.
%              Sup.r =  numeric vector containing row indexes  or cell
%                       array of strings or table containing supplementary
%                       rows. If indexes or cell array of strings are
%                       supplied, we assume that supplementary rows belong
%                       to contingency table N. For example, if Sup.r=[2 5]
%                       (that is Sup.r is a numeric vector which contains
%                       row indexes) we use rows 2 and 5 of the input
%                       contingency table as supplementary rows. For
%                       example, if Sup.r={'Junior-Managers'
%                       'Senior-Employees'} (that is Sup.r is a cell array
%                       of strings) we use rows named 'Junior-Managers' and
%                       'Senior-Employees' of the input contingency table
%                       as supplementary rows. Of course the length of
%                       Sup.r must be smaller than the number of rows of
%                       the contigencey matrix divided by 2. If, on the
%                       other hand, Sup.r is a table supplementary rows do
%                       not belong to N.
%               Sup.c = numeric vector containing column indexes or cell
%                       array of string containing names of the columns to
%                       use as supplementary columns, or table.
%                       If indexes or cell array of strings are supplied,
%                       we assume that supplementary columns belong to
%                       contingency table N.
%                       For example, if Sup.c=[2 3] (that is Sup.c is a
%                       numeric vector which contains column indexes) we use
%                       columns 2 and 3 of the input contingency table as
%                       supplementary columns.
%                       For example, if Sup.c={'Smokers' 'NonSmokers'}
%                       (that is Sup.c is a cell array of strings) we use
%                       columns of the contingency table labelled 'Smokers'
%                       and 'NonSmokers' of the input contingency table N
%                       as supplementary columns.
%                       Of course the length of Sup.c must be smaller than
%                       the number of columns of the contigencey matrix
%                       divided by 2.
%                       If, on the other hand, Sup.c is a table
%                       supplementary columns  do not belong to N.
%                       Example - 'Sup',Sup=struct; Sup.c={'c2' 'c4'}
%                       Data Types - struct
%                       REMARK: The default value of Sup is a missing value
%                       that is we assume that there are no supplementary
%                       rows or columns.
%      datamatrix   :  data matrix or contingency table. Boolean. If
%                       datamatrix is true the first input argument N is
%                       forced to be interpreted as a data matrix, else if
%                       the input argument is false N is treated as a
%                       contingency table. The default value of datamatrix
%                       is false, that is the procedure automatically
%                       considers N as a contingency table
%               Example - 'datamatrix',true
%               Data Types - logical
%       plots : Plot on the screen. Scalar or structure.
%               If plots = 1, a plot which shows  the Principal
%               coordinates of rows and columns is shown on the screen. If
%               plots is a structure it may contain the following fields:
%               plots.alpha = type of plot. scalar in the interval [0 1] or
%               a string identifyind the type of coordinates to use in the
%               plot.
%               If $plots.alpha='rowprincipal'$ the row points are in
%                   principal coordinates and the column coordinates are
%                   standard coordinates. Distances between row points are
%                   (approximated) chi-squared distances
%                   (row-metric-preserving). The position of the row points
%                   are at the weighted average of the column points.
%                   Note that 'rowprincipal' can also be specified setting
%                   plots.alpha=1.
%               If $plots.alpha='colprincipal'$, the column
%                   coordinates are referred to as principal coordinates
%                   and the row coordinates as standard coordinates.
%                   Distances between column points are (approximated)
%                   chi-squared distances (row-metric-preserving). The
%                   position of the column points are at the weighted
%                   average of the row points.
%                   Note that 'colwprincipal' can also be
%                   specified setting plots.alpha=0.
%               If $plots.alpha='symbiplot'$, the row and column coordinates
%                   are scaled similarly. The sum of weighted squared
%                   coordinates for each dimension is equal to the
%                   corresponding singular values.  These coordinates are often
%                   called symmetrical coordinates. This representation is
%                   particularly useful if one is primarily interested in
%                   the relationships between categories of row and column
%                   variables rather than in the distances among rows or
%                   among columns. 'symbiplot' can also be specified
%                   setting plots.alpha=0.5;
%               If $plots.alpha='bothprincipal'$, both the rows and columns
%                   are depicted in principal coordinates. Such a plot is
%                   often referred to as a symmetrical plot or French
%                   symemtrical model. Note that such a symmetrical plot
%                   does not provide a feasible solution in the sense that
%                   it does not approximate matrix
%                   $D_r^{-0.5}(P-rc')D_c^{-0.5}$.
%               If $plots.alpha='rowgab'$, rows are in principal coordinates
%                   and columns are in standard coordinates multiplied by
%                   the mass. This biplot has been suggested by Gabriel and
%                   Odoroff (1990).
%               If $plots.alpha='colgab'$, columns are in principal coordinates
%                   and rows are in standard coordinates multiplied by the
%                   mass. This biplot has been suggested by Gabriel and
%                   Odoroff (1990).
%               If $plots.alpha='rowgreen'$, rows are in principal
%                   coordinates and columns are in standard coordinates
%                   multiplied by square root of the mass.
%                  This biplot has been suggested by Greenacre and
%                  incorporates the contribution of points. In this
%                  display, points that contribute very little to the
%                  solution, are close to the center of the biplot and are
%                  relatively unimportant to the interpretation. This
%                  biplot is often referred as contribution biplot because
%                  it visually shows the most contributing points
%                  (Greenacre 2006b).
%               If $plots.alpha='colgreen'$, columns in principal coordinates
%                   and rows in standard coordinates multiplied by the
%                   square root of the mass.
%                  This biplot has been suggested by Greenacre and
%                  incorporates the contribution of points. In this
%                  display, points that contribute very little to the
%                  solution, are close to the center of the biplot and are
%                  relatively unimportant to the interpretation. This
%                  biplot is often referred as contribution biplot because
%                  it shows visually the most contributing points
%                  (Greenacre 2006b).
%               If $plots.alpha=scalar$ in the interval [0 1], row
%                   coordinates are given by $D_r^{-1/2} U \Gamma^\alpha$
%                   and column coordinates are given by $D_c^{-1/2} V
%                   \Gamma^{1-\alpha}$. Note that for any choice of $\alpha$
%                   the matrix product $ D_r^{-1/2} U \Gamma^\alpha (D_c^{-1/2} V
%                   \Gamma^{1-\alpha})^T$ optimally approximates matrix
%                   $D_r^{-0.5}(P-rc')D_c^{-0.5}$, in the sense that the
%                   sum of squared differences between $D_r^{1/2}
%                   D_r^{-1/2} U \Gamma^\alpha (D_c^{-1/2} V
%                   \Gamma^{1-\alpha})^T D_c^{1/2}$ and
%                   $D_r^{-0.5}(P-rc')D_c^{-0.5}$ is as small as possible.
%              plots.dim = vector with two elements which specifies which
%                   dimensions to show in the factor map. The default is to
%                   show to first two dimensione therefore plots.dim=[1 2]
%              plots.FontSize = scalar which specifies the font size or row
%                   (column) labels. The default value is 10.
%              plots.MarkerSize =  scalar which specifies the marker size
%                   or symbols associated with rows or columns. The default
%                   value is 10.
%                 Example - 'plots',1
%                 Data Types - scalar double or struct
%  dispresults :  display results on the screen. Boolean.
%                 If dispresults is true (default) it is possible to see on the
%                 screen all the summary results of the analysis.
%                 Example - 'dispresults',false
%                 Data Types - Boolean
%        d1    :  dimension to show on the horizontal axis. Positive
%                 integer. Positive integer in the range 1, 2, .., K which
%                 indicates the dimension to show on the x axis. The
%                 default value of d1 is 1.
%                 Example - 'd1',2
%                 Data Types - single | double
%        d2    :  dimension to show on the vertical axis. Positive
%                 integer. Positive integer in the range 1, 2, .., K which
%                 indicates the dimension to show on the y axis. The
%                 default value of d2 is 2.
%                 Example - 'd2',3
%                 Data Types - single | double
%
%  Output:
%
%  out :     A structure containing the following fields
%
%
% 		out.N         =   $I$-by-$J$-table containing contigency table
%                         referred to active rows (i.e. referred to the rows which
%                         participated to the fit).
%                         The $(i,j)$-th element is equal to $n_{ij}$,
%                         $i=1, 2, \ldots, I$ and $j=1, 2, \ldots, J$. The
%                         sum of the elements of out.P is $n$ (the grand
%                         total).
% 		out.P         =   $I$-by-$J$-table containing correspondence matrix
%                         (proportions). The $(i,j)$-th element is equal to
%                         $n_{ij}/n$, $i=1, 2, \ldots, I$ and $j=1, 2,
%                         \ldots, J$.  The sum of the elements of out.P is
%                         1.
% 		out.I         =   Number of active rows of contingency table.
% 		out.J         =   Number of active columns of contingency table.
% 		out.n         =   Grand total. out.n is equal to sum(sum(out.N)).
%                         This is the number of observations.
% 		out.r         =   vector of length $I$ containing row masses.
%                         \[
%                           r=(f_{1.},	f_{2.}, \ldots, f_{I.})'
%                         \]
%                         $r$ is also the centroid of column profiles.
% 		out.Dr        =   square matrix of size $I$ containing on the
%                         diagonal the row masses.  This is matrix $D_r$.
%                         \[
%                           D_r=diag(r)
%                         \]
% 		out.c         =   vector of length $J$ containing columns masses.
%                         \[
%                           c=(f_{.1},	f_{.2}, \ldots, f_{.J})'
%                         \]
%                         $c$ is also the centroid of row profiles.
% 		out.Dc        =   square matrix of size $J$ containing on the
%                         diagonal the column masses. This is matrix $D_c$.
%                         \[
%                           D_c=diag(c)
%                         \]
% 		out.ProfilesRows =  $I$-by-$J$-matrix containing row profiles.
%                           The $i,j$-th element of this matrix is given by
%                           $f_{ij}/f_{i.}=n_{ij}/n_{i.}$.
%                           Written in matrix form:
%                           \[
%                           ProfilesRows = D_r^{-1} \times P
%                           \]
% 		out.ProfilesCols =  $I$-by-$J$-matrix containing column profiles.
%                           The $i,j$-th element of this matrix is given by
%                           $f_{ij}/f_{.j}=n_{ij}/n_{.j}$.
%                           Written in matrix form:
%                           \[
%                           ProfilesCols = P \times D_c^{-1}
%                           \]
% 		out.K         =   scalar integer containing maximun number of
%                         dimensions. $K = \min(I-1,J-1)$.
% 		out.k         =   scalar integer containing number of retained
%                         dimensions.
% 		out.Residuals =   $I$-by-$J$-matrix containing standardized residuals.
%                         \[
%                           Residuals = D_r^{1/2}  (D_r^{-1}  P - r  c')  D_c^{-1/2} =
%                           D_r^{-1/2}  (P - r  c')  D_c^{-1/2}
%                         \]
%                         With the singular value decomposition (SVD) we
%                         obtain that:
%                          \[
%                           Residuals = U \Gamma V'
%                          \]
% 		out.TotalInertia =   scalar containing total inertia. Total inertia
%                           is equal (for example) to the sum of the
%                           squares of the elements of matrix
%                           out.Residuals.
% 		out.Chi2stat  =   scalar containing Chi-square statistic for the
%                           contingency table. $Chi2stat= TotalInertia \times n$.
% 		out.CramerV   =    scalar containing Cramer's $V$ index.
%                          \[
%                           V=\sqrt{Chi2stat/(n (\min(I,J)-1))}
%                          \]
%                          Cramer's index goes between 0 and 1.
% 		out.RowsPri   =    $I$-by-$K$ matrix containing principal coordinates
%                          of rows.
%                           \[
%                           RowsPri = D_r^{-1/2} \times U \times \Gamma;
%                           \]
% 		out.ColsPri   =    $J$-by-$K$ matrix containing Principal coordinates
%                          of columns.
%                           \[
%                           ColsPri = D_c^{-1/2} \times V \times \Gamma;
%                           \]
% 		out.RowsSta   =   $I$-by-$K$ matrix containing standard coordinates
%                          of rows.
%                           \[
%                           RowsSta = RowsPri \times \Gamma^{-1} = D_r^{-1/2} U
%                           \Gamma \Gamma^{-1}= D_r^{-1/2}  U
%                           \]
% 		out.ColsSta   =   $J$-by-$K$ matrix containing standard coordinates
%                          of columns.
%                           \[
%                           ColsSta = ColsPri \times \Gamma^{-1} = D_c^{-1/2} V
%                           \Gamma \Gamma^{-1}= D_c^{-1/2}  V
%                           \]
% 		out.RowsSym   =    $I$-by-$K$ matrix containing symmetrical coordinates
%                          of rows.
%                           \[
%                           RowsSym = D_r^{-1/2} \times U \times \Gamma^{1/2}
%                           \]
% 		out.ColsSym   =     $J$-by-$K$ matrix containing symmetrical coordinates
%                          of columns.
%                           \[
%                           ColsSym = D_c^{-1/2} \times V \times \Gamma^{1/2}
%                           \]
%                       Symetric plot represents the row and column
%                       profiles simultaneously in a common space
%                       (Bendixen, 2003). In this case, only the distance
%                       between row points or the distance between column
%                       points can be really interpreted.
%                       The distance between any row and column items is
%                       not meaningful! You can only make a general
%                       statements about the observed pattern. In order to
%                       interpret the distance between column and row
%                       points, the column profiles must be presented in
%                       row space or vice-versa. This type of map is called
%                       asymmetric biplot.
% out.InertiaRows     =    $I$-by-$2$ matrix containing absolute and relative
%                          contribution of each row to total inertia.
%                          The inertia of a point is the squared distance
%                          of point $d_i^2$ to the centroid. The absolute
%                          contribution of a point to total inertia is the
%                          inertia of the point multiplied by the point
%                          mass.
%                          1st column = absolute contribution of each row
%                          to TotalInertia. The sum of values of the first
%                          column is equal to TotalInertia;
%                          2nd column = relative contribution of each row
%                          to TotalInertia. The sum of the values of the
%                          second column is equal to 1.
% out.InertiaCols     =    $J$-by-$2$ matrix containing absolute and relative
%                          contribution of each column to total inertia.
%                          The inertia of a point is the squared distance
%                          of point $d_i^2$ to the centroid. The absolute
%                          contribution of a point to total inertia is the
%                          inertia of the point multiplied by the point
%                          mass.
%                          1st column = absolute contribution of each
%                          column to TotalInertia. The sum of values of the
%                          first column is equal to TotalInertia;
%                          2nd column = relative contribution of each
%                          column to TotalInertia. The sum of values of the
%                          second column is equal to 1.
%out.Point2InertiaRows  =  $I$-by-$K$ matrix  containing relative
%                          contributions of rows to inertia of the
%                          dimension. The inertia of first latent dimension
%                          is given by $\lambda_1=\gamma_{11}^2$. The
%                          inertia of second latent dimension is given by
%                          $\lambda_2=\gamma_{22}^2$ .... The sum of each
%                          column of matrix Point2InertiaRows is equal to
%                          1.
%                           Remark: the points  with the larger value of
%                           Point2Inertia are those which contribute the
%                           most to the definition of the dimension. If the
%                           row contributions were uniform, the expected
%                           value would be 1/size(contingeny_table,1) For a
%                           given dimension, any row with a contribution
%                           larger than this threshold could be considered
%                           as important in contributing to that dimension.
%out.Point2InertiaCols  =  $J$-by-$K$ matrix  Relative contributions of
%                          columns to inertia of the dimension. The sum of
%                          each column of matrix Point2InertiaCols is equal
%                          to 1.
% out.Dim2InertiaRows   =  $I$-by-$K$ matrix  containing relative
%                          contributions of latent dimensions to inertia of
%                          the row points.  These numbers can be
%                          interpreted as squared correlations and measures
%                          the degree of association between row points
%                          and a particular axis. The sum of
%                          each row of matrix Dim2InertiaRows is equal to
%                          1.
% out.Dim2InertiaCols   =  $J$-by-$K$ matrix  containing relative
%                          contributions of latent dimensions to inertia of
%                          the column points.  These numbers can be
%                          interpreted as squared correlations and measure
%                          the degree of association between columns points
%                          and a particular axis. The sum of each row of
%                          matrix Dim2InertiaCols is equal to 1.
% out.cumsumDim2InertiaRows   =  $I$-by-$K$ matrix  containing cumulative
%                          sum of the contributions of latent dimensions to
%                          inertia of the row points. These cumulative sums
%                          are equivalent to the communalities in PCA.
%                          The last column of matrix cumsumDim2InertiaRows
%                          is equal to 1.
% out.cumsumDim2InertiaCols   =  $J$-by-$K$ matrix  containing cumulative
%                          sum of the contributions of latent dimensions to
%                          inertia of the column points. These cumulative sums
%                          are equivalent to the communalities in PCA.
%                          The last column of matrix cumsumDim2InertiaCols
%                          is equal to 1.
% out.sqrtDim2InertiaRows = $I$-by-$K$ matrix  containing correlation of
%                           rows points with latent dimension axes. Similar
%                           to component loadings in PCA
% out.sqrtDim2InertiaCols = $I$-by-$K$ matrix  containing correlation of
%                           column points with latent dimension axes. Similar
%                           to component loadings in PCA.
% 		out.SupRowsN   =   matlab table containing contingency table referred
%                         to supplementary rows. If there are no
%                         supplementary rows this field is empty.
% 		out.SupColsN  =   matlab table containing contingency table related
%                         to supplementary columns. If there are no
%                         supplementary columns this field is empty.
% 		out.RowsPriSup    = Principal coordinates of supplementary rows
% 		out.RowsStaSup    = Standard coordinates of supplementary rows
%       out.RowsSymSup    = Symmetrical coordinates of supplementary
%                           rows
% 		out.ColsPriSup    = Principal coordinates of supplementary columns
% 		out.ColsStaSup    = Standard coordinates of of supplementary columns
%       out.ColsSymSup    = Symmetrical coordinates of supplementary
%                                columns
%       out.Summary       = $K$-times-4 table containing summary results
%                           for correpondence analysis.
%                           First column contains the singular values (the
%                           sum of the squared singular values is the total
%                           inertia).
%                           Second column contains the eigenvalues  (the
%                           sum of the eigenvalues is the total inertia).
%                           Third column contains the variance explained by
%                           each latent dimension. Fourth column contains
%                           the cumulative variance explained by each
%                           dimension.
%       out.OverviewRows = $I$-times-(k*3+2) table containing an overview
%                          of row points. More precisely if we suppose that $k=2$
%                          First column contains the row masses (vector
%                          $r$).
%                          Second column contains the scores of first dimension.
%                          Third column contains the scores of second dimension.
%                          Fourth column contains the inertia of each
%                          point, where inertia of point is the squared
%                          distance of point $d_i^2$ to the centroid.
%                          Fifth column contains the relative contribution
%                          of each point to the explanation of the inertia
%                          of the first dimension. The sum of the elements
%                          of this column is equal to 1.
%                          Sixth column contains the relative contribution
%                          of each point to the explanation of the inertia
%                          of the second dimension. The sum of the elements
%                          of this column is equal to 1.
%                          Seventh column contains the relative
%                          contribution of the first dimension to the
%                          explanation of the inertia of the point.
%                          Eight column contains the relative
%                          contribution of the second dimension to the
%                          explanation of the inertia of the point.
%       out.OverviewCols = $J$-times-(k*3+2) table containing an overview
%                          of row points. More precisely if we suppose that $k=2$
%                          First column contains the column masses (vector
%                          $c$).
%                          Second column contains the scores of first dimension.
%                          Third column contains the scores of second dimension.
%                          Fourth column contains the inertia of each
%                          point, where inertia of point is the squared
%                          distance of point $d_i^2$ to the centroid.
%                          Fifth column contains the relative contribution
%                          of each point to the explanation of the inertia
%                          of the first dimension. The sum of the elements
%                          of this column is equal to 1.
%                          Sixth column contains the relative contribution
%                          of each point to the explanation of the inertia
%                          of the second dimension. The sum of the elements
%                          of this column is equal to 1.
%                          Seventh column contains the relative
%                          contribution of the first dimension to the
%                          explanation of the inertia of the point.
%                          Eight column contains the relative
%                          contribution of the second dimension to the
%                          explanation of the inertia of the point.
%
% See also crosstab, rcontFS, CressieRead
%
% References:
%
% Benzecri, J.-P. (1992), Correspondence Analysis Handbook, New-York :
% Dekker
% Benzecri, J.-P. (1980), L’analyse des données tome 2 : l’analyse des
% correspondances, Paris : Bordas.
% Greenacre, M.J. (1993), Correspondence Analysis in Practice, London :
% Academic Press.
% Gabriel, K.R. and Odoroff, C. (1990). Biplots in biomedical research.
% Statistics in Medicine, 9, pp. 469-485.
% Greenacre, M.J. (1993) Biplots in correspondence Analysis, Journal of
% Applied Statistics, 20, pp. 251 - 269.
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('CorAna')">Link to the help function</a>
% Last modified 14-06-2016
%
% Last modified 31-05-2016
%
% Examples:

%{
    % CorAna with all the default options.
     load smoke

%? VERIFY MODIFICATION ?
    %N=crosstab(X(:,1),X(:,2));   
    N=crosstab(smoke.data(:,1),smoke.data(:,2));

    out=CorAna(N);
%}

%{
    %% CorAna with name pairs.
    % Input is the contingency table, labels for rows and columns are supplied
    load smoke

%? VERIFY MODIFICATION ?
    %N=crosstab(X(:,1),X(:,2));
    N=crosstab(smoke.data(:,1),smoke.data(:,2));
    %out=CorAna(N,'Lr',labels_rows,'Lc',labels_columns);
    labels_rows = unique(smoke.data(:,1))
    labels_columns = unique(smoke.data(:,2))
    out=CorAna(N,'Lr',labels_rows,'Lc',labels_columns);
%}

%% Beginning of code

% Check whether N is a contingency table or a n-by-p input dataset (in this
% last case the contigency table is built using the first tow columns of the
% input dataset).
if ~isempty(varargin)
    UserOptions=varargin(1:2:length(varargin));
    checkdatamatrix = strcmp(UserOptions,'datamatrix');
    if sum(checkdatamatrix)
        datamatrix = varargin{2*find(checkdatamatrix)};
    else
        datamatrix=false;
    end
else
    datamatrix=false;
end

% If input is a datamatrix it is necessary to construct the contingency
% table
if datamatrix == true
    [N,~,~,labels] =crosstab(N(:,1),N(:,2));
    [I,J]=size(N);
    % default labels for rows of contingency table
    Lr=labels(1:I,1);
    % default labels for columns of contingency table
    Lc=labels(1:J,2);
else
    [I,J]=size(N);
    % Size of N
    % default labels for rows of contingency table
    Lr=cellstr(strcat('r',num2str((1:I)')));
    % default labels for columns of contingency table
    Lc=cellstr(strcat('c',num2str((1:J)')));
end

% default value for supplementary units
Sup='';
k=2;
plots=1;

% Default font size for labels of rows or colums to add to the plot
FontSizedef=10;
MarkerSizedef=10;
dispresults=true;

% Dimensions to show in the plot. The default is to show the first two
% dimensions.
% d1= dimension to show in the x axis of correspondence analysis plot
d1=1;
% d2= dimension to show in the y axis of correspondence analysis plot
d2=2;

options=struct('k',k,'Sup',Sup,'plots',plots,'datamatrix',false,...
    'dispresults',dispresults,'d1',d1,'d2',d2);
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
    k   = options.k;
    Lr  = options.Lr;
    Lc  = options.Lc;
    Sup = options.Sup;
    plots=options.plots;
    dispresults=options.dispresults;
    d1  = options.d1;
    d2  = options.d2;
end

% Extract labels for rows and columns
if istable(N)
    Lc=N.Properties.VariableNames;
    Lr=N.Properties.RowNames;
    Ntable=N;
    N=table2array(N);
else
    if isempty(Lr)
        Lr=cellstr(num2str((1:I)'));
    else
        % Check that the length of Lr is equal to I
        if length(Lr)~=I
            error('Wrong length of row labels');
        end
    end
    
    if isempty(Lc)
        Lc=cellstr(num2str((1:J)'));
    else
        % Check that the length of Lc is equal to J
        if length(Lc)~=J
            error('Wrong length of column labels');
        end
    end
    Ntable=array2table(N,'RowNames',Lr,'VariableNames',Lc);
end

% Nred will contain the contingency table after removing supplementary rows
% and columns (if supplementary rows and columns belong to the table)
Nred      = N;
Nredtable = Ntable;

if ~isempty(Sup)
    labels=struct;
    
    % if Sup.r (Sup.c) is a cell or is character or it is a numeric vector,
    % then the supplementary rows (columns) belong to the the actual
    % contingency table, else if Sup.r (Sup.c) is a Matlab table then the
    % supplementary units do not belong to the actual contingency table N
    
    if isfield(Sup,'r')
        
        if iscell(Sup.r) || ischar(Sup.r)
            % find the indexes of the rows of matrix N to delete (rows to
            % use as supplementary rows)
            Indexesr=zeros(length(Sup.r),1);
            for i=1:length(Sup.r)
                if iscell(Sup.r)
                    Indexesr(i)=find(strcmp(Sup.r{i},Lr),1);
                else
                    Indexesr(i)=find(strcmp(Sup.r(i),Lr),1);
                end
            end
            
            labels.sr = Lr;
            labels.sr = labels.sr(Indexesr);
            % Delete the labels of contingency table associated to
            % supplementary rows
            Lr(Indexesr) = [];
            
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
        % out.labels.sr = labels.sr;
    else
        Indexesr='';
    end
    % end of part referred to labels for supplementary rows
    
    % beginning of part referred to labels for supplementary columns
    if isfield(Sup,'c')
        if iscell(Sup.c) || ischar(Sup.c)
            % find the indexes of the rows to delete (rows to use as
            % supplementary rows)
            Indexesc=zeros(length(Sup.c),1);
            for i=1:length(Sup.c)
                if iscell(Sup.c)
                    Indexesc(i)=find(strcmp(Sup.c{i},Lc),1);
                else
                    Indexesc(i)=find(strcmp(Sup.c(i),Lc),1);
                end
            end
            labels.sc=Lc;
            labels.sc = labels.sc(Indexesc); %#ok<STRNU>
            % Delete the labels of contingency table associated to
            % supplementary rows
            Lc(Indexesc)=[];
            
        elseif istable(Sup.c)
            Indexesc='';
            Nsupc=table2array(Sup.c);
            labels.sc=Sup.c.Properties.VariableNames; %#ok<STRNU>
        else
            
            Indexesc=Sup.c;
            if min(Indexesc)<1 || max(Indexesc)> size(N,2)
                error('FSDA:CorAna:wrongInput',['Numeric indexes of supplementary columns must be integers between 1 and ' num2str(size(N,2))])
            end
            labels.sc=Lc;
            labels.sc = labels.sc(Indexesc); %#ok<STRNU>
            % Delete the labels of contingency table associated to
            % supplementary rows
            Lc(Indexesc)=[];
            
        end
        %  out.labels.sc = labels.sc;
    else
        Indexesc='';
    end
    
    % if ~isempty(Indexesr) this means that supplementary rows belong to
    % matrix N
    if ~isempty(Indexesr)
        % Contingency table referred to supplementary rows.
        if ~isempty(Indexesc)
            Nsupr=N(Indexesr,Indexesc);
            Nsuprtable=Ntable(Indexesr,Indexesc);
        else
            Nsupr=N(Indexesr,:);
            Nsuprtable=Ntable(Indexesr,:);
        end
    else
        Nsuprtable='';
    end
    
    % if ~isempty(Indexesc) this means that supplementary columns belong to
    % matrix N
    if ~isempty(Indexesc)
        % Contingency table referred to supplementary columns.
        if ~isempty(Indexesr)
            Nsupc=N(Indexesr,Indexesc);
            Nsupctable=Ntable(Indexesr,Indexesc);
        else
            Nsupc=N(:,Indexesc);
            Nsupctable=Ntable(:,Indexesc);
        end
    else
        Nsupctable='';
    end
    
    % Delete the rows and columns of contingency table associated to
    % supplementary rows and/or associated to supplementary columns
    if ~isempty(Indexesr)
        Nred(Indexesr,:)=[];
        Nredtable(Indexesr,:)=[];
    end
    
    if ~isempty(Indexesc)
        Nred(:,Indexesc)=[];
        Nredtable(:,Indexesc)=[];
    end
end

% Store contingency table (in Matlab table format)
out.N = Nredtable;

[I,J]=size(Nred);

%vectors of ones
onesI = ones(I,1);
onesJ = ones(J,1);

%grand total
n=sum(sum(Nred));

% P = correspondence matrix  containing relative frequencies
P = (1/n) * Nred;

% Store P in table format
Ptable=array2table(P,'RowNames',Lr,'VariableNames',Lc);
out.P=Ptable;

out.I=I;        %number of active rows (excluding supplementary rows)
out.J=J;        %number of active columns (excluding supplementary columns)
out.n=n;        %grand total

% out.Lr=Lr;
% out.Lc=Lc;

% r= vector which contains row masses = centroids of the column profiles
r  = P * onesJ ;
Dr = diag(r);
out.r = r;          %row masses (vector)
out.Dr = Dr;        %row masses (diagonal matrix)

% Column masses = centroids of the row profiles r' * Dr^(-1) * P = 1' * P = c'
c  = (onesI' * P)';
Dc = diag(c);
out.c = c;          %column masses (vector); c is the centroid of row profiles
out.Dc = Dc;        %column masses (diagonal matrix)

%Rows profiles (equation 4.14)
ProfilesRows = Dr^(-1) * P;
out.ProfilesRows = ProfilesRows;

%Columns profiles  (equation 4.14)
ProfilesCols = P * Dc^(-1);
out.ProfilesCols = ProfilesCols;

% K = maximun number of dimensions
K = min(I-1,J-1);
out.K = K;

% k= number of retained dimensions
out.k = k;

% Standarized residuals
%
% Residuals   = Dr^(1/2) * (Dr^(-1) * P - Ir * c') * Dc^(-1/2) = Dr^(-1/2) * (P - r * c') * Dc^(-1/2);
% Residualsij = sqrt( (p_{ij} - r_ic_j)^2 /r_ic_j ) =(p_{ij} - r_ic_j)/sqrt(r_ic_j)
% Chi-square distances
Residuals     =  Dr^(-1/2) * (P - r * c') * Dc^(-1/2);
out.Residuals = Residuals; 

% SVD of Residuals = U*Gam*V'
[U,Gam,V] = svd(Residuals);
Gam = Gam(1:K,1:K);
U   = U(:,1:K);
V   = V(:,1:K);

% Total inertia
% TotalInertia = sum_i sum_j (pij - ricj)^2 / ricj = chis/n
TotalInertia     = sum(sum(Residuals.^2));
out.TotalInertia = TotalInertia;

% Chi-square statistic for the contingency table
%
% chi2     = n * sum_i sum_j (pij - ricj)^2 / ricj
% Chi2stat = onesI'* (Residuals.*Residuals) * onesJ * n;
Chi2stat     = n*TotalInertia;
out.Chi2stat = Chi2stat;

% Cramer's V
out.CramerV = sqrt(Chi2stat/(n*(min(I,J)-1)));

% Principal inertias
% eigenvalues of Residuals'*Residuals
% percentages of the total inertia

% cumsumTotalInertia = cumulative proportion of explained inertia
Gam2 = Gam.^2;
cumsumTotalInertia = cumsum(diag(Gam2))/TotalInertia;

% InertiaExplained is a matrix with 4 columnn.
% - First column contains the singular values (the sum of the squared
%   singular values is the total inertia)
% - Second column contains the eigenvalues (the sum of the eigenvalues is 
%   the total inertia)
% - Third column contains the variance explained by each latent dimension.
% - Fourth column contains the cumulative variance explained by each
%   dimension.
InertiaExplained=[diag(Gam) diag(Gam2) diag(Gam2 / TotalInertia) cumsumTotalInertia];

% Principal coordinates of rows  (alpha=1 for the rows) F=...
RowsPri     = Dr^(-1/2) * U*Gam;
out.RowsPri = RowsPri;

% Principal coordinates of columns G= (Dc^(-1/2)*V*Gamma)
ColsPri     = Dc^(-1/2) * V*Gam;
out.ColsPri = ColsPri;

% Standard coordinates of rows
% RowsSta (X) = RowsPri * Gam^(-1) = Dr^(-1/2) * U * Gam * Gam^(-1) = Dr^(-1/2) * U
RowsSta     = Dr^(-1/2) * U ;
out.RowsSta = RowsSta;

% Standard coordinates of columns
% ColsSta (Y) = ColsPri * Gam^(-1) = Dc^(-1/2) * V * Gam * Gam^(-1) = Dc^(-1/2) * V
ColsSta     = Dc^(-1/2) * V ;
out.ColsSta = ColsSta;

% Symmetrical coordinates of rows
RowsSym     = Dr^(-1/2) * U*Gam^(1/2);
out.RowsSym = RowsSym;

% Symmetrical coordinates of columns
ColsSym     = Dc^(-1/2) * V*Gam^(1/2);
out.ColsSym = ColsSym;

%Inertia of each row
InertiaRows = Dr*sum(RowsPri.*RowsPri,2);

% Relative inertia of each row
InertiaRows_relative=InertiaRows / TotalInertia;
out.InertiaRows = [InertiaRows InertiaRows_relative];

% Inertia of each column
% InertiaCols = Dc*(ColsPri.*ColsPri)*ones(K,1);
InertiaCols = Dc*sum(ColsPri.*ColsPri,2);

% Relative inertia of each column
InertiaCols_relative=InertiaCols / TotalInertia;
out.InertiaCols = [InertiaCols InertiaCols_relative];

% Contributions of rows to inertia of the dimension
% Absolute row contributions: contributions of points to axes.
Point2InertiaRows = Dr*RowsSta.*RowsSta;
out.Point2InertiaRows = Point2InertiaRows;

% Contributions of columns to inertia of the dimension
% Absolute column contributions: Contributions of points to axes.
Point2InertiaCols = Dc*ColsSta.*ColsSta;
out.Point2InertiaCols = Point2InertiaCols;

% Squared correlations of row points with axes
Dim2InertiaRows = diag(InertiaRows)^(-1)*Dr*(RowsPri.*RowsPri);
% Relative row contributions: contributions of axes (latent dimension) to
% points; squared correlations of rows with axes
out.Dim2InertiaRows=Dim2InertiaRows;

% Squared correlations of column points with axes. Relative columns
% contributions: contributions of axes to points; squared correlations of
% points with axes
Dim2InertiaCols = diag(InertiaCols)^(-1)*Dc*(ColsPri.*ColsPri);
out.Dim2InertiaCols=Dim2InertiaCols;

% Cumulative sum of quality of display for each row in the reduced space.
% The qualities are equivalent to the communalities in PCA
cumsumDim2InertiaRows = cumsum(Dim2InertiaRows,2);
out.cumsumDim2InertiaRows = cumsumDim2InertiaRows;

% Cumulative sum of quality of display for each colum  in the reduced space.
% The qualities are equivalent to the communalities in PCA
cumsumDim2InertiaCols = cumsum(Dim2InertiaCols,2);
out.cumsumDim2InertiaCols = cumsumDim2InertiaCols;

% Correlation of rows points with axes. Similar to component loadings in PCA
sqrtDim2InertiaRows = sign(RowsPri).*sqrt(Dim2InertiaRows);
out.sqrtDim2InertiaRows = sqrtDim2InertiaRows;

% Correlation of columns points with axes. Similar to component loadings in PCA
sqrtDim2InertiaCols = sign(ColsPri).*sqrt(Dim2InertiaCols);
out.sqrtDim2InertiaCols = sqrtDim2InertiaCols;

if exist('Sup','var')
    %Supplementary rows
    if isfield(Sup,'r')
        % Store table referred to supplementary rows
        out.SupRowsN = Nsuprtable;
        
        % The sum of each row of h must be equal to 1
        % h=Nsup(Indexesr,:)/(diag(sum(Nsup(Indexesr,:))));
        h=bsxfun(@rdivide,Nsupr,sum(Nsupr,2));
        
        RowsPriSup=h*ColsSta;             %Principal coordinates of supplementary rows
        RowsStaSup=h*ColsSta*Gam^(-1);    %Standard coordinates of supplementary rows
        RowsSymSup=h*ColsSta*Gam^(-1/2);  %Symmetrical coordinates of supplementary rows
        out.RowsPriSup=RowsPriSup;
        out.RowsStaSup=RowsStaSup;
        out.RowsSymSup=RowsSymSup;
        % rrc=(Gsup(:,1:k).*Gsup(:,1:k))/trace(Gsup(:,1:k)'*Gsup(:,1:k));
    end
    
    %Supplementary columns
    if isfield(Sup,'c')
        % Store table referred to supplementary columns
        out.SupColsN = Nsupctable;
        
        % The sum of each column of h must be equal to 1
        h=Nsupc/(diag(sum(Nsupc)));
        ColsPriSup=h'*RowsSta;                              %Principal coordinates of supplementary columns
        ColsStaSup=h'*RowsSta*Gam^(-1);                     %Standard coordinates of supplementary columns
        ColsSymSup=h'*RowsSta*Gam^(-1/2);                   %Symmetrical coordinates of supplementary columns
        out.ColsPriSup=ColsPriSup;
        out.ColsStaSup=ColsStaSup;
        out.ColsSymSup=ColsSymSup;
        %     rrc=(Gsup(:,1:k).*Gsup(:,1:k))/trace(Gsup(:,1:k)'*Gsup(:,1:k));
    end
end

if isstruct(plots) || plots==1
    FontName='Times';
    FontSizeAxisLabels=12;
    
    if isstruct(plots)
        % This anonymous function anables to extract the variable name to a
        % string
        ExtractVariableName=@(x) inputname(1);
        
        if isfield(plots,'alpha')
            if strcmp(plots.alpha,'rowprincipal')
                typeR='RowsPri'; % rows are in principal coordinates
                typeC='ColsSta'; % columns are in standard coordinates
                titl={'Rows principal coordinates, and column standard coordinates' , ...
                      '$\alpha=1$, $X=D_r^{-1/2}U\Gamma$ and $Y= D_c^{-1/2} V$'};
                
            elseif strcmp(plots.alpha,'colprincipal')
                typeR='RowsSta'; % rows are in standard coordinates
                typeC='ColsPri'; % columns are in principal coordinates
                titl={'Rows standard coordinates, and column principal coordinates' , ...
                      '$\alpha=0$, $X=D_r^{-1/2}U $ and $G= D_c^{-1/2} V \Gamma$'};
                
            elseif strcmp(plots.alpha,'symbiplot')
                % equivalent to alpha=0.5
                typeR='RowsSym';        % rows are in symmetrical coordinates
                typeC='RowsSym';        % columns are in symmetrical coordinates
                titl='Biplot symmetrical model $\alpha=0.5$ $X=D_r^{-1/2}U\Gamma^{1/2} $ and $Y= D_c^{-1/2} \Gamma V^{1/2}$';
                
            elseif strcmp(plots.alpha,'bothprincipal')
                typeR='RowsPri';        % rows are in principal coordinates
                typeC='ColsPri';        % columns are in principal coordinates
                titl={'French symmetrical model: rows and cols in principal coordinates.' , ...
                      'Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$'};
                
            elseif strcmp(plots.alpha,'rowgab')
                %  If plots.alpha='rowgab'  rows are in principal coordinates
                %  and columns are in standard coordinates multiplied by the
                %  mass.
                typeR='RowsPri';        % rows are in principal coordinates
                ColsStaDc=Dc*ColsSta;
                typeC=ExtractVariableName(ColsStaDc);
                titl={'Rows principal coordinates, and column standard coordinates times masses' , ...
                      '$\alpha=1$,$X=D_r^{-1/2}U\Gamma $ and $Y= D_c^{1/2} V$'};
                      
            elseif strcmp(plots.alpha,'colgab')
                % If plots.alpha='colgab'  columns are in principal coordinates
                % and rows are in standard coordinates multiplied by the
                % masses.
                RowsStaDr=Dr*RowsSta;
                typeR=ExtractVariableName(RowsStaDr);
                typeC='ColsPri';        % columns are in principal coordinates
                titl={'Rows standard coordinates multiplied by masses ' , ...
                      'and column principal coordinates $X=D_r^{-1/2} U$ and $Y= D_c^{-1/2} V \Gamma$'};
                
            elseif strcmp(plots.alpha,'rowgreen')
                %  If plots.alpha='rowgreen'  rows are in principal
                %  coordinates and columns are in standard coordinates
                %  multiplied by square root of the mass.
                typeR='RowsPri';        % rows are in principal coordinates
                ColsStaDcSqrt=(Dc^(1/2))*ColsSta;
                typeC= ExtractVariableName(ColsStaDcSqrt);
                titl={'Rows principal coordinates, and column standard coordinates ' , ...
                      'times sqrt of masses $X=D_r^{-1/2}U\Gamma $ and $Y= V$'};
                
            elseif strcmp(plots.alpha,'colgreen')
                %  If plots.alpha='colgreen' columns in principal coordinates
                %  and rows in standard coordinates multiplied by the square
                %  root of the mass.
                RowsStaDrSqrt=(sqrt(Dr))*RowsSta;
                typeR=ExtractVariableName(RowsStaDrSqrt);
                typeC='ColsPri';        % columns are in principal coordinates
                titl={'Rows standard coordinates times sqrt of masses,' ...
                      'and column principal coordinates, $X=U $ and $G= D_c^{-1/2} V \Gamma$'};
                
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
        else
            typeR='RowsPri';        % rows are in principal coordinates
            typeC='ColsPri';        % columns are in principal coordinates
            titl={'French symmetrical model: rows and cols in principal coordinates.' ...
                  'Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$'};
        end
        
        if isfield(plots,'FontSize')
            FontSize=plots.FontSize;
        else
            FontSize=FontSizedef;
        end
        if isfield(plots,'MarkerSize')
            MarkerSize=plots.MarkerSize;
        else
            MarkerSize=MarkerSizedef;
        end
        
    else
        typeR='RowsPri';        % rows are in principal coordinates
        typeC='ColsPri';        % columns are in principal coordinates
        titl={'French symmetrical model: rows and cols in principal coordinates.'...
              'Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$'};
        FontSize=FontSizedef;
        MarkerSize=MarkerSizedef;
        
    end
    symbolrows='o';
    symbolcols='^';
    symbolsuprows='o';
    symbolsupcols='^';
    % Color for symbols and text for rows points
    colorrows='b';
    % Color for symbols and text for  column rows
    colorcols='r';
    % Color for symbols and text for supplementary row points
    colorsuprows='b';
    % Color for symbols and text for supplementary column points
    colorsupcols='r';
    
    d1str=num2str(d1);
    d2str=num2str(d2);
    
    figure
    hold('on')
    % Plot row points
    propR=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolrows ,''',''Color'',''', colorrows , ...
        ''',''MarkerSize'',', num2str(MarkerSize)   ,'');
    
    eval(['plot(' typeR '(:,' d1str '),' typeR '(:,' d2str '),' propR ')'])
    
    % Plot column points
    propC=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolcols ,''',''Color'',''', colorcols, ...
        ''',''MarkerSize'',', num2str(MarkerSize)   ,'');
    
    eval(['plot(' typeC '(:,' d1str '),' typeC '(:,' d2str '),' propC ')'])
    
    % Add labels for row points and column points
    % addx = adds a small right horizontal displacement for labels
    addx=0.04;
    eval(['text(' typeR '(:,' d1str ')+' num2str(addx) ',' typeR '(:,' d2str '),Lr,''Interpreter'',''None'',''FontSize'',' num2str(FontSize) ',''Color'',''' colorrows ''')'])
    eval(['text(' typeC '(:,' d1str ')+' num2str(addx) ',' typeC '(:,' d2str '),Lc,''Interpreter'',''None'',''FontSize'',' num2str(FontSize) ',''Color'',''' colorcols ''')'])
    
    title(titl,'Interpreter','Latex');
    
    % Labels for axes
    xlabel(['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',InertiaExplained(d1,3)*100),'%)'],'FontName', FontName, 'FontSize', FontSizeAxisLabels);
    ylabel(['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',InertiaExplained(d2,3)*100),'%)'],'FontName', FontName, 'FontSize', FontSizeAxisLabels);
    
    % Add points and text associated to supplementary rows
    if isstruct(Sup) && isfield(Sup,'r')
        propsupR=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolsuprows ,''',''Color'',''', colorsuprows , ''',''MarkerFaceColor'',''', colorsuprows ,...
            ''',''MarkerSize'',', num2str(MarkerSize)   ,'');
        
        eval(['plot(' typeR 'Sup(:,d1),' typeR 'Sup(:,d2),' propsupR ')'])
        %         eval(['text(' typeR 'Sup(:,d1)+' num2str(addx) ',' typeR 'Sup(:,d2),''' labels.sr{:} ''',''Interpreter'',''None'',''FontSize'',' num2str(FontSize) ',''Color'',''' colorrows ''')'])
        eval(['text(' typeR 'Sup(:,d1)+' num2str(addx) ',' typeR 'Sup(:,d2),labels.sr,''Interpreter'',''None'',''FontSize'',' num2str(FontSize) ',''Color'',''' colorrows ''')'])
        
    end
    
    % Add points and text associated to supplementary columns
    if isstruct(Sup) && isfield(Sup,'c')
        propsupC=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolsupcols ,''',''Color'',''', colorsupcols , ''',''MarkerFaceColor'',''', colorsupcols ,...
            ''',''MarkerSize'',', num2str(MarkerSize)   ,'');
        
        eval(['plot(' typeC 'Sup(:,d1),' typeC 'Sup(:,d2),' propsupC ')'])
        % eval(['text(' typeC 'Sup(:,d1)+' num2str(addx) ',' typeC 'Sup(:,d2),''' labels.sc{:} ''',''Interpreter'',''None'',''FontSize'',' num2str(FontSize) ',''Color'',''' colorcols ''')'])
        eval(['text(' typeC 'Sup(:,d1)+' num2str(addx) ',' typeC 'Sup(:,d2),labels.sc,''Interpreter'',''None'',''FontSize'',' num2str(FontSize) ',''Color'',''' colorcols ''')'])
    end
    
    % Make axis equal and add cartesian axes
    axis(gca,'equal')
    axis(gca,'equal')
    vv=axis;
    line([vv(1);vv(2)],[0;0])
    line([0;0],[vv(3);vv(4)])
    
end
% Score Rows and ScoreCols respectively contain scores for rows and columns
ScoreRows=eval(typeR);
ScoreCols=eval(typeC);


% Explained is a matrix with 4 columnn.
% First column contains the singular values (the sum of the squared singular values is the
% total inertia)
% Second column contains the eigenvalues  (the sum of the eigenvalues is the
% total inertia)
% Third column contains the variance explained by each latent dimension.
% Fourth column contains the cumulative variance explained by each
% dimension.
ColNamesSummary={'Singular_value' 'Inertia' 'Accounted_for' 'Cumulative'};
RowNamesSummary=strcat(cellstr(repmat('dim_',K,1)), cellstr(num2str((1:K)')));
RowNamesSummary=regexprep(RowNamesSummary,' ','');

InertiaExplainedtable=array2table(InertiaExplained,'VariableNames',ColNamesSummary, ...);
    'RowNames',RowNamesSummary);

out.Summary = InertiaExplainedtable ;

ColNames={'Mass' 'Score_1' 'Score_2' 'Inertia' ,...
    'CntrbPnt2In_1' 'CntrbPnt2In_2' ...
    'CntrbDim2In_1' 'CntrbDim2In_2'};
OverviewRows=array2table([out.r ScoreRows(:,1:k) InertiaRows Point2InertiaRows(:,1:k) Dim2InertiaRows(:,1:k)],...
    'VariableNames',ColNames,'RowNames',Lr);
out.OverviewRows=OverviewRows;

OverviewCols=array2table([out.c ScoreCols(:,1:k) InertiaCols Point2InertiaCols(:,1:k) Dim2InertiaCols(:,1:k)],...
    'VariableNames',ColNames,'RowNames',Lc);
out.OverviewCols=OverviewCols;

if dispresults==true
    disp('Summary')
    disp(out.Summary)
    
    VarNamesforTab={'Scores', 'CntrbPnt2In' 'CntrbDim2In'};
    % CntrbPnt2In = relative contribution of points to explain total Inertia
    % CntrbDim2In = relative contribution of latent dimension to exaplin total
    %               Inertia
    disp('ROW POINTS')
    disp(['Results for dimension: ' d1str])
    Tabresults=array2table(eval(strcat('[',typeR,'(:,', d1str ,') Point2InertiaRows(:,', d1str ,')    Dim2InertiaRows(:,', d1str ,')         ]')));
    Tabresults.Properties.RowNames=Lr;
    Tabresults.Properties.VariableNames=VarNamesforTab;
    disp(Tabresults)
    
    disp(['Results for dimension: ' d2str])
    Tabresults=array2table(eval(strcat('[',typeR,'(:,', d2str ,') Point2InertiaRows(:,', d2str ,')    Dim2InertiaRows(:,', d2str ,')         ]')));
    Tabresults.Properties.RowNames=Lr;
    Tabresults.Properties.VariableNames=VarNamesforTab;
    disp(Tabresults)
    
    disp('COLUMN POINTS')
    disp(['Results for dimension: ' d1str])
    Tabresults=array2table(eval(strcat('[',typeC,'(:,', d1str ,') Point2InertiaCols(:,', d1str ,')    Dim2InertiaCols(:,', d1str ,')         ]')));
    Tabresults.Properties.RowNames=Lc;
    Tabresults.Properties.VariableNames=VarNamesforTab;
    disp(Tabresults)
    
    disp(['Results for dimension: ' d2str])
    Tabresults=array2table(eval(strcat('[',typeC,'(:,', d2str ,') Point2InertiaCols(:,', d2str ,')    Dim2InertiaCols(:,', d2str ,')         ]')));
    Tabresults.Properties.RowNames=Lc;
    Tabresults.Properties.VariableNames=VarNamesforTab;
    disp(Tabresults)
    disp('-----------------------------------------------------------')
    disp('Overview ROW POINTS')
    disp(OverviewRows)
    disp('Overview COLUMN POINTS')
    disp(OverviewCols)
    disp('-----------------------------------------------------------')
    disp('Legend')
    disp('CntrbPnt2In = relative contribution of points to explain total Inertia of the latent dimension')
    disp('              The sum of the numbers in a column is equal to 1')
    disp('CntrbDim2In = relative contribution of latent dimension to explain total Inertia of a point')
    disp('              CntrbDim2In_1+CntrbDim2In_2+...+CntrbDim2In_K=1')
    
    %TODO
    %     disp('SUPPLEMENTARY ROW POINTS')
    %     disp(['Results for dimension: ' d1str])
    %    if isstruct(Sup) && isfield(Sup,'r')
    %         propsupR=strcat('''LineStyle'',','''none''',',''Marker'',''', symbolsuprows ,''',''Color'',''', colorsuprows , ''',''MarkerFaceColor'',''', colorsuprows ,'''');
    %         eval(['plot(' typeR 'sup(:,d1),' typeR 'sup(:,d2),' propsupR ')'])
    %         eval(['text(' typeR 'sup(:,d1),' typeR 'sup(:,d2),labels.sr,''Interpreter'',''None'')'])
    %
    %
    %     disp(['Results for dimension: ' d1str])
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
    
    % Point2Inertia= relative contribution or each point to the inertia of the
    % dimension.
    % The points  with the larger value of Point2Inertia are those which contribute the
    % most to the definition of the dimension. If the row contributions
    % were uniform, the expected value would be 1/size(contingeny_table,1)
    % For a given dimension, any row with a
    % contribution larger than this threshold could be considered as
    % important in contributing to that dimension.
    
    % Dim2Inertia= contribution of dimension to the inertia of point
    % (where inertia of point is the squared distance of point d_i^2 to the
    % centroid).
    % The quality of representation of the rows on the factor map is called
    % the squared cosine (cos2) or the squared correlations.
    % Dim2Inertia measures the degree of association between points
    % (rows/columns) and a particular axis.
    % The values of the cos2 are forced to be between 0 and 1
    % If a row item is well represented by two dimensions, the sum of the
    % Dim2Inertia is close to one.
end

end
%FScategory:MULT-Categorical