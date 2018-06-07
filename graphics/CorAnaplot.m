function CorAnaplot(out,varargin)
%CorAnaplot draws the Correspondence Analysis (CA) graphs with confidence ellipses.
%
% The typical graphic in a correspondence analysis is to visualize the data
% in a two dimensional space using the first two extracted coordinates from
% both rows and columns. Although we could visualize the rows and the
% columns separately, the usual approach is to plot both in a single
% graphic to get an idea of the association between them.
%
%<a href="matlab: docsearchFS('CorAnaplot')">Link to the help function</a>
%
%  Required input arguments:
%
%  out :  Structure containing the output of function CorAna. Structure.
%               Structure containing the following fields.
%
% 		out.Lr         =  cell of length $I$ containing the labels of
%                         active rows (i.e. the rows which participated to
%                         the fit).
% 		out.Lc         =  cell of length $J$ containing the labels of
%                         active columns (i.e. the columns which participated to
%                         the fit).
% 		out.n         =   Grand total. out.n is equal to sum(sum(out.N)).
%                         This is the number of observations.
% 		out.Dr        =   Square matrix of size $I$ containing on the
%                         diagonal the row masses.  This is matrix $D_r$.
%                         \[
%                           D_r=diag(r)
%                         \]
% 		out.Dc        =   Square matrix of size $J$ containing on the
%                         diagonal the column masses. This is matrix $D_c$.
%                         \[
%                           D_c=diag(c)
%                         \]
% out.InertiaExplained =  matrix with 4 columnn.
%                         - First column contains the singular values (the
%                           sum of the squared singular values is the total
%                           inertia).
%                         - Second column contains the eigenvalues (the sum
%                           of the eigenvalues is the total inertia).
%                         - Third column contains the variance explained by
%                           each latent dimension.
%                         - Fourth column contains the cumulative variance
%                           explained by each dimension.
% 		out.LrSup         =  cell containing the labels of the
%                         supplementary rows (i.e. the rows whicg did not
%                         participate to the fit).
% 		out.LcSup         =  cell containing the labels of
%                         supplementary columns (i.e. the columns which did
%                         not participate to the fit).
% 		out.SupRowsN   =  matrix of size length(LrSup)-by-c
%                         referred to supplementary rows. If there are no
%                         supplementary rows this field is not present.
% 		out.SupColsN  =   matlab of size r-by-length(LcSup) referred to
%                         supplementary columns.
%                         If there are no supplementary columns this field
%                         is not present.
% 	out.RowsPriSup    =   Principal coordinates of supplementary rows.
%                         If there are no supplementary rows this field
%                         is not present.
% 	out.RowsStaSup    =   Standard coordinates of supplementary rows.
%                         If there are no supplementary rows this field
%                         is not present.
%   out.RowsSymSup    =   Symmetrical coordinates of supplementary rows.
%                         If there are no supplementary rows this field
%                         is not present.
% 	out.ColsPriSup    =   Principal coordinates of supplementary columns.
%                         If there are no supplementary columns this field
%                         is not present.
% 	out.ColsStaSup    =   Standard coordinates of of supplementary columns.
%                         If there are no supplementary columns this field
%                         is not present.
%   out.ColsSymSup    =   Symmetrical coordinates of supplementary columns.
%                         If there are no supplementary columns this field
%                         is not present.
%              Data Types - struct
%
%
%  Optional input arguments:
%
%       plots : Customize plot appearance. Scalar or structure.
%               If plots is not a structure, a plot which shows the Principal
%               coordinates of rows and columns is shown on the screen. If
%               plots is a structure it may contain the following fields:
%               plots.alpha = type of plot, scalar in the interval [0 1] or
%               a string identifying the type of coordinates to use in the
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
%                   chi-squared distances (column-metric-preserving). The
%                   position of the column points are at the weighted
%                   average of the row points.
%                   Note that 'colwprincipal' can also be
%                   specified setting plots.alpha=0.
%               If $plots.alpha='symbiplot'$, the row and column coordinates
%                   are scaled similarly. The sum of weighted squared
%                   coordinates for each dimension is equal to the
%                   corresponding singular values. These coordinates are often
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
%                   This biplot has been suggested by Greenacre and
%                   incorporates the contribution of points. In this
%                   display, points that contribute very little to the
%                   solution, are close to the center of the biplot and are
%                   relatively unimportant to the interpretation. This
%                   biplot is often referred as contribution biplot because
%                   it shows visually the most contributing points
%                   (Greenacre 2006b).
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
%              plots.FontSize = scalar which specifies the font size of row
%                   (column) labels. The default value is 10.
%              plots.MarkerSize = scalar which specifies the marker size
%                   of symbols associated with rows or columns. The default
%                   value is 10.
%              plots.SymbolRows= character which specifies the symbol to
%                   use for row points. If this field is not present the
%                   default symbols is 'o'.
%              plots.SymbolCols= character which specifies the symbol to
%                   use for column points. If this field is not present the
%                   default symbols is '^'.
%              plots.SymbolRowsSup= character which specifies the symbol to
%                   use for supplementary row points. If this field is not
%                   present the default symbols is 'o'.
%              plots.SymbolColsSup= character which specifies the symbol to
%                   use for supplementary column points. If this field is
%                   not present the default symbols is '^'.
%              plots.ColorRows= character which specifies the symbol to
%                   use for row points or RGB triplet. If this field is not
%                   present the default color is 'b'.
%              plots.ColorCols= character which specifies the symbol to
%                   use for column points or RGB triplet. If this field is
%                   not present the default color is 'r'.
%              plots.ColorRowsSup= character which specifies the symbol to
%                   use for row points or RGB triplet. If this field is not
%                   present the default color is 'b'.
%              plots.ColorColsSup= character which specifies the symbol to
%                   use for supplementary column points or RGB triplet. If
%                   this field is not present the default color is 'r'.
%              plots.MarkerFaceColorRows= character which specifies the
%                   marker fill color to use for active row points or RGB
%                   triplet. If this field is not present the default color
%                   is 'auto'.
%              plots.MarkerFaceColorCols= character which specifies the
%                   marker fill color to use for active column points or RGB
%                   triplet. If this field is not present the default color
%                   is 'auto'.
%              plots.MarkerFaceColorRowsSup= character which specifies the
%                   marker fill color to use for supplementary row points or RGB
%                   triplet. If this field is not present the default color
%                   is 'auto'.
%              plots.MarkerFaceColorColsSup= character which specifies the
%                   marker fill color to use for supplementary column points or RGB
%                   triplet. If this field is not present the default color
%                   is 'auto'.
%              Example - 'plots',plots=struct; plots.colorcols='k'
%              Data Types - double
%       addx : horizontal displacement for labels. Scalar. Amount of
%              horizontal displacement which has been put on the labels in the
%              plot. The defalut value of addx is 0.04.
%              Example - 'addx',0.01
%              Data Types - double
%changedimsign: change chosen dimension sign. Boolean vector of length 2.
%               Sometimes for better interpretability it is necessary to
%               change the sign of the coordinates for the chosen
%               dimension. If changedimsign(1) is true the sign of the
%               coordinates for first chosen dimension is changed. If
%               changedimsign(2) is true the sign of the coordinates for
%               first chosen dimension is changed. As default the
%               dimensions are the first and the second however, they can
%               be changed using option plots.dim. The defaul value of
%               changedimsign is [false false] that is the sign is not
%               changed.
%              Example - 'changedimsign', [true false]
%              Data Types - boolean
% confellipse : confidence ellipses around rows and/or columns points.
%               Scalar or struct.
%               If confellipse is 1, 90 per cent confidence ellipses are
%               drawn around each row and column point based on multinomial method.
%               If confellipse is a struct it may contain the followgin
%               fields.
%               confellipse.conflev= number in the interval (0 1) which
%                   defines the confidence level of each ellipse.
%               confellipse.method= cell which specifies the method(s) to use to
%                   compute confidence ellipses. Possible values are:
%                   {'multinomial'} = in this case the original contigencey
%                   table with the active elements is taken as a reference.
%                   Then new data tables are drawn in the following way:
%                   $r\times c$ values are drawn from a multinomial
%                   distribution with theoretical frequencies equals to
%                   $n_{ij}/n$.
%                   {'bootrows'} = the values are bootstrapped row by row:
%                   Given row i, $n_{i.}$ are extracted with repetition and
%                   a frequencey distribution is computed using classes
%                   $[0, n_{i1}]$,$[n_{i1}, n_{i1}+n_{i2}]$, $\ldots$
%                   $[\sum_{j=1}^{J-1} n_{ij}, \sum_{j=1}^{J} n_{ij}$.
%                   {'bootcols'} = the values are bootstrapped column by
%                   column. If  confellipse.method for example is
%                   {'bootrows' 'bootcols'} two ellipses are drawn to each
%                   point. In this case it is possible to appreciate the
%                   stability of both methods.
%               confellipse.nsimul=scalar which defines the number of
%                   contingency tables which have to ge generated. The
%                   default value of confellipse.nsimul is 1000. Thus
%                   nsimul new contingency tables are projected as
%                   supplementary rows and/or supplementary columns.
%               confellipse.selrows= vector which specifies for which row
%                   points it is necessary to draw the ellipses.
%                   confellipse.selRows either a boolean vector of length
%                   I containing a true in correspondence of the row
%                   elements for which the ellipse has to be drawn or a
%                   numeric vector which contains the indexes of the units
%                   which have to be drawn or a cell arrary containing the
%                   names of the rows for which the ellipse has to be drawn.
%                   For example if I=3 and the second row is
%                   called 'row2' in order to show just the confidence
%                   ellipse for this row it is possible to use
%                   confellipse.selRows=[false true false], or
%                   confellipse.selRows=2 or confellipse.selRows={'row2'}.
%               confellipse.selCols= vector which specifies for which
%                   column points it is necessary to draw the ellipses.
%                   confellipse.selCols can be either a boolean vector of
%                   length J containing a true in correspondence of the
%                   column elements for which the ellipse has to be drawn
%                   or a numeric vector which contains the indexes of the
%                   columns which have to be drawn or a cell arrary
%                   containing the names of the columns for which the
%                   ellipse has to be drawn. For example if J=3 and one the
%                   third column is called 'Col3' in order to show just the
%                   confidence ellipse for this element it is possible to
%                   use confellipse.selCols=[false false true], or
%                   confellipse.selCols=3 or confellipse.selCols={'Col3'}.
%               confellipse.AxesEllipse = boolean which specifies whether
%                   it is necessary to show the major axes of the ellipse.
%                   The default value of confellipse.AxesEllipse is true
%                   that is the axes are shown.
%               Example - 'confellipse', 0
%               Data Types - scalar or struct
%       xlimx   :   Min and Max of the x axis. Vector. Vector with two
%                   elements controlling minimum and maximum
%                   of the x axis.
%                   Example - 'xlimx',[-1 1]
%                   Data Types - double
%       ylimy   :   Min and Max of the y axis. Vector. Vector with two
%                   elements controlling minimum and
%                   maximum of the y axis.
%                   Example - 'ylimy',[0 1]
%                   Data Types - double
%        d1    :  Dimension to show on the horizontal axis. Positive
%                 integer. Positive integer in the range 1, 2, .., K which
%                 indicates the dimension to show on the x axis. The
%                 default value of d1 is 1.
%                 Example - 'd1',2
%                 Data Types - single | double
%        d2    :  Dimension to show on the vertical axis. Positive
%                 integer. Positive integer in the range 1, 2, .., K which
%                 indicates the dimension to show on the y axis. The
%                 default value of d2 is 2.
%                 Example - 'd2',3
%                 Data Types - single | double%
%           h : the axis handle of a figure where to send the CorAnaplot.
%               This can be used to host the CorAnaplot in a subplot of a
%               complex figure formed by different panels (for example a panel
%               with CorAnaplot from plots.alpha=0.2 and another
%               with CorAnaplot from plots.alpha=0.5).
%               Example -'h',h1 where h1=subplot(2,1,1)
%               Data Types - Axes object (supplied as a scalar)
% Output:
%
%
% See also: corAna
%
% References:
% Benzecri, J.-P. (1992), "Correspondence Analysis Handbook", New-York,
% Dekker.
% Benzecri, J.-P. (1980), "L'analyse des donnees tome 2: l'analyse des
% correspondances", Paris, Bordas.
% Greenacre, M.J. (1993), "Correspondence Analysis in Practice", London,
% Academic Press.
% Gabriel, K.R. and Odoroff, C. (1990), Biplots in biomedical research,
% "Statistics in Medicine", Vol. 9, pp. 469-485.
% Greenacre, M.J. (1993), Biplots in correspondence Analysis, "Journal of
% Applied Statistics", Vol. 20, pp. 251-269.
% Urbano, L.-S., van de Velden, M. and Kiers, H.A.L. (2009),
% CAR: A MATLAB Package to Compute Correspondence Analysis with Rotations,
% "Journal of Statistical Software", Vol. 31.
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('CorAnaplot')">Link to the help function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit

% Examples:
%


%{
    %% CorAnaplot with all the default options.
    N=[51	64	32	29	17	59	66	70;
        53	90	78	75	22	115	117	86;
        71	111	50	40	11	79	88	177;
        1	7	5	5	4	9	8	5;
        7	11	4	3	2	2	17	18;
        7	13	12	11	11	18	19	17;
        21	37	14	26	9	14	34	61;
        12	35	19	6	7	21	30	28;
        10	7	7	3	1	8	12	8;
        4	7	7	6	2	7	6	13;
        8	22	7	10	5	10	27	17;
        25	45	38	38	13	48	59	52;
        18	27	20	19	9	13	29	53;
        35	61	29	14	12	30	63	58;
        2	4	3	1	4	nan  nan	nan	  ;
        2	8	2	5	2	nan  nan	nan;
        1	5	4	6	3	nan  nan	nan;
        3	3	1	3	4	nan  nan	nan];
    % rowslab = cell containing row labels
    rowslab={'money','future','unemployment','circumstances',...
        'hard','economic','egoism','employment','finances',...
        'war','housing','fear','health','work','comfort','disagreement',...
        'world','to_live'};
    % colslab = cell containing column labels
    colslab={'unqualified','cep','bepc','high_school_diploma','university',...
        'thirty','fifty','more_fifty'};
    if verLessThan('matlab','8.2.0')==0
        tableN=array2table(N,'VariableNames',colslab,'RowNames',rowslab);
        % Extract just active rows
        Nactive=tableN(1:14,1:5);
        Nsupr=tableN(15:18,1:5);
        Nsupc=tableN(1:14,6:8);
        Sup=struct;
        Sup.r=Nsupr;
        Sup.c=Nsupc;
        % Compute Correspondence analysis
    else
        Nactive=N(1:14,1:5);
        Lr=rowslab(1:14);
        Lc=colslab(1:5);
        Sup=struct;
        Sup.r=N(15:end,1:5);
        Sup.Lr=rowslab(15:end);
        Sup.c=N(1:14,6:8);
        Sup.Lc=colslab(6:8);
        
        out=CorAna(Nactive,'Sup',Sup,'plots',0,'dispresults',false);
        % Rows and columns are showN in principal coordinates
        CorAnaplot(out)
    end
%}

%{
    %% CorAnaplot with personalized symbols.
    close all
    % Six-pointed star (hesagram) for supplementary rows
    SymbolRows='h';
    % Five-pointed star (pentagram) for supplementary rows
    SymbolRowsSup='p';
    % Color for active rows
    ColorRows='b';
    % Color for supplementary rows (dark blue)
    ColorRowsSup=[6 13 123]/255;
    % Blue fill color for active rows
    MarkerFaceColorRows='b';

    % Right-pointing triangle for active columns
    Symbolcols='^';
    % Six-pointed star (hexagram) for supplementary columns
    SymbolcolsSup='h';
    % Color for active columns
    ColorCols='r';
    % Red fill color for active rows
    MarkerFaceColorCols='r';

    % Color for supplementary columns (dark red)
    ColorColsSup=[128 0 0]/255;
    plots=struct;
    plots.SymbolRows=SymbolRows;
    plots.SymbolRowsSup=SymbolRowsSup;
    plots.ColorRows=ColorRows;
    plots.ColorRowsSup=ColorRowsSup;
    plots.MarkerFaceColorRows=MarkerFaceColorRows;

    plots.SymbolCols=Symbolcols;
    plots.SymbolColsSup=SymbolcolsSup;
    plots.ColorCols=ColorCols;
    plots.ColorColsSup=ColorColsSup;
    plots.MarkerFaceColorCols=MarkerFaceColorCols;

    % change the sign of the second dimension
    changedimsign=[false true];
    CorAnaplot(out,'plots',plots,'changedimsign',changedimsign)
%}

%{
    % CorAnaplot with personalized displacement.
    close all
    % No horizontal displacement for the labels.
    addx=0;
    CorAnaplot(out,'addx',addx)
%}

%{
    %% Correpondence analysis plot with selected ellipses.
    N=[51	64	32	29	17	59	66	70;
        53	90	78	75	22	115	117	86;
        71	111	50	40	11	79	88	177;
        1	7	5	5	4	9	8	5;
        7	11	4	3	2	2	17	18;
        7	13	12	11	11	18	19	17;
        21	37	14	26	9	14	34	61;
        12	35	19	6	7	21	30	28;
        10	7	7	3	1	8	12	8;
        4	7	7	6	2	7	6	13;
        8	22	7	10	5	10	27	17;
        25	45	38	38	13	48	59	52;
        18	27	20	19	9	13	29	53;
        35	61	29	14	12	30	63	58;
        2	4	3	1	4	nan  nan	nan	  ;
        2	8	2	5	2	nan  nan	nan;
        1	5	4	6	3	nan  nan	nan;
        3	3	1	3	4	nan  nan	nan];
    % rowslab = cell containing row labels
    rowslab={'money','future','unemployment','circumstances',...
        'hard','economic','egoism','employment','finances',...
        'war','housing','fear','health','work','comfort','disagreement',...
        'world','to_live'};
    % colslab = cell containing column labels
    colslab={'unqualified','cep','bepc','high_school_diploma','university',...
        'thirty','fifty','more_fifty'};
    if verLessThan('matlab','8.2.0')==0
        tableN=array2table(N,'VariableNames',colslab,'RowNames',rowslab);
        % Extract just active rows
        Nactive=tableN(1:14,1:5);
        Nsupr=tableN(15:18,1:5);
        Nsupc=tableN(1:14,6:8);
        Sup=struct;
        Sup.r=Nsupr;
        Sup.c=Nsupc;
    else
        Nactive=N(1:14,1:5);
        Lr=rowslab(1:14);
        Lc=colslab(1:5);
        Sup=struct;
        Sup.r=N(15:end,1:5);
        Sup.Lr=rowslab(15:end);
        Sup.c=N(1:14,6:8);
        Sup.Lc=colslab(6:8);
    end

    % Superimpose confidence ellipses for rows 2 and 4 and for column 3
    confellipse=struct;
    confellipse.selRows=[2 4];
    % Ellipse for column 3 using an interger
    confellipse.selCols=3;
    % Ellipse for column 3 using a Boolean vector
    confellipse.selCols=[ false false true false false];
    % confellipse.selcols={'c3'};
    % Use the 3 methods below in order to compute the confidence ellipses for
    % the selected rows and columns of the input contingency table
    confellipse.method={'multinomial' 'bootRows' 'bootCols'};
    % Set number of simulations
    confellipse.nsimul=500;
    % Set confidence interval
    confellipse.conflev=0.50;
    out=CorAna(Nactive,'Sup',Sup,'plots',0,'dispresults',false);

    % Draw correspondence analysis plot with requested confidence ellipses
    CorAnaplot(out,'plots',1,'confellipse',confellipse)
%}

%{
    % Correpondence analysis plot with ellipses only on column points.
    N=[51	64	32	29	17	59	66	70;
        53	90	78	75	22	115	117	86;
        71	111	50	40	11	79	88	177;
        1	7	5	5	4	9	8	5;
        7	11	4	3	2	2	17	18;
        7	13	12	11	11	18	19	17;
        21	37	14	26	9	14	34	61;
        12	35	19	6	7	21	30	28;
        10	7	7	3	1	8	12	8;
        4	7	7	6	2	7	6	13;
        8	22	7	10	5	10	27	17;
        25	45	38	38	13	48	59	52;
        18	27	20	19	9	13	29	53;
        35	61	29	14	12	30	63	58;
        2	4	3	1	4	nan  nan	nan	  ;
        2	8	2	5	2	nan  nan	nan;
        1	5	4	6	3	nan  nan	nan;
        3	3	1	3	4	nan  nan	nan];
    % rowslab = cell containing row labels
    rowslab={'money','future','unemployment','circumstances',...
        'hard','economic','egoism','employment','finances',...
        'war','housing','fear','health','work','comfort','disagreement',...
        'world','to_live'};
    % colslab = cell containing column labels
    colslab={'unqualified','cep','bepc','high_school_diploma','university',...
        'thirty','fifty','more_fifty'};
    if verLessThan('matlab','8.2.0')==0
        tableN=array2table(N,'VariableNames',colslab,'RowNames',rowslab);
        % Extract just active rows
        Nactive=tableN(1:14,1:5);
        Nsupr=tableN(15:18,1:5);
        Nsupc=tableN(1:14,6:8);
        Sup=struct;
        Sup.r=Nsupr;
        Sup.c=Nsupc;
    else
        Nactive=N(1:14,1:5);
        Lr=rowslab(1:14);
        Lc=colslab(1:5);
        Sup=struct;
        Sup.r=N(15:end,1:5);
        Sup.Lr=rowslab(15:end);
        Sup.c=N(1:14,6:8);
        Sup.Lc=colslab(6:8);
    end
    % Superimpose confidence ellipses 
    confellipse=struct;
    % No confdence ellipse for row points
    confellipse.selRows=[];
    % Ellipse for all the column points using a Boolean vector
    confellipse.selCols=[ true true true true true];

    % Compare methods 'multinomial' and 'bootCols'
    confellipse.method={'multinomial' 'bootCols'};
    % Set number of simulations
    confellipse.nsimul=10000;
    % Set confidence interval
    confellipse.conflev=0.90;
    out=CorAna(Nactive,'Sup',Sup,'plots',0,'dispresults',false);

    % Draw correspondence analysis plot with requested confidence ellipses
    CorAnaplot(out,'plots',1,'confellipse',confellipse)
%}

%{
    % Correpondence analysis plot using latent dimensions 3 and 4.
    N=[51	64	32	29	17	59	66	70;
        53	90	78	75	22	115	117	86;
        71	111	50	40	11	79	88	177;
        1	7	5	5	4	9	8	5;
        7	11	4	3	2	2	17	18;
        7	13	12	11	11	18	19	17;
        21	37	14	26	9	14	34	61;
        12	35	19	6	7	21	30	28;
        10	7	7	3	1	8	12	8;
        4	7	7	6	2	7	6	13;
        8	22	7	10	5	10	27	17;
        25	45	38	38	13	48	59	52;
        18	27	20	19	9	13	29	53;
        35	61	29	14	12	30	63	58;
        2	4	3	1	4	nan  nan	nan	  ;
        2	8	2	5	2	nan  nan	nan;
        1	5	4	6	3	nan  nan	nan;
        3	3	1	3	4	nan  nan	nan];
    % rowslab = cell containing row labels
    rowslab={'money','future','unemployment','circumstances',...
        'hard','economic','egoism','employment','finances',...
        'war','housing','fear','health','work','comfort','disagreement',...
        'world','to_live'};
    % colslab = cell containing column labels
    colslab={'unqualified','cep','bepc','high_school_diploma','university',...
        'thirty','fifty','more_fifty'};
    if verLessThan('matlab','8.2.0')==0
        tableN=array2table(N,'VariableNames',colslab,'RowNames',rowslab);
        % Extract just active rows
        Nactive=tableN(1:14,1:5);
        Nsupr=tableN(15:18,1:5);
        Nsupc=tableN(1:14,6:8);
        Sup=struct;
        Sup.r=Nsupr;
        Sup.c=Nsupc;
    else
        Nactive=N(1:14,1:5);
        Lr=rowslab(1:14);
        Lc=colslab(1:5);
        Sup=struct;
        Sup.r=N(15:end,1:5);
        Sup.Lr=rowslab(15:end);
        Sup.c=N(1:14,6:8);
        Sup.Lc=colslab(6:8);
    end
    % Superimpose ellipses 
    confellipse=struct;
    % Ellipse for the first 3 row points
    confellipse.selRows=1:3;
    % Ellipse for selected column points using a Boolean vector
    confellipse.selCols=[ false false true true false];

    % Compare methods 'multinomial' and 'bootCols'
    confellipse.method={'multinomial' 'bootCols'};
    % Set number of simulations
    confellipse.nsimul=10000;
    % Set confidence interval
    confellipse.conflev=0.90;
    d1=3;
    d2=4;
    out=CorAna(Nactive,'Sup',Sup,'plots',0,'dispresults',false,'d1',d1,'d2',d2);

    % Draw correspondence analysis plot with requested confidence ellipses
    CorAnaplot(out,'plots',1,'confellipse',confellipse,'d1',d1,'d2',d2)
%}

%{
    %% Correpondence analysis of the smoke data.
    % In this example we compare the results which are obtained using
    % option  plots.alpha='colprincipal'; (which implicitly implies
    % alpha=0) with those which come out imposing directly plots.alpha=0.
    load smoke
    X=smoke.data;
    [N,~,~,labels]=crosstab(X(:,1),X(:,2));
    [I,J]=size(N);
    labels_rows=labels(1:I,1);
    labels_columns=labels(1:J,2);

    out=CorAna(N,'Lr',labels_rows,'Lc',labels_columns,'plots',0,'dispresults',false);
    plots=struct;
    plots.alpha='rowgab';
    plots.alpha='colgab';
    plots.alpha='rowgreen';
    plots.alpha='colgreen';
    % Add confidence ellipses
    confellipse=1;
    plots.alpha='bothprincipal';
    plots.alpha='rowprincipal';
    plots.alpha='colprincipal';
    h1=subplot(1,2,1);
    CorAnaplot(out,'plots',plots,'confellipse',confellipse,'h',h1)
    h2=subplot(1,2,2);
    plots.alpha=0;
    CorAnaplot(out,'plots',plots,'confellipse',confellipse,'h',h2);
%}

%% Initialization
% Default font size for labels of rows or colums to add to the plot
FontSizedef=10;
MarkerSizedef=10;

options=struct('plots',1,'xlimx','','ylimy','','changedimsign',[false false],...
    'addx',0.04,'confellipse',0,'d1',1,'d2',2,'h','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:CorAnaplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 1
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

plots=options.plots;

Lr=out.Lr;
Lc=out.Lc;
LrSup=out.LrSup;
LcSup=out.LcSup;

h=options.h;

InertiaExplained=out.InertiaExplained;
Gam=diag(InertiaExplained(:,1));
Dr=out.Dr;
Dc=out.Dc;
n=out.n;

FontName='Times';
FontSizeAxisLabels=12;

d1  = options.d1;
d2  = options.d2;

dim=[d1 d2];



RowsSta=out.RowsSta;
RowsPri=out.RowsPri;
RowsSym=out.RowsSym;

ColsSta=out.ColsSta;
ColsPri=out.ColsPri;
ColsSym=out.ColsSym;

% Check if the sign of the chosen dimensions has to be changed
changedimsign=options.changedimsign;

for i=1:length(changedimsign)
    if changedimsign(i) == true
        RowsSta(:,dim(i))=-RowsSta(:,dim(i));
        RowsPri(:,dim(i))=-RowsPri(:,dim(i));
        RowsSym(:,dim(i))=-RowsSym(:,dim(i));
        ColsSta(:,dim(i))=-ColsSta(:,dim(i));
        ColsPri(:,dim(i))=-ColsPri(:,dim(i));
        ColsSym(:,dim(i))=-ColsSym(:,dim(i));
    end
end

if ~isempty(LrSup)
    RowsStaSup=out.RowsStaSup;
    RowsPriSup=out.RowsPriSup;
    RowsSymSup=out.RowsSymSup;
    
    for i=1:length(changedimsign)
        if changedimsign(i) == true
            RowsStaSup(:,dim(i))=-RowsStaSup(:,dim(i));
            RowsPriSup(:,dim(i))=-RowsPriSup(:,dim(i));
            RowsSymSup(:,dim(i))=-RowsSymSup(:,dim(i));
        end
    end
    SupRowsN=out.SupRowsN;
end

if ~isempty(LcSup)
    ColsStaSup=out.ColsStaSup;
    ColsPriSup=out.ColsPriSup;
    ColsSymSup=out.ColsSymSup;
    
    for i=1:length(changedimsign)
        if changedimsign(i) == true
            ColsStaSup(:,dim(i))=-ColsStaSup(:,dim(i));
            ColsPriSup(:,dim(i))=-ColsPriSup(:,dim(i));
            ColsSymSup(:,dim(i))=-ColsSymSup(:,dim(i));
        end
    end
    
    SupColsN=out.SupColsN;
end

if isstruct(plots)
    plotsdef=struct;
    plotsdef.alpha=1;
    plotsdef.FontSize=10;
    plotsdef.MarkerSize=10;
    plotsdef.SymbolRows='o';
    plotsdef.SymbolCols='^';
    plotsdef.SymbolRowsSup='o';
    plotsdef.SymbolColsSup='^';
    plotsdef.ColorRows='b';
    plotsdef.ColorRowsSup='b';
    plotsdef.ColorCols='r';
    plotsdef.ColorColsSup='r';
    plotsdef.MarkerFaceColorRows='auto';
    plotsdef.MarkerFaceColorCols='auto';
    plotsdef.MarkerFaceColorRowsSup='r';
    plotsdef.MarkerFaceColorColsSup='r';
    
    fld=fieldnames(plots);
    
    % Check if user options inside plots are valid options
    chkoptions(plotsdef,fld)
    
    % This anonymous function anables to extract the variable name to a
    % string
    ExtractVariableName=@(x) inputname(1);
    
    if isfield(plots,'alpha')
        if strcmp(plots.alpha,'rowprincipal')
            typeR='RowsPri'; % rows are in principal coordinates
            typeC='ColsSta'; % columns are in standard coordinates
            titl={'Rows principal coordinates, and columns standard coordinates' , ...
                '$\alpha=1$, $X=D_r^{-1/2}U\Gamma$ and $Y= D_c^{-1/2} V$'};
            typeRSup='RowsPriSup';
            typeCSup='ColsStaSup';
            
        elseif strcmp(plots.alpha,'colprincipal')
            typeR='RowsSta'; % rows are in standard coordinates
            typeC='ColsPri'; % columns are in principal coordinates
            titl={'Rows standard coordinates, and columns principal coordinates' , ...
                '$\alpha=0$, $X=D_r^{-1/2}U $ and $G= D_c^{-1/2} V \Gamma$'};
            typeRSup='RowsStaSup';
            typeCSup='ColsPriSup';
            
            
        elseif strcmp(plots.alpha,'symbiplot')
            % equivalent to alpha=0.5
            typeR='RowsSym';        % rows are in symmetrical coordinates
            typeC='ColsSym';        % columns are in symmetrical coordinates
            titl='Biplot symmetrical model $\alpha=0.5$ $X=D_r^{-1/2}U\Gamma^{1/2} $ and $Y= D_c^{-1/2} \Gamma V^{1/2}$';
            typeRSup='RowsSymSup';
            typeCSup='ColsSymSup';
            
            
        elseif strcmp(plots.alpha,'bothprincipal')
            typeR='RowsPri';        % rows are in principal coordinates
            typeC='ColsPri';        % columns are in principal coordinates
            titl={'French symmetrical model: rows and cols in principal coordinates.' , ...
                'Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$'};
            typeRSup='RowsPriSup';
            typeCSup='ColsPriSup';
            
            
        elseif strcmp(plots.alpha,'rowgab')
            %  If plots.alpha='rowgab'  rows are in principal coordinates
            %  and columns are in standard coordinates multiplied by the
            %  mass.
            typeR='RowsPri';        % rows are in principal coordinates
            ColsStaDc=Dc*ColsSta;
            typeC=ExtractVariableName(ColsStaDc);
            titl={'Rows principal coordinates, and columns standard coordinates times masses' , ...
                '$X=D_r^{-1/2}U\Gamma $ and $Y= D_c^{1/2} V$'};
            typeRSup='RowsPriSup';
            
            if ~isempty(LcSup)
                ColsStaDcSup=diag(sum(SupColsN,1)/n)*ColsStaSup;
                typeCSup=ExtractVariableName(ColsStaDcSup);
            end
            
            
        elseif strcmp(plots.alpha,'colgab')
            % If plots.alpha='colgab'  columns are in principal coordinates
            % and rows are in standard coordinates multiplied by the
            % masses.
            RowsStaDr=Dr*RowsSta;
            typeR=ExtractVariableName(RowsStaDr);
            typeC='ColsPri';        % columns are in principal coordinates
            titl={'Rows standard coordinates multiplied by masses ' , ...
                'and columns principal coordinates $X=D_r^{1/2} U$ and $Y= D_c^{-1/2} V \Gamma$'};
            
            if ~isempty(LrSup)
                RowsStaDrSup=diag(sum(SupRowsN,2)/n)*RowsStaSup;
                typeRSup=ExtractVariableName(RowsStaDrSup);
            end
            typeCSup='ColsPriSup';
            
        elseif strcmp(plots.alpha,'rowgreen')
            %  If plots.alpha='rowgreen'  rows are in principal
            %  coordinates and columns are in standard coordinates
            %  multiplied by square root of the mass.
            typeR='RowsPri';        % rows are in principal coordinates
            ColsStaDcSqrt=(Dc^(1/2))*ColsSta;
            typeC= ExtractVariableName(ColsStaDcSqrt);
            titl={'Rows principal coordinates, and column standard coordinates ' , ...
                'times sqrt of masses $X=D_r^{-1/2}U\Gamma $ and $Y= V$'};
            typeRSup='RowsPriSup';        % rows are in principal coordinates
            
            if ~isempty(LcSup)
                ColsStaDcSqrtSup=(diag(sum(SupColsN,1)/n))^(1/2)*ColsStaSup;
                typeCSup= ExtractVariableName(ColsStaDcSqrtSup);
            end
            
        elseif strcmp(plots.alpha,'colgreen')
            %  If plots.alpha='colgreen' columns in principal coordinates
            %  and rows in standard coordinates multiplied by the square
            %  root of the mass.
            RowsStaDrSqrt=(sqrt(Dr))*RowsSta;
            typeR=ExtractVariableName(RowsStaDrSqrt);
            typeC='ColsPri';        % columns are in principal coordinates
            titl={'Rows standard coordinates times sqrt of masses,' ...
                'and columns principal coordinates, $X=U $ and $G= D_c^{-1/2} V \Gamma$'};
            if ~isempty(LrSup)
                RowsStaDrSqrtSup=sqrt(diag(sum(SupRowsN,2)/n))*RowsStaSup;
                typeRSup=ExtractVariableName(RowsStaDrSqrtSup);
            end
            typeCSup='ColsPriSup';        % columns are in principal coordinates
        else
            if isnumeric(plots.alpha)
                if plots.alpha>=0 && plots.alpha<=1
                    RowsAlpha= RowsSta*Gam^plots.alpha;
                    ColsAlpha= ColsSta*Gam^(1-plots.alpha);
                    typeR=ExtractVariableName(RowsAlpha);
                    typeC=ExtractVariableName(ColsAlpha);
                    titl=['$\alpha='  num2str(plots.alpha) '\qquad   X=D_r^{-1/2} U \Gamma^{' num2str(plots.alpha) '}$ and $Y= D_c^{-1/2} V \Gamma^{1-'  num2str(plots.alpha) '}$'];
                    if ~isempty(LrSup)
                        hm=bsxfun(@rdivide,SupRowsN,sum(SupRowsN,2));
                        RowsAlphaSup= hm*ColsSta *Gam^(plots.alpha-1);
                        typeRSup=ExtractVariableName(RowsAlphaSup);
                    end
                    if ~isempty(LcSup)
                        hm=bsxfun(@rdivide,SupColsN,sum(SupColsN,1));
                        ColsAlphaSup= hm'*RowsSta* Gam^(-plots.alpha);
                        typeCSup=ExtractVariableName(ColsAlphaSup);
                    end
                else
                    error('FSDA:CorAnaplot:WrongInputOpt','Value of plots.alpha must lie in the interval [0 1]')
                end
            else
                listStrings={'rowprincipal'; 'colprincipal'; 'symbiplot'; 'bothprincipal'; 'rowgab'; 'colgab'; 'rowgreen'; 'colgreen'};
                warning('FSDA:CorAnaplot:WrongInputOpt',['Input string ''' plots.alpha ''' is  not found'])
                disp('Possible strings are')
                disp(listStrings)
                error('FSDA:CorAnaplot:WrongInputOpt','Please use one of the above strings')
            end
        end
        
    else
        typeR='RowsPri';        % rows are in principal coordinates
        typeC='ColsPri';        % columns are in principal coordinates
        titl={'French symmetrical model: rows and cols in principal coordinates.' ...
            'Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$'};
        typeRSup='RowsPriSup';
        typeCSup='ColsPriSup';
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
    if isfield(plots,'SymbolRows')
        SymbolRows=plots.SymbolRows;
    else
        SymbolRows='o';
    end
    
    if isfield(plots,'SymbolRowsSup')
        SymbolRowsSup=plots.SymbolRowsSup;
    else
        SymbolRowsSup='o';
    end
    
    if isfield(plots,'SymbolCols')
        SymbolCols=plots.SymbolCols;
    else
        SymbolCols='^';
    end
    
    if isfield(plots,'SymbolColsSup')
        SymbolColsSup=plots.SymbolColsSup;
    else
        SymbolColsSup='^';
    end
    
    if isfield(plots,'ColorRows')
        ColorRows=plots.ColorRows;
    else
        ColorRows='b';
    end
    
    if isfield(plots,'ColorRowsSup')
        ColorRowsSup=plots.ColorRowsSup;
    else
        ColorRowsSup='r';
    end
    
    if isfield(plots,'MarkerFaceColorRows')
        MarkerFaceColorRows=plots.MarkerFaceColorRows;
    else
        MarkerFaceColorRows='auto';
    end
    
    if isfield(plots,'MarkerFaceColorRowsSup')
        MarkerFaceColorRowsSup=plots.MarkerFaceColorRowsSup;
    else
        MarkerFaceColorRowsSup=ColorRowsSup;
    end
    
    if isfield(plots,'ColorCols')
        ColorCols=plots.ColorCols;
    else
        ColorCols='r';
    end
    if isfield(plots,'ColorColsSup')
        ColorColsSup=plots.ColorColsSup;
    else
        ColorColsSup='r';
    end
    
    if isfield(plots,'MarkerFaceColorCols')
        MarkerFaceColorCols=plots.MarkerFaceColorCols;
    else
        MarkerFaceColorCols='auto';
    end
    
    if isfield(plots,'MarkerFaceColorColsSup')
        MarkerFaceColorColsSup=plots.MarkerFaceColorColsSup;
    else
        MarkerFaceColorColsSup=ColorColsSup;
    end
    
else
    typeR='RowsPri';        % rows are in principal coordinates
    typeC='ColsPri';        % columns are in principal coordinates
    titl={'French symmetrical model: rows and cols in principal coordinates.'...
        'Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$'};
    FontSize=FontSizedef;
    MarkerSize=MarkerSizedef;
    typeRSup='RowsPriSup';
    typeCSup='ColsPriSup';
    SymbolRows='o';
    SymbolCols='^';
    SymbolRowsSup='o';
    SymbolColsSup='^';
    % Color for symbols and text for rows points
    ColorRows='b';
    % Color for symbols and text for supplementary row points
    ColorRowsSup='b';
    % Color for symbols and text for  active column points
    ColorCols='r';
    % Color for symbols and text for supplementary column points
    ColorColsSup='r';
    % Marker fill color for active rows
    MarkerFaceColorRows='auto';
    % Marker fill color for active columns
    MarkerFaceColorCols='auto';
    % Marker fill color for supplementary rows
    MarkerFaceColorRowsSup=ColorRowsSup;
    % Marker fill color for supplementary columns
    MarkerFaceColorColsSup=ColorColsSup;
end

Carows= eval(typeR);
Cacols= eval(typeC);

% Create the figure that will host the CorAnaplot
hfig = figure('Name', 'Correspondence analysis plot', 'NumberTitle', 'off',...
    'Tag','pl_CorAna');

% Get figure's axis
afig = axes('Parent',hfig);

hold('on')
% Plot row points
plot(afig,Carows(:,d1),Carows(:,d2),'LineStyle','none','Marker',SymbolRows ,'Color', ColorRows , ...
    'MarkerSize',MarkerSize,'MarkerFaceColor', MarkerFaceColorRows)

% Plot column points
plot(afig,Cacols(:,d1),Cacols(:,d2),'LineStyle','none','Marker',SymbolCols ,'Color', ColorCols , ...
    'MarkerSize',MarkerSize,'MarkerFaceColor', MarkerFaceColorCols)

% Add labels for row points and column points
% addx = adds a small right horizontal displacement for labels
addx=options.addx;

text(Carows(:,d1)+addx , Carows(:,d2), Lr,'Interpreter','None','FontSize',FontSize,'Color', ColorRows )
text(Cacols(:,d1)+addx , Cacols(:,d2), Lc,'Interpreter','None','FontSize',FontSize,'Color', ColorCols )


% Add points and text associated to supplementary rows
if ~isempty(LrSup)
    CarowsSup= eval(typeRSup);
    plot(afig,CarowsSup(:,d1),CarowsSup(:,d2),'LineStyle','none','Marker',SymbolRowsSup ,...
        'Color', ColorRowsSup , 'MarkerFaceColor', MarkerFaceColorRowsSup,'MarkerSize',MarkerSize)
    text(CarowsSup(:,d1)+addx , CarowsSup(:,d2), LrSup,'Interpreter','None','FontSize',FontSize,'Color', ColorRowsSup )
end

% Add points and text associated to supplementary columns
if ~isempty(LcSup)
    CacolsSup= eval(typeCSup);
    plot(afig,CacolsSup(:,d1),CacolsSup(:,d2),'LineStyle','none','Marker',SymbolColsSup ,...
        'Color', ColorColsSup , 'MarkerFaceColor', MarkerFaceColorColsSup, 'MarkerSize',MarkerSize)
    text(CacolsSup(:,d1)+addx , CacolsSup(:,d2), LcSup,'Interpreter','None','FontSize',FontSize,'Color', ColorColsSup )
end

% set the x and y axis
xlimx=options.xlimx;
ylimy=options.ylimy;
if ~isempty(xlimx)
    xlim(xlimx);
end
if ~isempty(ylimy)
    ylim(ylimy);
end

Narray=out.N;
n=out.n;
I=out.I;
J=out.J;

seq1I=1:I;
seq1J=1:J;

confellipse=options.confellipse;
if isstruct(confellipse) || confellipse ==1
    
    if isstruct(confellipse)
        confellipsedef = struct;
        confellipsedef.conflev=0.90;
        confellipsedef.method={'multinomial'};
        confellipsedef.nsimul=1000;
        confellipsedef.selRows=seq1I;
        confellipsedef.selCols=seq1J;
        confellipsedef.AxesEllipse=true;
        
        fld=fieldnames(confellipse);
        
        % Check if user options inside plots are valid options
        chkoptions(confellipsedef,fld)
        
        % confellipse.nsimul=scalar which defines the number of
        % contingency tables which have to ge generated.
        if isfield(confellipse,'nsimul')
            nsimul=confellipse.nsimul;
        else
            nsimul=1000;
        end
        
        % confellipse.method= cell which specifies the method which is used to
        % compute confidence ellipses.
        if isfield(confellipse,'method')
            method=confellipse.method;
            % Check that the supplied method is a cell
            if ~iscell(method)
                error('FSDA:CorAnaplot:WrongInputOpt','method of collipse.mthod must be a cell')
            end
        else
            method={'multinomial'};
        end
        
        % confellipse.conflev= number in the interval (0 1) which
        % defines the confidence level of each ellipse.
        if isfield(confellipse,'conflev')
            conflev=confellipse.conflev;
        else
            conflev=0.90;
        end
        
        % Find if it is necessary to shows the major axes of the ellipses
        if isfield(confellipse,'AxesEllipse')
            AxesEllipse= confellipse.AxesEllipse;
        else
            AxesEllipse=true;
        end
        
        % Find the rows and columns for which ellipses are requested
        if isfield(confellipse,'selRows')
            selRows=confellipse.selRows;
            if iscell(selRows)
                [~,selRows]=intersect(out.Lr,selRows);
            elseif islogical(selRows)
                selRows=seq1I(selRows);
            else
                % Make sure that the indexes are in the interval [1 I]
                if ~isempty(selRows) && (max(selRows)>I || min(selRows)<1)
                    error('FSDA:CorAnaplot:WrongInputOpt',['The indexes of the requested ellipses must be integers in the interval [1,' num2str(I)])
                end
            end
        else
            selRows=seq1I;
        end
        
        if isfield(confellipse,'selCols')
            selCols=confellipse.selCols;
            if iscell(selCols)
                [~,selCols]=intersect(out.Lc,selCols);
            elseif islogical(selCols)
                selCols=seq1J(selCols);
            else
                % Make sure that the indexes are in the interval [1 I]
                if ~isempty(selCols) &&  max(selCols)>J  || min(selCols)<1
                    error('FSDA:CorAnaplot:WrongInputOpt',['The indexes of the requested ellipses must be integers in the interval [1,' num2str(J)])
                end
            end
        else
            selCols=seq1J;
        end
        
    else
        conflev=0.90;
        method={'multinomial'};
        nsimul=1000;
        selRows=seq1I;
        selCols=seq1J;
        AxesEllipse=true;
    end
    
    % Compute matrices Arows and Acols
    % Arows = [ A_1
    %           A_2
    %           ...
    %           A_{nsimul}]
    % Acols = [ A_1 A_2 ...  A_{nsimul}]
    % where matrix A_r, r=1, 2, ..., nsimul is the I-by-J matrix which is
    % associated to the r-th simulated replicate of the original
    % contingency table
    
    otherR=127;
    otherC=110;
    
    s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz';
    numRands = length(s);
    %specify length of random string to generate
    sLength = 15;
    %generate random string
    randString = s( ceil(rand(nsimul*max(I,J),sLength)*numRands) );
    randString=cellstr(randString);
    Lrsup=randString(1:nsimul*I);
    Lcsup=randString(1:nsimul*J);
    
    % Old inefficient code to generate labels for supplementary rows and
    % supplementary columns
    %         % default labels for rows of contingency table
    %     Lrsup=cellstr(strcat('rsup',num2str((1:nsimul*I)')));
    %     % default labels for columns of contingency table
    %     Lcsup=cellstr(strcat('csup',num2str((1:nsimul*J)')));
    %     Lrsup=strrep(Lrsup,' ','');
    %     Lcsup=strrep(Lcsup,' ','');
    
    % Find method to use to draw the confidence ellipses
    % methodcount is a counter which counts the number of requested methods
    % for which a matching has been found
    methodcount=0;
    legall=cell(3,2);
    selmethods=false(3,2);
    
    if max(strcmp(method,{'multinomial'})) ==1
        methodcount=methodcount+1;
        A=mnrnd(n,Narray(:)/n,nsimul);
        % The sum of each row of matrix A is n (total number of observations)
        Atransp=A';
        Acols=reshape(Atransp,I,J*nsimul);
        
        % A1row=reshape(Atmp(:),nr*nrep,nc);
        % Reshape A' into a 3D array
        Atmp=reshape(Atransp,I,J,nsimul);
        Arows = permute(Atmp,[1 3 2]);
        Arows = reshape(Arows,[],size(Atmp,2),1);
        Sup=struct;
        Sup.r=Arows;
        Sup.c=Acols;
        
        Sup.Lr=Lrsup;
        Sup.Lc=Lcsup;
        outext=CorAna(Narray,'Sup',Sup,'plots',0,'dispresults',false);
        
        if isstruct(plots) &&  isfield(plots,'alpha')
            if strcmp(plots.alpha,'rowgab')
                %  If plots.alpha='rowgab'  rows are in principal coordinates
                %  and columns are in standard coordinates multiplied by the
                %  mass.
                % typeR='RowsPri';        % rows are in principal coordinates
                ColsStaDcSup=bsxfun(@times,(sum(Acols,1)/n)',outext.ColsStaSup);
                outext.ColsStaDcSup=ColsStaDcSup;
                
                
            elseif strcmp(plots.alpha,'colgab')
                % If plots.alpha='colgab'  columns are in principal coordinates
                % and rows are in standard coordinates multiplied by the
                % masses.
                typeC='ColsPri';        % columns are in principal coordinates
                
                %RowsStaDrSup=diag(sum(Arows,2)/n)*outext.RowsStaSup;
                RowsStaDrSup=bsxfun(@times,sum(Arows,2)/n,outext.RowsStaSup);
                outext.RowsStaDrSup=RowsStaDrSup;
                
            elseif strcmp(plots.alpha,'rowgreen')
                %  If plots.alpha='rowgreen'  rows are in principal
                %  coordinates and columns are in standard coordinates
                %  multiplied by square root of the mass.
                typeR='RowsPri';        % rows are in principal coordinates
                ColsStaDcSqrtSup=bsxfun(@times,((sum(Acols,1)/n)').^(1/2),outext.ColsStaSup);
                outext.ColsStaDcSqrtSup=ColsStaDcSqrtSup;
                
            elseif strcmp(plots.alpha,'colgreen')
                %  If plots.alpha='colgreen' columns in principal coordinates
                %  and rows in standard coordinates multiplied by the square
                %  root of the mass.
                % typeC='ColsPri';        % columns are in principal coordinates
                
                % RowsStaDrSqrtSup=sqrt(diag(sum(Arows,2)/n))*outext.RowsStaSup;
                RowsStaDrSqrtSup=bsxfun(@times,sqrt(sum(Arows,2)/n),outext.RowsStaSup);
                outext.RowsStaDrSqrtSup=RowsStaDrSqrtSup;
                
            elseif isnumeric(plots.alpha) && plots.alpha>=0 && plots.alpha<=1
                hm=bsxfun(@rdivide,Arows,sum(Arows,2));
                outext.RowsAlphaSup=hm*outext.ColsSta *Gam^(plots.alpha-1);
                hm=bsxfun(@rdivide,Acols,sum(Acols,1));
                outext.ColsAlphaSup= hm'*outext.RowsSta* Gam^(-plots.alpha);
            else
            end
        end
        
        Carowsext= eval(['outext.' typeR 'Sup']);
        Cacolsext= eval(['outext.'  typeC 'Sup']);
        
        hold('on')
        if ~isempty(selCols)
            for j=selCols
                se=j:J:(J*(nsimul-1)+j);
                EcoCols=Cacolsext(se,dim);
                me=nanmean(EcoCols);
                co=nancov(EcoCols);
                [~,hColsMultinomial]=ellipse(me,co,conflev,[255 0 0]/255,AxesEllipse);
            end
            legall{1,2}=hColsMultinomial;
            selmethods(1,2)=true;
        end
        
        if ~isempty(selRows)
            for i=selRows
                se=i:I:(I*(nsimul-1)+i);
                EcoRows=Carowsext(se,dim);
                me=nanmean(EcoRows);
                co=nancov(EcoRows);
                [~,hRowsMultinomial]=ellipse(me,co,conflev,[0 0 255]/255,AxesEllipse);
            end
            legall{1,1}=hRowsMultinomial;
            selmethods(1,1)=true;
        end
    end
    
    % BOOTSTRAP METHOD bootRows
    if max(strcmp(method,{'bootRows'}))==1
        methodcount=methodcount+1;
        Acols=zeros(I,J*nsimul);
        nidot=sum(Narray,2);
        for i=1:I
            % samp = random sample with replacement of ni(i)*nrep elements from ni(i)
            % elements
            samp=randsample(nidot(i),nidot(i)*nsimul,true);
            % reshape samp in order to have a 3D array with nidot(i) rows and nrep columns
            samp=reshape(samp,nidot(i),nsimul);
            % Setup intervals for the frequency distribution
            edges=(0:J)*0.1/(nidot(i))+[0 cumsum(Narray(i,:))];
            counts=histc(samp,edges(1:end));
            % First column of counts is
            counts=counts(1:end-1,:);
            Acols(i,:)=(counts(:))';
        end
        
        Atransp=reshape(Acols,I*J,nsimul);
        Atmp=reshape(Atransp,I,J,nsimul);
        Arows = permute(Atmp,[1 3 2]);
        Arows = reshape(Arows,[],size(Atmp,2),1);
        Sup=struct;
        Sup.r=Arows;
        Sup.c=Acols;
        
        Sup.Lr=Lrsup;
        Sup.Lc=Lcsup;
        outext=CorAna(Narray,'Sup',Sup,'plots',0,'dispresults',false); %#ok<NASGU>
        
        Carowsext= eval(['outext.' typeR 'Sup']);
        Cacolsext= eval(['outext.'  typeC 'Sup']);
        
        hold('on')
        if ~isempty(selCols)
            for j=selCols
                se=j:J:(J*(nsimul-1)+j);
                EcoCols=Cacolsext(se,dim);
                me=nanmean(EcoCols);
                co=nancov(EcoCols);
                [~,hColsBootRows]=ellipse(me,co,conflev,[255 otherC 0]/255,AxesEllipse);
            end
            legall{2,2}=hColsBootRows;
            selmethods(2,2)=true;
        end
        
        if ~isempty(selRows)
            for i=selRows
                se=i:I:(I*(nsimul-1)+i);
                EcoRows=Carowsext(se,dim);
                me=nanmean(EcoRows);
                co=nancov(EcoRows);
                [~,hRowsBootRows]=ellipse(me,co,conflev,[0 otherR 255]/255,AxesEllipse);
            end
            legall{2,1}=hRowsBootRows;
            selmethods(2,1)=true;
        end
    end
    
    % BOOTSTRAP METHOD bootCols
    if max(strcmp(method,{'bootCols'})) ==1
        methodcount=methodcount+1;
        Arows=zeros(I*nsimul,J);
        ndotj=sum(Narray,1);
        for j=1:J
            % samp = random sample with replacement of ni(i)*nrep elements from ni(i)
            % elements
            samp=randsample(ndotj(j),ndotj(j)*nsimul,true);
            % reshape samp in order to have a 3D array with nidot(i) rows and nrep columns
            samp=reshape(samp,ndotj(j),nsimul);
            % Setup intervals for the frequency distribution
            edges=((0:I)')*0.1/(ndotj(j))+[0;cumsum(Narray(:,j))];
            counts=histc(samp,edges(1:end));
            % First column of counts is
            counts=counts(1:end-1,:);
            Arows(:,j)=(counts(:))';
        end
        % Note that
        % sum(Arows(1:I,:),1) = sum(Narray,1)
        % sum(Arows(I+1:2*I,:),1) =  sum(Narray,1)
        % ...
        
        % Now from Arows compute Acols
        Arowst=Arows';
        
        % Create a 3D array with J rows I cols and nsimul replicates
        Acols=reshape(Arowst,J,I,nsimul);
        % For each slice of the 3D array compute the transpose
        % The reults if is I-by-J-by-nimsul 3D array
        Acols=permute(Acols,[2,1,3]);
        % Transform the 3D array into a 2D array with size I-by-J*simul
        Acols=reshape(Acols,I,J*nsimul);
        % Note that
        %  Acols(1:14,1:J) = Arows(1:I,1:J)
        % Acols(1:14,J+1:2*J) = Arows(I+1:2*I,1:J)
        
        
        Sup=struct;
        Sup.r=Arows;
        Sup.c=Acols;
        
        Sup.Lr =Lrsup;
        Sup.Lc=Lcsup;
        outext=CorAna(Narray,'Sup',Sup,'plots',0,'dispresults',false); %#ok<NASGU>
        
        Carowsext= eval(['outext.' typeR 'Sup']);
        Cacolsext= eval(['outext.'  typeC 'Sup']);
        
        hold('on')
        if ~isempty(selCols)
            for j=selCols
                se=j:J:(J*(nsimul-1)+j);
                EcoCols=Cacolsext(se,dim);
                me=nanmean(EcoCols);
                co=nancov(EcoCols);
                [~,hColsBootCols]=ellipse(me,co,conflev,[255 2*otherC 0]/255,AxesEllipse);
            end
            legall{3,2}=hColsBootCols;
            selmethods(3,2)=true;
        end
        
        if ~isempty(selRows)
            for i=selRows
                se=i:I:(I*(nsimul-1)+i);
                EcoRows=Carowsext(se,dim);
                me=nanmean(EcoRows);
                co=nancov(EcoRows);
                [~,hRowsBootCols]=ellipse(me,co,conflev,[0 2*otherR 255]/255,AxesEllipse);
            end
            legall{3,1}=hRowsBootCols;
            selmethods(3,1)=true;
        end
    end
    legstring={'rows MultinomialR', 'cols MultinomialC';
        'rows BootRows', 'cols BootRows';
        'rows BootCols' 'cols BootCols'};
    legallLEG=[legall{selmethods(:)}];
    
    if verLessThanFS(9.2)==0
        % hColsMultinomial hColsBootRows hColsBootCols
        legend(legallLEG,...
            legstring(selmethods(:)),'AutoUpdate','off')
    else
        legend(legallLEG,legstring(selmethods(:)))
    end
    
    if methodcount ==0
        warning('FSDA:CorAnaplot:WrongInputOpt','Valid methods not found in input cell confellipse.method')
        disp('Methods supplied are')
        disp(confellipse.method)
        disp('Possible strings are')
        listStrings={'multinomial', 'bootRows', 'bootCols'};
        disp(listStrings)
        error('FSDA:CorAnaplot:WrongInputOpt','Please use one of the above strings inside cell confellipse.method')
    elseif isstruct(confellipse) && methodcount < length(confellipse.method)
        warning('FSDA:CorAnaplot:WrongInputOpt','Some methods specified in input cell confellipse.method are not valid')
        listStrings={'multinomial', 'bootRows', 'bootCols'};
        disp('Valid methods are')
        disp(listStrings)
        disp('Supplied methods are')
        disp(confellipse.method)
    else
    end
    
    
    % Supply supplementary rows and units using table format
    %     A1supcol=array2table(Acols);
    %     sup.c=A1supcol;
    %     A1suprow=array2table(Arows);
    %     sup.r=A1suprow;
    %     Ntable=out.Ntable;
    %     outext=CorAna(Ntable,'Sup',sup,'plots',0);
end


if ~isempty(h)
    % Eventually send the CorAnaxplot into a different figure/subplot
    hfigh = get(h,'Parent');
    
    set(hfigh,'Name','Correspondene analysis plot','NumberTitle','off');
    set(h,'Tag','pl_subplot');
    copyobj(allchild(afig),h);
    pause(0.0000001);
    delete(hfig);
    
end

title(titl,'Interpreter','Latex');

% Labels for axes
xlab=['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',InertiaExplained(d1,3)*100),'%)'];
xlabel(xlab,'FontName', FontName, 'FontSize', FontSizeAxisLabels);
ylab=['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',InertiaExplained(d2,3)*100),'%)'];
ylabel(ylab,'FontName', FontName, 'FontSize', FontSizeAxisLabels);
% Make axes equal and add cartesian axes
axis(gca,'equal')
vv=axis;
line([vv(1);vv(2)],[0;0])
line([0;0],[vv(3);vv(4)])

end
%FScategory:VIS-Mult