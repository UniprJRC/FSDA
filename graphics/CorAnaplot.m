function CorAnaplot(out,varargin)
%CorAnaplot draws the Correspondence Analysis (CA) graphs.
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
%              plots.dim = vector with two elements which specifies which
%                   dimensions to show in the factor map. The default is to
%                   show the first two dimensions, therefore plots.dim=[1 2]
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
%               confellipse.method= character which specifies the method which is used to
%                   compute confidence ellipses. Possible values are:
%                   'multinomial' = in this case the original contigencey
%                   table with the active elements is taken as a reference.
%                   Then new data tables are drawn in the following way:
%                   $r\times c$ values are drawn from a multinomial
%                   distribution with theoretical frequencies equals to
%                   $n_{ij}/n$.
%                   'bootrows = the values are bootstrapped row by row:
%                   Given row i, $n_{i.}$ are extracted with repetition and
%                   a frequencey distribution is computed using classes
%                   $[0, n_{i1}]$,$[n_{i1}, n_{i1}+n_{i2}]$, $\ldots$
%                   $[\sum_{j=1}^{J-1} n_{ij}, \sum_{j=1}^{J} n_{ij}$.
%                   'bootcols = the values are bootstrapped column by
%                   column.
%               confellipse.nsimul=scalar which defines the number of
%                   contingency tables which have to ge generated. The
%                   default value of confellipse.nsimul is 1000. Thus
%                   nsimul new contingency tables are projected as
%                   supplementary rows and/or supplementary columns.
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
%
%
% Output:
%
%   plotopt : options which have been used to create the plot. Cell array
%               of strings. Store all options which have been used to
%               generate the plot inside cell plotopt.
%
% See also: corAna
%
% References:
%
%   Atkinson and Riani (2000), Robust Diagnostic Regression Analysis,
%   Springer Verlag, New York.
%
% Copyright 2008-2017.
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
    % rows and columns are show in princicipal coordinates
    CorAnaplot(out)
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


%% Initialization
% Default font size for labels of rows or colums to add to the plot
FontSizedef=10;
MarkerSizedef=10;

options=struct('plots',1,'xlimx','','ylimy','','changedimsign',[false false],...
    'addx',0.04);

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

% Find the dimensions to extract
if isstruct(plots)
    if isfield(plots,'dim')
        d1=plots.dim(1);
        d2=plots.dim(2);
    else
        d1=1;
        d2=2;
    end
else
    d1=1;
    d2=2;
end


Lr=out.Lr;
Lc=out.Lc;
LrSup=out.LrSup;
LcSup=out.LcSup;

InertiaExplained=out.InertiaExplained;
Gam=diag(InertiaExplained(:,1));
Dr=out.Dr;
Dc=out.Dc;
n=out.n;

FontName='Times';
FontSizeAxisLabels=12;

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
            
            ColsStaDcSqrtSup=(diag(sum(SupColsN,1)/n))^(1/2)*ColsStaSup;
            typeCSup= ExtractVariableName(ColsStaDcSqrtSup);
            
            
        elseif strcmp(plots.alpha,'colgreen')
            %  If plots.alpha='colgreen' columns in principal coordinates
            %  and rows in standard coordinates multiplied by the square
            %  root of the mass.
            RowsStaDrSqrt=(sqrt(Dr))*RowsSta;
            typeR=ExtractVariableName(RowsStaDrSqrt);
            typeC='ColsPri';        % columns are in principal coordinates
            titl={'Rows standard coordinates times sqrt of masses,' ...
                'and columns principal coordinates, $X=U $ and $G= D_c^{-1/2} V \Gamma$'};
            RowsStaDrSqrtSup=sqrt(diag(sum(SupRowsN,2)/n))*RowsStaSup;
            
            typeRSup=ExtractVariableName(RowsStaDrSqrtSup);
            typeCSup='ColsPriSup';        % columns are in principal coordinates
            
        else
            if isnumeric(plots.alpha)
                if plots.alpha>=0 && plots.alpha<=1
                    RowsAlpha= RowsSta*Gam^plots.alpha;
                    ColsAlpha= ColsSta*Gam^(1-plots.alpha);
                    typeR=ExtractVariableName(RowsAlpha);
                    typeC=ExtractVariableName(ColsAlpha);
                    titl=['$\alpha='  num2str(plots.alpha) '\qquad   X=D_r^{-1/2} U \Gamma^{' num2str(plots.alpha) '}$ and $Y= D_c^{-1/2} V \Gamma^{1-'  num2str(plots.alpha) '}$'];
                    
                    RowsAlphaSup= RowsPriSup *Gam^(1-plots.alpha);
                    ColsAlphaSup= ColsPriSup* Gam^(-plots.alpha);
                    typeRSup=ExtractVariableName(RowsAlphaSup);
                    typeCSup=ExtractVariableName(ColsAlphaSup);
                    
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

figure
hold('on')
% Plot row points
plot(Carows(:,d1),Carows(:,d2),'LineStyle','none','Marker',SymbolRows ,'Color', ColorRows , ...
    'MarkerSize',MarkerSize,'MarkerFaceColor', MarkerFaceColorRows)

% Plot column points
plot(Cacols(:,d1),Cacols(:,d2),'LineStyle','none','Marker',SymbolCols ,'Color', ColorCols , ...
    'MarkerSize',MarkerSize,'MarkerFaceColor', MarkerFaceColorCols)

% Add labels for row points and column points
% addx = adds a small right horizontal displacement for labels
addx=options.addx;

text(Carows(:,d1)+addx , Carows(:,d2), Lr,'Interpreter','None','FontSize',FontSize,'Color', ColorRows )
text(Cacols(:,d1)+addx , Cacols(:,d2), Lc,'Interpreter','None','FontSize',FontSize,'Color', ColorCols )

title(titl,'Interpreter','Latex');

% Labels for axes
xlabel(['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',InertiaExplained(d1,3)*100),'%)'],'FontName', FontName, 'FontSize', FontSizeAxisLabels);
ylabel(['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',InertiaExplained(d2,3)*100),'%)'],'FontName', FontName, 'FontSize', FontSizeAxisLabels);

% Add points and text associated to supplementary rows
if ~isempty(LrSup)
    CarowsSup= eval(typeRSup);
    plot(CarowsSup(:,d1),CarowsSup(:,d2),'LineStyle','none','Marker',SymbolRowsSup ,...
        'Color', ColorRowsSup , 'MarkerFaceColor', MarkerFaceColorRowsSup,'MarkerSize',MarkerSize)
    text(CarowsSup(:,d1)+addx , CarowsSup(:,d2), LrSup,'Interpreter','None','FontSize',FontSize,'Color', ColorRowsSup )
end

% Add points and text associated to supplementary columns
if ~isempty(LcSup)
    CacolsSup= eval(typeCSup);
    plot(CacolsSup(:,d1),CacolsSup(:,d2),'LineStyle','none','Marker',SymbolColsSup ,...
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

% Make axis equal and add cartesian axes
axis(gca,'equal')
axis(gca,'equal')
vv=axis;
line([vv(1);vv(2)],[0;0])
line([0;0],[vv(3);vv(4)])

end
