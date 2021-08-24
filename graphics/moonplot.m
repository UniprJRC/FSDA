function moonplot(out,varargin)
%moonplot draws the Correspondence Analysis (CA) moonplot.
%
% moonplot is an altenative way of visualizing the output of correspondence
% analysis. Row points are plotted in the traditional way. Column points
% are plotted equidistant from the origin, with their directions from the
% origin as in traditional correspondence analysis plots, and the
% information traditionally communicated by the distance of the points to
% the origin instead communicated by the size of the fonts of the labels.
% The closer row points are together, the more similar their positioning.
% The further a row point is from the middle, the more different it is from
% the “average”. Row points with very large masses are often in the middle
% because they define the average. Attributes (column points) are shown on
% the perimeter. The larger the font size of the attribute, the greater the
% level of discrimination of row points on the attribute (note:
% discrimination does not mean importance). The closer a row point to an
% attribute, the greater the association between the two.
% Moonplots are a better way to visualize brand maps than standard
% correspondence analysis outputs, which are often difficult to read
% correctly. The Moonplot resolves the key interpretation issues of
% correspondence analysis and is usually a better alternative.
%
%<a href="matlab: docsearchFS('moonplot')">Link to the help function</a>
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
%                   Note that 'colprincipal' can also be
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
%               If $plots.alpha='bothstandard'$, both the rows and columns
%                   are depicted in standard coordinates.
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
%                   labels. The default value is 10.
%              plots.FontSizeSup = scalar which specifies the font size of row
%                   labels of supplementary points. The default value is 10.
%              plots.MarkerSize = scalar which specifies the marker size
%                   of symbols associated with rows. The default
%                   value is 10.
%              plots.SymbolRows= character which specifies the symbol to
%                   use for row points. If this field is not present the
%                   default symbols is 'o'.
%              plots.SymbolRowsSup= character which specifies the symbol to
%                   use for supplementary row points. If this field is not
%                   present the default symbols is 'o'.
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
%              plots.MarkerFaceColorRowsSup= character which specifies the
%                   marker fill color to use for supplementary row points or RGB
%                   triplet. If this field is not present the default color
%                   is 'auto'.
%              Example - 'plots',plots=struct; plots.colorcols='k'
%              Data Types - double
%
%       addx : horizontal displacement for labels. Scalar. Amount of
%              horizontal displacement which has been put on the labels in the
%              plot. The defalut value of addx is 0.04.
%              Example - 'addx',0.01
%              Data Types - double
%
%       addy : vertical displacement for labels. Scalar. Amount of
%              vertical displacement which has been put on the labels in the
%              plot. The defalut value of addy is 0.
%              Example - 'addy',0.01
%              Data Types - double
%
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
%
%
%       xlimx   :   Min and Max of the x axis. Vector. Vector with two
%                   elements controlling minimum and maximum
%                   of the x axis.
%                   Example - 'xlimx',[-1 1]
%                   Data Types - double
%
%       ylimy   :   Min and Max of the y axis. Vector. Vector with two
%                   elements controlling minimum and
%                   maximum of the y axis.
%                   Example - 'ylimy',[0 1]
%                   Data Types - double
%
%        d1    :  Dimension to show on the horizontal axis. Positive
%                 integer. Positive integer in the range 1, 2, .., K which
%                 indicates the dimension to show on the x axis. The
%                 default value of d1 is 1.
%                 Example - 'd1',2
%                 Data Types - single | double
%
%        d2    :  Dimension to show on the vertical axis. Positive
%                 integer. Positive integer in the range 1, 2, .., K which
%                 indicates the dimension to show on the y axis. The
%                 default value of d2 is 2.
%                 Example - 'd2',3
%                 Data Types - single | double
%
%           h : the axis handle of a figure where to send the moonplot.
%               This can be used to host the moonplot in a subplot of a
%               complex figure formed by different panels (for example a panel
%               with moonplot from plots.alpha=0.2 and another
%               with moonplot from plots.alpha=0.5).
%               Example -'h',h1 where h1=subplot(2,1,1)
%               Data Types - Axes object (supplied as a scalar)
%
% Output:
%
%
% See also: corAnaplot, corAna, mcdCorAna
%
% References:
% Benzecri, J.-P. (1992), "Correspondence Analysis Handbook", New-York,
% Dekker.
% Benzecri, J.-P. (1980), "L'analyse des donnees tome 2: l'analyse des
% correspondances", Paris, Bordas.
% Greenacre, M.J. (1993), "Correspondence Analysis in Practice", London,
% Academic Press.
% Bock, T. (2011), Improving the display of correspondence analysis using moon
% plots, "International Journal of Market Research", Vol. 53, pp. 307-326.
% Gabriel, K.R. and Odoroff, C. (1990), Biplots in biomedical research,
% "Statistics in Medicine", Vol. 9, pp. 469-485.
% Greenacre, M.J. (1993), Biplots in correspondence Analysis, "Journal of
% Applied Statistics", Vol. 20, pp. 251-269.
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('moonplot')">Link to the help function</a>
%
%$LastChangedDate:: 2017-11-17 15:01:40 #$: Date of the last commit

% Examples:
%


%{
    %% moonplot with all the default options.
    % Prepare the data.
    load('csdPerceptions')
    N=csdPerceptions;
    out=CorAna(N,'plots',0);
    moonplot(out);
%}

%{
    %% moonplot with option changedimsign.
    % Prepare the data.
    load('mobilephone')
    N=mobilePhone;
    out=CorAna(N,'plots',0);
    % Use of option changedimsign
    moonplot(out,'changedimsign',[true false])
%}


%% Beginning of code

% Initialization
% Default font size for labels of rows or colums to add to the plot
FontSizedef=10;
MarkerSizedef=10;

options=struct('plots',1,'xlimx','','ylimy','','changedimsign',[false false],...
    'addy',0,'addx',0.04,'d1',1,'d2',2,'h','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:moonplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    plotsdef.SymbolRowsSup='o';
    plotsdef.ColorRows='b';
    plotsdef.ColorRowsSup='b';
    plotsdef.ColorCols='r';
    plotsdef.ColorColsSup='r';
    plotsdef.MarkerFaceColorRows='auto';
    plotsdef.MarkerFaceColorRowsSup='b';
    
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
            titl={'Mooplot starting from rows principal coordinates, and columns standard coordinates'};
            typeRSup='RowsPriSup';
            typeCSup='ColsStaSup';
            
        elseif strcmp(plots.alpha,'colprincipal')
            typeR='RowsSta'; % rows are in standard coordinates
            typeC='ColsPri'; % columns are in principal coordinates
            titl={'Mooplot starting from rows standard coordinates, and columns principal coordinates'};
            typeRSup='RowsStaSup';
            typeCSup='ColsPriSup';
            
            
        elseif strcmp(plots.alpha,'symbiplot')
            % equivalent to alpha=0.5
            typeR='RowsSym';        % rows are in symmetrical coordinates
            typeC='ColsSym';        % columns are in symmetrical coordinates
            titl='Moonplot starting from rows and columns in symmetrical coordinates';
            typeRSup='RowsSymSup';
            typeCSup='ColsSymSup';
            
            
        elseif strcmp(plots.alpha,'bothprincipal')
            typeR='RowsPri';        % rows are in principal coordinates
            typeC='ColsPri';        % columns are in principal coordinates
            titl={'Moonplot starting from rows and cols in principal coordinates.'};
            typeRSup='RowsPriSup';
            typeCSup='ColsPriSup';
            
        elseif strcmp(plots.alpha,'bothstandard')
            typeR='RowsSta';        % rows are in standard coordinates
            typeC='ColsSta';        % columns are in standard coordinates
            titl={'Moonplot starting from rows and cols in standard coordinates.'};
            typeRSup='RowsStaSup';
            typeCSup='ColsStaSup';
            
        elseif strcmp(plots.alpha,'rowgab')
            %  If plots.alpha='rowgab'  rows are in principal coordinates
            %  and columns are in standard coordinates multiplied by the
            %  mass.
            typeR='RowsPri';        % rows are in principal coordinates
            ColsStaDc=Dc*ColsSta;
            typeC=ExtractVariableName(ColsStaDc);
            titl={'Moonplot starting from rows principal coordinates, and columns standard coordinates times masses'};
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
            titl={'Moonplot starting from rows standard coordinates multiplied by masses ' , ...
                'and columns principal coordinates'};
            
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
            titl={'Mooplot starting from rows principal coordinates, and column standard coordinates ' , ...
                'times sqrt of masses'};
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
            titl={'Moonplot starting from rows standard coordinates times sqrt of masses,' ...
                'and columns principal coordinates'};
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
                    error('FSDA:moonplot:WrongInputOpt','Value of plots.alpha must lie in the interval [0 1]')
                end
            else
                listStrings={'rowprincipal'; 'colprincipal'; 'symbiplot'; 'bothstandard'; 'bothprincipal'; 'rowgab'; 'colgab'; 'rowgreen'; 'colgreen'};
                warning('FSDA:moonplot:WrongInputOpt',['Input string ''' plots.alpha ''' is  not found'])
                disp('Possible strings are')
                disp(listStrings)
                error('FSDA:moonplot:WrongInputOpt','Please use one of the above strings')
            end
        end
        
    else
        typeR='RowsSta';        % rows are in standard coordinates
        typeC='ColsSta';        % columns are in standard coordinates
        titl={'Moonplot starting from rows and cols in standard coordinates'};
        typeRSup='RowsStaSup';
        typeCSup='ColsStaSup';
    end
    
    
    if isfield(plots,'FontSizeSup')
        FontSizeSup=plots.FontSizeSup;
    else
        FontSizeSup=FontSizedef;
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
    
    
    if isfield(plots,'ColorRows')
        ColorRows=plots.ColorRows;
    else
        ColorRows='b';
    end
    
    if isfield(plots,'ColorRowsSup')
        ColorRowsSup=plots.ColorRowsSup;
    else
        ColorRowsSup='b';
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
    
    
else
    
    typeR='RowsPri';        % rows are in principal coordinates
    typeC='ColsPri';        % columns are in principal coordinates
    typeRSup='RowsPriSup';
    typeCSup='ColsPriSup';
    
    titl={'Moonplot starting from rows and cols in principal coordinates'};
    
    
    %     typeR='RowsPri';        % rows are in principal coordinates
    %     typeC='ColsPri';        % columns are in principal coordinates
    %     titl={'French symmetrical model: rows and cols in principal coordinates.'...
    %         'Plot of $X=D_r^{-1/2}U \Gamma$ and $Y= D_r^{-1/2} V \Gamma$'};
    %     typeRSup='RowsPriSup';
    %     typeCSup='ColsPriSup';
    
    FontSize=FontSizedef;
    FontSizeSup=FontSize;
    MarkerSize=MarkerSizedef;
    SymbolRows='o';
    SymbolRowsSup='o';
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
    % Marker fill color for supplementary rows
    MarkerFaceColorRowsSup=ColorRowsSup;
end

Carows= eval(typeR);
Cacols= eval(typeC);

% Carows=out.RowsSta;
% Cacols=out.ColsSta;

% Create the figure that will host the moonplot
hfig = figure('Name', 'Correspondence analysis plot', 'NumberTitle', 'off',...
    'Tag','pl_CorAna','Position', [100 100 600 600]);

% Get figure's axis
afig = axes('Parent',hfig);

% Add labels for row points and column points
% addx = adds a small horizontal displacement for labels
addx=options.addx;
% addy = adds a small vertical displacement for labels
addy=options.addy;

sum2=sum(Carows(:,dim).^2,2);
corrRows=(max(sum2))*1.10;
CarowsST=Carows(:,dim)/sqrt(corrRows);
theta=0:0.01:2*pi;
hold('on')
plot(sin(theta),cos(theta))

% Plot row points
plot(afig,CarowsST(:,d1),CarowsST(:,d2),'LineStyle','none','Marker',SymbolRows ,'Color', ColorRows , ...
    'MarkerSize',MarkerSize,'MarkerFaceColor', MarkerFaceColorRows)

% xlim([-1 1])
% ylim([-1 1])
% axis('square')
axis('off')
pbaspect([1 1 1])


% addx=0.02;
% addy=0.02;
text(CarowsST(:,1)+addx,CarowsST(:,2)+addy,Lr,'Color',ColorRows)



% Add points and text associated to supplementary rows
if ~isempty(LrSup)
    CarowsSup= eval(typeRSup);
    CarowsSupST=CarowsSup/sqrt(corrRows);
    plot(afig,CarowsSupST(:,d1),CarowsSupST(:,d2),'LineStyle','none','Marker',SymbolRowsSup ,...
        'Color', ColorRowsSup , 'MarkerFaceColor', MarkerFaceColorRowsSup,'MarkerSize',MarkerSize)
    text(CarowsSupST(:,d1)+addx , CarowsSupST(:,d2)+addy, LrSup,'Interpreter','None','FontSize',FontSizeSup,'Color', ColorRowsSup )
end

% Column points (labels around the circle)

ds=sqrt(sum(Cacols(:,dim).^2,2));
Fsize01=ds/max(ds);
% FontSize in the interval [7 18]
Fsize=7+Fsize01*11;
% The sum of squares of the rows of matrix CacolsST is 1
CacolsST=Cacols(:,dim)./ds;

rotat=atan2d(Cacols(:,dim(2)),Cacols(:,dim(1)));
for i=1:out.J
    if  rotat(i)>-90 && rotat(i)<90
        t=text(CacolsST(i,1),CacolsST(i,2),Lc(i),'FontSize',Fsize(i),'HorizontalAlignment','left','Color',ColorCols);
        t.Rotation=rotat(i);
    else
        t=text(CacolsST(i,1),CacolsST(i,2),Lc(i),'FontSize',Fsize(i),'HorizontalAlignment','right','Color',ColorCols);
        t.Rotation=rotat(i)-180;
    end
end

% Add points and text associated to supplementary columns
if ~isempty(LcSup)
    CacolsSup= eval(typeCSup);
    
    dsSup=sqrt(sum(CacolsSup(:,dim).^2,2));
    Fsize01Sup=dsSup/max(ds);
    % FontSize in the interval [7 18]
    FsizeSup=7+Fsize01Sup*11;
    % The sum of squares of the rows of matrix CacolsST is 1
    CacolsSupST=CacolsSup(:,dim)./ds;
    
    rotatSup=atan2d(CacolsSup(:,dim(2)),CacolsSup(:,dim(1)));
    
    for i=1:length(LcSup)
        if  rotatSup(i)>-90 && rotatSup(i)<90
            t=text(CacolsSupSupST(i,1),CacolsSupST(i,2),LcSup(i),'FontSize',FsizeSup(i),'HorizontalAlignment','left','Color',ColorColsSup);
            t.Rotation=rotat(i);
        else
            t=text(CacolsSupST(i,1),CacolsSupST(i,2),LcSup(i),'FontSize',FsizeSup(i),'HorizontalAlignment','right','Color',ColorColsSup);
            t.Rotation=rotat(i)-180;
        end
    end
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


if ~isempty(h)
    % Eventually send the moonplot into a different figure/subplot
    hfigh = get(h,'Parent');
    
    set(hfigh,'Name','Moonplot','NumberTitle','off');
    set(h,'Tag','pl_subplot');
    copyobj(allchild(afig),h);
    pause(0.0000001);
    delete(hfig);
    
end



title(titl,'Interpreter','Latex');


% Labels for axes
% xlab=['Dimension ',sprintf('%2.0f',d1),' (',sprintf('%5.1f',InertiaExplained(d1,3)*100),'%)'];
% xlabel(xlab,'FontName', FontName, 'FontSize', FontSizeAxisLabels);
% ylab=['Dimension ',sprintf('%2.0f',d2),' (',sprintf('%5.1f',InertiaExplained(d2,3)*100),'%)'];
% ylabel(ylab,'FontName', FontName, 'FontSize', FontSizeAxisLabels);
% % Make axes equal and add cartesian axes
% axis(gca,'equal')
vv=axis;
line([vv(1);vv(2)],[0;0],'LineStyle','--')
line([0;0],[vv(3);vv(4)],'LineStyle','--')

% Add cross in the origin of the axes
text(0,0,'+','FontSize',20,'HorizontalAlignment','center')

end
%FScategory:VIS-Mult