function out = boxplotb(Y,varargin)
%boxplotb computes a bivariate boxplot
%
%<a href="matlab: docsearchFS('boxplotb')">Link to the help function</a>
%
% Required input arguments:
%
%             Y : Observations. Matrix. n x 2 data matrix: n observations
%               and 2 variables. Rows of Y represent observations, and
%               columns represent variables.
%
% Optional input arguments:
%
%           coeff : expansion factor. Scalar.
%                  Coefficient which enables us to pass from
%                 a contour which contains 50% of the data (hinge) to a contour
%                 which contains a prespecified portion of the data.
%                 Table below (taken from Zani, Riani and Corbellini, 1998,
%                 CSDA) shows the coefficients which must be used to obtain
%                 a theoretical threshold of 75, 90, 95 or 99 per cent in
%                 presence of normally distributed data:
%
%                   confidence level 0.75 -> coefficient 0.43;
%                   confidence level 0.90 -> coefficient 0.83;
%                   confidence level 0.95 -> coefficient 1.13;
%                   confidence level 0.99 -> coefficient 1.68.
%
%                   Example - 'coeff',1.68
%                   Data Types - double
%                 Remark: The default value of coeff is 1.68, that is 99%
%                 confidence level contours are produced.
%
% strictlyinside: additional peeling. Scalar. If strictlyinside=1 an
%                 additional convex hull is done on the 50% hull in order
%                 to increase the robustness properties of the method. In
%                 fact there may in general be some loss of robustness in
%                 small samples due to the use of peeling, therefore if we
%                 suspect to be in presence of a considerable propotion of
%                 outliers it may be necessary to do an additional peeling.
%                 The default value of strictlyinside is 0.
%                   Example - 'strictlyinside',1
%                   Data Types - double
%
%       plots   : graphical output. missing value | scalar | structure.
%                 This options specifies whether it
%                 is necessary to produce the bivariate boxplot on the
%                 screen.
%                 If plots is a missing value or is a scalar equal to 0 no
%                 plot is produced.
%                 If plots is a scalar equal to 1 (default) the bivariate
%                 boxplot with the outliers labelled is produced.
%                 If plots is a structure it may contain the following fields:
%                    plots.ylim = vector with two elements controlling minimum and maximum
%                       on the y axis. Default value is '' (automatic
%                       scale).
%                    plots.xlim = vector with two elements controlling minimum and maximum
%                       on the x axis. Default value is '' (automatic
%                       scale).
%                    plots.labeladd = If this option is '1', the outliers in the
%                       spm are labelled with the unit row index. The
%                       default value is labeladd='1', i.e. the row numbers are
%                       added.
%                    plots.InnerColor = a three element vector which specifies the
%                       color in RGB format to fill the inner contour
%                       (hinge). The default value of InnerColor is
%                       InnerColor=[168/255 150/255 255/255].
%                    plots.OuterColor = a three element vector which specifies the
%                       color in RGB format to fill the outer contour
%                       (fence). The default value of OuterColor is
%                       OuterColor=[210/255 203/255 255/255].
%                   Example - 'plots',1
%                   Data Types - double
%
%        resolution : resolution to use. Scalar. Resolution which must be
%                     used to produce the inner and outer spline.
%                     The default value of resolution is 1000, that is the
%                     splines are plotted on the screen using
%                     1000-by-(number of vertices of the inner hull) points.
%                   Example - 'resolution',5000
%                   Data Types - double
%
% Remark:       The user should only give the input arguments that have to
%               change their default value.
%               The name of the input arguments needs to be followed by
%               their value. The order of the input arguments is of no
%               importance.
% Output:
%
%         out:   structure which contains the following fields
%
%
%  out.outliers = vector containing the list of the units which lie outside the
%             outer contour.
%             REMARK: if no unit lies outside the outer spline outliers is a
%             Empty matrix: 0-by-1
%
%         out.cent = 2 x 1 vector containing the coordinates
%                of the robust centroid.
%                cent[1] = x coordinate;
%                cent[2] = y coordinate.
%
%          out.Spl = r-by-4 matrix containing the coordinates
%                of the inner and outer spline. r (rows of matrix Spl) is
%                approximately equal to the number of vertices of the inner hull
%                multiplied by the resolution which is used.
%                The first two columns refer to the (x,y) coordinates of
%                the inner spline.
%                The last  two columns refer to the (x,y) coordinates  of the
%                outer spline.
%
%
% See also convhull.m, FSM.m
%
% References:
%
% Zani, S., Riani M. and Cerioli A. (1998), Robust bivariate boxplots
% and multiple outlier detection, "Computational Statistics and Data
% Analysis", Vol. 28, pp. 257-270.
% Corbellini A., Riani M. and Atkinson A.C. (2015), Discussion of the
% paper 'Multivariate Functional Outlier Detection' by Hubert, Rousseeuw
% and Segaert, "Statistical Methods and Applications".
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('boxplotb')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% boxplotb with all default options.
    % Bivariate boxplot of the writing data at time t=5.
    % This example reproduces Figure 1 of Corbellini, Riani and Atkinson,
    % 2015, Statistical Methods and Applications
    close all
    X=load('writingdata.txt');
    out=boxplotb(X);
    xlabel('horizontal coordinate')
    ylabel('vertical coordinate')
    title('Bivariate boxplot of the writing data at time $t=5$','Interpreter','Latex')
%}

%{
    %% boxplotb with optional arguments.
    % Bivariate boxplot of the stars data
    % This example reproduces Figure 4 of Zani Riani and Corbellini
    close all
    X=load('stars.txt');
    out=boxplotb(X,'strictlyinside',1);
    xlabel('Log effective surface temperature')
    ylabel('Log light intensity')
%}

%{
    % Bivariate boxplot of the brain data.
    % This example reproduces Figure 4 of Zani Riani and Corbellini
    close all
    X=load('bodybrain.txt');
    X=log10(X);
    out=boxplotb(X);
    xlabel('Log (to the base 10) body weight')
    ylabel('Log (to the base 10) brain weight')
    title('Bivariate boxplot of Log brain weight and Log body weight for 28 animals')
%}

%{
    % Bivariate boxplot of the stars data.
    % Now we change the colors of the inner and outer contour to white
    % In this example we explore the various graphical options
    close all
    X=load('stars.txt');
    plots=struct;
    plots.InnerColor=[0 0 0]+1; % remove the color for the hinge
    plots.OuterColor=[0 0 0]+1; % remove the color for the fence
    plots.labeladd=0; % do not include the labels for the outliers
    plots.xlim=[min(X(:,1)) max(X(:,1))];  % tight xlim
    plots.ylim=[min(X(:,2)) max(X(:,2))];  % tight ylim
    out=boxplotb(X,'strictlyinside',1,'plots',plots);
    xlabel('Log effective surface temperature')
    ylabel('Log light intensity')
%}

%{
    % Bivariate boxplot of two variables of Emilia Romagna data.
    % This example reproduces Figure 2 of Zani Riani and Corbellini
    close all
    load('emilia2001')
    Y=emilia2001{:,:};
    % Extract the variables y1 and y3
    % y1= Percentage of infant population (that is the percentage of
    % population aged less than 10)
    % y3 = % of single member (one component) families
    X=Y(:,[1 3]);
    % In order to reproduce exactly Figure 2 of Zani, Riani and Corbellini
    % (1998), CSDA, we remove municipalities with a percentage of single
    % members greater than 45%
    X=X(X(:,2)<45,:);
    out=boxplotb(X,'strictlyinside',1);
    xlabel('y1=Percentage of infant population')
    ylabel('y3 = Percentage of single member families')
%}

%% Beginning of code

% Input parameters checking
% Extract size of the data
n=size(Y,1);
% seq = sequence from 1 to n which enables us to identify observations
seq=(1:n)';

if nargin<1
    error('FSDA:boxplotv:missingInputs','Initial data matrix is missing');
end

resolution=1000;
coeff=1.68;
plo=1;
strictlyinside=0;

options=struct('resolution',resolution,'coeff',coeff,...
    'plots',plo,'strictlyinside',strictlyinside);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:boxplotb:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


if nargin > 1
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    coeff=options.coeff;
    resolution=options.resolution;
    plo=options.plots;
    strictlyinside=options.strictlyinside;
end

Ywithseq = [Y seq];

% Do convex hull peeling until the first hull which contains not more than
% 50% of the points is obtained
% Remark: n-size(Ywithseq,1) = length of the points which have been
% eliminated up to each step of the loop; The first time that the removed
% points is greater or equal than n/2 procedure stops. In other words, loop
% continues while n-size(Ywithseq,1)<(n/2)
while n-length(Ywithseq)<(n/2)
    % disp(n-length(Ywithseq))

    K = convhull(Ywithseq(:,1),Ywithseq(:,2));

    % sel = vector of true row coordinates of the units forming the current
    % convex hull
    sel=Ywithseq(K(1:end-1),3);

    % Rows of matrix Ywithseq which belong to the current convex hull
    % are set to zero
    Ywithseq(K(1:end-1),:)=0;

    % Redefine matrix Ywithseq deleting the rows which are equal to zero
    Ywithseq = Ywithseq(any(Ywithseq,2),:);
end

% if strictyinside==1 do an additional peeling
if strictlyinside==1
    K = convhull(Ywithseq(:,1),Ywithseq(:,2));

    % sel = vector of true row coordinates of the units forming the current
    % convex hull
    sel=Ywithseq(K(1:end-1),3);
end


% FiftyPerCentHull = largest convex hull which contains not more than 50%
% of the data. To this convex hull we must superimpose the spline curve
FiftyPerCentHull=Y(sel,:);

% Insplx,Insply are respectively the x and y coordinates of the B-spline
% which is superimposed on the 50% convex hull
[Insplx,Insply]=bsplineFS(FiftyPerCentHull(:,1),FiftyPerCentHull(:,2), resolution);
Inspl=[Insplx Insply];

% cent = robust centroid. Mean of the units inside the 50% hull
cent=mean(Ywithseq(:,1:2),1);

% Construct outer contour
% Ospl = l x 2 matrix containing the coordinates
%       of the outer spline @
%     Ospl=(Ispl-cent')*(coeff+1)+cent';
Ospl=bsxfun(@minus,Inspl,cent)*(1+coeff);
Ospl=bsxfun(@plus,Ospl,cent);

% Find the points which are outside the outer contour
[in]=inpolygon(Y(:,1),Y(:,2),Ospl(:,1),Ospl(:,2));
outliers=seq(in==0);

% outcor is the old function
% [outlierschk]=outcor(Y(:,1),Y(:,2),Ospl(:,1),Ospl(:,2),cent(1:2));

if isstruct(plo) || (~isstruct(plo) && plo~=0)
    if isstruct(plo)

        fplo=fieldnames(plo);

        d=find(strcmp('xlim',fplo));
        if d>0
            xlimx=plo.xlim;
        else
            xlimx='';
        end

        d=find(strcmp('ylim',fplo));
        if d>0
            ylimy=plo.ylim;
        else
            ylimy='';
        end


        d=find(strcmp('InnerColor',fplo));
        if d>0
            InnerColor=plo.InnerColor;
        else

            InnerColor=[168/255 150/255 255/255];
        end

        d=find(strcmp('OuterColor',fplo));
        if d>0
            OuterColor=plo.OuterColor;
        else
            OuterColor=[210/255 203/255 255/255];
        end

        d=find(strcmp('labeladd',fplo));
        if d>0
            labeladd=plo.labeladd;
        else
            labeladd='1';
        end


    else
        xlimx='';
        ylimy='';
        labeladd='1';
        InnerColor=[168/255 150/255 255/255];
        OuterColor=[210/255 203/255 255/255];
    end


    hold('on')

    % set RGB colors
    hh1=fill(Ospl(:,1),Ospl(:,2),OuterColor);
    hh1.Annotation.LegendInformation.IconDisplayStyle = 'off';

    hh2=fill(Inspl(:,1),Inspl(:,2),InnerColor);
    hh2.Annotation.LegendInformation.IconDisplayStyle = 'off';


    sel=setdiff(seq,outliers);
    % Plot the data excluding the outliers
    hh6=plot(Y(sel,1),Y(sel,2),'ko','MarkerFaceColor','k','Markersize',4);
    hh6.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % Plot the outliers
    hh7=plot(Y(outliers,1),Y(outliers,2),'*r');
    if ~isempty(hh7)
        hh7.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    
    % Plot coordinates of the 50% hull
    % plot(FiftyPerCentHull(:,1),FiftyPerCentHull(:,2),'k');

    % Plot inner spline defining the hinge (50% of the data)
    hh3=plot(Inspl(:,1),Inspl(:,2),'Color','k');
    hh3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Plot coordinates of the outer spline curve
    hh4=plot(Ospl(:,1),Ospl(:,2),'Color','k');
    hh4.Annotation.LegendInformation.IconDisplayStyle = 'off';

    if strcmp(labeladd,'1')
        % Plot text associated with the outliers
        if sum(isnan(outliers))==0
            text(Y(outliers,1),Y(outliers,2),cellstr(num2str(outliers)))
        end
    end

    % Plot robust centroid
    hh5=plot(cent(1),cent(2),'r+','MarkerSize',14,'LineWidth',3);
    hh5.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Set the x and y scale
    if ~isempty(xlimx)
        xlim(xlimx)
    end

    % Set the x and y scale
    if ~isempty(ylimy)
        ylim(ylimy)
    end

end



% Store robust centroid
out.cent=cent;

% Store list of units falling outside the outer contour
out.outliers=outliers;

% Store coordinates of inner and outer spline
out.Spl=[Inspl Ospl];
% fcns = localfunctions;
% display(fcns)

end




% local functions
function [X,Y]= bsplineFS(x,y, resolution)

% mesh points overlapping:
% needed to close the path of the bicubic spline
x=[x; x(1:3)];
y=[y; y(1:3)];

dims=length(x)*resolution;

X=zeros(dims,1);
Y=zeros(dims,1);

% spline resolution:
% number of spline segments between each consecutive
% mesh points, where mesh point are the vertexes of the convex hull
if nargin<3
    resolution=3000;
end
z=0;
n=length(x);

k=1;

for i = 2:n-2
    % x
    xa = x(i - 1);
    xb = x(i);
    xc = x(i + 1);
    xd = x(i + 2);
    % y
    ya = y(i - 1);
    yb = y(i);
    yc = y(i + 1);
    yd = y(i + 2);
    % spline coefs
    a3 = ( - xa + 3 * (xb - xc) + xd)/6;
    b3 = ( - ya + 3 * (yb - yc) + yd)/6;
    a2 = (xa - 2 * xb + xc)/2;
    b2 = (ya - 2 * yb + yc)/2;
    a1 = (xc - xa)/2;
    b1 = (yc - ya)/2;
    a0 = (xa + 4 * xb + xc)/6;
    b0 = (ya + 4 * yb + yc)/6;
    for j = 1:resolution
        t=j / resolution;
        z=k;

        X(k)=(((a3 * t + a2) * t + a1) * t + a0);
        Y(k)=(((b3 * t + b2) * t + b1) * t + b0);
        k=k+1;

    end
end
% X=X(2:z - 30);
% Y=Y(2:z - 30);

X=X(1:z);
Y=Y(1:z);

end

% % local function
% function [pointIndex]=outcor(x,y,splx,sply,cent)
%
% centrx=cent(1);
% centry=cent(2);
%
% % sensitivity
% eps=0.001;
% % points length
% nn=length(x);
% X =zeros(nn,1);
% % spline length
% ll = length(splx);
%
% % accumulation counters
% k = 0;
% n = 1;
%
%     xa = centrx;
%     ya = centry;
%
% while(n < nn)
%     xb = x(n);
%     yb = y(n);
%     j = 2;
%     while(j < ll)
%         xc = splx(j - 1);
%         xd = splx(j);
%         yc = sply(j - 1);
%         yd = sply(j);
%         den = (xb - xa) * (yd - yc) - (yb - ya) * (xd - xc);
%         r = ((ya - yc) * (xd - xc) - (xa - xc) * (yd - yc))/den;
%         s = ((ya - yc) * (xb - xa) - (xa - xc) * (yb - ya))/den;
%         if((s >= 0 && s <= 1) && (r >= 0 && r <= 1))
%             xk = xa + r * (xb - xa);
%             yk = ya + r * (yb - ya);
%             dp1 = sqrt((xk-xa)*(xk-xa)+(yk-ya)*(yk-ya));
%             dp2 = sqrt((xb-xa)*(xb-xa)+(yb-ya)*(yb-ya));
%             dp = sqrt((dp1-dp2)*(dp1-dp2));
%             if (dp > eps)
%                 % the point lies outside the spline curve
%                 X(n) = 1;
%                 k=k+1;
%             else
%                 % the point lies inside the spline curve
%                 X(n) = 0;
%                 k=k+1;
%             end
%         end
%         j=j+1;
%     end
%     n=n+1;
% end
% % find the indexes of outlying points
% pointIndex=find(X==1);
% end

%FScategory:VIS-Mult
