function [fre]=unibiv(Y,varargin)
%unibiv has the purpose of detecting univariate and bivariate outliers
%
%<a href="matlab: docsearchFS('unibiv')">Link to the help function</a>
%
% Required input arguments:
%
% Y :           Input data. Matrix. 
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
%
% Optional input arguments:
%
%           rf  :  It specifies the confidence
%                  level of the robust bivariate ellipses. Scalar. 0<rf<1.
%                  The default value is 0.95 that is the outer contour in presence of normality
%                  for each ellipse should leave outside 5% of the values.
%                 Example - 'rf',0.99 
%                 Data Types - double
%      robscale :   It specifies the indexes to use
%                   to compute the dispersion of each variable and the
%                   correlation among each pair of variables. Scalar.
%                   robscale=1 (default): the program uses the median correlation
%                   and the MAD as estimate of the dispersion of each variable; 
%                   robscale=2: the correlation coefficient among ranks is used
%                   (Spearman's rho) and the MAD as estimate of the dispersion
%                   of each variable; 
%                   robscale=3: the correlation coefficient is based on Kendall's tau b
%                   and the MAD as estimate of the dispersion of each
%                   variable; 
%                   robscale=4: tetracoric correlation coefficient is used and the MAD
%                   as estimate of the dispersion of each variable; 
%                   otherwise the correlation and the dispersion of the variables are
%                   computed using the traditional (non robust) formulae
%                   around the univariate medians.
%                 Example - 'robscale',2 
%                 Data Types - double
%         plots :   It specify whether it is necessary to produce a plot
%                   with univariate standardized boxplots on the
%                   main diagonal and bivariate confidence ellipses out of
%                   the main diagonal. Scalar. If plots is equal to 1 a plot
%                   which contains univariate standardized boxplots on the
%                   main diagonal and bivariate confidence ellipses out of
%                   the main diagonal is produced on the screen. If plots is
%                   <> 1 no plot is produced. As default no plot is
%                   produced.
%                 Example - 'plots',2 
%                 Data Types - double
%       textlab : It controls the labels in the plots. Scalar. If textlab=1 and
%                   plots=1 the labels associated
%                   to the units which are univariate outliers or which are
%                   outside the confidence levels of the contours are
%                   displayed on the screen.
%                 Example - 'textlab',0
%                 Data Types - double
%       tag     :   It identifies the handle of the plot which
%                   is about to be created. Character. The default is to use tag
%                   'pl_unibiv'. Notice that if the program finds a plot which
%                   has a tag equal to the one specified by the user, then
%                   the output of the new plot overwrites the existing one
%                   in the same window else a new window is created.
%                 Example - 'tag','new_tag'
%                 Data Types - char
%       madcoef :   Coefficient which is used to scale MAD
%                   coefficient to have a robust estimate of dispersion. Scalar. The
%                   default is 1.4815 so that 1.4815*MAD(N(0,1))=1. 
%                 Example - 'madcoef',2 
%                 Data Types - double
%                   Remark: if mad =median(y-median(y))=0 then the interquartile
%                   range is used. If also the interquartile range is 0
%                   than the MD (mean absolute deviation) is used.  In
%                   other words MD=mean(abs(y-mean(Y))
%
%
% Output:
%
%   fre  :  n x 4 matrix which contains details about the univariate and
%           bivariate outliers.
%           1st col = index of the units; 
%           2nd col = number of times unit has been declared
%           univariate outliers; 
%           3rd col = number of times unit has been declared
%           bivariate outlier; 
%           4th col = pseudo MD as sum of bivariate MD.
%
%
% See also: FSMmmd
%
% References:
%
%       Riani, M., Zani S. (1997). An iterative method for the detection of
%       multivariate outliers, Metron, vol. LV, pp. 101-117.
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('unibiv')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % unibiv with all default options.
    % Run this code to see the output shown in the help file
    n=500;
    p=5;
    randn('state', 123456);
    Y=randn(n,p);
    [out]=unibiv(Y);
%}

%{
    %% unibiv with optional arguments.
    % Stack loss data.
    Y=load('stack_loss.txt');
    % Show robust confidence ellipses
    out=unibiv(Y,'plots',1,'textlab',1);
%}

[n,v]=size(Y);

% Default confidence level for bivariate ellipses
rfdef=0.95;

options=struct('rf',rfdef,'plots',0,'textlab',0,...
    'tag','pl_unibiv','robscale',1,'madcoef',1.4815);

%% Input parameters checking

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:unibiv:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end


% Write in structure 'options' the options chosen by the user
if nargin > 1
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

rf=options.rf;
plo=options.plots;
textlab=options.textlab;
robscale=options.robscale;

if rf<=0 || rf>=1;
    error('FSDA:unibiv:Wrongrf','The confidence threshold must be greater than 0 and smaller than 1')
end;

if plo==1
    % Create a figure to host the plot or clear the existing one
    h=findobj('-depth',1,'tag',options.tag);
    if (~isempty(h))
        clf(h);
        figure(h)
        axes;
    else
        % create a new figure
        h=figure;
    end
    set(h,'tag',options.tag);
    set(h,'Name', 'Univariate boxplots and robust bivariate confidence ellipses', 'NumberTitle', 'off');
    hold('all');
end

% bivT contains the frequency distribution of biv. outliers
bivT=zeros(n,1);
% univT contains the frequency distribution of univ. outliers
univT=zeros(n,1);
% MDbiv= vector which contains the sum of the bivaraite Mahalanobis distances for each unit
MDbiv=bivT;
% madcoef=1;
madcoef=options.madcoef;

seq=(1:n)';

for il=1:v;      % il is linked to the rows
    for jl=il:v;    % jl is linked to columns
        
        if il==jl;
            
            % Ys = vector which contains standardized data
            y=Y(:,il);
            Ty=median(y);
            
            mady=mad(y,1);
            
            if mady>0
                Ys=(y-Ty)/(madcoef*mady);
            else
                iqry=iqr(y);
                if iqry >0
                    mady=1.3490*iqr(y)/0.6745;
                    Ys=(y-Ty)/(madcoef*mady);
                else
                    mady=1.2533*mad(y)/0.6745;
                    Ys=(y-Ty)/(madcoef*mady);
                end
            end
            
            
            if robscale>4;
                % Sy is the unrobust standard deviation of y
                Sy=sqrt((y-Ty)'*(y-Ty)/(n-1));
                Ys=(y-Ty)/Sy;
            end
            
            % datax x add a sequence to standardized data
            datax=[seq Ys];
            
            % quan = 1 x 3 vector which contins 1% quartile median and 3rd
            % quartile
            quan = quantile(Ys,[.25 .50 .75]);
            % di= interquartile difference
            di=quan(3)-quan(1);
            % uq=upper truncation point
            uq=quan(3)+1.5*di;
            % lq=lower truncation point
            lq=quan(1)-1.5*di;
            % outy is a (l+1) x 2 matrix. the first column contains
            % the indexes of the units declared univariate
            % outliers, the second columns gives the standardized
            % values of the outliers
            
            outy=datax(datax(:,2)>uq | datax(:,2)<lq  | abs(datax(:,2))>3 ,:);
            
            
            if ~isempty(outy)
                % Increase by 1 the frequencey distribution of
                % univariate outliers in vector univT
                univT(outy(:,1))=univT(outy(:,1))+1;
            end
            
            
            % plotting part
            % produce a vertical boxplot @
            if plo==1;
                axis(axis);
                hbuf=subplot(v,v,jl+(il-1)*v);
                abuf=get(hbuf,'Position');
                h=boxplot(Ys);
                set(hbuf,'Position',abuf);
                hold('on');
                curaxuni=gca;
                if textlab==1
                    h=findobj(h,'tag','Outliers');
                    outcooy=sort(get(h,'Ydata'))';
                    outcoox=get(h,'Xdata')';
                    outy=sortrows(outy,2);
                    text(0.1+outcoox,outcooy,num2str(outy(:,1)));
                end
                
            end
            % end of univariate part
            
        else
            
            % beginning of bivariate part
            
            % Tx is the coordinate of the  Median
            x=Y(:,jl);
            Tx=median(x);
            madx=mad(x,1);
            
            if madx>0
                Xs=(x-Tx)/(madcoef*madx);
            else
                if robscale==1
                    warning('FSDA:unibiv:Median0','Median is 0 therefore robut correlation is computed using ranks')
                    robscale=2;
                end
                
                iqrx=iqr(x);
                if iqrx>0
                    madx=1.3490*iqrx/0.6745;
                    Xs=(x-Tx)/(madcoef*madx);
                    warning('FSDA:unibiv:IqrUsed','Interquartile range is used to scale the data')
                else
                    warning('FSDA:unibiv:MadUsed','Mean absolute deviation is used to scale the data')
                    madx=1.2533*mad(x)/0.6745;
                    Xs=(x-Tx)/(madcoef*madx);
                end
            end
            
            
            if robscale==1;
                us=abs(Xs+Ys);
                vs=abs(Xs-Ys);
                mus=(median(us))^2;
                mvs=(median(vs))^2;
                r=(mus-mvs)/(mus+mvs);
                if isnan(r)
                    r=0;
                end
                
            elseif robscale==2;
                % r is computed using ranks
                r=corr(x,y,'type','Spearman');
            elseif robscale==3;
                % r is based on the linear correlation
                % between the "concordances" sign(x(i)-x(j))*sign(y(i)-y(j)), i<j, with
                % an adjustment for ties.  This is often referred to as Kendall's tau-b.
                r=corr(x,y,'type','Kendall');
            elseif robscale==4;
                % r is based on the tetracoric correlation
                r=sum(sign((x-Tx).*(y-Ty)))/n;
            else
                Sx=sqrt((x-Tx)'*(x-Tx)/(n-1));      % Sx is the unrobust standard deviation of x
                
                Xs=(x-Tx)/Sx;   % standardization of x
                Ys=(y-Ty)/Sy;   % standardization of y
                r=Xs'*Ys/(n-1); % r= unrobust correlation 
            end
            
            % Spearman's rho and Kendall's tau and tetracoric correlation
            % are discrete-valued statistics, and
            % their distributions have positive probability at 1 and -1.
            % If their value is equal to 1 or -1 then artificially put the value equal
            % to 0.999 or -0.999
            if r==1;
                r=0.999;
            elseif r==-1;
                r=-0.999;
            end
            
            % Now we calculate Mahalanobis distances with centroid defined
            % by medians
            E=sqrt((Xs.^2+Ys.^2-2*r*Xs.*Ys)/(1-r^2));
            
            % Create the vector of pseudoMahalanobis distances (based on the sum
            % of each bivariate projection)
            MDbiv=MDbiv+E;
            
            %  Em  is the median  of vector E
            Em=median(E);
            emax= (2*(n-1))/(n-2)* finv(rf,2,n-2);
            emax=sqrt(emax);
            
            
            % hinge=ellipse containing 50% of the values
            R1=Em*sqrt((1+r)/2);
            R2=Em*sqrt((1-r)/2);
            th=0:(2*pi/3600):2*pi;
            The1=R1*cos(th);
            The2=R2*sin(th);
            Xx=(The1+The2);
            Yy=(The1-The2);
            
            % fence=ellipe containing (1-\alpha)% of the values
            R1max=emax*sqrt((1+r)/2);
            R2max=emax*sqrt((1-r)/2);
            The1max=R1max*cos(th);
            The2max=R2max*sin(th);
            Xmax=(The1max+The2max);
            Ymax=(The1max-The2max);
            
            f=[Xs Ys];
            
            % rotation of the coordinates
            cost=emax^2;
            a=1/(cost*(1-r^2));
            c=a;
            b=-r*a;
            
            aut=eig([a b;b c]);
            aut=sort(aut,'descend');
            costh=b/sqrt((c-aut(2))^2+b^2);
            sinth=(aut(2)-c)/sqrt((c-aut(2))^2+b^2);
            
            % M is the orthogonal matrix which enables the rotation of the axes
            M=[costh sinth; -sinth costh];
            
            % fou = fires of the ellipse
            fuo=sqrt(abs(1/aut(1)-1/aut(2)));
            rot=[Xmax;Ymax];
            
            % new1 = 2 x n matrix which contains the coordinates of the rotated ellipse
            new1=M*rot;
            % xnew = n x 1 vector which contains x coord. of rotated points
            xnew=new1(1,:)';
            
            % ynew = n x 1 vector which contains y coord. of rotated points
            ynew=new1(2,:)';
            
            % new2 = 2 x n matrix which contains the coordinates of the rotated points
            new2=M*(f');
            
            
            % ch is the fixed distance of each point lying on the ellipse
            i=10;
            ch=sqrt((xnew(i)+fuo)^2+ynew(i)^2)+sqrt((xnew(i)-fuo)^2+ynew(i)^2);
            
            % biv is the 1 x n vector which contains potential bivariate outliers
            biv=zeros(1,n);
            
            for jk=1:n;
                if sqrt((new2(1,jk)+fuo)^2+new2(2,jk)^2)+sqrt((new2(1,jk)-fuo)^2+new2(2,jk)^2)>ch;
                    biv(1,jk)=biv(1,jk)+1;
                else
                end
            end
            % bivT contains cumulative distribution of bivariate outliers
            bivT=bivT+biv';
            
            % the following lines plot the hinge together with the fence
            if plo==1;
                axis(axis);
                subplot(v,v,jl+(il-1)*v);
                hold('on');
                plot(Xx,Yy,Xmax,Ymax);
                % Add the points to the plot
                plot(Xs,Ys,'o')
                
                curaxbiv=gca;
                yl=get(curaxbiv,'Ylim');
                set(curaxuni,'Ylim',yl);
                
                % if textlab =1 labels are plotted for the units outside
                % the confidence ellipses
                if textlab==1;
                    text(Xs(biv==1),Ys(biv==1),num2str(seq(biv==1)));
                end
            end
            
            
        end; % endif of il=jl
        
    end % endif of jl=1:v
    
end % endif of il=1:v

fre=[seq univT bivT MDbiv];

if plo==1;
    hold('off')
end
end

%FScategory:MULT-Multivariate