function [mmdinv] = FSMinvmmd(mmd,v,varargin)
%FSMinvmmd converts values of minimum Mahalanobis distance into confidence levels
%
%<a href="matlab: docsearchFS('FSMinvmmd')">Link to the help function</a>
%
%  Required input arguments:
%
% mmd :          Distances. Matrix. n-m0 x 2 matrix. 
%                1st col = fwd search index; 
%                2nd col = minimum Mahalanobis distance.
%                Data Types - single | double
% v :            Number of variables. Scalar. 
%                Number of variables of the underlying dataset. 
%                Data Types - single | double
%
%  Optional input arguments:
%
%        n:     It specifies the size of the sample. Scalar. If it is not specified
%               it is set equal to mmd(end,1)+1.
%                 Example - 'n',5
%                 Data Types - double
%   plots :     Plot on the screen. Scalar or structure.
%               If plots = 1, a plot which shows the
%               confidence level of mmd in each step is shown on the
%               screen. Three horizontal lines associated respectively with
%               values 0.01, 0.5 and 0.99  are added to the plot. 
%               If plots is a structure, it may contain the following fields: 
%                   plots.conflev = vector containing horizontal lines associated
%                       with confidence levels; 
%                   plots.conflevlab = scalar if it is equal 1 labels associated with
%                       horizontal lines are shown on the screen; 
%                   plots.xlim = minimum and maximum on the x axis; 
%                   plots.ylim = minimum and maximum on the y axis; 
%                   plots.LineWidth = Line width of the trajectory of mmd in
%                   normal coordinates; 
%                   plots.LineStyle = Line style of the
%                   trajectory of mle of transformation parameters; 
%                   plots.LineWidthEnv = Line width of the horizontal lines; 
%                   plots.Tag = tag of the plot (default is pl_mmdinv); 
%                   plots.FontSize = font size of the text labels which identify
%                   the trajectories
%                 Example - 'plots',1
%                 Data Types - double
%
%  Output:
%
%   mmdinv:     confidence levels plotted in normal coordinates. 
%               (n-m0) x 3 matrix (same rows of input matrix mmd).
%               It contains information about requested
%               confidence levels plotted in normal coordinates.
%               1st col = fwd search index from m0 to n-1; 
%               2nd col = confidence level of each value of mmd; 
%               3rd col = confidence level in normal coordinates. 
%                    50 per cent conf level becomes norminv(0.50)=0; 
%                    99 per cent conf level becomes norminv(0.99)=2.33.
%
%
% More About:
%
%Let $d^2_i(m )$ and  $d_{\mbox{min}}(m)$ be respectively the deletion
%distance for unit $i$ based
% on a subset of size $m$ and $d_{\mbox{min}}(m)$ the min. Mahalanobis
% distance in the forward search at step m. Testing for outliers requires a
% reference distribution for $d^2_i(m )$ in and hence for
% $d_{\mbox{min}}(m)$ in (\ref{min}). When $\Sigma$ is estimated from all
% $n$ observations, the squared statistics  have an $F$ distribution.
% However, the estimate $\hat{\Sigma}(m)$ in the search uses the central
% $m$ out of $n$ observations, so that the variability is underestimated.
% The consistency factor $c(m,n)$ given below
% \[
%  c(m,n)=\frac{n}{m} C_{v+2} \{C_{v}^{-1} (m/n) \}
% \]
% where $C_r$ is the c.d.f. of the $\chi^2$ distribution on $r$ degrees of
% freedom, allows for estimation from this truncated distribution,  providing an
%  approximately unbiased estimate of $\Sigma$.
% 
% We can treat the distribution of the rescaled deletion Mahalanobis distance
% $c(m,n)d_{\mbox{min}}^2(m)$ as a squared deletion distance on $m-1$
% degrees of freedom, whose distribution is (Atkinson Riani and Cerioli,
% 2004; pp. 43-44)
% \begin{equation}\label{F}
%   \frac{m^2-1}{m(m-v)} F_{v,m-v},
% \end{equation}
% 
% The distribution of the rescaled min Mahalanobis distance
%  $c(m,n) d_{\mbox{min}}^2(m)$
% of a subset of size $m$ constructed in such a way that
%  the centroid and covariance matrix of the subset are taken using the
%  units having the $m$ smallest Mahalanobis distances can be treated as
%  the distribution of the $(m+1)$th order statistic from ($F_{v,m-v}$).
% 
% The results of order statistics $Y_{(1)}$, $Y_{(2)}$, $\cdots$,  $Y_{(n)}$ from a sample of size $n$ from a distribution with CDF $G(y)$, state that
% \begin{equation}
% \label{orderstat}
%  P\{Y_{(m+1)} \le y \}  =  P \left\{ F_{2(n-m),2(m+1)} >  \frac{1-G(y)}{G(y)} \times
%  \frac{m+1}{n-m} \right\}
%  \end{equation}
% Given that in our case $G(y)$ is the CDF of the  $F_{v,m-v}$ we can rewrite this equation  as
% 
% \begin{eqnarray*}
%   &&  P\{d_{\mbox{ min}}^2(m)  \leq \widehat{ d_{\mbox{min}}^2(m)} \} = \\
%    && 1-  F_{2(n-m),2(m+1)} \left( \left( \frac{1}{  F_{v,m-v} \left( \frac{m(m-v)}{m^2-1 } c(m,n) d_{\mbox{min}}^2(m) \right)  }-1 \right) \frac{m+1}{n-m} \right)
% \end{eqnarray*}
% where $F_{a,b}(y)$ is the CDF of the $F$ distribution with $a$ and $b$ degrees of freedom evaluated in $y$.
% The value of the min. Mahalanobis distance transformed in normal coordinates computed by this routine is
% nothing but
% \[
% \Phi^{-1} \left(  P\left\{ d_{\mbox{min}}^2(m)  \leq \widehat{ d_{\mbox{min}}^2(m)} \right\}    \right)
% \]
% where $\Phi^{-1}$ is the inverse of the CDF of the standard normal distribution.
%
% See also FSMenvmmd, FSM.m, FSMeda.m
%
% References:
%
%   Atkinson, A.C. and Riani, M. (2006). Distribution theory and
%   simulations for tests of outliers in regression. Journal of
%   Computational and Graphical Statistics, Vol. 15, pp. 460-476
%   Riani, M. and Atkinson, A.C. (2007). Fast calibrations of the forward
%   search for testing multiple outliers in regression, Advances in Data
%   Analysis and Classification, Vol. 1, pp. 123-141.
%
% Copyright 2008-2015
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSMinvmmd')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
    % FSMinvmmd with all default options.
    % After creating 99 per cent confidence envelopes based on 1000
    % observations and 5 variables are created, their confidence level 
    % is calculated with FSMinvmmd.
      v=5;
      mmdenv=FSMenvmmd(1000,v,'prob',0.99);
      mmdinv=FSMinvmmd(mmdenv,v);
    % mmdinv is a matrix which in the second colum contains
    % all values equal to 0.99.
     
%}

%{
    % FSMinvmmd with optional arguments.
    % Example of finding confidence level of mmd. Forgery Swiss Banknotes data. 
    load('swiss_banknotes');
    Y=swiss_banknotes.data;
    Y=Y(101:200,:);
    % The line below shows the plot of mmd
    [out]=FSM(Y,'plots',1);

    % The line below transforms the values of mmd into observed confidence
    % levels and shows the output in a plot in normal coordinates using all
    % default options
    plots=struct;
    plots.conflev=[0.01 0.5 0.99 0.999 0.9999 0.99999];
    mmdinv=FSMinvmmd(out.mmd,size(Y,2),'plots',plots);
%}

%{
    %% Resuperimposing envelopes and normal coordinates.
    % Comparison of resuperimposing envelopes using mmd coordinates and normal
    % coordinates. Forgery Swiss Banknotes data. 
    load('swiss_banknotes');
    Y=swiss_banknotes.data;
    Y=Y(101:200,:);
    % The line below shows the plot of mmd
    [out]=FSM(Y,'plots',2);

    n0=83:86;
    quantplo=[0.01 0.5 0.99 0.999 0.9999 0.99999];
    ninv=norminv(quantplo);
    lwdenv=2;
    supn0=max(n0);

    ij=0;
    for jn0=n0;
        ij=ij+1;
        MMDinv = FSMinvmmd(out.mmd,size(Y,2),'n',jn0);
        % Resuperimposed envelope in normal coordinates
        subplot(2,2,ij)
        plot(MMDinv(:,1),MMDinv(:,3),'LineWidth',2)
        xlim([out.mmd(1,1) supn0])
        v=axis;
        line(v(1:2)',[ninv;ninv],'color','g','LineWidth',lwdenv,'LineStyle','--','Tag','env');
        text(v(1)*ones(length(quantplo),1),ninv',strcat(num2str(100*quantplo'),'%'));
        title(['Resuperimposed envelope n=' num2str(jn0)]);
    end
%}


%% Input parameters checks
nn=mmd(end,1)+1; % This line is introduced for publishFS
options=struct('n',nn,'plots','');

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:FSMinvmmd:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)
end

if nargin>2
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end

n=options.n;

if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSMinvmmd:Wrongn','n must be scalar non empty and non missing!!');
end

if ~isscalar(v) || isempty(n) || isnan(v)
    error('FSDA:FSMinvmmd:Wrongv','v must be scalar non empty and non missing!!!');
end


%% Find confidence level of each value of mmd

% mm = column vector which contains fwd search index.
mm =mmd(:,1);

% Compute variance of the truncated normal distribution.
% mm/n is the percentage of observations inside subset.
% corr=v*((mm+1)./mm).*(mm-1)./(mm-v);

a=chi2inv(mm/n,v);
corr=(n./mm).*(chi2cdf(a,v+2));

invmmd=(mmd(:,2).^2).*corr.*(mm./(mm+1)).*((mm-v)./(v*(mm-1)));

% regression qq=-mm-1+(mm+1)./(2*tcdf(mmd(:,2).*sqrt(corr),mm-v)-1);
qq=-mm-1+(mm+1)./(fcdf(invmmd,v,mm-v));

qq=qq./(n-mm);

% mmd= 1-confidence level
mmd=fcdf(qq,2*(n-mm),2*(mm+1));
% mmdt = confidence level
mmdt=1-mmd;

% Confidence level in normal coordinates
% One may wonder why mmdncoo cannot be computed as follows
% mmdncoo=norminv(mmt);
% The reason is that while norminv(1-\epsilon)=inf but
% -norminv(\epsilon) is exact where \epsilon is a very very small number
mmdncoord=-norminv(mmd);

mmdinv=[mm mmdt mmdncoord];

%% Plotting part

if ~isempty(options.plots)
    
    plots=options.plots;
    
    if isstruct(plots)
        
        fplots=fieldnames(plots);
        
        d=find(strcmp('xlim',fplots));
        if d>0
            xlimx=plots.xlim;
        else
            xlimx='';
        end
        
        
        d=find(strcmp('LineWidth',fplots));
        if d>0
            LineWidth=plots.LineWidth;
        else
            LineWidth=2;
        end
        
        % LineWidthEnv = line width of the horizontal lines associated
        % with required confidence envelopes
        % the plot of monitoring of MLE of transformation parameters or the
        % horizontal lines associated with the asymptotic confidence levels
        % for the likelihood ratio test
        d=find(strcmp('LineWidthEnv',fplots));
        if d>0
            LineWidthEnv=plots.LineWidthEnv;
        else
            LineWidthEnv=1;
        end
        
        % Specify the line type for the trajectory of mmd
        % in normal coordinate
        d=find(strcmp('LineStyle',fplots));
        if d>0
            LineStyle=plots.LineStyle;
        else
            LineStyle={'-'};
        end
        
        % Horizontal lines associated with confidence level
        d=find(strcmp('conflev',fplots));
        if d>0
            conflev=plots.conflev;
        else
            conflev=[0.01 0.5 0.99];
        end
        
        d=find(strcmp('ylim',fplots));
        if d>0
            ylimy=plots.ylim;
        else
            ylimy=[min([mmdncoord;norminv(min(conflev))]) max([mmdncoord;norminv(max(conflev))])];
        end
        
        d=find(strcmp('conflevlab',fplots));
        if d>0
            conflevlab=plots.conflevlab;
        else
            conflevlab=1;
        end
        
        d=find(strcmp('Tag',fplots));
        if d>0
            tag=plots.Tag;
        else
            tag='pl_mmdinv';
        end
        
        % Font size for the number of the axes and the titles of the
        % axes
        d=find(strcmp('FontSize',fplots));
        if d>0
            FontSize=plots.FontSize;
        else
            FontSize=12;
        end
        
        % Font size for text messages associated with the
        % labels of the horizontal lines
        d=find(strcmp('FontSizeLab',fplots));
        if d>0
            FontSizeLab=plots.FontSizeLab;
        else
            FontSizeLab=12;
        end
        
    else
        
        xlimx=[mmdinv(1,1) mmdinv(end,1)+1];
        ylimy=[min([mmdncoord;norminv(0.009)]) max([mmdncoord;norminv(0.991)])];
        LineWidth=2;
        LineWidthEnv=1;
        tag='pl_mmdinv';
        FontSize=12;
        FontSizeLab=12;
        LineStyle={'-'};
        conflev=[0.01 0.5 0.99];
        conflevlab=1;
    end
    
    % Specify where to send the output of the monitoring of mmd
    % in normal coordinates
    hmle=findobj('-depth',1,'tag',tag);
    if (~isempty(hmle))
        clf(hmle);
        figure(hmle)
        axes;
    else
        figure;
        set(gcf,'Name','mmd in normal coordinates');
    end
    
    plot1=plot(mmdinv(:,1),mmdncoord,'LineWidth',LineWidth);
    
    % Specify the line style of the trajectory of mmd in normal coordinates
    set(plot1,{'LineStyle'},LineStyle);
    
    
    set(gcf,'Tag',tag)
    
    if ~isempty(xlimx)
        xlim(xlimx);
    end
    
    if ~isempty(ylimy)
        ylim(ylimy);
    end
    
    % Add horizontal lines associated with confidenve levels
    v=axis;
    ninv=norminv(conflev);
    for i=1:length(conflev)
        if conflev(i)==0.5
            col=[1 0.69 0.39];
        elseif conflev(i)>=0.01 &&  conflev(i)<=0.99
            col='b';
        else
            col='r';
        end
        line(v(1:2)',[ninv(i);ninv(i)],'color',col,'LineWidth',LineWidthEnv,'LineStyle','--','Tag','env');
    end
    
    % add text label associated to horizontal confidence levels
    if conflevlab==1
        text(v(1)*ones(length(ninv),1),ninv+0.2,strcat(num2str(100*conflev'),'%'),...
            'FontSize',FontSizeLab,'HorizontalAlignment','Left');
    end
    
    xlabel('Subset size m','FontSize',FontSize);
    ylabel('mmd in normal coordinates','FontSize',FontSize);
    set(gca,'FontSize',FontSize)
    
end

end

%FScategory:MULT-Multivariate

