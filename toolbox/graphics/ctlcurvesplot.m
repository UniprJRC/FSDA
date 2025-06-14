function outCTLnew=ctlcurvesplot(outCTL,varargin)
%ctlcurvesplot plots the output of routine ctlcurves
%
%<a href="matlab: docsearchFS('ctlcurvesplot')">Link to the help function</a>
%
%   ctlcurvesplot takes as input the output of function ctlcurves (that is
%   a series of matrices which contain the values of the CTL bands, the set
%   of tentative solutions found using the intersection of the bands, and
%   the likelihood ratio tests for testing k vs k+1 given alpha.
%   The plot enables interaction in the sense that, if option
%   databrush has been activated, it is possible to click on a point in the
%   plot of the ctlcurves and to see the associated classification in the
%   scatter plot matrix.
%
%  Required input arguments:
%
%           outCTL : output structure produced by function ctlcurves. Structure.
%                It contains the following fields.
%             outCTL.Mu = cell of size length(kk)-by-length(alpha)
%                       containing the estimate of the centroids for each
%                       value of k and each value of alpha. More precisely,
%                       suppose kk=1:4 and alpha=[0 0.05 0.1], out.Mu{2,3}
%                       is a matrix with two rows and v columns containing
%                       the estimates of the centroids obtained when
%                       alpha=0.1.
%         outCTL.Sigma = cell of size length(kk)-by-length(alpha)
%                       containing the estimate of the covariance matrices
%                       for each value of k and each value of alpha. More
%                       precisely, suppose kk=1:4 and alpha=[0 0.05 0.1],
%                       out.Sigma{2,3} is a 3D  array of size v-by-v-by-2
%                       containing the estimates of the covariance matrices
%                       obtained when alpha=0.1.
%          outCTL.Pi   = cell of size length(kk)-by-length(alpha)
%                       containing the estimate of the group proportions
%                       for each value of k and each value of alpha. More
%                       precisely, suppose kk=1:4 and alpha=[0 0.05 0.1],
%                       out.Pi{2,3} is a 3D  array of size v-by-v-by-2
%                       containing the estimates of the covariance matrices
%                       obtained when alpha=0.1.
%          outCTL.IDX   = cell of size length(kk)-by-length(alpha)
%                       containing the final assignment for each value of k
%                       and each value of alpha. More precisely, suppose
%                       kk=1:4 and alpha=[0 0.05 0.1], out.IDX{2,3} is a
%                       vector of length(n) containing the containinig the
%                       assignment of each unit obtained when alpha=0.1.
%                       Elements equal to zero denote unassigned units.
%         outCTL.CTL    = matrix of size length(kk)-by-length(alpha)
%                       containing the values of the trimmed likelihood
%                       curves for each value of k and each value of alpha.
%                       This output (and all the other output below which
%                       start with CTL) are present only if input option
%                       bands is true or is a struct. All the other fields
%                       belo
%    outCTL.CTLbands    = 3D array of size
%                       length(kk)-by-length(alpha)-by-nsimul containing
%                       the nsimul replicates of out.CTL.
%      outCTL.CTLlikLB    =  matrix of size length(kk)-by-length(alpha)
%                       containing the lower confidence bands of the
%                       trimmed likelihood curves for each value of k and
%                       each value of alpha.
%      outCTL.CTLlikUB    =  matrix of size length(kk)-by-length(alpha)
%                       containing the upper confidence bands of the
%                       trimmed likelihood curves for each value of k and
%                       each value of alpha. T
%      outCTL.CTLlik050    =  matrix of size length(kk)-by-length(alpha)
%                       containing the central confidence bands of the
%                       trimmed likelihood curves for each value of k and
%                       each value of alpha.
%    outCTL.CTLtentSol  = matrix with size m-by-4. Details of the ordered
%                          solutions where there was intersection between
%                          two consecutive trimmed likelihood curves. First
%                          column contains the value of k, second column
%                          the value of alpha, third column the index
%                          associated to the best value of alpha, fourth
%                          colum index associated with the best value of
%                          kk.
%    outCTL.CTLoptimalAlpha = scalar, optimal value of trimming.
%     outCTL.CTLoptimalK = scalar, optimal number of clusters, stored
%                        as a positive integer value.
%   outCTL.CTLoptimalIDX  = n-by-1 vector containing assignment of each unit to
%                       each of the k groups in correspodence of
%                       OptimalAlpha and OptimalK. Cluster names are
%                       integer numbers from 1 to k. 0 indicates trimmed
%                       observations. The fields which follow which start
%                       with LRT refer to the likilhood ratio test
%
%     outCTL.LRTpval =  table with size length(kk)-1-times-length(alpha)
%                           which stores the relative frequency in which
%                           the Likelihood ratio test is greater than the
%                           corresponding bootstrap test.
%                        as a positive integer value.
%     outCTL.LRTtentSol  = matrix with size m-by-8. Details of the ordered
%                         solutions using the likelihood ratio test.
%                         First column (index): the index number of the
%                         solution.
%                         Second column (k): the value of k.
%                         Third column (alpha): the value of alpha.
%                         Fourth column (Truesol) contains 1 if the
%                         p-value beyond the threshold is always above $k$
%                         shown in the second column for all $k^* >k$ else
%                         it contains 0.
%                         Fifth and sixth columns (kindex and
%                         alphaindex) contain the index numbers of optimal
%                         input values kk and alpha.
%                         Seventh column (kbestGivenalpha) contains 1 in
%                         correspondence of the best solution for each
%                         value of k (given alpha).
%                         Eight column (ProperSize) contains 1 if
%                         the solution which has been found has a minimum
%                         group size which is greater than n*max(alpha).
%   outCTL.LRTtentSolt  = table with the same size of LRTtentSol containing
%                        the same information of array  out.TentSolLR in
%                        table format.
%    outCTL.LRTtentSolIDX = matrix with size n-by-size(out.LRTtentSol,1) with
%                        the allocation associated with the tentative
%                        solutions found in out.LRTtentSol. First column
%                        refers to solution in row 1 of out.LRTtentSol ...
%    outCTL.LRToptimalAlpha = scalar, optimal value of trimming.
%    outCTL.LRToptimalK = scalar, optimal number of clusters, stored
%                        as a positive integer value.
%    outCTL.LRToptimalIDX  = n-by-r vector containing assignment of each unit to
%                       each of the k groups in correspondence of
%                       OptimalAlpha and OptimalK. Cluster names are
%                       integer numbers from 1 to k. 0 indicates trimmed
%                       observations. The fields which follow which start
%                       with LRT refer to the likelihood ratio test.
%              outCTL.kk = vector containing the values of k (number of
%                       components) which have been considered. This  vector
%                       is equal to input optional argument kk if kk had been
%                       specified else it is equal to 1:5.
%             outCTL.alpha = vector containing the values of the trimming
%                       level which have been considered. This
%                       vector is equal to input optional argument alpha.
%            outCTL.tclustPAR = a structure  containing all the parameters
%                       used by tclust (restrfactor, mixt, restrtype.
%                       refsteps, equalweights, reftol, cshape, nsamp,
%                       nsampExtra, outliersFromUniform, nsimulExtra and
%                       usepriorSolExtra).
%             outCTL.Y  = Original data matrix Y. The field is present if
%                       option 
%                 Data Types - struct
%
%  Optional input arguments:
%
%    crit : criterion for sgnificance. Scalar in the in interval (0 1) or
%               empty value.
%                Scalar which defines the p-value threshold to define a
%                solution as significant and to reject the null hypothesis
%                of smaller groups in the likelihood ratio test. For
%                example when crit is 0.05 when we test $k$ against $k+1$
%                given a certain value of $\alpha$, if the empirical
%                p-value is greater than crit we accept the null hypothesis
%                of $k$ groups. The default value of crit is empty, that is
%                we do not change the tentative solutions found by
%                producere ctlcurves.
%                Example - 'crit',0.05
%                Data Types - double  
%
%     thresh    :   threshold which defines where to put NaN in the
%                   out.pvalLRtest matrix. Scalar in the interval [0 1].
%                   The default value is 0.05. In other words the
%                   subsequent values in a particular column of
%                   out.pvalLRtest which follow a number breater than 0.05
%                   will be set to NaN.
%                   Example - 'thresh',0.10
%                   Data Types - char
%
%    tagCtl     :   Personalized tag for CTL curves plot. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag 'pl_Ctl' for the classification
%                   likelihood curves plot with bands.
%                   Note that if the program finds a plot which has a tag
%                   equal to the one specified by the user, then the output
%                   of the new plot overwrites the existing one in the same
%                   window else a new window is created.
%                   Example - 'tag','myplot'
%                   Data Types - char
%
% tagPortofino  :   Personalized tag for Portofino plot. String. String which identifies the
%                   handle of the plot which is about to be created. The
%                   default is to use tag 'pl_Portofino' for the Portofino plot.
%                   Note that if the program finds a plot which has a tag
%                   equal to the one specified by the user, then the output
%                   of the new plot overwrites the existing one in the same
%                   window else a new window is created.
%                   Example - 'tag','myplot1'
%                   Data Types - char
%
%         conflev : confidence level of the bands. Empty value (default) or scalar.
%                   Scalar in the interval (0 1) which contains
%                   the confidence level of the bands.
%                   The default is to use the confidence level taken from
%                   ctlcurves.m (that is 50 per cent confidence
%                   level).
%                   Example - 'conflev',0.9
%                   Data Types - double
%
%   datatooltip :   interactive clicking. Empty value (default) or
%                   structure. The default is datatooltip=''.
%                   If datatooltip = 1, the user can select with the
%                   mouse a solution in order to
%                   have the following information:
%                   1) value of k which has been selected
%                   2) value of alpha which has been selected
%                   3) frequency distribution of the associated
%                   classification
%                   If datatooltip is a structure it may contain the
%                   following the fields
%                   datatooltip.DisplayStyle = Determines how the data
%                       cursor displays. datatip | window.
%                       - datatip displays data cursor
%                       information in a small yellow text box attached to a
%                       black square marker at a data point you interactively
%                       select.
%                       - window displays data cursor information for the
%                       data point you interactively select in a floating
%                       window within the figure.
%                   datatooltip.SnapToDataVertex=  Specifies whether the
%                       data cursor snaps to the nearest data value or is
%                       located at the actual pointer position.  on | off.
%                       - on data cursor snaps to the nearest data value
%                       - off data cursor is located at the actual pointer
%                       position.
%                   (see the MATLAB function datacursormode or the examples
%                   below). Default values are datatooltip.DisplayStyle =
%                   'Window' and datatooltip.SnapToDataVertex = 'on'.
%                   Example - 'datatooltip',''
%                   Data Types - scalar double or struct
%
%    databrush  :   interactive mouse brushing. Empty value, scalar or structure.
%                   If databrush is an empty value (default), no brushing
%                   is done. The activation of this option (databrush is a
%                   scalar or a structure) enables the user to select a set
%                   of values of ctl curves in the current plot and to see
%                   the corresponding classification highlighted in the
%                   scatter plot matrix (spm). If spm does not exist it is
%                   automatically created. Please, note that the window
%                   style of the other figures is set equal to that which
%                   contains the ctl plot. In other words, if the ctl plot
%                   is docked all the other figures will be docked too.
%                   DATABRUSH IS A SCALAR. If databrush is a scalar the
%                   default selection tool is a rectangular brush and it is
%                   possible to brush only once (that is persist='').
%                   DATABRUSH IS A STRUCTURE. If databrush is a structure,
%                   it is possible to use all optional arguments of
%                   function selectdataFS.m and the following optional
%                   arguments: - databrush.persist = repeated brushing
%                   enabled. Persist is an empty value or a scalar
%                     containing the strings 'on' or 'off'.
%                     The default value of persist is '', that is brushing
%                     is allowed only once.
%                     If persist is 'on' or 'off' brushing can be done as
%                     many time as the user requires.
%                     If persist='on' then the corresponding spmplot of the
%                     solutions currently brushed are added to those
%                     previously brushed.
%                     If persist='off' every time a new brush is performed
%                     scatter plot matrices corresponding to previously
%                     brushed solutions are removed.
%                   - dispopt = string which controls how to fill the
%                      diagonals in the scatterplot matrix of the brushed
%                      solutions. Set dispopt to 'hist' (default) to plot
%                      histograms, or 'box' to plot boxplots.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
%         nameY  : variable labels. Cell array. Cell array of strings
%                   containing the labels of the
%                   variables. As default value, the labels which are added
%                   are Y1, ..., Yv.
%                   Example - 'nameY',{'myY1', 'myY2'}
%                   Data Types - cell
%
%  Output:
%
%
%  outCTLnew    : structure. This structure is exactly equal to input
%                 structure outCTL unless optional input option crit is
%                 specfied. In this case using the new value of crit
%                 starting from the input table of p-values outCTL.LRTpval
%                 we can obtain a new set of tentative solutions
%                 (outCTLnew.LRTtentSol outCTLnew.LRTtentSolt) and
%                 associated classification (outCTLnew.LRTtentSolIDX), a
%                 new optimal value of k (outCTLnew.LRToptimalK) a new
%                 optimal alpha (outCTLnew.LRToptimalAlpha) a new optional
%                 classification (outCTLnew.LRToptimalIDX). In what follows
%                 we shows show just the fields which might change
%   outCTLnew.LRTtentSol  = matrix with size m-by-8. Details of the ordered
%                         solutions using the likelihood ratio test.
%                         First column (index): the index number of the
%                         solution.
%                         Second column (k): the value of k.
%                         Third column (alpha): the value of alpha.
%                         Fourth column (Truesol) contains 1 if the
%                         p-value beyond the threshold is always above $k$
%                         shown in the second column for all $k^* >k$ else
%                         it contains 0.
%                         Fifth and sixth columns (kindex and
%                         alphaindex) contain the index numbers of optimal
%                         input values kk and alpha.
%                         Seventh column (kbestGivenalpha) contains 1 in
%                         correspondence of the best solution for each
%                         value of k (given alpha).
%                         Eight column (ProperSize) contains 1 if
%                         the solution which has been found has a minimum
%                         group size which is greater than n*max(alpha).
% outCTLnew.LRTtentSolt  = table with the same size of LRTtentSol containing
%                        the same information of array  out.TentSolLR in
%                        table format.
%  outCTLnew.LRTtentSolIDX = matrix with size n-by-size(out.LRTtentSol,1) with
%                        the allocation associated with the tentative
%                        solutions found in out.LRTtentSol. First column
%                        refers to solution in row 1 of outCTLnew.LRTtentSol ...
%  outCTLnew.LRToptimalAlpha = scalar, optimal value of trimming.
%  outCTLnew.LRToptimalK = scalar, optimal number of clusters, stored
%                        as a positive integer value.
%  outCTLnew.LRToptimalIDX  = n-by-r vector containing assignment of each unit to
%                       each of the k groups in correspondence of
%                       OptimalAlpha and OptimalK. Cluster names are
%                       integer numbers from 1 to k. 0 indicates trimmed
%                       observations. The fields which follow which start
%                       with LRT refer to the likelihood ratio test.

%
% See also: ctlcurves, tclust
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('ctlcurvesplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% ctlcurves and Portofino plot for the Gyeser data with all default options.
    Y=load('geyser2.txt');
    outCTl=ctlcurves(Y,'plots',false,'nsamp',20)
    ctlcurvesplot(outCTl);
%}

%{
    % Interactive_example
    %   Example of the use of option databrush.
    %   (brushing is done only once using a rectangular selection tool)
    %  Use the geyser data.
    Y=load('geyser2.txt');
    databrush=struct;
    databrush.selectionmode='Rect';
    outCTl=ctlcurves(Y,'plots',false,'nsamp',20)
    ctlcurvesplot(outCTl,'databrush',databrush);
%}

%{
    % Interactive_example
    %   Example of the use of option databrush.
    %   (brushing is persistent)
    %  Use the geyser data.
    Y=load('geyser2.txt');
    databrush=struct;
    databrush.persist='on';
    outCTl=ctlcurves(Y,'plots',false,'nsamp',20)
    ctlcurvesplot(outCTl,'databrush',databrush);
%}


%% Beginning of code

% no brushsing as default
databrush='';
datatooltip=1;
tagCtl='tagCtl';
tagPortofino='tagPortofino';
thresh=0.05;
conflev=[];
crit=[];

if nargin>1
    options=struct('datatooltip',datatooltip, 'tagCtl',tagCtl,...
        'tagPortofino',tagPortofino,'conflev',conflev,...
        'databrush', databrush,'nameY','','thresh',thresh,'crit',crit);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)

        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:ctlcurvesplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end

        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:ctlcurvesplot:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end


    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    datatooltip=options.datatooltip;
    databrush=options.databrush;
    tagCtl=options.tagCtl;
    tagPortofino=options.tagPortofino;
    thresh=options.thresh;
    conflev=options.conflev;
    crit=options.crit;
end

% Extract the values of k (number of groups)
kk=outCTL.kk;
alpha=outCTL.alpha(:);
lkk=length(kk);
lalpha=length(alpha);

if isfield(outCTL,'likLB')
    ComputeBands=true;
    if isempty(conflev)
        likLB=outCTL.likLB;
        likUB=outCTL.likUB;
        lik050=outCTL.lik050;
    else
        BandsCTL=outCTL.BandsCTL;
        likLB = zeros(lkk,lalpha);
        likUB =likLB;
        gamma=(1-conflev)/2;
        for k=1:length(kk)  % loop for different values of k (number of groups)
            parfor j=1:lalpha
                likLB(k,j) = quantile(BandsCTL(k,j,:), gamma);
                likUB(k,j) = quantile(BandsCTL(k,j,:), 1- gamma);
            end
        end
        lik050=outCTL.lik050;
        IDX=outCTL.IDX;
        % Recall routine which computes the best tentative solutions using
        % confev just supplied
        [TentSol]=findOptimalSolutions(likUB,likLB,lik050,IDX,alpha,lkk);
        disp(['New set of tentative solutions using conflev=' num2str(conflev)])
        disp(TentSol)
    end

else
    ComputeBands=false;
end

CTLVal=outCTL.CTL;

% Global variables to declare is datatooltip is not empty
if ~isempty(datatooltip)
    hTarget=[];
    hTargetlwd=[];
    hTargetcol=[];
end

%% Prepate the figure to display the ctlplot
% Create a figure to host the plot or clear the existing one
if isempty(tagCtl)
    hCTL=figure;
else
    hCTL=findobj('-depth',1,'tag',tagCtl);
    if (~isempty(hCTL))
        for i=1:length(hCTL)
            clf(hCTL(i));
        end
        figure(hCTL(1));
        axes;
    else
        hCTL=figure;
    end
end
set(hCTL,'Name', 'ctlcurves plot', 'NumberTitle', 'off');
hold('all');

if ComputeBands==true
    linetype1 = repmat({'-.','-.','-.','-.','-.'},1,10);
    color = repmat({'r','g','b','c','k'},1,10);
    LineWidth = 1;
    hold('on')
    for i = 1:length(kk)
        plot(alpha,likLB(i,:), 'LineStyle',linetype1{i}, 'Color', color{i}, 'LineWidth', LineWidth)
        % plot(alphaTrim,lik050(i,:), 'LineStyle',linetype{i}, 'Color', color{i}, 'LineWidth', LineWidth)
        text(alpha(end),lik050(i,end),[' k = ' num2str(kk(i))],'FontSize',16, 'Color', color{i})
        plot(alpha,likUB(i,:), 'LineStyle',linetype1{i}, 'Color', color{i}, 'LineWidth', LineWidth)
    end
else
    plot(alpha', CTLVal')
    one=ones(length(alpha),1);
    for i = 1:length(kk)
        text(alpha',CTLVal(i,:)',num2str(kk(i)*one),'FontSize', 14)
    end
end
xlabel('Trimming level alpha')
ylabel('Log likelihood')
set(gca,'XTick',alpha);
title('Classification Trimmed Likelihood Curves')

if isempty(tagCtl)
    set(gcf,'Tag','pl_Ctl');
else
    set(gcf,'Tag',tagCtl);
end

plot1=gcf;
% Datatooltip mode (call to function ICplotLbl)
if ~isempty(datatooltip)
    PrepareDatatooltip(outCTL)
end



tbootGTtobs=outCTL.LRTpval{:,:};


%% Create Portofino plot
if isfield(outCTL,'pvalLRtest')

    if isempty(tagPortofino)
        hPorto=figure;
    else
        hPorto=findobj('-depth',1,'tag',tagPortofino);
        if (~isempty(hPorto))
            clf(hPorto);
            figure(hPorto);
            axes;
        else
            hPorto=figure;
        end
    end
    set(hPorto,'Name', 'Portofino plot', 'NumberTitle', 'off');
    hold('all');

    tbootGTtobsf=tbootGTtobs;

    for j=1:lalpha
        nextj=false;
        for i=1:lkk-1
            if tbootGTtobs(i,j)>thresh
                tbootGTtobsf(i+1:end,j)=NaN;
                nextj=true;
            end
            if nextj==true
                break
            end
        end
    end

    tbootGTtobsf(tbootGTtobsf==0)=NaN;
    plot(alpha,tbootGTtobsf','x','LineWidth',3)
    % plot(alpha,a','LineWidth',3)
    area(alpha,tbootGTtobsf','LineWidth',3,'FaceColor',[0.9290 0.6940 0.1250],'FaceAlpha',0.3) % )
    xlabel('Level of trimming')
    ylabel('p-value of k groups')
    ylim([0 1.05])
    for i=1:size(tbootGTtobsf,1)
        FindFirstNotNaN=find(~isnan(tbootGTtobsf(i,:)),1);
        text(alpha(FindFirstNotNaN)+0.005,tbootGTtobsf(i,FindFirstNotNaN),['k=' num2str(i)])
    end
    set(gca,'XTick',alpha)
    set(gca,'XDir','reverse')

    if isempty(tagPortofino)
        set(gcf,'Tag','pl_Portofino');
    else
        set(gcf,'Tag',tagPortofino);
    end
    title('Portofino plot')
    if ~isempty(datatooltip)
        outCTL.tbootGTtobsf=[tbootGTtobsf;NaN(1,lalpha)];
        PrepareDatatooltip(outCTL)
    end

end

outCTLnew=outCTL;

if ~isempty(crit)

    PiVal=outCTL.Pi;
    MuVal=outCTL.Mu;
    SigmaVal=outCTL.Sigma;
    IDX=outCTL.IDX;
    alphaTrim=outCTL.alpha;

    Y=outCTL.Y;
    n=size(Y,1);
    LRTtentSol=zeros(lkk-1,8);
    LRTtentSolIDX=zeros(n,lkk-1);

    % Recall from outCTL.tclustPAR all the parameters used by tclust
    tclustPAR=outCTL.tclustPAR;
    outliersFromUniform=tclustPAR.outliersFromUniform;
    nsimulExtra=tclustPAR.nsimulExtra;
    usepriorSolExtra=tclustPAR.usepriorSolExtra;
    restrfactor=tclustPAR.restrfactor;
    mixt=tclustPAR.mixt;
    restr=tclustPAR.restrtype;
    refsteps=tclustPAR.refsteps;
    equalweights=tclustPAR.equalweights;
    reftol=tclustPAR.reftol;
    cshape=tclustPAR.cshape;
    nsampExtra=tclustPAR.nsampExtra;


    ij=1;

    for j=1:lalpha
        % For each value of alpha soli and indexSpuriousSolution are reset
        indexSpuriousSolution=true;
        soli=0;
        increment=1;
        while indexSpuriousSolution==true

            soli=find(tbootGTtobs((soli+increment):end,j)>=crit,1)+soli+increment-1;
            if ~isempty(soli)

                % Store solution number, value of k, value of alpha
                idxij=IDX{soli,j};
                LRTtentSolIDX(:,ij)=idxij;
                tabij=tabulate(idxij(idxij>0));
                if min(tabij(:,2))<n*max(alphaTrim)
                    properSize=false;
                else
                    properSize=true;
                end
                LRTtentSol(ij,[1:3 5:6 8])=[ij kk(soli) alphaTrim(j) soli j properSize];

                if soli<lkk

                    % sum(tbootGTtobs(soli+1:end,j)>=crit)==lkk-soli
                    if all(tbootGTtobs(soli+1:end,j)>=crit)
                        LRTtentSol(ij,4)=1;
                        indexSpuriousSolution=false;
                    else
                        LRTtentSol(ij,4)=0;
                        increment=find(tbootGTtobs(soli+1:end,j)<crit,1);

                    end
                    % Check whether this solution had been previously found
                    % that is if for the same value of k  we had already
                    % obtained the same solution
                    if ij>1
                        if LRTtentSol(ij,2) == LRTtentSol(ij-1,2) && LRTtentSol(ij,4)==LRTtentSol(ij-1,4)
                            LRTtentSol(ij,:)=[];
                            LRTtentSolIDX(:,ij)=[];
                            ij=ij-1;
                        end
                    end
                else
                    LRTtentSol(ij,4)=1;
                    indexSpuriousSolution=false;
                end

                % idxLR(:,ij)=IDX{i,soli};
                ij=ij+1;
            else
                indexSpuriousSolution=false;
            end
        end
    end

    if size(LRTtentSol,1)>=1
        LRTtentSol=LRTtentSol(1:ij-1,:);
        LRTtentSolIDX=LRTtentSolIDX(:,1:ij-1);

        numsol=(1:size(LRTtentSol,1))';
        nsoleti="Sol"+numsol;

        % FOR EACH VALUE OF ALPHA COMPUTE ADDITIONAL LIK RATIO TESTS TO FIND
        % UNIQUE BEST SOLUTION
        % For example suppose that for a particular alpha we found that there
        % are the potential solutions k=3, k=6 and k=10 then LR test between
        % k=3 and k=6 is performed, and the best among these two values of k
        % is stored and is later tested against k=10
        % The best value among these 3 candidates gets a value of 1 in the last
        % column of matrix LRTtentSolt
        %         nsimulExtra=50;
        %         nsampExtra=5000;
        %         usepriorSolExtra=false;

        for j=1:lalpha
            alphaTrimj=alphaTrim(j);
            MultSolalpha=find(LRTtentSol(:,3)==alphaTrimj);
            kTentative=LRTtentSol(MultSolalpha,2);
            kTentativepos=LRTtentSol(MultSolalpha,5);

            if length(MultSolalpha)>1

                indbestk=1;
                for jjj=2:length(MultSolalpha)

                    % kpos = position of first value of k to use in the additional LR test
                    kpos=kTentativepos(jjj-1);
                    % kH1pos = position larger value of k to use in the additional LR test
                    kH1pos=kTentativepos(jjj);
                    % k and kH1= true values of k  to use in the additional
                    % LR test
                    k=kTentative(jjj-1);
                    kH1=kTentative(jjj);

                    ktrue = length(PiVal{kpos, j});
                    Mutrue = MuVal{kpos, j};
                    Mutrue=Mutrue(1:ktrue,:);
                    Sigmatrue = SigmaVal{kpos, j};
                    Sigmatrue = Sigmatrue(:,:,1:ktrue);
                    Pitrue=PiVal{kpos, j};
                    alphaTrimj=alphaTrim(j);
                    ngood=round(n*(1-alphaTrimj));
                    nout=n-ngood;
                    idxkj=IDX{kpos,j};

                    % if outliersFromUniform is false the outliers in the replicate
                    % samples are the units which have been trimmed
                    if outliersFromUniform == false
                        Yadd=Y(idxkj==0,:);
                    else
                        Yadd=[];
                    end
                    tboot=zeros(nsimulExtra,1);


                    if usepriorSolExtra ==true
                        idxkH1=IDX{kH1pos,j};
                    else
                        idxkj=[];
                        idxkH1=[];
                    end
                    parfor (zz = 1:nsimulExtra)

                        % Alternative instruction using traditional for
                        % instead of parfor
                        % for zz = 1:nsimulExtra
                        if outliersFromUniform == true
                            [Ysim]=simdataset(ngood, Pitrue, Mutrue, Sigmatrue,'noiseunits', nout);
                            if size(Ysim,1)<n
                                Yadd1=repmat(Ysim(end,:),n-size(Ysim,1),1);
                                Ysim=[Ysim;Yadd1];
                            end
                        else
                            [Ysim]=simdataset(ngood, Pitrue, Mutrue, Sigmatrue,'noiseunits', 0);
                            Ysim=[Ysim;Yadd];
                        end

                        outtcSIM=tclust(Ysim,k,alphaTrimj,restrfactor,'plots',0,'msg',0,'mixt',mixt, ...
                            'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                            'reftol',reftol,'cshape',cshape,...
                            'priorSol',idxkj,'nsamp',nsampExtra);

                        outtcSIMkH1=tclust(Ysim,kH1,alphaTrimj,restrfactor,'plots',0,'msg',0,'mixt',mixt, ...
                            'restrtype',restr,'nocheck',1,'refsteps',refsteps,'equalweights',equalweights,...
                            'reftol',reftol,'cshape',cshape, ...
                            'priorSol',idxkH1,'nsamp',nsampExtra);

                        tboot(zz)=outtcSIMkH1.obj-outtcSIM.obj;

                    end


                    tobs= outCTL.CTL(kH1pos,j)-outCTL.CTL(kpos,j);

                    % tboot = values of the difference between target function using kH1
                    % and k groups for data simulated assuming k groups
                    addLRtest=sum(tboot>tobs)/nsimulExtra;

                    if addLRtest>crit
                        kTentative(jjj)=kTentative(jjj-1);
                    else
                        indbestk=indbestk+1;
                    end
                end
                LRTtentSol(MultSolalpha(indbestk),7)=1;
            elseif length(MultSolalpha) ==1
                LRTtentSol(MultSolalpha,7)=1;
            else
            end
        end

        LRTtentSolt=array2table(LRTtentSol,'RowNames',nsoleti, ...
            'VariableNames',{'index' 'k' 'alpha' 'Truesol' 'kindex' 'alphaindex' 'kbestGivenalpha' 'ProperSize'});

    else
        LRTtentSol=[];
        LRTtentSolt=[];
        LRTtentSolIDX=[];
    end

    outCTLnew.LRTtentSol=LRTtentSol;
    outCTLnew.LRTtentSolt=LRTtentSolt;
    outCTLnew.LRTtentSolIDX=LRTtentSolIDX;

    % Find  LRToptimalK, LRToptimalAlpha and LRToptimalIDX
    if ~isempty(LRTtentSol)
        % Just in case there are multiple solutions take as optimal the one
        % with the smallest alpha independently of group size (if
        % there is at least one solution which satisfies the minimum group
        % size)
        if size(LRTtentSol,1)>1
            if any(LRTtentSol(:,8))
                % indexBestSolLRT=find(LRTtentSol(:,7)==1 & LRTtentSol(:,8)==1);
                indexBestSolLRT=find(LRTtentSol(:,7)==1);
            else
                indexBestSolLRT=1;
            end
        else
            indexBestSolLRT=1;
        end
        LRToptimalK=LRTtentSol(indexBestSolLRT,2);
        LRToptimalAlpha=LRTtentSol(indexBestSolLRT,3);
        LRToptimalIDX=LRTtentSolIDX(:,indexBestSolLRT);
    else
        LRToptimalK=[];
        LRToptimalAlpha=[];
        LRToptimalIDX=[];
    end
    outCTLnew.LRToptimalK=LRToptimalK;
    outCTLnew.LRToptimalAlpha=LRToptimalAlpha;
    outCTLnew.LRToptimalIDX=LRToptimalIDX;
end


%% Interactivity through databrush
if isstruct(databrush)
    fdatabrush=fieldnames(databrush);

    % persist option
    dpers=find(strcmp('persist',fdatabrush));
    if dpers>0
        persist=databrush.persist;
        databrush=rmfield(databrush,'persist');
        fdatabrush=fieldnames(databrush);

        ColorOrd=[1 0 0;0 1 1; 1 0 1; 1 1 0; 0 0 0; 0 1 0; 0 0 1];
        ColorOrd=repmat(ColorOrd,4,1);
    else
        persist='';
        ColorOrd=[1 0 0];
    end

    % FlagColor option Initialize colors for the brushing option: default
    % colors are blue (unbrushed unit) and red (brushed units)

else
    persist='';
    ColorOrd=[1 0 0];
end

%% Brush mode (call to function selectdataFS)
if ~isempty(databrush) || isstruct(databrush)

    % Define marker type
    styp={'+';'o';'*';'x';'s';'d';'^';'v';'>';'<';'p';'h';'.'};
    styp=repmat(styp,ceil(lkk/length(styp)),1);

    IDX=outCTL.IDX;
    Y=outCTL.Y;


    if isstruct(databrush)
        % Choose what to put on the main diagonal of the spm
        d=find(strcmp('dispopt',fieldnames(databrush)),1);
        if  isempty(d)
            dispopt='hist';
        else
            dispopt=databrush.dispopt;
            databrush=rmfield(databrush,'dispopt');
            fdatabrush=fieldnames(databrush);
        end


        % transform the input structure databrush into a cell array
        cv=[fdatabrush struct2cell(databrush)]';
        sele=cv(:)';
        % Add the FlagSize of the brushed points if it has not been
        % previously specified by the user
        d=find(strcmp('FlagSize',fdatabrush));
        if d==0
            sele=[sele 'FlagSize' '3'];
        end
    else
        sele={'selectionmode' 'Rect' };
        dispopt='hist';
    end

    % Set the labels of the axes in the spmplot which shows the brushed solutions.
    d=find(strcmp('nameY',fieldnames(outCTL)),1);
    if  isempty(d)
        v=size(Y,2);
        p1=1:v;
        nameY=cellstr(num2str(p1','y%d'));
    else
        nameY=outCTL.nameY;
    end

    plo=struct;
    plo.nameY=nameY;

    % some local variables
    but=0; brushcum=[]; ij=1;

    sele=[sele 'FlagColor' ColorOrd(ij,:) 'FlagMarker' char(styp(ij+1))];

    % loop brushing
    while but<=1

        figure(plot1);

        % Remark: function selectdataFS cannot be used on the current
        % figure if the "selection mode" or the "zoom tool" are on. Setting
        % the plotedit mode initially to on and then to off, has the effect
        % to deselect plotedit mode.
        plotedit on
        plotedit off

        if strcmp(persist,'on')
            % add to cell sele option FlagColor (color of selection) and
            % FlagMarker (symbol to be used for selection)

            if ij>1
                chkexist=find(strcmp('FlagColor',sele)==1);
                sele{chkexist+1}=ColorOrd(ij,:);
                sele{chkexist+3}=char(styp(ij+1));
            end
        end


        % CALL to selectdataFS
        disp('Select a region to brush in the ctlcurves plot');
        [pl,xselect,yselect] = aux.selectdataFS(sele{:});

        % exit if the ctlcurves figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (pl == -999)
            return
        end

        if ~isempty(cell2mat(pl))

            if iscell(xselect)

                for z=1:length(xselect)
                    if z==1
                        xselect_all=xselect{z};
                        yselect_all=yselect{z};
                    else
                        xselect_all=cat(1,xselect_all,xselect{z} ) ;
                        yselect_all=cat(1,yselect_all,yselect{z} ) ;
                    end
                end
            else    %xselect is a vector if there is only one selection.
                xselect_all=xselect;
                yselect_all=yselect;
            end

            % Find unique values of alpha and k selected
            % Fist column = alpha
            % Second column = k
            % Third column = row index of selected solution
            % Fourth column = column index of selected solution
            alphak=NaN(50,4);

            ijk=1;
            seqalpha=1:lalpha;
            seqk=1:lkk;
            for i=1:length(xselect_all)
                % selalphai = index associated with the value of alpha
                % which has been selected
                selalphai=seqalpha(alpha==xselect_all(i));
                % alphai = value of alpha which has been selected
                alphai=alpha(selalphai);

                selk1=likUB(:,selalphai)==yselect_all(i);
                selk2=likLB(:,selalphai)==yselect_all(i);
                % selki =index associated with the value of k which has been
                % selected
                selki=seqk(or(selk1,selk2));
                % ki = value of k which has been selected
                ki=kk(selki);
                alphak(ijk,:)=[ki alphai, selki, selalphai];
                ijk=ijk+1;
            end
            % Find the unique combinations of alpha and k
            alphak=alphak(1:ijk-1,:);
            alphak=unique(alphak,'sorted','rows');

            nbrush=1;
            disp('Values  k and alpha selected');
            disp(array2table(alphak(:,1:2),'VariableNames',{'k','Trimming level alpha',}));

        else
            disp('Wrong selection: Try again');
            disp('Select a region to brush in the CTL curves plot');
            figure(plot1);
            nbrush='';
        end

        %% For each brushing operation, do the following:
        if ~isempty(nbrush)
            % brushcum = - the list of selected observations in all
            % iterations if
            %   persist=on
            % - the list of selected observations in the current iteration
            %   if persist=off
            if strcmp(persist,'on')
                brushcum=unique([brushcum; nbrush]);
            else
                brushcum=nbrush;
            end

            for r=1:size(alphak,1)
                Iwithk=alphak(r,3);
                Jwitha=alphak(r,4);
                % It is necessary to put inside the tag the word group in
                % order to let spmplot understand that we are not dearling
                % with brushed units but simply with groups.

                figure('Tag','pl_spmCTL');
                set(gcf,'WindowStyle',get(plot1,'WindowStyle'));
                spmplot(Y,IDX{Iwithk,Jwitha},plo,dispopt);
                title(['CTL' num2str(likLB(Iwithk,Jwitha),'%10.2f') '-' num2str(likUB(Iwithk,Jwitha),'%10.2f') ...
                    ' k='  num2str(kk(Iwithk)) ' alpha='  num2str(alpha(Jwitha)) ])
            end

            %% - check condition to exit from the brush mode
            % If the option persistent is not equal off or on than get out
            % of the loop
            if strcmp('on',persist) || strcmp('off',persist)

                if strcmp('on',persist)
                    ij=ij+1;
                    % but=1 makes that previous highlightments in other
                    % figures are not deleted
                    but=1;
                end

                % Before waitforbuttonpress:
                % - the CTL curves plot is highlighted again
                figure(plot1);
                % - and a function to be executed on figure close is set
                set(gcf,'CloseRequestFcn',@aux.closereqFS);

                % Lay down the plots before continuing
                position(plot1);
                disp('Highlight the CTL plot then: click on it to continue brushing or press a keyboard key to stop');
                ss=aux.waitforbuttonpressFS;
                disp('------------------------');

                % After waitforbuttonpress:
                % - the standard MATLAB function to be executed on figure
                %   close is recovered
                set(gcf,'CloseRequestFcn','closereq');

                if  strcmp('off',persist)
                    Open_spm = findobj(0, 'type', 'figure','tag','pl_spm');
                    delete(Open_spm);
                end

                Open_mal = findobj(0, 'type', 'figure','tag','tagCtl');
                if isempty(Open_mal)  % User closed the main brushing window
                    Open_spm = findobj(0, 'type', 'figure','tag','pl_spm');

                    if ~isempty(Open_spm)
                        delete(Open_spm);
                    end % spmplot  is deleted
                    delete(get(0,'CurrentFigure')); % deletes Figure if still one left open
                end

                % - and the 'but' variable is set if keyboard key was
                % pressed
                if ss==1
                    but=2;
                end
            else
                but=2;
            end % close loop associated with persist 'on' or persist 'off'
            position(plot1)
        end % for each brushing operation do ...
    end % close loop associated with but (loop brushing)
else
    % Apply cascade to existing plots in case databrush is not invoked
    position(0);

end



    function PrepareDatatooltip(IC)
        try
            % chkgpu=gpuDeviceCount;
            % datacursormode on;
            hdt = datacursormode;
            set(hdt,'Enable','on');
            % If options.datatooltip is not a struct then use our default options
            if ~isstruct(datatooltip)
                set(hdt,'DisplayStyle','window','SnapToDataVertex','on');
            else
                % options.datatooltip contains a structure where the user can set
                % the properties of the data cursor
                set(hdt,datatooltip);
            end

            LineColor=[1 0 0];
            % Declare a custom datatooltip update function to display additional
            % information about the selected unit
            set(hdt,'UpdateFcn',{@ICplotLbl,IC,LineColor});
        catch
            disp('No graphical device, interactive datatooltip not enabled')
        end
    end

    function output_txt = ICplotLbl(~,event_obj,IC,~)
        % ICplotLbl provides information about
        %
        % Required input arguments:
        %
        %       obj =   Currently not used (empty but necessary)
        % event_obj =   Handle to event object
        %               (event_obj=graphics.datatipevent)
        %               Remark: the first two arguments are implicit in the
        %               sense that these arguments are automatically passed
        %               to the function when it executes.
        %       out =   a structure containing the following fields
        %               IC = a matrix containing the IC
        %
        % Output:
        %
        %   output_txt=  Datatip text (string or string cell array) which
        %                informs about the value of k and c which have been
        %                selected and the frequency distribution of the
        %                associated classification
        %
        %
        %
        % Written by FSDA team

        % find the plot the user has clicked on
        titl=get(gca,'title');
        titl=titl.String;

        if strcmp(titl,'Portofino plot')
            ICselB=IC.tbootGTtobsf;
            ICselU=[];
        else
            if isfield(IC,'likLB') % bands have been computed
                ICselB=IC.likLB;
                ICselU=IC.likUB;
            else
                ICselB=IC.CTL;
                ICselU=[];
            end
        end

        ICIDXsel=IC.IDX;

        if ~isempty(hTarget) && isvalid(hTarget)
            % set old line width and old color for old selection
            %             if strcmp(titl,'Portofino plot') && isfield(struct(hTarget),'FaceColor')
            %                 set(hTarget,'LineWidth',hTargetlwd,'FaceColor',hTargetcol);
            %             else
            %                 set(hTarget,'LineWidth',hTargetlwd,'Color',hTargetcol);
            %             end

            try
                set(hTarget,'LineWidth',hTargetlwd,'FaceColor',hTargetcol);
            catch
                set(hTarget,'LineWidth',hTargetlwd,'Color',hTargetcol);
            end
        else
        end

        % Store line width and color of selected trajectory
        % Notice that changing event_obj.Target in subsequent lines seems
        % to affect also hTarget
        hTarget=event_obj.Target;
        hTargetlwd=get(hTarget,'LineWidth');
        if strcmp(titl,'Portofino plot')
            hTargetcol=get(hTarget,'FaceColor');
            % Increase Line width and keep line color
            set(hTarget,'LineWidth',hTargetlwd+1.5,'FaceColor',hTargetcol);
        else
            hTargetcol=get(hTarget,'Color');
            % Increase Line width and keep line color
            set(hTarget,'LineWidth',hTargetlwd+1.5,'Color',hTargetcol);
        end


        pos = get(event_obj,'Position');

        % x and y, plot coordinates of the mouse
        x1 = pos(1); y1 = pos(2);

        % Find index to retrieve value of k (number of groups) and value of c (restriction factor)
        % Consider that find return the
        % linear indexing of matrix xydata

        % Linear indexing is transformed into normal indexing using
        % function ind2sub row and column contain the column and row
        % indexed of the observation which has been selected with the mouse
        idx = find(ICselB == y1,1);
        if isempty(idx)
            idx = find(ICselU == y1,1);
            [row,col] = ind2sub(size(ICselU),idx);
        else
            [row,col] = ind2sub(size(ICselB),idx);
        end



        if isempty(row)
            output_txt{1}=['no ctl x,y' num2str(x1) '' num2str(y1)] ;
        else
            output_txt=cell(4,1);

            % output_txt is what it is shown on the screen
            output_txt(1,1) = {[ titl ' equal to: ',num2str(y1,4)]};

            % Add information about k and alpha
            output_txt{2,1} = ['k= ' num2str(row) ', alpha=' num2str(IC.alpha(col))];


            output_txt{3,1} = 'Classification';

            % Add information about the corresponding frequency
            % distribution  of associated classification
            clas=tabulate(ICIDXsel{row,col});
            output_txt{4,1} = num2str(clas);

            set(0,'ShowHiddenHandles','on');    % Show hidden handles
            hText = findobj('Type','text','Tag','DataTipMarker');
            set(hText,'Interpreter','latex');
        end
    end

end


function [TentSol,kfin,alphafin,idxOptimal]=findOptimalSolutions(likUB,likLB,lik050,IDX,alphaTrim,lkk)
conv=0;
% First column of TentSol will contain the value of k while the second
% column the associated trimming level
% Third column is the index of the elment of alpha containing the best
% optimal trimming level
TentSol=zeros(lkk-1,3);
jj=1;

for j = 1:lkk-1
    alphaBest='';
    for jalpha =1:size(likUB,2)
        if (likUB(j,jalpha) > likLB(j+1,jalpha) &&  lik050(j+1,jalpha)> lik050(j,jalpha)) || ...
                (likLB(j,jalpha) < likUB(j+1,jalpha) &&  lik050(j+1,jalpha) < lik050(j,jalpha))
            % Find best trimming level
            alphaBest = alphaTrim(jalpha);
            break
        else
        end

    end

    if ~isempty(alphaBest)
        TentSol(jj,:)=[j alphaBest, jalpha];
        jj=jj+1;
    end
end
if TentSol(1,1)>0
    conv=1;
    TentSol=TentSol(1:jj-1,:);
    kfin=TentSol(1,1);
    jalpha=TentSol(1,3);
    alphafin=TentSol(1,2);
else
    TentSol=NaN;
end


if conv == 1
    idxOptimal=IDX{kfin,jalpha};
else
    disp('No intersection among the curves has been found for the selected trimming levels and number of groups')
    disp('Please increase k or alpha')
    alphafin=max(alphaTrim);
    kfin=max(kk);
    idxOptimal = IDX{end, end};
end
end
%FScategory:VIS-Clu


