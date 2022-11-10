function h  = carbikeplotGPCM(RelSol,varargin)
%carbikeplot produces the carbike plot to find best relevant clustering solutions
%
%<a href="matlab: docsearchFS('carbikeplotGPCM')">Link to the help function</a>
%
%   carbikeplotGPCM takes as input the output of function tclustICsolGPCM
%   (that is a structure containing the best relevant solutions) and
%   produces the car-bike plot. This plot provides a concise summary of the
%   best relevant solutions. This plot shows on the horizontal axis the
%   value of $c_{det}$ and $c_{shw}$ restriction factors and on the
%   vertical axis the value of $k$. For each solution we draw two rectangle
%   (associated with  $c_{det}$ and $c_{shw}$) which are respectively
%   referred to interval of values for which the solution is best and
%   stable and a horizontal line which departs from the rectangle for the
%   values of $c_{det}$ ($c_{shw}$) in which the solution is only stable.
%   Finally, for the best value of $c_{det}$ ($c_{shw}$) associated to the
%   solution, we show a circle with a  number indicating the ranked
%   solution among those which are not spurious. This plot has been
%   baptized ``car-bike'', because the first best solutions (in general 2
%   or 3) are generally best and stable for a large number of values of $c$
%   and therefore will have large rectangles. In addition, these solutions
%   are likely to be stable for additional values of $c_det$ ($c_{shw}$)
%   and therefore are likely to have horizontal lines departing from the
%   rectangles (from here the name ``cars''). Finally, local minor
%   solutions (which are associated with particular values of $c_{det}$
%   ($c_{shw}$) and $k$) do not generally present rectangles or lines and
%   are shown with circles (from here the name ``bikes'')
%
%  Required input arguments:
%
%           RelSol : Relevant solutions produced by function tclustICsolGPCM. Structure.
%                It contains the following fields:
%
%   RelSol.MIXMIXbs = cell of size NumberOfBestSolutions-times-8 which contains
%                the details of the best solutions for MIXMIX (BIC).
%                Each row refers to a solution.  The information which is
%                stored in the columns is as follows.
%                1st col = scalar, value of k for which solution takes place;
%                2nd col = scalar, value of cdet for which solution takes place;
%                3rd col = row vector of length d which contains the values
%                   of cdet for which the solution is uniformly better.
%                4th col = row vector of length d+r which contains the
%                   values of cdet for which the solution is considered
%                   stable (i.e. for which the value of the adjusted Rand
%                   index, or the adjusted Fowlkes and Mallows index) does
%                   not go below the threshold defined in input option
%                   ThreshRandIndex).
%                5th col = string which contains 'true' or 'spurious'. The
%                   solution is labelled spurious if the value of the
%                   adjusted Rand index with the previous solutions is
%                   greater than ThreshRandIndex.
%                6th col = scalar, value of cshw for which solution takes
%                   place.
%                7th col = row vector of length d which contains the values
%                   of cshw  for which the solution is uniformly better.
%                8th col = row vector of length d+r which contains the
%                   values of cshw for which the solution is considered
%                   stable (i.e. for which the value of the adjusted Rand
%                   index, or the adjusted Fowlkes and Mallows index) does
%                   not go below the threshold defined in input option
%                   ThreshRandIndex).
%               Remark: field out.MIXMIXbs is present only if input option
%               'whichIC' is 'ALL' or 'whichIC' is 'MIXMIX'.
%
%  RelSol.MIXMIXbsari =  matrix of adjusted Rand indexes (or Fowlkes and Mallows
%               indexes) associated with the best
%                solutions for MIXMIX. Matrix of size
%                NumberOfBestSolutions-times-NumberOfBestSolutions whose
%                i,j-th entry contains the adjusted Rand index between
%                classification produced by solution i and solution j,
%                $i,j=1, 2, \ldots, NumberOfBestSolutions$.
%               Remark: field out.MIXMIXbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'MIXMIX'.
%
%  RelSol.MIXCLAbs = this output has the same structure as out.MIXMIXbs but
%               it is referred to MIXCLA.
%               Remark: field out.MIXCLAbs is present only if 'whichIC' is
%               'ALL' or 'whichIC' is 'MIXCLA'.
%
% RelSol.MIXCLAbsari = this output has the same structure as out.MIXMIXbs but
%               it is referred to MIXCLA.
%               Remark: field out.MIXCLAbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'MIXCLA'.
%
%  RelSol.CLACLAbs = this output has the same structure as out.MIXMIXbs but
%               it is referred to CLACLA.
%               Remark: field out.CLACLAbs is present only if 'whichIC' is
%               'ALL' or 'whichIC' is 'CLACLA'.
%
% RelSol.CLACLAbsari = this output has the same structure as out.MIXMIXbs but
%               it is referred to CLACLA.
%               Remark: field out.MIXCLAbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'CLACLA'
%
%          RelSol.kk = vector containing the values of k (number of
%                   components) which have been considered.
%
%      RelSol.ccdet = vector containing the values of cdet (values of the
%                   restriction factor for determinants) which have been
%                   considered. 
%
%      RelSol.ccshw = vector containing the values of cshw (values of the
%                   restriction factor for shape elements inside each
%                   group) which have been considered. 
%
%          out.alpha = scalar containing the value of $\alpha$ (trimming
%                   level) which has been considered. 
%          Data Types - struct
%
%  Optional input arguments:
%
%
% SpuriousSolutions  :  Include or nor spurious solutions. Boolean. As
%                       default spurios solutions are not included into the
%                       plot.
%                 Example - 'SpuriousSolutions',false
%                 Data Types - single | double
%
%
%  Output:
%
%         h:   graphics handle to the plot. Graphics handle. Graphics
%               handle which is produced on the screen.
%
%
%
% See also: carbikeplot, tclustICgpcm, tclust, tclustICsolGPCM
%
% References:
%
% Cerioli, A. Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani, M. (2017),
% Finding the Number of Groups in Model-Based Clustering via Constrained
% Likelihoods, "Journal of Computational and Graphical Statistics", pp. 404-416,
% https://doi.org/10.1080/10618600.2017.1390469
%
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('carbikeplotGPCM')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Simulated data with 3 components.
    % Data generation
    rng('default') % Reinitialize the random number generator to its startup configuration
    rng(10);
    ktrue=3;
    % n = number of observations
    n=150;
    % v= number of dimensions
    v=2;
    % Imposed average overlap and restriction factor.
    BarOmega=0.04;
    restrfact=5;
    outMS=MixSim(ktrue,v,'BarOmega',BarOmega, 'restrfactor',restrfact);
    % data generation given centroids and cov matrices
    [Y,id]=simdataset(n, outMS.Pi, outMS.Mu, outMS.S);
    % Number of subsets to extract
    nsamp=100;
    % Computation of information criterion using MIXMIX
    outICmixt=tclustICgpcm(Y,'plots',0,'nsamp',nsamp,'kk',1:4);
    % Specify number of solutions
    NumberOfBestSolutions=3;
    % Extract the best solutions using as Information criterion MIXMIX
    [outMIXMIX]=tclustICsolGPCM(outICmixt,'whichIC','MIXMIX','plots',0,'NumberOfBestSolutions',NumberOfBestSolutions);
    carbikeplotGPCM(outMIXMIX);
%}

%{
    %% car-bike plot for the geyser data.
    Y=load('geyser2.txt');
    out=tclustICgpcm(Y,'cleanpool',false,'plots',0,'alpha',0.1,'nsamp',100,'kk',2:4);

    % Find the best solutions using as Information criterion MIXMIX
    disp('Best solutions using MIXMIX')
    [outMIXMIX]=tclustICsolGPCM(out,'whichIC','MIXMIX','plots',0,'NumberOfBestSolutions',3);
    % Produce the car-bike plot
    carbikeplotGPCM(outMIXMIX)
%}



%% Beginning of code

if ~isstruct(RelSol)
    error('FSDA:carbikeplot:WrongInput','First input argument must be a structure.');
end

SpuriousSolutions=false;

if nargin>1
    options=struct('SpuriousSolutions',SpuriousSolutions);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:carbikeplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:tclustBICsol:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end
    
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    SpuriousSolutions=options.SpuriousSolutions;
    
end


% Extract the values of k (number of groups)
kk=RelSol.kk;
% Extract the values of c_det and c_shw (values of restriction factors)
ccdet=RelSol.ccdet;
ccshw=RelSol.ccshw;

if isfield(RelSol,'CLACLAbs')
    ICbs=RelSol.CLACLAbs;
elseif isfield(RelSol,'MIXMIXbs')
    ICbs=RelSol.MIXMIXbs;
elseif isfield(RelSol,'MIXCLAbs')
    ICbs=RelSol.MIXCLAbs;
else
    error('FSDA:carbikeplotGPCM:WrongInput','A field of input structure RelSol must be ''MIXMIXbs'' , ''MIXCLAbs'',  ''CLACLAbs''')
end

cc=unique([ccshw(:); ccdet(:)]);
lcc=length(cc);
  

ngroups=max(kk);

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);

maxcOralpha=lcc;
axis([0   maxcOralpha+1 0 ngroups+1]);
grid off;
set(axes1,'YTick',0:ngroups,'XTick',1:(lcc+2));
    xlabel('Restriction factor c_{det} or c_{shw}')
ylabel('Number of groups k')


a=cell(lcc,1);
    a(:)={'c='};
legstr=strcat(a, cellstr(num2str(cc)));
legstr=strrep(legstr,' ','');
legstr=[legstr; {' '}];
set(axes1,'XtickLabel',legstr)

numsol=size(ICbs,1);

area = zeros(1,numsol*2);
hr   = zeros(1,numsol*2);
ij=1;
for i=1:numsol
    
    if strcmp(ICbs{i,5},'true') || SpuriousSolutions == true
        kbest=ICbs{i,1};
        cdetBest=find(cc==ICbs{i,2});
        
        if isempty(ICbs{i,3})
            minindcdet=cdetBest;
            maxindcdet=cdetBest;
        else
            minindcdet=find(cc==min(ICbs{i,3}));
            maxindcdet=find(cc==max(ICbs{i,3}));
        end
        
        if isempty(ICbs{i,4})
            minindstablecdet=cdetBest;
            maxindstablecdet=cdetBest;
        else
            minindstablecdet=find(cc==min(ICbs{i,4}));
            maxindstablecdet=find(cc==max(ICbs{i,4}));
        end
        
        area(ij) = ((maxindcdet-minindcdet)*(0.5*(1- i/numsol))) / (numsol*maxindcdet);
        % rectangle
        hrect=0.25*(1- i/numsol)+eps;
        hr(ij)   = rectangle('position',[minindcdet kbest maxindcdet-minindcdet+eps hrect],'facecolor','c','Curvature',[0.2 0.2]);
        % circle associated to the best c
        rectangle('position',[cdetBest-0/2 kbest 0.5/2 0.5/2],'facecolor','w','Curvature',[1 1])
        % line associated with the values of stable c
        minl=min([minindcdet minindstablecdet]);
        rectangle('position',[minl kbest max([maxindcdet maxindstablecdet])-minl+eps eps],'facecolor','w');
        
        text(cdetBest,kbest,[' ' num2str(i)],'HorizontalAlignment','left','FontSize',13,'VerticalAlignment','bottom');
        
        % end of part referred to cdet
        % beginning of part referred to cshw
        cshwBest=find(cc==ICbs{i,6});
        
        if isempty(ICbs{i,7})
            minindcshw=cshwBest;
            maxindcshw=cshwBest;
        else
            minindcshw=find(cc==min(ICbs{i,7}));
            maxindcshw=find(cc==max(ICbs{i,7}));
        end
        
        if isempty(ICbs{i,8})
            minindstablecshw=cshwBest;
            maxindstablecshw=cshwBest;
        else
            minindstablecshw=find(cc==min(ICbs{i,8}));
            maxindstablecshw=find(cc==max(ICbs{i,8}));
        end
        ij=ij+1;
        vdisp=0.1;
        area(ij) = ((maxindcshw-minindcshw)*(0.5*(1- i/numsol))) / (numsol*maxindcshw);
        % rectangle
        hr(ij)   = rectangle('position',[minindcshw kbest+hrect+vdisp maxindcshw-minindcshw+eps 0.25*(1- i/numsol)+eps],'facecolor','g','Curvature',[0.2 0.2]);
        % circle associated to the best c
        rectangle('position',[cshwBest-0/4 kbest+hrect+vdisp 0.5/2 0.5/2],'facecolor','w','Curvature',[1 1])
        minl=min([minindcshw minindstablecshw]);
        % lines out of the car
        rectangle('position',[minl kbest+hrect+vdisp max([maxindcshw maxindstablecshw])-minl+eps eps],'facecolor','w');
        ij=ij+1;
        
        
        % soltruen=sum(strcmp(ICbs(1:i,5),'true'));
        text(cshwBest,kbest+hrect+vdisp,[' ' num2str(i)],'HorizontalAlignment','left','FontSize',13,'VerticalAlignment','bottom');
        
        % text(cshwORalphabest,kbest+hrect*2+vdisp,num2str(i),'HorizontalAlignment','center','FontSize',15,'VerticalAlignment','middle');
    else
        ij=ij+2;
    end
end
hold('on')
b1=bar(2,0,'cyan');
b2=bar(3,0,'green');
legend([b1 b2],'c_{det}','c_{shw}','Location','southeast')

% A = rescaleFS(nanmean(abs(area),1),1,0);
% ivalid = find(area>0);
% colormapres = num2cell(colormap([zeros(numel(ivalid),1) , A(ivalid)' , ones(numel(ivalid),1)]),2);
% set(hr(ivalid),{'facecolor'},colormapres);
set(gca,'FontSize',16)
box('on');

h=gcf;
end
%FScategory:VIS-Clu
