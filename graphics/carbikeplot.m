function h  = carbikeplot(RelSol,varargin)
%carbikeplot produces the carbike plot to find best relevant clustering solutions
%
%<a href="matlab: docsearchFS('carbikeplot')">Link to the help function</a>
%
%   carbikeplot takes as input the output of function tclustICsol (that is
%   a structure containing the best relevant solutions) and produces the
%   car-bike plot. This plot provides a concise summary of the best
%   relevant solutions. This plot shows on the horizontal axis the value of
%   $c$ and on the vertical axis the value of $k$. For each solution we
%   draw a rectangle for the interval of values for which the solution is
%   best and stable and a horizontal line which departs from the rectangle
%   for the values of $c$ in which the solution is only stable. Finally,
%   for the best value of $c$ associated to the solution, we show a circle
%   with two numbers, the first number indicates the ranked solution among
%   those which are not spurious and the second one the ranked number
%   including the spurious solutions. This plot has been baptized
%   ``car-bike'', because the first best solutions (in general 2 or 3) are
%   generally best and stable for a large number of values of $c$ and
%   therefore will have large rectangles. In addition, these solutions are
%   likely to be stable for additional values of $c$ and therefore are
%   likely to have horizontal lines departing from the rectangles (from
%   here the name ``cars''). Finally, local minor solutions (which are
%   associated with particular values of $c$ and $k$) do not generally
%   present rectangles or lines and are shown with circles (from here the
%   name ``bikes'')
%
%  Required input arguments:
%
%           RelSol : Relevant solutions produced by function tclustICsol. Structure.
%                It contains the following fields:
%
%   RelSol.MIXMIXbs = cell of size NumberOfBestSolutions-times-5 which contains
%                the details of the best solutions for MIXMIX (BIC).
%                Each row refers to a solution.  The information which is
%                stored in the columns is as follows.
%                1st col = scalar, value of k for which solution takes place;
%                2nd col = scalar, value of c for which solution takes place;
%                3rd col = row vector of length d which contains the values
%                   of c for which the solution is uniformly better.
%                4th col = row vector of length d+r which contains the
%                   values of c for which the solution is considered stable
%                   (i.e. for which the value of the adjusted Rand index
%                   (or the adjusted Fowlkes and Mallows index)
%                   does not go below the threshold defined in input option
%                   ThreshRandIndex).
%                5th col = string which contains 'true' or 'spurious'. The
%                   solution is labelled spurious if the value of the
%                   adjusted Rand index with the previous solutions is
%                   greater than ThreshRandIndex.
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
% See also: tclustICsol, tclustIC, tclust, tclustregIC, tclustreg, 
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
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('carbikeplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Car-bike plot for simulated data.
    % Generate the data
    restrfact=5;
    rng('default') % Reinitialize the random number generator to its startup configuration
    rng(20000);
    ktrue=3;
    % n = number of observations
    n=150;
    % v= number of dimensions
    v=2;
    % Imposed average overlap
    BarOmega=0.04;

    out=MixSim(ktrue,v,'BarOmega',BarOmega, 'restrfactor',restrfact);
    % data generation given centroids and cov matrices
    [Y,id]=simdataset(n, out.Pi, out.Mu, out.S);

    nsamp=1000;

    % Computation of information criterion
    out=tclustIC(Y,'cleanpool',false,'plots',0,'nsamp',nsamp);

    % Computation of the best solutions
    % Plot first 3 best solutions using as Information criterion CLACLA
    disp('Best 9 solutions using CLACLA')
    ThreshRandIndex=0.8;
    NumberOfBestSolutions=9;
    [outCLACLA]=tclustICsol(out,'whichIC','CLACLA','plots',0,'NumberOfBestSolutions',NumberOfBestSolutions,'ThreshRandIndex',ThreshRandIndex);
    % Car-bike plot to show what are the most relevant solutions
    carbikeplot(outCLACLA);

%}

%{
    %% car-bike plot for the geyser data.
    Y=load('geyser2.txt');
    out=tclustIC(Y,'cleanpool',false,'plots',0,'alpha',0.1);

    % Find the best solutions using as Information criterion MIXMIX
    disp('Best solutions using MIXMIX')
    [outMIXMIX]=tclustICsol(out,'whichIC','MIXMIX','plots',0,'NumberOfBestSolutions',6);
    % Produce the car-bike plot
    carbikeplot(outMIXMIX);
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
% Extract the values of c (number of components)
cc=RelSol.cc;

if isfield(RelSol,'CLACLAbs')
    ICbs=RelSol.CLACLAbs;
elseif isfield(RelSol,'MIXMIXbs')
    ICbs=RelSol.MIXMIXbs;
elseif isfield(RelSol,'MIXCLAbs')
    ICbs=RelSol.MIXCLAbs;
else
    error('FSDA:carbikeplot:WrongInput','A field of input structure RelSol must be ''MIXMIXbs'' , ''MIXCLAbs'',  ''CLACLAbs''')
end


lcc=length(cc);
ngroups=max(kk);

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);

maxc=length(cc);
axis([0   maxc+1 0 ngroups+1]);
grid off;
set(axes1,'YTick',0:ngroups,'XTick',1:(lcc+2));
xlabel('Restriction factor c')
ylabel('Number of groups k')


a=cell(lcc,1);
a(:)={'c='};
if isrow(cc)
    legstr=strcat(a, cellstr(num2str(cc')));
else
    legstr=strcat(a, cellstr(num2str(cc)));
end
legstr=strrep(legstr,' ','');
legstr=[legstr; {' '}];
set(axes1,'XtickLabel',legstr)

numsol=size(ICbs,1);

for i=1:numsol
    if strcmp(ICbs{i,end},'true') || SpuriousSolutions == true
        kbest=ICbs{i,1};
        cbest=find(cc==ICbs{i,2});
        
        if isempty(ICbs{i,3})
            minindc=cbest;
            maxindc=cbest;
        else
            minindc=find(cc==min(ICbs{i,3}));
            maxindc=find(cc==max(ICbs{i,3}));
        end
        
        if isempty(ICbs{i,4})
            minindstablec=cbest;
            maxindstablec=cbest;
        else
            minindstablec=find(cc==min(ICbs{i,4}));
            maxindstablec=find(cc==max(ICbs{i,4}));
        end
        
        rectangle('position',[minindc kbest maxindc-minindc+eps 0.5*(1- i/numsol)+eps],'facecolor','c');
        rectangle('position',[cbest-0.25 kbest 0.5 0.5],'facecolor','w','Curvature',[1 1])
        minl=min([minindc minindstablec]);
        rectangle('position',[minl kbest max([maxindc maxindstablec])-minl+eps eps],'facecolor','w');
        
        soltruen=sum(strcmp(ICbs(1:i,end),'true'));
        
        text(cbest,kbest+0.25,[num2str(soltruen) ',' num2str(i)],'HorizontalAlignment','center','FontSize',16,'VerticalAlignment','middle');
    end
end
set(gca,'FontSize',16)
h=gcf;
end
%FScategory:VIS-Clu
