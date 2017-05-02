function out  = tclustICsol(IC,varargin)
%tclustICsol extracts a set of best relevant solutions
%
%<a href="matlab: docsearchFS('tclustICsol')">Link to the help function</a>
%
%   tclustICsol takes as input the output of function tclustIC (that is a
%   series of matrices which contain the values of the information criteria
%   BIC/ICL/CLA for different values of k and c) and extracts the first
%   best solutions. Two solutions are considered equivalent if the value of
%   the adjusted Rand index (or the adjusted Fowlkes and Mallows index) is
%   above a certain threshold. For each tentative solution the program
%   checks the adjacent values of c for which the solution is stable. A
%   matrix with adjusted Rand indexes is given for the extracted solutions.
%
%  Required input arguments:
%
%           IC : Information criterion to use. Structure. It contains
%                the following fields.
%                IC.CLACLA = matrix of size length(kk)-by-length(cc)
%                   containinig the values of the penalized
%                   classification likelihood (CLA).
%                   This field is linked with out.IDXCLA.
%                IC.IDXCLA = cell of size length(kk)-by-length(cc).
%                   Each element of the cell is a vector of length n
%                   containinig the assignment of each unit using the
%                   classification model.
%                Remark: fields CLACLA and IDXCLA are linked together.
%                   CLACLA and IDXCLA are compulsory just if optional input
%                   argument 'whichIC' is 'CLACLA' or 'ALL'.
%                IC.MIXMIX = matrix of size length(kk)-by-length(cc)
%                   containinig the value of the penalized
%                   mixture likelihood (BIC). This field is linked with
%                   out.IDXMIX.
%                IC.MIXCLA = matrix of size length(kk)-times length(cc)
%                   containinig the value of the ICL. This field is linked
%                   with out.IDXMIX.
%                IC.IDXMIX = cell of size length(kk)-times length(cc).
%                   Each element of the cell is a vector of length n
%                   containinig the assignment of each unit using the
%                   mixture model.
%                Remark 1: fields MIXMIX and IDXMIX are linked together.
%                   MIXMIX and IDXMIX are compulsory just if optional input
%                   argument 'whichIC' is 'CLACLA' or 'ALL'.
%                Remark 2: fields MIXCLA and IDXMIX are linked together.
%                   MIXCLA and IDXMIX are compulsory just if optional input
%                   argument 'whichIC' is 'MIXCLA' or 'ALL'.
%                IC.kk = vector containing the values of k (number of
%                   components) which have been considered.
%                IC.cc = vector containing the values of c (values of the
%                   restriction factor) which have been considered.
%                IC.Y =  original n-times-v data matrix on which the IC
%                   (Information criterion) has
%                    been computed
%                 Data Types - struct
%
%  Optional input arguments:
%
%
% NumberOfBestSolutions: number of solutions to consider. Scalar integer
%                       greater than 0. Number of best solutions to
%                       extract from BIC/ICL matrix. The default value
%                       of NumberOfBestSolutions is 5
%                       Example - 'NumberOfBestSolutions',5
%                       Data Types - int16 | int32 | single | double
%
% ThreshRandIndex     : threshold to identify spurious solutions. Positive
%                       scalar between 0 and 1. Scalar which specifies the
%                       threshold of the adjusted Rnd index to use to
%                       consider two solutions as equivalent. The default
%                       value of ThreshRandIndex is 0.7
%                       Example - 'ThreshRandIndex',0.8
%                       Data Types - single | double
%
%   whichIC  : character which specifies the information criterion to use
%               to extract best solutions. Character.
%               Possible values for whichIC are:
%               'CLACLA' = in this case best solutions are referred to
%                   the classification likelihood.
%               'MIXMIX' = in this case in this case best solutions are
%                   referred to the mixture likelihood (BIC).
%               'MIXCLA' = in this case in this case best solutions are
%                   referred to ICL.
%                 'ALL'  = in this case best solutions both three solutions
%                          using classification and mixture likelihood are
%                          produced. In output structure out all the three
%                          matrices out.MIXMIXbs, out.CLACLAbs and
%                          out.MIXCLAbs are given.
%               The default value of 'whichIC' is 'ALL'
%                 Example - 'whichIC','ALL'
%                 Data Types - character
%
%       plots : plots of best solutions on the screen. Scalar. It specifies
%               whether to plot on the screen the best solutions which have
%               been found.
%                 Example - 'plots',1
%                 Data Types - single | double
%
%       msg  :  Message on the screen. Scalar. Scalar which controls
%               whether to display or not messages about code execution.
%               The default value of msg is 0, that is no message is
%               displayed on the screen.
%                 Example - 'msg',1
%                 Data Types - single | double
%
%     Rand   :  Index to use to compare partitions. Scalar. If Rand =1
%               (default) the adjusted Rand index is used, else the
%               adjusted Fowlkes and Mallows index is used
%                 Example - 'Rand',1
%                 Data Types - single | double
%
%  Output:
%
%         out:   structure which contains the following fields:
%
%   out.MIXMIXbs = cell of size NumberOfBestSolutions-times-5 which contains
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
%  out.MIXMIXbsari =  matrix of adjusted Rand indexes (or Fowlkes and Mallows
%               indexes) associated with the best
%                solutions for MIXMIX. Matrix of size
%                NumberOfBestSolutions-times-NumberOfBestSolutions whose
%                i,j-th entry contains the adjusted Rand index between
%                classification produced by solution i and solution j,
%                $i,j=1, 2, \ldots, NumberOfBestSolutions$.
%               Remark: field out.MIXMIXbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'MIXMIX'.
%
% out.ARIMIX = Matrix of adjusted Rand indexes between two consecutive value of c.
%                 Matrix of size k-by-length(cc)-1. The first column
%                 contains the ARI indexes between 
%                 with cc(2) and cc(1) given k. The second column contains
%                 the the ARI indexes  between cc(3) and cc(2) given k.
%                 This output is also present in table format (see below)
%               Remark: field ARIMIX is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'MIXMIX' or 'MIXLCA'
% out.ARIMIXtable = Table with the same meaning of matrix ARIMIX above.
%                 A Matlab table has also been been given to faciliate the
%                 interpretation of the rows and columns. The Rownames of
%                 this table correspond to the values of k which are used
%                 and the colNames of this table contain in a dynamic way
%                 the two values of c which are considered. For example if
%                 the first two values of c are c=3 and c=7, the first
%                 column name of this table is c3_v_c7 to denote that the
%                 entry of this column are the ARI indexes between c=3 and
%                 c=7
%               Remark: field ARIMIXtable is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'MIXMIX' or 'MIXLCA'
%  out.MIXCLAbs = this output has the same structure as out.MIXMIXbs but
%               it is referred to MIXCLA.
%               Remark: field out.MIXCLAbs is present only if 'whichIC' is
%               'ALL' or 'whichIC' is 'MIXCLA'.
%
% out.MIXCLAbsari = this output has the same structure as out.MIXMIXbs but
%               it is referred to MIXCLA.
%               Remark: field out.MIXCLAbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'MIXCLA'.
%
%  out.CLACLAbs = this output has the same structure as out.MIXMIXbs but
%               it is referred to CLACLA.
%               Remark: field out.CLACLAbs is present only if 'whichIC' is
%               'ALL' or 'whichIC' is 'CLACLA'.
%
% out.CLACLAbsari = this output has the same structure as out.MIXMIXbs but
%               it is referred to CLACLA.
%               Remark: field out.MIXCLAbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'CLACLA'
%
% out.ARICLA = Matrix of adjusted Rand indexes between two consecutive value of c.
%                 Matrix of size k-by-length(cc)-1. The first column
%                 contains the ARI indexes between 
%                 with cc(2) and cc(1) given k. The second column contains
%                 the the ARI indexes  between cc(3) and cc(2) given k.
%                 This output is also present in table format (see below)
%               Remark: field ARICLA is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'CLACLA'
% out.ARICLAtable = Table with the same meaning of matrixo CLACLAari above.
%                 A Matlab table has also been been given to faciliate the
%                 interpretation of the rows and columns. The Rownames of
%                 this table correspond to the values of k which are used
%                 and the colNames of this table contain in a dynamic way
%                 the two values of c which are considered. For example if
%                 the first two values of c are c=3 and c=7, the first
%                 column name of this table is c3_v_c7 to denote that the
%                 entry of this column are the ARI indexes between c=3 and
%                 c=7
%               Remark: field ARICLAtable is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'CLACLA'
%           out.kk = vector containing the values of k (number of
%                   components) which have been considered. This  vector
%                   is equal to input optional argument kk if kk had been
%                   specified else it is equal to 1:5.
%
%          out.cc = vector containing the values of c (values of the
%                   restriction factor) which have been considered. This
%                   vector is equal to input optional argument cc if cc had
%                   been specified else it is equal to [1, 2, 4, 8, 16, 32,
%                   64, 128].
%
% See also: tclustIC, tclust, carbikeplot
%
% References:
%
% A. Cerioli, L.A. Garcia-Escudero, A. Mayo-Iscar and M. Riani (2016).
% Finding the Number of Groups in Model-Based Clustering via Constrained
% Likelihoods, submitted.
%
% L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of
% Classification 2:193-218
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustICsol')">Link to the help function</a>
% Last modified 31-05-2016

% Examples:

%{
    %% Plot of first two best solutions for Geyser data.
    Y=load('geyser2.txt');
    out=tclustIC(Y,'cleanpool',false,'plots',0,'alpha',0.1);

    % Plot first two best solutions using as Information criterion MIXMIX
    disp('Best solutions using MIXMIX')
    [outMIXMIX]=tclustICsol(out,'whichIC','MIXMIX','plots',1,'NumberOfBestSolutions',2);
    disp(outMIXMIX.MIXMIXbs)
%}

%{
    %% Simulated data: compare first 3 best solutions using MIXMIX and CLACLA.
    % Data generation
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

    % Computation of information criterion
    out=tclustIC(Y,'cleanpool',false,'plots',0,'nsamp',200);
    % Plot first 3 best solutions using as Information criterion MIXMIX
    disp('Best 3 solutions using MIXMIX')
    [outMIXMIX]=tclustICsol(out,'whichIC','MIXMIX','plots',1,'NumberOfBestSolutions',3);
    disp(outMIXMIX.MIXMIXbs)
    disp('Best 3 solutions using CLACLA')
    [outCLACLA]=tclustICsol(out,'whichIC','CLACLA','plots',1,'NumberOfBestSolutions',3);
    disp(outCLACLA.CLACLAbs)
%}



%{
    % An example with input option kk.
    Y=load('geyser2.txt');
    out=tclustIC(Y,'cleanpool',false,'plots',1,'alpha',0.1,'whichIC','CLACLA','kk',[2 3 4 6])
    [outCLACLCA]=tclustICsol(out,'whichIC','CLACLA','plots',1,'NumberOfBestSolutions',3);

%}

%{
    % Comparison between the use of Rand index and FM index.
    Y=load('geyser2.txt');
    out=tclustIC(Y,'cleanpool',false,'plots',1,'alpha',0.1,'whichIC','CLACLA')
    [outCLACLCA]=tclustICsol(out,'whichIC','CLACLA','plots',0,'NumberOfBestSolutions',5,'Rand',1);
    disp('Matrix of adjusted Rand indexes among the first 5 best solutions')
    disp(outCLACLCA.CLACLAbsari)
    [outCLACLCA]=tclustICsol(out,'whichIC','CLACLA','plots',0,'NumberOfBestSolutions',5,'Rand',0);
    disp('Matrix of adjusted Fowlkes and Mallows indexes among the first 5 best solutions')
    disp(outCLACLCA.CLACLAbsari)

%}



%% Beginning of code

% Default number of best solutions to consider
NumberOfBestSolutions=5;
% Specify information criterion for which you want to see best solutions
whichIC='ALL';
% Default value of Adjust Rand index above which two solutions are
% considered as equivalent
ThreshRandIndex=0.7;
% Plots on the screen
plots=1;
% Message about code execution
msg=0;

if nargin>1
    options=struct('whichIC',whichIC,'plots',plots,...
        'NumberOfBestSolutions',NumberOfBestSolutions, 'ThreshRandIndex', ThreshRandIndex,'msg',msg,'Rand',1);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustICsol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    
    
    NumberOfBestSolutions=options.NumberOfBestSolutions;
    ThreshRandIndex=options.ThreshRandIndex;
    whichIC=options.whichIC;
    plots=options.plots;
    msg=options.msg;
    Rand=options.Rand;
end


% Extract the values of k (number of groups)
kk=IC.kk;
lkk=length(kk);
% Extract the values of c (number of components)
cc=IC.cc;
lcc=length(cc);

if strcmp(whichIC,'ALL')
    typeIC=3;
elseif strcmp(whichIC,'MIXMIX')
    typeIC=2;
elseif strcmp(whichIC,'MIXCLA')
    typeIC=1;
elseif strcmp(whichIC,'CLACLA')
    typeIC=0;
else
    warning('FSDA:tclustICsol:WrongIC','Supplied string for whichIC is not supported.')
    error('FSDA:tclustICsol:WrongIC','Specified information criterion is not supported: possible values are ''MIXMIX'' , ''MIXCLA'',  ''CLACLA'', ''ALL''')
end


if typeIC>0
    ARIMIX=ones(lkk,lcc);
    IDXMIX=IC.IDXMIX;
end

if typeIC==0 || typeIC==3
    ARICLA=ones(lkk,lcc);
    
    IDXCLA=IC.IDXCLA;
end

for k=1:length(kk) % loop for different values of k (number of groups)
    if kk(k) ~= 1 % just do the computations or ARI or FM indexes when k is different from 1
        for j=2:length(cc) % loop through the different values of c
            if typeIC>0
                if Rand==1
                    ARIMIX(k,j)=RandIndexFS(IDXMIX{k,j-1},IDXMIX{k,j});
                else
                    ARIMIX(k,j)=FowlkesMallowsIndex(IDXMIX{k,j-1},IDXMIX{k,j});
                end
                
            end
            if typeIC==0 || typeIC==3
                if Rand==1
                    ARICLA(k,j)=RandIndexFS(IDXCLA{k,j-1},IDXCLA{k,j});
                else
                    ARICLA(k,j)=FowlkesMallowsIndex(IDXCLA{k,j-1},IDXCLA{k,j});
                end
                
            end
        end
    end
end

% Prepare rownames and colsnames for table which will contain
% in the rows the number of groups from
rownamesARI=strcat(cellstr(repmat('k=',length(kk),1)), cellstr(num2str(kk')));

lc1=length(cc)-1;
cup=strcat(cellstr(repmat('c',lc1,1)), cellstr(num2str(cc(2:end)')));
clow=strcat(cellstr(repmat('c',lc1,1)), cellstr(num2str(cc(1:end-1)')));
cuplow=strcat(cup,repmat('_vs_',lc1,1),clow);
colnamesARI=regexprep(cuplow,' ','');



out=struct;

if typeIC==2 || typeIC==3
    [MIXMIXbs,MIXMIXbsari]=findBestSolutions(IC.MIXMIX,ARIMIX,IC.IDXMIX,kk,cc,NumberOfBestSolutions,ThreshRandIndex,msg);
    out.MIXMIXbs=MIXMIXbs;
    out.MIXMIXbsari=MIXMIXbsari;
    
    % Store matrix which contains in the columns the details of the
    % classification
    MIXMIXbsIDX=plotBestSolutions(IC.Y,IC.IDXMIX,MIXMIXbs,kk,cc,'MIXMIX',plots);
    out.MIXMIXbsIDX=MIXMIXbsIDX;
    out.ARIMIX=ARIMIX(:,2:end);
    ARIMIXtable=array2table(ARIMIX(:,2:end),'RowNames',rownamesARI,'VariableNames',colnamesARI);
    out.ARIMIXtable=ARIMIXtable;
    
end

if typeIC==1 || typeIC==3
    [MIXCLAbs,MIXCLAbsari]=findBestSolutions(IC.MIXCLA,ARIMIX,IC.IDXMIX,kk,cc,NumberOfBestSolutions,ThreshRandIndex,msg);
    out.MIXCLAbs=MIXCLAbs;
    out.MIXCLAbsari=MIXCLAbsari;
    % Store matrix which contains in the columns the details of the
    % classification
    MIXCLAbsIDX=plotBestSolutions(IC.Y,IC.IDXMIX,MIXCLAbs,kk,cc,'MIXCLA',plots);
    out.MIXCLAbsIDX=MIXCLAbsIDX;
    out.ARIMIX=ARIMIX(:,2:end);
    ARIMIXtable=array2table(ARIMIX(:,2:end),'RowNames',rownamesARI,'VariableNames',colnamesARI);
    out.ARIMIXtable=ARIMIXtable;
end

if typeIC==0 || typeIC==3
    [CLACLAbs,CLACLAbsari]=findBestSolutions(IC.CLACLA,ARICLA,IC.IDXCLA,kk,cc,NumberOfBestSolutions,ThreshRandIndex,msg);
    out.CLACLAbs=CLACLAbs;
    out.CLACLAbsari=CLACLAbsari;
    % Store matrix which contains in the columns the details of the
    % classification
    CLACLAbsIDX=plotBestSolutions(IC.Y,IC.IDXCLA,CLACLAbs,kk,cc,'CLACLA',plots);
    out.CLACLAbsIDX=CLACLAbsIDX;
    out.ARICLA=ARICLA(:,2:end);
    ARICLAtable=array2table(ARICLA(:,2:end),'RowNames',rownamesARI,'VariableNames',colnamesARI);
    out.ARICLAtable=ARICLAtable;
    
end

% Store values of c and k which have been used
out.cc=cc;
out.kk=kk;

end


function [Bestsols,ARIbest]  = findBestSolutions(PENloglik,ARI,IDX,kk,cc,NumberOfBestSolutions,ThreshRandIndex,msg)
% PENloglik = matrix of size length(kk)-times-length(cc) containing penalized log likelihood
% ARI = matrix of size length(kk)-times-length(cc) containing Adjusted R
%       indexes between two consecutive values of c for fixed k
% IDX = cell of size length(kk)-times-length(cc) containing the indication
%       where units have been classified for each k (rows) and c (columns)
% kk  = vector containing values of k which have to be considered
% cc  = vector containing values of c which have to be considered
ARI=ARI';

% The rows of ARI refer to the values of c. The columns of ARI refer to k

% this instruction should not be necessary anymore
% ARI(2:end,1)=1;

Xcmod=PENloglik';
lcc=length(cc);
seqcc=1:lcc;
seqkk=1:length(kk);
Bestsols=cell(NumberOfBestSolutions,5);
Bestsols{1,5}='true';
endofloop=0;
NumberOfExistingSolutions=NumberOfBestSolutions;
for z=1:NumberOfBestSolutions
    
    % valmin= mimimum for IC each value of k
    [valmin,indmin]=min(Xcmod,[],1);
    
    % minBICk identifies position of k where there is the optimal solution
    [~,minBICk]=min(valmin);
    % minBICc position of the best value of c where there is the optimal
    % solution
    minBICc=indmin(minBICk);
    
    if min(valmin)<Inf
        
        if msg==1
            disp(['Best solution number: ' num2str(z)])
            disp(['k=' num2str(kk(minBICk))])
            disp(['c=' num2str(cc(minBICc))])
        end
        % overall minimum is in col minBICk
        % overall minimum is in row minBICc
        Xcmod(indmin(minBICk),minBICk)=Inf;
        cwithbestsol=nan(length(cc),1);
        cwithbestsol(minBICc)=1;
        % Store value of k in column 1 of Bestsols
        Bestsols{z,1}=kk(minBICk);
        % Store value of c in column 2 of Bestsols
        Bestsols{z,2}=cc(minBICc);
        if msg==1
            disp('Find for which adjacent value of c (and fixed k) best solution extends to')
        end
        
        % Find overall minimum of matrix Xcmod after excluding the column
        % associated with the value of k where minBICk lies
        XcmodWithoutBestk=Xcmod;
        XcmodWithoutBestk(:,kk(minBICk))=[];
        % minICconstr= min value of IC excluding the values involving column
        % minBICk
        minICconstr=min(min(XcmodWithoutBestk));
        cctoadd=zeros(length(cc),1);
        
        candcabove=seqcc(seqcc>minBICc);
        if ~isempty(candcabove)
            for r=1:length(candcabove)
                posctoadd=candcabove(r);
                if ARI(posctoadd,minBICk)>ThreshRandIndex && Xcmod(posctoadd,minBICk)<minICconstr
                    Xcmod(posctoadd,minBICk)=Inf;
                    cctoadd(posctoadd)=1;
                else
                    break
                end
            end
        end
        
        candcbelow=seqcc(seqcc<minBICc);
        if ~isempty(candcbelow)
            for r=length(candcbelow):-1:1
                posctoadd=candcbelow(r);
                if ARI(posctoadd+1,minBICk)>ThreshRandIndex && Xcmod(posctoadd,minBICk)<minICconstr
                    Xcmod(posctoadd,minBICk)=Inf;
                    cctoadd(posctoadd)=1;
                else
                    break
                end
            end
        end
        
        
        cwithbestsol(cctoadd==1)=1;
        % Store the values of c associated to the best solution
        Bestsols{z,3}=cc(cwithbestsol==1);
        
        % The interval of values of c for which the solution is uniformly
        % better has been found, however, now we have to make sure that we
        % do not consider anymore the solutions for the same k which have a
        % value of R index adjacent to those which have already been found
        % greater than a certain threshold
        % Remark. If minBICk =1
        
        
        if minBICk ==1
            Bestsols{z,4}=cc;
            Xcmod(:,minBICk)=Inf;
        else
            intc=seqcc(cwithbestsol==1);
            outc=cc(cwithbestsol~=1);
            cctoadd=zeros(length(cc),1);
            
            if ~isempty(outc)
                % candcbelow = indexes of the set of candidate values for c which are smaller than
                % those associated with best solutions found so far. For
                % example, if candc =(1 2) it means that we check whether for the same value of k
                % solution with cc(2) and cc(1) have a R index greater than a certain threshold.
                % If it is the case this means that these solutions do not have to be considered anymore
                candcbelow=seqcc(seqcc<min(intc));
                
                if ~isempty(candcbelow)
                    for r=length(candcbelow):-1:1
                        posctoadd=candcbelow(r);
                        if ARI(posctoadd+1,minBICk)>ThreshRandIndex
                            Xcmod(posctoadd,minBICk)=Inf;
                            cctoadd(posctoadd)=1;
                        else
                            break
                        end
                    end
                end
                
                candcabove=seqcc(seqcc>max(intc));
                if ~isempty(candcabove)
                    for r=1:length(candcabove)
                        posctoadd=candcabove(r);
                        if ARI(posctoadd,minBICk)>ThreshRandIndex
                            Xcmod(posctoadd,minBICk)=Inf;
                            cctoadd(posctoadd)=1;
                        else
                            break
                        end
                    end
                end
                
                Bestsols{z,4}=cc(cctoadd==1);
            end
        end
        % get out of the loop because in this case there is a new candidate
        % solution (with the same k because the R index was smaller than a
        % certain threshold, or with a different k)
        
        % Before getting out of the loop, check if the solution
        % which has just been found, has a Rand index greater than a
        % certain threshold with those which have already been
        % found
        if z>1
            idxcurrentz=IDX{seqkk(kk==Bestsols{z,1}),seqcc(cc==Bestsols{z,2})};
            
            for  j=1:z-1
                idxpreviousz=IDX{seqkk(kk==Bestsols{j,1}),seqcc(cc==Bestsols{j,2})};
                if RandIndexFS(idxpreviousz,idxcurrentz)>ThreshRandIndex  && strcmp(Bestsols{j,5},'true')==1
                    Bestsols{z,5}='spurious';
                    break
                else
                    Bestsols{z,5}='true';
                end
            end
        end
        
        
    else % if Xcmod is full of inf than get out of the loop
        NumberOfExistingSolutions=z-1;
        Bestsols=Bestsols(1:NumberOfExistingSolutions,:);
        endofloop=1;
        break
    end
end

if endofloop==1
    if msg==1
        disp(['There are at most ' num2str(z) ' different solutions'])
    end
    %  break
end


%% Find matrix of ARI for the z solutions which have been found
ARIbest=zeros(NumberOfExistingSolutions,NumberOfExistingSolutions);
for i=1:NumberOfExistingSolutions
    for  j=1:NumberOfExistingSolutions
        idxi=IDX{seqkk(kk==Bestsols{i,1}),seqcc(cc==Bestsols{i,2})};
        idxj=IDX{seqkk(kk==Bestsols{j,1}),seqcc(cc==Bestsols{j,2})};
        ARIbest(i,j)=RandIndexFS(idxi,idxj);
    end
end
end


function IDXout=plotBestSolutions(Y,IDX,Bestsols,kk,cc,lab,plots)
seqkk=1:length(kk);
seqcc=1:length(cc);
nbestsol=size(Bestsols,1);
IDXout=zeros(size(Y,1),nbestsol);
for i=1:nbestsol
    IDXi=IDX{seqkk(kk==Bestsols{i,1}),seqcc(cc==Bestsols{i,2})};
    IDXout(:,i)=IDXi;
    if plots==1
        figure
        spmplot(Y,IDXi,[],'box');
        detsol=[lab ': solution ' num2str(i)  ': k=' num2str(Bestsols{i,1}) ' c=' num2str(Bestsols{i,2})];
        bestsol=[' Best in c ' num2str(min(Bestsols{i,3})) '-' num2str(max(Bestsols{i,3})) ];
        Bestsolc=[Bestsols{i, 3} Bestsols{i,4}];
        stabsol=[' Stable in c ' num2str(min(Bestsolc)) '-' num2str(max(Bestsolc)) ];
        title([detsol bestsol stabsol ' Sol:' Bestsols{i,5}])
    end
end

if plots==1
    cascade
end

end
%FScategory:CLUS-RobClaMULT