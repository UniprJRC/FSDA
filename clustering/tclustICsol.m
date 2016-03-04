function out  = tclustICsol(IC,varargin)
%tclustICsol extracts a set of best relevant solutions
%
%<a href="matlab: docsearchFS('tclustICsol')">Link to the help function</a>
%
%   tclustICsol takes as input the output of function tclustIC (that is a
%   series of matrices which contain the values of the information criteria
%   BIC/ICL/CLA for different values of k and c) and extracts the first
%   best solutions. Two solutions are considered equivalent if the value of
%   the adjusted Rand index is above a certain threshold. For each
%   tentative solution the program checks the adjacent values of c for
%   which the solution is stable. A matrix with adjusted Rand indexes is
%   given for the extracted solutions.
%
%  Required input arguments:
%
%           IC : Information criterion to use. Structure. It contains
%                the following fields.
%                IC.CLACLA = matrix of size length(kk)-times length(cc)
%                   containinig the values of the penalized
%                   classification likelihood (CLA).
%                   This field is linked with out.IDXCLA.
%                IC.IDXCLA = cell of size length(kk)-times length(cc).
%                   Each element of the cell is a vector of length n
%                   containinig the assignment of each unit using the
%                   classification model.
%                Remark: fields CLACLA and IDXCLA are linked together.
%                   CLACLA and IDXCLA are compulsory just if optional input
%                   argument 'whichIC' is 'CLACLA' or 'ALL'.
%                IC.MIXMIX = matrix of size length(kk)-times
%                   length(cc) containinig the value of the penalized
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
%                   does not go below the threshold defined in input option
%                   ThreshRandIndex).
%                5th col = string which contains 'true' or 'spurious'. The
%                   solution is labelled spurious if the value of the
%                   adjusted Rand index with the previous solutions is
%                   greater than ThreshRandIndex.
%               Remark: field out.MIXMIXbs is present only if input option
%               'whichIC' is 'ALL' or 'whichIC' is 'MIXMIX'.
%
%  out.MIXMIXbsari =  matrix of adjusted Rand indexes associated with the best
%                solutions for MIXMIX. Matrix of size
%                NumberOfBestSolutions-times-NumberOfBestSolutions whose
%                i,j-th entry contains the adjusted Rand index between
%                classification produced by solution i and solution j,
%                $i,j=1, 2, \ldots, NumberOfBestSolutions$.
%               Remark: field out.MIXMIXbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'MIXMIX'.
%
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
%
% See also: tclustIC, tclust
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
% Copyright 2008-2015.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustICsol')">Link to the help function</a>
% Last modified 06-Feb-2015

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
    [outCLACLCA]=tclustICsol(out,'whichIC','CLACLA','plots',1,'NumberOfBestSolutions',3);
    disp(outCLACLCA.CLACLAbs)
%}

%{
    % An example with input option kk.
    Y=load('geyser2.txt');
    out=tclustIC(Y,'cleanpool',false,'plots',1,'alpha',0.1,'whichIC','CLACLA','kk',[2 3 4 6])
    [outCLACLCA]=tclustICsol(out,'whichIC','CLACLA','plots',1,'NumberOfBestSolutions',3);

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
        'NumberOfBestSolutions',NumberOfBestSolutions, 'ThreshRandIndex', ThreshRandIndex,'msg',msg);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustBICsol:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
    
    
    NumberOfBestSolutions=options.NumberOfBestSolutions;
    ThreshRandIndex=options.ThreshRandIndex;
    whichIC=options.whichIC;
    plots=options.plots;
    msg=options.msg;
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
    warning('Supplied string for whichIC is not supported.')
    error('FSDA:tclustICsol:WrongIC','Specified information criterion is not supported: possible values are ''MIXMIX'' , ''MIXCLA'',  ''CLACLA'', ''ALL''')
end


if typeIC>0
    ARIMIX=zeros(lkk,lcc);
    IDXMIX=IC.IDXMIX;
end

if typeIC==0 || typeIC==3
    ARICLA=zeros(lkk,lcc);
    IDXCLA=IC.IDXCLA;
end

for k=1:length(kk) % loop for different values of k (number of groups)
    for j=2:length(cc) % loop through the different values of c
        if typeIC>0
            ARIMIX(k,j)=RandIndexFS(IDXMIX{k,j-1},IDXMIX{k,j});
        end
        if typeIC==0 || typeIC==3
            ARICLA(k,j)=RandIndexFS(IDXCLA{k,j-1},IDXCLA{k,j});
        end
    end
end


if typeIC==2 || typeIC==3
    [MIXMIXbs,MIXMIXbsari]=findBestSolutions(IC.MIXMIX,ARIMIX,IC.IDXMIX,kk,cc,NumberOfBestSolutions,ThreshRandIndex,msg);
    out.MIXMIXbs=MIXMIXbs;
    out.MIXMIXbsari=MIXMIXbsari;
    
    % Store matrix which contains in the columns the details of the
    % classification
    MIXMIXbsIDX=plotBestSolutions(IC.Y,IC.IDXMIX,MIXMIXbs,kk,cc,'MIXMIX',plots);
    out.MIXMIXbsIDX=MIXMIXbsIDX;
end

if typeIC==1 || typeIC==3
    [MIXCLAbs,MIXCLAbsari]=findBestSolutions(IC.MIXCLA,ARIMIX,IC.IDXMIX,kk,cc,NumberOfBestSolutions,ThreshRandIndex,msg);
    out.MIXCLAbs=MIXCLAbs;
    out.MIXCLAbsari=MIXCLAbsari;
    % Store matrix which contains in the columns the details of the
    % classification
    MIXCLAbsIDX=plotBestSolutions(IC.Y,IC.IDXMIX,MIXCLAbs,kk,cc,'MIXCLA',plots);
    out.MIXCLAbsIDX=MIXCLAbsIDX;
    
end

if typeIC==0 || typeIC==3
    [CLACLAbs,CLACLAbsari]=findBestSolutions(IC.CLACLA,ARICLA,IC.IDXCLA,kk,cc,NumberOfBestSolutions,ThreshRandIndex,msg);
    out.CLACLAbs=CLACLAbs;
    out.CLACLAbsari=CLACLAbsari;
    % Store matrix which contains in the columns the details of the
    % classification
    CLACLAbsIDX=plotBestSolutions(IC.Y,IC.IDXCLA,CLACLAbs,kk,cc,'CLACLA',plots);
    out.CLACLAbsIDX=CLACLAbsIDX;
end

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
ARI(2:end,1)=1;
Xcmod=PENloglik';
seqcc=1:length(cc);
seqkk=1:length(kk);
Bestsols=cell(NumberOfBestSolutions,5);
Bestsols{1,5}='true';
endofloop=0;
for z=1:NumberOfBestSolutions
    [valmin,indmin]=min(Xcmod);
    
    % indminall identifies for which k there is the optimal solution
    [~,minBICk]=min(valmin);
    
    minBICc=indmin(minBICk);
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
    % Store value of k
    Bestsols{z,1}=kk(minBICk);
    % Store value of c
    Bestsols{z,2}=cc(minBICc);
    if msg==1
        disp('Find for which values of c best solution extends to')
    end
    for m=1:1000
        [valmin,indmin]=min(Xcmod);
        if min(valmin)<Inf
            % indminall identifies for which k there is the optimal solution
            [~,minBICknew]=min(valmin);
            % disp(['k=' num2str(kk(minBICknew))])
            minBICcnew=indmin(minBICknew);
            twoc=[minBICcnew minBICc];
            
            if minBICknew==minBICk && min(ARI(min(twoc)+1:max(twoc),minBICk))>ThreshRandIndex
                if msg ==1
                    disp(['c=' num2str(cc(minBICcnew))])
                end
                Xcmod(minBICcnew,minBICk)=Inf;
                % replace old minBICc with minBICcnew
                cwithbestsol(minBICcnew)=1;
                minBICc=minBICcnew;
                % as concerns k, it is the same therefore there is nothing to do
            else
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
                        % those assocaited with best solutions found so far. For
                        % example, if candc =(1 2) it means that we check whether for the same value of k
                        % solution with cc(2) and cc(1) have a R index greater than a certain threshold.
                        % If it is the case this means that these solutions do not have to be considered anymore
                        candcbelow=seqcc(seqcc<min(intc));
                        if ~isempty(candcbelow)
                            posctoadd=min(intc);
                            for r=1:length(candcbelow)
                                if ARI(posctoadd,minBICk)>ThreshRandIndex
                                    Xcmod(posctoadd-1,minBICk)=Inf;
                                    cctoadd(posctoadd-1)=1;
                                    posctoadd=posctoadd-1;
                                else
                                    break
                                end
                            end
                        end
                        
                        candcabove=seqcc(seqcc>max(intc));
                        if ~isempty(candcabove)
                            posctoadd=max(intc)+1;
                            for r=1:length(candcabove)
                                if ARI(posctoadd,minBICk)>ThreshRandIndex
                                    Xcmod(posctoadd,minBICk)=Inf;
                                    cctoadd(posctoadd)=1;
                                    posctoadd=posctoadd+1;
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
                    for  j=1:z-1
                        idxpreviousz=IDX{seqkk(kk==Bestsols{j,1}),seqcc(cc==Bestsols{j,2})};
                        idxcurrentz=IDX{seqkk(kk==Bestsols{z,1}),seqcc(cc==Bestsols{z,2})};
                        if RandIndexFS(idxpreviousz,idxcurrentz)>ThreshRandIndex
                            Bestsols{z,5}='spurious';
                            break
                        else
                            Bestsols{z,5}='true';
                        end
                    end
                end
                
                break
            end
            
            % Check if the solution which has just been found has a
            % Rand index greater than a certain threshold with those which have
            % already been found
            if z>1
                for  j=1:z-1
                    idxpreviousz=IDX{seqkk(kk==Bestsols{j,1}),seqcc(cc==Bestsols{j,2})};
                    idxcurrentz=IDX{seqkk(kk==Bestsols{z,1}),seqcc(cc==Bestsols{z,2})};
                    if RandIndexFS(idxpreviousz,idxcurrentz)>ThreshRandIndex
                        Bestsols{z,5}='spurious';
                        break
                    else
                        Bestsols{z,5}='true';
                    end
                end
            end
        else % if Xcmod is full of inf than get out of the loop
            Bestsols=Bestsols(1:z,:);
            endofloop=1;
            break
        end
    end
    if endofloop==1
        if msg==1
            disp(['There are at most ' num2str(z) ' different solutions'])
        end
        break
    end
end

%% Find matrix of ARI for the z solutions which have been found
ARIbest=zeros(z,z);
for i=1:z
    for  j=1:z
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
cascade

end
