function out  = tclustICsolGPCM(IC,varargin)
%tclustICsolGPCM extracts a set of best relevant solutions from 3D array computed using function tclustICgpcm
%
%<a href="matlab: docsearchFS('tclustICsolGPCM')">Link to the help function</a>
%
%   tclustICsolGPCM takes as input the output of function tclustICgpcm that
%   is a series of matrices which contain the values of the information
%   criteria BIC/ICL/CLA for different values of $k$ and $c_{det}$ and
%   $c_{shw}$ (for fixed trimming level $\alpha$) and extracts the first
%   best solutions. Two solutions are considered equivalent if the value of
%   the adjusted Rand index (or the adjusted Fowlkes and Mallows index) is
%   above a certain threshold. For each tentative solution the program
%   checks the adjacent values of $c_{det}$ and $c_{shw}$ for which the
%   solution is stable. A matrix with adjusted Rand indexes is given for
%   the extracted solutions.
%
%  Required input arguments:
%
%           IC : Information criterion to use. Structure. It contains
%                the following fields.
%                IC.CLACLA = 3D array of size
%                   length(kk)-by-length(cdet)-by-length(cshw) containinig
%                   the values of the penalized
%                   classification likelihood (CLA).
%                   This field is linked with IC.IDXCLA.
%                IC.IDXCLA = 3D array of size
%                   length(kk)-by-length(cdet)-by-length(csshw).
%                   Each element of the cell is a vector of length n
%                   containinig the assignment of each unit using the
%                   classification model.
%                Remark: fields CLACLA and IDXCLA are linked together.
%                   CLACLA and IDXCLA are compulsory just if optional input
%                   argument 'whichIC' is 'CLACLA'.
%                IC.MIXMIX = 3D array of size length(kk)-by-length(cdet)-by-length(cshw) 
%                   containinig the value of the penalized
%                   mixture likelihood (BIC). This field is linked with
%                   IC.IDXMIX.
%                IC.MIXCLA = 3D array of size
%                   length(kk)-by-length(cdet)-by-length(cshw) containinig
%                   the value of the ICL. This field is linked with
%                   IC.IDXMIX.
%                IC.IDXMIX = 3D cell of size
%                   length(kk)-by-length(cdet)-by-length(cshw).
%                   Each element of the cell is a vector of length n
%                   containinig the assignment of each unit using the
%                   mixture model.
%                Remark 1: fields MIXMIX and IDXMIX are linked together.
%                   MIXMIX and IDXMIX are compulsory just if optional input
%                   argument 'whichIC' is 'CLACLA'.
%                Remark 2: fields MIXCLA and IDXMIX are linked together.
%                   MIXCLA and IDXMIX are compulsory just if optional input
%                   argument 'whichIC' is 'MIXCLA'.
%                IC.kk = vector containing the values of k (number of
%                   components) which have been considered.
%                IC.ccdet = vector containing the values of cdet (values of the
%                   restriction factor for ratio of determinants) which
%                   have been considered.
%                IC.ccshw = vector containing the values of cshw (values of the
%                   restriction factor for ratio of elements of shape
%                   matrices inside each group) which have been considered.
%                IC.alpha = scalar containing the values of trimming level
%                   which has been considered.
%                IC.Y =  original n-times-v data matrix on which the IC
%                   (Information criterion) has
%                    been computed. This input option is present only if IC
%                    comes from tclustIC.
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
%               'MIXMIX' = in this case best solutions are
%                   referred to the mixture likelihood (BIC).
%               'MIXCLA' = in this case best solutions are
%                   referred to ICL.
%               The default value of 'whichIC' is 'MIXMIX'
%                 Example - 'whichIC','CLACLA'
%                 Data Types - character
%
%       plots : plots of best solutions on the screen. Scalar. It specifies
%               whether to plot on the screen the best solutions which have
%               been found.
%                 Example - 'plots',1
%                 Data Types - single | double
%
% SpuriousSolutions  :  Include or nor spurious solutions in the plot. Boolean.
%                       As default spurios solutions are shown in the plot.
%                 Example - 'SpuriousSolutions',false
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
%   out.MIXMIXbs = cell of size NumberOfBestSolutions-times-8 which contains
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
%               'whichIC' is 'MIXMIX'.
%
%  out.MIXMIXbsari =  matrix of adjusted Rand indexes (or Fowlkes and Mallows
%               indexes) associated with the best
%                solutions for MIXMIX. Matrix of size
%                NumberOfBestSolutions-times-NumberOfBestSolutions whose
%                i,j-th entry contains the adjusted Rand (or Fowlkes and Mallows) index between
%                classification produced by solution i and solution j,
%                $i,j=1, 2, \ldots, NumberOfBestSolutions$.
%               Remark: field out.MIXMIXbsari is present only if 'whichIC'
%               is 'MIXMIX'.
%
%
%  out.MIXCLAbs = this output has the same structure as out.MIXMIXbs but
%               it is referred to MIXCLA.
%               Remark: field out.MIXCLAbs is present only if 'whichIC' is
%               'MIXCLA'.
%
% out.MIXCLAbsari = this output has the same structure as out.MIXMIXbs but
%               it is referred to MIXCLA.
%               Remark: field out.MIXCLAbsari is present only if 'whichIC'
%               is 'MIXCLA'.
%
%  out.CLACLAbs = this output has the same structure as out.MIXMIXbs but
%               it is referred to CLACLA.
%               Remark: field out.CLACLAbs is present only if 'whichIC' is
%               'CLACLA'.
%
% out.CLACLAbsari = this output has the same structure as out.MIXMIXbs but
%               it is referred to CLACLA.
%               Remark: field out.MIXCLAbsari is present only if 'whichIC'
%               is 'ALL' or 'whichIC' is 'CLACLA'
%
%
% out.MIXCLAbsIDX = matrix of dimension n-by-NumberOfBestSolutions
%               containing the allocations for MIXCLA associated with the best
%               NumberOfBestSolutions. This field is present only if 'whichIC'
%               is 'MIXCLA'.
%
% out.MIXMIXbsIDX = matrix of dimension n-by-NumberOfBestSolutions
%               containing the allocations for MIXMIX associated with the best
%               NumberOfBestSolutions. This field is present only if 'whichIC'
%               is 'MIXMIX'.
%
% out.CLACLAbsIDX = matrix of dimension n-by-NumberOfBestSolutions
%               containing the allocations for CLACLA associated with the
%               best NumberOfBestSolutions. This field is present only if
%               'whichIC' is 'CLACLA'.
%
%
%           out.kk = vector containing the values of k (number of
%                   components) which have been considered. This  vector
%                   is equal to input optional argument kk if kk had been
%                   specified else it is equal to 1:5.
%
%          out.ccdet = vector containing the values of cdet (values of the
%                   restriction factor for determinants) which have been
%                   considered. This vector is equal to input argument
%                   IC.cdet.
%
%          out.ccshw = vector containing the values of cshw (values of the
%                   restriction factor for shape elements inside each
%                   group) which have been considered. This vector is equal
%                   to input argument IC.cshw.
%
%          out.alpha = scalar containing the value of $\alpha$ (trimming
%                   level) which have been considered. This
%                   output is equal to input argument IC.alpha.
%
% See also: tclustICgpcm, tclust, carbikeplot
%
% References:
%
% Cerioli, A., Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2017),
% Finding the Number of Groups in Model-Based Clustering via Constrained
% Likelihoods, "Journal of Computational and Graphical Statistics", pp. 404-416,
% https://doi.org/10.1080/10618600.2017.1390469
%
% Hubert L. and Arabie P. (1985), Comparing Partitions, "Journal of
% Classification", Vol. 2, pp. 193-218.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('tclustICsolGPCM')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
   %% Plot of first two best solutions for Geyser data.
    Y=load('geyser2.txt');
    % nsamp=30 to reduce computational time
    outIC=tclustICgpcm(Y,'cleanpool',false,'plots',0,'alpha',0.1,'nsamp',30);

    % Plot first two best solutions using as Information criterion MIXMIX
    disp('Best solutions using MIXMIX')
    [out]=tclustICsolGPCM(outIC,'whichIC','MIXMIX','plots',1,'NumberOfBestSolutions',2);
    disp(out.MIXMIXbs)
%}

%{
    % Simulated data: compare first 2 best solutions using MIXMIX and CLACLA.
    % Data generation
    restrfact=5;
    rng('default') % Reinitialize the random number generator to its startup configuration
    rng(10000);
    ktrue=3;
    % n = number of observations
    n=150;
    % v= number of dimensions
    v=2;
    % Imposed average overlap
    BarOmega=0.04;
    
    outMS=MixSim(ktrue,v,'BarOmega',BarOmega, 'restrfactor',restrfact);
    % data generation given centroids and cov matrices
    [Y,id]=simdataset(n, outMS.Pi, outMS.Mu, outMS.S);
    % Specify number of solutions
    NumberOfBestSolutions=2;
    % Number of subsets to extract
    nsamp=100;
    % Computation of information criterion using MIXMIX
    outICmixt=tclustICgpcm(Y,'plots',0,'nsamp',nsamp);
    % Plot first 2 best solutions using as Information criterion MIXMIX
    disp('Best 2 solutions using MIXMIX')
    [outMIXMIX]=tclustICsolGPCM(outICmixt,'whichIC','MIXMIX','plots',1,'NumberOfBestSolutions',NumberOfBestSolutions);
    disp(outMIXMIX.MIXMIXbs)
     % Computation of information criterion using CLACLA
    outICcla=tclustICgpcm(Y,'whichIC','CLACLA','plots',0,'nsamp',nsamp);
    [outCLACLA]=tclustICsolGPCM(outICcla,'whichIC','CLACLA','plots',1,'NumberOfBestSolutions',NumberOfBestSolutions);
    disp('Best 2 solutions using CLACLA')
    disp(outCLACLA.CLACLAbs)
%}



%{
   % An example with input options Rand and kk.
    Y=load('geyser2.txt');
    nsamp=100;
    pa=struct;
    pa.cdet=[2 4];
    pa.shw=[8 16 32];
    kk=[2 3 4 6];
    out=tclustICgpcm(Y,'pa',pa,'cleanpool',false,'plots',0,'alpha',0.1,'whichIC','CLACLA','kk',kk,'nsamp',nsamp);
    [outCLACLA]=tclustICsolGPCM(out,'whichIC','CLACLA','plots',1,'NumberOfBestSolutions',3,'Rand',0);
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
% Include or not the plot of spurious solutions.
SpuriousSolutions = true;

% Message about code execution
msg=0;
% Rand =1 implies the adjusted Rand index else FM index is used
Rand = 1;
if nargin>1
    options=struct('whichIC',whichIC,'plots',plots,'SpuriousSolutions',SpuriousSolutions,...
        'NumberOfBestSolutions',NumberOfBestSolutions, ...
        'ThreshRandIndex', ThreshRandIndex,'msg',msg,'Rand',Rand);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        
        
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:tclustICsolGPCM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        
        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:tclustICsolGPCM:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
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
    SpuriousSolutions=options.SpuriousSolutions;
    msg=options.msg;
    Rand=options.Rand;
end


% Extract the values of k (number of groups)
kk=IC.kk;

% Extract the values of restriction factors (number of components)
ccdet=IC.ccdet;
ccshw=IC.ccshw;
% Extract the values of alpha (trimming level)
alpha=IC.alpha;

% Prepare rownames and colsnames for table which will contain
% in the rows the number of groups from
% rownamesARI=strcat(cellstr(repmat('k=',length(kk),1)), cellstr(num2str(kk')));

if strcmp(whichIC,'MIXMIX')
    typeIC=2;
elseif strcmp(whichIC,'MIXCLA')
    typeIC=1;
elseif strcmp(whichIC,'CLACLA')
    typeIC=0;
else
    warning('FSDA:tclustICsolGPCM:WrongIC','Supplied string for whichIC is not supported.')
    error('FSDA:tclustICsolGPCM:WrongIC','Specified information criterion is not supported: possible values are ''MIXMIX'' , ''MIXCLA'',  ''CLACLA''')
end

out=struct;

if typeIC==2
    
    % best solutions in terms of cdet
    PenLik=permute(IC.MIXMIX,[2 1 3]);
    IDX=IC.IDXMIX;
    [cdetMIXMIXbs,cdetMIXMIXbsari]=findBestSolutions(PenLik,IDX,kk,...
        ccdet,ccshw,NumberOfBestSolutions,ThreshRandIndex,msg,Rand);
    
    % best solutions in terms of cshw
    PenLik=permute(IC.MIXMIX,[3 1 2]);
    IDX=permute(IC.IDXMIX,[1 3 2]);
    [cshwMIXMIXbs,~]=findBestSolutions(PenLik,IDX,kk,...
        ccshw,ccdet,NumberOfBestSolutions,ThreshRandIndex,msg,Rand);
    rSolcdet=cell2mat(cdetMIXMIXbs(:,[1 2 6]));
    rSolcshw=cell2mat(cshwMIXMIXbs(:,[1 6 2]));
    % intersection between the solutions of rSolcdet and rSolcshw
    [~,icdet,icshw]=intersect(rSolcdet,rSolcshw,'stable','rows');
    
    MIXMIXbsfin=[cdetMIXMIXbs(icdet,:) cshwMIXMIXbs(icshw,3:4)];
    MIXMIXbsarifin=cdetMIXMIXbsari(icdet,icdet);
    
    out.MIXMIXbs=MIXMIXbsfin;
    out.MIXMIXbsari=MIXMIXbsarifin;
    
    % Store matrix which contains in the columns the details of the
    % classification
    MIXMIXbsIDX=plotBestSolutionsMult(IC.Y,IC.IDXMIX,MIXMIXbsfin,kk,...
        ccdet,ccshw,'MIXMIX',plots,SpuriousSolutions);
    
    out.MIXMIXbsIDX=MIXMIXbsIDX;
    
elseif typeIC==1
    
    % best solutions in terms of cdet
    PenLik=permute(IC.MIXCLA,[2 1 3]);
    IDX=IC.IDXMIX;
    [cdetMIXMIXbs,cdetMIXMIXbsari]=findBestSolutions(PenLik,IDX,kk,...
        ccdet,ccshw,NumberOfBestSolutions,ThreshRandIndex,msg,Rand);
    
    % best solutions in terms of cshw
    PenLik=permute(IC.MIXMIX,[3 1 2]);
    IDX=permute(IC.IDXMIX,[1 3 2]);
    [cshwMIXMIXbs,~]=findBestSolutions(PenLik,IDX,kk,...
        ccshw,ccdet,NumberOfBestSolutions,ThreshRandIndex,msg,Rand);
    rSolcdet=cell2mat(cdetMIXMIXbs(:,[1 2 6]));
    rSolcshw=cell2mat(cshwMIXMIXbs(:,[1 6 2]));
    % intersection between the solutions of rSolcdet and rSolcshw
    [~,icdet,icshw]=intersect(rSolcdet,rSolcshw,'stable','rows');
    
    MIXMIXbsfin=[cdetMIXMIXbs(icdet,:) cshwMIXMIXbs(icshw,3:4)];
    MIXMIXbsarifin=cdetMIXMIXbsari(icdet,icdet);
    
    out.MIXCLAbs=MIXMIXbsfin;
    out.MIXCLAbsari=MIXMIXbsarifin;
    
    % Store matrix which contains in the columns the details of the
    % classification
    MIXMIXbsIDX=plotBestSolutionsMult(IC.Y,IC.IDXMIX,MIXMIXbsfin,kk,...
        ccdet,ccshw,'MIXMIX',plots,SpuriousSolutions);
    
    out.MIXCLAbsIDX=MIXMIXbsIDX;
    
    
else %  typeIC==0
    % best solutions in terms of cdet
    PenLik=permute(IC.CLACLA,[2 1 3]);
    IDX=IC.IDXCLA;
    [cdetCLACLAbs,cdetCLACLAbsari]=findBestSolutions(PenLik,IDX,kk,...
        ccdet,ccshw,NumberOfBestSolutions,ThreshRandIndex,msg,Rand);
    
    % best solutions in terms of cshw
    PenLik=permute(IC.CLACLA,[3 1 2]);
    IDX=permute(IC.IDXCLA,[1 3 2]);
    [cshwCLACLAbs,~]=findBestSolutions(PenLik,IDX,kk,...
        ccshw,ccdet,NumberOfBestSolutions,ThreshRandIndex,msg,Rand);
    rSolcdet=cell2mat(cdetCLACLAbs(:,[1 2 6]));
    rSolcshw=cell2mat(cshwCLACLAbs(:,[1 6 2]));
    % intersection between the solutions of rSolcdet and rSolcshw
    [~,icdet,icshw]=intersect(rSolcdet,rSolcshw,'stable','rows');
    
    CLACLAbsfin=[cdetCLACLAbs(icdet,:) cshwCLACLAbs(icshw,3:4)];
    CLACLAbsarifin=cdetCLACLAbsari(icdet,icdet);
    
    out.CLACLAbs=CLACLAbsfin;
    out.CLACLAbsari=CLACLAbsarifin;
    
    % Store matrix which contains in the columns the details of the
    % classification
    CLACLAbsIDX=plotBestSolutionsMult(IC.Y,IC.IDXCLA,CLACLAbsfin,kk,...
        ccdet,ccshw,'MIXMIX',plots,SpuriousSolutions);
    
    out.CLACLAbsIDX=CLACLAbsIDX;
end

% Store values of ccdet, ccshw and k which have been used
out.ccdet=ccdet;
out.ccshw=ccshw;
out.kk=kk;
out.alpha=alpha;
end



function [ARI]=findARI(IDX,kk, Rand)
[lkk,lcc]=size(IDX);

ARI=ones(lkk,lcc);

for k=1:length(kk) % loop for different values of k (number of groups)
    if kk(k) ~= 1 % just do the computations or ARI or FM indexes when k is different from 1
        for j=2:lcc % loop through the different values of cdet
            
            idjm1=IDX{k,j-1};
            idj=IDX{k,j};
            idjm1(idjm1<0)=0;
            idj(idj<0)=0;
            if Rand==1
                ARI(k,j)=RandIndexFS(idjm1,idj,0);
            else
                ARI(k,j)=FowlkesMallowsIndex(idjm1,idj,0);
            end
        end
    end
end

end

function [Bestsols,ARIbest]  = findBestSolutions(PENloglik,IDX,kk,cc,ccother,NumberOfBestSolutions,ThreshRandIndex,msg,Rand)
% PENloglik = 3D array of size length(cc)-by-length(kk)-by-length(ccother)
% containing penalized log likelihood.
% if cc is ccdet ccother is cshw and viceversa.
% IDX = cell of size length(kk)-by-length(cc)-by-length(ccother) containing the indication
%       where units have been classified for each k (rows) and cc (columns)
%       and ccother (third dimensione)
% kk  = vector containing values of k which have to be considered
% cc  = vector containing values of ccdet or ccshw which have to be
%   considered. If cc contains ccdet than ccother contains the values of
%   ccshw and viceversa.
% ccother  = vector containing values of ccdet or ccshw which have to be
%   considered. If ccother contains ccdet than cc contains the values of
%   ccshw and viceversa.
% NumberOfBestSolutions = scalar, maximum number of solutions to consider.
% ThreshRandIndex = threshold for ARI (FM) index.
% msg = scalar display message about best solution
% Rand = ARI or Fowlkes and Mallows index

lcc=length(cc);
seqcc=1:lcc;
lccother=length(ccother);
seqccother=1:lccother;

seqkk=1:length(kk);
% 1st col = value of best k
% cols 2-4 details of values of cdet (cshw)
% col 5 spurious or true
% cols 6-8 details of values of shw (cdet)
Bestsols=cell(NumberOfBestSolutions,5);
Bestsols{1,5}='true';
endofloop=0;

NumberOfExistingSolutions=NumberOfBestSolutions;
for z=1:NumberOfBestSolutions
    
    % valmin= mimimum for IC each value of k
    [valmin,~]=min(PENloglik,[],1);
    
    % minBICk identifies position of k where there is the optimal solution
    [minlist,~]=min(valmin,[],2);
    % minBICc position of the best value of c where there is the optimal
    
    % minBICk has dim 1x1xlength(ccother)
    % minBICcshw index of ccother for which there is the minimum
    [~,minBICcshw]=min(minlist,[],3);
    
    Bestsols{z,6}=ccother(minBICcshw);
    
    % valmin= mimimum for IC each value of k
    [valmin,indmin]=min(PENloglik(:,:,minBICcshw),[],1);
    
    % minBICk identifies position of k where there is the optimal solution
    [~,minBICk]=min(valmin);
    % minBICc position of the best value of cdet where there is the optimal
    % solution
    minBICc=indmin(minBICk);
    
    if min(valmin)<Inf
        
        if msg==1
            disp(['Best solution number: ' num2str(z)])
            disp(['k=' num2str(kk(minBICk))])
            disp(['cdet=' num2str(cc(minBICc))    ' shw='  num2str(ccother(minBICcshw))])
        end
        % overall minimum is in col minBICk
        % overall minimum is in row minBICc
        PENloglik(indmin(minBICk),minBICk,:)=Inf;
        % Xcshw(indmincshw(minBICkcshw),minBICkcshw,:)=Inf;
        
        cwithbestsol=nan(length(cc),1);
        cwithbestsol(minBICc)=1;
        % Store value of k in column 1 of Bestsols
        Bestsols{z,1}=kk(minBICk);
        % Store value of c in column 2 of Bestsols
        Bestsols{z,2}=cc(minBICc);
        if msg==1
            disp('Find for which adjacent value of c (cdet and cshw) and fixed k best solution extends to')
        end
        
        % PENloglik is 3D
        % 1st dim values of k
        % 2nd dim values of cc
        % 3rd dim values of ccother
        
        % Find overall minimum of matrix Xcmod after excluding the column
        % associated with the value of k where minBICk lies
        % Columns  of Xcmod refer to k
        XcmodWithoutBestk=PENloglik;
        XcmodWithoutBestk(:,minBICk,:)=[];
        
        
        % minICconstr= min value of IC excluding the values involving column
        % minBICk
        minICconstr=min(min(min(XcmodWithoutBestk)));
        cctoadd=zeros(length(cc),1);
        
        % The rows of ARI refer to the values of cc. The columns of ARI refer to k
        % ARI = matrix of size length(kk)-times-length(cc) containing Adjusted R
        %       indexes between two consecutive values of c (cdet or cshw) for fixed k
        
        ARI=findARI(IDX(:,:,minBICcshw),kk, Rand);
        ARI=ARI';
        % Loop over c (cdet or cshw)
        candcabove=seqcc(seqcc>minBICc);
        if ~isempty(candcabove)
            for r=1:length(candcabove)
                posctoadd=candcabove(r);
                if ARI(posctoadd,minBICk)>ThreshRandIndex && PENloglik(posctoadd,minBICk,minBICcshw)<minICconstr
                    PENloglik(posctoadd,minBICk,:)=Inf;
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
                if ARI(posctoadd+1,minBICk)>ThreshRandIndex && PENloglik(posctoadd,minBICk,minBICcshw)<minICconstr
                    PENloglik(posctoadd,minBICk,:)=Inf;
                    cctoadd(posctoadd)=1;
                else
                    break
                end
            end
        end
        
        
        cwithbestsol(cctoadd==1)=1;
        % Store the values of c (cdet or cshw)  associated to the best solution
        Bestsols{z,3}=cc(cwithbestsol==1);
        
        % The interval of values of cdet (cshw) for which the solution is uniformly
        % better has been found, however, now we have to make sure that we
        % do not consider anymore the solutions for the same k which have a
        % value of R index adjacent to those which have already been found
        % greater than a certain threshold
        
        
        if minBICk ==1
            Bestsols{z,4}=cc;
            PENloglik(:,minBICk,:)=Inf;
        else
            intc=seqcc(cwithbestsol==1);
            outc=cc(cwithbestsol~=1);
            cctoadd=zeros(length(cc),1);
            
            if ~isempty(outc)
                % candcbelow = indexes of the set of candidate values for c which are smaller than
                % those associated with best solutions found so far. For
                % example, if candc =(1 2) it means that we check whether for the same value of k
                % solution with cc(2) and cc(1) have an ARI index greater than a certain threshold.
                % If it is the case this means that these solutions do not have to be considered anymore
                candcbelow=seqcc(seqcc<min(intc));
                
                if ~isempty(candcbelow)
                    for r=length(candcbelow):-1:1
                        posctoadd=candcbelow(r);
                        if ARI(posctoadd+1,minBICk)>ThreshRandIndex
                            PENloglik(posctoadd,minBICk,:)=Inf;
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
                            PENloglik(posctoadd,minBICk,:)=Inf;
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
            idxcurrentz=IDX{seqkk(kk==Bestsols{z,1}),seqcc(cc==Bestsols{z,2}), seqccother(minBICcshw)};
            
            for  j=1:z-1
                idxpreviousz=IDX{seqkk(kk==Bestsols{j,1}),seqcc(cc==Bestsols{j,2}),seqccother(ccother==Bestsols{j,6}) };
                if RandIndexFS(idxpreviousz,idxcurrentz,0)>ThreshRandIndex  && strcmp(Bestsols{j,5},'true')==1
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
        idxi=IDX{seqkk(kk==Bestsols{i,1}),seqcc(cc==Bestsols{i,2}), seqccother(ccother==Bestsols{i,6})};
        idxj=IDX{seqkk(kk==Bestsols{j,1}),seqcc(cc==Bestsols{j,2}), seqccother(ccother==Bestsols{j,6})};
        idxi(idxi<0)=0;
        idxj(idxj<0)=0;
        ARIbest(i,j)=RandIndexFS(idxi,idxj,0);
    end
end
end


function IDXout=plotBestSolutionsMult(Y,IDX,Bestsols,kk,ccdet,ccshw,lab,plots,SpuriousSolutions)
seqkk=1:length(kk);
seqccdet=1:length(ccdet);
seqccshw=1:length(ccshw);

nbestsol=size(Bestsols,1);
IDXout=zeros(size(Y,1),nbestsol);
for i=1:nbestsol
    IDXi=IDX{seqkk(kk==Bestsols{i,1}),seqccdet(ccdet==Bestsols{i,2}), seqccshw(ccshw==Bestsols{i,6})};
    IDXout(:,i)=IDXi;
    if (plots==1 &&  SpuriousSolutions ==true) || (plots==1 && strcmp(Bestsols{i,5},'true'))
        figure
        spmplot(Y,IDXi,1,'box');
        detsol=[lab ': sol ' num2str(i)  ': k=' num2str(Bestsols{i,1}) ' c_{det}=' num2str(Bestsols{i,2})    ' c_{shw}=' num2str(Bestsols{i,6})  ];
        bestsol=[' Best in c_{det} ' num2str(min(Bestsols{i,3})) '-' num2str(max(Bestsols{i,3})) ];
        Bestsolsi3=Bestsols{i, 3};
        Bestsolsi4=Bestsols{i,4};
        Bestsolc=[Bestsolsi3(:); Bestsolsi4(:)];
        stabsol=[' Stable in c_{det} ' num2str(min(Bestsolc)) '-' num2str(max(Bestsolc)) ];
        
        bestsol1=[' Best in c_{shw} ' num2str(min(Bestsols{i,7})) '-' num2str(max(Bestsols{i,7})) ];
        Bestsolsi7=Bestsols{i, 7};
        Bestsolsi8=Bestsols{i,8};
        Bestsolcshw=[Bestsolsi7(:); Bestsolsi8(:)];
        stabsolcshw=[' Stable in c_{shw} ' num2str(min(Bestsolcshw)) '-' num2str(max(Bestsolcshw)) ];
        
        title([detsol bestsol stabsol  bestsol1 stabsolcshw ' Sol:' Bestsols{i,5}])
    end
end

if plots==1
    cascade
end
end


%FScategory:CLUS-RobClaMULT