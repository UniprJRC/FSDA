function [Seqsaligned, WrngAlignment]= multialign2ref(refSeqs, Seqs2align, varargin)
%multialign2ref performs multialignment to reference sequences
%
%<a href="matlab: docsearchFS('multialign2ref')">Link to the help function</a>
%
% This function performs multialignment of set of Amino or nucleotide
% sequences not aligned (Seqs2align) to a set of reference sequences
% already aligned (refSeqs). In order to achieve this target function
% multialign of the BioInformatics toolbox is called several times with all
% the default options. Note that if the alignment of a sequence inside
% Seqs2align changes the reference sequences we call again function
% multialign using the 'Name',Value, GapOpen',100000 in order to allow the
% possibility of creating gaps or adding additional characters in the
% reference sequence. If the new call to multialign:
% 1) creates gap(s) in the reference sequences (refSeqs) we delete the
% characters corresponding to the gaps in the aligned sequence
% and store this information inside boolean variable usedGap of output
% table SeqsMultiAligned.
% 2) creates additional characters at the end of the reference sequences
% (refSeqs) we delete the additional characters and store this information
% inside boolean vector usedDeletion of output table SeqsMultiAligned.
%
% As a result of this process the number of sites (characters) in the
% reference sequence is kept fixed and all the aligned sequences will have
% the same number of characters of the reference sequences. We also store
% the information about the sequences which could not be aligned inside
% second output argument WrngAlignment. So far with the hundreds of millions
% of sequences we have aligned this case never took place.
%
% Required input arguments:
%
%     refSeqs : Reference sequences already aligned. Vector of structures.
%               Vector of structures of length a with the fields 
%               refSeqs.Sequence =for the residues and 
%               refSeqs.Header = (or 'refSeqs.Name') for the labels.
%                 Data Types - struct
%
%  Seqs2align : Sequences which have to be aligned. Vector of structures.
%               Vector of structures of length n with the fields 
%               Seqs2align.Sequence = for the residues and
%               Seqs2align.Header = or (or 'Seqs2align.Name') for the labels.
%                 Data Types - struct
%
% Optional input arguments:
%
%  UseParallel  : use or not the parallel computing toolbox.
%                Boolean or positive integer. If UseParallel is true
%                the parallel computing toolbox is used.
%                 Example - 'UseParallelValue','false'
%                 Data Types - boolean or positive integer
%
% ScoringMatrix : scoring method to use for the alignment. Character
%               vector or string. This option specifies the scoring method
%               to use for the alignment. For further details about this
%               option see option ScoringMatrixValue of the multialignment
%               function of the BioInformatics toolbox.
%                 Example - 'ScoringMatrix','BLOSUM30'
%                 Data Types - character or string
%
%
% NumberSeqsEachIter : Number of sequences to add for each iteration.
%                      Positive integer.
%                This option controls the number of sequences to add for
%                each iteration. The default is to add 25 additional
%                sequences for each iteration before calling multialign.
%                 Example - 'NumberSeqsEachIter', 50
%                 Data Types - single or double
%
%  verbose     : Show estimated time to complete the process. Boolean.
%                If verbose is true (default) the estimated time to perform
%                the alignment is shown on the screen.
%                 Example - 'verbose', false
%                 Data Types - logical
%
%
% Output:
%
% Seqsaligned : Sequences which have been aligned. Vector of structures.
%               Vector of structures of the same length n of input
%               Seqs2align with the following 4 fields.
%               Seqsaligned.Sequence =  n aligned sequences
%               Seqsaligned.Header = containing the labels  (this field did
%                   not change from input Seqs2align.Header).
%               Seqsaligned.usedGap =  boolean containing true for the
%                   sequence in which in order to compute the alignment it
%                   was necessary to modify option 'GapOpen', to 100000
%               Seqsaligned.usedDeletion =  boolean containing true when
%                   the sequences to align contained a number of
%                   characters greater than those of the reference
%                   sequences.
%                 Data Types - array of struct
%
% WrngAlignment: sequences which could not be aligned. Vector.
%                Vector containing the numbers of sequences for which
%                usedGap and usedDeletion was not sufficient to produce
%                the alignment. For example if WrngAlignment=[400 800] the
%                it was not possible to align sequences 400 and 800. If
%                WrngAlignment is en empty value all sequences could be
%                aligned.
%                 Data Types - vector with natural numbers
%
% See also multialign,
%
% References:
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('multialign2ref')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of multialign2ref without optional arguments.
    % Load fastafile containing original covid and other sequences
    Seqs2align = fastaread("X01sel.txt");
    
    % Load fasta file containing the 5 covid variants
    % (Alpha, Beta, Delta, Gamma, Omicron)
    variants = fastaread('Variants.txt');
    
    % Original covid sequence and variants make up the reference sequences
    refSequences=[Seqs2align(1); variants];
    % Perform multialignment on the reference sequences
    refSequences=multialign(refSequences);
    
    % Remove initial covid sequence from Seqs2align
    Seqs2align=[Seqs2align(2:101)];
    
    %Call of multialgin2ref with all default arguments
    Seqsaligned=multialign2ref(refSequences,Seqs2align);
%}

%{
    % Example without using the parallel computing toolbox.
    % Load fastafile containing original covid and other sequences
    Seqs2align = fastaread("X01sel.txt");
    
    % Load fasta file containing the 5 covid variants
    % (Alpha, Beta, Delta, Gamma, Omicron)
    variants = fastaread('Variants.txt');
    
    % Original covid sequence and variants make up the reference sequences
    refSequences=[Seqs2align(1); variants];
    % Perform multialignment on the reference sequences
    refSequences=multialign(refSequences);
    
    % Remove initial covid sequence from Seqs2align
    Seqs2align=[Seqs2align(2:50)];

    % Call multialign2ref with option 'UseParallel',false
    Seqsaligned=multialign2ref(refSequences,Seqs2align,'UseParallel',false);
%}

%{
    %% Example of two output arguments.
    % Load fastafile containing original covid and other sequences
    Seqs2align = fastaread("X01sel.txt");
    
    % Load fasta file containing the 5 covid variants
    % (Alpha, Beta, Delta, Gamma, Omicron)
    variants = fastaread('Variants.txt');
    
    % Original covid sequence and variants make up the reference sequences
    refSequences=[Seqs2align(1); variants];
    % Perform multialignment on the reference sequences
    refSequences=multialign(refSequences);
    
    % Remove initial covid sequence from Seqs2align
    Seqs2align=[Seqs2align(2:101)];
    
    %Call of multialign2ref with all default arguments
    [Seqsaligned,WrngAlignment]=multialign2ref(refSequences,Seqs2align);
%}

%% Beginning of code

UseParallel=true;
ScoringMatrix=[];
NumberSeqsEachIter=25;
verbose=true;

options=struct('UseParallel',UseParallel, ...
    'ScoringMatrix',ScoringMatrix,'NumberSeqsEachIter',NumberSeqsEachIter, ...
    'verbose',verbose);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:multialign2ref:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    aux.chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

UseParallel=options.UseParallel;
ScoringMatrix=options.ScoringMatrix;
if isempty(ScoringMatrix)
    defaultScoringMatrix=true;
else
    defaultScoringMatrix=false;
end

NumberSeqsEachIter=options.NumberSeqsEachIter;
verbose=options.verbose;

% Seqsaligned is the struct which will contain the sequences which are
% aligned, it is initialized with Seqs2align
Seqsaligned=Seqs2align;


% nSeqs2align = number of sequences to align
nSeqs2align=length(Seqs2align);
% nrefSeqs = number of reference sequences (already aligned)
nrefSeqs=length(refSeqs);

seq100 = 100*(1:1:ceil(nSeqs2align/100));
seq100boo=false(nSeqs2align,1);
seq100boo(seq100)=true;

WrngAlignment=zeros(1000,1);
usedGap=zeros(1000,1);
usedDeletion=usedGap;
ijWrngAlignment=0;
ijusedGap=0; ijusedDeletion=0;


%% Alignment loop
steplength=NumberSeqsEachIter;
if verbose==true
    tStart=tic;
end
for i=1:(floor(nSeqs2align/steplength)+1)
    if i<=floor(nSeqs2align/steplength)
        sel=((i-1)*steplength+1):(i*steplength);
    else
        sel=(steplength*floor(nSeqs2align/steplength)+1):nSeqs2align;
    end

    if  i==3 && verbose == true 
        disp('Process started')
        disp(datetime('now'))
          b=toc(tStart);
        totest=b*nSeqs2align/(steplength*3);
        if totest>60
            disp(['Total estimated time (minutes)=' num2str(totest/60)])
        elseif  totest>3600
            disp(['Total estimated time (hours)=' num2str(totest/3600)])
        else
            disp(['Total estimated time (seconds)=' num2str(totest)])
        end
    end

    if defaultScoringMatrix ==true
        p1 = multialign([refSeqs; Seqs2align(sel)] ,'UseParallel',UseParallel);
    else
        p1 = multialign([refSeqs; Seqs2align(sel)] ,'UseParallel',UseParallel,'ScoringMatrix',ScoringMatrix);
    end


    if ~isequal(refSeqs,p1(1:nrefSeqs))
        % Do alignment one by one because refSeqs have changed
        for j=1:steplength
            if defaultScoringMatrix ==true
                p1j = multialign([refSeqs; Seqs2align(sel(j))] ,'UseParallel',UseParallel);
            else
                p1j = multialign([refSeqs; Seqs2align(sel(j))] ,'UseParallel',UseParallel,'ScoringMatrix',ScoringMatrix);
            end

            if ~isequal(refSeqs,p1j(1:nrefSeqs))
                % Retry the multialign using option 'GapOpen',100000
                if defaultScoringMatrix ==true
                    p1j = multialign([refSeqs; Seqs2align(sel(j))] ,'UseParallel',true,'GapOpen',100000);
                else
                    p1j = multialign([refSeqs; Seqs2align(sel(j))] ,'UseParallel',true,'GapOpen',100000,'ScoringMatrix',ScoringMatrix);
                end

                if ~isequal(refSeqs,p1j(1:nrefSeqs))
                    % If the original sequences have been modified
                    % we delete the extra columns in the aligned sequences
                    % seqs = cell array of size 1x(nrefSeqs+1)
                    seqs = {p1j(:).Sequence};
                    % seqs1= cell array of size nrefSeqs+1x1
                    seqs1 = strtrim(seqs(:));     % trim sequences
                    % seqs2= character array of size 7-by-unknown_length
                    seqs2=char(seqs1);
                    % Find the columns of seqs2 which have to be deleted (ie those which have -
                    % in the nrefSeqs original sequences)
                    cols2delete=sum(seqs2(1:nrefSeqs,:)==45,1)==6;
                    seqs3=seqs2;
                    seqs3(:,cols2delete)=[];
                    p2j=p1j;
                    for jj=1:(nrefSeqs+1)
                        p2j(jj).Sequence=seqs3(jj,:);
                    end

                    if isequal(refSeqs,p2j(1:nrefSeqs))
                        Seqsaligned(sel(j))=p2j(end);
                        % Seqsaligned.usedDeletion(sel(j))=true;
                        ijusedDeletion=ijusedDeletion+1;
                        usedDeletion(ijusedDeletion)=sel(j);
                    else

                        disp(['Wrong multialignemnt in row=' num2str(sel(j))])
                        ijWrngAlignment=ijWrngAlignment+1;
                        WrngAlignment(ijWrngAlignment)=sel(j);
                    end
                else
                    % In this case the original sequences did not change
                    % thanks to option 'GapOpen',100000
                    ijusedGap=ijusedGap+1;
                    usedGap(ijusedGap)=sel(j);
                    % Seqsaligned.usedGap(sel(j))=true;
                    Seqsaligned(sel(j))=p1j(end);
                end
            else
                % In this case the original reference squences
                % did not change
                Seqsaligned(sel(j))=p1j(end);
            end
        end

    else
        % In this case the original reference squences
        % did not change
        Seqsaligned(sel)=p1((nrefSeqs+1):end);
    end

    if  seq100boo(i) == true
        disp(['i=' int2str(i)]);
    end
    % disp(i)
end
%
usedGap=usedGap(1:ijusedGap);
usedDeletion=usedDeletion(1:ijusedDeletion);
WrngAlignment=WrngAlignment(1:ijWrngAlignment);

% Add fields usedGap and usedDeletion to Seqsaligned
[Seqsaligned.usedGap] = deal(false);
[Seqsaligned.usedDeletion] = deal(false);
% Insert true in correspondence of usedDeletion and usedGap
% of the corresponding fields
[Seqsaligned(usedDeletion).usedDeletion] = deal(true);
[Seqsaligned(usedGap).usedGap] = deal(true);

end

%FScategory:UTIGEN

