function [ROsim,ROsimT]=RatcliffeObesrhelp(STR, varargin)
%RatcliffeObesrhelp computes the Ratcliffe Obershelp algorithm to measure the similarity between two texts
%
%<a href="matlab: docsearchFS('RatcliffeObesrhelp')">Link to the help function</a>
%
% The Ratcliffe Obershelp algorithm (or Gestalt pattern matching) is a
% technique used to measure the similarity between two texts. It was
% developed by Paul Ratcliffe and Peter Obershelp in the 1970s. The
% algorithm works by comparing the words used in both texts and calculating
% a similarity score based on the number of matching words and their
% respective positions within the texts. The algorithm is often used in
% plagiarism detection, document clustering, and information retrieval. It
% is a popular tool for analyzing large amounts of text data and
% identifying patterns and similarities between documents. For additional
% information see the "More about" section in this documents.
%
%
% Required input arguments:
%
%    STR:       Strings or files to compare or list of extensions.
%               String array or cell array of strings of characters or table.
%               Vector or table containting, or the list of strings to
%               compare or the list of files or the list of extensions of
%               the files which have to be used for camparison. The default
%               is to assume that the rows of STR contain the strings to
%               compare. For example, if STR=["first string", "Second
%               string"], then STR(1) is compared with STR(2).
%               On the other hand, if STR=["filename1.m"
%               "fn.txt" "f3.prn"] and optional input argument
%               comparisonType="file", then the contents of files named
%               "filename1.m",  "fn.txt" and "f3.prn" are compared. Finally
%               if STR=["*.m", "*.mlx"] and comparisonType="ext" then all
%               files in the current folder which have extensions .m and
%               .mlx are compared. Note that when ".mlx" extension is used
%               the files are preliminarly transformed into .m format in
%               order to perform the comparison.
%       Data type: string array | character vector | cell array of character vectors | cell array of string arrays | table
%
% Optional input arguments:
%
% comparisonType:  string which specifies what are the elements inside
%                  required input argument STR. Possibile values of
%                  comparisonType are "str" or "file" or "ext".
%                  If comparisonType  = "str" then the strings
%                  inside STR identify texts that have to be compared.
%                  If comparisonType="file" then the elements inside
%                  STR contain the names of the files which have to be
%                  compared.
%                  If comparisonType="ext" then the elements inside
%                  STR contain the extensions of the files which have to be
%                  compared. For example if STR=["*.m" "*.txt"] all files
%                  in the current folder which have extension .m or
%                  extension .txt are used in the comparison.
%                       Example - 'comparisonType',"file"
%                       Data Types - string scalar
%
% remStopWords: remove Stop Words from the strings.
%                If remStopWords is true the following 3 functions of the Text
%                Analytics toolbox are called tokenize, removeStopWords and
%                joinWords, in order to remove stop words. The default
%                value of remStopWords is false and the comparison does not
%                remove stop words.
%                       Example - 'remStopWords',true
%                       Data Types - logical
%
%  Output:
%
%         ROsim:   array with Ratcliffe Obershelp similarity measures.
%                   Array of size n-by-n where n is the numbers of strings
%                   (files) which have been compared. The elements on the
%                   main diagonal are equal to 1.
%         ROsimT:   table with Ratcliffe Obershelp similarity measures.
%                   Table of size n-by-n where n is the numbers of strings
%                   (files) which have been compared. The elements on the
%                   main diagonal are equal to 1. The RowNames
%                   (VariableNames) of ROsimT are equal to the names of the
%                   files which have been selected, if optional input
%                   option comparisonType is equal to "file" or "ext", or
%                   are equal to the RowNames of STR if STR is a table. In
%                   all the other cases the  RowNames
%                   (VariableNames) of ROsimT are equal to "str1", ...,
%                   "strn", where n is the number of elements of STR.
%
%
% More About:
%
% The similarity of two strings $S_1$ and $S_2$ in the Gestalt pattern
% matching (Ratcliff/Obershelp pattern recognition) is determined by twice
% the number of matching characters $K_m$ divided by the total number of
% characters of both strings ($|S_1| + |S_2|$). The matching characters
% $K_m$ are defined as some longest common substring plus recursively the
% number of matching characters in the non-matching regions on both sides
% of the longest common substring In symbols:
%
% \[
%   S_{RO} = \frac{ 2K_m }{ |S_1| + |S_2| }
% \]
%
%  $0 \leq S_{RO} \leq 1$. The value of 1 stands for the complete match of
%  the two strings, whereas the value of 0 means there is no match and not
%  even one common letter.
%
% Note that the similarity score is not a similarity index, because it is
% is not symmetric (not commutative) in the sense that
% \[
%  S_{RO}(S_1, S_2) \ne S_{RO}(S_2, S_1)
% \]
%
%
% See also: tkmeans, tclustIC, tclusteda
%
% References:
%
%   Ratcliff, J. W., Metzener, D. (1988). "Pattern Matching: The Gestalt Approach".
%   Dr. Dobb's Journal, 46.
%   [https://www.drdobbs.com/database/pattern-matching-the-gestalt-approach/184407970?pgno=5]
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('RatcliffeObesrhelp')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %%  Example where input is a set of strings
    STR=["hello"; "hello world"];
    ROsim=RatcliffeObesrhelp(STR);
    disp('Matrix of similarity similarities')
    disp(ROsim)
%}

%{
    %% Another example where input is a set of strings.
    S1="GESTALT PATTERN MATCHING";
    S2="GESTALT PRACTICE";
    STR=[S1;S2];
    % RatcliffeObesrhelp is called with two output arguments
    [ROsim,ROsimT]=RatcliffeObesrhelp(STR);
    disp('Show table of similarity indexes')
    disp(ROsimT)
    disp('Note that the index is not commutative (not symmetric)')
%}

%% Beginning of code

% As default we assume that the elements of STR contain the strings to
% compare
comparisonType="str";
compTypeNumeric=1;
remStopWords=false;

if nargin > 1

    % the 'options' struct in this implementation is provided by the for loop
    options=struct('comparisonType',comparisonType,'remStopWords',remStopWords);

    for i=1:2:(length(varargin)-1)
        options.(varargin{i})=varargin{i+1};
    end


    comparisonType=options.comparisonType;
    remStopWords=options.remStopWords;

    % Make sure that comparisonType contains
    % either "str" or "file" or "ext".
    if comparisonType=="str"
        compTypeNumeric=1;
    elseif comparisonType=="file"
        compTypeNumeric=2;
    elseif comparisonType=="ext"
        compTypeNumeric=3;
    else
        warning('FSDA:RatcliffeObesrhelp:WrongFamily','comparisonType which has been chosen is not supported')
        error('FSDA:RatcliffeObesrhelp:WrongFamily','Supported values are "str" or "file" or "ext"')
    end
end

if istable(STR)
    nameRows=STR.Properties.RowNames;
    STR=STR{:,1};
else
    nameRows='';
end

listOfFilesToDelete=[];

% The user has supplied the list of file extensions inside STR
if compTypeNumeric ==3

    % If "*.mlx" is an extenmsion inside STR it is necessary to convert mlx
    % files to m format
    if ismember("*.mlx",STR)

        % List all .mlx files in the current folder
        files = dir('*.mlx');
        fileNames = {files.name};


        % Convert every .mlx file into .m format
        for i = 1:length(fileNames)
            % name = NameOfFile withtout extension
            [~, name, ~] = fileparts(fileNames{i});

            % Name of destination file
            mFileName = [name, '.m'];

            if exist(mFileName,'file')~=2
                % convert .mlx into formato .m
                matlab.internal.liveeditor.openAndConvert(fileNames{i}, mFileName);

                fprintf('Convert %s in temporary file %s\n', fileNames{i}, mFileName);
                listOfFilesToDelete=[listOfFilesToDelete; string(mFileName)]; %#ok<AGROW>
            else

                %  fprintf('Già Convertito %s in %s\n', fileNames{i}, mFileName);
            end
        end

        % If mlx extension was present make sure that the "*.m" is added to STR
        STR=unique([STR(:); "*.m"]);
        % Remove "*.mlx" because already converted in .m format
        STR(STR=="*.mlx")=[];
    end


    fileNames=[];
    for i=1:length(STR)
        % list all files with required extension in the current folder
        files = dir(char(STR(i)));
        fNames = string({files.name});
        fileNames=[fileNames;fNames(:)]; %#ok<AGROW>
    end
    numComparisons = length(fileNames);
    if isempty(nameRows)
        nameRows=fileNames;
    end

    % elements of STR are FileNames
elseif compTypeNumeric==2
    fileNames=STR;
    numComparisons = length(STR);
    if isempty(nameRows)
        nameRows=fileNames;
    end

else
    numComparisons = length(STR);
    if isempty(nameRows)
        nameRows="str"+(1:numComparisons)';
    end
end



%% Compute similarity matrix (not that the matrix is not symmetric)

% Initialize Ratcliffe Obershelp similarity matrix
ROsim = zeros(numComparisons);

% Compute similarity between each pair of strings (or of files)
for i = 1:numComparisons
    for j = 1:numComparisons

        if compTypeNumeric>1
            % In this case the inputs are the file names
            str1 = fileread(fileNames{i});
            str2 = fileread(fileNames{j});
        else
            % In this case the inputs are the strings inside STR
            str1=STR(i);
            str2=STR(j);
        end

        if remStopWords ==true
            % Step 1: tokenize str1 and Remove stop words
            docSTR1 = tokenizedDocument(str1);
            cleanDoc1 = removeStopWords(docSTR1);
            % Step 2: Convert back to string
            str1 = joinWords(cleanDoc1);
            % Step 1: tokenize str2 and Remove stop words
            docSTR2 = tokenizedDocument(str2);
            cleanDoc2 = removeStopWords(docSTR2);
            % Step 2: Convert back to string
            str2 = joinWords(cleanDoc2);
        end

        % Find index of similarity between str1 and str2
        similarity = ratcliffobershelp(str1, str2);

        % fills the similarity matrix
        ROsim(i, j) = similarity;
    end
end

% Delete the temporary .m files which had been converted from .mlx
if ~isempty(listOfFilesToDelete)
    [~,nameMformat]=fileparts(nameRows);
    [~,nameMLXformat]=fileparts(listOfFilesToDelete);
    [~,ia]=intersect(nameMformat,nameMLXformat);
    nameRows(ia)=nameMLXformat+".mlx";
    for j=1:length(listOfFilesToDelete)
        fileMLXj=listOfFilesToDelete(j);
        disp(['Deleting temporary file '  char(fileMLXj)])
        delete(listOfFilesToDelete(j))
    end
end

% Create the associated table
ROsimT=array2table(ROsim,"RowNames",nameRows,"VariableNames",nameRows);

end

% Begin of inner functions

function similarity = ratcliffobershelp(str1, str2)
% RATCLIFF_OBERSHELP computes the similarity between two strings
% using the Ratcliff/Obershelp algorithm.
% see: https://yassineelkhal.medium.com/the-complete-guide-to-string-similarity-algorithms-1290ad07c6b7
% for implementation and background literature


if nargin < 2
    error('FSDA:ratcliffobershelp:WrgInput','Two input strings are required.')
end

% Convert strings to char arrays if they are not already
if ~ischar(str1)
    str1 = char(str1);
end
if ~ischar(str2)
    str2 = char(str2);
end

% Km = number of matching characters
% The matching characters are defined as some longest common substring
% plus recursively the number of matching characters in the non-matching
% regions on both sides of the longest common substring
Km = NumberOfMatchingCharacters(str1, str2);

% Calculate the Ratcliff Obershelp measure of similarity
similarity = 2*Km/(length(str1)+length(str2));

end

function score = NumberOfMatchingCharacters(s1, s2)
% case 1: empty string(s)
if isempty(s1) || isempty(s2)
    score = 0;
    return
end

% case 2: identical strings
if strcmp(s1, s2)
    score = length(s1);
    return
end

% Find the longest common substring (lcs)
%  The algorithm focuses on finding and aligning the longest common
%  substrings, which helps in recognizing significant similarities even in
%  the presence of small differences.


[lcs, start1, start2] = longest_common_substring(s1, s2);

% Stop recursion when no more common substrings are found
if isempty(lcs)
    score = 0;
    return
end

score = length(lcs);
% Recursively compute the score
% By recursively applying the algorithm to the segments of the strings, it
% accounts for multiple overlapping matches, making it robust against minor
% variations.
score = score + NumberOfMatchingCharacters(s1(1:start1-1), s2(1:start2-1)) ...
    + NumberOfMatchingCharacters(s1(start1+length(lcs):end), s2(start2+length(lcs):end));
end

function [lcs, start1, start2] = longest_common_substring(s1, s2)
% Initialize variables
lcs = '';
start1 = 0;
start2 = 0;
len1 = length(s1);
len2 = length(s2);
max_len = 0;

% Create an array to store lengths of longest common suffixes
lcsuff = zeros(len1+1, len2+1);

% Build the table in bottom-up fashion
for i = 1:len1
    for j = 1:len2
        if s1(i) == s2(j)
            lcsuff(i+1, j+1) = lcsuff(i, j) + 1;
            if lcsuff(i+1, j+1) > max_len
                max_len = lcsuff(i+1, j+1);
                lcs = s1(i-max_len+1:i);
                start1 = i - max_len + 1;
                start2 = j - max_len + 1;
            end
        end
    end
end
end
