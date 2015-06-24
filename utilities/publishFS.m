function publishFS(file,varargin)
%Enables to create automatic HELP FILES from structured .m function files
%
% Required input arguments:
%
%    file:         MATLAB File. String. Full or partial
%                  path of the MATLAB file for which Structured Matlab
%                  HTML help has to be created
%                  Example-'myfile.m'
% The .m file (which must be located on the MATLAB path or on the currect
% folder) must satisfy the following characteristics to be correctly
% processed.
%
% 1) The row below the row which starts with function .... must contain the
% description of the purpose of the .m file. Remark: generally the row
% which starts with function .... is the first row of an .m file
% 2) String 'Required input arguments:' must be present. The lines below
% this string must contain the description of the compulsory input
% arguments. Each argument must have  the name (a series of spaces from 0
% to 10) symbol ':' and then the description. The format of the description
% is as follows:
% The first sentence after symbol ':' is the title of the input argument
% and in the HTML file it will appear in bold face in the same line of the
% input argument. (this is the short description of the optional input
% argument). The second sentence after symbol ':' describes the objects
% (for example, scalar, vector, 3D array) and in the HTML file will appear
% in the second row. These first two rows will always be visible. What
% starts with the third sentence after symbol ':' is the detailed
% description of that particular input argument and in the HTML file it
% will be visible just if the user clicks on any point in the first two
% lines or the user clicks on the option expand all. The last line may
% start with the words "Data Types:" and contains the specification of a
% particular input argument (e.g. Data Types: single | double). For
% example, suppose that the .m routine which has to be processed has two
% required input arguments which are respectively called y and X, then the
% accepted format is as follows.
%
%               Required input arguments:
%
%                y:         Response variable. Vector. Response variable,
%                           specified as a vector of length n, where n is
%                           the number of observations. Each entry in y is
%                           the response for the corresponding row of X.
%                           Data Types: single | double
%              X :          Predictor variables. Matrix of explanatory
%                           variables (also called 'regressors') of
%                           dimension n x (p-1) where p denotes the number
%                           of explanatory variables including the
%                           intercept. Rows of X represent observations,
%                           and columns represent variables. By default,
%                           there is a constant term in the model, unless
%                           you explicitly remove it using input option
%                           intercept, so do not include a column of 1s in
%                           X.
%                           Data Types: single | double
%
% IMPORTANT NOTICE: if an input argument is a structure (publishFS
% automatically checks if the input argument contains the word "structure"
% then the fields of the structure will be automatically included into a
% HTML table). In this case the fields of the structure are identified as
% the lines which contain "a series of spaces" "a_word" "a series
% of spaces followed by symbol '='". For example suppose the an input
% option is called bayes and object bayes is a structure with field names
% beta0, R, tau0 and n0, the accepted format is as follows.
%
%    bayes      : prior information. Structure.
%                       It contains the following fields
%               beta0=  p-times-1 vector containing prior mean of \beta
%               R    =  p-times-p positive definite matrix which can be
%                       interpreted as X0'X0 where X0 is a n0 x p matrix
%                       coming from previous experiments (assuming that the
%                       intercept is included in the model
%               tau0 = scalar. Prior estimate of
%                       \[ \tau=1/ \sigma^2 =a_0/b_0 \]
%               n0   = scalar. Sometimes it helps to think of the prior
%                      information as coming from n0 previous experiments.
%
% 2) String 'Optional input arguments:' must be present even if there are
% no optional arguments. publishFS, in order to understand what are the names
% of the optional input arguments scans the rows below the string 'Optional
% input arguments:' and identifies the lines which contain the optional
% arguments as those which contain "a series of spaces" "a_word" "a series
% of spaces followed by symbol ':'". The first sentence after symbol ':' is
% the title of the optional input argument and in the HTML file it will
% appear in the same row of the name of the optional input argument (this
% is the short description of the optional input argument). The second
% sentence after symbol ':' describes the objects (for example, scalar,
% vector, 3D array) and in the HTML file will appear
% in the second row. These first two rows will always be visible. What
% starts with the third sentence after symbol ':' is the detailed
% description of that particular optional input argument and in the HTML
% file it will be visible just if the user clicks on any point in the first
% two lines or the user clicks on the option expand all.
% The last two lines of each optional input argument MUST start with the
% words 'Example -' and 'Data Types -' followed by a string without spaces
% which specifies a possible use of the option and the type of data. For
% example, suppose that the first two optional arguments are called
% respecively 'intercept' and 'h', then the
% accepted format is as follows
%
%               Optional input arguments:
%
%               intercept :  Indicator for constant term. Scalar.
%                           If 1, a model with constant term will be fitted
%                           (default),
%                           if 0, no constant term will be included.
%                           Example - 'intercept',1
%                           Data Types - double
%                       h : The number of observations that have determined the least
%                             trimmed squares estimator. Scalar.
%                             Example - 'h',round(n*0,75)
%                             Data Types - double
%
%
%  IMPORTANT NOTICE: given that options are identified as those which have
%  symbol "%" followed by "a series of spaces" then "a word" then "a series
%  of spaces" then symbol ":", each line inside the description does not
%  have to start as follows "%   ANYWORD   :" because the parser will
%  wrongly identify "ANYWORD" as an additional optional input argument. The
%  only once exception to this rule is the word "%  REMARK :". However, if
%  there is a remark, it must be put at the very end of the description of
%  the optional input argument. At the very end means after the rows
%   Example and Data Types
%
% 4) String 'Output:' must be present. The lines after string 'Output:'
% must contain the list of the output arguments using the same rules
% described above for the optional arguments. In this case, however, the
% identification of the ouptut arguments is easier because they are
% extracted directly from the first line of the file (e.g. if the first
% line is as follows
% function [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,bsb,varargin) then the 5
% output arguments are immediately known to the parser).
% In the case of output argument publishFS checks if the first 50
% characters contain the words "which contains" or "containing" e.g.:
%
%              mdr:         n -init x 2
%                           matrix containing the
%                           monitoring of minimum deletion residual.
%                           1st col = fwd search
%                           ........
%                Un:        (n-init) x 11 Matrix
%                           which contains the unit(s) included in the
%                           subset at each step of the search.
%                           ...........
%
% In this case what is after the strings "which contains" or "containing"
% will appear in bold face as the title of the output argument. What is
% before the strings "which contains" or "containing" will appear in the
% second row.
%
% For example, the above lines will be processed as follows:
%
%      mdr   -  Monitoring of minimum deletion residual
%      n -init x 2 matrix
%
%      Un    - unit(s) included in the subset at each step of the search.
%      (n-init) x 11 Matrix
%
% If in the HTML file the user clicks on them the expdanded description
% (that is what starts after the second full stop will appear).
%
% Alternatively, if the first 50 characters of each output argument do not
% contain the strings "which contains" or "containing" the following
% convention is used. The first sentence after symbol ":" is assumed
% to be the title of the output argument and in the HTML file it will
% appear in bold face in the same line of the name of output
% argument. The second sentence (words between first and second full stop)
% will appear in the second row. The third sentence is the full description
% of the output argument. For example, suppose that the output of a
% procedure contains the objects mdr and Un, the accepted format
% is as follows.
%
%              %  Output:
%              mdr:         Minimum deletion residual. Matrix.  n -init x 2
%                           matrix which contains the
%                           monitoring of minimum deletion residual.
%                           1st col = ...
%              Un:          Units included. Matrix. (n-init) x 11 Matrix
%                           which contains the unit(s) included in the
%                           subset at each step of the search.
%                           REMARK: in every step ....
%
% The above lines will be processed as follows:
%
%      mdr   -  Minimum deletion residual
%      Matrix
%
%      Un    - Units included.
%      Matrix
%
% If in the HTML file the user clicks on them the expdanded description
% (that is what starts after the second full stop will appear).
%
% IMPORTANT NOTICE: Similarly to what happend for each input argument, if
% an output argument is a structure, publishFS automatically checks if it
% contains the words "structure" and "field". In this case, the fields of
% the structure will be automatically included into a HTML table. The
% fields of the structure are identified as the lines which contain "a
% series of spaces" "name_of_output_structure.a_word" "a series of spaces
% followed by symbol '='". For example suppose that the output of a
% procedure is an object called out which is a structure with two fields
% out.rew and out.beta, an accepted format is as follows
%
%                          %  Output:
%
%                          out :     A structure containing the following fields
%
%                                    out.rew  = Scalar if out.rew=1 all
%                                               subsequent output refers to
%                                               reweighted else no
%                                               reweighting is done.
%                                    out.beta = Vector of beta LTS (LMS)
%                                               coefficient estimates,
%                                               including the intercept
%                                               when options.intercept=1.
%                                               out.beta=[intercept
%                                               slopes].
%
%
% PLEASE REMEMBER THAT THE FIELDS of an output instance HAVE TO CONTAIN THE
% "=" SIGN AND NOT THE ":" SIGN
%
% REMARK: If there is the string REMARK after the description of the last
% field of the structure, all the words after REMARK are put outside and
% below the HTML table
%
% If the description of a particular output has the string "which contains"
% or "containing",  as follows
%
%              mdr:          n -init x 2 matrix which contains the
%                           monitoring of minimum deletion residual at each
%                           step of the forward search.
%                           1st col = fwd search index (from init to n-1).
%                           2nd col = minimum deletion residual.
%
%publishFS will try to put what comes before the string "which
%contains" or "containing" inside the subtitle (second row) of the each
%ouptut argument in the HTML file. For example, the example above in the
%HTML file will be processed as follows:
%                mdr —Monitoring of minimum deletion residual at each step of the forward search.
%                n -init -by- 2 matrix
% If, in the HTML file, the user clicks on the first line,
%               "mdr —Monitoring..."
%the expanded description will automatically appear
%
% 5) A line which starts with string 'See also:' must be present. Linked m
% files must be separated by symbol ",". For example, suppose that files
% FSRBmdr.m and FSR.m have connections with the current file, then an
% accepted format is
%
%                   See also: FSRBmdr, FSR.m
%
%
% 6) A line which starts with string 'References:' must be present.
% The year of each reference must be enclosed in round parenthesis.
% PublishFS decides that a new reference starts if a new line contains
% symbol "(" + "a sequence of 4 or 5 characthers identifying the year
% because the reference can be for example 2003 or 2003a" + symbol ")"
% For example, an acceptable format for the two references below is
%
%
%                 Chaloner and Brant (1988). A Bayesian Approach to Outlier
%                 Detection and Residual Analysis, Biometrika, Vol 75 pp.
%                 651-659.
%                 Riani M., Corbellini A., Atkinson A.C. (2015), Very
%                 Robust Bayesian Regression for Fraud Detection, submitted
%
% 7) All the examples associated with the file which has to be processed
% must be enclosed inside Percent-braces (comments blocks, i.e. smbols %{
% and %} ). The first sentence identifies the title of the comment which
% will appear in the HTML file.
% IMPORTANT NOTICE: if the comment has to be executed, the first line
% associated with the title must start with two "%%" symbols instead of just
% one "%" symbol. The examples in the first positions will appear in
% the HTML file under the caption "Examples" while the latest will appear
% under the caption "Related Examples". More precisely, if the output of a
% procedure contains k outputs and some optional input arguments the first
% k+1 comment blocks will appear in the HTML file under "Examples".
% First comment block is associated with the call of the procedure with
% just one output and all default input arguments
% Second comment block is associated with the call of the procedure with
% just one output and with some optional input arguments
% Third comment block is associated with the call of the procedure with
% two output arguments
% ...
% k+1 comment block is associated with the call of the procedure with
% k output arguments
% k+2 comment block is the first which in the HTML file will appear under
% the heading "Related Examples
% For example, suppose that the first example of procedure FSRmdr has to be
% executed and its output must be included into the HTML file, then the accepted
% format is as follows ("please notice the two symbols %% in the first row")
%
%
%                 %{
%                     %% FSRmdr with all default options.
%                     % Compute minimum deletion residual.
%                     % Monitor minimum deletion residual in each step of the forward search.
%                     % Common part to all examples: load fishery dataset.
%                      load('fishery');
%                      y=fishery.data(:,2);
%                      X=fishery.data(:,1);
%                      % Find starting subset
%                      [out]=LXS(y,X,'nsamp',10000);
%                      [mdr] = FSRmdr(y,X,out.bs);
%                      plot(mdr(:,1),mdr(:,2))
%                      title('Monitoring of minimum deletion residual')
%                 %}
%
%                 %{
%                     % FSRmdr with optional arguments.
%                     % Choose step to start monitoring.
%                     % Compute minimum deletion residual and start monitoring it from step
%                     % 60.
%                     [mdr] = FSRmdr(y,X,out.bs,'init',60);
%                 %}
%
% 8) If a procedure contains varargout then a section string
%               Optional Output:
% must be present. For example suppose there is a function called mcd which
% has the following sintax:
%
%                 function [RAW,REW,varargout] = mcd(Y,varargin)
%
% then at the end of the output argument the format must be as follows
%
%                       Optional Output:
%
%                                 C     : matrix of the indices of the
%                                         subsamples extracted for
%                                         computing the estimate
%
% 9) If the .m file contains string  "More About:" a particular section
% called "More about" in the HTML file will be created (just before See
% Also).
% 10) If the .m file contains string 'Acknowledgements:' then a particular
% section named "Acknowledgements" will be created just above the
% references.
%
% GENERAL REMARKS:
%
% ------------------------------------------------------------------------
% REMARK1: if symbol % is wanted and it is not a simple comment delimiter, it
% must be replaced by words "per cent". For example, string "50% envelope"
% must become "50 per cent" envelope.
% ---------------------------------------------------
% REMARK2:
% If there is just one output argument it can be without square brackets.
% Among the input elements of a procedure the number of spaces between
% them is not important. For example
% "y,X,varargin" or "y, X   ,  varargin"   are both fine
% --------------------------------------------------
% REMARK 3: publishFS uses javascript matlab-highlighter.min.js in order to
% automatically color the examples in the HTML file
% --------------------------------------------------
% REMARK 4: publishFS uses MathJax javascript in order to interpret the
% equations inside the .m file written in Latex style. For in line
% equations both symbols $ $ and \( \) are accepted.
% For single line equations symbols \[ \] must be used
% ---------------------------------------------------
% REMARK 5: if there are not enough examples in the .m file the procedure
% still runs but a warning will be automatically produced
% -----------------------------------------------------
% REMARK 6: the output file to be correctly viewed must be located in a
% folder which contains a subfolder named includesFS containing the files
% present inside (home FSDA)\helpfiles\FSDA\includesFS. Therefore if the
% the directory which contains the output file is not located inside (home
% FSDA)\helpfiles\FSDA subfolder includesFS must be copied into the current
% folder
% -----------------------------------------------------
% REMARK 7: strings are HTML formatted as follows. Every time there is
% symbol ". one_or_more_space symbol_of carriage_return" or
% ". one_or_more_space symbol_of carriage_return" the parser adds HTML
% string '</p>/<p>' after just symbol "."  or symbol ":".
% This is done using subfunction named formatHTML at the end of this file.
% subfunction formatHTMLwithMATHJAX is even more general because it applies
% function formatHTML just to the parts of the input string which are not
% enclosed inside symbols '\[ \]'.
%
%
% Optional input arguments:
%
%   Display : Level of display. String.
%             'off' or 'none' displays no output.
%             'iter' displays a series of messages on the screen about
%             the execution of the different section of the input .m file
%             'iter-detailed' displays a series of messages on the screen not only about
%             the execution of the different section of the input .m file,
%             but also about cells containing information about the required input arguments,
%             optional input arguments, and output arguments
%             Example - 'Display','none'
%             Data Types - string
% outputDir : Output folder. String.
%             Output folder to which the HTML document is saved, specified
%             as the comma-separated pair consisting of 'outputDir' and the
%             full path. You must specify the full path as a string, for
%             example 'C:\myPublishedOutput'.
%             The default value, '', specifies the (FSDA root)\helpfiles\FSDA
%             path.
%             Remark - outputDir must be a valid path.
%             Example - 'outputDir','C:\myPublishedOutput'
%             Data Types - string
% imagesDir : Output folder of png images. String.
%             Output folder to which the images attached to the HTML
%             document are saved, specified as the comma-separated pair
%             consisting of 'outputDir' and the full path. You must specify
%             the full path as a string, for example
%             'C:\myPublishedOutput'.
%             The default value, '', specifies the
%             "(FSDA root)\helpfiles\FSDA\images" path.
%             Remark: if imageDir is not specified but outputDir is
%             specified images will be save into the same folder of the
%             HTML output file
%             Remark - imagesDir must be a valid path.
%             Example - 'outputDir','C:\myPublishedOutput'
%             Data Types - string
% evalCode :  Option to run code. Logical. Option to evaluate code of the
%             examples in the input .m files enclosed in tags %{ %} whise
%             first line starts with symbols %% and include the MATLAB
%             output in the HTML file. The iamges will be save in subfolder
%             iamges_help of the outputDir. The default value of evalCode
%             is true.
%             Example - 'evalCode','false'
%             Data Types - Boolean
%
%
% See also: publish
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('publishFS')">Link to the help function</a>
% Last modified 06-Feb-2015

% Examples:

%{
  % Create file FSRmdr.html starting from file FSRmdr.
  % Display detailed information about the Input, Output and Optional arguments
  publishFS('FSRmdr','evalCode',false,'Display','iter-detailed')
%}

%{
  % Create HTML file with embedded images in folder C:\tmp
    publishFS('FSRmdr','evalCode',true,'outputDir','C:\tmp')
%}

%% Beginning of code

if ~ischar(file)
    error('FSDA:publishFS:WrongInputOpt','input must be a string containing the name of the file.');
end

evalCode=true;
Display='none';

% Write output file in subfolder \(FSDAroot)\helpfiles\FSDA
FileWithFullPath=which('docsearchFS.m');
[pathstr]=fileparts(FileWithFullPath);
outputDir=[pathstr '\helpfiles\FSDA'];
imagesDir=[pathstr '\helpfiles\FSDA\images'];

if nargin>1
    options=struct('evalCode',evalCode,'Display',Display,'outputDir',outputDir);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:regressB:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
        
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin);
            options.(varargin{i})=varargin{i+1};
        end
        
    end
    
    evalCode=options.evalCode;
    outputDir=options.outputDir;
    Display=options.Display;
    checkimageDir = strcmp(UserOptions,'imageDir');
    checkoutputDir = strcmp(UserOptions,'outputDir');
    
    if sum(checkimageDir)==0 && sum(checkoutputDir)==1
        imagesDir=outputDir;
    elseif sum(checkimageDir)==1 && sum(checkoutputDir)==1
        imagesDir=options.imageDir;
    else
    end
end

%% Open input .m file put it in astring a do preliminary operations
FileWithFullPath=which(file);
[pathstrcf,name,ext]=fileparts(FileWithFullPath);

if ~strcmp('.m',ext)
    error('FSDA:publishFS:WrongFileExt','Input file must have m extension')
end

if isempty(pathstrcf)
    error('FSDA:publishFS:WrongFile','SourceNotFound');
end

filename=FileWithFullPath;
% f = fopen(filename);

fileID = fopen(char(filename), 'r+');


% Insert the file into fstring
fstring=fscanf(fileID,'%c');
% Replace < and > symbols with HTML code
%fstring=regexprep(fstring,'[^%]<','&lt;');
fstring=regexprep(fstring,'<','&lt;');
fstring=regexprep(fstring,'>','&gt;');


%-----------------
%% CREATE Name-Value Pair Arguments SECTION OF HTML FILE
inselOpt=regexp(fstring,'Optional input arguments:');
if isempty(inselOpt)
    disp('Please check HTML input file')
    error('FSDA:missInps','HTML file does not contain ''Optional input arguments:'' string')
end

% substring to search start from Optional input arguments:
fstringselOpt=fstring(inselOpt(1):end);

endpoint=regexp(fstringselOpt,'Output:');
if isempty(endpoint)
    disp('Please check HTML input file')
    error('FSDA:missOuts','HTML file does not contain ''Output:'' string')
end
fstringselOpt=fstringselOpt(1:endpoint-2);

% Find any string which
% begins with % sign then
% contains a series of white space which can go from 0 to 15 then
% contains any single word
% a series of white spaces which can go from 0 to 10 then
% character : then
% a seris of white spaces
% The inipoint of the following two regular xpessions is the same but
% however we want to exclude the lines where symbol : is followed by a
% series of white spaces and then by a carriage return because these are
% input arguments but simply are the beginning of a list
[~,fint]=regexp(fstringselOpt,'%\s{0,15}\w*\s{0,10}:');
[ini,~]=regexp(fstringselOpt,'%\s{0,15}\w*\s{0,10}:\s{0,10}\w');
fin=fint(1:length(ini));

% listOptArgs = list which contains all optional arguments
% The first column will contain the names (just one word)
% The second column will contain the title of the option (the first
% sentence which finishes with a full stop sign)
% The third column will contain the type (the second sentence which
% finishes with a comma or full stop sign)
% The fourth column will contain the long description. What starts with the
% third sentence
% The fifth column will contain the example what starts just after
% string  Example -
% The sixth column will contain the example what starts just after
% string  Data Types -
listOptArgs=cell(length(ini),6);

ij=1;
for i=1:length(ini)
    % fin(i)-1 because character ':' does not have to be extracted
    opti=fstringselOpt(ini(i):fin(i)-1);
    % Remove from string descri all '% signs
    posPercentageSigns=regexp(opti,'%');
    opti(posPercentageSigns)=[];
    % Remove from string opti leading and trailing white spaces
    opti=strtrim(opti);
    % Check if optional argument is the string rEmArK (written in a case
    % insensitive way)
    
    CheckIfRemark=regexp(opti,'remark','match','ignorecase');
    if ~isempty(CheckIfRemark)
        if i<length(ini)
            descradd=fstringselOpt(ini(i):ini(i+1));
        else
            descradd=fstringselOpt(ini(i):end);
        end
        
        % Remove from string descradd all '% signs
        posPercentageSigns=regexp(descradd,'%');
        descradd(posPercentageSigns)=[];
        descradd=strtrim(descradd);
        listOptArgs{ij-1,4}=[listOptArgs{ij-1,4} descradd];
        % listOptArgs{i-1,2}=[listOptArgs{i-1,2}
    else
        % Store name in the first column of listOptArgs
        listOptArgs{ij,1}=opti;
        % Store short description in the 3nd col of listOptArgs
        if i<length(ini)
            descrtosplit=fstringselOpt(fin(i)+1:ini(i+1)-1);
        else
            descrtosplit=fstringselOpt(fin(i)+1:end);
        end
        
        % Remove from string descrtosplit all '% signs
        posPercentageSigns=regexp(descrtosplit,'%');
        descrtosplit(posPercentageSigns)=[];
        
        [inifullstops]=regexp(descrtosplit,'\.');
        if isempty(inifullstops)
            error('FSDA:publishFS:WrongInp',['Sentence''' descrtosplit '''must contain at least two full stops'])
            % error('Wrong input')
        end
        descrtitle=strtrim(descrtosplit(1:inifullstops(1)-1));
        listOptArgs{ij,2}=descrtitle;
        
        try
            descrtype=strtrim(descrtosplit(inifullstops(1)+1:inifullstops(2)-1));
        catch
            error('FSDA:publishFS:WrongInp',['Option: ' listOptArgs{ij,1} '\nSentence\n''' strtrim(descrtosplit) '''\nmust contain at least two full stops'])
        end
        
        
        listOptArgs{ij,3}=descrtype;
        
        try
            descrlong=strtrim(descrtosplit(inifullstops(2)+1:end));
        catch
            error('FSDA:publishFS:WrongInp',['Sentence''' descrtosplit '''must contain at least two full stops'])
        end
        
        % Check if descr long contains
        % Example - and Data types -
        
        CheckExample=regexp(descrlong,'Example -','once');
        if ~isempty(CheckExample)
            Datatypes=regexp(descrlong,'Data Types -','once');
            listOptArgs{ij,4}=descrlong(1:CheckExample-1);
            
            % The first word of example code must be embedded around tags <code> </code>
            examplecode=descrlong(CheckExample+10:Datatypes-1);
            posspace=regexp(examplecode,'      ');
            examplecode=['<code>' examplecode(1:posspace-1) '</code>' examplecode(posspace:end)];
            listOptArgs{ij,5}=strtrim(examplecode);
            listOptArgs{ij,6}=descrlong(Datatypes+13:end);
        else
            listOptArgs{ij,4}=descrlong;
        end
        
        ij=ij+1;
    end
    
    %     if strcmp(opti,'
    %     listOptArgs{i}=opti;
end
listOptArgs=listOptArgs(1:ij-1,:);

if strcmp(Display,'iter-detailed')
    disp('Detailed information about Optional arguments')
    disp(listOptArgs)
end



%------------------

% Define iniTable and cloTable which are respectively the beginning and the
% end of the table which will contains the fields of the structure argument
% if present in the input or output parameters of the procedure
iniTable=[sprintf(['<table border="2" cellpadding="4" cellspacing="0" class="body">\r'...
    '<colgroup>\r']) ...
    '<col width="50%"><col width="50%">'...
    sprintf(['\r</colgroup>\r'...
    '<thead>\r'...
    '<tr valign="top">\r'...
    '<th valign="top">Value</th>\r'...
    '<th valign="top">Description</th>\r'...
    '</tr>\r'...
    '</thead>'])];
cloTable='</table>';


%% Add title
beforetitl=['<!DOCTYPE HTML> \r'  ...
    '<html itemscope="" xmlns="http://www.w3.org/1999/xhtml">\r' ...
    '<head>\r' ...
    '<title>\r'];
aftertitle='</title>';
titl=sprintf([beforetitl    name  aftertitle]);

%% Add purpose of the file (extract what is in the second row of .m file)
beforemetacontent=['<meta content="refpage" name="chunktype">\r' ...
    '<meta content="function:' name '" itemprop="refentity" name="refentity">\r'...
    '<meta content="text/javascript" http-equiv="Content-Script-Type">\r'...
    '<meta content="fcn" itemprop="pagetype" name="toctype">\r'...
    '<meta content="ref/function" itemprop="infotype" name="infotype" />\r'...
    '<meta content="'];
[startIndex] = regexp(fstring,'%');
% startIndex(2)-3 because there is also the carriage return
purpose=fstring(startIndex(1)+1:startIndex(2)-3);
aftermetacontent=['." itemprop="description" name="description" />\r'...
    '<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />\r'...
    '<meta http-equiv="X-UA-Compatible" content="IE=EmulateIE7" />\r'...
    '<script type="text/x-mathjax-config">\r'...
    'MathJax.Hub.Config({\r'...
    'extensions: ["tex2jax.js"],\r'...
    'jax: ["input/TeX","output/HTML-CSS"],\r'...
    'menuSettings: {zoom: "Double-Click", zscale: "300"},\r'...
    'tex2jax: {inlineMath: [["$","$"],["\\\\(","\\\\)"]]},\r'...
    'MathMenu: {showRenderer: false},\r'...
    '"HTML-CSS": {\r'...
    'availableFonts: ["TeX"],\r'...
    'preferredFont: "TeX",\r'...
    'imageFont: null\r'...
    '}\r'...
    '});\r'...
    '</script>\r'...
    '<script type="text/javascript" src="includesFS/Mathjax/MathJax.js"></script>\r'...
    '<script src="includesFS/matlab-highlighter.min.js"></script>\r'...
    '<link href="includesFS/matlab-highlighter.css" rel="stylesheet" type="text/css">\r'...
    '<script src="includesFS/jquery-latest.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/l10n.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/docscripts.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.toc.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.imagescaling.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.imageanimation.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/bottom.js" type="text/javascript"></script>\r'...
    '<link href="includesFS/reset.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/960.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/site5.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center_installed.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includes/product/css/doc_center_print.css" media="print" rel="stylesheet" type="text/css">\r'...
    '</head>\r'...
    '<body  onload="highlightMATLABCode();">\r'];

metacontent=sprintf([beforemetacontent purpose aftermetacontent]);
% Insert navigation bar on top of the page
% necessary to insert sprintf later because there is symbol % in 100%
insnav=['<table border="0" cellpadding="0" cellspacing="0" class="nav" width="100%">'...
    sprintf(['\r<tr valign="top">\r'...
    '<td align="left" width="20"><a href=" ">\r'...
    '<img align="bottom" alt="" border="0" src="images_help/b_prev.gif"></a>\r'...
    '</td>\r'...
    '<td align="left">left</td>\r'...
    '<td>&nbsp;</td>\r'...
    '<td align="right"> right</td>\r'...
    '<td align="right" width="20"><a href=" ">\r'...
    '<img align="bottom" alt="score" border="0" src="images_help/b_next.gif"></a></td>\r'...
    '</tr>\r'...
    '</table>'])];

%% CONTAINER + SINTAX SECTION
% Di seguito qullo che c'è in doc_center.css
%{
/* TOC styles */
.site_container { padding-right:30px; }
.site_container.site_toc_closed { margin-left:175px; }
.site_container.site_toc_opened { margin-left:300px; }
%}

%{
doc_center_installed.css nella parte che segue controlla il margine
verticale
/* TOC */
.site_container.site_toc_closed { margin-left:600px; } QUIQUI
.toc_pane { padding-top:0px; }

.toc_container_wrapper { position:fixed; top:15px; }

/* Fixed Search Box and Breadcrumbs */
#search_crumb_container { padding:15px 0px 10px; background:#fff; float:left; position:fixed; top:0px; z-index:1001; }
.content_container {
	padding-top: 10px;      QUIQUI
}
%}

% inisitecont=sprintf(['<div class="site_container site_toc_opened">\r'...
inisitecont=sprintf(['<div class="site_container site_toc_closed">\r'...
    '<div class="page_container">\r'...
    '<div class="content_frame">\r'...
    '<div id="content_container" class="content_container">\r'...
    '<section id="doc_center_content">\r'...
    '<div class="function_ref">\r'...
    '<div itemprop="content">\r'...
    '<h1 itemprop="title">'  name '</h1>\r'...
    '<div class="doc_topic_desc">\r'...
    '<div class="purpose_container">\r'...
    '<p itemprop="purpose">\r'...
    purpose '\r'...
    '</p>\r'...
    '<div class="switch">\r'...
    '<a id="expandAllPage" href="javascript:void(0);">\r'...
    'expand all in page</a></div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '<div class="ref_sect">\r']);

[gendescini]=regexp(fstring,'</a>','once');
[gendescfin]=regexp(fstring,'Required input arguments','once');
gendesc=fstring(gendescini+4:gendescfin-1);
posPercentageSigns=regexp(gendesc,'%');
gendesc(posPercentageSigns)=[];
% Remove from string descri leading and trailing white spaces
gendesc=strtrim(gendesc);

if length(gendesc)<3
    htmlsitecont='';
else
    htmlsitecont=['<p>' gendesc '</p>'];
end



finsitecont=sprintf(['<h2 id="syntax">Syntax</h2>\r'...
    '<div class="syntax_signature">\r'...
    '<div class="syntax_signature_module">\r'...
    '<ul>\r']);

sitecont=[insnav inisitecont htmlsitecont finsitecont];

%% Create Sintax section of HTML file (referred to the loop of input arguments)
% Find point where first input argument starts

% Now find the number of required input
% arguments and if there are varargin
% and the number of required output
% arguments
% Find the string containing output arguments
%outargs= [out1, out2, out3]
[startIndexEq] = regexp(fstring,'=');

[startIndex] = regexp(fstring(1:startIndexEq(1)),'[');
[endIndex] = regexp(fstring(1:startIndexEq(1)),']');
% if startindex is empty it means there is a single output which is not
% enclosed in square brackets
if isempty(startIndex)
    outargs=['[' strtrim(fstring(9:startIndexEq(1)-1)) ']'];
else
    outargs=fstring(startIndex(1):endIndex(1));
end


% Find number of output arguments
% nargout = number of commas in string  outargs= [out1, out2, out3] +1
[commasOut] = regexp(outargs,',');

nargout=length(commasOut)+1;
if isempty(commasOut)
    commasOut=length(outargs);
end

% Required input arguments
% Find number of compulasory input arguments
% nargin = number of commas in string  functionname(inp1,inp2,inp3, ...)
[startIndexInp] = regexp(fstring,'(');
[endIndexInp] = regexp(fstring,')');
% inputargs =  (inp1,inp2,inp3, ...)
InputArgs=fstring(startIndexInp(1):endIndexInp(1));
% Check if inputargs contains the string varargin
[OptArgsVarargin]=regexp(InputArgs,'varargin');


[commasIn] = regexp(InputArgs,',');
j=1;

% nTOTargin= total number of input arguments (requested + optional),
% excluding name-value pairs arguments

if isempty(OptArgsVarargin)
    % Write all required input arguments in cell listargins
    nTOTargin=length(commasIn)+1;
else
    nTOTargin=length(commasIn);
end

listargins=cell(nTOTargin,1);
for i=1:nTOTargin
    if nTOTargin>1
        if i==1
            inpi=InputArgs(2:commasIn(i)-1);
        elseif i==nTOTargin && isempty(OptArgsVarargin)
            inpi=InputArgs(commasIn(i-1)+1:end-1);
            
        else
            inpi=InputArgs(commasIn(i-1)+1:commasIn(i)-1);
        end
    else
        if isempty(OptArgsVarargin) % if there are no optional input arguments
            inpi=InputArgs(2:end-1);
        else
            inpi=InputArgs(2:commasIn(i)-1);
        end
    end
    inpi=strtrim(inpi);
    listargins{i}=inpi;
end

% Check if among the input arguments there are explicitly declared inputs
% which are optionals. That is check if the intersection between the first
% columns of cell listOptArgs and vector listargins is empty
[OptArgsWithoutNameValue,PosOpt]=intersect(listargins(:,1),listOptArgs(:,1));
nOPTargin=length(OptArgsWithoutNameValue);
% nREQargin = number of required input arguments
nREQargin=nTOTargin-nOPTargin;

sintax=cell(nargout+1+nOPTargin,1);

if nOPTargin>0
    for j=1:nOPTargin;
        sintax{j}=[outargs(2:commasOut(1)-1) '=' name InputArgs(1:commasIn(PosOpt-1)-1) ')'];
    end
    j=j+1;
else
    j=1;
end


if isempty(OptArgsVarargin)
    sintax{j}=[outargs(2:commasOut(1)-1) '=' name InputArgs];
    j=j+1;
else
    % In this case there is also Name Value line
    strinputarg=strtrim(InputArgs(1:OptArgsVarargin-1));
    if strcmp(strinputarg(end),',')
        strinputarg=strinputarg(1:end-1);
    end
    
    %sintax{1}=[outargs(2:commasOut(1)-1) '=' name strtrim(inputargs(1:optargs1-2)) ')'];
    %sintax{2}=[outargs(2:commasOut(1)-1) '=' name strtrim(inputargs(1:optargs1-2)) ',Name,Value)'];
    
    sintax{j}=[outargs(2:commasOut(1)-1) '=' name strinputarg ')'];
    sintax{j+1}=[outargs(2:commasOut(1)-1) '=' name strinputarg ',Name,Value)'];
    
    j=j+2;
end
if j>1
    sintax=sintax(1:j-1);
end
%% Create Sintax section of HTML file (referred to the loop of output arguments)
if nargout>1
    for i=2:nargout
        if i<nargout
            sintax{j}=[outargs(1:commasOut(i)-1) ']=' name '(___)'];
        else
            sintax{j}=[outargs '=' name '(___)'];
        end
        j=j+1;
    end
end

sintaxhtml='';

for j=1:length(sintax)
    sintaxhtml= [sintaxhtml '<li><code class="synopsis">'  sintax{j} '</code>' ...
        '<span class="syntax_example">'...
        '<a class="intrnllnk" href="' name '.html#ex' num2str(j) '">' ...
        'example</a></span></li>\r'];
end
sintaxhtml=sprintf(sintaxhtml);

%  										<li><code class="synopsis">idx = kmeans(X,k)</code>
%  										<span class="syntax_example">
%  										<a class="intrnllnk" href="kmeans.html#ex1">
%  										example</a></span></li>


sintaxclose=sprintf(['</ul>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r']);


%% CREATE DESCRIPTION SECTION OF HTML FILE
inidescription=sprintf(['	<div class="ref_sect" itemprop="content">\r'...
    '							<h2 id="description">Description</h2>\r'...
    '							<div class="descriptions">\r'...
    '								<div class="description_module">\r']);

descriptionhtml='';
if strcmp(Display,'iter-detailed')
    disp('Examples')
    disp(sintax)
    disp('---------------')
end

% start of example j
[startIndexEx] = regexp(fstring,'%\{');
[endIndexEx] = regexp(fstring,'%\}');

%listEx = list which contains the examples (associated to sintax)
% First column= title, second column detailed description. Third column
% code
% Fourth colum dummy variable which indicates whether the example must be
% executed or not) If 1 example is executed
listEx=cell(length(sintax),4);

for j=1:length(sintax)
    
    descriptionini=sprintf(['<div class="description_element">\r'...
        '	<p class="syntax_example">\r'...
        '	<a class="intrnllnk" href="' name '.html#ex'  num2str(j) '">\r'...
        '	example</a></p>\r'...
        '	<p><span itemprop="syntax"><code>\r']);
    
    %For each element of input and output argument a hypertext link is added
    sintaxj=sintax{j};
    % Locate in expression [out1,out2,...]=namefunc(inp1,inp2,...) the
    % position of equal sign
    [startIndex] = regexp(sintaxj,'=');
    outs=sintaxj(1:startIndex-1);
    
    
    commaspos=regexp(outs,',');
    if isempty(commaspos);
        noutel=1;
    else
        noutel=length(commaspos)+1;
    end
    if j==length(sintax)
        % Write in cell listargouts the list of output arguments
        listargouts=cell(noutel,1);
    end
    outstring='';
    if noutel>1
        for i=1:noutel
            if i==1
                outi=['[' outs(2:commaspos(i))];
                outstring=[outstring sprintf(['[' '<a class="intrnllnk" href="#outputarg_' strtrim(outi(2:end-1)) '"><code>' outi(2:end-1) '</code></a>,\r'])];
                if j==length(sintax)
                    listargouts{i}=strtrim(outi(2:end-1));
                end
            elseif i==noutel
                outi=outs(commaspos(i-1)+1:end);
                outstring=[outstring sprintf(['<a class="intrnllnk" href="#outputarg_' strtrim(outi(1:end-1)) '"><code>' outi(1:end-1) '</code></a>]\r'])];
                if j==length(sintax)
                    listargouts{i}=strtrim(outi(1:end-1));
                end
            else
                outi=outs(commaspos(i-1)+1:commaspos(i));
                outstring=[outstring sprintf(['<a class="intrnllnk" href="#outputarg_' strtrim(outi(1:end-1)) '"><code>' outi(1:end-1) '</code></a>,\r'])];
                if j==length(sintax)
                    listargouts{i}=strtrim(outi(1:end-1));
                end
            end
            
        end
    else
        outi=strtrim(outs);
        outstring=sprintf(['<a class="intrnllnk" href="#outputarg_' outi '"><code>' outi '</code></a>\r']);
        if j==length(sintax)
            % TOCHECK
            % listargouts{j}=outi;
            listargouts=cell(1,1);
            listargouts{1}=outi;
        end
    end
    
    % Locate in  expression [out1,out2,...]=namefunc(inp1,inp2,...) the
    % position of open parenthesis sign
    [startIndex] = regexp(sintaxj,'(');
    inps=sintaxj(startIndex+1:end-1);
    commaspos=regexp(inps,',');
    if isempty(commaspos);
        ninpel=1;
    else
        ninpel=length(commaspos)+1;
    end
    inpstring='';
    if ninpel>1
        for i=1:ninpel
            if i==1
                inpi=inps(1:commaspos(i));
                inpi=strtrim(inpi);
            elseif i==ninpel
                inpi=[inps(commaspos(i-1)+1:end) ' '];
            else
                inpi=inps(commaspos(i-1)+1:commaspos(i));
                inpi=strtrim(inpi);
            end
            
            if (strcmp(inpi,'Name,') + strcmp(inpi,'Value'))>0
                inpstring=[inpstring sprintf('<a class="intrnllnk" href="#namevaluepairarguments"><code>Name, Value</code></a>\r')]; %#ok<*AGROW>
                break
            elseif  strcmp(inpi,'___') ==1
                inpstring=[inpstring sprintf([inpi '\r'])]; %#ok<*AGROW>
            else
                if i<ninpel
                    inpstring=[inpstring sprintf(['<a class="intrnllnk" href="#inputarg_' strtrim(inpi(1:end-1)) '"><code>' inpi(1:end-1) '</code></a>,\r'])]; %#ok<*AGROW>
                else
                    inpstring=[inpstring sprintf(['<a class="intrnllnk" href="#inputarg_' strtrim(inpi(1:end-1)) '"><code>' inpi(1:end-1) '</code></a>\r'])]; %#ok<*AGROW>
                end
                
            end
        end
    else
        inpi=inps;
        if strcmp(inpi,'___') ==1
            inpstring=sprintf([inpi '\r']);
        else
            inpstring=sprintf(['<a class="intrnllnk" href="#inputarg_' inpi '"><code>' inpi '</code></a>\r']);
        end
    end
    
    description=[outstring '=' name '(' inpstring ')'];
    %     description=sprintf(['<a class="intrnllnk" href="#outputarg_idx"><code>idx</code></a>\r'...
    %         '	= kmeans(<a class="intrnllnk" href="#inputarg_X"><code>X</code></a>\r'...
    %         '	,\r'...
    %         '   <a class="intrnllnk" href="#inputarg_k"><code>k</code></a>)\r']);
    
    %---------
    try
        stri=fstring(startIndexEx(j)+2:endIndexEx(j)-1);
    catch
        %disp(stri)
        warning('FSDA:wrongEx','This file does not contain enough examples, please add them!')
        stri='EXAMPLS TO ADD';
    end
    
    % What is before the first full stop is the title.
    % What is after the second full stop is the description
    % The first line which does not contain symbol % is the beginning of the
    % code
    [endtitle] = regexp(stri,'\.\s{1,3}','once');
    strititle=stri(1:endtitle);
    % Remove percentage signs if present.
    posPercentageSigns=regexp(strititle,'%');
    % Insert in fourth column of listEx information on the fact that the
    % example must be executed  (if there are two consecutive %% then it must be
    % executed)
    if length(posPercentageSigns)>1 && posPercentageSigns(2)-posPercentageSigns(1)==1
        listEx{j,4}=1;
    else
        listEx{j,4}=0;
    end
    
    strititle(posPercentageSigns)=[];
    % Remove from string strititle leading and trailing white spaces
    strititle=strtrim(strititle);
    % Store title
    listEx{j,1}=strititle;
    
    % Find point where description ends
    inicr=regexp(stri,'\r');
    % This is the first line which does not contain symbol %
    for jj=1:length(inicr)-1
        strtest=stri(inicr(jj):inicr(jj+1));
        if isempty(regexp(strtest,'%','once'));
            break
        end
    end
    findescriptionEx=inicr(jj);
    strdescrEx=stri(endtitle+1:findescriptionEx);
    % remove % signs from strdescrEx
    posPercentageSigns=regexp(strdescrEx,'%');
    strdescrEx(posPercentageSigns)=[];
    
    listEx{j,2}=strdescrEx;
    % listEx{j,3}=stri(findescriptionEx+1:end);
    
    StringWithLTandGT=stri(findescriptionEx+1:end);
    StringWithoutLTandGT=strrep(StringWithLTandGT,'<','&lt;');
    StringWithoutLTandGT=strrep(StringWithoutLTandGT,'>','&gt;');
    listEx{j,3}=StringWithoutLTandGT;
    
    %---------
    
    descriptionend=[sprintf('</code></span>\r') ...
        listEx{j,1} sprintf(['</p>\r'...
        '</div>'])];
    description=[descriptionini description descriptionend];
    descriptionhtml= [descriptionhtml description];
end


closedescription=sprintf(['								</div>\r'...
    '							</div>\r'...
    '							<div class="clear">\r'...
    '							</div>\r'...
    '						</div>']);

description=[inidescription descriptionhtml closedescription];

if strcmp(Display,'iter-detailed')
    disp('Detailed information about all the examples')
    disp(listEx)
end

%% CREATE EXAMPLES SECTION OF HTML FILE

% imgtemplate = iamge to include for the examples which can be executed
imgtemplate='<img alt="" src="images_help/M.gif" style="width: 12px; height: 12px"> ';

% the examples which are inside %{   %} are put here.
% The first sentence which end with a full stop is the title of the example
iniexamples=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<div class="examples">\r'...
    '<h2 id="examples">Examples</h2>\r'... % start of expandable examples
    '<div id="expandableExamples" class="expandableContent">\r'...
    '<p class="switch"><a class="expandAllLink"' ...
    'href="javascript:void(0);">expand all</a>' ...
    '</p>']);

exampleshtml='';
for j=1:length(sintax)
    
    if listEx{j,4}==1
        addimg=imgtemplate;
    else
        addimg='';
    end
    
    exampleshtml=[exampleshtml  sprintf(['<div id="example_' num2str(j) '" class="example_module expandableContent">\r'...
        '<div id="ex' num2str(j) '">\r'...
        '</div>\r'...
        '<h3 class="expand"><span>\r'...
        '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        '<span class="example_title">']) addimg listEx{j,1} sprintf(['</span></a></span></h3>\r'...
        '<div class="collapse">\r'...
        '<p>']) listEx{j,2} sprintf(['<div class="programlisting">\r'...
        '<div class="codeinput"><pre class="matlab-code">\r']) ...
        listEx{j,3} sprintf(['</pre></div>\r'...
        '</div>\r'])...
        sprintf(['</p>\r'...
        '</div>\r'...
        '</div>\r'])]; % close div id="example_j"
end

%% CREATE RELATED EXAMPLES SECTION OF HTML FILE
if length(startIndexEx)>length(sintax)
    lsintax=length(sintax);
    % Fourth column of listextraEX contains flag 1 or 0 depending on the
    % fact that the example must be executed or not
    listExtraEx=cell(length(startIndexEx)-lsintax,4);
    
    
    for j=1:size(listExtraEx,1)
        stri=fstring(startIndexEx(j+lsintax)+2:endIndexEx(j+lsintax)-1);
        % What is before the first full stop is the title.
        % What is after the second full stop is the description
        % The first line which does not contain symbol % is the beginning of the
        % code
        [endtitle] = regexp(stri,'\.','once');
        strititle=stri(1:endtitle);
        % Remove percentage signs if present.
        posPercentageSigns=regexp(strititle,'%');
        
        if length(posPercentageSigns)>1 && posPercentageSigns(2)-posPercentageSigns(1)==1
            listExtraEx{j,4}=1;
        else
            listExtraEx{j,4}=0;
        end
        
        
        strititle(posPercentageSigns)=[];
        % Remove from string strititle leading and trailing white spaces
        strititle=strtrim(strititle);
        % Store title
        listExtraEx{j,1}=strititle;
        
        % Find point where description ends
        inicr=regexp(stri,'\r');
        % This is the first line which does not contain symbol %
        for jj=1:length(inicr)-1
            strtest=stri(inicr(jj):inicr(jj+1));
            if isempty(regexp(strtest,'%','once'));
                break
            end
        end
        findescriptionEx=inicr(jj);
        strdescrEx=stri(endtitle+1:findescriptionEx);
        % remove % signs from strdescrEx
        posPercentageSigns=regexp(strdescrEx,'%');
        strdescrEx(posPercentageSigns)=[];
        
        listExtraEx{j,2}=strdescrEx;
        % Replace symbols < and > with &lt; and  &gt;
        StringWithLTandGT=stri(findescriptionEx+1:end);
        StringWithoutLTandGT=strrep(StringWithLTandGT,'<','&lt;');
        StringWithoutLTandGT=strrep(StringWithoutLTandGT,'>','&gt;');
        
        % listExtraEx{j,3}=stri(findescriptionEx+1:end);
        listExtraEx{j,3}=StringWithoutLTandGT;
    end
    
    if strcmp(Display,'iter-detailed')
        disp('Detailed information about all the Extra examples')
        disp(listExtraEx)
    end
else
    listExtraEx='';
end

closeexamples=sprintf(['</div>\r'... % close div id="expandableExamples
    '<p> </p>\r']);
% Related examples are below
iniRelatedExamples='';
RelatedExamples='';
if length(startIndexEx)>length(sintax)
    iniRelatedExamples=sprintf('<h3 class="bottom_ruled">Related Examples</h3>\r');
    
    for j=1:size(listExtraEx,1)
        
        if listExtraEx{j,4}==1
            addimg=imgtemplate;
        else
            addimg='';
        end
        
        RelatedExamples=[RelatedExamples  sprintf(['<div id="example_' num2str(j) '" class="example_module expandableContent">\r'...
            '<div id="ex' num2str(j) '">\r'...
            '</div>\r'...
            '<h3 class="expand"><span>\r'...
            '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            '<span class="example_title">']) addimg listExtraEx{j,1} sprintf(['</span></a></span></h3>\r'...
            '<div class="collapse">\r'...
            '<p>']) listExtraEx{j,2} sprintf(['<div class="programlisting">\r'...
            '<div class="codeinput"><pre class="matlab-code">\r']) ...
            listExtraEx{j,3} sprintf(['</pre></div>\r'...
            '</div>\r'])...
            sprintf(['</p>\r'...
            '</div>\r'...
            '</div>\r'])]; % close div id="example_j"
    end
    
end

closeallex=sprintf(['</div>\r'... % div class="examples"
    '</div>']);	 % close class="ref_sect

examples=[iniexamples exampleshtml closeexamples iniRelatedExamples...
    RelatedExamples closeallex];

%% CREATE REQUIRED NPUT ARGUMENTS SECTION OF HTML file
iniReqInputArgs=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<h2 id="inputs">Input Arguments</h2>\r'...
    '<div class="expandableContent">\r'...
    '<div class="arguments">\r'...
    '<div class="input_argument_container">\r'...
    '<p class="switch">\r'...
    '<a class="expandAllLink" href="javascript:void(0);">\r'...
    'expand all</a></p>']);


reqargs='';
%% Create listInptArgs and related HTML code
% listInpArgs = list which contains all required input arguments
% The first column will contain the names (just one word)
% The second column will contain the title of the input argument (the first
% sentence which finishes with a full stop sign)
% The third column will contain the type (the second sentence which
% finishes with a comma or full stop sign), e.g. scalar, matrix ...
% The fourth column will contain the long description. What starts with the
% third sentence
% The fifth column will contain the example what starts just after
% string  Example - (This is necessary just if the input argument is
% optional)
% The sixth column will contain the example what starts just after
% string  Data Types -

listInpArgs=cell(length(nTOTargin),6);

for i=1:nTOTargin
    
    % Name of the input argument (just one word)
    inpi=listargins{i};
    listInpArgs{i,1}=inpi;
    
    
    insel=regexp(fstring,'Required input arguments:');
    if isempty(insel)
        disp('Please check .m input file')
        error('FSDA:missInps','.m file does not contain ''Required input arguments:'' string')
    end
    
    % substring to search start from Required input arguments:
    fstringsel=fstring(insel(1):end);
    
    inipoint=regexp(fstringsel,['%\s{0,10}' listargins{i} '\s{0,10}:']);
    
    % The endpoint of the substring is See also or the next optional input argument
    if i <nREQargin
        endpoint=regexp(fstringsel,['%\s{0,10}' listargins{i+1} '\s{0,10}:']);
    elseif i==nREQargin
        endpoint=regexp(fstringsel,'Optional input arguments:');
        if isempty(endpoint)
            disp('Please check .m input file')
            error('FSDA:missOuts','.m file does not contain ''Optional input arguments:'' string')
        end
    elseif i<nTOTargin
        endpoint=regexp(fstringsel,['%\s{0,10}' listargins{i+1} '\s{0,10}:']);
    else
        endpoint=regexp(fstringsel,'Output:');
        if isempty(endpoint)
            disp('Please check .m input file')
            error('FSDA:missOuts','.m file does not contain ''Output:'' string')
        end
    end
    % DescrInputToSplit = string which contains all the information about the i-th input
    % argument (excluding xxxx :)
    % fstringseltmp = string which contains all the information about the i-th input
    % argument (including  xxxx :)
    fstringseltmp=fstringsel(inipoint(1)+1:endpoint(1)-1);
    inipoint=regexp(fstringseltmp,':');
    
    
    DescrInputToSplit=fstringseltmp(inipoint(1)+1:end);
    %    DescrInputToSplit=fstringseltmp((inipoint(1)+length(listargins{i})+2):endpoint(1)-1);
    
    % Remove from string descri all '% signs
    posPercentageSigns=regexp(DescrInputToSplit,'%');
    DescrInputToSplit(posPercentageSigns)=[];
    % Remove from string descri leading and trailing white spaces
    DescrInputToSplit=strtrim(DescrInputToSplit);
    %------------------
    % Add an artificial space character at the end just in case sentence
    % terminates with a full stop followed by no white space character
    % because in the next regexp we search for full stops followed by one
    % up to three spaces. At least one space is necessary otherwise we
    % misinterpret number as 0.234 (in this last case the full stop is a
    % decimal separatorand not the end of a sentence)
    DescrInputToSplit=[DescrInputToSplit ' '];
    
    [inifullstops]=regexp(DescrInputToSplit,'\.[\s1-3]');
    if isempty(inifullstops)
        warning('FSDA:publishFS:WrongInp',['Input option: ''' inpi '''\n Sentence''' DescrInputToSplit '''must contain at least two full stops'])
        % error('Wrong input')
    end
    shortdesc=strtrim(DescrInputToSplit(1:inifullstops(1)-1));
    % Store title of the i-th input argument
    %     % remove sign : if present at the beginning of the sentence
    if strcmp(shortdesc(1),':')
        shortdesc=strtrim(shortdesc(2:end));
    end
    listInpArgs{i,2}= shortdesc;
    
    try
        descrtype=strtrim(DescrInputToSplit(inifullstops(1)+1:inifullstops(2)-1));
    catch
        % warning('FSDA:publishFS:WrongInp',['Input: ' listInpArgs{i,1}])
        errmsg=['Input argument ' listInpArgs{i,1} ' Sentence ''' DescrInputToSplit ''' must contain at least two full stops followed by a white space'];
        error('FSDA:publishFS:WrongInp',errmsg)
    end
    
    listInpArgs{i,3}=descrtype;
    
    try
        descrlong=strtrim(DescrInputToSplit(inifullstops(2)+1:end));
    catch
        error('FSDA:publishFS:WrongInp',['Sentence''' DescrInputToSplit '''must contain at least two full stops'])
    end
    
    
    
    if i<=nREQargin
        
        Datatypes=regexp(descrlong,'Data Types -','once');
        if ~isempty(Datatypes)
            listInpArgs{i,6}=descrlong(Datatypes+13:end);
            
            descrlong=descrlong(1:Datatypes-1);
            descrlongHTML=formatHTML(descrlong);
            listInpArgs{i,4}=descrlongHTML;
        else
            
            descrlongHTML=formatHTML(descrlong);
            
            listInpArgs{i,4}=descrlongHTML;
            warning('FSDA:publishFS:MissingDataType',['Input argument ''' inpi ''' does not contain DataType line, by default string  ''single| double'' has been added'])
            
            listInpArgs{i,6}='single| double';
        end
        
        
        reqargs=[reqargs sprintf(['<div class="expandableContent">\r'...
            ' <div id="inputarg_' inpi '" class="clearfix">\r'...
            ' </div>\r'...
            ' <h3 id="input_argument_' inpi '" class="expand">\r'...
            ' <span>\r'...
            ' <a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            ' <span class="argument_name"><code>' inpi '</code> &#8212; ']) listInpArgs{i,2} sprintf([' </span> \r'...  % &#8212; = long dash
            ' </a><span class="example_desc">']) listInpArgs{i,3} sprintf(['</span></span></h3>\r'...
            ' <div class="collapse">\r'...
            ' <p>']) listInpArgs{i,4} sprintf(['</p>\r'...
            ' <p class="datatypelist"><strong>\r'...
            ' Data Types: </strong><code>' listInpArgs{i,6}  '</code></p>\r'...
            ' </div>\r'...
            ' </div>\r'])];
    else
        
        
        % Check if descrlong contains
        % Example - and Data types -
        
        CheckExample=regexp(descrlong,'Example -','once');
        if ~isempty(CheckExample)
            Datatypes=regexp(descrlong,'Data Types -','once');
            descrlonginp=descrlong(1:CheckExample-1);
            descrlongHTML=formatHTML(descrlonginp);
            listInpArgs{i,4}=descrlongHTML;
            
            % The first word of example code must be embedded around tags <code> </code>
            examplecode=descrlong(CheckExample+10:Datatypes-1);
            posspace=regexp(examplecode,'      ');
            examplecode=['<code>' examplecode(1:posspace-1) '</code>' examplecode(posspace:end)];
            listInpArgs{i,5}=strtrim(examplecode);
            listInpArgs{i,6}=descrlong(Datatypes+13:end);
            
            
        else
            listInpArgs{i,4}=descrlong;
            warning('FSDA:publishFS:MissingExample',['Optional input argument ''' inpi ''' does not contain an Example'])
            
        end
        
        
        reqargs=[reqargs sprintf(['<div class="expandableContent">\r'...
            ' <div id="inputarg_' inpi '" class="clearfix">\r'...
            ' </div>\r'...
            ' <h3 id="input_argument_' inpi '" class="expand">\r'...
            ' <span>\r'...
            ' <a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            ' <span class="argument_name"><code>' inpi '</code> &#8212; ']) listInpArgs{i,2} sprintf([' </span> \r'...  % &#8212; = long dash
            ' </a><span class="example_desc">']) listInpArgs{i,3} sprintf(['</span></span></h3>\r'...
            ' <div class="collapse">\r'...
            ' <p>']) listInpArgs{i,4} sprintf(['</p>\r'...
            '	<p class="description_valueexample">\r'...
            '       <strong>Example: </strong>' listInpArgs{i,6} '</p>\r'...
            ' <p class="datatypelist"><strong>\r'...
            ' Data Types: </strong><code>' listInpArgs{i,5}  '</code></p>\r'...
            ' </div>\r'...
            ' </div>\r'])];
        
        
    end
    if i==nREQargin && i<nTOTargin
        reqargs = [reqargs sprintf(['<div id="optionalarguments" class="clearfix">\r'...
            '</div>\r' ...
            '<h3 id="namevaluepairs" class="bottom_ruled">\r'...
            'Optional Arguments</h3>'])];
    end
end


%% CREATE Optional Arguments SECTION OF HTML FILE (excluding Name-Value pair)


if strcmp(Display,'iter-detailed')
    disp('Detailed information about Input arguments')
    disp(listInpArgs)
end


%% CREATE Name-Value Pair Arguments SECTION OF HTML FILE
% codewithexample=['''Distance'',''cosine'',''Replicates'',10,' ...
%     '''Options'',statset(''UseParallel'',1)'];
if isempty(OptArgsVarargin)
    OptArgsNameValueHeading='';
    OptArgsNameValue='';
    
else
    codewithexample='';
    for i=1:size(listOptArgs,1)
        
        NamVali=listOptArgs{i,5};
        if isempty(NamVali)
            error('FSDA:missingex',['Optional input argument  ' listOptArgs{i,1} ...
                ' does not seem to contain an example (or alternatively string remark has not been put at the end)'])
        end
        % Add as example only those which do finish with </code>, that is just
        % those which do not contain exaplanations
        if strcmp('</code>',NamVali(end-6:end))
            if i<size(listOptArgs,1)
                codewithexample=[codewithexample NamVali ',' ];
            else
                codewithexample=[codewithexample NamVali];
            end
        end
    end
    OptArgsNameValueHeading=sprintf(['<div id="namevaluepairarguments" class="clearfix">\r'...
        '</div>\r' ...
        '<h3 id="namevaluepairs" class="bottom_ruled">\r'...
        'Name-Value Pair Arguments</h3>\r'...
        '<div class="namevaluecontainer">\r'...
        '<p>Specify optional comma-separated pairs of <code>Name,Value</code> arguments.\r'...
        ' <code>Name</code> is the argument name and <code>Value</code>\r'...
        ' is the corresponding value. <code>Name</code> must appear \r'...
        ' inside single quotes (<code>'' ''</code>). \r'...
        ' You can specify several name and value pair arguments in any order as <code> \r'...
        ' Name1,Value1,...,NameN,ValueN</code>.</p> \r'...
        ' <span class="example_desc"><strong>Example:\r'...
        '</strong><code>' codewithexample '</code>\r'...
        '</span></div>']);
    
    
    OptArgsNameValue='';
    for i=1:size(listOptArgs,1);
        nameoptarg=listOptArgs{i,1};
        titloptarg=listOptArgs{i,2};
        shortdesc=listOptArgs{i,3};
        if strcmp(shortdesc,'Structure') && ~isempty(strfind(listOptArgs{i,4},'field'))
            longdesc=listOptArgs{i,4};
            
            [inistructfield,finstructfield]=regexp(longdesc,'\s{8,18}\w*\s{0,8}=');
            posREMARK=regexp(longdesc,'REMARK','once');
            % If there is the string REMARK it means that all that which is after
            % remark is a general statement which does not have to be contained in the
            % table with the associated fields
            if ~isempty(posREMARK)
                boo=inistructfield<posREMARK;
                inistructfield=inistructfield(boo);
                finstructfield=finstructfield(boo);
            end
            % Insert all fields of inside listStruArgs
            listOutArgs=cell(length(inistructfield),2);
            for k=1:length(inistructfield)
                % fin(i)-1 because character ':' does not have to be extracted
                StructFieldName=longdesc(inistructfield(k):finstructfield(k)-1);
                % Remove from string opti leading and trailing white spaces
                StructFieldName=strtrim(StructFieldName);
                
                % Store name in the first column
                listOutArgs{k,1}=StructFieldName;
                % Store short description in the 3nd col of listOptArgs
                if k<length(inistructfield)
                    StructFieldCont=longdesc(finstructfield(k)+1:inistructfield(k+1));
                else
                    if ~isempty(posREMARK)
                        StructFieldCont=longdesc(finstructfield(k)+1:posREMARK-1);
                    else
                        StructFieldCont=longdesc(finstructfield(k)+1:end);
                    end
                end
                % Store name in the first column
                listOutArgs{k,2}=StructFieldCont;
            end
            
            
            
            Tablehtml='';
            for k=1:length(inistructfield)
                Tablehtml=[Tablehtml sprintf(['<tr valign="top">\r'...
                    '<td><code>' listOutArgs{k,1} '</code></td>\r'...
                    '<td>\r'...
                    '<p>']) listOutArgs{k,2} sprintf(['</p>\r'...
                    '</td>\r'...
                    '</tr>'])];
            end
            
            % Add the Remark after the table, if it is present
            if ~isempty(posREMARK)
                descrREMARK=longdesc(posREMARK:end);
                descrREMARKHTML=formatHTML(descrREMARK);
                
                longdescription=[iniTable Tablehtml cloTable '<p>' descrREMARKHTML '</p>'];
            else
                longdescription=[iniTable Tablehtml cloTable];
            end
        else
            longdescriptionHTML=formatHTML(listOptArgs{i,4});
            
            longdescription=longdescriptionHTML;
        end
        % datatype = type of data for that particular option
        %     examplecode=['''Display'',''final'''];
        %     datatype='char';
        
        examplecode=listOptArgs{i,5};
        datatype=listOptArgs{i,6};
        
        OptArgsNameValue=[OptArgsNameValue sprintf(['<div class="expandableContent">\r'...
            '<div id="inputarg_Display" class="clearfix">\r'...
            '</div>\r'...
            '<h3 id="input_argument_namevalue_display" class="expand">\r'...
            '<span>\r'...
            '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            '<span class="argument_name"><code>' nameoptarg  '</code> \r'...
            '&#8212;']) titloptarg sprintf(['</span></a><span class="example_desc">' shortdesc  '</span></span></h3>\r'...
            '<div class="collapse">\r'...
            '	<p>']) longdescription sprintf(['</p>\r'...
            '	<p class="description_valueexample">\r'...
            '       <strong>Example: </strong>' examplecode '</p>\r'...
            '	<p class="datatypelist"><strong>Data Types: </strong><code>' datatype '</code></p>\r'...
            '</div>\r'...
            '</div>'])];
    end
    % CLOSE OPT ARGS NAME VALUE
end

closeinputargs=sprintf(['</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r']);
InputArgs=[iniReqInputArgs reqargs  OptArgsNameValueHeading OptArgsNameValue closeinputargs];


%% CREATE OUTPUT ARGUMENTS SECTION OF HTML FILE

inioutargs=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<h2>Output Arguments</h2>\r'...
    '<div class="expandableContent">\r'...
    '<div class="arguments">\r'...
    '<div class="output_argument_container">\r'...
    '<p class="switch">\r'...
    '<a class="expandAllLink" href="javascript:void(0);">expand all</a></p>']);

% outargs = strings which contains output arguements (includeing [])
% nargout = number of output arguments
outargshtml='';

% check if the last element of listargouts is varargout
% listargouts
%  Optional Output:


for i=1:nargout
    
    % listargouts is a cell which contains the list of output arguments
    outi=listargouts{i};
    
    outsel=regexp(fstring,'Output:');
    if isempty(outsel)
        disp('Please check HTML input file')
        error('FSDA:missOuts','HTML file does not contain ''Output:'' string')
    end
    
    % substring to search starting from Output:
    fstringsel=fstring(outsel(1):end);
    
    % The initial point of the string is 'listargouts{i}' is there is just
    % one output else is string 'listargouts{i} :' is there is more than
    % one output and this is not varargout
    % else if there is varargour the initialpoint is the string
    % "Optional Output:"
    if length(listargouts)==1
        inipoint=regexp(fstringsel,listargouts{i});
    elseif  i<length(listargouts)
        inipoint=regexp(fstringsel,[listargouts{i} '\s{0,7}:']);
    else
        if strcmp(listargouts{end},'varargout') ==0
            inipoint=regexp(fstringsel,[listargouts{i} '\s{0,7}:']);
        else
            inipoint=regexp(fstringsel,'Optional Output:')+8;
        end
    end
    
    if isempty(inipoint)
        error('FSDA:missingOuts',['Output argument ' listargouts{i} ' has not been found'])
    end
    
    % The endpoint of the substring is 'more About'. or sSee also or the next output argument
    if i <nargout-1
        endpoint=regexp(fstringsel,[listargouts{i+1} '\s{0,7}:']);
    elseif i==nargout-1
        
        if strcmp(listargouts{end},'varargout') ==0
            endpoint=regexp(fstringsel,[listargouts{i+1} '\s{0,7}:']);
        else
            % In this case there are also optional arguments
            endpoint=regexp(fstringsel,'Optional Output:');
        end
        
    else
        
        inipointSeeAlso=regexp(fstringsel,'See also','once');
        
        if isempty(inipointSeeAlso)
            disp('Please check .m input file')
            error('FSDA:missOuts','Input .m file does not contain ''See also:'' string')
        end
        
        inipointMoreAbout=regexp(fstringsel,'More About:','once');
        if ~isempty(inipointMoreAbout);
            MoreAbout=fstringsel(inipointMoreAbout+15:inipointSeeAlso-1);
            posPercentageSigns=regexp(MoreAbout,'%');
            MoreAbout(posPercentageSigns)=[];
            MoreAboutHTML=formatHTMLwithMATHJAX(MoreAbout);
        else
            MoreAboutHTML='';
            inipointMoreAbout=Inf;
        end
        
        endpoint=min(inipointSeeAlso,inipointMoreAbout);
        
    end
    
    % descri = string which contains the description of i-th output
    % argument
    try
        descrioutput=fstringsel((inipoint(1)+length(listargouts{i})+2):endpoint(1)-1);
        if isempty(descrioutput)
            initmp=inipoint(1);
            disp('Starting point of the description')
            disp([fstringsel(initmp:initmp+50) '.....'])
            disp('Final point of the description')
            endtmp=endpoint(1);
            disp(['....' fstringsel(endtmp-50:endtmp)]);
            disp(['FSDA:WrongOut','Could not process correctly output argument ' listargouts{i}])
        end
    catch
        disp(['FSDA:WrongOut','Could not process correctly output argument ' listargouts{i}])
    end
    
    % Remove from string descri all '% signs
    posPercentageSigns=regexp(descrioutput,'%');
    descrioutput(posPercentageSigns)=[];
    % Remove from string descri leading and trailing white spaces
    descrioutput=strtrim(descrioutput);
    if strcmp(descrioutput(1),':')
        descrioutput=strtrim(descrioutput(2:end));
    end
    
    % Check if the output is a structure. If this is the case
    checkifstructure=regexp(descrioutput,[outi '\.\w'],'once');
    
    
    if ~isempty(checkifstructure)
        [ini,fin]=regexp(descrioutput,['\s{0,8}' outi '\.\w*\s{0,8}=']);
        if isempty(ini)
            disp('Probably ":" symbols  must be replaced with "=" symbols in out description')
            error('FSDA:MissingArg',['Parser cannot find string \n''' outi '.xxxx'' = \n for output structure ' outi])
        else
        end
        
        listOutArgs=cell(length(ini),2);
        
        for k=1:length(ini)
            % fin(i)-1 because character ':' does not have to be extracted
            StructFieldName=descrioutput(ini(k):fin(k)-1);
            % Remove from string opti leading and trailing white spaces
            StructFieldName=strtrim(StructFieldName);
            StructFieldName=StructFieldName(length(outi)+2:end);
            
            % Store name in the first column
            listOutArgs{k,1}=StructFieldName;
            % Store short description in the 3nd col of listOptArgs
            if k<length(ini)
                StructFieldCont=descrioutput(fin(k)+1:ini(k+1)-1);
            else
                StructFieldCont=descrioutput(fin(k)+1:end);
            end
            
            % Store name in the first column
            listOutArgs{k,2}=StructFieldCont;
        end
        
        % rowtodel = vector which contains the duplicate rows of
        % listStruArgs which have to be deleted
        inisel=1:size(listOutArgs,1);
        
        % Check if cell listStruArgs contains duplicates in the first column
        for j=2:size(listOutArgs,1)
            if strcmp(listOutArgs{j,1},listOutArgs{j-1,1})
                listOutArgs{j-1,2}=[listOutArgs{j-1,2} listOutArgs{j,1} listOutArgs{j,2}];
                inisel(j)=999;
            end
        end
        
        if strcmp(Display,'iter-detailed')
            disp('Detailed information about Output arguments')
            disp(listOutArgs)
        end
        
        % remove from inisel the rows equal to 999 (that is the rows which
        % correspond to duplicated arguments)
        inisel(inisel==999)=[];
        
        Tablehtml='';
        for k=inisel % length(ini)
            
            descrlong=listOutArgs{k,2};
            
            descrlongHTML=formatHTML(descrlong);
            
            % listOutArgs{k,2}
            Tablehtml=[Tablehtml sprintf(['<tr valign="top">\r'...
                '<td><code>' listOutArgs{k,1} '</code></td>\r'...
                '<td>\r'...
                '<p>']) descrlongHTML  sprintf(['</p>\r'...
                '</td>\r'...
                '</tr>'])];
        end
        
        descrioutput=[iniTable Tablehtml cloTable];
        
        preamble='A structure containing the following fields:';
        
        outargshtml=[outargshtml sprintf(['<div class="expandableContent">\r'...
            '<div id="outputarg_' outi '" class="clearfix">\r'...
            '</div>\r'...
            '<h3 id="output_argument_' outi '" class="expand">\r'...
            '<span>\r'...
            '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            '<span class="argument_name"><code>' outi '</code> &#8212; description</span></a>\r'...
            '<span class="example_desc">' preamble '</span></span></h3>\r'...
            '<div class="collapse">\r'...
            '<p>']) descrioutput sprintf(['</p>\r'...
            '</div>\r'...
            '</div>'])];
        
    else
        % Check if string descrioutput contains the words 'which contains' or
        % 'containing'; in the first 'numcharcontains'
        poswhichcontains=regexp(descrioutput,'which contains');
        poscontaining=regexp(descrioutput,'containing');
        numcharcontains=50;
        
        if ~isempty(poswhichcontains) && poswhichcontains(1)<numcharcontains
            preamble=descrioutput(1:poswhichcontains(1)-1);
            descrioutput=descrioutput(poswhichcontains(1)+14:end);
            % Remove word the at the beginning of the sentence and starts with
            % uppercase
            StartsWithThe=regexp(descrioutput,'the');
            if ~isempty(StartsWithThe)
                if StartsWithThe(1)<4
                    descrioutput=descrioutput(StartsWithThe(1)+4:end);
                end
            end
            descrioutput=strtrim(descrioutput);
            descrioutput=[upper(descrioutput(1)) descrioutput(2:end)];
        elseif ~isempty(poscontaining) && poscontaining(1)<numcharcontains
            preamble=descrioutput(1:poscontaining(1)-1);
            descrioutput=descrioutput(poscontaining(1)+10:end);
            % Remove word the at the beginning of the sentence and starts with
            % uppercase
            StartsWithThe=regexp(descrioutput,'the');
            if StartsWithThe(1)<4
                descrioutput=descrioutput(StartsWithThe(1)+4:end);
            end
            descrioutput=strtrim(descrioutput);
            descrioutput=[upper(descrioutput(1)) descrioutput(2:end)];
        else
            posfullstops=regexp(descrioutput,'\.[1-3\s]');
            if length(posfullstops)<2
                warn1=[' Wrong format for ouptut argument ' outi '\n'];
                error('FSDA:publishFS:WrongOut',[ warn1 'Sentence ''' descrioutput ''' must contain at least two full stops'])
                
            else
                preamble=descrioutput(posfullstops(1)+1:posfullstops(2)-1);
                descrioutput=[descrioutput(1:posfullstops(1)) descrioutput(posfullstops(2)+1:end)];
            end
        end
        
        % From
        posfullstop=regexp(descrioutput,'\.', 'once');
        if ~isempty(posfullstop)
            descroutputtitl=descrioutput(1:posfullstop);
            if length(descrioutput)> posfullstop
                descrioutput=descrioutput(posfullstop+1:end);
            else
                descrioutput='';
            end
        else
            descroutputtitl='FULL STOP MISSING IN THE OUTPUT DESCRIPTION';
        end
        
        % transform x with by and write in italic the dimensions of the
        % matrices
        if ~strcmp(preamble,'TOWRITE')
            
            outvect=regexp(preamble,'\wector', 'once');
            if ~isempty(outvect)
                beforepreamble=preamble(1:outvect-1);
                beforepreamble=strrep(beforepreamble, 'x', '-by-');
                preamble=['<code>' beforepreamble  '</code>' preamble(outvect:end)];
            end
            
            
            outmat=regexp(preamble,'\watrix', 'once');
            if ~isempty(outmat)
                beforepreamble=preamble(1:outmat-1);
                beforepreamble=strrep(beforepreamble, 'x', '-by-');
                preamble=['<code>' beforepreamble '</code>' preamble(outmat:end)];
            end
        end
        
        preamble=strtrim(preamble);
        
        outargshtml=[outargshtml sprintf(['<div class="expandableContent">\r'...
            '<div id="outputarg_' outi '" class="clearfix">\r'...
            '</div>\r'...
            '<h3 id="output_argument_' outi '" class="expand">\r'...
            '<span>\r'...
            '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            '<span class="argument_name"><code>' outi '</code> &#8212;']) descroutputtitl   sprintf(['</span></a>\r'...
            '<span class="example_desc">']) preamble sprintf(['</span></span></h3>\r'...
            '<div class="collapse">\r'...
            '<p>']) descrioutput sprintf(['</p>\r'...
            '</div>\r'...
            '</div>'])];
        
    end
    
end

closeoutargs=sprintf(['</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>']);


outargs=[inioutargs outargshtml closeoutargs];

%% CREATE MORE ABOUT SECTION

if ~isempty(MoreAboutHTML)
    Moreabout=[sprintf(['<div class="moreabout ref_sect">\r'...
        '<h2 id="moreabout">More About</h2>\r'...
        '<div class="expandableContent">\r'...
        '<p class="switch">\r'...
        '<a class="expandAllLink" href="javascript:void(0);">\r'...
        'expand all</a></p>\r'...
        '<div class="expandableContent" itemprop="content">\r'...
        '<h3 class="expand"><span>\r'...
        '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        '<span>Methodological Details </span></a></span></h3>\r'...
        '<div class="collapse">\r'...
        '<p>']) MoreAboutHTML sprintf(['</p>\r'...
        '</div>\r'...
        '</div>\r'...
        '</div>\r'...
        '</div>'])];
else
    Moreabout='';
end

%% REFERENCES
inipointAcknowledgements=regexp(fstring,'Acknowledgements:');
if isempty(inipointAcknowledgements);
    Acknowledgements='';
end

iniref=regexp(fstring,'References:');
if isempty(iniref)
    warning('FSDA:missInp',['File ' name '.m does not contain string ''References:'''])
    % refsargs={''};
    References='';
    
else
    
    
    
    inipointCopyright=regexp(fstring,'Copyright');
    
    
    if ~isempty(inipointAcknowledgements);
        Acknowledgements=fstring(inipointAcknowledgements+18:inipointCopyright-1);
        posPercentageSigns=regexp(Acknowledgements,'%');
        Acknowledgements(posPercentageSigns)=[];
        
    else
        Acknowledgements='';
        inipointAcknowledgements=Inf;
    end
    
    endref=min(inipointAcknowledgements,inipointCopyright);
    
    % stringsel = block of test which contains the references
    fstringsel=fstring(iniref(1)+1:endref(1)-1);
    
    
    % Now we must try to infer how many references there are, that is where
    % each reference starts and ends
    % refsargs is a cell which contains in each row the references
    refsargs=cell(10,1);
    ij=1;
    findnewline=regexp(fstringsel,'\n');
    begref=0;
    endref=0;
    begreftoaddtmp='';
    for i=1:length(findnewline)-1
        % Find candidate for beginning of a reference
        candiniref=fstringsel(findnewline(i):findnewline(i+1));
        % findref=regexp(candiniref,'\(....\)','once');
        % Find 4 or 5 characters inside parenthesis. Note that 4 or 5
        % because the year can be of type (2002) or for example (2002b)
        findref=regexp(candiniref,'\(.{4,5}\)','once');
        
        if ~isempty(findref) && begref==0
            begreftoadd=findnewline(i);
            begref=1;
        elseif ~isempty(findref) && begref==1
            endreftoadd=findnewline(i)-1;
            begreftoaddtmp=findnewline(i);
            endref=1;
        else
        end
        
        if (begref==1 && endref ==1) || (begref==1 && endref==0 && i==length(findnewline)-1)
            if i<length(findnewline)-1
                ref2add=fstringsel(begreftoadd:endreftoadd);
            else
                ref2add=fstringsel(begreftoadd:findnewline(i));
            end
            
            % Remove % characters and white spaces
            posPercentageSigns=regexp(ref2add,'%');
            ref2add(posPercentageSigns)=[];
            ref2add=strtrim(ref2add);
            refsargs{ij}=ref2add;
            ij=ij+1;
            % begref=0;
            endref=0;
            begreftoadd=begreftoaddtmp;
        end
    end
    
    % Now check if there is a final open reference
    refsargs=refsargs(1:ij-1);
    
    
    Referenceshtml='';
    iniReferences=sprintf(['<div class="ref_sect" itemprop="content">\r'...
        '<div class="bibliography">\r'...
        '<h2 id="references">References</h2> \r']);
    
    for i=1:length(refsargs)
        Referenceshtml=sprintf([Referenceshtml  '<div><p>' refsargs{i} '</p></div>\r']);
    end
    Referencesclose=sprintf(['</div>\r'...
        '</div>']);
    References=[iniReferences Referenceshtml Referencesclose];
end


%% ACKNOWLEDGEMENTS
if ~isempty(Acknowledgements)
    iniAcknowledgements=sprintf(['<div class="ref_sect" itemprop="content">\r'...
        '<div class="bibliography">\r'...
        '<h2 id="references">Acknowledgements</h2> \r']);
    
    Acknowledgementshtml=sprintf(['<div><p>' Acknowledgements '</p></div>\r']);
    
    Ack=[iniAcknowledgements Acknowledgementshtml Referencesclose];
else
    Ack='';
end


%% SEE ALSO

iniSeealso=sprintf(['<div class="ref_sect">\r'...
    '<h2>See Also</h2>\r'...
    '<p>\r']);

iniref=regexp(fstring,'See also','once');
endref=regexp(fstring(iniref:end),'%','once');
seealsostr=fstring(iniref+8:iniref+endref-2);

% Remove : character and % character
posColonSign=regexp(seealsostr,':');
seealsostr(posColonSign)=[];
posPercentageSigns=regexp(seealsostr,'%');
seealsostr(posPercentageSigns)=[];
seealsostr=strtrim(seealsostr);

% count number of see also
poscommas=regexp(seealsostr,',');
nseealso=length(poscommas)+1;

Seealsohtml='';
for i=1:nseealso
    if nseealso==1;
        Seealsoitem= seealsostr(1:end);
    else
        if i==nseealso
            Seealsoitem=seealsostr(poscommas(i-1)+1:end);
        elseif i==1
            Seealsoitem=seealsostr(1:poscommas(i)-1);
        else
            Seealsoitem=seealsostr(poscommas(i-1)+1:poscommas(i)-1);
        end
    end
    % Remove .m if present at the end of the reference
    if ~isempty(Seealsoitem) &&  strcmp(Seealsoitem(end-1:end),'.m')
        Seealsoitem=Seealsoitem(1:end-2);
    end
    
    Seealsohtml=[Seealsohtml sprintf(['<span itemprop="seealso">\r'...
        '<a href="' Seealsoitem '.html" itemprop="url">\r'...
        '<span itemprop="name"><code>' Seealsoitem '</code></span></a></span>\r'])];
    
    if i < nseealso
        Seealsohtml=[Seealsohtml ' | '];
    end
end

closeSeealso=sprintf('</div>\r');

Seealso=[iniSeealso Seealsohtml closeSeealso];
% Seealso='';

%% CLOSE TAGS SECTION

clos=sprintf(['<h1>Automatically generated by routine publishFS</h1>\r'...
    '</div>\r'...
    '</section>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r']);


insbarra=sprintf(['<script type ="text/javascript" language="javaScript">\r'...
    'document.write(barra);\r'...
    '</script>']);


closbody=sprintf(['</body>\r'...
    '</html>']);

fclose(fileID);

outstring=([titl metacontent sitecont sintaxhtml sintaxclose description  ....
    examples InputArgs outargs Moreabout References Ack Seealso clos insnav insbarra closbody]);

file1ID=fopen([outputDir '\' name 'tmp.html'],'w');

if file1ID==-1
    outputDir=strrep(outputDir,'\','\\');
    errmsg= [' Path ' outputDir '\\' name '.html does not exist or output file '  name '.html is not writable'];
    error('FSDA:publishFS:WrngOutFolder',errmsg)
end

if evalCode==true
    %% EXECUTE THE EXAMPLES WHICH START WITH SYMBOLS %%
    % Create a temporary file with all the examples which must be executed
    % ExToExec= string which contains the examples which must be executed
    ExToExec='';
    numexToExec=0;
    for i=1:size(listEx,1)
        if listEx{i,4}==1
            ExToExec=[ExToExec '%% Ex' num2str(i) listEx{i,3}];
            numexToExec=numexToExec+1;
        end
    end
    
    numextraexToExec=0;
    if ~isempty(listExtraEx)
        for i=1:size(listExtraEx,1)
            if listExtraEx{i,4}==1
                ExToExec=[ExToExec '%% Ex' num2str(i) listExtraEx{i,3}];
                numextraexToExec=numextraexToExec+1;
            end
        end
    end
    
    if numextraexToExec+numexToExec>0
        % tmp .file containing all the .m examples will be created. It will be
        % created in subfolder tmp of helpfiles and then automatically removed.
        % This subfolder will be added to put for this session
        nametmp=[name 'tmp.m'];
        fullPathToScript=[pathstr '\helpfiles\FSDA\tmp\' nametmp];
        filetmp=fopen(fullPathToScript,'w');
        addpath([pathstr '\helpfiles\FSDA\tmp'])
        addpath([pathstr '\utilities\privateFS'])
        
        % Replace < and > HTML symbols with < and >
        ExToExec=strrep(ExToExec,'&lt;','<');
        ExToExec=strrep(ExToExec,'&gt;','>');
        
        fprintf(filetmp,'%s',ExToExec);
        fclose(filetmp);
        
        options=struct;
        options = supplyDefaultOptions(options);
        options.codeToEvaluate=[name 'tmp'];
        options.createThumbnail=0;
        [dom,cellBoundaries] = m2mxdom(ExToExec);
        prefix=[name 'tmp'];
        % file='C:\Users\MarcoAW\D\matlab\FSDA\examples\tmp.m';
        dom = evalmxdom(fullPathToScript,dom,cellBoundaries,prefix,imagesDir,outputDir,options);
        %
        ext='html';
        AbsoluteFilename = fullfile(outputDir,[prefix '.' ext]);
        [xResultURI]=xslt(dom,options.stylesheet,AbsoluteFilename);
        
        % Now remove the temporary .m file with the examples which had been created
        delete(fullPathToScript)
        
        % load html output in a string and extract the parts which are required
        fileHTML = fopen(xResultURI(7:end), 'r+');
        % Insert the file into fstring
        fstringHTML=fscanf(fileHTML,'%c');
        
        totex=numexToExec+numextraexToExec;
        texttoadd=cell(totex,1);
        
        fHTML=regexp(fstringHTML,'<h2>Ex');
        if isempty(fHTML)
            fHTML=regexp(fstringHTML,'<pre class="codeoutput">','once');
        end
        % If fHTML is still empty it means that the ouptu only generates images
        if isempty(fHTML)
            fHTML=regexp(fstringHTML,'<img vspace','once');
        end
        
        for j=1:totex
            if j<totex && totex>1
                fcand=fstringHTML(fHTML(j):fHTML(j+1)-1);
            else
                fendHTML=regexp(fstringHTML,'<p class="footer">','once');
                fcand=fstringHTML(fHTML(end):fendHTML-1);
            end
            
            % in fcand search the two following strings
            fcode=regexp(fcand,'<pre class="codeoutput">','once');
            if isempty(fcode)
                fcode=Inf;
            end
            fimg=regexp(fcand,'<img','once');
            if isempty(fimg)
                fimg=Inf;
            end
            if min(fcode,fimg)<Inf
                texttoadd{j}=fcand(min(fcode,fimg):end);
            end
        end
        
        % Now insert the strings which have been stored in cell texttoadd in the
        % appropriate position of outstring
        a=cell2mat(listEx(:,4));
        seqa=1:length(a);
        sel=seqa(a==1);
        
        for i=1:length(sel)
            % Process string listEx{i,1}
            listExi=listEx{sel(i),1};
            % If there are signs $ ^ [ ] replace them with \$ and \^ \[ \]
            listExi=strrep(listExi,'$','\$');
            listExi=strrep(listExi,'^','\^');
            listExi=strrep(listExi,'[','\[');
            listExi=strrep(listExi,']','\]');
            listExi=strrep(listExi,'(','\(');
            listExi=strrep(listExi,')','\)');
            
            iniout=regexp(outstring,listExi);
            if length(iniout)<2
                errmsg= [' Title of example \n''' listExi '''\n could not be found \n'...
                    'Probably because the string contains special characters\n' ...
                    'which cannot be interpreted by MATLAB function regexp'];
                error('FSDA:publishFS:WrngOutFolder',errmsg)
            end
            
            finout=regexp(outstring,'</pre>');
            finout=finout(finout>iniout(2));
            % outstring(finout:finout+11)
            % inclplint = point where output of the example must be included
            inclpoint=finout(1)+18;
            % incl= string which contains the output of the code
            incl=texttoadd{i};
            outstring=[outstring(1:inclpoint) incl outstring(inclpoint+1:end)];
        end
        
        if ~isempty(listExtraEx)
            a=cell2mat(listExtraEx(:,4));
            seqa=1:length(a);
            sel=seqa(a==1);
            
            for i=1:length(sel)
                % Process string listEx{i,1}
                listExi=listExtraEx{sel(i),1};
                % If there are signs $ ^ [ ] replace them with \$ and \^ \[ \]
                listExi=strrep(listExi,'$','\$');
                listExi=strrep(listExi,'^','\^');
                listExi=strrep(listExi,'[','\[');
                listExi=strrep(listExi,']','\]');
                listExi=strrep(listExi,'(','\(');
                listExi=strrep(listExi,')','\)');
                listExi=strrep(listExi,'.','\.');
                
                iniout=regexp(outstring,listExi);
                if length(iniout)>2
                    disp(['Duplicate name for: ' listExi ' found'])
                    warning('FSDA:WrongArg','There are examples with the same title, please use a title which is unique')
                elseif isempty(iniout)
                    errmsg= [' Title of example \n''' listExi '''\n could not be found \n'...
                        'Probably because the string contains special characters\n' ...
                        'which cannot be interpreted by MATLAB function regexp'];
                    error('FSDA:publishFS:WrngEx',errmsg)
                else
                end
                
                iniout=iniout(1);
                
                finout=regexp(outstring,'</pre>');
                finout=finout(finout>iniout);
                % outstring(finout:finout+11)
                % inclplint = point where output of the example must be included
                inclpoint=finout(1)+18;
                % incl= string which contains the output of the code
                incl=texttoadd{i+numexToExec};
                outstring=[outstring(1:inclpoint) incl outstring(inclpoint+1:end)];
            end
        end
        
        close all
        
        % Remove folder which ahve temporarily added to path
        rmpath([pathstr '\helpfiles\FSDA\tmp'])
        rmpath([pathstr '\utilities\privateFS'])
    end
end

%% WRITE string outstring into the final HTML file
fprintf(file1ID,'%s',outstring);
fclose('all');



end

% This inner function  has the purpose of add symbols </p> <p> every time
% a full stop is followed by a series of space and then a carriage return.
function descrlongHTML=formatHTML(descrlong)
newlinewithFullStop=regexp(descrlong,'[\.\s1-50]\r');
newlinewithColon=regexp(descrlong,'[\:\s1-50]\r');
newl=sort([newlinewithColon newlinewithFullStop]);
if ~isempty(newl)
    descrlongHTML=['<p>' descrlong(1:newl(1))];
    if length(newl)==1
        descrlongHTML=[descrlongHTML '</p> <p>' descrlong(newl(1)+1:end)];
    else
        for j=1:(length(newl)-1)
            descrlongHTML=[descrlongHTML '</p> <p>' descrlong(newl(j)+1:newl(j+1))];
        end
        descrlongHTML=[descrlongHTML descrlong(newl(j+1):end)];
    end
    descrlongHTML=[descrlongHTML '</p>'];
else
    descrlongHTML=descrlong;
end
end

function StringHTML=formatHTMLwithMATHJAX(inputSring)

% Check if symbols \[ \] are present
% If this is the case it is necessary to split inputSring into
% the text_part and the Mathjax_part and apply HTML format just
% to the complementary of the MathJax part
iniMathJax=regexp(inputSring,'\\\[');
finMathJax=regexp(inputSring,'\\\]');

if ~isempty(iniMathJax)
    MoreA=formatHTML(inputSring(1:iniMathJax-1));
    for k=1:length(iniMathJax)
        MoreA=[MoreA inputSring(iniMathJax(k):finMathJax(k)+1)];
        if k==length(iniMathJax)
            MoreA=[MoreA formatHTML(inputSring(finMathJax(k)+2:end))];
        else
            MoreA=[MoreA formatHTML(inputSring(finMathJax(k)+2:iniMathJax(k+1)))];
        end
    end
    StringHTML=MoreA;
else
    % In this case there are not latex formulae so just apply
    % routine formatHTML
    StringHTML=formatHTML(inputSring);
end

end

