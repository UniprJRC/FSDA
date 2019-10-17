function out=mreadFS(file,varargin)
%Enables to create a structure with InputArgs/OptArgs/OutArgs ... from .m function files
%
%<a href="matlab: docsearchFS('mreadFS')">Link to the help function</a>
%
%   mreadFS creates a structure from .m fucntion file. To understand the
%   characteristics the .m file must have, please see the "More About" section
%
% Required input arguments:
%
%    file:         MATLAB File. String. Full or partial
%                  path of the MATLAB file for which a structure containing
%                  all thre relevant help elements has to be created.
%                  Example-'myfile.m'
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
%
% Output:
%
%  The output consists of a structure 'out' containing the following fields:
% out.title     = title of HTML file. String.
%                   String to be included in HTML tag title
% out.purpose   = purpose of the routine. String.
%                 String extracting what is in the second row of the
%                 original input file.
%out.description= short description of the file. String
%                 If this field is not empty it is included in the HTML file
%                 at the beginning of the section "Description"
% out.InpArgs   = Required and Optional input arguments. Cell.
%                 Cell of size k-by-7 containing information about required
%                 and optional input argument.
%                 1st column = name
%                 2nd column = short description
%                 3rd column = type indication (Scalar, matrix, ...)
%                 4th column = string containing long description. If the
%                 input argument is a struct, columns 2 and 3 will be empty and all
%                 information about the fields of the structure will be
%                 included in the 4th column.
%                 5th column = example (if present)
%                 6th column = Data type
%                 7th column = string contaning '1' if the argument is
%                 required and '0' if it is optional
%                 8th column = this column is empty unless 3rd column
%                 contains the word structure. If this is the case the 8th
%                 column will contain a cell with two columns containing
%                 the Values/Description of the associated structure.
% out.OptArgs   = Optional input arguments specified as name/values pairs. Cell.
%                 Cell of size r-by-7 containing information about
%                 and optional input argument specified as name/values pairs.
%                 1st column = name.
%                 2nd column = short description.
%                 3rd column = type indication (Scalar, matrix, strucutre, ...).
%                 4th column = string containing long description.
%                 5th column = example (if present).
%                 6th column = Data type (ex. Single |Double).
%                 7th column = this column is empty unless 3rd column
%                 contains the word structure. If this is the case the 7th
%                 column will contain a cell with two columns containing
%                 the Values/Description of the associated structure.
% out.OutArgs   = Required and Optional (varaargout) output arguments. Cell.
%                 Cell of size k-by-5 containing information about output
%                 and varargout output argument.
%                 1st column = name.
%                 2nd column = short description.
%                 3rd column = type indication (Scalar, matrix, strucutre, ...).
%                 4th column = string containing long description.
%                 5th column =  this column is empty unless output argument
%                 is a structure. If an output argument is a struct,
%                 columns 2, 3 and 4 will be empty and all information
%                 about the fields of the structure will be included in the
%                 5th column.
% out.MoreAbout = More About. String. String containing what in the HTML
%                 file will appear under the section "More About".
%out.Acknowledgements = Acknowledgements. String. String containing what in the HTML
%                 file will appear under the section "Acknowledgements".
%out.References = References. cell. Cell of length r containing the
%                 references.
%   out.SeeAlso = References. cell. Cell of length s containing the
%                 references to linked files.
%        out.Ex = Examples. cell. Cell of length t containing the
%                 examples.
%                 First column= title of the example;
%                 Second column = detailed description;
%                 Third column = MATLAB code;
%                 Fourth column = dummy variable which indicates whether
%                 the example must be executed or not) If 1 example is executed
%    out.ExtraEx= Extra Examples. cell. Cell of length u containing the u
%                 extra examples.
%                 First column= title of the example;
%                 Second column = detailed description;
%                 Third column = MATLAB code;
%                 Fourth column = dummy variable which indicates whether
%                 the example must be executed or not) If 1 example is executed
%   out.laste   = object of class MException which provides information
%                 about last error in executing the examples. If all
%                 examples run without errors laste is an empty value;
%out.InpArgsMisMatch = cell of size k-by-3 which in presence of name/pairs
%                 optional arguments enables to understand which are the
%                 optional arguments which are described but are not used
%                 inside the file and vice versa. More precisely, the first
%                 column contains the list of the options for which there
%                 is a mismatch. The second column contains 1 if the option
%                 was described. The third column contains 1 if the option
%                 was effectively used. Of course, the sum of columns two
%                 and three is always 1.
%                 For example if InpArgsMisMatch is equal to
%                         []    'Options described'    'Options used'
%                 'nomes'       [                0]    [           1]
%                 'refsteps'    [                0]    [           1]
%                 'reftol'      [                0]    [           1]
%                 it means that options 'nomes', 'refstesps' and 'reftol'
%                 have not been described, but are used inside the .m file
%out.OutArgsStructMisMatch = cell of size r-by-3 which in presence of output
%                 arguments which are structures enables to highlight
%                 the fields of the structures which are described but
%                 are not used inside the file and vice versa. More precisely,
%                 the first column contains the list of the fields for
%                 which there is a mismatch. The second column contains 1 if
%                 the field was described. The third column contains 1 is
%                 the field was effectively used. Of course the sum of
%                 columns two and three is always 1. For example if
%                 OutArgsStructMisMatch is
%                   []    'Fields described'    'Fields used'
%                   'out'                    []               []
%                   'Y'      [               1]    [          0]
%                it means that field Y of output structure out has been
%                described but has not been produced by the output of the
%                .m file
%                 It is important to remark that this output is
%                 present only if option evalCode is 1 because publishFS
%                 takes the output of the examples which are run to check
%                 if the output structure contains all fields
%                 which are described. If the output of the function is
%                 called out, then publishFS tries to find if the examples
%                 which were run have produced a structure which contains
%                 the name out. Therefore, if an example produced a
%                 structure named outMM or Myout then publishFS checks the
%                 fields present in outMM or Myout to detect eventual
%                 mismatches.
%
%
%
% More About:
%
%         The .m file (which must be located on the MATLAB path or on the currect
%         folder) must satisfy the following characteristics to be correctly
%         processed.
%
%         [1] The row below the row which starts with 'function ....' must contain the
%         description of the purpose of the .m file.
%         Remark: generally the row which starts with 'function ....' is the first
%         row of an .m file.
%         [2] String 'Required input arguments:' must be present. The lines below
%         this string must contain the description of the compulsory input
%         arguments. Each argument must have the name (a series of spaces from 0
%         to 10) symbol ':' and then the description. The format of the description
%         is as follows:
%         The first sentence after symbol ':' is the title of the input argument
%         and in the HTML file it will appear in bold face in the same line of the
%         input argument (this is the short description of the required input
%         argument).
%         The second sentence after symbol ':' describes the objects
%         (for example, scalar, vector, 3D array) and in the HTML file will appear
%         in the second row.
%         These first two rows will always be visible.
%         What starts with the third sentence after symbol ':' is the detailed
%         description of that particular input argument and in the HTML file it
%         will be visible just if the user clicks on any point in the first two
%         lines or the user clicks on the option expand all.
%         The last line may start with the words "Data Types:" and contains the
%         specification of a particular input argument (e.g. Data Types: single |
%         double). For example, suppose that the .m routine which has to be
%         processed has two required input arguments which are respectively called
%         y and X, then the accepted format is as follows.
%
%                       Required input arguments:
%
%                        y:         Response variable. Vector. Response variable,
%                                   specified as a vector of length n, where n is
%                                   the number of observations. Each entry in y is
%                                   the response for the corresponding row of X.
%                                   Data Types: single | double.
%                      X :          Predictor variables. Matrix of explanatory
%                                   variables (also called 'regressors') of
%                                   dimension n x (p-1) where p denotes the number
%                                   of explanatory variables including the
%                                   intercept. Rows of X represent observations,
%                                   and columns represent variables. By default,
%                                   there is a constant term in the model, unless
%                                   you explicitly remove it using input option
%                                   intercept, so do not include a column of 1s in
%                                   X.
%                                   Data Types: single | double.
%
%         IMPORTANT NOTICE: if an input argument is a structure (publishFS
%         automatically checks if the input argument contains the word "structure"
%         then the fields of the structure will be automatically included into a
%         HTML table). In this case the fields of the structure are identified as
%         the lines which contain "a series of spaces" "a_word" "a series
%         of spaces followed by symbol '='". For example suppose the an input
%         option is called bayes and object bayes is a structure with field names
%         beta0, R, tau0 and n0, the accepted format is as follows.
%
%            bayes      : prior information. Structure.
%                               It contains the following fields.
%                   out.beta0=  p-times-1 vector containing prior mean of \beta.
%                   out.R    =  p-times-p positive definite matrix which can be
%                               interpreted as X0'X0 where X0 is a n0 x p matrix
%                               coming from previous experiments (assuming that the
%                               intercept is included in the model.
%                     out.tau0 = scalar. Prior estimate of tau.
%                     out.n0   = scalar. Sometimes it helps to think of the prior
%                              information as coming from n0 previous experiments.
%
%
%
%         [3] If the input .m file between the row which starts with
%          <a href="matlab: docsearchFS(.....
%          and the row with the string
%          "Required input arguments:"
%          contains a series of sentences, they will be automatically included in
%          the HTML output just below the description.
%
%          An example is given below:
%                 'function [out , varargout]  = tclust(Y,k,alpha,restrfactor,varargin);
%                 tclust computes trimmed clustering with restricitons on
%                 the eigenvalues.
%
%                 <a href="matlab: docsearchFS('tclust')">Link to the help function</a>;
%
%                   tclust partitions the points in the n-by-v data matrix Y into k
%                   clusters.  This partition minimizes the trimmed sum, over all
%                   clusters, of the within-cluster sums of
%                   point-to-cluster-centroid distances....
%
%                  Required input arguments:'.
%
%         [4] String 'Optional input arguments:' must be present even if there are
%         no optional arguments.
%         publishFS, in order to understand what are the names of the optional
%         input arguments scans the rows below the string "Optional input
%         arguments:" and identifies the lines which contain the optional arguments
%         as those which contain "a series of spaces" "a_word" "a series of spaces
%         followed by symbol ':'".
%         The first sentence after symbol ':' is the title of the optional input
%         argument and in the HTML file it will appear in the same row of the name
%         of the optional input argument (this is the short description of the
%         optional input argument).
%         The second sentence after symbol ':' describes the objects (for example,
%         scalar, vector, 3D array) and in the HTML file will appear in the second
%         row. These first two rows will always be visible.
%         What starts with the third sentence after symbol ':' is the detailed
%         description of that particular optional input argument and in the HTML
%         file it will be visible just if the user clicks on any point in the first
%         two lines or the user clicks on the option expand all.
%         The last two lines of each optional input argument MUST start with the
%         words 'Example -' and 'Data Types -' followed by a string without spaces
%         which specifies a possible use of the option and the type of data. For
%         example, suppose that the first two optional arguments are called
%         respecively 'intercept' and 'h', then the
%         accepted format is as follows:
%
%
%                       Optional input arguments:
%
%                       intercept :  Indicator for constant term. Scalar.
%                                   If 1, a model with constant term will be fitted
%                                   (default),
%                                   else no constant term will be included.
%                                   Example - 'intercept',1.
%                                   Data Types - double.
%                               h : The number of observations that have determined the least
%                                     trimmed squares estimator. Scalar.
%                                     Example - 'h',round(n*0,75).
%                                     Data Types - double.
%
%
%
%
%          IMPORTANT NOTICE: given that options are identified as those which have
%          symbol "%" followed by "a series of spaces" then "a word" then "a series
%          of spaces" then symbol ":", each line inside the description does not
%          have to start as follows "%   ANYWORD   :" because the parser will
%          wrongly identify "ANYWORD" as an additional optional input argument. The
%          only once exception to this rule is the word "%  REMARK :". However, if
%          there is a remark, it must be put at the very end of the description of
%          the optional input argument. At the very end means after the rows
%           Example and Data Types.
%
%         [5] String 'Output:' must be present.
%         The lines after string 'Output:'
%         must contain the list of the output arguments using the same rules
%         described above for the optional arguments. In this case, however, the
%         identification of the ouptut arguments is easier because they are
%         extracted directly from the first line of the file (e.g. if the first
%         line is as follows
%         function [mdr,Un,BB,Bols,S2] = FSRmdr(y,X,bsb,varargin) then the 5
%         output arguments are immediately known to the parser).
%         In the case of output argument publishFS checks if the first 50
%         characters contain the words "which contains" or "containing" e.g.:
%
%                      mdr:         n -init x 2
%                                   matrix containing the
%                                   monitoring of minimum deletion residual.
%                                   1st col = fwd search
%                                   ........
%                        Un:        (n-init) x 11 Matrix
%                                   which contains the unit(s) included in the
%                                   subset at each step of the search.
%                                   ...........
%
%         In this case what is after the strings "which contains" or "containing"
%         will appear in bold face as the title of the output argument. What is
%         before the strings "which contains" or "containing" will appear in the
%         second row.
%
%         For example, the above lines will be processed as follows:
%
%              mdr   -  Monitoring of minimum deletion residual;
%              n -init x 2 matrix.
%
%              Un    - unit(s) included in the subset at each step of the search.
%              (n-init) x 11 Matrix which contains the unit(s) included in the
%              subset at each step of the search.
%
%         If in the HTML file the user clicks on them the expdanded description
%         (that is what starts after the second full stop will appear).
%
%         Alternatively, if the first 50 characters of each output argument do not
%         contain the strings "which contains" or "containing" the following
%         convention is used. The first sentence after symbol ":" is assumed
%         to be the title of the output argument and in the HTML file it will
%         appear in bold face in the same line of the name of output
%         argument. The second sentence (words between first and second full stop)
%         will appear in the second row. The third sentence is the full description
%         of the output argument. For example, suppose that the output of a
%         procedure contains the objects mdr and Un, the accepted format
%         is as follows.
%               Output:
%
%                      mdr:         Minimum deletion residual. Matrix.  n -init x 2
%                                   matrix which contains the
%                                   monitoring of minimum deletion residual.
%                                   1st col = ...
%                      Un:          Units included. Matrix. (n-init) x 11 Matrix
%                                   which contains the unit(s) included in the
%                                   subset at each step of the search.
%                                   REMARK: in every step ....
%
%         The above lines will be processed as follows:
%
%              mdr   -  Minimum deletion residual. Matrix.
%                       n -init x 2 matrix which contains the
%                       monitoring of minimum deletion residual.
%                       1st col = ...
%
%              Un    - Units included. Matrix.
%                       (n-init) x 11 Matrix which contains the unit(s)
%                       included in the subset at each step of the search.
%                       REMARK: in every step ....
%
%         If in the HTML file the user clicks on them the expdanded description
%         (that is what starts after the second full stop will appear).
%
%         IMPORTANT NOTICE: Similarly to what happend for each input argument, if
%         an output argument is a structure, publishFS automatically checks if it
%         contains the words "structure" and "field". In this case, the fields of
%         the structure will be automatically included into a HTML table. The
%         fields of the structure are identified as the lines which contain "a
%         series of spaces" "name_of_output_structure.a_word" "a series of spaces
%         followed by symbol '='". For example suppose that the output of a
%         procedure is an object called out which is a structure with two fields
%         out.rew and out.beta, an accepted format is as follows:
%
%                                  %  Output:
%
%                                  out :     A structure containing the following fields:
%                                            out.rew  = Scalar if out.rew=1 all
%                                                       subsequent output refers to
%                                                       reweighted else no
%                                                       reweighting is done.
%                                            out.beta = Vector of beta LTS (LMS)
%                                                       coefficient estimates,
%                                                       including the intercept
%                                                       when options.intercept=1.
%                                                       out.beta=[intercept
%                                                       slopes].
%
%
%         PLEASE REMEMBER THAT THE FIELDS of an output instance HAVE TO CONTAIN THE
%         "=" SIGN AND NOT THE ":" SIGN.
%
%         REMARK: If there is the string REMARK after the description of the last
%         field of the structure, all the words after REMARK are put outside and
%         below the HTML table.
%
%         If the description of a particular output has the string "which contains"
%         or "containing",  as follows
%
%                      mdr:          n -init x 2 matrix which contains the
%                                   monitoring of minimum deletion residual at each
%                                   step of the forward search.
%                                   1st col = fwd search index (from init to n-1).
%                                   2nd col = minimum deletion residual.
%
%         publishFS will try to put what comes before the string "which
%         contains" or "containing" inside the subtitle (second row) of the each
%         ouptut argument in the HTML file. For example, the example above in the
%         HTML file will be processed as follows:
%                        mdr -Monitoring of minimum deletion residual at each step of the forward search.
%                        n -init -by- 2 matrix.
%         If, in the HTML file, the user clicks on the first line:
%                       "mdr -Monitoring...";
%         the expanded description will automatically appear.
%
%         [6] A line which starts with string "See also:" must be present.
%         Linked m files must be separated by symbol ",". For example, suppose that
%         files FSRBmdr.m and FSR.m have connections with the current file, then an
%         accepted format is         See also: FSRBmdr.m, FSR.m.
%
%         [7] A line which starts with string "References:" must be present.
%         The year of each reference must be enclosed in round parenthesis.
%         PublishFS decides that a new reference starts if a new line contains
%         symbol "(" + "a sequence of 4 or 5 characthers identifying the year
%         because the reference can be for example 2003 or 2003a" + symbol ")"
%         For example, an acceptable format for the two references below is:
%
%
%                         Chaloner and Brant (1988). A Bayesian Approach to Outlier
%                         Detection and Residual Analysis, Biometrika, Vol 75
%                         pp. 651-659.
%                         Riani M., Corbellini A., Atkinson A.C. (2015), Very
%                         Robust Bayesian Regression for Fraud Detection,
%                         submitted.
%
%         [8] All the examples associated with the file which has to be processed
%         must be enclosed inside Percent-braces (comments blocks, i.e.
%         symbols "%{" and "%}" ).
%         The first sentence identifies the title of the comment which
%         will appear in the HTML file.
%         IMPORTANT NOTICE: if the comment has to be executed, the first line
%         associated with the title must start with two "%%" symbols instead of just
%         one "%" symbol. The examples in the first positions will appear in
%         the HTML file under the caption "Examples" while the latest will appear
%         under the caption "Related Examples". More precisely, if the output of a
%         procedure contains k outputs and some optional input arguments the first
%         k+1 comment blocks will appear in the HTML file under "Examples".
%         First comment block is associated with the call of the procedure with
%         just one output and all default input arguments.
%         Second comment block is associated with the call of the procedure with
%         just one output and with some optional input arguments.
%         Third comment block is associated with the call of the procedure with
%         two output arguments.
%         ...
%         k+1 comment block is associated with the call of the procedure with
%         k output arguments.
%         k+2 comment block is the first which in the HTML file will appear under
%         the heading "Related Examples".
%         For example, suppose that the first example of procedure FSRmdr has to be
%         executed and its output must be included into the HTML file, then the accepted
%         format is as follows ("please notice the two symbols "%%" in the
%         first row").
%
%
%                         "%{".
%                             "%% FSRmdr with all default options".
%                             "% Compute minimum deletion residual".
%                             "% Monitor minimum deletion residual in each step of the forward search".
%                             "% Common part to all examples: load fishery dataset".
%                              load('fishery');
%                              y=fishery.data(:,2);
%                              X=fishery.data(:,1);
%                              "% Find starting subset"
%                              [out]=LXS(y,X,'nsamp',10000);
%                              [mdr] = FSRmdr(y,X,out.bs);
%                              plot(mdr(:,1),mdr(:,2))
%                              title('Monitoring of minimum deletion residual')
%                         "%}".
%
%                         "%{".
%                             "% FSRmdr with optional arguments".
%                             "% Choose step to start monitoring".
%                             "% Compute minimum deletion residual and
%                              start monitoring it from step 60".
%                             [mdr] = FSRmdr(y,X,out.bs,'init',60).
%                         "%}".
%
%          In order to understand where the example finish and the MATLAB code
%          starts publishFS checks if one of the following strings
%         is present
%         "Input parameters checking",
%         "Beginning of code",
%         "nargin".
%         If this is the case the search of "comments blocks signs"
%         (i.e. symbols  "%{" and "%}") is limited to the first character prior
%         to the detection of one of the previous strings. This
%         modification has been added in order to make sure that there are
%         no additional block signs within matlab code.
%
%         [9] If a procedure contains varargout, then the string:
%                       "Optional Output:"
%         must be present. For example suppose there is a function called mcd which
%         has the following sintax:
%
%                         function [RAW,REW,varargout] = mcd(Y,varargin).
%
%         Then at the end of the output argument the format must be as follows:
%
%                               Optional Output:
%
%                                         C     : matrix of the indices of the
%                                                 subsamples extracted for
%                                                 computing the estimate.
%
%
%--------------------------------------------------------------------------
%
%         GENERAL REMARKS:
%
%         -----------------------------------------------------------------------:.
%         REMARK1: if symbol "%" is wanted and it is not a simple comment delimiter, it
%         must be replaced by words "per cent". For example, string "50% envelope"
%         must become "50 per cent" envelope.
%         -----------------------------------------------------------------------:.
%         REMARK2: If there is just one output argument it can be without square
%         brackets. Among the input elements of a procedure the number of spaces
%         between them is not important. For example
%         "y,X,varargin" or "y, X   ,  varargin"   are both fine.
%         -----------------------------------------------------------------------:.
%         REMARK 4: publishFS uses MathJax javascript in order to interpret the
%         equations inside the .m file written in Latex style. For in line
%         equations both symbols dollar dollar  and backslash(  backslash) are accepted.
%         For single line equations symbols backslash[ backslash] must be used.
%         -----------------------------------------------------------------------:.
%         REMARK 5: if there are not enough examples in the .m file the procedure
%         still runs but a warning will be automatically produced.
%         -----------------------------------------------------------------------:.
%
% See also: publishFS.m
%
%
% References:
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mreadFS')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
  % Create output structure out starting from file FSRmdr.
  out=mreadFS('FSRmdr')
%}

%{
  % Create output structure out starting from file FSRmdr and
  % display detailed information about the Input, Output and Optional
  % arguments.
  out=publishFS('FSRmdr','Display','iter-detailed')
%}



%% Beginning of code

if ~ischar(file)
    error('FSDA:mreadFS:WrongInputOpt','input must be a string containing the name of the file.');
end

evalCode=true;
Display='none';


if nargin>1
    options=struct('Display',Display);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:regressB:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
        
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
        
    end
    
    Display=options.Display;
    
end

%% Open input .m file, put it in a string and do a series of preliminary operations
FileWithFullPath=which(file);
[pathstrcf,name,ext]=fileparts(FileWithFullPath);

if isempty(pathstrcf)
    error('FSDA:publishFS:WrongFile','SourceNotFound');
end

if ~strcmp('.m',ext)
    error('FSDA:publishFS:WrongFileExt','Input file must have m extension')
end


filename=FileWithFullPath;
% f = fopen(filename);

fileID = fopen(char(filename), 'r+');


% Insert the file into fstring
fstring=fscanf(fileID,'%c');

if strcmp(Display,'iter-detailed')
    disp('Check that file contains appropriate reference to itself inside docsearchFS')
end
%linkHTML must a vector with length equal to 2
linkHTML=regexp(fstring,['docsearchFS\(''' name '''\)']);

% Replace < and > symbols with HTML code
%fstring=regexprep(fstring,'[^%]<','&lt;');
fstring=regexprep(fstring,'<','&lt;');
fstring=regexprep(fstring,'>','&gt;');
% replace if present symbol ü with its HTML code
fstring=regexprep(fstring,'ü','&uuml;');



%-----------------
%% CREATE Name-Value Pair Arguments SECTION OF HTML FILE
inselOpt=regexp(fstring,'Optional input arguments:');
if isempty(inselOpt)
    disp('Please check source .m input file')
    error('FSDA:WrgInpFile',['Input file ''' file '.m'' does not contain ''Optional input arguments:'' string'])
end

% substring to search start from Optional input arguments:
fstringselOpt=fstring(inselOpt(1):end);

endpoint=regexp(fstringselOpt,'Output:');
if isempty(endpoint)
    disp('Please check .m input file')
    error('FSDA:missOuts','Input .m file does not contain ''Output:'' string')
end
fstringselOpt=fstringselOpt(1:endpoint-2);

% Find any string which
% begins with % sign then
% contains a series of white space which can go from 0 to 15 then
% contains any single word
% a series of white spaces which can go from 0 to 10 then
% character : then
% a seris of white spaces
% The inipoint of the following two regular expressions is the same but
% however we want to exclude the lines where symbol : is followed by a
% series of white spaces and then by a carriage return because these are
% input arguments but simply are the beginning of a list
[iniA,finA]=regexp(fstringselOpt,'%\s*\w*\s*:'); % '%\s{0,15}\w*\s{0,10}:'
[inichk,~]=regexp(fstringselOpt,'%\s*\w*\s*:\s*\w'); %'%\s{0,15}\w*\s{0,10}:\s{0,10}\w'
% Select all rows of iniA and finA in which the elements of iniA are equal
% to those of inichk
if length(iniA)>length(inichk)
    [~,ia]=intersect(iniA,inichk);
    ini=iniA(ia);
    fin=finA(ia);
else
    ini=iniA;
    fin=finA;
end


% listOptArgs = list which contains all optional arguments (7 columns)
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
% string  Data Types - (i.e.: char, double....)
% If the third column contains the string struct then the content of
% value/description of the struct input argument will appear in the seventh
% column else the seventh column will be empty
listOptArgs=cell(length(ini),7);

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
        % OLD
        % posPercentageSigns=regexp(descradd,'%');
        % descradd(posPercentageSigns)=[];
        
        posPercentageSigns=regexp(descradd,'\n%');
        descradd(posPercentageSigns+1)=[];
        % Remove First character % if remained
        if ~isempty(descradd) && strcmp(descradd(1),'%')
            descradd=descradd(2:end);
        end
        
        % descradd=strtrim(descradd);
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
        % OLD
        % posPercentageSigns=regexp(descrtosplit,'%');
        % descrtosplit(posPercentageSigns)=[];
        
        
        posPercentageSigns=regexp(descrtosplit,'\n%');
        descrtosplit(posPercentageSigns+1)=[];
        
        [inifullstops]=regexp(descrtosplit,'\.');
        if isempty(inifullstops)
            error('FSDA:publishFS:WrongInp',['Sentence''' descrtosplit '''must contain at least two full stops'])
            % error('Wrong input')
        end
        descrtitle=strtrim(descrtosplit(1:inifullstops(1)));
        listOptArgs{ij,2}=descrtitle;
        
        try
            descrtype=strtrim(descrtosplit(inifullstops(1)+1:inifullstops(2)));
        catch
            %errmsg=[' Options found by the parser are:\n'  listOptArgs{:,1}];
            %disp(errmsg)
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
        
        CheckExample=regexp(descrlong,'Example -');
        
        if ~isempty(CheckExample)
            % Just in case take only the very last example
            CheckExample=CheckExample(end);
            Datatypes=regexp(descrlong,'Data Types -','once');
            
            
            % The first word of example code must be embedded around tags <code> </code>
            examplecode=descrlong(CheckExample+9:Datatypes-1);
            % Store string containing the examples
            listOptArgs{ij,5}=examplecode;
            listOptArgs{ij,6}=descrlong(Datatypes+13:end);
            descrlong=descrlong(1:CheckExample-1);
        end
        
        if ~isempty(strfind(listOptArgs{ij,3},'tructure'))
            [listStructureArgs]=formatHTMLstructure(descrlong,listOptArgs{ij,1});
            listOptArgs{ij,7}=listStructureArgs;
        else
            listOptArgs{ij,4}=descrlong;
        end
        
        
        ij=ij+1;
    end
    
end
listOptArgs=listOptArgs(1:ij-1,:);

% Check if the last row, column 4 of list listOptArgs contains the sentence
%' Remark:      The user should only give .....'
% Given that this sentence is very generic and not applied to the last
% optional input argument, if it is present it is deleted
if ~isempty(listOptArgs)
    Checklastremark=listOptArgs{end,4};
    DelTheUser=regexp(Checklastremark,'Remark\s*:\s*The user','once','start','ignorecase');
    if ~isempty(DelTheUser)
        
        listOptArgsIns=Checklastremark(1:DelTheUser-1);
        % Remove extra % signs if present in the string
        posPercentageSigns=regexp(listOptArgsIns,'\s%\s*');
        listOptArgsIns(posPercentageSigns+1)=[];
        listOptArgs{end,4}=listOptArgsIns;
    end
end
if strcmp(Display,'iter-detailed')
    disp('Detailed information about Optional arguments')
    disp(listOptArgs)
end



[startIndex] = regexp(fstring,'%');
% startIndex(2)-3 because there is also the carriage return
purpose=fstring(startIndex(1)+1:startIndex(2)-3);




[gendescini]=regexp(fstring,'/a','once');
[gendescfin]=regexp(fstring,'Required input arguments','once');
gendesc=fstring(gendescini+6:gendescfin-1);
posPercentageSigns=regexp(gendesc,'%');
gendesc(posPercentageSigns)=[];
% Remove from string descri leading and trailing white spaces
gendesc=strtrim(gendesc);


%% Create Sintax section of HTML file (referred to the loop of input arguments)
% Find point where first input argument starts

% Now find the number of required input
% arguments and if there are varargin
% and the number of required output
% arguments
% Find the string containing output arguments
%outargs= [out1, out2, out3]
[startIndexEq] = regexp(fstring,'=');

if startIndexEq(1)<regexp(fstring,'%','once')
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
else
    outargs='';
    nargout=0;
end


% Required input arguments
% Find number of compulasory input arguments
% nargin = number of commas in string  functionname(inp1,inp2,inp3, ...)
[startIndexInp] = regexp(fstring,'(');
[endIndexInp] = regexp(fstring,')');
% inputargs =  (inp1,inp2,inp3, ...)
InputArgs=fstring(startIndexInp(1):endIndexInp(1));
% Check if inputargs contains the string varargin but not the string
% vvargin
[OptArgsVarargin]=regexp(InputArgs,'varargin');

[OptArgsvvarargin]=regexp(InputArgs,'vvarargin');
if ~isempty(OptArgsvvarargin)
    OptArgsVarargin=[];
end

[commasIn] = regexp(InputArgs,',');
j=1;

% nTOTargin= total number of input arguments (requested + optional),
% excluding name-value pairs arguments

if isempty(OptArgsVarargin)
    % Write all required input arguments in cell listargins
    if endIndexInp(1)-startIndexInp(1)>1
        nTOTargin=length(commasIn)+1;
    else
        nTOTargin=0;
    end
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
% which are optionals. That is, check if the intersection between the first
% columns of cell listOptArgs and vector listargins is empty
[OptArgsWithoutNameValue,PosOpt]=intersect(listargins(:,1),listOptArgs(:,1));
nOPTargin=length(OptArgsWithoutNameValue);
% nREQargin = number of required input arguments
nREQargin=nTOTargin-nOPTargin;
PosOpt=sort(PosOpt);

sintax=cell(nargout+1+nOPTargin,1);

if nOPTargin>0
    for j=1:nOPTargin
        if ~isempty(outargs)
            sintax{j}=[outargs(2:commasOut(1)-1) '=' name InputArgs(1:commasIn(PosOpt(j)-1)-1) ')'];
        else
            if ~isempty(commasIn)
                sintax{j}=[name InputArgs(1:commasIn(PosOpt(j)-1)-1) ')'];
            else
                sintax{j}=name;
            end
        end
        
    end
    j=j+1;
else
    j=1;
end


if isempty(OptArgsVarargin)
    if ~isempty(outargs)
        sintax{j}=[outargs(2:commasOut(1)-1) '=' name InputArgs];
    else
        sintax{j}=[name InputArgs];
    end
    j=j+1;
else
    % In this case there is also Name Value line
    strinputarg=strtrim(InputArgs(1:OptArgsVarargin-1));
    if strcmp(strinputarg(end),',')
        strinputarg=strinputarg(1:end-1);
    end
    
    if ~isempty(outargs)
        sintax{j}=[outargs(2:commasOut(1)-1) '=' name strinputarg ')'];
        if length(strinputarg)>1
            sintax{j+1}=[outargs(2:commasOut(1)-1) '=' name strinputarg ',Name,Value)'];
        else
            % just in case function has no compulsary input argument then
            % the comma before 'Name.value' is unnecessary
            sintax{j+1}=[outargs(2:commasOut(1)-1) '=' name strinputarg 'Name,Value)'];
        end
        
    else
        sintax{j}=[name strinputarg ')'];
        sintax{j+1}=[name strinputarg ',Name,Value)'];
    end
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



%% CREATE DESCRIPTION SECTION OF HTML FILE

if strcmp(Display,'iter-detailed')
    disp('Examples')
    disp(sintax)
    disp('---------------')
end

% start of example j
[startIndexEx] = regexp(fstring,'\n%\{[\s1-20]');
[endIndexEx] = regexp(fstring,'\n%\}[\s1-20]');

% Try to understand where preamble finishes and MATLAB code starts
EndOfExtry1=regexp(fstring,'[^"]Input parameters checking');
EndOfExtry2=regexp(fstring,'[^"]Beginning of code');
EndOfExtry3=regexp(fstring,'[^"]nargin');
EndOfExtry=min([EndOfExtry3 EndOfExtry1 EndOfExtry2]);
if ~isempty(EndOfExtry)
    startIndexEx=startIndexEx(startIndexEx<EndOfExtry);
    endIndexEx=endIndexEx(endIndexEx<EndOfExtry);
end

%listEx = list which contains the examples (associated to sintax)
% First column= title, second column detailed description. Third column
% code
% Fourth column is a dummy variable which indicates whether the example must be
% executed or not). If 1 example is executed
listEx=cell(length(sintax),4);

for j=1:length(sintax)
    
    
    %For each element of input and output argument a hypertext link is added
    sintaxj=sintax{j};
    % Locate in expression [out1,out2,...]=namefunc(inp1,inp2,...) the
    % position of equal sign
    [startIndex] = regexp(sintaxj,'=');
    
    outs=sintaxj(1:startIndex-1);
    
    if ~isempty(outs)
        commaspos=regexp(outs,',');
        if isempty(commaspos)
            noutel=1;
        else
            noutel=length(commaspos)+1;
        end
        if j==length(sintax)
            % Write in cell listargouts the list of output arguments
            listargouts=cell(noutel,1);
        end
        if noutel>1
            for i=1:noutel
                if i==1
                    outi=['[' outs(2:commaspos(i))];
                    if j==length(sintax)
                        listargouts{i}=strtrim(outi(2:end-1));
                    end
                elseif i==noutel
                    outi=outs(commaspos(i-1)+1:end);
                    if j==length(sintax)
                        listargouts{i}=strtrim(outi(1:end-1));
                    end
                else
                    outi=outs(commaspos(i-1)+1:commaspos(i));
                    if j==length(sintax)
                        listargouts{i}=strtrim(outi(1:end-1));
                    end
                end
                
            end
        else
            outi=strtrim(outs);
            if j==length(sintax)
                listargouts=cell(1,1);
                listargouts{1}=outi;
            end
        end
    else
        listargouts='';
    end
    
    
    try
        stri=fstring(startIndexEx(j)+3:endIndexEx(j)-1);
    catch
        %disp(stri)
        warning('FSDA:wrongEx','This file does not contain enough examples, please add them!')
        stri='EXAMPLES TO ADD';
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
    if isempty(inicr) && strcmp(stri,'EXAMPLES TO ADD')~=1
        disp('String below seems to be without carriage return')
        disp('------------------------------------------------')
        disp(stri)
        errmsg=['Carriage return could not be found in the example section \n'...
            'Probably file has been created using Linux'];
        error('FSDA:wrongdelimiter',errmsg)
    end
    
    % This is the first line which does not contain symbol %
    for jj=1:length(inicr)-1
        strtest=stri(inicr(jj):inicr(jj+1));
        if isempty(regexp(strtest,'%','once'))
            break
        end
    end
    findescriptionEx=inicr(jj);
    strdescrEx=stri(endtitle+1:findescriptionEx);
    
    posPercentageSigns=regexp(strdescrEx,'\D%');
    strdescrEx(posPercentageSigns+1)=[];
    
    listEx{j,2}=strdescrEx;
    % listEx{j,3}=stri(findescriptionEx+1:end);
    
    StringWithLTandGT=stri(findescriptionEx+1:end);
    StringWithoutLTandGT=strrep(StringWithLTandGT,'<','&lt;');
    StringWithoutLTandGT=strrep(StringWithoutLTandGT,'>','&gt;');
    listEx{j,3}=StringWithoutLTandGT;
    
    %---------
    
end

% Now check whether the first column of cell listEx contains the string interactive_example
NumOfInterEx=1;
for i=1:size(listEx,1)
    [StartInteractive,EndInteractive]=regexp(listEx{i,1},'[Ii]nteractive_example.');
    if ~isempty(StartInteractive)
        StringToReplace=listEx{i,1};
        listEx{i,1}=['<i>Interactive example ' num2str(NumOfInterEx)  '.</i>' StringToReplace(EndInteractive+1:end)];
        NumOfInterEx=NumOfInterEx+1;
    end
end

if strcmp(Display,'iter-detailed')
    disp('Detailed information about all the examples')
    disp(listEx)
end


%% CREATE EXAMPLES SECTION OF HTML FILE



%% CREATE RELATED EXAMPLES SECTION OF HTML FILE
if length(startIndexEx)>length(sintax)
    lsintax=length(sintax);
    % Fourth column of listextraEX contains flag 1 or 0 depending on the
    % fact that the example must be executed or not
    listExtraEx=cell(length(startIndexEx)-lsintax,4);
    
    
    for j=1:size(listExtraEx,1)
        stri=fstring(startIndexEx(j+lsintax)+3:endIndexEx(j+lsintax)-1);
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
        for jj=1:length(inicr)-1
            strtest=stri(inicr(jj):inicr(jj+1));
            
            % break when you find the first line which does not contain symbol %
            %             if isempty(regexp(strtest,'%','once'));
            %                 break
            %             end
            
            % NEW code: break when you find the first line which does not start with symbol %
            strtest=strtrim(strtest);
            if  ~isempty(strtest) && ~strcmp(strtest(1),'%')
                break
            end
        end
        findescriptionEx=inicr(jj);
        strdescrEx=stri(endtitle+1:findescriptionEx);
        % remove % signs from strdescrEx
        % OLD
        %posPercentageSigns=regexp(strdescrEx,'%');
        %strdescrEx(posPercentageSigns)=[];
        
        posPercentageSigns=regexp(strdescrEx,'\D%');
        strdescrEx(posPercentageSigns+1)=[];
        
        
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
    
    % Now check whether the first column of cell listExtraEx contains the string interactive_example
    for i=1:size(listExtraEx,1)
        [StartInteractive,EndInteractive]=regexp(listExtraEx{i,1},'[Ii]nteractive_example.');
        if ~isempty(StartInteractive)
            StringToReplace=listExtraEx{i,1};
            listExtraEx{i,1}=['<i>Interactive example ' num2str(NumOfInterEx)  '.</i>' StringToReplace(EndInteractive+1:end)];
            NumOfInterEx=NumOfInterEx+1;
        end
    end
else
    listExtraEx='';
end





%% Create listInptArgs and related HTML code
% listInpArgs = list which contains all input arguments (required or not)
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
% The seventh col will contain '1' or '0' according to the fact that the
% associated argument is required or optional
% The eigth col will contain a cell with the value/desctiption
% arguments of the structure

listInpArgs=cell(nTOTargin,8);
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
    
    % OLD WAY OF FINDING inipoint
    % inipoint=regexp(fstringsel,['%\s{0,10}' listargins{i} '\s{0,10}:']);
    inipoint=regexp(fstringsel,['%\s*' listargins{i} '\s*:']);
    if isempty(inipoint)
        error('FSDA:missInpDescr',['No description has been found for compulsory input argument '  listargins{i} ])
        
    else
    end
    
    % The endpoint of the substring is "See also" or the next optional input argument
    if i <nREQargin
        endpoint=regexp(fstringsel,['%\s*' listargins{i+1} '\s*:']);
        if isempty(endpoint)
            disp('Please check .m input file')
            error('FSDA:missInps',['Input .m file does not contain the description' ...
                ' for input argument '''  listargins{i+1} ''''])
        end
        
    elseif i==nREQargin
        endpoint=regexp(fstringsel,'Optional input arguments:');
        if isempty(endpoint)
            disp('Please check .m input file')
            error('FSDA:missOuts','.m file does not contain ''Optional input arguments:'' string')
        end
    elseif i<nTOTargin
        endpoint=regexp(fstringsel,['%\s*' listargins{i+1} '\s*:']);
        if isempty(endpoint)
            disp('Please check .m input file')
            error('FSDA:missInps',['Input .m file does not contain the description' ...
                ' for input argument '''  listargins{i+1} ''''])
        end
        
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
    DescrInputToSplit=[DescrInputToSplit ' ']; %#ok<AGROW>
    
    [inifullstops]=regexp(DescrInputToSplit,'\.[\s1-3]');
    if isempty(inifullstops)
        error('FSDA:publishFS:WrongInp',['Input option: ''' inpi '''\n Sentence\n''' DescrInputToSplit '''\nmust contain at least two full stops'])
        % error('Wrong input')
    end
    shortdesc=strtrim(DescrInputToSplit(1:inifullstops(1)));
    % Store title of the i-th input argument
    %     % remove sign : if present at the beginning of the sentence
    if strcmp(shortdesc(1),':')
        shortdesc=strtrim(shortdesc(2:end));
    end
    listInpArgs{i,2}= shortdesc;
    
    try
        descrtype=strtrim(DescrInputToSplit(inifullstops(1)+1:inifullstops(2)));
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
    
    % Check if the input is a structure with fields. In this case it is
    % necessary to create a table
    
    if ~isempty(strfind(listInpArgs{i,3},'tructure')) && ~isempty(strfind(descrlong,'field'))
        
        Datatypes=regexp(descrlong,'Data Types -','once');
        if ~isempty(Datatypes)
            listInpArgs{i,6}=descrlong(Datatypes+13:end);
            
            descrlong=descrlong(1:Datatypes-1);
        else
            warning('FSDA:publishFS:MissingDataType',['Input argument ''' inpi ''' does not contain DataType line, by default string  ''single| double'' has been added'])
            listInpArgs{i,6}='single| double';
        end
        
        Examplesini=regexp(descrlong,'Example -','once');
        if ~isempty(Examplesini)
            listInpArgs{i,5}=descrlong(Examplesini:end);
            descrlong=descrlong(1:Examplesini-1);
        else
            listInpArgs{i,5}='';
            %TODO5
        end
        % Extract what comes before the description of the fields of the
        % structure
        [iniA]=regexp(descrlong,[ inpi '\.\w*\s*=']);
        if ~isempty(iniA)
            
            descrlongini=descrlong(1:iniA(1)-1);
            listInpArgs{i,4}=descrlongini;
            descrlong=descrlong(iniA(1):end);
        end
        % Include the fields of the structure in the 8th column
        [listStructureArgs]=formatHTMLstructure(descrlong,inpi);
        listInpArgs{i,8}=listStructureArgs;
        
    else
        
        
        if i<=nREQargin
            
            Datatypes=regexp(descrlong,'Data Types -','once');
            if ~isempty(Datatypes)
                listInpArgs{i,6}=descrlong(Datatypes+13:end);
                
                descrlong=descrlong(1:Datatypes-1);
                
                
                listInpArgs{i,4}=descrlong;
            else
                
                listInpArgs{i,4}=descrlong;
                if strcmp(Display,'iter-detailed')
                    warning('FSDA:publishFS:MissingDataType',['Input argument ''' inpi ''' does not contain DataType line, by default string  ''single| double'' has been added'])
                end
                
                listInpArgs{i,6}='single| double';
            end
        else
            
            % Check if descrlong contains
            % Example - and Data types -
            
            CheckExample=regexp(descrlong,'Example -','once');
            if ~isempty(CheckExample)
                Datatypes=regexp(descrlong,'Data Types -','once');
                descrlonginp=descrlong(1:CheckExample-1);
                listInpArgs{i,4}=descrlonginp;
                
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
            
        end
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
    listOptArgs='';
end

%% CREATE OUTPUT ARGUMENTS SECTION OF HTML FILE



% check if the last element of listargouts is varargout
% listargouts
%  Optional Output:
outsel=regexp(fstring,'Output:');
if isempty(outsel)
    disp('Please check HTML input file')
    error('FSDA:missOuts','HTML file does not contain ''Output:'' string')
end
% substring to search starting from Output:
fstringsel=fstring(outsel(1):end);

% cell which will contain the details of output arguments
listOutArgs=cell(length(listargouts),5);
listOutArgs(:,1)=listargouts;

if nargout>0
    for i=1:nargout
        
        % listargouts is a cell which contains the list of output arguments
        outi=listargouts{i};
        
        
        % The initial point of the string is 'listargouts{i}' is there is just
        % one output else is string 'listargouts{i} :' is there is more than
        % one output and this is not varargout
        % else if there is varargout the initialpoint is the string
        % "Optional Output:"
        if length(listargouts)==1 && strcmp(listargouts,'varargout') ==0
            inipoint=regexp(fstringsel,listargouts{i});
        elseif  i<length(listargouts)
            inipoint=regexp(fstringsel,[listargouts{i} '\s{0,11}:']);
        else
            if strcmp(listargouts{end},'varargout') ==0
                inipoint=regexp(fstringsel,[listargouts{i} '\s{0,11}:']);
            else
                inipoint=regexp(fstringsel,'Optional Output:')+8;
            end
        end
        
        if isempty(inipoint)
            error('FSDA:missingOuts',['Output argument ' listargouts{i} ' has not been found'])
        end
        
        % The endpoint of the substring is 'more About'. or See also or the next output argument
        if i <nargout-1
            endpoint=regexp(fstringsel,[listargouts{i+1} '\s{0,7}:']);
        elseif i==nargout-1
            
            if strcmp(listargouts{end},'varargout') ==0
                endpoint=regexp(fstringsel,[listargouts{i+1} '\s{0,7}:']);
                if isempty(endpoint)
                    % warning('FSDA:wrongOutDescription',)
                    errmsg=['Error in processing output argument  ''' listargouts{i} '''\n' ...
                        'Parser could not find string:  ''' listargouts{i+1} '       :''\n' ...
                        'Endpoint for the description of output argument ''' listargouts{i} '''not found'];
                    error('FSDA:wrongOutDescription',errmsg)
                end
            else
                % In this case there are also optional arguments
                endpoint=regexp(fstringsel,'Optional [Oo]utput:');
                if isempty(endpoint)
                    error('FSDA:missOuts','varagout is present but input .m file does not contain ''Optional Output:'' string')
                end
            end
            
        else
            inipointSeeAlso=regexp(fstringsel,'%\s*See also','once');
            % fstringsel1=fstringsel(iniref+1:end);
            
            % inipointSeeAlso=regexp(fstringsel1,'See also','once');
            
            if isempty(inipointSeeAlso)
                disp('Please check .m input file')
                error('FSDA:missOuts','Input .m file does not contain ''See also:'' string')
            end
            
            inipointMoreAbout=regexp(fstringsel,'More About:','once');
            if ~isempty(inipointMoreAbout)
                MoreAbout=fstringsel(inipointMoreAbout+15:inipointSeeAlso-1);
                posPercentageSigns=regexp(MoreAbout,'[^"^%]%')+1;
                
                MoreAbout(posPercentageSigns)=[];
                MoreAboutHTML=MoreAbout;
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
                disp(['Starting point of the description of ouput argument' '' listargouts{i} ''])
                disp([fstringsel(initmp:initmp+50) '.....'])
                disp(['Final point of the description of ouput argument ''' listargouts{i} ''])
                endtmp=endpoint(1);
                disp(['....' fstringsel(endtmp-100:endtmp-1)]);
                disp(['FSDA:WrongOut','Could not process correctly output argument ' listargouts{i}])
                disp('Probably the list of the output arguments has been inverted')
            end
        catch
            disp(['FSDA:WrongOut: ','Could not process correctly output argument ' listargouts{i}])
            
            disp('Below is the wrong description that parser publishFS has extracted')
            starterr=(inipoint(1)+length(listargouts{i})+2);
            disp(fstringsel(starterr:starterr+300))
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
            
            [listStructureArgs]=formatHTMLstructure(descrioutput,outi);
            listOutArgs{i,5}=listStructureArgs;
            
            
            
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
                    error('FSDA:publishFS:WrongOut',[ warn1 'Sentence ''' descrioutput ''' must contain at least \n two full stops each followed by a description'])
                    
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
            
            listOutArgs{i,2}=preamble;
            listOutArgs{i,3}=descroutputtitl;
            listOutArgs{i,4}=descrioutput;
            
        end
        
    end
else
    
    inipointSeeAlso=regexp(fstringsel,'See also','once');
    
    if isempty(inipointSeeAlso)
        disp('Please check .m input file')
        error('FSDA:missOuts','Input .m file does not contain ''See also:'' string')
    end
    
    
    
    inipointMoreAbout=regexp(fstringsel,'More About:','once');
    if ~isempty(inipointMoreAbout)
        MoreAbout=fstringsel(inipointMoreAbout+15:inipointSeeAlso-1);
        posPercentageSigns=regexp(MoreAbout,'%');
        MoreAbout(posPercentageSigns)=[];
        MoreAboutHTML=MoreAbout;
    else
        MoreAbout='';
        MoreAboutHTML='';
    end
    % Remove from string descri leading and trailing white spaces
    
end



%% CREATE MORE ABOUT SECTION
if isempty(MoreAboutHTML)
    MoreAbout='';
end


%% REFERENCES
inipointAcknowledgements=regexp(fstring,'%\s*Acknowledgements:');
if isempty(inipointAcknowledgements)
    Acknowledgements='';
end

iniref=regexp(fstring,'%\s*References:');
if isempty(iniref)
    warning('FSDA:missInp',['File ' name '.m does not contain string ''References:'''])
    % refsargs={''};
    refsargs='';
else
    
    inipointCopyright=regexp(fstring,'Copyright');
    
    
    if ~isempty(inipointAcknowledgements)
        Acknowledgements=fstring(inipointAcknowledgements+19:inipointCopyright-1);
        posPercentageSigns=regexp(Acknowledgements,'%');
        Acknowledgements(posPercentageSigns)=[];
        
    else
        Acknowledgements='';
        inipointAcknowledgements=Inf;
    end
    
    endref=min([inipointAcknowledgements(:);inipointCopyright(:)]);
    
    if isempty(endref)
        disp('Please check .m input file')
        error('FSDA:missOuts','Input .m file does not contain ''Copyright'' string')
    end
    
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
    
    
end



%% ACKNOWLEDGEMENTS



%% SEE ALSO

iniref=regexp(fstring,'%\s*See also','once');
fstring=fstring(iniref+1:end);
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

% Store in cell listSeeAlso, See also items
listSeeAlso=cell(nseealso,1);

for i=1:nseealso
    if nseealso==1
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
    if ~isempty(Seealsoitem)
        Seealsoitem=strtrim(Seealsoitem);
        str=which(Seealsoitem);
        
        % Store See also item
        listSeeAlso{i}=Seealsoitem;
        
        if isempty(str)
            error('FSDA:publishFS:WrngSeeAlso',['Wrong reference in "See Also:" cannot find a reference to ' Seealsoitem ]);
        else
            % Check if the reference is towards a file present in the FSDA toolbox
            FSDAtoolboxfile=regexpi(str,'FSDA', 'once');
            
            % Create the string which contains the 'destination' of the hyperlink
            % DestHyperLink is the 'destination' of the hyperlink
            if ~isempty(FSDAtoolboxfile) % reference is towards a function of the FSDA toolbox
            else % reference is towards a MATLAB function or a function of another toolbox
                pathdocroot=docroot;
                % Find path of .html documentation file
                pathExtHelpFile=findFile(pathdocroot,'InclFiles',[Seealsoitem '.html']);
                
                if isempty(pathExtHelpFile)
                    error('FSDA:publishFS:WrngSeeAlso',['cannot find a reference to doc file ' Seealsoitem '.html']);
                end
                
                % replace '\' with '/'
                % DestHyperLink=['matlab:web(fullfile(docroot,''' addSubPath '.html''))'];
            end
        end
        
    end
    
end

% Seealso='';

%% CLOSE TAGS SECTION



%% Create output structure
out=struct;
% save title
out.titl=name;
% save purpose
out.purpose=purpose;
% Save description
out.description=gendesc;
% Save input arguments (required +optional)
listInpArgs(:,7)={'0'};
listInpArgs(1:nREQargin,7)={'1'};
out.InpArgs=listInpArgs;


out.OptArgs=listOptArgs;
% Save output arguments
out.OutArgs=listOutArgs;
% Save more about
out.MoreAbout=MoreAbout;
% Save Acknowledgements
out.Acknowledgements=Acknowledgements;
% Store references
out.References=refsargs;
% listSeeAlso
out.SeeAlso=listSeeAlso;
% Store cell containing examples
out.Ex=listEx;
% Store cell contaning extra examples
out.ExtraEx=listExtraEx;



%% Check if all name pairs arguments are commented inside the HTML
% Check that (if varagin is present all optional arguments contained in
% "options=struct('................');" are present in the column of
% out.OptArgs or inside the first column of listOptArgs
OptMisMatch=cell(1,3);
OptMisMatch{1,2}='Options described';
OptMisMatch{1,3}='Options used';

if ~isempty(OptArgsVarargin)
    posoptionsini=regexp(fstring,'options\s*=\s*struct(');
    if isempty(posoptionsini)
        warning('FSDA:publishFS:WrongUseVarArgInt','varargin is used but structure containing the options of varargin has not been defined inside .m file')
    else
        posoptionsini=posoptionsini(1);
        posoptionsfin=regexp(fstring,');');
        posoptionsfin=posoptionsfin(posoptionsfin>posoptionsini);
        posoptionsfin=posoptionsfin(1);
        NamePairs=fstring(posoptionsini:posoptionsfin-1);
        findFirstRoundBracket=regexp(NamePairs,'(');
        NamePairs=NamePairs(findFirstRoundBracket+1:end);
        
        % Remove (if present) "..." symbols, carriage returns and white spaces
        % from, string NamePairs
        posThreeFullStops=regexp(NamePairs,'\.');
        NamePairs(posThreeFullStops)=[];
        NamePairs=strtrim(NamePairs);
        posCR=regexp(NamePairs,'\n');
        NamePairs(posCR)=[];
        posLF=regexp(NamePairs,'\r');
        NamePairs(posLF)=[];
        posWS=regexp(NamePairs,'\s');
        NamePairs(posWS)=[];
        
        
        % count number of commas (it must be an odd number). Check that this is the
        % case
        poscommas=regexp(NamePairs,',');
        if length(poscommas)/2==floor(length(poscommas)/2)
            error('FSDA:publishFS:WrongInputOpt','Name pairs must be in pairs: something wrong')
        else
            
            numberOptArgs=(length(poscommas)+1)/2;
            OptArgsUsed=cell(numberOptArgs,1);
            ij=1;
            for j=1:2:length(poscommas)
                if j>1
                    OptArgsUsed{ij}=NamePairs(poscommas(j-1)+2:poscommas(j)-2);
                else
                    % regexp(NamePairs,'\');
                    OptArgsUsed{ij}=NamePairs(2:poscommas(j)-2);
                end
                ij=ij+1;
            end
        end
        % Now compare listOptArgs (described options) with OptArgsUsed
        % (effectively used options)
        if strcmp(Display,'iter-detailed')
            disp('Check if Name/Pairs optional input arguments are documented')
        end
        OptMisMatch=CompareDescribedUsed(listOptArgs(:,1),OptArgsUsed);
        
        %     if ~isequal(OptArgsDescribed,OptArgsUsed)
        %         warning('Options described are different from Option effectively used')
        %         lused=length(OptArgsUsed);
        %         ldesc=length(OptArgsDescribed);
        %         if ldesc<lused
        %             OptArgsDescribed=[OptArgsDescribed; cell(lused-ldesc,1)];
        %         elseif ldesc>lused
        %             OptArgsUsed=[OptArgsUsed; cell(ldesc-lused,1)];
        %         end
        %         disp('Options described --  Options used')
        %         disp([OptArgsDescribed OptArgsUsed])
        %     end
    end
end
out.InpArgsMisMatch=OptMisMatch;

%% Check if all fields of output arguments which are struct are commented inside the HTML file
% In other words, suppose that one of the output arguments is called out
% and that the following fields have been commented
% out.a, out.b and out.c
% the code below checks that fields out.a out.b and out.c are effectively
% used inside the code.
if evalCode ==1
    if strcmp(Display,'iter-detailed')
        disp('Check if all fields of output arguments which are structure are documented')
    end
    OutArgsMisMatch=cell(100,3);
    OutArgsMisMatch{1,2}='Fields described';
    OutArgsMisMatch{1,3}='Fields used';
    ik=1;
    
    for i=1:size(listOutArgs,1)
        % If the fourth column of listOutArgs is a cell then we are dealing
        % with a structure
        if iscell(listOutArgs{i,5})
            cellOut=listOutArgs{i,5};
            OutputDescribed=cellOut(:,1);
            listouti=listOutArgs{i,1};
            
            % Find at least one example which is executed inside
            % listExtraEx or listEx
            if ~isempty(listExtraEx)
                boo=cell2mat(listExtraEx(:,4))==1;
            else
                boo=0;
            end
            if sum(boo)>0
                listExtraExSel=listExtraEx(boo,:);
                WhereToSearch=listExtraExSel{1,3};
            else
                boo=cell2mat(listEx(:,4))==1;
                if sum(boo)>0
                    listExSel=listEx(boo,:);
                    WhereToSearch=listExSel{1,3};
                else
                    WhereToSearch='';
                end
            end
            
            if ~isempty(WhereToSearch)
                % Remove from string WhereToSearch symbols [ or ]
                posSquarePar=regexp(WhereToSearch,'[');
                WhereToSearch(posSquarePar)=[];
                posSquarePar=regexp(WhereToSearch,']');
                WhereToSearch(posSquarePar)=[];
                
                findBegFileNameString=regexp(WhereToSearch,['\w*=\s*' name]);
                findEndFileNameString=regexp(WhereToSearch,['=\s*' name]);
                if size(listOutArgs,1)==1
                    try
                        varToSearch=strtrim(WhereToSearch(findBegFileNameString(1):findEndFileNameString(1)-1));
                    catch
                        varToSearch='';
                    end
                else
                    varToSearch='';
                end
                
                % Check if variable varToSearch/listouti exists in main workspace
                if ~isempty(varToSearch)
                    ChkExistinWS=evalin('base', ['who(''' varToSearch ''')']);
                else
                    ChkExistinWS=evalin('base', ['who(''' listouti ''')']);
                end
                
                % If it is empty check if there are variables in the main workspace whose name contains string listouti
                if isempty(ChkExistinWS)
                    ChkExistinWS=evalin('base', ['who(''*' listouti '*'')']);
                    % If there is more than 1 instance in output just take the
                    % first one
                    if ~isempty(ChkExistinWS)
                        NameToSearchinWS=ChkExistinWS{1,1};
                    else
                        NameToSearchinWS=listouti;
                    end
                else
                    NameToSearchinWS=listouti;
                end
                
                try
                    [~]=evalin('base', NameToSearchinWS);
                    OutputProduced = fieldnames(evalin('base', NameToSearchinWS));
                catch
                    warning('FSDA:publishFS:WrongOut',['In the examples which were executed output argument containing string ' listouti ' which is of class struct has not been found'])
                    OutputProduced='';
                end
                if ~isempty(OutputProduced)
                    if strcmp(Display,'iter-detailed')
                        disp(['Analysis of output argument: ''' listouti ''''])
                    end
                    OutiMisMatch=CompareDescribedUsed(OutputDescribed,OutputProduced);
                    if size(OutiMisMatch,1)>1
                        OutArgsMisMatch{ik+1,1}=listouti;
                        OutArgsMisMatch(ik+2:ik+size(OutiMisMatch,1),:)=OutiMisMatch(2:end,:);
                        ik=ik+size(OutiMisMatch,1);
                    end
                end
            else
                warnmsg=['No example has been executed therefore'...
                    ' coherance control in output structure is impossible'];
                warning('FSDA:publishFS:WrongOut',warnmsg);
            end
        end
    end
    out.OutArgsStructMisMatch=OutArgsMisMatch(1:ik,:);
else
    out.OutArgsStructMisMatch='';
end
if length(linkHTML)~=2
    out.linkHTMLMisMatch=true;
else
    out.linkHTMLMisMatch=false;
end


end





function [listStructureArgs]=formatHTMLstructure(descriinput,StructureName)

[iniA,finA]=regexp(descriinput,['\s{0,8}' StructureName '\.\w*\s*=']);
%[iniA,finA]=regexp(descriinput,['\s{0,8}' StructureName '\.\w*\(\w*\)\s*=']);

% Sometimes field names are repeated and therefore it is necessary to find
% just the first instance. THe following lines do that
FieldNamesStruct=cell(length(iniA),1);
for k=1:length(iniA)
    FieldNamesStruct{k}=strtrim(descriinput(iniA(k):finA(k)-1));
end

[~,seluni]=unique(FieldNamesStruct,'stable');
ini=iniA(seluni);
fin=finA(seluni);



if isempty(ini)
    disp('Probably ":" symbols  must be replaced with "=" symbols in out description')
    error('FSDA:MissingArg',['Parser cannot find string \n''' StructureName '.xxxx'' = \n for output structure ' StructureName])
else
end


posREMARK=regexp(descriinput,'\r\s*\r\s*remark','once','ignorecase');
% If there is the string REMARK it means that all that which is after
% remark is a general statement which does not have to be contained in the
% table with the associated fields
if ~isempty(posREMARK)
    boo=ini<posREMARK;
    ini=ini(boo);
    fin=fin(boo);
end


% listStructureArgs = cells which will contains the name/values
% pairs of the structure
% fieldnames will go into the first column of the table
% the content of the field names will fo into the second column of
% the table
if ~isempty(ini)
    listStructureArgs=cell(length(ini),2);
    
    for k=1:length(ini)
        % fin(i)-1 because character '=' does not have to be extracted
        StructFieldName=descriinput(ini(k):fin(k)-1);
        % Remove from string opti leading and trailing white spaces
        StructFieldName=strtrim(StructFieldName);
        StructFieldName=StructFieldName(length(StructureName)+2:end);
        
        % Store name in the first column
        listStructureArgs{k,1}=StructFieldName;
        % Store field content in the 2nd col of listStructureArgs
        if k<length(ini)
            StructFieldCont=descriinput(fin(k)+1:ini(k+1)-1);
        else
            
            if ~isempty(posREMARK)
                StructFieldCont=descriinput(fin(k)+1:posREMARK-1);
            else
                StructFieldCont=descriinput(fin(k)+1:end);
            end
        end
        
        % Store content in the second column
        StructFieldCont=strtrim(StructFieldCont);
        if strcmp(StructFieldCont(end),'-')
            StructFieldCont=StructFieldCont(1:end-1);
        end
        listStructureArgs{k,2}=StructFieldCont;
    end
    
    % rowtodel = vector which contains the duplicate rows of
    % listStruArgs which have to be deleted
    inisel=1:size(listStructureArgs,1);
    
    % Check if cell listStruArgs contains duplicates in the first column
    for j=2:size(listStructureArgs,1)
        if strcmp(listStructureArgs{j,1},listStructureArgs{j-1,1})
            listStructureArgs{j-1,2}=[listStructureArgs{j-1,2} listStructureArgs{j,1} listStructureArgs{j,2}];
            inisel(j)=999;
        end
    end
    
    
    % Add the Remark after the table, if it is present
    
    % use capital letter for the first word.
else
    listStructureArgs='';
end

end


function OptMisMatch=CompareDescribedUsed(OptArgsDescribed,OptArgsUsed)

OptArgsDescribed=sort(OptArgsDescribed);
OptArgsUsed=sort(OptArgsUsed);


[~,ia,ib] = setxor(OptArgsDescribed,OptArgsUsed);
% lused=length(OptArgsUsed);
% ldesc=length(OptArgsDescribed);
OptMisMatch=cell(length(ia)+length(ib)+1,3);
OptMisMatch{1,2}='Options described';
OptMisMatch{1,3}='Options used';
ij=1;
if ~isempty(ia)
    disp('Elements described but not used')
    % OptMisMatch{ij+1:ij+length(ia),1}=OptArgsDescribed{ia};
    OptMisMatch(ij+1:ij+length(ia),1)=OptArgsDescribed(ia);
    OptMisMatch(ij+1:ij+length(ia),2:3)=num2cell([true(length(ia),1),false(length(ia),1)]);
    ij=ij+length(ia);
    disp(OptArgsDescribed(ia))
end
if ~isempty(ib)
    disp('Elements used but not described')
    OptMisMatch(ij+1:ij+length(ib),1)=OptArgsUsed(ib);
    OptMisMatch(ij+1:ij+length(ib),2:3)=num2cell([false(length(ib),1),true(length(ib),1)]);
    disp(OptArgsUsed(ib))
end
end





%FScategory:UTIHELP