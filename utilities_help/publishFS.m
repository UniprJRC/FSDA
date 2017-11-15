function out=publishFS(file,varargin)
%publishFS enables to create automatic HELP FILES from structured .m function files
%
%<a href="matlab: docsearchFS('publishFS')">Link to the help function</a>
%
%   publishFS creates HTML files from structured .m file. To understand the
%   characteristics the .m file must have, please see the "More About" section
%
% Required input arguments:
%
%    file:         MATLAB File. String. Full or partial
%                  path of the MATLAB file for which Structured Matlab
%                  HTML help has to be created
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
% outputDir : Output folder. String.
%             Output folder to which the HTML document is saved, specified
%             as the comma-separated pair consisting of 'outputDir' and the
%             full path. You must specify the full path as a string, for
%             example 'C:\PublishedOutput'.
%             The default value, '', specifies the (FSDA root)\helpfiles\FSDA
%             path.
%             Remark - outputDir must be a valid path.
%             Example - 'outputDir','C:'
%             Data Types - string
% imagesDir : Output folder of png images. String.
%             Output folder to which the images attached to the HTML
%             document are saved, specified as the comma-separated pair
%             consisting of 'outputDir' and the full path. You must specify
%             the full path as a string, for example
%             'C:\PublishedOutput'.
%             The default value, '', specifies the
%             "(FSDA root)\helpfiles\FSDA\images" path.
%             Remark - if imageDir is not specified but outputDir is
%             specified images will be save into the same folder of the
%             HTML output file
%             Remark - imagesDir must be a valid path.
%             Example - 'imagesDir','C:'
%             Data Types - string
% evalCode :  Option to run code. Logical. Option to evaluate code of the
%             examples in the input .m files enclosed in tags "%{" "%}" whose
%             first line starts with symbols "%%".
%             If evalcode is true the code associated with the examples
%             which start with symbols '%%' will be run and the output will
%             be automatically included into the HTML output file. The
%             images will be saved in subfolder images_help of the
%             outputDir. The default value of evalCode is true.
%             Example - 'evalCode',false
%             Data Types - logical
% write2file: Option to write HTML file. Logical. Option which specifies
%             whether HTML file must be created or if just structure out
%             must be created. The default value of write2file is true,
%             that is html file is created
%             Example - 'write2file',false
%             Data Types - logical
%
% Output:
%
%  The output consists of a structure 'out' containing the following fields:
% out.title     = title of HTML file. String.
%                   String to be included in HTML tag title
% out.purpose   = purpose of the routine. String.
%                 String forming second row of output HTML file
%                 String forming second row of output HTML file
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
%                 input argument is a struct, all
%                 information about the fields of the structure will be
%                 included in the 8th column.
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
% out.OutArgs   = Required and Optional (varargout) output arguments. Cell.
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
%                 the example must be executed or not. If it is equal to 1
%                 example is executed.
%    out.ExtraEx= Extra Examples. cell. Cell of length u containing the u
%                 extra examples.
%                 First column= title of the example;
%                 Second column = detailed description;
%                 Third column = MATLAB code;
%                 Fourth column = dummy variable which indicates whether
%                 the extra example must be executed or not. If it is equal to 1
%                 extra  example is executed.
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
%                 it means that options 'nomes', 'refsteps' and 'reftol'
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
%         The .m file (which must be located on the MATLAB path or on the current
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
%         The last line may start with the words "Data Types -" and
%         contains the specification of a particular input argument 
%         (e.g. Data Types - single | double). For example, suppose that the .m
%         routine which has to be processed has two required input
%         arguments which are respectively called y and X, then the
%         accepted format is as follows.
%
%                       Required input arguments:
%
%                        y:         Response variable. Vector. Response variable,
%                                   specified as a vector of length n, where n is
%                                   the number of observations. Each entry in y is
%                                   the response for the corresponding row of X.
%                                   Data Types - single | double.
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
%                                   Data Types - single | double.
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
%                                   (default), else no constant term will be included.
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
%         IMPORTANT NOTICE: Similarly to what happened for each input argument, if
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
%         [10] If the .m file contains string  "More About:" a particular section
%         called "More about" in the HTML file will be created
%         (just before See Also).
%         [11] If the .m file contains string "Acknowledgements:" then a particular
%         section named "Acknowledgements" will be created just above the
%         references.
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
%         REMARK 3: publishFS uses javascript matlab-highlighter.min.js in order to
%         automatically color the examples in the HTML file.
%         -----------------------------------------------------------------------:.
%         REMARK 4: publishFS uses MathJax javascript in order to interpret the
%         equations inside the .m file written in Latex style. For in line
%         equations both symbols dollar dollar  and backslash(  backslash) are accepted.
%         For single line equations symbols backslash[ backslash] must be used.
%         -----------------------------------------------------------------------:.
%         REMARK 5: if there are not enough examples in the .m file the procedure
%         still runs but a warning will be automatically produced.
%         -----------------------------------------------------------------------:.
%         REMARK 6: the output file to be correctly viewed must be located in a
%         folder which contains a subfolder named includesFS containing the files
%         present inside
%         (home FSDA) filesep helpfiles filesep FSDA filesep includesFS.
%         Therefore if the
%         the directory which contains the output file is not located inside
%         (home FSDA) filesep helpfiles filesep FSDA subfolder
%         includesFS must be copied into the current folder.
%         -----------------------------------------------------------------------:.
%         REMARK 7: strings are HTML formatted as follows. Every time there is
%         symbol ". one_or_more_space symbol_of carriage_return" or
%         ": one_or_more_space symbol_of carriage_return" the parser adds HTML
%         string '</p>/<p>' after just symbol "."  or symbol ":".
%         This is done using subfunction named formatHTML at the end of this file.
%         subfunction formatHTMLwithMATHJAX is even more general because it applies
%         function formatHTML just to the parts of the input string which are not
%         enclosed inside symbols 'backslash[ backslash]'.
%         -----------------------------------------------------------------------:.
%         REMARK 8: parser automatically puts an hyperlink every time in the text
%         there is something which starts with "http": or every time there is a
%         reference to a .m file. For example the sentence: "More details can be
%         found in routine ncx2mixtcdf.m" is converted as follows.
%         "More details can be found in routine
%         <a href="ncx2mixtcdf.html">ncx2mixtcdf</a>".
%         Similarly, a sentence such as:
%         The full paper can be downloaded from "http://www.mysite.org".
%         is converted as follows:
%         "The full paper can be donwloaded from
%         <a href="http://www.mysite.org"> http://www.mysite.org </a>".
%         -----------------------------------------------------------------------:.
%
% See also: publish
%
%
% References:
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('publishFS')">Link to the help function</a>
%
% Last modified 14-06-2016

% Examples:


%{
  % Create file FSRmdr.html starting from file FSRmdr.
  out=publishFS('FSRmdr','evalCode',false,'Display','iter-detailed')
%}

%{
  % Option Display.
  % Create file FSRmdr.html starting from file FSRmdr and
  % display detailed information about the Input, Output and Optional
  % arguments.
  out=publishFS('FSRmdr','evalCode',false,'Display','iter-detailed')
%}

%{
  % Option outputDir.
  % Create HTML file with embedded images in folder C:\tmp.
  out=publishFS('FSRmdr','evalCode',true,'outputDir','C:\tmp')
%}

%% Beginning of code

if ~ischar(file)
    error('FSDA:publishFS:WrongInputOpt','input must be a string containing the name of the file.');
end

evalCode=true;
write2file=true;
Display='none';

% % Use file separator of current operating system
% % \ = Windows
% % / = Unix
fsep=filesep;

% Write output file in subfolder \(FSDAroot)\helpfiles\FSDA
FileWithFullPath=which('docsearchFS.m');
[pathFSDAstr]=fileparts(FileWithFullPath);

outputDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA'];
imagesDir=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA' fsep 'images'];

%  pathstr=docroot;
%  outputDir=[pathstr fsep 'FSDA'];
%  imagesDir=[pathstr fsep 'FSDA' fsep 'images'];

% Check matlab version. The help system has changed from 2012b and jar
% files have been deleted therefore, given that in the see also part parser
% checks whether a corresponding funciton xxx.html exists, we must take
% into account that this function will not be found in releases older tha
% 2012b
matlabversion=verLessThan('matlab','8.1.0');

if nargin>1
    options=struct('evalCode',evalCode,'Display',Display,'outputDir',outputDir,'imagesDir',imagesDir,'write2file',true);
    
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
    
    evalCode=options.evalCode;
    write2file=options.write2file;
    outputDir=options.outputDir;
    Display=options.Display;
    checkimageDir = strcmp(UserOptions,'imagesDir');
    checkoutputDir = strcmp(UserOptions,'outputDir');
    
    if sum(checkimageDir)==0 && sum(checkoutputDir)==1
        imagesDir=outputDir;
    elseif sum(checkimageDir)==1 && sum(checkoutputDir)==1
        imagesDir=options.imagesDir;
    else
    end
    
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
    disp('Please check HTML input file')
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

% fin=fint(1:length(ini));
% [ini,fin]=regexp(fstringselOpt,'%\s*\w*\s*:\s*\w');
%
% [iniCR,finCR]=regexp(fstringselOpt,'%\s*\w*\s*:\s*\r');


% listOptArgs = list which contains all optional arguments (8 columns)
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
listOptArgs=cell(length(ini),8);

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
            listOptArgs{ij,4}=descrlong(1:CheckExample-1);
            
            % The first word of example code must be embedded around tags <code> </code>
            examplecode=descrlong(CheckExample+9:Datatypes-1);
            % Store string containing the examples before adding
            % HTML strings '<code>'  '</code>
            listOptArgs{ij,7}=examplecode;
            
            posspace=regexp(examplecode,'      ');
            examplecode=['<code>' examplecode(1:posspace-1) '</code>' examplecode(posspace:end)];
            listOptArgs{ij,5}=strtrim(examplecode);
            listOptArgs{ij,6}=descrlong(Datatypes+13:end);
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
    '<link href="includesFS/bootstrap.min.css" rel="stylesheet" type="text/css">'...
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
    '<script type="text/javascript" src="includesFS/jquery-latest.js"></script>\r'...
    '<script>\r'...
    '$(document).ready(function(){\r'...
    '      $("#divtop").load("includesFS/top.html");\r'...
    '      $("#divbottom").load("includesFS/bottom.html");\r'...
    '});\r'...
    '</script>\r'...
    '<script src="includesFS/jquery-latest.js" type="text/javascript"></script>\r'...
    '<link href="includesFS/site6.css?201505100807" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/site6_lg.css?201505100807" media="screen and (min-width: 1200px)" rel="stylesheet">\r'...
    '<link href="includesFS/site6_md.css?201505100807" media="screen and (min-width: 992px) and (max-width: 1199px)" rel="stylesheet">\r'...
    '<link href="includesFS/site6_sm+xs.css?201505100807" media="screen and (max-width: 991px)" rel="stylesheet">\r'...
    '<link href="includesFS/site6_sm.css?201505100807" media="screen and (min-width: 768px) and (max-width: 991px)" rel="stylesheet">\r'...
    '<link href="includesFS/site6_xs.css?201505100807" media="screen and (max-width: 767px)" rel="stylesheet">\r'...
    '<link href="includesFS/site6_offcanvas.css?201505100807" rel="stylesheet" type="text/css">\r'...
    '<script src="includesFS/l10n.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/SHAREDdocscripts.js"></script>\r'...
    '<script src="includesFS/PRODUCTdocscripts.js"></script>\r'...
    '<script src="includesFS/mw.toc.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.imagescaling.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.imageanimation.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/jquery.highlight.js"></script>\r'...
    '<script src="includesFS/bootstrap.min.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/global.js"></script>\r'...   ' enables scrolling
    '<script src="includesFS/bottom.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/suggest.js" type="text/javascript"></script>\r'... % for search engine
    '<script src="includesFS/underscore-min.js"></script>\r'...                 % for search engine
    '<link href="includesFS/reset.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/960.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center_installed.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center_installed.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includes/product/css/doc_center_print.css" media="print" rel="stylesheet" type="text/css">\r'...
    '</head>\r'...
    '<body  onload="highlightMATLABCode();" id="responsive_offcanvas">\r'];


searchenginestring=[' ' ' '];


metacontent=sprintf([beforemetacontent purpose aftermetacontent searchenginestring]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Insert navigation bar on top of the page
% necessary to insert sprintf later because there is symbol % in 100%
% insnav=['<table border="0" cellpadding="0" cellspacing="0" class="nav" width="100%">'...
%     sprintf(['\r<tr valign="top">\r'...
%     '<td align="left" width="20"><a href=" ">\r'...
%     '<img align="bottom" alt="" border="0" src="images_help/b_prev.gif"></a>\r'...
%     '</td>\r'...
%     '<td align="left">left</td>\r'...
%     '<td>&nbsp;</td>\r'...
%     '<td align="right"> right</td>\r'...
%     '<td align="right" width="20"><a href=" ">\r'...
%     '<img align="bottom" alt="score" border="0" src="images_help/b_next.gif"></a></td>\r'...
%     '</tr>\r'...
%     '</table>'])];

insnav=sprintf('<div id="divtop"></div>');



%% CONTAINER + SINTAX SECTION

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

[gendescini]=regexp(fstring,'/a','once');
[gendescfin]=regexp(fstring,'Required input arguments','once');
gendesc=fstring(gendescini+6:gendescfin-1);
posPercentageSigns=regexp(gendesc,'%');
gendesc(posPercentageSigns)=[];
% Remove from string descri leading and trailing white spaces
gendesc=strtrim(gendesc);

if length(gendesc)<3
    htmlsitecont='';
    gendesc='';
else
    htmlsitecont=['<p>' formatHTMLwithMATHJAX(gendesc) '</p>'];
end

finsitecont=sprintf(['<h2 id="syntax">Syntax</h2>\r'...
    '<div class="syntax_signature">\r'...
    '<div class="syntax_signature_module">\r'...
    '<ul>\r']);

sitecont=[insnav inisitecont finsitecont];

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
    
    %sintax{1}=[outargs(2:commasOut(1)-1) '=' name strtrim(inputargs(1:optargs1-2)) ')'];
    %sintax{2}=[outargs(2:commasOut(1)-1) '=' name strtrim(inputargs(1:optargs1-2)) ',Name,Value)'];
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

sintaxhtml='';

for j=1:length(sintax)
    sintaxhtml= [sintaxhtml '<li><code class="synopsis">'  sintax{j} '</code>' ...
        '<span class="syntax_example">'...
        '<a class="intrnllnk" href="' name '.html#Example_' num2str(j) '">' ...
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
    '							<h2 id="Description">Description</h2>\r'...
    '							<div class="descriptions">\r'...
    '								<div class="description_module">\r']);


descriptionhtml=htmlsitecont;
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
    
    descriptionini=sprintf(['<div class="description_element">\r'...
        '	<p class="syntax_example">\r'...
        '	<a class="intrnllnk" href="' name '.html#Example_'  num2str(j) '">\r'...
        '	example</a></p>\r'...
        '	<p><span itemprop="syntax"><code>\r']);
    
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
                listargouts=cell(1,1);
                listargouts{1}=outi;
            end
        end
    else
        listargouts='';
        outstring='';
    end
    
    % Locate in  expression [out1,out2,...]=namefunc(inp1,inp2,...) the
    % position of open parenthesis sign
    [startIndex] = regexp(sintaxj,'(');
    inps=sintaxj(startIndex+1:end-1);
    commaspos=regexp(inps,',');
    if isempty(commaspos)
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
        elseif isempty(inpi)
            inpstring='';
        else
            inpstring=sprintf(['<a class="intrnllnk" href="#inputarg_' inpi '"><code>' inpi '</code></a>\r']);
        end
    end
    
    if ~isempty(outs)
        description=[outstring '=' name '(' strtrim(inpstring) ')'];
    else
        description=[name '(' strtrim(inpstring) ')'];
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
    % remove % signs from strdescrEx
    % OLD
    %posPercentageSigns=regexp(strdescrEx,'%');
    % strdescrEx(posPercentageSigns)=[];
    
    posPercentageSigns=regexp(strdescrEx,'\D%');
    strdescrEx(posPercentageSigns+1)=[];
    
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

% imgtemplate = image to include for the examples which can be executed
imgtemplate='<img alt="" src="images_help/M.gif" style="width: 12px; height: 12px"> ';

% the examples which are inside %{   %} are put here.
% The first sentence which ends with a full stop is the title of the example
iniexamples=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<div class="examples">\r'...
    '<h2 id="Examples">Examples</h2>\r'... % start of expandable examples
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
        '<div id="Example_' num2str(j) '">\r'...
        '</div>\r'...
        '<h3 class="expand"><span>\r'...
        '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        '<span class="example_title"> <li>']) addimg listEx{j,1} sprintf(['</li></span></a></span></h3>\r'...
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


closeexamples=sprintf(['</div>\r'... % close div id="expandableExamples
    '<p> </p>\r']);
% Related examples are below
iniRelatedExamples='';
RelatedExamples='';
if length(startIndexEx)>length(sintax)
    iniRelatedExamples=sprintf('<h3 id="ExtraExamples" class="bottom_ruled">Related Examples</h3>\r');
    
    for j=1:size(listExtraEx,1)
        
        if listExtraEx{j,4}==1
            addimg=imgtemplate;
        else
            addimg='';
        end
        
        RelatedExamples=[RelatedExamples  sprintf(['<div id="ExtraExample_' num2str(j) '" class="example_module expandableContent">\r'...
            '<div id="ExtraExample_' num2str(j) '">\r'...
            '</div>\r'...
            '<h3 class="expand"><span>\r'...
            '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            '<span class="example_title"><li>']) addimg listExtraEx{j,1} sprintf(['</li></span></a></span></h3>\r'...
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

%% CREATE REQUIRED INPUT ARGUMENTS SECTION OF HTML file
iniReqInputArgs=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<h2 id="Inputs">Input Arguments</h2>\r'...
    '<div class="expandableContent">\r'...
    '<div class="arguments">\r'...
    '<div class="input_argument_container">\r'...
    '<p class="switch">\r'...
    '<a class="expandAllLink" href="javascript:void(0);">\r'...
    'expand all</a></p>']);


reqargs='';
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
    DescrInputToSplit=[DescrInputToSplit ' '];
    
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
        
        [descrlongHTML,listStructureArgs]=formatHTMLstructure(descrlong,inpi);
        listInpArgs{i,4}=descrlongHTML;
        listInpArgs{i,7}=descrlong;
        
        listInpArgs{i,8}=listStructureArgs;
        
        jins=6;
    else
        
        
        if i<=nREQargin
            
            Datatypes=regexp(descrlong,'Data Types -','once');
            if ~isempty(Datatypes)
                listInpArgs{i,6}=descrlong(Datatypes+13:end);
                
                descrlong=descrlong(1:Datatypes-1);
                
                
                descrlongHTML=formatHTMLwithMATHJAX(descrlong);
                listInpArgs{i,4}=descrlongHTML;
                listInpArgs{i,7}=descrlong;
            else
                
                descrlongHTML=formatHTMLwithMATHJAX(descrlong);
                
                listInpArgs{i,4}=descrlongHTML;
                listInpArgs{i,7}=descrlong;
                if strcmp(Display,'iter-detailed')
                    warning('FSDA:publishFS:MissingDataType',['Input argument ''' inpi ''' does not contain DataType line, by default string  ''single| double'' has been added'])
                end
                
                listInpArgs{i,6}='single| double';
            end
            jins=6;
        else
            
            % Check if descrlong contains
            % Example - and Data types -
            
            CheckExample=regexp(descrlong,'Example -','once');
            if ~isempty(CheckExample)
                Datatypes=regexp(descrlong,'Data Types -','once');
                descrlonginp=descrlong(1:CheckExample-1);
                descrlongHTML=formatHTML(descrlonginp);
                listInpArgs{i,4}=descrlongHTML;
                listInpArgs{i,7}=descrlonginp;
                
                % The first word of example code must be embedded around tags <code> </code>
                examplecode=descrlong(CheckExample+10:Datatypes-1);
                posspace=regexp(examplecode,'      ');
                examplecode=['<code>' examplecode(1:posspace-1) '</code>' examplecode(posspace:end)];
                listInpArgs{i,5}=strtrim(examplecode);
                listInpArgs{i,6}=descrlong(Datatypes+13:end);
                jins=6;
                
            else
                listInpArgs{i,4}=descrlong;
                listInpArgs{i,7}=descrlong;
                warning('FSDA:publishFS:MissingExample',['Optional input argument ''' inpi ''' does not contain an Example'])
                jins=5;
                
            end
            
        end
    end
    
    % Add the row Example: ... just if it is present
    if ~isempty(listInpArgs{i,5})
        example=[sprintf(['	<p class="description_valueexample">\r'...
            '       <strong>Example: </strong>']) listInpArgs{i,5} sprintf('</p>\r')];
    else
        example='';
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
        ' <p>']) listInpArgs{i,4} sprintf('</p>\r') ...
        example ...
        sprintf([' <p class="datatypelist"><strong>\r'...
        ' Data Types: </strong><code>' listInpArgs{i,jins}  '</code></p>\r'...
        ' </div>\r'...
        ' </div>\r'])];
    
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
    listOptArgs='';
else
    codewithexample='';
    for i=1:size(listOptArgs,1)
        
        NamVali=listOptArgs{i,5};
        if isempty(NamVali)
            error('FSDA:missingex',['Optional input argument  ' listOptArgs{i,1} ...
                ' does not seem to contain an example (or alternatively string remark has not been put at the end)'])
        end
        % Add as example only those which do finish with </code>, that is
        % just those which do not contain explanations
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
        '<h3 id="NameValuePairs" class="bottom_ruled">\r'...
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
    for i=1:size(listOptArgs,1)
        nameoptarg=listOptArgs{i,1};
        
        titloptarg=listOptArgs{i,2};
        
        shortdesc=listOptArgs{i,3};
        % Remove carriage return if they are present in string shortdesc
        shortdesc=(regexprep(shortdesc,'\r\n|\n|\r',''));
        % replace string ' or ' with ' |'
        shortdesc=strrep(shortdesc,' or ',' | ');
        % Write in lower case first letter of shortdesc
        if length(shortdesc)>1
            shortdesc=[lower(shortdesc(1)) shortdesc(2:end)];
        end
        
        % Check if optional input argument is a structure
        % find just "tructure" and not "structure" or "Structure" because
        % the search is case sensitive
        if ~isempty(strfind(shortdesc,'tructure')) && ~isempty(strfind(listOptArgs{i,4},'field'))
            longdesc=listOptArgs{i,4};
            
            %    [inistructfield,finstructfield]=regexp(longdesc,'\s{8,18}\w*\s{0,8}=');
            [longdescription,listStructureArgs]=formatHTMLstructure(longdesc,nameoptarg);
            listOptArgs{i,8}=listStructureArgs;
            
        else
            longdescriptionHTML=formatHTMLwithMATHJAX(listOptArgs{i,4});
            
            longdescription=longdescriptionHTML;
        end
        % datatype = type of data for that particular option
        %     examplecode=['''Display'',''final'''];
        %     datatype='char';
        
        examplecode=listOptArgs{i,5};
        datatype=listOptArgs{i,6};
        
        OptArgsNameValue=[OptArgsNameValue sprintf(['<div class="expandableContent">\r'...
            '<div id="inputarg_' listOptArgs{i,1} '" class="clearfix">\r'...
            '</div>\r'...
            '<h3 id="input_argument_namevalue_' listOptArgs{i,1} '" class="expand">\r'...
            '<span>\r'...
            '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
            '<span class="argument_name"><code>' nameoptarg  '</code> \r'...
            '&#8212;']) titloptarg sprintf('</span></a><span class="example_desc">') shortdesc  sprintf(['</span></span></h3>\r'...
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
    '<h2 id="OutputArgs" >Output Arguments</h2>\r'...
    '<div class="expandableContent">\r'...
    '<div class="arguments">\r'...
    '<div class="output_argument_container">\r'...
    '<p class="switch">\r'...
    '<a class="expandAllLink" href="javascript:void(0);">expand all</a></p>']);

% outargs = strings which contains output arguments (including [])
% nargout = number of output arguments
outargshtml='';

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
if ~isempty(listargouts)
    listOutArgs(:,1)=listargouts;
end
endpoint=1;
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
            error('FSDA:publishFS:missingOuts',['Output argument ' listargouts{i} ' has not been found'])
        end
        
        % Just in case inipoint has more than one element take just the
        % first but before make sure that you take the first among those
        % which have not been considered yet (that is among those which are
        % greater than endpoint). Suppose for example that first output
        % argument is called Abk and the second output argument is called
        % bk, you want to avoid that string 'bk :' is found inside 'Abk   :'. 
        % Please notice that endpoint has been initialized with number  at the beginning of this loop.
        % Also, please notice that inipoint can have more than element
        % when, for example, say output is out and then
        % string out(:,1) ... is found. 
        inipoint=inipoint(inipoint>=endpoint(1));
        inipoint=inipoint(1);
        
        % The endpoint of the substring is 'more About'. or See also or the next output argument
        if i <nargout-1
            % Note that the endpoint must be searched from position
            % inipoint+length(listargouts{i}) of fstringsel because in
            % order to avaoid cases in which the first output argument is
            % for example ABk and the second output argument is Bk
            endpoint=inipoint+length(listargouts{i})-1+regexp(fstringsel(inipoint+length(listargouts{i}):end),[listargouts{i+1} '\s{0,7}:']);
        elseif i==nargout-1
            
            if strcmp(listargouts{end},'varargout') ==0
                % Note that also in this case the endpoint must be searched from position
                % inipoint+length(listargouts{i}) of fstringsel
                endpoint=inipoint+length(listargouts{i})-1+regexp(fstringsel(inipoint+length(listargouts{i}):end),[listargouts{i+1} '\s{0,7}:']);
                % endpoint=regexp(fstringsel,[listargouts{i+1} '\s{0,7}:']);
                if isempty(endpoint)
                    % warning('FSDA:wrongOutDescription',)
                    errmsg=['Error in processing output argument  ''' listargouts{i} '''\n' ...
                        'Parser could not find string:  ''' listargouts{i+1} '       :''\n' ...
                        'Endpoint for the description of output argument ''' listargouts{i} '''not found'];
                    error('FSDA:wrongOutDescription',errmsg)
                end
                % just in case endpoint has more than one element
                endpoint=endpoint(1);
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
                MoreAboutHTML=formatHTMLwithMATHJAX(MoreAbout);
                MoreAboutHTML=formatHTMLwithList(MoreAboutHTML);
                
            else
                MoreAboutHTML='';
                inipointMoreAbout=Inf;
            end
            
            endpoint=min(inipointSeeAlso,inipointMoreAbout);
            %             if isempty(endpoint)
            %                error('Strings See also and More about are not found inside the file')
            %             end
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
            
            [descrioutput,listStructureArgs]=formatHTMLstructure(descrioutput,outi);
            listOutArgs{i,5}=listStructureArgs;
            % preamble='A structure containing the following fields:';
            preamble='Structure';
            
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
                '<p>']) formatHTMLwithMATHJAX(descrioutput) sprintf(['</p>\r'...
                '</div>\r'...
                '</div>'])];
            
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
        MoreAboutHTML=formatHTMLwithMATHJAX(MoreAbout);
        descrioutput=fstringsel(9:inipointMoreAbout-1);
        
    else
        MoreAbout='';
        MoreAboutHTML='';
        descrioutput=fstringsel(9:inipointSeeAlso-1);
    end
    posPercentageSigns=regexp(descrioutput,'%');
    descrioutput(posPercentageSigns)=[];
    % Remove from string descri leading and trailing white spaces
    descrioutput=strtrim(descrioutput);
    
    outargshtml=formatHTMLwithMATHJAX(descrioutput);
    
end

closeoutargs=sprintf(['</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>']);


outargs=[inioutargs outargshtml closeoutargs];


%% CREATE MORE ABOUT SECTION

if ~isempty(MoreAboutHTML)
    MoreaboutHTMLwithdiv=[sprintf(['<div class="moreabout ref_sect">\r'...
        '<h2 id="MoreAbout">More About</h2>\r'...
        '<div class="expandableContent">\r'...
        '<p class="switch">\r'...
        '<a class="expandAllLink" href="javascript:void(0);">\r'...
        'expand all</a></p>\r'...
        '<div class="expandableContent" itemprop="content">\r'...
        '<h3 class="expand"><span>\r'...
        '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        '<span>Additional Details</span></a></span></h3>\r'...
        '<div class="expand">\r'...
        '<p>']) MoreAboutHTML sprintf(['</p>\r'...
        '</div>\r'...
        '</div>\r'...
        '</div>\r'...
        '</div>'])];
else
    MoreAbout='';
    MoreaboutHTMLwithdiv='';
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
    References='';
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
    
    
    Referenceshtml='';
    iniReferences=sprintf(['<div class="ref_sect" itemprop="content">\r'...
        '<div class="bibliography">\r'...
        '<h2 id="References">References</h2> \r']);
    
    for i=1:length(refsargs)
        Referenceshtml=sprintf([Referenceshtml  '<div><p>' formatHTMLwithMATHJAX(refsargs{i}) '</p></div>\r']);
    end
    Referencesclose=sprintf(['</div>\r'...
        '</div>']);
    References=[iniReferences Referenceshtml Referencesclose];
end



%% ACKNOWLEDGEMENTS
if ~isempty(Acknowledgements)
    iniAcknowledgements=sprintf(['<div class="ref_sect" itemprop="content">\r'...
        '<div class="bibliography">\r'...
        '<h2 id="Acknowledgements">Acknowledgements</h2> \r']);
    
    Acknowledgementshtml=sprintf(['<div><p>' formatHTMLwithMATHJAX(Acknowledgements) '</p></div>\r']);
    
    Ack=[iniAcknowledgements Acknowledgementshtml Referencesclose];
else
    Ack='';
end


%% SEE ALSO

iniSeealso=sprintf(['<div class="ref_sect">\r'...
    '<h2 id="SeeAlso">See Also</h2>\r'...
    '<p>\r']);

% See also:'\:\s*\r

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

Seealsohtml='';
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
                DestHyperLink=[Seealsoitem '.html'];
            else % reference is towards a MATLAB function or a function of another toolbox
                pathdocroot=docroot;
                % Find path of .html documentation file
                pathExtHelpFile=findFile(pathdocroot,'InclFiles',[Seealsoitem '.html']);
                
                if matlabversion==0
                    if isempty(pathExtHelpFile)
                        error('FSDA:publishFS:WrngSeeAlso',['cannot find a reference to doc file ' Seealsoitem '.html']);
                    end
                    pathExtHelpFile=char(pathExtHelpFile{1});
                    addSubPath=pathExtHelpFile(length(pathdocroot)+2:end);
                    
                    % replace '\' with '/'
                    addSubPath=strrep(addSubPath,'\','/') ;
                    % DestHyperLink=['matlab:web(fullfile(docroot,''' addSubPath '.html''))'];
                    DestHyperLink=['matlab:web(fullfile(docroot,''' addSubPath '''))'];
                else
                    DestHyperLink='';
                end
                
            end
        end
        
        Seealsohtml=[Seealsohtml sprintf(['<span itemprop="seealso">\r'...
            '<a href="' DestHyperLink '" itemprop="url">\r'...
            '<span itemprop="name"><code>' Seealsoitem '</code></span></a></span>\r'])];
    end
    
    if i < nseealso
        Seealsohtml=[Seealsohtml ' | '];
    end
end

closeSeealso=sprintf('</div>\r');

Seealso=[iniSeealso Seealsohtml closeSeealso];
% Seealso='';

%% CLOSE TAGS SECTION

clos=sprintf(['<i>This page has been automatically generated by our routine <a href="publishFS.html">publishFS</a></i>\r'...
    '</div>\r'...
    '</section>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r']);


% insbarra=sprintf(['<script type ="text/javascript" language="javaScript">\r'...
%     'document.write(barra);\r'...
%     '</script>']);
insbarra=sprintf('<div id="divbottom"></div>');


closbody=sprintf(['</body>\r'...
    '</html>']);

fclose(fileID);


% site5 is not present anymore in 2015b
%     '<link href="includesFS/site5.css" rel="stylesheet" type="text/css">\r'...


%% Create output structure
out=struct;
% save title
out.titl=name;
% save purpose
out.purpose=purpose;
% Save description
out.description=gendesc;
% Save input arguments (required +optional)
listInpArgs(:,4)=listInpArgs(:,7);
listInpArgs(:,7)={'0'};
listInpArgs(1:nREQargin,7)={'1'};
out.InpArgs=listInpArgs;

% Save optional input arguments of the kind name/pairs
if size(listOptArgs,2)==8
    listOptArgs(:,5)=listOptArgs(:,7);
    listOptArgs=listOptArgs(:,[1:6 8]);
end
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

%% CREATE ON THIS PAGE SECTION WHICH WILL APPEAR IN THE LEFT PANEL
OnThisPageini=['<div id="doc_header_spacer" class="header"></div>\r' ...
'<div class="sticky_header_container includes_subnav">\r'...
  '<div class="section_header level_3">\r'...
    '<div class="container-fluid">\r'...
      '<div class="row" id="mobile_search_row">\r'...
        '<div class="col-xs-12 col-sm-6 col-sm-push-6 col-md-5 col-md-push-7" id="mobile_search">\r'...
          '<div class="search_nested_content_container">\r'...
            '<form id="docsearch_form" name="docsearch_form" method="get" data-release="R2017a" data-language="en" action="../templates/searchresults.html">\r'...
              '<div class="input-group tokenized_search_field">\r'...
                '<label class="sr-only">Search Help</label>\r'...
                '<input type="text" class="form-control conjoined_search" autocomplete="off" name="qdoc" placeholder="Search Help" id="docsearch">\r'...
                '<div class="input-group-btn">\r'...
                  '<button type="submit" name="submitsearch" id="submitsearch" class="btn icon-search btn_search_adjacent btn_search icon_16" tabindex="-1"></button>\r'...
                '</div>\r'...
              '</div>\r'...
            '</form>\r'...
          '</div>\r'...
          '<button class="btn icon-remove btn_search pull-right icon_32 visible-xs" data-toggle="collapse" href="#mobile_search" aria-expanded="false" aria-controls="mobile_search"></button>\r'...
        '</div>\r'...
        '<div class="col-sm-6 col-sm-pull-6 col-md-7 col-md-pull-5" id="section_header_title">\r'...
         '<div class="section_header_content">\r'...
          '<div class="section_header_title">\r'...
              '<h1><a href="../documentation-center.html">Documentation</a></h1>\r'...
            '</div>\r'...
          '</div>\r'...
        '</div>\r'...
        '<div class="visible-xs" id="search_actuator">\r'...
          '<button class="btn icon-search btn_search pull-right icon_16" data-toggle="collapse" href="#mobile_search" aria-expanded="false" aria-controls="mobile_search"></button>\r'...
        '</div>\r'...
      '</div>\r'...
      '<!--END.CLASS row--> \r'...
    '</div>\r'...
    '<!--END.CLASS container-fluid--> \r'...
  '</div>\r'...
  '<!--END.CLASS section_header level_3-->\r'...
  '<div class="horizontal_nav">\r'...
    '<div class="horizontal_nav_container">\r'...
      '<div class="offcanvas_actuator" data-toggle="offcanvas" data-target="#sidebar" id="nav_toggle">\r'...
        '<button type="button" class="btn"><span class="sr-only">Toggle navigation</span><span class="icon-menu icon_24"></span></button>\r'...
        '<span class="offcanvas_actuator_label"></span><span class="offcanvas_actuator_close"></span> </div>\r'...
      '<div class="offcanvas_horizontal_nav">\r'...
        '<div class="container-fluid">\r'...
          '<div class="row">\r'...
            '<div class="col-md-12 hidden-xs hidden-sm"></div>\r'...
          '</div>\r'...
        '</div>\r'...
      '</div>\r'...
    '</div>\r'...
  '</div>\r'...
'</div>\r'...
'<!--END.CLASS sticky_header_container-->\r'...
'<div class="row-offcanvas row-offcanvas-left">\r'...
'<div id="sidebar" class="sidebar-offcanvas" role="navigation">\r'...
 '<nav class="offcanvas_nav" role="navigation">\r'...
    '<ul class="nav_breadcrumb" xmlns:atict="http://www.arbortext.com/namespace/atict">\r'...
      '<li itemprop="breadcrumb" itemscope="" itemtype="http://www.data-vocabulary.org/Breadcrumb"> <a href="../documentation-center.html" itemprop="url"> <span itemprop="title">All Products</span></a></li>\r'...
    '</ul>\r'...
    '<ul class="nav_disambiguation" xmlns:atict="http://www.arbortext.com/namespace/atict">\r'...
      '<li class="product"><a href="index.html">Flexible Statistics and Data Analysis (FSDA)</a>\r'...
        '<div class="dropdown"> <span id="dropdownMenu1" class="icon-arrow-down icon_16" data-toggle="dropdown"> </span>\r'...
          '<ul aria-labelledby="dropdownMenu1" class="dropdown-menu dropdown-menu-right" role="menu">\r'...
            '<li role="presentation"> <a href="examples.html" role="menuitem" tabindex="-1">Examples</a></li>\r'...
            '<li role="presentation"> <a href="function-alpha.html" role="menuitem" tabindex="-1"> Functions and Other Reference</a></li>\r'...
            '<li role="presentation"> <a href="release-notes.html" role="menuitem" tabindex="-1"> Release Notes</a></li>\r'...
            '<li role="presentation"> <a href="tutorials.html" role="menuitem" tabindex="-1"> Tutorials</a></li>\r'...
          '</ul>\r'...
        '</div>\r'...
      '</li>\r'...
    '</ul>\r'...
    '<ul xmlns:atict="http://www.arbortext.com/namespace/atict" class="nav_disambiguation">\r'...
    ' <li><a href="index.html">Flexible Statistics and Data Analysis Toolbox</a> </li>\r'...
     ' <li><a href="function-cate.html" itemprop="url">Functions</a></li>\r'...
    '</ul>\r'...
    '<!-- -->\r'...
    '<ul class="nav_scrollspy nav" xmlns:atict="http://www.arbortext.com/namespace/atict">\r'...
    '  <li class="nav_scrollspy_function"> <a href="#responsive_offcanvas">' name '</a></li>\r'...
     ' <li id="SSPY810-refentry" class="nav_scrollspy_title">On this page</li>\r'...
      '<li><a class="intrnllnk" href="#syntax">Syntax</a></li>\r'...
      '<li><a class="intrnllnk" href="#Description">Description</a></li>'];



if ~isempty(listEx)
    Ini=['<li><a class="intrnllnk" href="#Examples">Examples</a>\r'...
        '<ul>\r'];
    Core='';
    for i=1:size(listEx,1)
        Core=[Core '<li><a class="intrnllnk" href="#Example_' num2str(i) '">' listEx{i,1} '</a></li>\r'];
    end
    Fin=['</ul>\r'...
        '</li>'];
    Core=strrep(Core,'\h','\\h');
    Core=strrep(Core,'\(','\\(');
    Core=strrep(Core,'\)','\\)');
    Core=strrep(Core,'\b','\\b');
    Core=strrep(Core,'\lam','\\lam');
    
    OnThisPageExamples=([Ini Core Fin]);
else
    OnThisPageExamples='';
end

if ~isempty(listExtraEx)
    Ini=['<li><a class="intrnllnk" href="#ExtraExamples">Extra Examples</a>\r'...
        '<ul>\r'];
    Core='';
    for i=1:size(listExtraEx,1)
        Core=[Core '<li><a class="intrnllnk" href="#ExtraExample_' num2str(i) '">' listExtraEx{i,1} '</a></li>\r'];
    end
    Fin=['</ul>\r'...
        '</li>'];
    Core=strrep(Core,'\h','\\h');
    Core=strrep(Core,'\(','\\(');
    Core=strrep(Core,'\)','\\)');
    Core=strrep(Core,'\b','\\b');
    Core=strrep(Core,'\lam','\\lam');
    
    
    OnThisPageExtraExamples=([Ini Core Fin]);
else
    OnThisPageExtraExamples='';
end

if ~isempty(listInpArgs)
    Ini=['<li><a class="intrnllnk" href="#Inputs">Input Arguments</a>\r'...
        '<ul>\r'];
    Core='';
    for i=1:size(listInpArgs,1)
        Core=[Core '<li><a class="intrnllnk" href="#input_argument_' listInpArgs{i,1} '">' listInpArgs{i,1} '</a></li>\r'];
    end
    Fin=['</ul>\r'...
        '</li>'];
    
    OnThisPageInputArguments=([Ini Core Fin]);
else
    OnThisPageInputArguments='';
end

if ~isempty(listOptArgs)
    Ini=['<li><a class="intrnllnk" href="#NameValuePairs">Name-Value Pair Arguments</a>\r'...
        '<ul>\r'];
    Core='';
    for i=1:size(listOptArgs,1)
        Core=[Core '<li><a class="intrnllnk" href="#input_argument_namevalue_' listOptArgs{i,1} '">' listOptArgs{i,1} '</a></li>\r'];
    end
    Fin=['</ul>\r'...
        '</li>'];
    
    OnThisPageNameValuePairs=([Ini Core Fin]);
else
    OnThisPageNameValuePairs='';
end


if ~isempty(listOutArgs)
    Ini=['<li><a class="intrnllnk" href="#OutputArgs">Output Arguments</a>\r'...
        '<ul>\r'];
    Core='';
    for i=1:size(listOutArgs,1)
        Core=[Core '<li><a class="intrnllnk" href="#output_argument_' listOutArgs{i,1} '">' listOutArgs{i,1} '</a></li>\r'];
    end
    Fin=['</ul>\r'...
        '</li>'];
    OnThisPageOutputArgs=([Ini Core Fin]);
else
    OnThisPageOutputArgs='';
end

if ~isempty(MoreAbout)
    OnThisPageMoreAbout=['<li><a class="intrnllnk" href="#MoreAbout">More About</a>\r'...
        '</li>'];
else
    OnThisPageMoreAbout='';
end


if ~isempty(References)
    OnThisPageReferences=['<li><a class="intrnllnk" href="#References">References</a>\r'...
        '</li>'];
else
    OnThisPageReferences='';
end

if ~isempty(Acknowledgements)
    OnThisPageAcknowledgements=['<li><a class="intrnllnk" href="#Acknowledgements">Acknowledgements</a>\r'...
        '</li>'];
else
    OnThisPageAcknowledgements='';
end


OnThisPageSeeAlso=['<li><a class="intrnllnk" href="#SeeAlso">See Also</a>\r'...
    '</li>'];


OnThisPagefin=['</ul>\r'...
    '</nav>\r'...
    '<script src="includesFS/offcanvas.js"></script>\r'...
    '</div>'];
% <!--END.CLASS sidebar-offcanvas-->
new2015b=[OnThisPageini OnThisPageExamples OnThisPageExtraExamples OnThisPageInputArguments ...
    OnThisPageNameValuePairs OnThisPageOutputArgs OnThisPageMoreAbout  ...
    OnThisPageAcknowledgements OnThisPageReferences OnThisPageSeeAlso OnThisPagefin];
metacontent2015b=[metacontent sprintf(new2015b)];

%% Write all pieces in a HTML file
outstring=([titl metacontent2015b sitecont sintaxhtml sintaxclose description  ....
    examples InputArgs outargs MoreaboutHTMLwithdiv References Ack Seealso clos  insbarra closbody]);
%insnav before insbarra has been deleted
%insnav


if write2file
    file1ID=fopen([outputDir fsep name '.html'],'w');
    
    if file1ID==-1
        
        if ismac || isunix
            errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
        elseif ispc
            outputDir=strrep(outputDir,'\','\\');
            errmsg= [' Path ' outputDir '\\' name '.html does not exist or output file '  name '.html is not writable'];
        else
            errmsg= [' Path ' outputDir '/' name '.html does not exist or output file '  name '.html is not writable'];
        end
        
        error('FSDA:publishFS:WrngOutFolder',errmsg);
        
    end
end

%% EXECUTE THE EXAMPLES WHICH START WITH SYMBOLS %%
if evalCode==true
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
                ExToExec=[ExToExec '%% ExExtra' num2str(i) listExtraEx{i,3}];
                numextraexToExec=numextraexToExec+1;
            end
        end
    end
    
    if numextraexToExec+numexToExec>0
        % tmp .file containing all the .m examples will be created. It will be
        % created in subfolder tmp of helpfiles and then automatically removed.
        % This subfolder will be added to put for this session
        nametmp=[name 'tmp.m'];
        % nametmp=[name '.m'];
        
        fullPathToScript=[outputDir fsep 'tmp' fsep nametmp];
        % fullPathToScript=[pathstr fsep 'helpfiles' fsep 'FSDA' fsep 'tmp' fsep nametmp];
        
        
        filetmp=fopen(fullPathToScript,'w');
        addpath([outputDir fsep 'tmp'])
        % addpath([pathstr fsep 'helpfiles' fsep 'FSDA' fsep 'tmp'])
        addpath([pathFSDAstr fsep 'utilities' fsep 'privateFS'])
        
        %        addpath([pathstr fileseparator '\helpfiles\FSDA\tmp'])
        %        addpath([pathstr '\utilities\privateFS'])
        
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
        prefix=name;
        
        
        % file='C:\Users\MarcoAW\D\matlab\FSDA\examples\tmp.m';
        [dom,laste] = evalmxdom(fullPathToScript,dom,cellBoundaries,prefix,imagesDir,outputDir,options);
        %
        drawnow;
        
        ext='html';
        AbsoluteFilename = fullfile(outputDir,[prefix '.' ext]);
        [xResultURI]=xslt(dom,options.stylesheet,AbsoluteFilename);
        drawnow;
        
        % Now remove the temporary .m file with the examples which had been created
        delete(fullPathToScript)
        
        % load html output in a string and extract the parts which are required
        if ismac || isunix
            fileHTML = fopen(xResultURI(6:end), 'r+');
        elseif ispc
            fileHTML = fopen(xResultURI(7:end), 'r+');
        else
            fileHTML = fopen(xResultURI(6:end), 'r+');
            disp('Cannot recognize platform: I use unix as default')
        end
        
        % Insert the file into fstring
        fstringHTML=fscanf(fileHTML,'%c');
        
        totex=numexToExec+numextraexToExec;
        texttoadd=cell(totex,1);
        
        fHTML=regexp(fstringHTML,'<h2>Ex');
        if isempty(fHTML)
            fHTML=regexp(fstringHTML,'<pre class="codeoutput">','once');
        end
        % If fHTML is still empty it means that the ouptut only generates images
        if isempty(fHTML)
            fHTML=regexp(fstringHTML,'<img vspace','once');
        end
        
        % if fHTML is still empty search codeinput
        if isempty(fHTML)
            fHTML=regexp(fstringHTML,'<pre class="codeinput">','once');
            %  if fHTML is still empty produce an error
            if isempty(fHTML)
                errmsg='Parser could not run an example';
                error('FSDA:publishFS:WrngOutFolder',errmsg)
            end
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
        
        ij=1;
        
        for i=1:length(sel)
            % Process string listEx{i,1}
            listExi=listEx{sel(i),1};
            
            % Add symbol \ before special characters in string listExi
            % otherwise regexp will not  find listExi inside
            % outstring
            listExi=SpecialCharacters(listExi);
            
            
            iniout=regexp(outstring,listExi);
            
            
            if length(iniout)<2
                errmsg= [' Title of example \n''' listExi '''\n could not be found \n'...
                    'Probably because the string contains special characters\n' ...
                    'which cannot be interpreted by MATLAB function regexp'];
                error('FSDA:publishFS:WrngOutFolder',errmsg)
            end
            
            finout=regexp(outstring,'</pre>');
            % finout=finout(finout>iniout(2));
            
            % Make sure that string '</pre>' which has been found is after
            % the i-th instance of M.gif
            posextoinclude=regexp(outstring,'M\.gif');
            finout=finout(finout>posextoinclude(ij));
            ij=ij+1;
            
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
                
                % Add symbol \ in before special characters in string listExi
                % otherwise regexp will not  find listExi inside
                % outstring
                listExi=SpecialCharacters(listExi);
                
                iniout=regexp(outstring,listExi);
                if length(iniout)>2
                    disp(['Duplicate name for: ' listExi ' found'])
                    warning('FSDA:WrongArg','There are examples with the same title, please use a title which is unique')
                    warning('FSDA:WrongArg','Or alternatively there is an empty line before the start of the example')
                elseif isempty(iniout)
                    errmsg= [' Title of example \n''' listExi '''\n could not be found \n'...
                        'Probably because the string contains special characters\n' ...
                        'which cannot be interpreted by MATLAB function regexp'];
                    error('FSDA:publishFS:WrngEx',errmsg)
                else
                end
                
                % iniout=iniout(1);
                
                finout=regexp(outstring,'</pre>');
                % finout=finout(finout>iniout);
                
                posextoinclude=regexp(outstring,'M\.gif');
                finout=finout(finout>posextoinclude(ij));
                ij=ij+1;
                
                % outstring(finout:finout+11)
                % inclplint = point where output of the example must be included
                inclpoint=finout(1)+18;
                % incl= string which contains the output of the code
                incl=texttoadd{i+numexToExec};
                outstring=[outstring(1:inclpoint) incl outstring(inclpoint+1:end)];
            end
        end
        
        close all
        
        % Remove folder which had temporarily added to path
        %         rmpath([pathstr '\helpfiles\FSDA\tmp'])
        %         rmpath([pathstr '\utilities\privateFS'])
        
        %rmpath([pathstr fsep 'helpfiles' fsep 'FSDA' fsep 'tmp'])
        
        
        rmpath([outputDir fsep 'tmp'])
        rmpath([pathFSDAstr fsep 'utilities' fsep 'privateFS'])
    else
        laste='';
    end
else
    laste='';
end

out.laste=laste;

%% WRITE string outstring into final HTML file
if write2file
    fprintf(file1ID,'%s',outstring);
    fclose('all');
end

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


% descrlongHTMLwithref = this inner function  has the purpose of adding symbols </p> <p> every
% time a full stop, colon or semicolon symbol followed by a series of space
% and then a carriage return.
function descrHTTPwithref=formatHTML(descrlong)
% Find all lines which terminate with full stop
newlinewithFullStop=regexp(descrlong,'\.\s*\r');
% Do not consider the line which terminate with full stop but are preeceed
% by "one space" and then "pp."
[ppwithfullstop]=regexp(descrlong,'\spp\.\s*\r')+3;
newlinewithFullStop=setdiff(newlinewithFullStop,ppwithfullstop);

% Find all lines which terminate with symbol : or with symbol ;
newlinewithColon=regexp(descrlong,'\:\s*\r');
newlinewithSemiColon=regexp(descrlong,'\;\s*\r');
newl=sort([newlinewithColon newlinewithSemiColon newlinewithFullStop]);
if ~isempty(newl)
    descrlongHTML=['<p>' descrlong(1:newl(1))];
    if length(newl)==1
        descrlongHTML=[descrlongHTML '</p> <p>' descrlong(newl(1)+1:end)];
    else
        for j=1:(length(newl)-1)
            descrlongHTML=[descrlongHTML '</p> <p>' descrlong(newl(j)+1:newl(j+1))];
        end
        descrlongHTML=[descrlongHTML '</p><p>' descrlong(newl(j+1)+1:end)];
    end
    descrlongHTML=[descrlongHTML '</p>'];
else
    descrlongHTML=descrlong;
end
% put a hypertext link to all words which end with .m
[IniRefFilem,FinRefFilem]=regexp(descrlongHTML,'\w*\.m[\s\.,]');
% Make sure that the .m string is not standalone, that is make sure
% that the .m string is preceeded by some characters.
boo=(FinRefFilem-IniRefFilem)>2;
IniRefFilem=IniRefFilem(boo);
FinRefFilem=FinRefFilem(boo);

if ~isempty(IniRefFilem)
    FinRefFilem=FinRefFilem-1;
    descrlongHTMLwithref='';
    for i=1:length(IniRefFilem)
        
        namewithoutHTML=descrlongHTML(IniRefFilem(i):FinRefFilem(i)-2);
        namewithHTML=[namewithoutHTML '.html'];
        if i==1 && length(IniRefFilem)==1
            descrlongHTMLwithref=[descrlongHTMLwithref descrlongHTML(1:IniRefFilem(i)-1) ...
                '<a href="' namewithHTML '">' namewithoutHTML '</a>'...
                descrlongHTML(FinRefFilem(i)+1:end)];
        elseif i==1
            descrlongHTMLwithref=[descrlongHTMLwithref descrlongHTML(1:IniRefFilem(i)-1) ...
                '<a href="' namewithHTML '">' namewithoutHTML '</a>'];
            
        elseif i==length(IniRefFilem)
            descrlongHTMLwithref=[descrlongHTMLwithref descrlongHTML(FinRefFilem(i-1)+1:IniRefFilem(i)-1) ...
                '<a href="' namewithHTML '">' namewithoutHTML '</a>'...
                descrlongHTML(FinRefFilem(i)+1:end)];
        else
            descrlongHTMLwithref=[descrlongHTMLwithref descrlongHTML(FinRefFilem(i-1)+1:IniRefFilem(i)-1) ...
                '<a href="' namewithHTML '">' namewithoutHTML '</a>'];
        end
        
        %
    end
else
    descrlongHTMLwithref=descrlongHTML;
end

% put a hypertext link to all words which start with http
descrHTTP=descrlongHTMLwithref;
[IniRefhttp]=regexp(descrHTTP,'[^"]http');

if ~isempty(IniRefhttp)
    FinRefhttp=zeros(length(IniRefhttp),1);
    descrHTTPwithref='';
    for i=1:length(IniRefhttp)
        
        descrHTTPsel=descrHTTP(IniRefhttp(i)+1:end);
        findfirstspace=regexp(descrHTTPsel,'\s','once');
        if isempty(findfirstspace)
            FinRefhttp(i)=IniRefhttp(i)+length(descrHTTPsel)+1;
        else
            FinRefhttp(i)=IniRefhttp(i)+findfirstspace;
        end
        
        
        namehttp=strtrim(descrHTTP(IniRefhttp(i):FinRefhttp(i)-1));
        % Make sure you do not select string </p> at the end
        if strcmp(namehttp(end-3:end),'</p>')
            namehttp=namehttp(1:end-4);
        end
        
        
        if strcmp(namehttp(end),'.') || strcmp(namehttp(end),',')
            namehttp=namehttp(1:end-1);
        end
        
        FinRefhttp(i)=IniRefhttp(i)+length(namehttp);
        
        if i==1 && length(IniRefhttp)==1
            descrHTTPwithref=[descrHTTPwithref descrHTTP(1:IniRefhttp(i)-1) ...
                '<a href="' namehttp '">' namehttp '</a>'...
                descrHTTP(FinRefhttp(i)+1:end)];
        elseif i==1
            descrHTTPwithref=[descrHTTPwithref descrHTTP(1:IniRefhttp(i)-1) ...
                '<a href="' namehttp '">' namehttp '</a>'];
            
        elseif i==length(IniRefhttp)
            descrHTTPwithref=[descrHTTPwithref descrHTTP(FinRefhttp(i-1)+1:IniRefhttp(i)-1) ...
                '<a href="' namehttp '">' namehttp '</a>'...
                descrHTTP(FinRefhttp(i):end)];
        else
            descrHTTPwithref=[descrHTTPwithref descrlongHTML(FinRefhttp(i-1)+1:IniRefhttp(i)-1) ...
                '<a href="' namehttp '">' namehttp '</a>'];
        end
    end
else
    descrHTTPwithref=descrHTTP;
end
end

function StringHTML=formatHTMLwithMATHJAX(inputString)

% Check if symbols \[ \] are present
% If this is the case it is necessary to split inputString into
% the text_part and the Mathjax_part and apply HTML format just
% to the complementary of the MathJax part
iniMathJax=regexp(inputString,'\\\[');
finMathJax=regexp(inputString,'\\\]');



if ~isempty(iniMathJax) || ~isempty(finMathJax)
    
    if length(iniMathJax) ~= length(finMathJax)
        disp('Latex error in string:')
        disp('---------')
        disp(inputString)
        disp('---------')
        error('FSDA:wrongLatex','There is a non matching math symbol in the LaTeX equation')
    end
    
    MoreA=formatHTML(inputString(1:iniMathJax-1));
    for k=1:length(iniMathJax)
        MoreA=[MoreA inputString(iniMathJax(k):finMathJax(k)+1)];
        if k==length(iniMathJax)
            MoreA=[MoreA formatHTML(inputString(finMathJax(k)+2:end))];
        else
            MoreA=[MoreA formatHTML(inputString(finMathJax(k)+2:iniMathJax(k+1)-1))];
        end
    end
    StringHTML=MoreA;
else
    % In this case there are not latex formulae so just apply
    % routine formatHTML
    StringHTML=formatHTML(inputString);
end

end

function [descrioutput,listStructureArgs]=formatHTMLstructure(descriinput,StructureName)
iniTable=[sprintf(['<table border="2" cellpadding="4" cellspacing="0" class="body">\r'...
    '<colgroup>\r']) ...
    '<col width="20%"><col width="80%">'...
    sprintf(['\r</colgroup>\r'...
    '<thead>\r'...
    '<tr valign="top">\r'...
    '<th valign="top">Value</th>\r'...
    '<th valign="top">Description</th>\r'...
    '</tr>\r'...
    '</thead>'])];
cloTable='</table>';

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

%[~,finchk]=regexp(descriinput,['\n\s*' StructureName '\.\w*\s*=']);
% Select all rows of iniA and finA in which the elements of finA are equal
% to those of finchk
% if length(finA)>length(finchk)
%     [~,ia]=intersect(finA,finchk);
%     ini=iniA(ia);
%     fin=finA(ia);
% else
%     ini=iniA;
%     fin=finA;
% end

if isempty(ini)
    disp('Probably ":" symbols  must be replaced with "=" symbols in out description')
    error('FSDA:MissingArg',['Parser cannot find string \n''' StructureName '.xxxx'' = \n for structure ' StructureName])
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
    
    %         if strcmp(Display,'iter-detailed')
    %             disp('Detailed information about Output arguments')
    %             disp(listOutArgs)
    %         end
    
    % remove from inisel the rows equal to 999 (that is the rows which
    % correspond to duplicated arguments)
    inisel(inisel==999)=[];
    
    Tablehtml='';
    for k=inisel % length(ini)
        
        descrlong=listStructureArgs{k,2};
        
        descrlongHTML=formatHTMLwithMATHJAX(descrlong);
        
        % listOutArgs{k,2}
        Tablehtml=[Tablehtml sprintf(['<tr valign="top">\r'...
            '<td><code>' listStructureArgs{k,1} '</code></td>\r'...
            '<td>\r'...
            '<p>']) descrlongHTML  sprintf(['</p>\r'...
            '</td>\r'...
            '</tr>'])];
    end
    
    preamble=descriinput(1:ini(1)-1);
    preambleHTML=formatHTMLwithMATHJAX(preamble);
    
    % Add the Remark after the table, if it is present
    if ~isempty(posREMARK)
        descrREMARK=descriinput(posREMARK:end);
        descrREMARKHTML=formatHTMLwithMATHJAX(descrREMARK);
        
        descrioutput=[preambleHTML iniTable Tablehtml cloTable '<p>' descrREMARKHTML '</p>'];
    else
        descrioutput=[preambleHTML iniTable Tablehtml cloTable];
    end
    
    % use capital letter for the first word.
    descrioutput=[upper(descrioutput(1)) descrioutput(2:end)];
else
    descrioutput=formatHTMLwithMATHJAX(descriinput);
    listStructureArgs='';
end

end


function StringwithSpecialCharacters=SpecialCharacters(StringwithSpecialCharacters)
% Add symbol \ before special characters in string StringwithSpecialCharacters
% otherwise regexp will not find this string in
% outstring

% If there are signs $ ^ [ ] replace them with \$ and \^ \[ \]
StringwithSpecialCharacters=strrep(StringwithSpecialCharacters,'$','\$');
StringwithSpecialCharacters=strrep(StringwithSpecialCharacters,'^','\^');
StringwithSpecialCharacters=strrep(StringwithSpecialCharacters,'[','\[');
StringwithSpecialCharacters=strrep(StringwithSpecialCharacters,']','\]');
StringwithSpecialCharacters=strrep(StringwithSpecialCharacters,'(','\(');
StringwithSpecialCharacters=strrep(StringwithSpecialCharacters,')','\)');
StringwithSpecialCharacters=strrep(StringwithSpecialCharacters,'.','\.');
end


function descrlongHTML=formatHTMLwithList(descrlong)
newl=regexp(descrlong,'\[\d*\]');
if ~isempty(newl)
    descrlongHTML=[descrlong(1:newl(1)-3) '<hr>[<b>1</b>] '];
    if length(newl)==1
        descrlongHTML=[descrlongHTML descrlong(newl(1)+3:end)];
    else
        for j=1:(length(newl)-1)
            descrlongHTML=[descrlongHTML descrlong(newl(j)+4:newl(j+1)-3) '<hr>[<b>' num2str(j+1) '</b>] ' ];
        end
        descrlongHTML=[descrlongHTML descrlong(newl(j+1)+4:end)];
    end
    descrlongHTML=[descrlongHTML '<hr>'];
else
    descrlongHTML=descrlong;
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

% function descrlongHTML=formatHTMLwithList(descrlong)
% newl=regexp(descrlong,'\[\d*\]');
% if ~isempty(newl)
%     descrlongHTML=[descrlong(1:newl(1)-3) '<ol><li>'];
%     if length(newl)==1
%         descrlongHTML=[descrlongHTML descrlong(newl(1)+3:end)];
%     else
%         for j=1:(length(newl)-1)
%             descrlongHTML=[descrlongHTML descrlong(newl(j)+4:newl(j+1)-3) '</li><li>'];
%         end
%         descrlongHTML=[descrlongHTML descrlong(newl(j+1)+4:end)];
%     end
%     descrlongHTML=[descrlongHTML '</li></ol>'];
% else
%     descrlongHTML=descrlong;
% end
% end

%FScategory:UTIHELP