function [outstring,laste]=htmlwriteFS(IPS,varargin)
%htmlwriteFS is an obsolete function which will be removed in future releases. Use publishFS.m instead.
%
%<a href="matlab: docsearchFS('htmlwriteFS')">Link to the help function</a>
%
%   htmlwriteFS creates HTML files a MATLAB structure created
%   with function mreadFS.m. Note that htmlwriteFS and mreadFS.m are
%   obsolete and will be removed in future releases. USE publishFS.m
%   instead.
%
% Required input arguments:
%
%   IPS:  Specific MATLAB strutcure. Structure. A structure created with
%                   file mreadFS.m or modified with the GUI.
%                   containing the following fields:
% IPS.title     = title of HTML file. String.
%                   String to be included in HTML tag title
% IPS.purpose   = purpose of the routine. String.
%                 String forming second row of output HTML file
%IPS.description= short description of the file. String
%                 If this field is not empty it is included in the HTML file
%                 at the beginning of the section "Description"
% IPS.InpArgs   = Required and Optional input arguments. Cell.
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
% IPS.OptArgs   = Optional input arguments specified as name/values pairs. Cell.
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
% IPS.OutArgs   = Required and Optional (varaargout) output arguments. Cell.
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
% IPS.MoreAbout = More About. String. String containing what in the HTML
%                 file will appear under the section "More About".
%IPS.Acknowledgements = Acknowledgements. String. String containing what in the HTML
%                 file will appear under the section "Acknowledgements".
%IPS.References = References. cell. Cell of length r containing the
%                 references.
%   IPS.SeeAlso = References. cell. Cell of length s containing the
%                 references to linked files.
%        IPS.Ex = Examples. cell. Cell of length t containing the
%                 examples.
%                 First column= title of the example;
%                 Second column = detailed description;
%                 Third column = MATLAB code;
%                 Fourth column = dummy variable which indicates whether
%                 the example must be executed or not) If 1 example is executed
%    IPS.ExtraEx= Extra Examples. cell. Cell of length u containing the u
%                 extra examples.
%                 First column= title of the example;
%                 Second column = detailed description;
%                 Third column = MATLAB code;
%                 Fourth column = dummy variable which indicates whether
%                 the example must be executed or not) If 1 example is executed
%             Data Types - struct
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
%             images will be save in subfolder images_help of the
%             outputDir. The default value of evalCode is true.
%             Example - 'evalCode','false'
%             Data Types - Boolean
% write2file: Option to write HTML file. Logical. Option which specifies
%             whether HTML file must be created or if just string outstring
%             must be created. The default value of write2file is true,
%             that is html file is created
%             Example - 'write2file','false'
%             Data Types - Boolean
%ErrWrngSeeAlso: Option to check links in the see also part. Logical.
%            If ErrWrngSeeAlso is true publishFS checks whether the strings
%            inside see also are valid files and puts an hyperlink to the
%            corresponding HTML file. If publishFS cannot find the files
%            exits the procedure with an error. If ErrWrngSeeAlso is false
%            no check is done and empty links are produced. Use
%            ErrWrngSeeAlso set to false if the purpose is just to check
%            the code (e.g. in external environment like TRAVIS) and not to
%            buid the help system. Default value of ErrWrngSeeAlso is true.
%             Example - 'ErrWrngSeeAlso',false
%             Data Types - logical
%
% Output:
%
% outstring : string wchich contains the processed HRML file. String.
%             String containing parsed html file.
%  laste   : Information about errors. MException class. Object of class
%            MException which provides information about last error in
%            executing the examples. If the procedure runs without errors
%            laste is an empty value;
%
%
% More About:
%  This function uses routine strjoin and can be used just by those which
%  have at least MATLAB R2013a
%
% See also: publishFS
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('htmlwriteFS')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    % htmlwriteFS with all the default options.
    % Create file FSRmdr.html starting from file FSRmdr.
    NameFile='FSRmdr.m';
    out=publishFS(NameFile,'ErrWrngSeeAlso',false,'evalCode',false);
%}

%{
    % htmlwriteFS with display option.
    % Create file FSRmdr.html starting from file FSRmdr and
    % display detailed information about the Input, Output and Optional
    % arguments.
    NameFile='FSRmdr.m';
    out=publishFS(NameFile,'Display','iter-detailed','ErrWrngSeeAlso',false,'evalCode',false);
%}



%% Beginning of code

if ~isstruct(IPS)
    error('FSDA:htmlwriteFS:WrongInput','input must be a MATLAB struct.');
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

ErrWrngSeeAlso=true;

if nargin>1
    options=struct('evalCode',evalCode,'Display',Display,'outputDir',outputDir,...
        'imagesDir',imagesDir,'write2file',true,'ErrWrngSeeAlso',ErrWrngSeeAlso);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:htmlwriteFS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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
    ErrWrngSeeAlso=options.ErrWrngSeeAlso;
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

if strcmp(Display,'iter-detailed')
    disp('Check that file contains appropriate reference to itself inside docsearchFS')
end
%linkHTML must a vector with length equal to 2


%% Extract all relevant items from input structure IPS

name=IPS.titl;
purpose=IPS.purpose;
gendesc=IPS.description;
ErrWrngSeeAlso=IPS.ErrWrngSeeAlso;

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
    '      $("#div002").load("includesFS/top.html");\r'...
    '      $("#div001").load("includesFS/bottom.html");\r'...
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


searchenginestring=['<div class="section_header level_3">\r'...
    '<div class="container-fluid">\r'...
    '<div class="row" id="mobile_search_row">\r'...
    '<div class="col-xs-12 col-sm-6 col-sm-push-6 col-md-5 col-md-push-7" id="mobile_search">\r'...
    '<div class="search_nested_content_container">\r'...
    '<form id="docsearch_form" name="docsearch_form" method="get" data-release="R2016a" data-language="en" action="../templates/searchresults.html">\r'...
    '<div class="input-group tokenized_search_field">\r'...
    '<label class="sr-only">Search Documentation</label>\r'...
    '<input type="text" class="form-control conjoined_search" autocomplete="off" name="qdoc" placeholder="Search Documentation" id="docsearch">\r'...
    '<div class="input-group-btn">\r'...
    '<button type="submit" name="submitsearch" id="submitsearch" class="btn icon-search btn_search_adjacent btn_search icon_16" tabindex="-1"></button>\r'...
    '</div>\r'...
    '</div>\r'...
    '</form>\r'...
    '</div>\r'...
    '<button class="btn icon-remove btn_search pull-right icon_32 visible-xs" data-toggle="collapse" href="#mobile_search" aria-expanded="false" aria-controls="mobile_search"></button></div>\r'...
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
    '</div><!--END.CLASS row-->\r'...
    '</div><!--END.CLASS container-fluid-->\r'...
    '</div>\r'];


metacontent=sprintf([beforemetacontent purpose aftermetacontent searchenginestring]);


insnav=sprintf('<div id="div002"></div>');



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


if length(gendesc)<3
    htmlsitecont='';
else
    htmlsitecont=['<p>' formatHTMLwithMATHJAX(gendesc) '</p>'];
end

finsitecont=sprintf(['<h2 id="syntax">Syntax</h2>\r'...
    '<div class="syntax_signature">\r'...
    '<div class="syntax_signature_module">\r'...
    '<ul>\r']);

sitecont=[insnav inisitecont finsitecont];

%% Create Sintax section of HTML file (referred to the loop of input arguments)
% Find number of required input arguments
nREQargin=sum(str2num(cell2mat(IPS.InpArgs(:,7)))); %#ok<ST2NM>
% Find number of optional input arguments (not varargin)
nTOTargin=size(IPS.InpArgs,1);
nOPTargin=nTOTargin-nREQargin;
if size(IPS.OutArgs,1)>0
    nargout=size(IPS.OutArgs,1);
else
    outargs='';
    nargout=0;
end

OutArgsNames=IPS.OutArgs(:,1);
InputArgs=IPS.InpArgs(:,1);

sintax=cell(nargout+1+nOPTargin,1);

if nOPTargin>0
    [commasIn] = regexp(InputArgs,',');
    for j=1:nOPTargin
        if ~isempty(OutArgsNames)
            sintax{j}=[OutArgsNames{1} '=' name InputArgs{j} ')'];
        else
            if ~isempty(commasIn)
                sintax{j}=[name InputArgs{j} ')'];
            else
                sintax{j}=name;
            end
        end
        
    end
    j=j+1;
else
    j=1;
end

OptArgsVarargin=size(IPS.OptArgs,1);


if OptArgsVarargin==0
    if ~isempty(outargs)
        sintax{j}=[OutArgsNames{j} '=' name InputArgs];
    else
        sintax{j}=[name InputArgs];
    end
    j=j+1;
else
    strinputarg=['(' strjoin(InputArgs,',')];
    
    
    if ~isempty(OutArgsNames)
        sintax{j}=[OutArgsNames{j} '=' name strinputarg ')'];
        if length(strinputarg)>1
            sintax{j+1}=[OutArgsNames{j} '=' name strinputarg ',Name,Value)'];
        else
            % just in case function has no compulsary input argument then
            % the comma before 'Name.value' is unnecessary
            sintax{j+1}=[OutArgsNames{j} '=' name strinputarg 'Name,Value)'];
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
            outargs=['[' strjoin(OutArgsNames(1:i),',') ']'];
            
            sintax{j}=[outargs ']=' name '(___)'];
        else
            outargs=['[' strjoin(OutArgsNames,',') ']'];
            
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


%listEx = list which contains the examples (associated to sintax)
% First column= title, second column detailed description. Third column
% code
% Fourth column is a dummy variable which indicates whether the example must be
% executed or not). If 1 example is executed
listEx=IPS.Ex;

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
% Fourth column of listextraEX contains flag 1 or 0 depending on the
% fact that the example must be executed or not
listExtraEx=IPS.ExtraEx;



if strcmp(Display,'iter-detailed')
    disp('Detailed information about all the Extra examples')
    disp(listExtraEx)
end

% Now check whether the first column of cell listExtraEx contains the string interactive_example
%     for i=1:size(listExtraEx,1)
%         [StartInteractive,EndInteractive]=regexp(listExtraEx{i,1},'[Ii]nteractive_example.');
%         if ~isempty(StartInteractive)
%             StringToReplace=listExtraEx{i,1};
%                 listExtraEx{i,1}=['<i>Interactive example ' num2str(NumOfInterEx)  '.</i>' StringToReplace(EndInteractive+1:end)];
%             NumOfInterEx=NumOfInterEx+1;
%         end
%     end


closeexamples=sprintf(['</div>\r'... % close div id="expandableExamples
    '<p> </p>\r']);
% Related examples are below
iniRelatedExamples='';
RelatedExamples='';
if ~isempty(listExtraEx)
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

listInpArgs=IPS.InpArgs;
for i=1:nTOTargin
    
    % Add the row Example: ... just if it is present
    if ~isempty(listInpArgs{i,5})
        example=[sprintf(['	<p class="description_valueexample">\r'...
            '       <strong>Example: </strong>']) listInpArgs{i,5} sprintf('</p>\r')];
    else
        example='';
    end
    
    if ~isempty(strfind(listInpArgs{i,3},'tructure'))
        FormattedInpArg=formatHTMLtable(listInpArgs{i,8});
        if ~isempty(listInpArgs{i,4})
            FormattedInpArg=[listInpArgs{i,4} '<p>'  FormattedInpArg '</p>'];
        end
    else
        FormattedInpArg=listInpArgs{i,4};
    end
    reqargs=[reqargs sprintf(['<div class="expandableContent">\r'...
        ' <div id="inputarg_' listInpArgs{i,1} '" class="clearfix">\r'...
        ' </div>\r'...
        ' <h3 id="input_argument_' listInpArgs{i,1} '" class="expand">\r'...
        ' <span>\r'...
        ' <a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        ' <span class="argument_name"><code>' listInpArgs{i,1} '</code> &#8212; ']) listInpArgs{i,2} sprintf([' </span> \r'...  % &#8212; = long dash
        ' </a><span class="example_desc">']) listInpArgs{i,3} sprintf(['</span></span></h3>\r'...
        ' <div class="collapse">\r'...
        ' <p>']) FormattedInpArg sprintf('</p>\r') ...
        example ...
        sprintf([' <p class="datatypelist"><strong>\r'...
        ' Data Types: </strong><code>' listInpArgs{i,6}  '</code></p>\r'...
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

listOptArgs=IPS.OptArgs;
if isempty(listOptArgs)
    OptArgsNameValueHeading='';
    OptArgsNameValue='';
else
    codewithexample='';
    for i=1:size(listOptArgs,1)
        
        NamVali=listOptArgs{i,5};
        
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
        if ~isempty(listOptArgs{i,7})
            longdesc=listOptArgs{i,7};
            
            %    [inistructfield,finstructfield]=regexp(longdesc,'\s{8,18}\w*\s{0,8}=');
            [longdescription]=formatHTMLtable(longdesc);
            if ~isempty(listOptArgs{i,4})
                longdescription=[longdescription '<p>' listOptArgs{i,4} '</p>'];
            end
        else
            longdescription=formatHTMLwithMATHJAX(listOptArgs{i,4});
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


listOutArgs=IPS.OutArgs;
listargouts=listOutArgs(:,1);

if ~isempty(listargouts)
    for i=1:nargout
        
        % listargouts is a cell which contains the list of output arguments
        outi=listargouts{i};
        
        
        % The initial point of the string is 'listargouts{i}' is there is just
        % one output else is string 'listargouts{i} :' is there is more than
        % one output and this is not varargout
        % else if there is varargout the initialpoint is the string
        % "Optional Output:"
        
        
        
        
        % descri = string which contains the description of i-th output
        
        
        % Check if the output is a structure. If this is the case
        if ~cellfun(@isempty,listOutArgs(i,5))
            descrioutput=listOutArgs{i,5};
            [descrioutput]=formatHTMLtable(descrioutput);
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
            preamble=listOutArgs{i,2};
            descroutputtitl=listOutArgs{i,3};
            descrioutput=listOutArgs{i,4};
            
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
end

closeoutargs=sprintf(['</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>']);


outargs=[inioutargs outargshtml closeoutargs];


%% CREATE MORE ABOUT SECTION
MoreAbout=IPS.MoreAbout;
MoreAboutHTML=formatHTMLwithMATHJAX(MoreAbout);
MoreAboutHTML=formatHTMLwithList(MoreAboutHTML);

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
        '<span>Additional Details </span></a></span></h3>\r'...
        '<div class="collapse">\r'...
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
Acknowledgements=IPS.Acknowledgements;

refsargs=IPS.References;


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



Seealsohtml='';
% Store in cell listSeeAlso, See also items
listSeeAlso=IPS.SeeAlso;
nseealso=length(listSeeAlso);

for i=1:nseealso
    
    
    % Store See also item
    Seealsoitem=listSeeAlso{i};
    
    str=which(Seealsoitem);
    
    if ErrWrngSeeAlso == true
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
                
                if isempty(pathExtHelpFile)
                    error('FSDA:publishFS:WrngSeeAlso',['cannot find a reference to doc file ' Seealsoitem '.html']);
                end
                pathExtHelpFile=char(pathExtHelpFile{1});
                addSubPath=pathExtHelpFile(length(pathdocroot)+2:end);
                
                % replace '\' with '/'
                addSubPath=strrep(addSubPath,'\','/') ;
                % DestHyperLink=['matlab:web(fullfile(docroot,''' addSubPath '.html''))'];
                DestHyperLink=['matlab:web(fullfile(docroot,''' addSubPath '''))'];
            end
        end
    else
        DestHyperLink='';
    end
    
    
    Seealsohtml=[Seealsohtml sprintf(['<span itemprop="seealso">\r'...
        '<a href="' DestHyperLink '" itemprop="url">\r'...
        '<span itemprop="name"><code>' Seealsoitem '</code></span></a></span>\r'])];
    if i < nseealso
        Seealsohtml=[Seealsohtml ' | '];
    end
end

closeSeealso=sprintf('</div>\r');

Seealso=[iniSeealso Seealsohtml closeSeealso];


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
insbarra=sprintf('<div id="div001"></div>');


closbody=sprintf(['</body>\r'...
    '</html>']);


%% CREATE ON THIS PAGE SECTION WHICH WILL APPEAR IN THE LEFT PANEL
OnThisPageini=['<!--START.CLASS sticky_header_container-->\r'...
    '<div class="sticky_header_container includes_subnav">\r'...
    '<div class="horizontal_nav_container">\r'...
    '<div id="nav_toggle" class="offcanvas_actuator" data-target="#sidebar" data-toggle="offcanvas">\r'...
    '<button class="btn" type="button"><span class="sr-only">Toggle navigation</span><span class="icon-menu icon_32"></span>\r'...
    '</button><span class="offcanvas_actuator_label"></span>\r'...
    '<span class="offcanvas_actuator_close"></span></div>\r'...
    '<div class="offcanvas_horizontal_nav">\r'...
    '<div class="container-fluid">\r'...
    '<div class="row">\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '<div class="row-offcanvas row-offcanvas-left">\r'...
    '<div id="sidebar" class="sidebar-offcanvas" role="navigation">\r'...
    '<nav class="offcanvas_nav" role="navigation">\r'...
    '<ul class="nav_breadcrumb" xmlns:atict="http://www.arbortext.com/namespace/atict">\r'...
    '<li itemprop="breadcrumb" itemscope="" itemtype="http://www.data-vocabulary.org/Breadcrumb">\r'...
    '<a href="../documentation-center.html" itemprop="url">\r'...
    '<span itemprop="title">All Products</span></a></li>\r'...
    '</ul>\r'...
    '<ul class="nav_disambiguation" xmlns:atict="http://www.arbortext.com/namespace/atict">\r'...
    '<li class="product"><a href="index.html">Flexible Statistics and \r'...
    'Data Analysis (FSDA)</a>\r'...
    '<div class="dropdown">\r'...
    '<span id="dropdownMenu1" class="icon-arrow-down icon_16" data-toggle="dropdown">\r'...
    '</span>\r'...
    '<ul aria-labelledby="dropdownMenu1" class="dropdown-menu dropdown-menu-right" role="menu">\r'...
    '<li role="presentation">\r'...
    '<a href="examples.html" role="menuitem" tabindex="-1">Examples</a></li>\r'...
    '<li role="presentation">\r'...
    '<a href="functionlist.html" role="menuitem" tabindex="-1">\r'...
    'Functions and Other Reference</a></li>\r'...
    '<li role="presentation">\r'...
    '<a href="release-notes.html" role="menuitem" tabindex="-1">\r'...
    'Release Notes</a></li>\r'...
    '<li role="presentation">\r'...
    '<a href="tutorials.html" role="menuitem" tabindex="-1">\r'...
    'Tutorials</a></li>\r'...
    '</ul>\r'...
    '</div>\r'...
    '</li>\r'...
    '</ul>\r'...
    '<ul xmlns:atict="http://www.arbortext.com/namespace/atict" class="nav_disambiguation">\r'...
    '<li><a href="index.html">Flexible Statistics and Data Analysis Toolbox</a>\r'...
    '</li>\r'...
    '<li><a href="function-cate.html" itemprop="url">Functions</a></li>\r'...
    '</ul>\r'...
    '<!-- -->\r'...
    '<ul class="nav_scrollspy nav" xmlns:atict="http://www.arbortext.com/namespace/atict">\r'...
    '<li class="nav_scrollspy_function">\r'...
    '<a href="#responsive_offcanvas">' name '</a></li>\r'...
    '<li id="SSPY810-refentry" class="nav_scrollspy_title">On this page</li>\r'...
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
%
%
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
        
        ext='html';
        AbsoluteFilename = fullfile(outputDir,[prefix '.' ext]);
        [xResultURI]=xslt(dom,options.stylesheet,AbsoluteFilename);
        
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


%% WRITE string outstring into final HTML file
if write2file
    fprintf(file1ID,'%s',outstring);
    fclose('all');
end


end


% This inner function  has the purpose of adding symbols </p> <p> every
% time a full stop, colon or semicolo symbol followed by a series of space
% and then a carriage return.
% descrlongHTMLwithref

function descrHTTPwithref=formatHTML(descrlong)
newlinewithFullStop=regexp(descrlong,'\.\s*\r');
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
% If this is the case it is necessary to split inputSring into
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

function [listStructureArgsFormatted]=formatHTMLtable(listStructureArgs)
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




% listStructureArgs = cells which will contains the name/values
% pairs of the structure
% fieldnames will go into the first column of the table
% the content of the field names will fo into the second column of
% the table


% rowtodel = vector which contains the duplicate rows of
% listStruArgs which have to be deleted
inisel=1:size(listStructureArgs,1);

preamble='Structure which contains the following fields';
preambleHTML=formatHTMLwithMATHJAX(preamble);


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


%     % Add the Remark after the table, if it is present
%     if ~isempty(posREMARK)
%         descrREMARK=descriinput(posREMARK:end);
%         descrREMARKHTML=formatHTMLwithMATHJAX(descrREMARK);
%
%         descrioutput=[preambleHTML iniTable Tablehtml cloTable '<p>' descrREMARKHTML '</p>'];
%     else
listStructureArgsFormatted=[preambleHTML iniTable Tablehtml cloTable];
%     end


end


function StringwithSpecialCharacters=SpecialCharacters(StringwithSpecialCharacters)
% Add symbol \ in before special characters in string StringwithSpecialCharacters
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


%FScategory:UTIHELP