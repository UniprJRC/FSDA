function publishFStmp(file)
% publishFS('Score.m')
% fullPathToScript = locateFile(file);
% if isempty(fullPathToScript)
%     error('FSDA:publishFS:WrongFile','SourceNotFound');
% end

FileWithFullPath=which(file);
[pathstr,name,ext]=fileparts(FileWithFullPath);

if ~strcmp('.m',ext)
    error('FSDA:publishFS:WrongFileExt','Input file must have m extension')
end

if isempty(pathstr)
    error('FSDA:publishFS:WrongFile','SourceNotFound');
end

filename=FileWithFullPath;
% f = fopen(filename);

fileID = fopen(char(filename), 'r+');


% Insert the file into fstring
fstring=fscanf(fileID,'%c');

%% Add title
beforetitl=['<!DOCTYPE HTML> \r'  ...
    '<html itemscope="" xmlns="http://www.w3.org/1999/xhtml">\r' ...
    '<head>\r' ...
    '<title>\r'];
aftertitle='</title>';
titl=sprintf([beforetitl    name  aftertitle]);

%% Add purpose of the file (what is in the second row of .m file)
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
    '<script src="includesFS/jquery-latest.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/l10n.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/docscripts.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.toc.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.imagescaling.js" type="text/javascript"></script>\r'...
    '<script src="includesFS/mw.imageanimation.js" type="text/javascript"></script>\r'...
    '<link href="includesFS/reset.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/960.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/site5.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includesFS/doc_center_installed.css" rel="stylesheet" type="text/css">\r'...
    '<link href="includes/product/css/doc_center_print.css" media="print" rel="stylesheet" type="text/css">\r'...
    '</head>\r'...
    '<body>\r'];
metacontent=sprintf([beforemetacontent purpose aftermetacontent]);

%% CONTAINER + SINTAX SECTION
inisitecont=sprintf(['<div class="site_container site_toc_opened">\r'...
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

% Required input arguments:

finsitecont=sprintf(['<h2 id="syntax">Syntax</h2>\r'...
    '<div class="syntax_signature">\r'...
    '<div class="syntax_signature_module">\r'...
    '<ul>\r']);

sitecont=[inisitecont htmlsitecont finsitecont];

% Find point where first argument are declared

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

% Find number of compulasory input arguments
% nargin = number of commas in string  functionname(inp1,inp2,inp3, ...)
[startIndexInp] = regexp(fstring,'(');
[endIndexInp] = regexp(fstring,')');
% inputargs =  (inp1,inp2,inp3, ...)
inputargs=fstring(startIndexInp(1):endIndexInp(1));
% Check if inputargs contains the string varargin
[optargs1]=regexp(inputargs,'varargin');
[commasIn] = regexp(inputargs,',');
if isempty(optargs1)
    sintax=cell(nargout,1);
    sintax{1}=[outargs(2:commasOut(1)-1) '=' name inputargs];
    j=2;
    nREQargin=length(commasIn)+1; % nREQargin = number of required input args
else
    % In this case there is also Name Value line
    sintax=cell(nargout+1,1);
    strinputarg=strtrim(inputargs(1:optargs1-1));
    if strcmp(strinputarg(end),',')
        strinputarg=strinputarg(1:end-1);
    end
    
    %sintax{1}=[outargs(2:commasOut(1)-1) '=' name strtrim(inputargs(1:optargs1-2)) ')'];
    %sintax{2}=[outargs(2:commasOut(1)-1) '=' name strtrim(inputargs(1:optargs1-2)) ',Name,Value)'];
    
    sintax{1}=[outargs(2:commasOut(1)-1) '=' name strinputarg ')'];
    sintax{2}=[outargs(2:commasOut(1)-1) '=' name strinputarg ',Name,Value)'];
    
    j=3;
    nREQargin=length(commasIn); % nREQargin = number of required input args
end



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


%% DESCRIPTION SECTION



inidescription=sprintf(['	<div class="ref_sect" itemprop="content">\r'...
    '							<h2 id="description">Description</h2>\r'...
    '							<div class="descriptions">\r'...
    '								<div class="description_module">\r']);

descriptionhtml='';

disp('Examples')
disp(sintax)
disp('---------------')

% start of example j
[startIndexEx] = regexp(fstring,'%\{');
[endIndexEx] = regexp(fstring,'%\}');

%listEx = list which contains the examples (associated to sintax)
% First column= title, second column detailed description. Third column
% code
listEx=cell(length(sintax),3);

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
    stri=fstring(startIndexEx(j)+2:endIndexEx(j)-1);
    % What is before the first full stop is the title.
    % What is after the second full stop is the description
    % The first line which does not contain symbol % is the beginning of the
    % code
    [endtitle] = regexp(stri,'\.','once');
    strititle=stri(1:endtitle);
    % Remove percentage signs if present.
    posPercentageSigns=regexp(strititle,'%');
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
    listEx{j,3}=stri(findescriptionEx+1:end);
    
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

%% EXAMPLES
% the examples which are inside %{   %} are put here.
% The first sentence which end with a full stop is the title of the example
iniexamples=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<div class="examples">\r'...
    '<h2 id="examples">Examples</h2>\r'... % start of expandable examples
    '<div id="expandableExamples" class="expandableContent">\r'...
    '<p class="switch"><a class="expandAllLink"' ...
    'href="javascript:void(0);">expand all</a>' ...
    '</p>']);





% for j=1:length(sintax)
%
% end



exampleshtml='';
for j=1:length(sintax)
    exampleshtml=[exampleshtml  sprintf(['<div id="example_' num2str(j) '" class="example_module expandableContent">\r'...
        '<div id="ex' num2str(j) '">\r'...
        '</div>\r'...
        '<h3 class="expand"><span>\r'...
        '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        '<span class="example_title">']) listEx{j,1} sprintf(['</span></a></span></h3>\r'...
        '<div class="collapse">\r'...
        '<p>']) listEx{j,2} sprintf(['<div class="programlisting">\r'...
        '<div class="codeinput"><pre>\r']) ...
        listEx{j,3} sprintf(['</pre></div>\r'...
        '</div>\r'])...
        sprintf(['</p>\r'...
        '</div>\r'...
        '</div>\r'])]; % close div id="example_j"
end

if length(startIndexEx)>length(sintax)
    lsintax=length(sintax);
    listextraEX=cell(length(startIndexEx)-lsintax,3);
    
    
    for j=1:size(listextraEX,1)
        stri=fstring(startIndexEx(j+lsintax)+2:endIndexEx(j+lsintax)-1);
        % What is before the first full stop is the title.
        % What is after the second full stop is the description
        % The first line which does not contain symbol % is the beginning of the
        % code
        [endtitle] = regexp(stri,'\.','once');
        strititle=stri(1:endtitle);
        % Remove percentage signs if present.
        posPercentageSigns=regexp(strititle,'%');
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
        listExtraEx{j,3}=stri(findescriptionEx+1:end);
    end
else
    listExtraEx='';
end

closeexamples=sprintf(['</div>\r'... % close div id="expandableExamples
    '<p> </p>\r']);
% Related examples are below
iniRelatedExamples=sprintf('<h3 class="bottom_ruled">Related Examples</h3>\r');

RelatedExamples='';
for j=1:size(listextraEX,1)
    RelatedExamples=[RelatedExamples  sprintf(['<div id="example_' num2str(j) '" class="example_module expandableContent">\r'...
        '<div id="ex' num2str(j) '">\r'...
        '</div>\r'...
        '<h3 class="expand"><span>\r'...
        '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        '<span class="example_title">']) listExtraEx{j,1} sprintf(['</span></a></span></h3>\r'...
        '<div class="collapse">\r'...
        '<p>']) listExtraEx{j,2} sprintf(['<div class="programlisting">\r'...
        '<div class="codeinput"><pre>\r']) ...
        listExtraEx{j,3} sprintf(['</pre></div>\r'...
        '</div>\r'])...
        sprintf(['</p>\r'...
        '</div>\r'...
        '</div>\r'])]; % close div id="example_j"
end

% RelatedExamples=sprintf(['<ul>\r'...
%     '<li>rel1</li>\r'...
%     '<li>rel2 </li>\r'...
%     '<li>rel3 </li>\r'...
%     '</ul>\r']);

closeallex=sprintf(['</div>\r'... % div class="examples"
    '</div>']);	 % close class="ref_sect

examples=[iniexamples exampleshtml closeexamples iniRelatedExamples...
    RelatedExamples closeallex];

%% INPUT ARGUMENTS
iniinputargs=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<h2 id="inputs">Input Arguments</h2>\r'...
    '<div class="expandableContent">\r'...
    '<div class="arguments">\r'...
    '<div class="input_argument_container">\r'...
    '<p class="switch">\r'...
    '<a class="expandAllLink" href="javascript:void(0);">\r'...
    'expand all</a></p>']);

%       inputargs
%[optargs1,optargs2]=regexp(inputargs,'varargin');
%[commasIn] = regexp(inputargs,',');
%if isempty(optargs1)

% in string inpi there will the required input argument


% Write all required input arguments in cell listargins
listargins=cell(nREQargin,1);
for i=1:nREQargin
    if nREQargin>1
        if i==1
            inpi=inputargs(2:commasIn(i)-1);
        else
            inpi=inputargs(commasIn(i-1)+1:commasIn(i)-1);
        end
    else
        if isempty(optargs1) % if there are no optional input arguments
            inpi=inputargs(2:end-1);
        else
            inpi=inputargs(2:commasIn(i)-1);
        end
    end
    inpi=strtrim(inpi);
    listargins{i}=inpi;
end

reqargs='';
for i=1:nREQargin
    
    inpi=listargins{i};
    
    
    insel=regexp(fstring,'Required input arguments:');
    if isempty(insel)
        disp('Please check HTML input file')
        error('FSDA:missInps','HTML file does not contain ''Required input arguments:'' string')
    end
    
    % substring to search start from Required input arguments:
    fstringsel=fstring(insel(1):end);
    inipoint=regexp(fstringsel,[listargins{i} '\s{0,10}:']);
    
    % The endpoint of the substring is See also or the next output argument
    if i <nREQargin
        endpoint=regexp(fstringsel,[listargins{i+1} '\s{0,10}:']);
    else
        endpoint=regexp(fstringsel,'Optional input arguments:');
        if isempty(endpoint)
            disp('Please check HTML input file')
            error('FSDA:missOuts','HTML file does not contain ''Optional input arguments:'' string')
        end
    end
    % descri = string which contains the description of i-th input
    % argument
    descriinput=fstringsel((inipoint(1)+length(listargins{i})+2):endpoint(1)-1);
    
    % Remove from string descri all '% signs
    posPercentageSigns=regexp(descriinput,'%');
    descriinput(posPercentageSigns)=[];
    % Remove from string descri leading and trailing white spaces
    descriinput=strtrim(descriinput);
    % what is before the first comma or the first full stop is the
    % preamble, the reset in the description
    posfirstcomma=regexp(descriinput,',','once');
    posfirstfullstop=regexp(descriinput,'\.','once');
    sep=min([posfirstcomma posfirstfullstop]);
    if isempty(sep)
        warning('FSDA:MISSdescr',['In the description of ' listargins{i} ' String ''' descriinput ''' does not contain symbols '','' or ''.'''])
        preamble=descriinput;
    else
        preamble=descriinput(1:sep-1);
    end
    preamble=strtrim(preamble);
    % remove sign : if present at the beginning of the sentence
    if strcmp(preamble(1),':')
        preamble=strtrim(preamble(2:end));
    end
    descriinput=descriinput(sep+1:end);
    % Start with upper case character
    descriinput=strtrim(descriinput);
    if ~isempty(descriinput)
        descriinput=[upper(descriinput(1)) descriinput(2:end)];
    end
    
    % disp(inpi)
    reqargs=[reqargs sprintf(['<div class="expandableContent">\r'...
        ' <div id="inputarg_' inpi '" class="clearfix">\r'...
        ' </div>\r'...
        ' <h3 id="input_argument_' inpi '" class="expand">\r'...
        ' <span>\r'...
        ' <a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        ' <span class="argument_name"><code>' inpi '</code> &#8212; \r'...  % &#8212; = long dash
        ' Data</span></a><span class="example_desc">' preamble '</span></span></h3>\r'...
        ' <div class="collapse">\r'...
        ' <p>']) descriinput sprintf(['</p>\r'...
        ' <p class="datatypelist"><strong>\r'...
        ' Data Types: </strong><code>single</code> \r'...
        ' | <code>double</code></p>\r'...
        ' </div>\r'...
        ' </div>\r'])];
    
end


% reqargs=sprintf(['<div class="expandableContent">\r'...
% 			'<div id="inputarg_X" class="clearfix">\r'...
% 			'</div>\r'...
% 			'<h3 id="input_argument_X" class="expand">\r'...
% 			'<span>\r'...
% 			'<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
% 			'<span class="argument_name"><code>X</code> &#8212; \r'...  % &#8212; = long dash
% 			'Data</span></a><span class="example_desc">numeric\r'...
% 			'matrix</span></span></h3>\r'...
% 			'<div class="collapse">\r'...
% 				'<p>TOWRITE</p>\r'...
% 				'<p class="datatypelist"><strong>\r'...
% 				'Data Types: </strong><code>single</code> \r'...
% 				'| <code>double</code></p>\r'...
% 			'</div>\r'...
% 		'</div>\r']);


%% OPTIONAL INPUT ARGUMENTS PART

insel=regexp(fstring,'Optional input arguments:');
if isempty(insel)
    disp('Please check HTML input file')
    error('FSDA:missInps','HTML file does not contain ''Optional input arguments:'' string')
end

% substring to search start from Optional input arguments:
fstringsel=fstring(insel(1):end);

endpoint=regexp(fstringsel,'Output:');
if isempty(endpoint)
    disp('Please check HTML input file')
    error('FSDA:missOuts','HTML file does not contain ''Output:'' string')
end
fstringsel=fstringsel(1:endpoint-2);

% Find any string which
% begins with % sign then
% contains a series of white space which can go from 0 to 20 then
% contains any single word
% a series of white spaces which can go from 0 to 10 then
% character :
% [ini,fin]=regexp(fstringsel,'%\s{0,20}\w*\s{0,10}:');
[ini,fin]=regexp(fstringsel,'%\s{0,15}\w*\s{0,10}:');

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

listOptArgs=cell(length(ini),4);

ij=1;
for i=1:length(ini)
    % fin(i)-1 because character ':' does not have to be extracted
    opti=fstringsel(ini(i):fin(i)-1);
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
            descradd=fstringsel(ini(i):ini(i+1));
        else
            descradd=fstringsel(ini(i):end);
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
            descrtosplit=fstringsel(fin(i)+1:ini(i+1)-1);
        else
            % TOCHECK
            % descrtosplit=fstringsel(fin(i)+1:ini(i+1)-1);
            descrtosplit=fstringsel(fin(i)+1:end);
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
            jjjj=1;
        end
        
        
        listOptArgs{ij,3}=descrtype;
        
        try
            descrlong=strtrim(descrtosplit(inifullstops(2)+1:end));
        catch
            ddd=1;
        end
        
        % Check if descr long contains
        % Example - and Data types -
        
        CheckExample=regexp(descrlong,'Example -','once');
        if ~isempty(CheckExample)
            Datatypes=regexp(descrlong,'Data Types -','once');
            listOptArgs{ij,4}=descrlong(1:CheckExample-1);
            
            % The first word of example code must be embedded around tags <code> </code>
            examplecode=descrlong(CheckExample+10:Datatypes-1);
            posspace=regexp(examplecode,' ');
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

disp('Optional arguments found')
disp(listOptArgs(:,1))

% codewithexample=['''Distance'',''cosine'',''Replicates'',10,' ...
%     '''Options'',statset(''UseParallel'',1)'];
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
optargs=sprintf(['<div id="namevaluepairarguments" class="clearfix">\r'...
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


optargsexp='';
for i=1:size(listOptArgs,1);
    longdescription=listOptArgs{i,4};
    nameoptarg=listOptArgs{i,1};
    shortdesc=listOptArgs{i,3};
    titloptarg=listOptArgs{i,2};
    
    % datatype = type of data for that particular option
    %     examplecode=['''Display'',''final'''];
    %     datatype='char';
    
    examplecode=listOptArgs{i,5};
    datatype=listOptArgs{i,6};
    
    optargsexp=[optargsexp sprintf(['<div class="expandableContent">\r'...
        '<div id="inputarg_Display" class="clearfix">\r'...
        '</div>\r'...
        '<h3 id="input_argument_namevalue_display" class="expand">\r'...
        '<span>\r'...
        '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
        '<span class="argument_name"><code>' nameoptarg  '</code> \r'...
        '&#8212;' titloptarg '</span></a><span class="example_desc">' shortdesc  '</span></span></h3>\r'...
        '<div class="collapse">\r'...
        '	<p>']) longdescription sprintf(['</p>\r'...
        '	<p class="description_valueexample">\r'...
        '       <strong>Example: </strong>' examplecode '</p>\r'...
        '	<p class="datatypelist"><strong>Data Types: </strong><code>' datatype '</code></p>\r'...
        '</div>\r'...
        '</div>'])];
end
% CLOSE OPT ARGS

closeoptargs=sprintf(['</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r']);
inputargs=[iniinputargs reqargs  optargs optargsexp closeoptargs];


%% OUTPUT ARGUMENTS

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
    %     if nargout>1
    %             if i==1
    %                 outi=outs(2:commasOut(i)-1);
    %             elseif i==nargout
    %                 outi=outs(commasOut(i-1)+1:end-1);
    %             else
    %                 outi=outs(commasOut(i-1)+1:commasOut(i)-1);
    %             end
    %     else
    %         outi=outs;
    %     end
    
    % listargouts is a cell which contains the list of output arguments
    outi=listargouts{i};
    
    outsel=regexp(fstring,'Output:');
    if isempty(outsel)
        disp('Please check HTML input file')
        error('FSDA:missOuts','HTML file does not contain ''Output:'' string')
    end
    
    % substring to search start from Output:
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
    
    % The endpoint of the substring is See also or the next output argument
    if i <nargout-1
        endpoint=regexp(fstringsel,listargouts{i+1});
    elseif i==nargout-1
        
        if strcmp(listargouts{end},'varargout') ==0
            endpoint=regexp(fstringsel,listargouts{i+1});
        else
            % In this case there are also optional arguments
            endpoint=regexp(fstringsel,'Optional Output:');
        end
        
    else
        endpoint=regexp(fstringsel,'See also');
        if isempty(endpoint)
            disp('Please check HTML input file')
            error('FSDA:missOuts','HTML file does not contain ''See also:'' string')
        end
    end
    
    
    
    % descri = string which contains the description of i-th output
    % argument
    try
        descrioutput=fstringsel((inipoint(1)+length(listargouts{i})+2):endpoint(1)-1);
    catch
        disp(['FSDA:WrongOut','Could not process correctly output argument ' listargouts{i}])
        ddd=1;
    end
    
    % Remove from string descri all '% signs
    posPercentageSigns=regexp(descrioutput,'%');
    descrioutput(posPercentageSigns)=[];
    % Remove from string descri leading and trailing white spaces
    descrioutput=strtrim(descrioutput);
    
    % Check if the output is a structure. If this is the case
    checkifstructure=regexp(descrioutput,[outi '\.'],'once');
    
    
    if ~isempty(checkifstructure)
        [ini,fin]=regexp(descrioutput,['\s{0,8}' outi '\.\w*\s{0,8}=']);
        if isempty(ini)
            disp('Probably ":" symbols  must be replaced with "=" symbols in out description')
            error('FSDA:MissingArg',['Parser cannot find string' outi '.xxxx = for output structure'])
        else
        end
        
        listStruArgs=cell(length(ini),2);
        
        for k=1:length(ini)
            % fin(i)-1 because character ':' does not have to be extracted
            StructFieldName=descrioutput(ini(k):fin(k)-1);
            % Remove from string opti leading and trailing white spaces
            StructFieldName=strtrim(StructFieldName);
            StructFieldName=StructFieldName(length(outi)+2:end);
            
            % Store name in the first column
            listStruArgs{k,1}=StructFieldName;
            % Store short description in the 3nd col of listOptArgs
            if k<length(ini)
                StructFieldCont=descrioutput(fin(k)+1:ini(k+1)-1);
            else
                StructFieldCont=descrioutput(fin(k)+1:end);
            end
            
            % Store name in the first column
            listStruArgs{k,2}=StructFieldCont;
        end
        
        % rowtodel = vector which contains the duplicate rows of
        % listStruArgs which have to be deleted
        inisel=1:size(listStruArgs,1);
        
        % Check if cell listStruArgs contains duplicates in the first column
        for j=2:size(listStruArgs,1)
            if strcmp(listStruArgs{j,1},listStruArgs{j-1,1})
                listStruArgs{j-1,2}=[listStruArgs{j-1,2} listStruArgs{j,1} listStruArgs{j,2}];
                inisel(j)=999;
            end
        end
        % remove from inisel the rows equal to 999 (that is the rows which
        % correspond to duplicated arguments)
        inisel(inisel==999)=[];
        
        iniTable=sprintf(['<table border="2" cellpadding="4" cellspacing="0" class="body">\r'...
            '<colgroup>\r'...
            '<col width="50%"><col width="50%">\r'...
            '</colgroup>\r'...
            '<thead>\r'...
            '<tr valign="top">\r'...
            '<th valign="top">Value</th>\r'...
            '<th valign="top">Description</th>\r'...
            '</tr>\r'...
            '</thead>']);
        
        Tablehtml='';
        for k=inisel % length(ini)
            Tablehtml=[Tablehtml sprintf(['<tr valign="top">\r'...
                '<td><code>' listStruArgs{k,1} '</code></td>\r'...
                '<td>\r'...
                '<p>']) listStruArgs{k,2} sprintf(['</p>\r'...
                '</td>\r'...
                '</tr>'])];
        end
        
        cloTable='</table>';
        
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
        % 'containing';
        poswhichcontains=regexp(descrioutput,'which contains');
        if ~isempty(poswhichcontains) && poswhichcontains(1)<120
            preamble=descrioutput(1:poswhichcontains(1)-1);
            descrioutput=descrioutput(poswhichcontains(1)+14:end);
            % Remove word the at the beginning of the sentence and starts with
            % uppercase
            StartsWithThe=regexp(descrioutput,'the');
            if StartsWithThe(1)<4
                descrioutput=descrioutput(StartsWithThe(1)+4:end);
            end
            descrioutput=strtrim(descrioutput);
            descrioutput=[upper(descrioutput(1)) descrioutput(2:end)];
        else
            poscontaining=regexp(descrioutput,'containing');
            if ~isempty(poscontaining) && poscontaining(1)<120
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
                preamble='TOWRITE';
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

%% MORE ABOUT SECTION

Moreabout=sprintf(['<div class="moreabout ref_sect">\r'...
    '<h2 id="moreabout">More About</h2>\r'...
    '<div class="expandableContent">\r'...
    '<p class="switch">\r'...
    '<a class="expandAllLink" href="javascript:void(0);">\r'...
    'expand all</a></p>\r'...
    '<div class="expandableContent" itemprop="content">\r'...
    '<h3 class="expand"><span>\r'...
    '<a href="javascript:void(0);" style="display: block;" title="Expand/Collapse">\r'...
    '<span>ADDITIONAL EXAMPLE</span></a></span></h3>\r'...
    '<div class="collapse">\r'...
    '<p>describe describe </p>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>']);

%% REFERENCES

iniref=regexp(fstring,'References:');
if isempty(iniref)
    warning('FSDA:missInp',['File ' name '.m does not contain string ''References:'''])
    refsargs={''};
else
    
    endref=regexp(fstring,'Copyright');
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

Referenceshtml='';
iniReferences=sprintf(['<div class="ref_sect" itemprop="content">\r'...
    '<div class="bibliography">\r'...
    '<h2 id="references">References</h2 \r']);

for i=1:length(refsargs)
    Referenceshtml=sprintf([Referenceshtml  '<div><p>' refsargs{i} '</p></div>\r']);
end
Referencesclose=sprintf(['</div>\r'...
    '</div>']);
References=[iniReferences Referenceshtml Referencesclose];

%% SEE ALSO

iniSeealso=sprintf(['<div class="ref_sect">\r'...
    '<h2>See Also</h2>\r'...
    '<p>\r']);

iniref=regexp(fstring,'See also','once');
endref=regexp(fstring,'References','once');
seealsostr=fstring(iniref+8:endref-1);

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
        Seealsoitem= seealsostr;
    else
        if i==nseealso
            Seealsoitem=seealsostr(poscommas(i-1)+1:end);
        elseif i==1
            Seealsoitem=seealsostr(1:poscommas(i)-1);
        else
            Seealsoitem=seealsostr(poscommas(i-1)+1:poscommas(i)-1);
        end
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

clos=sprintf(['<h1>Automatically generated by FSDA parser</h1>\r'...
    '</div>\r'...
    '</section>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</div>\r'...
    '</body>'...
    '</html>']);




fclose(fileID);

outstring=([titl metacontent sitecont sintaxhtml sintaxclose description  ....
    examples inputargs outargs Moreabout References Seealso clos]);


file1ID=fopen([name 'tmp.html'],'w');
fprintf(file1ID,'%s',outstring);
fclose('all');


end
