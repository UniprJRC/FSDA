function strFormatted = wraptextFS(str, varargin)
% wraptextFS formats long strings into wrapped text of specified width.
%
%
%<a href="matlab: docsearchFS('wraptextFS')">Link to the help function</a>
%
%
% This function not only does text wrapping, but also enables us:
% 1) to control left margin of the text
% 2) to control the maximum width of the text or the right margin;
% 3) to add a (comment) sign at the beginning of each row of the wrapped text; 
% 4) to indent the first line of the text.
%
%
%  Required input arguments:
%
%       txt : Input text. Character vector. String which has to be analysed 
%             and formatted. 
%
%  Optional input arguments:
%
%  startcolumn :  Left margin of the text. Scalar (non negative integer).
%               This option controls the left margin of the text. The
%               default value of startcolumn is 1.
%               Example - 'startcolumn',10 
%               Data Types - double
%    endcolumn :  Right margin of the text. Scalar (non negative integer).
%               This option controls the right margin of the text.
%               The default value of endcolumn is 65.
%               Example - 'endcolumn',50 
%               Data Types - double
%    width :  width of the text. Scalar (non negative integer).
%               This option controls the width of the text.
%               The default value of width is 65.
%               Example - 'width',50 
%               Data Types - double
%               Remark: it is necessary just to give two values among,
%               width, startcolumn and endcolumn because the third is
%               automatically determined
% firstline :  indentation for first line. Boolean. if firstline is true
%              then the first line starts in column 3 and not in column
%              startcolumn, while the text in all the other columns starts
%              as specified by option startcolumn
%               The default value of firstline is false
%               Example - 'firstline',true 
%               Data Types - Boolean
%   comment :  specify whether text is a Maltab comment. Boolean. if
%              comment is true then the first character in each row will be
%              the percentage sign (comment symbol in Matlab). The deafult
%              value of comment is false
%               Example - 'comment',true 
%               Data Types - Boolean
%
%  Output: 
%
%   newTxt : Output text. Characted vector. Formatted string.
%            Text starts in column specified by option startcolumn, the
%            first line may have an indentation, and length of the text in
%            each row cannot exceed the prespecified width
%
%
% See also: strtrim, removeExtraSpacesLF
%
%
% References:
%
% Acknowledgements: 
%
%  This file had been inspired by function wraptext written by  Chad A. Greene of the University of Texas
% https://www.mathworks.com/matlabcentral/fileexchange/53176-wraptext
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('wraptextFS')">Link to the help function</a>
% Last modified 14-06-2016
%

% Examples

%{
    %% wraptextFS with all default options.
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    Newstr=wraptextFS(str);
%}

%{
    % start text in column 3 and put percentage sign at the beginning of
    % each line.
    Newstr=wraptextFS(str,'comment',true,'startcolumn',3)
%}

%{
    % Example of specification of startcolumn and text width.
    % Start text in column 5, the maximum text widh is 40.
    Newstr=wraptextFS(str,'comment',true,'startcolumn',10,'width',40)
%}

%{
    % Add an indentation for first line.
    Newstr=wraptextFS(str,'comment',true,'startcolumn',10,'width',40,'firstline',true)
%}


%{ 
   % Use the width of command window.
    startcolumn=10;
    cms = get(0,'CommandWindowSize');
    width = cms(1)-10;
    Newstr=wraptextFS(str,'comment',false,'startcolumn',startcolumn,'width',width)
%}
    
%% Input parameters checking

if nargin < 1
    error('FSDA:wraptextFS:missingInputs', ...
        'wraptextFS requires at least one input')
end

% Check the the input is a string
assert(ischar(str)==1,'Input str must be a string.')


% Set default parameters
width=75;
startcolumn=1;
endcolumn=75;
firstline=false;
comment=false;

% Write in structure 'options' the options chosen by the user
if nargin > 1
    options=struct('width',width,'startcolumn',startcolumn,'endcolumn',endcolumn,...
        'firstline',firstline,'comment',comment);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    
    % If the user has specified both startcolumn and endcolumn
    if sum(strcmp(UserOptions,'startcolumn')) && sum(strcmp(UserOptions,'endcolumn'))
        startcolumn=options.startcolumn;
        endcolumn=options.endcolumn;
        width=endcolumn-startcolumn+1;
        
        % If the user has specified both width and endcolumn
    elseif sum(strcmp(UserOptions,'width')) && sum(strcmp(UserOptions,'endcolumn'))
        startcolumn=options.endcolumn-width;
        width=options.width;
        
        % the user has specified both startcolumn and width
    elseif sum(strcmp(UserOptions,'width')) && sum(strcmp(UserOptions,'startcolumn'))
        
        startcolumn=options.startcolumn;
        width=options.width;
        
        % the user has specified only width
    elseif sum(strcmp(UserOptions,'width'))
        width=options.width;
        
        % the user has specified only startcolumn
    elseif sum(strcmp(UserOptions,'startcolumn'))
        startcolumn=options.startcolumn;
        
        % the user has specified only endcolumn
    elseif sum(strcmp(UserOptions,'endcolumn'))
        startcolumn=  options.startcolumn-options.width-1;
    else
        startcolumn=options.startcolumn;
        width=options.width;
    end
    
    firstline=options.firstline;
    comment=options.comment;
    %     if endcolumn-startcolumn+1 ~=width
    %         error('wrong numebr of cols specified')
    %     end
end

if comment
    commentsign='%';
else
    commentsign=[];
end


% Add a space at the end of the string
str=[str ' '];

% PosSpaces is the vector which contains the position of spaces inside
% input string

% str=regexprep(str,'\\[','\x0A\x0A\\[\x0A');
% str=regexprep(str,'\\]','\x0A\\]\x0A\x0A');


% Find the position of carriage returns
PosLF = regexp(str,'\x0A');
if ~isempty(PosLF) && PosLF(end)==length(str)-1
    str(end-1)=[];
end

PosSpaces = regexp(str,' ');


% Throw error if any words are longer than specified width:
if any(diff(PosSpaces)>width)
    error('FSDA:wraptextFS:TooSmallWidth','Some words in the input string are longer than the specified width. Either increase the width allowance or break up long words with a hyphen followd by a space.')
end

k = 0; % counter
i = 1;

% PosLinBreaks is the vector which contains the position of line breaks
% inside input string
PosLinBreaks=zeros(length(PosSpaces),1);

%  CellStackedStrings is the cell which contains the original string. Each
%  row refers to a line. It is initialized with length(PosSpaces), however,
%  its length will be much shorter
CellStackedStrings=cell(length(PosSpaces),1);


leftMargin= {repmat(' ',1,startcolumn-1-comment)};

while k <PosSpaces(end)
    if k==0 && firstline
        PosLinBreaks(i) = PosSpaces(find(PosSpaces<=k+width+startcolumn,1,'last'));
    else
        PosLinBreaks(i) = PosSpaces(find(PosSpaces<=k+width,1,'last'));
    end
    
    k = PosLinBreaks(i);
    if i>1
        strsel=str(PosLinBreaks(i-1)+1:PosLinBreaks(i));
        findLFinstrsel=regexp(strsel,'\x0A', 'once');
        % if carriage return is present and is not the penultimate
        % character (penultimate because we artificially added a space at
        % the end of the input string)
        if ~isempty(findLFinstrsel) && findLFinstrsel<length(strsel)-1
            PosLinBreaks(i)=PosLinBreaks(i-1)+findLFinstrsel;
            k = PosLinBreaks(i);
            CellStackedStrings{i}=(str(PosLinBreaks(i-1)+1:PosLinBreaks(i)));
        else
            CellStackedStrings{i}=strsel;
        end
    else
        strsel=str(1:PosLinBreaks(i));
        
        findLFinstrsel=regexp(strsel,'\x0A', 'once');
        if ~isempty(findLFinstrsel) && findLFinstrsel<length(strsel)-1
            PosLinBreaks(i)=findLFinstrsel;
            k = PosLinBreaks(i);
            strsel=str(1:PosLinBreaks(i));
            CellStackedStrings{i}=strsel;
        else
            CellStackedStrings{i}=strsel;
        end
    end
    i = i+1;
end

CellStackedStrings=CellStackedStrings(1:i-1);

% First line indentation
% Insert at the end of each row symbol 0A123 which will be replaced by
% carriage return symbol \x0A
if firstline==1
    str=strjoin(strcat(commentsign,leftMargin,strtrim(CellStackedStrings(2:end)),'0A123'),'');
    str=[commentsign CellStackedStrings{1} '0A123' str];
else
    str=strjoin(strcat(commentsign,leftMargin,strtrim(CellStackedStrings),'0A123'),'');
end

% strFormatted = output formatted string
strFormatted=regexprep(str,'0A123','\x0A');

end

