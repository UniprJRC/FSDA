function strFormatted = wraptextFS(str, varargin)
% wraptextFS formats long strings into wrapped text of specified width.
%
%
%<a href="matlab: docsearchFS('wraptextFS')">Link to the help function</a>
%
%
% This function not only does text wrapping, but also enables us:
% 1) to control left margin of the text;
% 2) to control the maximum width of the text or the right margin;
% 3) to add a (comment) sign at the beginning of each row of the wrapped text;
% 4) to indent the first line of the text.
% 5) to personalize comments, and left margin for comments. 
% This function uses routine strjoin and therefore can be used just by
% those who have a version of MATLAB>=2013a
%
%  Required input arguments:
%
%       str : Input text. Character vector. String which has to be analysed
%             and formatted.
%
%  Optional input arguments:
%
%  startcolumn :  Left margin of the text. Scalar (non negative integer).
%               This option controls the left margin of the text. The
%               default value of startcolumn is 1.
%               Example - 'startcolumn',10
%               Data Types - double
%
%    endcolumn :  Right margin of the text. Scalar (non negative integer).
%               This option controls the right margin of the text.
%               The default value of endcolumn is 75.
%               Example - 'endcolumn',50
%               Data Types - double
%
%    width :  width of the text. Scalar (non negative integer).
%               This option controls the width of the text.
%               The default value of width is 65.
%               Example - 'width',50
%               Data Types - double
%               Remark: it is necessary just to give two values among,
%               width, startcolumn and endcolumn because the third is
%               automatically determined
%
% firstline :  indentation for first line. Boolean. If firstline is true
%              then the first line starts in column 3 and not in column
%              startcolumn, while the text in all the other columns starts
%              as specified by option startcolumn.
%               The default value of firstline is false
%               Example - 'firstline',true
%               Data Types - Boolean
%
%   comment :  specify whether text is a Matlab comment. Boolean or structure. If
%              comment is true then the first character in each row will be
%              the percentage sign (comment symbol in Matlab). The default
%              value of comment is false. If comment is a structure it is
%              possible to personalize the symbol to put in from of each
%              row, and the left margin of the comment symbol. More
%              precisely, if comment is a structure it may contain the
%              following fields:
%               comment.commentsign = character(s) to be put at the beginning of
%               each row. String which identifies comment sign.
%               comment.startcolumn = starting column to include
%               commentsign.
%               Example - 'comment',true
%               Data Types - Boolean
%
%      code :  specify whether text is a Matlab code (with comments).
%              Boolean. If code than the extra space on the left is not
%              trimmed. The default value of code is false. Option code
%              must be set to true when we have to translate to .m file
%              code which contains wanted indentation.
%               Example - 'code',false
%               Data Types - Boolean
%
%
%  Output:
%
%   strFormatted : Output text. Character. Formatted string.
%                   Text starts in column specified by option startcolumn,
%                   the first line may have an indentation, and length of
%                   the text in each row cannot exceed the prespecified
%                   width.
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
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('wraptextFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples

%{
    %% wraptextFS with all default options.
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    if verLessThan('matlab','8.1') ==1
        warning('This function uses routine strjoin and works just with Matlab >=2013a')
    else
        Newstr=wraptextFS(str)
    end
%}

%{
    % start text in column 3 and put percentage sign at the beginning of
    % each line.
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    if verLessThan('matlab','8.1') ==1
        warning('This function uses routine strjoin and works just with Matlab >=2013a')
    else
        Newstr=wraptextFS(str,'comment',true,'startcolumn',3)
    end
%}

%{
    % Example of specification of startcolumn and text width.
    % Start text in column 5, the maximum text width is 40.
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    if verLessThan('matlab','8.1') ==1
        warning('This function uses routine strjoin and works just with Matlab >=2013a')
    else
        Newstr=wraptextFS(str,'comment',true,'startcolumn',10,'width',40)
    end
%}

%{
    % Add an indentation for first line.
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    if verLessThan('matlab','8.1') ==1
        warning('This function uses routine strjoin and works just with Matlab >=2013a')
    else
        Newstr=wraptextFS(str,'comment',true,'startcolumn',10,'width',40,'firstline',true)
    end
%}


%{
   % Use the width of command window.
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    if verLessThan('matlab','8.1') ==1
        warning('This function uses routine strjoin and works just with Matlab >=2013a')
    else
        startcolumn=10;
        cms = get(0,'CommandWindowSize');
        width = cms(1)-10;
        Newstr=wraptextFS(str,'comment',false,'startcolumn',startcolumn,'width',width)
    end
%}

%{
    % Example of input option comment supplied as structure.
    % Symbol '$$$' is included at the beginning of each row in column 5.
    % The width of the text is 60 and starts in column 12.
    comment=struct;
    comment.commentsign='$$$';
    comment.startcolumn=5;
    startcolumn=12;
    width=60;
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    if verLessThan('matlab','8.1') ==1
        warning('This function uses routine strjoin and works just with Matlab >=2013a')
    else
       Newstr=wraptextFS(str,'comment',comment,'startcolumn',startcolumn,'width',width);
    end
%}

%{ 
    % Example of use of option code.
    str='Paene insularum, Sirmio, insularumque ocelle, quascumque in liquentibus stagnis marique vasto fert uterque Neptunus, quam te libenter quamque laetus inviso, vix mi ipse credens Thuniam atque Bithunos liquisse campos et videre te in tuto. o quid solutis est beatius curis, cum mens onus reponit, ac peregrino labore fessi venimus larem ad nostrum, desideratoque acquiescimus lecto? hoc est quod unum est pro laboribus tantis. salve, o venusta Sirmio, atque ero gaude gaudente, vosque, o Lydiae lacus undae, ridete quidquid est domi cachinnorum.';
    if verLessThan('matlab','8.1') ==1
        warning('This function uses routine strjoin and works just with Matlab >=2013a')
    else
        out=xmlreadFS('tclust');
        ii=2;
        startcolumnEx=5;
        endcolumn=60;
        Ex=out.Ex;
        comment=struct;
        comment.commentsign='%';
        comment.startcolumn=startcolumnEx;
        endcolumnEx=60;
        i=2; jj=3;
        Exi=strtrim(Ex{i,jj});
        Exi{3,1}='%                        This is an example with extra spaces on the left which are trimmed';
        Exi{10,1}='for i=1:10';
        Exi{11,1}='    disp(i)';  
        % In this case the extra space on the left is wanted and it is not deleted
        Exi{12,1}='end';

        Eximod=Exi;
        for ii=1:size(Exi,1)
            % We must check whether it is comment or not
            % If it is a comment if the first character is symbol %
            Exii=Exi{ii,1};
            if ~isempty(Exii)
                if strcmp(Exii(1),'%')
                    % In this case strtrim is invoked inside wraptextFS (code is
                    % false)
                    descriFormatted=wraptextFS( Exii(2:end),'startcolumn',startcolumnEx,'endcolumn',endcolumn,'firstline',false,'comment',comment,'code',false);
                else
                    descriFormatted=wraptextFS( Exii,'startcolumn',startcolumnEx,'endcolumn',endcolumnEx,'firstline',false,'comment',false,'code',true);
                end
            end
            Eximod{ii,1}=descriFormatted;
        end

        % Before formatting
        disp(Exi)
        % After formatting
        disp([Eximod{:}])
    end

%}

%% Beginning of code 

% Input parameters checking


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
code=false;

% Write in structure 'options' the options chosen by the user
if nargin > 1
    options=struct('width',width,'startcolumn',startcolumn,'endcolumn',endcolumn,...
        'firstline',firstline,'comment',comment,'code',false);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSM:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Assign to each option structure argument the corresponding value
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
    code=options.code;
end

if isstruct(comment)
    fcomment=fieldnames(comment);
    
    
    d=find(strcmp('commentsign',fcomment));
    if d>0
        commentsign= comment.commentsign;
    else
        commentsign='%';
    end
    
    % Check if field startcolumn is present
    d=find(strcmp('startcolumn',fcomment));
    if d>0
        startcolumnsymbol= comment.startcolumn;
        if startcolumnsymbol>1
            commentsign=[repmat(' ',1,startcolumnsymbol-1) commentsign];
        end
    end
    
else
    
    if comment
        commentsign='%';
    else
        commentsign=[];
    end
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
% Include in the position of the spaces also the line feeds otherwise they
% will be overlooked
PosSpaces= sort([PosLF PosSpaces]);

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

% Left margin for text
leftMargin= {repmat(' ',1,startcolumn-1-length(commentsign))};

% Main code section: the loop below populates CellStackedStrings{i} with
% lines of desired length
while k <PosSpaces(end)
    if k==0 && firstline
        PosLinBreaks(i) = PosSpaces(find(PosSpaces<=k+width+startcolumn,1,'last'));
    else
        Posps=PosSpaces(find(PosSpaces<=k+width,1,'last'));
        if ~isempty(Posps)
            
            PosLinBreaks(i) = Posps;
        else
            PosLinBreaks(i)=length(str);
        end
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
            % deblank just in case the last hex character of strsel is 0A
            % type dec2hex(strsel) to see hexadecimal numbers of strsel
            CellStackedStrings{i}=deblank(strsel);
        else
            CellStackedStrings{i}=deblank(strsel);
        end
    end
    i = i+1;
end

CellStackedStrings=CellStackedStrings(1:i-1);

% First line indentation
% Insert at the end of each row symbol 0A123 which will be replaced by
% carriage return symbol \x0A
if firstline==1
    if code
        str=strjoin(strcat(commentsign,leftMargin,(CellStackedStrings(2:end)),'0A123'),'');
    else
        str=strjoin(strcat(commentsign,leftMargin,strtrim(CellStackedStrings(2:end)),'0A123'),'');
    end
    str=[commentsign CellStackedStrings{1} '0A123' str];
else
    if code
        str=strjoin(strcat(commentsign,leftMargin,(CellStackedStrings),'0A123'),'');
    else
        str=strjoin(strcat(commentsign,leftMargin,strtrim(CellStackedStrings),'0A123'),'');
    end
end

% strFormatted = output formatted string
strFormatted=regexprep(str,'0A123','\x0A');

end
%FScategory:UTIGEN
