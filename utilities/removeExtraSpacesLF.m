function newTxt = removeExtraSpacesLF(txt)
%removeExtraSpacesLF removes extra spaces and selected carriage returns from input string
%
%
%<a href="matlab: docsearchFS('removeExtraSpacesLF')">Link to the help function</a>
%
%
%    Given an input string possibly containing a series of carriage returns (CR) and white spaces,
%    removeExtraSpacesLF removes all carriage returns except those when:
%           1) symbol ';'  is followed by one or more spaces and a CR;
%           2) symbol ':'  is followed by one or more spaces and a CR;
%           3) symbol '.'  is followed by one or more spaces and a CR;
%           4) symbol '\[' is followed by one or more spaces and a CR;
%           5) symbol '\]' is preceded by one or more spaces and a CR.
%           6) symbol '\\' is preceded by one or more spaces and a CR.
%
%
%  Required input arguments:
%
%       txt : Input text. Character vector. String which has to be analysed.
%
%  Optional input arguments:
%
%
%  Output:
%
%   newTxt : Output text. Character. String without unwanted carriage
%            returns and extra spaces, as in cases 1-5 above.
%
%
% See also: strtrim, deblank
%
%
% References:
%
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('removeExtraSpacesLF')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples

%{
    % Create a string with unnecesseray line feeds and without text justification.
    % Create an input string containing a series of unwanted features.
    % The input string is extracted from the head of the FSDA function tclust.m.
    FileWithFullPath=which('tclust');
    filename=FileWithFullPath;
    fileID = fopen(char(filename), 'r');
    fstring=fscanf(fileID,'%c');
    % Starting and ending lines
    aa=regexp(fstring,'tclust partitions the points') ;
    bb=regexp(fstring,'constrained, Mahalanobis distances\.');
    % String
    str=fstring(aa:bb+35);

    % Remove from string all percentage signs
    posPercentageSigns=regexp(str,'%');
    str(posPercentageSigns)=[];

    str=[str 'x0ANow some wanted line feeds: x0A first item;   x0A   second item.'];
    str=regexprep(str,'x0A','\x0A');

    % Remove unnecessary spaces and extra line feeds from str but just keep
    % line break if before there was ':' or ';' or '.'
    a=removeExtraSpacesLF(str)
%}

%{
    % Create a string with a series of Latex equations.
    % The input string is extracted from the FSDA function tclust.m.
    FileWithFullPath=which('tclust');
    filename=FileWithFullPath;
    fileID = fopen(char(filename), 'r');
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    aa=regexp(fstring,'\\\[','once') ;
    bb=regexp(fstring,'\\\]','once');
    str=fstring(aa-145:bb+560);
    % Remove from string descri all percentage signs
    posPercentageSigns=regexp(str,'%');
    str(posPercentageSigns)=[];
    
    % str is the input string containing a series of Latex equations
    a=removeExtraSpacesLF(str)
%}

%% Beginning of code

if ~isempty(txt)
    % search for older Macs carriage returns (OS-9 and earlier) and replace
    % them with line feeds (unix and OS X and later)
    txt=regexprep(txt,'\x0D','\x0A');
    
    % Find position of wanted line feed
    [~,goodLF]=regexp(txt,'[\:\;\.\\\\]\s*[\n\r]');
    % [goodLF1st,goodLF1]=regexp(txt,'\\\[\s*[\n\r]');
    % Find 'line feed' then 'a sequence of spaces' then '\[' then 'another
    % sequence of spaces' then 'another line feed'
    [openLatexEqIni,openLatexEqEnd]=regexp(txt,'\x0A\s*\\\[\s*\x0A');
    
    [goodLF2,~]=regexp(txt,'[\n\r]\s*\\\]');
    [~,goodLF3]=regexp(txt,'\\\]\s*[\n\r]');
    
    % [goodLF2,goodLF3]=regexp(txt,'\x0A\s*\\\]\s*\x0A');
    
    goodLF=[goodLF openLatexEqEnd goodLF2 goodLF3 openLatexEqIni];
    
    allLF=regexp(txt,'[\n\r]');
    LFtoremove=setdiff(allLF,goodLF);
    % Remove unwanted line feeds
    txt(LFtoremove)=[];
    
    % Find the position of one or more spaces
    [a,b]=regexp(txt,'\s*');
    %[a,b]=regexp(stri,'\s*\x0A\s*');
    
    for i=length(a):-1:1
        % Check if stri(a(i):b(i)) contains a carriage return
        LFstri=regexp(txt(a(i):b(i)),'[\n\r]', 'once');
        if ~isempty(LFstri)
            % If there is a carriage return it must be left as is
            sel=setdiff(a(i):b(i),LFstri+a(i)-1);
            txt(sel)=[];
        else
            % we remove all unnecessary spaces
            txt(a(i)+1:b(i))=[];
        end
    end
    % Remove also unnecessary spaces at the beginning (if they are still
    % present)
    newTxt=strtrim(txt);
else
    newTxt=[];
end
end
%FScategory:UTIGEN
