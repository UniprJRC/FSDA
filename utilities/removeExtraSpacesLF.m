function newTxt = removeExtraSpacesLF(txt)
%removeExtraSpacesLF removes extra spaces and selected carriage returns from input string
%
%
%<a href="matlab: docsearchFS('removeExtraSpacesLF')">Link to the help function</a>
%
%
%    given an input string containing a series of carriage returns (CR) and
%    white spaces, it removes all carriage returns apart the cases below:
%           1) symbol ';'  is followed by one or more space and a CR; 
%           2) symbol ':'  is followed by one or more space and a CR;
%           2) symbol '.'  is followed by one or more space and a CR. 
%
%
%  Required input arguments:
%
%       txt : Input text. Character. String which has to be analysed. 
%
%  Optional input arguments:
%
%
%  Output: 
%
%   newTxt : Output text. Characted. String without carriage returns (apart
%            from those preceeded by symbols ':', ';' and '.') and without
%            extra spaces.
%
%
% See also: strtrim
%
%
% References:
%
%
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('removeExtraSpacesLF')">Link to the help function</a>
% Last modified 14-06-2016
%

% Examples

%{
    % Create a string with unnecesseray line feeds and without text justification.
    % Create an input string containing a series of unwanted features.
    FileWithFullPath=which('tclust');
    filename=FileWithFullPath;
    fileID = fopen(char(filename), 'r');
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    aa=regexp(fstring,'tclust partitions the points') ;
    bb=regexp(fstring,'constrained, Mahalanobis distances\.');
    str=fstring(aa:bb+35);
    % Remove from string descri all percentage signs
    posPercentageSigns=regexp(str,'%');
    str(posPercentageSigns)=[];
    str=[str 'x0ANow some wanted line feeds: x0A first item;   x0A   second item.'];
    str=regexprep(str,'x0A','\x0A');
    % Remove unnecessary spaces and extra line feeds from str but just keep
    % line break if before there was ':' or ';' or '.'
    a=removeExtraSpacesLF(str)
%}

%% Beginning of code

% Find position of wanted line feed
[~,goodLF]=regexp(txt,'[\:\;\.]\s*[\n\r]');
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
newTxt=txt;
end

