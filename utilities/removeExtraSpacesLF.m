function newstr = removeExtraSpacesLF(stri)
%removeExtraSpacesLF removes extra spaces and selected carriage returns from input string
%

%{
    % Create a string with  unnecesseray line feeds and without text justification
    FileWithFullPath=which('tclust');
    filename=FileWithFullPath;
    fileID = fopen(char(filename), 'r');
    % Insert the file into fstring
    fstring=fscanf(fileID,'%c');
    aa=regexp(fstring,'tclust partitions the points') ;
    bb=regexp(fstring,'constrained, Mahalanobis distances\.');
    str=fstring(aa:bb+35);
    % Remove from string descri all '% signs
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
[~,goodLF]=regexp(stri,'[\:\;\.]\s*[\n\r]');
allLF=regexp(stri,'[\n\r]');
LFtoremove=setdiff(allLF,goodLF);
% Remove unwanted line feeds
stri(LFtoremove)=[];

% Find the position of one or more spaces
[a,b]=regexp(stri,'\s*');
%[a,b]=regexp(stri,'\s*\x0A\s*');

for i=length(a):-1:1
    % Check if stri(a(i):b(i)) contains a carriage return
    LFstri=regexp(stri(a(i):b(i)),'[\n\r]', 'once');
    if ~isempty(LFstri)
        % If there is a carriage return it must be left as is
        sel=setdiff(a(i):b(i),LFstri+a(i)-1);
        stri(sel)=[];
    else
        % we remove all unnecessary spaces
        stri(a(i)+1:b(i))=[];
    end
end
newstr=stri;
end

