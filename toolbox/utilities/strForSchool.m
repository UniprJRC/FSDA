function [str] = strForSchool(header, corpus, footer, classes)

%% Beginning of code

[n,p]=size(corpus);

% classes is the string array to put in the second column of the table
if nargin<4
    classes='';
end

% Option which controls how to display the numbers inside the GUI.
% If rightAlignment is true numbers are right aligned else they are
% centered.
rightAlignment=true;

if rightAlignment==true
    FirstLine='$\begin{array}{rrrrrrrr} \\ $\boldmath{$';
else
    FirstLine='$\matrix{ $\boldmath{$';
end

for i=1:length(header)
    if i<length(header)
        FirstLine=[FirstLine header{i} '$}$ & $\boldmath{$'];
    else
        if rightAlignment==true
            FirstLine=[FirstLine header{i} '$}$ \\'];
        else
            FirstLine=[FirstLine header{i} '$}$ \cr'];
        end
    end
end

TableCorpus='';
if size(corpus)<=13
    for i=1:n
        for j=1:p
            corpusij=corpus(i,j);
            if isnan(corpusij)
                numij='-';
            else
                numij=num2str(corpusij);
            end

            if j==2 && isstring(classes)
                TableCorpus=[TableCorpus ' ' classes{i}  ' &'];

            elseif j<size(corpus,2)
                TableCorpus=[TableCorpus ' ' numij ' &'];
            else
                if rightAlignment==true
                    TableCorpus=[TableCorpus ' ' numij ' \\'];
                else
                    TableCorpus=[TableCorpus ' ' numij ' \cr'];
                end
            end
        end
    end
else
    for i=1:6
        for j=1:p
            numij=num2str(corpus(i,j));

            if j==2 && isstring(classes)
                TableCorpus=[TableCorpus ' ' classes{i}  ' &'];

            elseif j<p
                TableCorpus=[TableCorpus ' ' numij ' &'];
            else
                if rightAlignment==true
                    TableCorpus=[TableCorpus ' ' numij ' \\'];
                else
                    TableCorpus=[TableCorpus ' ' numij ' \cr'];
                end
            end
        end
    end

    % Insert simbol ... for the lines which are not shown
    for j=1:p
        numij='...';
        if j<p
            TableCorpus=[TableCorpus ' ' numij ' &'];
        else
            if rightAlignment==true
                TableCorpus=[TableCorpus ' ' numij ' \\'];
            else
                TableCorpus=[TableCorpus ' ' numij ' \cr'];
            end

        end
    end

    for i=n-5:n
        for j=1:p
            numij=num2str(corpus(i,j));

            if j==2 && isstring(classes)
                TableCorpus=[TableCorpus ' ' classes{i}  ' &'];

            elseif j<size(corpus,2)
                TableCorpus=[TableCorpus ' ' numij ' &'];
            else
                if rightAlignment==true
                    TableCorpus=[TableCorpus ' ' numij ' \\'];
                else
                    TableCorpus=[TableCorpus ' ' numij ' \cr'];
                end
            end
        end
    end
end

LastLine='';
for i=1:length(header)
    if isnan(footer(i))
        numi='-';
    else
        numi=num2str(footer(i));
    end
    if i<length(header)
        LastLine=[LastLine ' \bf ' numi ' &']; %#ok<*AGROW>
    else

        if rightAlignment==true
            LastLine=[LastLine  ' \bf ' numi ' \\ \end{array}$'];
        else
            LastLine=[LastLine  ' \bf ' numi '}$'];
        end
    end
end
str=[FirstLine TableCorpus LastLine];

if length(str)>1100
    % note that there is a maximum of string size of 1200 characters for the
    % LaTeX interpreter
    if rightAlignment==true
        startIndex=regexp(str,'\\');
    else
        startIndex=regexp(str,'\cr');
    end
    a=ceil(length(startIndex)/2);
    b=a+1;

    if rightAlignment==true
        toinsert=['\\ ' repmat('& ',1,length(header)-1) '\\ '];
    else
        toinsert=['\cr ' repmat('&',1,length(header)-1) ' '];
    end
    strbis=str;
    while length(strbis) >1100 % 1200
            strbis=[str(1:startIndex(a)-2) toinsert str(startIndex(b)-1:end)];
        a=a-1;
        b=b+1;
    end
    str=strbis;
end
end

