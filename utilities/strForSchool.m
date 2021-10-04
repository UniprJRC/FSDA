function [str] = strForSchool(header, corpus, footer, classes)
[n,p]=size(corpus);

% classes is the string array to put in the second column of the table
if nargin<4
    classes='';
end

FirstLine='$\matrix{ $\boldmath{$';
for i=1:length(header)
    if i<length(header)
        FirstLine=[FirstLine header{i} '$}$ & $\boldmath{$'];
    else
        FirstLine=[FirstLine header{i} '$}$ \cr'];
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
                TableCorpus=[TableCorpus ' ' numij ' \cr'];
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
                TableCorpus=[TableCorpus ' ' numij ' \cr'];
            end
        end
    end
    
    % Insert simbol ... for the lines which are not shown
    for j=1:p
        numij='...';
        if j<p
            TableCorpus=[TableCorpus ' ' numij ' &'];
        else
            TableCorpus=[TableCorpus ' ' numij ' \cr'];
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
                TableCorpus=[TableCorpus ' ' numij ' \cr'];
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
        LastLine=[LastLine  ' \bf ' numi '}$'];
    end
end
str=[FirstLine TableCorpus LastLine];

end

