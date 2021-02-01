function outDiff=getWebMatlabFiles(outputOFHtmlHelpFileWeb)
% Creates a list of html files which are inside outputOFHtmlHelpFileWeb  
% Copyright 2008-2021.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit


s=ls(outputOFHtmlHelpFileWeb);
ri=size(s,1);
ListofFiles2=cell(1,ri-2);
j=1;
for i=3:ri
    tmp=strsplit(strip(s(i,1:end)),'.');
    if size(tmp,2) >= 2
        if (strcmp(tmp{1,2},'html'))
            % disp([num2str(i), '\t']);
            ListofFiles2{1,j}=[tmp{1,1} '.html'];
            j=j+1;
        end
    end
end
outDiff=ListofFiles2(1:j-1);
end
