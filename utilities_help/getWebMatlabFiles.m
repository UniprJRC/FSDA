function outDiff=getWebMatlabFiles()

s=ls('C:\FSDA\helpfiles\FSDA');
[ri,co]=size(s);
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
