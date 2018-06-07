function References=extractRefs(OUT, FilesIncluded)
% Creates a table of references from input struct OUT

ri=size(OUT,1);
References=cell(5000,3);
j=1;
for ii=1:ri
    fileRefs=OUT{ii,1}.References;
    refsRows=size(fileRefs,1);
    for jj=1:refsRows
        References{j,1}=fileRefs(jj);
        References(j:j+refsRows,2)=FilesIncluded(ii,1);
        References(j:j+refsRows,3)=FilesIncluded(ii,9);
        
        j=j+1;
    end
end

References=References(1:j-1,:);
% end
end