function [outstring]=mwriteFS(out)
%create an XML file starting from input structure
%
% Required input arguments:
%
%    out:     structure created by publishFS or by GUI which contains
%

startpos=20;
startcolumn=22;
startcolumn1=25;
endcolumn=75;
startcolumnEx=6;
endcolumnEx=1000;
comment=struct;
comment.commentsign='%';
comment.startcolumn=startcolumnEx;
widthEx=65;

[FSDAroot]=fileparts(which('docsearchFS.m'));
fsep=filesep;
fullpathFileName=[FSDAroot fsep 'helpfiles' fsep 'XML' fsep 'tmpHELP.m'];

% load outtowritemp
% fileID = fopen('tmpHELP.m','w');
fileID = fopen(fullpathFileName,'w');

fprintf(fileID,'%s',['%' out.purpose{:}]);
% fileID = fopen('tmpHELP.m','a');
fprintf(fileID,'\n%s','%');
fprintf(fileID,'\n%s',['%<a href="matlab: docsearchFS(''' out.titl ''')">Link to the help function</a>']);
fprintf(fileID,'\n%s','%');
% fprintf(out.description)

fprintf(fileID,'\n');
for i=1:size(out.description,1)
    % fprintf(fileID,'\n%s','% ');
    % fprintf(fileID,'\n');
    if ~isempty(out.description{i})
        descri=out.description{i};
        descriFormatted=wraptextFS(descri,'startcolumn',3,'comment',true);
        % fprintf(fileID, '%s',[descriFormatted,'\n']);
        fprintf(fileID,'%s',descriFormatted);
    else
        fprintf(fileID,'%%\n');
    end
end
fprintf(fileID,'%%\n');
fprintf(fileID,'%s\n','%  Required Input Arguments:');
fprintf(fileID,'%%\n');

InpArgs=out.InpArgs;
for i=1:size(InpArgs,1)
    
    nameInpArgs=InpArgs{i,1};
    
    % Short description
    ShortDesc=strtrim(InpArgs{i,2});
    descri=[repmat(' ',1,startpos-2-length(nameInpArgs)) nameInpArgs ': ' ShortDesc];
    descriFormatted=wraptextFS(descri,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',true,'comment',true);
    fprintf(fileID,'%s',descriFormatted);
    
    % type indication (Scalar, matrix, ...)
    ShortDesc=strtrim(InpArgs{i,3});
    descriFormatted=wraptextFS(ShortDesc,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
    fprintf(fileID,'%s',descriFormatted);
    
    % Long description
    if ~cellfun(@isempty,InpArgs(i,4)) && ~isempty(strtrim(InpArgs{i,4}))
        LongDesc=strtrim(InpArgs{i,4});
        descriFormatted=wraptextFS(LongDesc,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
        fprintf(fileID,'%s',descriFormatted);
    end
    
    % Data type
    if ~cellfun(@isempty,InpArgs(i,6)) && ~isempty(strtrim(InpArgs{i,6}))
        DataType=InpArgs{i,6};
        
        ntimes=20;
        leftMargin= {[repmat(' ',1,ntimes) 'Data Types - ']};
        % newstr=strjoin(strcat('%%',leftMargin,strtrim(DataType),'\n'),'');
        
        newstr=strjoin(strcat('%',leftMargin,strtrim(DataType),'0A123'),'');
        newstr=regexprep(newstr,'0A123','\x0A');
        
        fprintf(fileID,'%s',newstr);
    end
end

fprintf(fileID,'%%\n');
fprintf(fileID,'%s\n','%  Optional input arguments:');
fprintf(fileID,'%%\n');


OptArgs=out.OptArgs;
for i=1:size(OptArgs,1)
    
    nameInpArgs=OptArgs{i,1};
    
    % Short description
    ShortDesc=OptArgs{i,2};
    descri=[repmat(' ',1,startpos-2-length(nameInpArgs)) nameInpArgs ': ' ShortDesc];
    descriFormatted=wraptextFS(descri,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',true,'comment',true);
    
    fprintf(fileID,'%s',descriFormatted);
    
    % type indication (Scalar, matrix, ...)
    ShortDesc=strtrim(OptArgs{i,3});
    %    descriFormatted=wraptext2(ShortDesc,numWrappedCols,0);
    descriFormatted=wraptextFS(ShortDesc,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
    
    fprintf(fileID,'%s',descriFormatted);
    
    % Long description
    if ~cellfun(@isempty,OptArgs(i,4)) && ~isempty(strtrim(OptArgs{i,4}))
        LongDesc=strtrim(OptArgs{i,4});
        % Check if string ReMarK: is present inside LongDesc
        CheckRemarkString=regexpi(LongDesc,'Remark\:');
        if ~isempty(CheckRemarkString)
            AddRemark=LongDesc(CheckRemarkString(1):end);
            LongDesc=LongDesc(1:CheckRemarkString(1)-1);
        else
            AddRemark='';
        end
        
        % descriFormatted=wraptext2(LongDesc,numWrappedCols,0);
        descriFormatted=wraptextFS(LongDesc,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
        fprintf(fileID,'%s',descriFormatted);
        
        
    end
    
    % Examples
    if ~cellfun(@isempty,OptArgs(i,5)) && ~isempty(strtrim(OptArgs{i,5}))
        Example=OptArgs{i,5};
        
        leftMargin= {[repmat(' ',1,startpos) 'Example - ']};
        
        newstr=strjoin(strcat('%',leftMargin,strtrim(Example),'0A123'),'');
        newstr=regexprep(newstr,'0A123','\x0A');
        
        %  newstr=strjoin(strcat('%%',leftMargin,strtrim(DataType),'\n'),'');
        fprintf(fileID,'%s',newstr);
    end
    
    % Data type
    if ~cellfun(@isempty,OptArgs(i,6)) && ~isempty(strtrim(OptArgs{i,6}))
        DataType=OptArgs{i,6};
        
        leftMargin= {[repmat(' ',1,startpos) 'Data Types - ']};
        
        newstr=strjoin(strcat('%',leftMargin,strtrim(DataType),'0A123'),'');
        newstr=regexprep(newstr,'0A123','\x0A');
        
        %  newstr=strjoin(strcat('%%',leftMargin,strtrim(DataType),'\n'),'');
        fprintf(fileID,'%s',newstr);
    end
    
    if ~isempty(AddRemark)
        fprintf(fileID,'%%\n');
        AddRemarkFormatted=wraptextFS(AddRemark,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
        fprintf(fileID,'%s',AddRemarkFormatted);
        fprintf(fileID,'%%\n');
    end
    
end

%% Output arguments section
fprintf(fileID,'%%\n');
fprintf(fileID,'%s\n','%   Output:');
fprintf(fileID,'%%\n');

OutArgs=out.OutArgs;
StartColOutArg=19;

for i=1:size(OutArgs,1)
    
    
    if ~cellfun(@isempty,OutArgs(i,2)) && ~isempty(strtrim(OutArgs{i,2}))
        
        if strcmp(OutArgs{i,1},'varargout')
            fprintf(fileID,'%%\n');
            fprintf(fileID,'%s\n','%    Optional Output:');
            fprintf(fileID,'%%\n');
            
            % type indication (Scalar, matrix, ...)
            ShortDesc=OutArgs{i,3};
            descriFormatted=wraptextFS(ShortDesc,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
            
            % descriFormatted=wraptext2(ShortDesc,numWrappedCols,0);
            fprintf(fileID,'%s',descriFormatted);
            
            % Long description
            if ~cellfun(@isempty,OutArgs(i,4)) && ~isempty(strtrim(OutArgs{i,4}))
                LongDesc=OutArgs{i,4};
                % descriFormatted=wraptext2(LongDesc,numWrappedCols,0);
                descriFormatted=wraptextFS(LongDesc,'startcolumn',startcolumn1,'endcolumn',endcolumn,'firstline',false,'comment',true);
                
                fprintf(fileID,'%s',descriFormatted);
            end
        else
            nameOutArgs=OutArgs{i,1};
            % Short description
            ShortDesc=OutArgs{i,3};
            descri=[repmat(' ',1,startpos-2-length(nameOutArgs)) nameOutArgs ': ' ShortDesc];
            % descriFormatted=wraptext2(descri,numWrappedCols);
            descriFormatted=wraptextFS(descri,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',true,'comment',true);
            
            fprintf(fileID,'%s',descriFormatted);
            
            % type indication (Scalar, matrix, ...)
            ShortDesc=OutArgs{i,2};
            % Make sure that the last character of ShortDesc is the full
            % stop and add it it is not present add at the end
            if ~isempty(ShortDesc)
                if regexp(ShortDesc,'\x0A')==length(ShortDesc)
                    ShortDesc=deblank(ShortDesc(1:end-1));
                end
                if ~strcmp(ShortDesc(end),'\.')
                    ShortDesc=[ShortDesc '.']; %#ok<AGROW>
                end
            end
            
            % descriFormatted=wraptext2(ShortDesc,numWrappedCols,0);
            descriFormatted=wraptextFS(ShortDesc,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
            
            fprintf(fileID,'%s',descriFormatted);
            
            % Long description
            if ~cellfun(@isempty,OutArgs(i,4)) && ~isempty(strtrim(OutArgs{i,4}))
                LongDesc=OutArgs{i,4};
                % descriFormatted=wraptext2(LongDesc,numWrappedCols,0);
                descriFormatted=wraptextFS(LongDesc,'startcolumn',startcolumn,'endcolumn',endcolumn,'firstline',false,'comment',true);
                
                fprintf(fileID,'%s',descriFormatted);
            end
        end
    else
        fprintf(fileID,'%%\n');
        nameOutArgs=OutArgs{i,1};
        spacestoadd=repmat(' ',1,StartColOutArg-length(nameOutArgs));
        fprintf(fileID,'%s',['%' spacestoadd nameOutArgs ': structure which contains the following fields']);
        fprintf(fileID,'\n%%\n');
        
        outTable=OutArgs{i,5};
        for ii=1:size(outTable,1)
            nameTab=outTable{ii,1};
            % Short description
            DescTab=outTable{ii,2};
            descri=[repmat(' ',1,StartColOutArg-2-length(nameTab)-length(nameOutArgs)) nameOutArgs '.' nameTab '= ' DescTab];
            % descriFormatted=wraptext2(descri,75-StartColOutArg-1);
            descriFormatted=wraptextFS(descri,'startcolumn',startcolumn1,'endcolumn',endcolumn,'firstline',true,'comment',true);
            
            fprintf(fileID,'%s',descriFormatted);
        end
    end
end

%% More About Section
if ~isempty(out.MoreAbout)
    fprintf(fileID,'%%\n');
    fprintf(fileID,'%s\n','%   More About:');
    fprintf(fileID,'%%\n');
    for i=1:size(out.MoreAbout,1)
        % fprintf(fileID,'\n%s','% ');
        % fprintf(fileID,'\n');
        if ~isempty(out.MoreAbout{i})
            descri=out.MoreAbout{i};
            % descriFormatted=wraptext1(descri,75);
            descriFormatted=wraptextFS(descri,'startcolumn',3,'comment',true);
            
            % fprintf(fileID, '%s',[descriFormatted,'\n']);
            fprintf(fileID,'%s',descriFormatted);
        else
            fprintf(fileID,'%%\n');
        end
    end
end
%% See also section
fprintf(fileID,'%%\n');
fprintf(fileID,'%s','%  See also: ');

newstr=strjoin(out.SeeAlso,', ');

fprintf(fileID,'%s',newstr);
fprintf(fileID,'\n%%\n');

%% References Section
if ~isempty(out.References)
    fprintf(fileID,'%%\n');
    fprintf(fileID,'%s\n','%   References:');
    fprintf(fileID,'%%\n');
    for i=1:size(out.References,1)
        % fprintf(fileID,'\n%s','% ');
        % fprintf(fileID,'\n');
        if ~isempty(out.References{i})
            descri=out.References{i};
            % descriFormatted=wraptext1(descri,75);
            descriFormatted=wraptextFS(descri,'startcolumn',3,'comment',true);
            
            % fprintf(fileID, '%s',[descriFormatted,'\n']);
            fprintf(fileID,'%s',descriFormatted);
        else
            fprintf(fileID,'%%\n');
        end
    end
end
out.LastModified='05-06-2016';
fprintf(fileID,'%%\n');
fprintf(fileID,'%%\n');
fprintf(fileID,'%s\n','%  Copyright 2008-2016.');
fprintf(fileID,'%s\n','%  Written by FSDA team');
fprintf(fileID,'%%\n');
fprintf(fileID,'%s',['%<a href="matlab: docsearchFS(''' out.titl ''')">Link to the help function</a>']);
fprintf(fileID,'\n%s\n',['% Last Modified ' out.LastModified]);
fprintf(fileID,'\n');
fprintf(fileID,'%s\n','% Examples:');
fprintf(fileID,'\n');

Ex=out.Ex;

for i=1:size(Ex,1)
    fprintf(fileID,'%s\n','%{');
    
    ExText=strtrim(Ex{i,1});
    if str2double(Ex{i,4})==1
        Commentsign='%% ';
    else
        Commentsign='% ';
    end
    
    descriFormatted=wraptextFS([Commentsign ExText],'startcolumn',startcolumnEx,'endcolumn',endcolumnEx,'firstline',false,'comment',false);
    fprintf(fileID,'%s',descriFormatted);
    
    for jj=2:3
        if ~isempty(Ex{i,jj})
            Exi=strtrim(Ex{i,jj});
            for ii=1:size(Exi,1)
                if jj==2
                    descriFormatted=wraptextFS( Exi{ii,1},'startcolumn',startcolumnEx+3,'endcolumn',endcolumn,'firstline',false,'comment',comment);
                else
                    % if j==3 we must check whether it is comment or not
                    % it is a comment if the first character is symbol %
                    Exii=Exi{ii,1};
                    if ~isempty(Exii)
                        if strcmp(Exii(1),'%')
                            % descriFormatted=wraptextFS( Exii(2:end),'startcolumn',startcolumnEx,'endcolumn',endcolumn,'firstline',false,'comment',true,'code',false);
                            descriFormatted=wraptextFS( Exii(2:end),'startcolumn',startcolumnEx,'endcolumn',endcolumn,'firstline',false,'comment',comment);
                        else
                            descriFormatted=wraptextFS( Exii,'startcolumn',startcolumnEx,'endcolumn',endcolumnEx,'firstline',false,'comment',false,'code',true);
                        end
                    end
                end
                fprintf(fileID,'%s',descriFormatted);
            end
        end
    end
    fprintf(fileID,'%s\n\n','%}');
end

Ex=out.ExtraEx;

for i=1:size(Ex,1)
    fprintf(fileID,'%s\n','%{');
    
    ExText=strtrim(Ex{i,1});
    if str2double(Ex{i,4})==1
        Commentsign='%% ';
    else
        Commentsign='% ';
    end
    
    descriFormatted=wraptextFS([Commentsign ExText],'startcolumn',startcolumnEx,'endcolumn',endcolumnEx,'firstline',false,'comment',false);
    fprintf(fileID,'%s',descriFormatted);
    
    for jj=2:3
        if ~isempty(Ex{i,jj})
            Exi=strtrim(Ex{i,jj});
            for ii=1:size(Exi,1)
                if jj==2
                    descriFormatted=wraptextFS( Exi{ii,1},'startcolumn',startcolumnEx+3,'endcolumn',endcolumn,'firstline',false,'comment',comment);
                else
                    % if j==3 we must check whether it is comment or not
                    % it is a comment if the first character is symbol %
                    
                    Exii=Exi{ii,1};
                    if ~isempty(Exii)
                        if strcmp(Exii(1),'%')
                            % descriFormatted=wraptextFS( Exii(2:end),'startcolumn',startcolumnEx,'endcolumn',endcolumn,'firstline',false,'comment',true,'code',false);
                            descriFormatted=wraptextFS( Exii(2:end),'startcolumn',startcolumnEx,'endcolumn',endcolumn,'firstline',false,'comment',comment);
                        else
                            descriFormatted=wraptextFS( Exii,'startcolumn',startcolumnEx,'endcolumn',endcolumnEx,'firstline',false,'comment',false,'code',true);
                        end
                    end
                end
                fprintf(fileID,'%s',descriFormatted);
            end
        end
    end
    fprintf(fileID,'%s\n\n','%}');
end

fclose(fileID);
end
