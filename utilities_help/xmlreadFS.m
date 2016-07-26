function out=xmlreadFS(FileName)
% read xml file and converts it to a structure with dedicated fields

%% Beginning of code
[FSDAroot]=fileparts(which('docsearchFS.m'));

% Remove .m extension if it exists
[~,Name] = fileparts(FileName);

fsep=filesep;
fullpathFileName=[FSDAroot fsep 'helpfiles' fsep 'XML' fsep Name '.xml'];

dom=parseXML(fullpathFileName);

out=struct;
% Extract title
out.titl=dom.Children(2).Children(2).Data;
% Extract purpose
out.purpose=dom.Children(4).Children(2).Data;
% Extract description
if length(dom.Children(6).Children)>1
    out.description=dom.Children(6).Children(2).Data;
else
    out.description='';
end

% Extract InpArgs
% Create cell which will contain InpArgs
index=8;
InpArgs=cell((length(dom.Children(index).Children)-2)/2,8);
for i=1:size(InpArgs,1)
    for j=1:8
        %InpArgs{i,j}=dom.Children(index).Children(2*i+1).Children(2*j).Children.Data;
        
        if j==8 && size(dom.Children(index).Children(2*i+1).Children(2*j).Children,2)>1
            itemCell=dom.Children(index).Children(2*i+1).Children(2*j).Children;
            % In this case the fourth column of OutArgs is a cell which
            % contains name (1st column) and description (2nd column)
            lCellItem=(length(itemCell)-1)/2;
            CellItem=cell(lCellItem,2);
            for ii=1:lCellItem
                for jj=1:2
                    CellItem{ii,jj}=itemCell(ii*2).Children(jj*2).Children.Data;
                end
            end
            data=CellItem;
%         elseif size(dom.Children(index).Children(2*i+1).Children(2*j).Children,2)==0
        else
            data=dom.Children(index).Children(2*i+1).Children(2*j).Children.Data;
        end
        
        if ~isempty(strtrim(data))
            InpArgs{i,j}=data;
        end
    end
end
out.InpArgs=InpArgs;


% Extract OptArgs
% Create cell which will contain OptArgs
index=10;
OptArgs=cell((length(dom.Children(index).Children)-2)/2,7);
for i=1:size(OptArgs,1)
    for j=1:7
        if j==7 && size(dom.Children(index).Children(2*i+1).Children(2*j).Children,2)>1
            itemCell=dom.Children(index).Children(2*i+1).Children(2*j).Children;
            % In this case the fourth column of OutArgs is a cell which
            % contains name (1st column) and description (2nd column)
            lCellItem=(length(itemCell)-1)/2;
            CellItem=cell(lCellItem,2);
            for ii=1:lCellItem
                for jj=1:2
                    CellItem{ii,jj}=itemCell(ii*2).Children(jj*2).Children.Data;
                end
            end
            data=CellItem;
        else
            data=dom.Children(index).Children(2*i+1).Children(2*j).Children.Data;
        end
        if ~isempty(strtrim(data))
            OptArgs{i,j}=data;
        end
        
    end
end
out.OptArgs=OptArgs;

% Extract OutArgs
% Create cell which will contain OutArgs
index=12;
OutArgs=cell((length(dom.Children(index).Children)-2)/2,4);
for i=1:size(OutArgs,1)
    for j=1:4
        if j==4 && isempty(OutArgs{i,2}) &&  isempty(OutArgs{i,3})
            itemCell=dom.Children(index).Children(2*i+1).Children(2*j).Children;
            % In this case the fourth column of OutArgs is a cell which
            % contains name (1st column) and description (2nd column)
            lCellItem=(length(itemCell)-1)/2;
            CellItem=cell(lCellItem,2);
            for ii=1:lCellItem
                for jj=1:2
                    CellItem{ii,jj}=itemCell(ii*2).Children(jj*2).Children.Data;
                end
            end
            data=CellItem;
        else
            data=dom.Children(index).Children(2*i+1).Children(2*j).Children.Data;
        end
        if ~isempty(strtrim(data))
            OutArgs{i,j}=data;
        end
    end
end
out.OutArgs=OutArgs;

% Extract MoreAbout
if length(dom.Children(14).Children)>1
    out.MoreAbout=dom.Children(14).Children(2).Data;
else
    out.MoreAbout='';
end

% Extract Acknowledgements
if length(dom.Children(16).Children)>1
    out.Acknowledgements=dom.Children(16).Children(2).Data;
else
    out.Acknowledgements='';
end

% Extract References
lReferences=(length(dom.Children(18).Children)-2)/2;
References=cell(lReferences,1);
for i=1:lReferences
    References{i}=dom.Children(18).Children(2*i+1).Children.Data;
end
out.References=References;

% Extract SeeAlso
lSeeAlso=(length(dom.Children(20).Children)-2)/2;
SeeAlso=cell(lSeeAlso,1);
for i=1:lSeeAlso
    SeeAlso{i}=dom.Children(20).Children(2*i+1).Children.Data;
end
out.SeeAlso=SeeAlso;

% Extract Ex
% Create cell which will contain Examples
index=22;
Ex=cell((length(dom.Children(index).Children)-2)/2,4);
for i=1:size(Ex,1)
    for j=1:4
        if ~isempty(dom.Children(index).Children(2*i+1).Children(2*j).Children)
            Ex{i,j}=dom.Children(index).Children(2*i+1).Children(2*j).Children.Data;
        end
        
    end
end
out.Ex=Ex;

% Extract ExtraEx
% Create cell which will contain ExtraExamples
indExtraEx=24;
ExtraEx=cell((length(dom.Children(indExtraEx).Children)-2)/2,4);
for i=1:size(ExtraEx,1)
    for j=1:4
        if ~isempty(dom.Children(indExtraEx).Children(2*i+1).Children(2*j).Children)
            ExtraEx{i,j}=dom.Children(indExtraEx).Children(2*i+1).Children(2*j).Children.Data;
        end
    end
end
out.ExtraEx=ExtraEx;
end