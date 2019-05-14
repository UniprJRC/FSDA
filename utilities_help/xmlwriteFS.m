function [docNode,docNodechr]=xmlwriteFS(out, varargin)
%create an XML file starting from input structure
%
% Required input arguments:
%
%    out:     structure created by publishFS or by GUI which contains
%                   the  following fields
%                   out.titl =........
%                   out.purpose =........
%
%
% Optional input arguments:
%
% write2file: Option to write HTML file. Logical. Option which specifies
%             whether XML file must be created or if just structure docNode
%             must be created. The default value of write2file is true,
%             that is xml file is overwritten in the following path
%             (main root of FSDA) filesep helpfiles filesep XML
%             where filesep means "\" or "/" according to the operating
%             system
%             Example - 'write2file','false'
%
% Output:
%
%    docNode:    Document Object Model node, as defined by the World Wide Web
%               consortium.
% docNodechr:    Character vector that contains the serialized DOM node as
%               it appears in an XML file.
%
% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit


%% Beginning of code
write2file=true;

if nargin>1
    options=struct('write2file',write2file);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:regressB:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
        
        % Write in structure 'options' the options chosen by the user
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
        
    end
    
    write2file=options.write2file;
end

%% Beginning of code
StructArgsNames={'Value' 'Description'};

docNode = com.mathworks.xml.XMLUtils.createDocument('HelpXML');

% Create file name section
Preamble=docNode.createElement('Title');
Preamble.appendChild(docNode.createComment('This is simply the filename'));
Preambletxt=docNode.createTextNode(out.titl);
Preamble.appendChild(Preambletxt);
docNode.getDocumentElement.appendChild(Preamble);

% Create purpose section
Preamble=docNode.createElement('Purpose');
Preamble.appendChild(docNode.createComment('This is the second line of the .m file'));
if iscell(out.purpose)
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.purpose{:}));
else
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.purpose));
end

Preamble.appendChild(Preambletxt);
docNode.getDocumentElement.appendChild(Preamble);

% Create description section
Preamble=docNode.createElement('Description');
Preamble.appendChild(docNode.createComment('Description section'));
% Preambletxt=docNode.createTextNode(out.description{:});
if iscell(out.description)
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.description{:}));
else
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.description));
end

Preamble.appendChild(Preambletxt);
docNode.getDocumentElement.appendChild(Preamble);
% [{'\n';'\n'} out.description]

%Create InpArgs section
InpArgsNames={'Name' 'ShortDesc' 'TypeInd' 'LongDesc' 'Example' 'DataType' 'ReqArg' 'Struct'};
InpArgs = docNode.createElement('InpArgs');
InpArgs.appendChild(docNode.createComment('REQUIRED INPUT ARGUMENT SECTION'));

for i=1:size(out.InpArgs,1)
    entry_node = docNode.createElement('Item');
    for j=1:length(InpArgsNames)
        name_node = docNode.createElement(InpArgsNames{j});
        if ~isempty(out.InpArgs{i,j})
            if iscell(out.InpArgs{i,j})
                InpArgsCell=out.InpArgs{i,j};
                for ii=1:size(InpArgsCell,1)
                    if ~isempty(InpArgsCell{ii,1}) || ~isempty(InpArgsCell{ii,2})
                        entry_nodeCell = docNode.createElement('ItemCell');
                        for jj=1:2
                            name_nodeCell = docNode.createElement(StructArgsNames{jj});
                            name_nodeCell.appendChild(docNode.createTextNode((InpArgsCell{ii,jj})));
                            entry_nodeCell.appendChild(name_nodeCell);
                        end
                        name_node.appendChild(entry_nodeCell);
                    end
                end
            else
                name_node.appendChild(docNode.createTextNode((out.InpArgs{i,j})));
            end
        else
            name_node.appendChild(docNode.createTextNode(' '));
        end
        
        
        entry_node.appendChild(name_node);
        InpArgs.appendChild(entry_node);
    end
end
docNode.getDocumentElement.appendChild(InpArgs);

%Create OptArgs section
OptArgsNames={'Name' 'ShortDesc' 'TypeInd' 'LongDesc' 'Example' 'DataType' 'Struct'};
OptArgs = docNode.createElement('OptArgs');
OptArgs.appendChild(docNode.createComment('OPTIONAL (NAME/PAIRS) INPUT ARGUMENT SECTION'));
for i=1:size(out.OptArgs,1)
    entry_node = docNode.createElement('Item');
    for j=1:length(OptArgsNames)
        name_node = docNode.createElement(OptArgsNames{j});
        if ~isempty(out.OptArgs{i,j})
            if iscell(out.OptArgs{i,j})
                OptArgsCell=out.OptArgs{i,j};
                for ii=1:size(OptArgsCell,1)
                    if ~isempty(OptArgsCell{ii,1}) || ~isempty(OptArgsCell{ii,2})
                        entry_nodeCell = docNode.createElement('ItemCell');
                        for jj=1:2
                            name_nodeCell = docNode.createElement(StructArgsNames{jj});
                            name_nodeCell.appendChild(docNode.createTextNode((OptArgsCell{ii,jj})));
                            entry_nodeCell.appendChild(name_nodeCell);
                        end
                        name_node.appendChild(entry_nodeCell);
                    end
                end
            else
                name_node.appendChild(docNode.createTextNode((out.OptArgs{i,j})));
            end
        else
            name_node.appendChild(docNode.createTextNode(' '));
        end
        
        entry_node.appendChild(name_node);
        OptArgs.appendChild(entry_node);
    end
end
docNode.getDocumentElement.appendChild(OptArgs);

%Create OutArgs section
OutArgsNames={'Name' 'ShortDesc' 'TypeInd' 'LongDesc' 'Structure'};

OutArgs = docNode.createElement('OutArgs');
OutArgs.appendChild(docNode.createComment('OUTPUT ARGUMENT SECTION'));

for i=1:size(out.OutArgs,1)
    entry_node = docNode.createElement('Item');
    for j=1:length(OutArgsNames)
        name_node = docNode.createElement(OutArgsNames{j});
        if ~isempty(out.OutArgs{i,j})
            if iscell(out.OutArgs{i,j})
                OutArgsCell=out.OutArgs{i,j};
                for ii=1:size(OutArgsCell,1)
                    entry_nodeCell = docNode.createElement('ItemCell');
                    for jj=1:2
                        name_nodeCell = docNode.createElement(StructArgsNames{jj});
                        name_nodeCell.appendChild(docNode.createTextNode((OutArgsCell{ii,jj})));
                        entry_nodeCell.appendChild(name_nodeCell);
                    end
                    name_node.appendChild(entry_nodeCell);
                end
            else
                name_node.appendChild(docNode.createTextNode((out.OutArgs{i,j})));
            end
        else
            name_node.appendChild(docNode.createTextNode(' '));
        end
        entry_node.appendChild(name_node);
        OutArgs.appendChild(entry_node);
    end
end
docNode.getDocumentElement.appendChild(OutArgs);

%% Create More About section
Preamble=docNode.createElement('MoreAbout');
Preamble.appendChild(docNode.createComment('MORE ABOUT SECTION'));
if iscell(out.MoreAbout)
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.MoreAbout{:}));
else
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.MoreAbout));
end
Preamble.appendChild(Preambletxt);
docNode.getDocumentElement.appendChild(Preamble);

%% Create Acknowledgements section
Preamble=docNode.createElement('Acknowledgements');
Preamble.appendChild(docNode.createComment('ACKNOWLEDGEMENTS SECTION'));
if iscell(out.Acknowledgements)
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.Acknowledgements{:}));
else
    Preambletxt=docNode.createTextNode(sprintf('%s\n',out.Acknowledgements));
end

Preamble.appendChild(Preambletxt);
docNode.getDocumentElement.appendChild(Preamble);

%% Create References section
References = docNode.createElement('References');
References.appendChild(docNode.createComment('REFERENCES SECTION'));
for i=1:length(out.References)
    entry_node = docNode.createElement('Item');
    entry_node.appendChild(docNode.createTextNode((out.References{i})));
    References.appendChild(entry_node);
end
docNode.getDocumentElement.appendChild(References);

%% Create SeeAlso section
SeeAlso = docNode.createElement('SeeAlso');
SeeAlso.appendChild(docNode.createComment('SEE ALSO SECTION'));

for i=1:length(out.SeeAlso)
    if ~isempty(out.SeeAlso{i})
        entry_node = docNode.createElement('Item');
        entry_node.appendChild(docNode.createTextNode((out.SeeAlso{i})));
        SeeAlso.appendChild(entry_node);
    end
end
docNode.getDocumentElement.appendChild(SeeAlso);


%% Create Example section
ExNames={'Title' 'Desc' 'MATLABcode' 'Exec'};
Ex = docNode.createElement('Ex');
Ex.appendChild(docNode.createComment('EXAMPLES SECTION'));
for i=1:size(out.Ex,1)
    entry_node = docNode.createElement('Item');
    for j=1:length(ExNames)
        name_node = docNode.createElement(ExNames{j});
        if ~isempty(out.Ex{i,j})
            if j==4
                name_node.appendChild(docNode.createTextNode(num2str(out.Ex{i,j})));
                
            else
                %                 if iscell(out.Ex{i,j}) && j~=1
                %                     name_node.appendChild(docNode.createTextNode(sprintf('%s\n',(out.Ex{i,j}{:}))));
                %                 elseif j~=1
                %                     name_node.appendChild(docNode.createTextNode(sprintf('%s\n',out.Ex{i,j})));
                %                 else
                %                     name_node.appendChild(docNode.createTextNode(out.Ex{i,j}));
                %                 end
                if iscell(out.Ex{i,j})
                    ExCell=out.Ex{i,j};
                    
                    for ii=1:size(ExCell,1)
                        entry_nodeCell = docNode.createElement('ItemCell');
                        try
                        entry_nodeCell.appendChild(docNode.createTextNode(deblank(ExCell{ii})));
                        catch
                            ddd=1;
                        end
                        name_node.appendChild(entry_nodeCell);
                    end
                else
                    name_node.appendChild(docNode.createTextNode(out.Ex{i,j}));
                end
            end
        else
            name_node.appendChild(docNode.createTextNode(' '));
        end
        entry_node.appendChild(name_node);
        Ex.appendChild(entry_node);
    end
end
docNode.getDocumentElement.appendChild(Ex);

%% Create Extra Example section
ExtraExNames={'Title' 'Desc' 'MATLABcode' 'Exec'};
ExtraEx = docNode.createElement('ExtraEx');
ExtraEx.appendChild(docNode.createComment('EXTRA EXAMPLES SECTION'));
for i=1:size(out.ExtraEx,1)
    entry_node = docNode.createElement('Item');
    for j=1:length(ExtraExNames)
        name_node = docNode.createElement(ExtraExNames{j});
        if ~isempty(out.ExtraEx{i,j})
            if j==4
                name_node.appendChild(docNode.createTextNode(num2str(out.ExtraEx{i,j})));
            else
                %                 if iscell(out.ExtraEx{i,j}) && j~=1
                %                     name_node.appendChild(docNode.createTextNode(sprintf('%s\n',(out.ExtraEx{i,j}{:}))));
                %                 elseif j~=1
                %                     name_node.appendChild(docNode.createTextNode(sprintf('%s\n',out.ExtraEx{i,j})));
                %                 else
                %                     name_node.appendChild(docNode.createTextNode(out.ExtraEx{i,j}));
                %                 end
                if iscell(out.ExtraEx{i,j})
                    ExCell=out.ExtraEx{i,j};
                    for ii=1:size(ExCell,1)
                        entry_nodeCell = docNode.createElement('ItemCell');
                        try
                            entry_nodeCell.appendChild(docNode.createTextNode(deblank(ExCell{ii})));
                        catch
                            ddd=1;
                        end
                        name_node.appendChild(entry_nodeCell);
                    end
                else
                    name_node.appendChild(docNode.createTextNode(out.ExtraEx{i,j}));
                end
                
                % name_node.appendChild(docNode.createTextNode(out.ExtraEx{i,j}));
            end
        else
            name_node.appendChild(docNode.createTextNode(' '));
        end
        entry_node.appendChild(name_node);
        ExtraEx.appendChild(entry_node);
    end
end
docNode.getDocumentElement.appendChild(ExtraEx);

[FSDAroot]=fileparts(which('docsearchFS.m'));

Name=out.titl;

fsep=filesep;
fullpathFileName=[FSDAroot fsep 'helpfiles' fsep 'XML' fsep strtrim(Name) '.xml'];
if write2file
    xmlwrite(fullpathFileName,docNode);
end
docNodechr=xmlwrite(docNode);
end

% http://stackoverflow.com/questions/35510368/python-remove-xd-from-xml
