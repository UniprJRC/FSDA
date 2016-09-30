function [docNode,docNodechr]=xmlcreateFS(FileName, varargin)
%create an XML file passing through publishFS
%
% Required input arguments:
%
%    FileName:     MATLAB File. String. Full or partial
%                  path of the MATLAB file for which Structured XML
%                  help has to be created
%                  Example-'myfile.m'
%
%
% Optional input arguments:
%
% write2file: Option to write XML file. Logical. Option which specifies
%             whether XML file must be created or if just structure docNode
%             must be created. The default value of write2file is true,
%             that is xml file is created in the following path
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


% out=publishFS('tclust','evalCode',false,'write2file',false);
out=publishFS(FileName,'evalCode',false,'write2file',false);

if ~isempty(out.description)
    out.description=removeExtraSpacesLF(out.description);
end

for j=2:4
    Inpi=out.InpArgs(:,j);
    for i=1:length(Inpi)
        if ~isempty(Inpi{i})
            Inpi{i}=removeExtraSpacesLF(Inpi{i});
        end
    end
    out.InpArgs(:,j)=Inpi;
end


for j=2:4
    Inpi=out.OptArgs(:,j);
    for i=1:length(Inpi)
        if ~isempty(Inpi{i})
            Inpi{i}=removeExtraSpacesLF(Inpi{i});
        end
    end
    out.OptArgs(:,j)=Inpi;
end

%%%%%%%%%%%%%%%%%
OutArgs=out.OutArgs;
for i=1:size(OutArgs,1)
    
    
    if ~cellfun(@isempty,OutArgs(i,2)) && ~isempty(strtrim(OutArgs{i,2}))
        % Short description
        for j=2:4
            Inpi=OutArgs{i,j};
            if ~isempty(Inpi)
                Inpi=removeExtraSpacesLF(Inpi);
            end
            out.OutArgs{i,j}=Inpi;
        end
    else
        
        % Table containing the details of the ith output arg which is a
        % structure
        outTable=OutArgs{i,5};
        for ii=1:size(outTable,1)
            % Short description
            DescTab=outTable{ii,2};
            outTable{ii,2}=removeExtraSpacesLF(DescTab);
        end
        out.OutArgs{i,5}=outTable;
    end
    
end

if ~isempty(out.MoreAbout)
    out.MoreAbout=removeExtraSpacesLF(out.MoreAbout);
end


Inpi=out.References;
for i=1:length(Inpi)
    if ~isempty(Inpi{i})
        Inpi{i}=removeExtraSpacesLF(Inpi{i});
    end
end
out.References=Inpi;

%% Clean examples part

for j=1:3
    
    Ex3=out.Ex(:,j);
    
    for i=1:length(Ex3)
        
        if j<=2
            str=removeExtraSpacesLF(Ex3{i});
        else
            % Given that this is code just remove extra consecutive line
            % feeds
            str=regexprep(Ex3{i},'\x0D','\x0A');
            str=regexprep(str,'\x0A*','\x0A');
            checkfirstLF=regexp(str,'\x0A');
            if checkfirstLF(1)==1
                str=str(2:end);
            end
        end
        if isempty(strtrim(str))
            out.Ex{i,j}=[];
        else
            PosLinBreaks = [regexp(str,'\x0A') length(str)];
            if j<=2
                CellStackedStrings=cell(length(PosLinBreaks),1);
                lP=length(PosLinBreaks);
            else
                CellStackedStrings=cell(length(PosLinBreaks)-1,1);
                lP=length(PosLinBreaks)-1;
            end
            
            for ii=1:lP
                if ii>1
                    strsel=str(PosLinBreaks(ii-1)+1:PosLinBreaks(ii));
                    % findLFinstrsel=regexp((strsel),'\n', 'once');
                    CellStackedStrings{ii}=strtrim(strsel);
                else
                    strsel=str(1:PosLinBreaks(ii));
                    CellStackedStrings{ii}=strtrim(strsel);
                end
            end
            
            % if j==1 given that the title is just one sentence we store
            % the content of the cell. In all the other cases we store the
            % cell itself
            if j==1
                out.Ex{i,j}=CellStackedStrings{:};
            else
                out.Ex{i,j}=CellStackedStrings; % str;
            end
        end
    end
end

%% Clean Extra examples part

for j=1:3
    
    Ex3=out.ExtraEx(:,j);
    
    for i=1:length(Ex3)
        if j<=2
            str=removeExtraSpacesLF(Ex3{i});
        else
            % Given that this is code just remove extra consecutive line
            % feeds
            str=regexprep(Ex3{i},'\x0D','\x0A');
            str=regexprep(str,'\x0A*','\x0A');
            checkfirstLF=regexp(str,'\x0A');
            if checkfirstLF(1)==1
                str=str(2:end);
            end
        end
        if isempty(strtrim(str))
            out.ExtraEx{i,j}=[];
        else
            
            PosLinBreaks = [regexp(str,'\x0A') length(str)];
            if j<=2
                CellStackedStrings=cell(length(PosLinBreaks),1);
                lP=length(PosLinBreaks);
            else
                CellStackedStrings=cell(length(PosLinBreaks)-1,1);
                lP=length(PosLinBreaks)-1;
            end
            
            for ii=1:lP
                if ii>1
                    strsel=str(PosLinBreaks(ii-1)+1:PosLinBreaks(ii));
                    % findLFinstrsel=regexp((strsel),'\n', 'once');
                    CellStackedStrings{ii}=strtrim(strsel);
                else
                    strsel=str(1:PosLinBreaks(ii));
                    CellStackedStrings{ii}=strtrim(strsel);
                end
            end
            
            % if j==1 given that the title is just one sentence we store
            % the content of the cell. In all the other cases we store the
            % cell itself
            if j==1
                out.ExtraEx{i,j}=CellStackedStrings{:};
            else
                out.ExtraEx{i,j}=CellStackedStrings; % str;
            end
        end
    end
end


% Now after cleaning out structure write to xml
[docNode,docNodechr]=xmlwriteFS(out,'write2file',write2file);

end

