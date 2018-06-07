function CreateFSDAbibliography()
% CreateFSDAbibliography creates the bibliograpy html page
%
%
% The purpose of this function is to automatize the creation of the
% bibliograpy page contained in the HTML documentation
% 
%
%
% Required input arguments:
%
%  Output:
%
%         bibliography.html:   file which contains all the bibliograpy
%         entries contained in the "References:" section of MATLAB source files 
%
%          
%
%
% Copyright 2008-2018.
% Written by FSDA team
%


%% STEP 1 create personalized contents files
FileName='addFSDA2path';
FullPath=which(FileName);
%Navigate to the main folder of FSDA
FSDAroot=fileparts(FullPath);
cd(FSDAroot)
% Specify subfolders of main folders for which contents file has to be
% created
InclDir={'graphics' 'regression' 'multivariate' 'clustering' 'combinatorial' ...
    'examples' 'utilities' 'utilities_stat' 'utilities_help'};
ExclDir={'privateFS'  'datasets'};
% Create list of folders which must have the personalized contents file
list = findDir(FSDAroot,'InclDir',InclDir,'ExclDir',ExclDir);
% Crete personalized contents file for main folder of FSDA
% and required subfolders.
[FilesIncluded,FilesExcluded]=makecontentsfileFS('dirpath',list, ...
    'FilterFileContent','%FScategory:','force',false,'printOutputCell','Contents.m');
disp('List of files which have been excluded (with path)')
disp(FilesExcluded(:,[1 9]))

%% STEP 2: create HTML for all files filtered using makecontentsFilesFS
[~,OUT]=publishFSallFiles(FilesIncluded, 'evalCode','false');

cits=extractRefs(OUT, FilesIncluded);

citst=sortrows(cell2table(cits));
%% vvvv
[~,ia,~]=unique(cits(:,1));
uniquesTab=cell2table(cits(ia,:));

%% STEP 3: create 'bibliography.html' file for all uniques references inside all 
%  MATLAB source files.

html='';

for ii = 2:size(uniquesTab,1)
    
    refEntry=char(uniquesTab{ii,1});
    refText=refEntry;
    
    % write in italic what is inside double quotes
    dbquotesVect=regexp(refText,'\"');
    if (~isempty(dbquotesVect))
        tmpText=[ refText(1:dbquotesVect(1)-1) '<em>' ...
            refText(dbquotesVect(1)+1:dbquotesVect(2)-1) '</em>' ...
            refText(dbquotesVect(2)+1:end)];
        refText=tmpText;
    end
    
    % delete every citation that is inside square brackets
    % also including brackets (inside html help files only!)
    quadStart=regexp(refText,'\[','once');
    if (~isempty(quadStart))
        quadEnd=regexp(refText,'\]','once');
        refText=[ refText(1:quadStart-1) refText(quadEnd+1:end)];
    end
       
    html = [html '<p>' refText '</p>' newline]; %#ok<AGROW>
    
    
end


% header and footer of the html reference file

topRefPage=sprintf(['<!DOCTYPE HTML>'...
    '<html itemscope="" itemtype="http://www.mathworks.com/help/schema/MathWorksDocPage" xmlns="http://www.w3.org/1999/xhtml">' ...
    '<script src="includesFS/headJSweb.js" type="text/javascript"></script>'...
    '<!--headJS loads all required javascripts--> <script type="text/javascript">'...
    'document.write(headJS); </script> <!--Insert title of the page--> <title>FSDA main page</title>'...
    '<!--Beginning of body--> <body id="responsive_offcanvas"> <div id="doc_header_spacer" class="header">7</div>'...
    '<!--Include serch engine--> <script type="text/javascript"> document.write(engine);'...
    ' </script> <!--Include left bar menu--> <div class="row-offcanvas row-offcanvas-left">'...
    '<script type="text/javascript"> 	document.write(lbar); 	</script>'...
    '<!--Include divs before FSDA text--> 	<script type="text/javascript">'...
    'document.write(divsbefore); </script> 	<!--BEGINNING OF FSDA TEXT--> 	<!-- --> <!-- --> 	<!-- -->'...
    '<h1 class="r2016a" itemprop="title content">Bibliography</h1>']);


bottomRefPage=sprintf(['	<!-- --> <!-- --> <!-- --> 	<!--END OF FSDA TEXT-->' ...
    '<!--Include divs after FSDA text--> 	<script language="javaScript" type="text/javascript">' ...
    'document.write(divsafter); </script> 	<!--Include fixed text at the bottom of the page-->' ...
    '<script language="javaScript" type="text/javascript"> document.write(barra);' ...
    '</script> 	<!--END.CLASS body_trail_container--> 	<!--close row-offcanvas--></div>' ...
    '<!--close class="row-offcanvas row-offcanvas-left" --> </body> </html> ']);

% build the reference file and write it to disk
RefPage=[topRefPage html bottomRefPage];

fileID = fopen('bibliography.html','w');
fprintf(fileID, RefPage);
fclose(fileID);


end

function References=refConv1(OUT, FilesIncluded)
% Creates a table of references

ri=size(OUT,1);
References=cell(2000,3);
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