function [dom,cellBoundaries] = m2mxdom(originalCode)
%M2MXDOM  Converts codepad-style m-code into a Document Object Model.
%   M2MXDOM(TXT) parses the char array TXT and returns the contents as
%   a cell script DOM.
 
% Copyright 1984-2012 The MathWorks, Inc.
 
%#ok<*AGROW> The length of chunkList and paragraphText are not known in advance.
 
% Normalize line endings to Unix-style.
code = regexprep(originalCode,'\r\n?','\n');
newLine = sprintf('\n');
 
% Trim trailing whitespace.
code = regexprep(code,'[ \t]+(\n|$)','\n');
 
% Exactly one newline at the end of the file.
code = regexprep(code,'(.)\n*$','$1\n');
 
% Find the cells.
cellLocations = com.mathworks.widgets.text.mcode.cell.CellUtils.getCellLocations(code);
if isempty(cellLocations) || (cellLocations(1) ~= 1)
    % Add the implicit outer cell.
    cellLocations = [1; sum(code==newLine); cellLocations];
end
cellStarts = cellLocations(1:2:end);
cellEnds = cellLocations(2:2:end);
cellBreaks = unique([cellStarts' cellEnds'+1]);
cellBoundaries = [cellStarts cellEnds];
 
dom = createDom;
ignoreMarker = 'XXX_IGNORE_THIS_LINE_XXX';
code = ignoreStuff(code,dom,ignoreMarker);

% Use the new parser. 
chunkList = findChunks(code,cellBreaks,cellEnds,newLine);
addChunksToDom(chunkList,dom,ignoreMarker);

finalizeDom(dom,originalCode)
end
 
%===============================================================================
function dom = createDom
% Now create the new DOM
dom = com.mathworks.xml.XMLUtils.createDocument('mscript');
dom.getDocumentElement.setAttribute('xmlns:mwsh', ...
    'http://www.mathworks.com/namespace/mcode/v1/syntaxhighlight.dtd')
 
% Add version.
newNode = dom.createElement('version');
matlabVersion = ver('MATLAB');
newTextNode = dom.createTextNode(matlabVersion(1).Version);
newNode.appendChild(newTextNode);
dom.getFirstChild.appendChild(newNode);

% Add release.
newNode = dom.createElement('release');
release = version('-release');
newTextNode = dom.createTextNode(release);
newNode.appendChild(newTextNode);
dom.getFirstChild.appendChild(newNode);
 
% Add date.
newNode = dom.createElement('date');
newTextNode = dom.createTextNode(datestr(now,29));
newNode.appendChild(newTextNode);
dom.getFirstChild.appendChild(newNode);
end
 
%===============================================================================
function code = ignoreStuff(code,dom,ignoreMarker)
% Exclude some lines from the published output.
    function excludeLine(pattern,iKeep,nodeName)
    paddedPattern = ['(?<=\n|^)' pattern '(?:\n)'];
    paddedReplace = ['% ' ignoreMarker '\n'];
    match = regexp(code,paddedPattern,'tokens','once');
    if ~isempty(match)
        code = regexprep(code,paddedPattern,paddedReplace);
        if (nargin > 1)
            node = dom.createElement(nodeName);
            dom.getDocumentElement.appendChild(node);
            node.appendChild(dom.createTextNode(match{iKeep}));
        end
    end
    end
excludeLine( ...
    '(%[ \t]*)(Copyright [\d-]+ [^\n]*The MathWorks, Inc\.)', ...
    2, ...
    'copyright')
excludeLine( ...
    '(%[ \t]*)(\$(?:Revision|Date):.*?\$[ \t]*)+', ...
    2, ...
    'revision')
excludeLine('([ \t]*)displayEndOfDemoMessage\(mfilename\)')
end
 
%===============================================================================
% New parser
function chunkList = findChunks(code,cellBreaks,cellEnds,newLine)
% Initialize variables for loop.
chunkList = [];
returns = find(code==newLine);
% Iterate over each cell, possibly busting it up into smaller chunks.
for iCodeBreaks = 1:numel(cellBreaks)-1
    % Find all cell breaks..
    cellStartLine = cellBreaks(iCodeBreaks);
    cellEndLine = cellBreaks(iCodeBreaks+1)-1;
    if cellStartLine == 1
        chunkStartIndex = 1;
    else
        chunkStartIndex = returns(cellStartLine-1)+1;
    end
    chunkEndIndex = returns(cellEndLine);
    chunk = code(chunkStartIndex:chunkEndIndex);
    nextChunk = newChunk();
    nextChunk.originalCode = chunk;
    chunkList = [chunkList nextChunk];
    chunkList(end).outputTargets = find(cellEndLine == cellEnds);
end
end
 
%===============================================================================
% New parser
function addChunksToDom(chunkList,dom,ignoreMarker)
 
% Loop over each chunk in the code.
codeCount = 0;
cellNumber = 0;
for n = 1:length(chunkList)
    frag = com.mathworks.publishparser.PublishParser.getDomFragment( ...
        dom, ...
        regexprep(chunkList(n).originalCode,['% ' ignoreMarker],''), ...
        cellNumber, ...
        codeCount);
    if frag.hasChildNodes
        for cellOutputTarget = chunkList(n).outputTargets'
            cellOutputTargetNode = dom.createElement('cellOutputTarget');
            newTextNode = dom.createTextNode(num2str(cellOutputTarget));
            cellOutputTargetNode.appendChild(newTextNode);
            frag.getLastChild().appendChild(cellOutputTargetNode);
        end
    end
    dom.getDocumentElement.appendChild(frag);
    codeCount = dom.getElementsByTagName('mcode-count').getLength();
    cellNumber = dom.getElementsByTagName('count').getLength();
end
end
 
%===============================================================================
function finalizeDom(dom,originalCode)
 
% Tag the first cell if it is an "Overview".
cellList = dom.getFirstChild.getElementsByTagName('cell');
if (cellList.getLength > 1) && ...
        (cellList.item(0).getElementsByTagName('mcode').getLength == 0)
    cellList.item(0).setAttribute('style','overview')
    % A title in the "Overview" cell is the document title.
    firstStepTitle = cellList.item(0).getElementsByTagName('steptitle');
    if (firstStepTitle.getLength == 1)
        firstStepTitle.item(0).setAttribute('style','document')
    end
end
 
% Potentially tag the first steptitle as the document title.
if (dom.getElementsByTagName('steptitle').getLength == 1)
    firstStepTitle = cellList.item(0).getElementsByTagName('steptitle');
    if (firstStepTitle.getLength == 1)
        firstStepTitle.item(0).setAttribute('style','document')
    end
end
 
% Save the virgin code in a node.
originalCodeNode = dom.createElement('originalCode');
dom.getDocumentElement.appendChild(originalCodeNode);
originalCodeNode.appendChild(dom.createTextNode(originalCode));
 
end
 
%===============================================================================
function chunk = newChunk()
chunk = [];
chunk.title = '';
chunk.text = {};
chunk.code = '';
chunk.outputTargets = [];
end
