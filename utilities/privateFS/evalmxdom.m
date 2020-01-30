function [dom, laste] = evalmxdom(f,dom,cellBoundaries,imagePrefix,imageDir,outputDir,options)
%EVALMXDOM   Evaluate cellscript Document Object Model, generating inline images.
%   dom = evaldom(dom,imagePrefix,imageDir,options)
%   imagePrefix is the prefix that will be used for the image files.

% Copyright 1984-2014 The MathWorks, Inc.
%
% Copyright 2008-2019.
%
%$LastChangedDate::                      $: Date of the last commit

%#ok<*AGROW> (nodeList needs to grow)


%% Beginning of code 

% Provide a clean, white figure.
if options.useNewFigure
    figureName = 'Figure';
    myFigure = figure( ...
        'IntegerHandle','off','NumberTitle','off','Name',figureName);
else
    myFigure = [];
end

% Run the code.
file = which(f);
[data,text,laste] = instrumentAndRun(file,cellBoundaries,imageDir,imagePrefix,options);
createThumbnail(options,data)

% Make paths relative when the imageDir is the same as the outputDir.
for i = find(strncmp(outputDir,data.pictureList,numel(outputDir)))
    data.pictureList{i} = data.pictureList{i}(numel(outputDir)+2:end);
end

% Cleanup the provided figure, if it is still around.
close(myFigure(ishandle(myFigure)))

% Handle errors.
if ~isempty(laste)
    if options.catchError
        if options.displayError
            disp(formatError(laste,true))
            beep
        end
    else
        rethrow(laste)
    end
end

% Dump all the new data into the DOM.
populateDom(dom,data,text,laste)

end

function [data,text,laste] = instrumentAndRun(file,cellBoundaries,imageDir,imagePrefix,options)

% Add conditional breakpoints.
originalDbstatus = dbstatus(file,'-completenames');
resetDbstatusObj = onCleanup(@()restoreDbstatus(file, originalDbstatus));
    function restoreDbstatus(file, originalDbstatus)
        safeDbclear(file)
        dbstop(originalDbstatus,'-completenames')
    end
endCondition = setConditionalBreakpoints(file,cellBoundaries);

% Initialize publishSnapshot.
data = [];
data.baseImageName = fullfile(imageDir,imagePrefix);
data.options = options;
data.marker = 'WDEavRCxr';

resetSnapnowObj = onCleanup(@()snapnow('set',[]));
snapnow('set',data)

% For consistency, set the published output to 80 columns, regardless of 
% the current Command Window width.
originalEightyColumns = feature('EightyColumns');
resetEightyColumnsObj = onCleanup(...
    @()feature('EightyColumns',originalEightyColumns) );
feature('EightyColumns',1)

% Be sure we don't capture the code itself in the output if echo is on.
echoState = get(0,'echo');
resetEchoObj = onCleanup(@()set(0,'echo',echoState));
echo off

% Run the command.
[~,tempVar] = fileparts(tempname);
resetAppdataObj = onCleanup(@()rmappdata(0,tempVar));
setappdata(0,tempVar,options.codeToEvaluate)

ret = char(10);
laste = [];
evalstr = [ ...
    'feature(''hotlinks'',0);' ret ...
    'try' ret ...
    'evalin(''base'',getappdata(0,''' tempVar '''));' ret ...
    endCondition ';' ret ...
    'catch laste, end'];

warnState = warning('off','backtrace');
resetWarnStateObj = onCleanup(@()warning(warnState));
[text] = evalc(evalstr);
delete(resetWarnStateObj); % To restore warning state now.

% Get the latest info out of publishSnapshot.
data = snapnow('get',data);
if isempty(data)
    m = message('MATLAB:publish:RecursivePublish');
    error(m);
end

% Restore to the original state
delete(resetAppdataObj);
delete(resetSnapnowObj);
delete(resetDbstatusObj);
delete(resetEightyColumnsObj);
delete(resetEchoObj);

end

function endCondition = setConditionalBreakpoints(file, cellBoundaries)

% Ensure there are no breakpoints in the file, for findLandingLines.
safeDbclear(file)

% Build a Map of the line's conditions.  
% Ensure there is an "Inf" entry for beyond-the-end.
snapnowArgs = containers.Map({Inf}, {{}});
for iCell = 1:size(cellBoundaries, 1)
    iStartLine = findLandingLine(file, cellBoundaries(iCell, 1));
    iEndLine   = findLandingLine(file, cellBoundaries(iCell, 2) + 1);
    
    if iStartLine == iEndLine
        % Nothing; the cell has no executable line.
    else
        addArgsToMap(iStartLine, 'beginCell', iCell);
        addArgsToMap(  iEndLine,   'endCell', iCell);
    end
end

% Make a structure suitable for DBSTOP.
[dbStruct, endCondition] = dbStructFromMap;

% Set the conditional breakpoints; DBSTOP crashes if dbStruct is "empty".
if ~isempty(dbStruct.line)
    dbstop(dbStruct,'-completenames')
end

    function addArgsToMap(iLine, beginOrEnd, iCell)
        if isKey(snapnowArgs, iLine)
            oldArgs = snapnowArgs(iLine); 
        else
            oldArgs = {};
        end
            
        if isinf(iLine) && strcmp(beginOrEnd, 'endCell')
            snapnowArgs(iLine) = [{beginOrEnd iCell} oldArgs];
        else
            snapnowArgs(iLine) = [oldArgs {beginOrEnd iCell}];
        end
    end

    function [dbStruct, endCondition] = dbStructFromMap        
        % Pull out the end condition.
        endConditionArgs = snapnowArgs(Inf);
        if numel(endConditionArgs) > 1            
            endConditionArgs{end} = 0-endConditionArgs{end};
        end
        endCondition = convertArgsToString(endConditionArgs);
        remove(snapnowArgs, Inf);
        
        % Form breakpoint structure.
        dbStruct.name = file;
        dbStruct.file = file;
        dbStruct.cond = [];
        dbStruct.identifier = {};
        dbStruct.line = cell2mat(keys(snapnowArgs));
        dbStruct.anonymous = zeros(size(dbStruct.line));
        dbStruct.expression = cellfun( ...
            @convertArgsToString, values(snapnowArgs), ...
            'UniformOutput', false);
    end

    function str = convertArgsToString(argList)
        if isempty(argList)
            str = '';
        else
            str = sprintf(', ''%s'', %d', argList{:});
            str = sprintf('snapnow(%s);', str(3:end));
        end
    end
end

function populateDom(dom,data,text,laste)

% Create nodeList.
nodeList = [];

% Populate images into nodeList.
for i = 1:numel(data.pictureList)
    imgNode = dom.createElement('img');
    imgNode.setAttribute('src',data.pictureList{i});
    nodeList(end+1).node = imgNode;
    nodeList(end).cell = data.placeList(2*i-1);
    nodeList(end).count = data.placeList(2*i);
end

% Apply any "backspace" characters.
while true
    codeOutputSize = numel(text);
    text = regexprep(text,('([^\b])\b'),'');
    if codeOutputSize == numel(text)
        break
    end
end

% Populate M-code output into nodeList.
[cellTextNumbers,cellTextTexts,cellTextCounts] = textparse(text,data.marker);
for i = 1:numel(cellTextNumbers)
    cell = cellTextNumbers(i);
    text = cellTextTexts{i};
    count = cellTextCounts(i);
    if ~isempty(text)
        mcodeOutputNode = dom.createElement('mcodeoutput');
        mcodeOutputTextNode = dom.createTextNode(text);
        mcodeOutputNode.appendChild(mcodeOutputTextNode);
        mcodeOutputNode.setAttribute('class', 'codeoutput');
        nodeList(end+1).node = mcodeOutputNode;
        nodeList(end).cell = cell;
        nodeList(end).count = count;
    end
end

% Populate error into nodeList.
if ~isempty(laste)
    mcodeOutputNode = dom.createElement('mcodeoutput');
    mcodeOutputTextNode = dom.createTextNode(formatError(laste,false));
    mcodeOutputNode.appendChild(mcodeOutputTextNode);
    mcodeOutputNode.setAttribute('class', 'codeoutput error');
    nodeList(end+1).node = mcodeOutputNode;
    % TODO Flag as error node so it will be red.
    nodeList(end).cell = data.lastGo;
    nodeList(end).count = data.counter();
end

% Put nodeList into cells.
if ~isempty(nodeList)
    [~,ind] = sort([nodeList.count]);
    nodeList = nodeList(ind);
    cellOutputTargetList = dom.getElementsByTagName('cellOutputTarget');
    for i = 1:numel(nodeList)
        for iCellOutputTargetList = 1:cellOutputTargetList.getLength
            cellOutputTarget = cellOutputTargetList.item(iCellOutputTargetList-1);
            if str2double(char(cellOutputTarget.getTextContent)) == nodeList(i).cell
                cellOutputTarget.getParentNode.appendChild(nodeList(i).node);
                break
            end
        end
    end
end

end

function createThumbnail(options,data)
% Create thumbnail.

if ~options.createThumbnail
    return
end

for iPictureList = numel(data.pictureList):-1:1
    picture = data.pictureList{iPictureList};
    [~,~,ext] = fileparts(picture);
    if any(strcmp(ext,{'.png','.jpg','.tif','.gif','.bmp'}))
        [myFrame.cdata,myFrame.colormap] = robustImread(picture);
        internal.matlab.publish.writeImage( ...
            [data.baseImageName '.png'],'png',myFrame, ...
            options.maxThumbnailHeight,options.maxThumbnailWidth);
        break
    end
end

end

function [X,map] = robustImread(picture)
% In case of filesystem lag, where the file we just wrote out is not yet
% available for reading, try a few times.
nTries = 10;
for iTries = 1:nTries
    try
        [X,map] = imread(picture);
        break
    catch e
        warning(pm('CannotReadForThumbnail',picture,nTries))
        missingErr = 'MATLAB:imagesci:imread:fileDoesNotExist';
        if ~strcmp(e.identifier,missingErr) || (iTries == nTries)
            rethrow(e)
        end
        pause(1)
    end
end
end

function landingLine = findLandingLine(file,targetLine)
% Probe to see where the breakpoint wants to land.
% Precondition: this file has no existing breakpoints.
try
    dbstop(file,num2str(targetLine),'-completenames')
catch anError %#ok<NASGU> Only assigned to a variable to prevent changing LASTERROR.
    % This errors if there aren't any executible lines left or if there is
    % any sort of parse error.
end
tempDbstatus = dbstatus(file,'-completenames');
safeDbclear(file)
if isempty(tempDbstatus) || isempty(tempDbstatus.line)
    landingLine = Inf;
else
    landingLine = tempDbstatus.line;
end

end

function safeDbclear(file)
try
    dbclear('in',file,'-completenames')
catch anError %#ok<NASGU> Only assigned to a variable to prevent changing LASTERROR.
    % This errors when the file contains a parse error, which corresponds
    % to a variety of error IDs.
end
end

function [cellTextNumbers,cellTextTexts,cellTextCounts] = textparse( ...
    text,marker, ...
    cellTextNumbers,cellTextTexts,cellTextCounts, ...
    enclosingCell,enclosingCount)

% Initialize varables to hold return data.
if nargin < 3
    cellTextNumbers = [];
    cellTextTexts = {};
    cellTextCounts = [];
end

% Loop and recurse until all has been parsed.
while ~isempty(text)

    % Find the first cell start marker in TEXT.
    [tMatches,mStartStart,mStartEnd] = ...
        regexp(text,[marker '(\d+)A(\d+)X'],'tokens','start','end','once');

    % If there are no start markers left, push the text on the list.
    if isempty(tMatches)
        % Trim any remaining markers.
        text = regexprep(text,[marker '\d+[AZ]\d+X'],'');
        % Add to the document.
        if exist('enclosingCell','var')
            cellTextNumbers(end+1) = enclosingCell;
            cellTextTexts{end+1} = text;
            cellTextCounts(end+1) = enclosingCount;
        else
            % This shouldn't happen, so display for debugging purposes.
            disp(text)
        end
        break
    end
    
    % Extract the cell and count for this start marker.
    startMarker = text(mStartStart:mStartEnd);
    iCell = str2double(tMatches(1));
    iCount = str2double(tMatches(2));

    % Locate the matching close maker for this cell.
    closeMarker = sprintf('(?<=%s.*)%s%iZ(\\d+)X',startMarker,marker,iCell);
    [match,mEndStart,mEndEnd] = regexp( ...
        text,closeMarker, ...
        'match','start','end','once');
    if isempty(match)
        % We should only get here if the code we're publishing errors.
        mEndStart = numel(text)+1;
        mEndEnd = numel(text);
    end
    
    % Recurse on the contents of this match and excise it from TEXT.
    [cellTextNumbers,cellTextTexts,cellTextCounts] = textparse( ...
        text(mStartEnd+1:mEndStart-1),marker, ...
        cellTextNumbers,cellTextTexts,cellTextCounts, ...
        iCell,iCount);
    text(mStartStart:mEndEnd) = [];

end

end

function m = formatError(laste,hotlinks)
    % Get the stack trace message from the error.
    origHotlinks = feature('hotlinks',hotlinks);
    m = getReport(laste);
    feature('hotlinks',origHotlinks)
    iRet = find(m==sprintf('\n'));
    
    % Find where the publishing stack begins in the error message.
    instrumentAndRunPos = strfind(m,'evalmxdom>instrumentAndRun');
    if isempty(instrumentAndRunPos)
        % This shouldn't happen, but in the case of misbehaved errors, just
        % show the whole thing.
        stackStart = [];
    else
        stackStart = max(iRet(iRet < instrumentAndRunPos(end)));
    end
    
    % Find the end of the publishing stack in the error message.
    causePos = strfind(m,getCauseString);
    if isempty(causePos)
        stackEnd = numel(m);
    else
        stackEnd = max(iRet(iRet < causePos(end)))-1;
    end
    
    % Trim out the publishing stack from the error.
    m(stackStart:stackEnd) = []; 
    
    % Trim trailing newline.
    if ~isempty(causePos)
        m(end) = [];
    end
end

function s = getCauseString
    M1 = MException('foo:bar','');
    M2 = MException('foo2:bar2','');
    M2 = M2.addCause(M1);
    s = strtrim(getReport(M2));
end

function m = pm(id,varargin)
m = message(['MATLAB:publish:' id],varargin{:});
end
