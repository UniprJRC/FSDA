function options = supplyDefaultOptions(options)
% Supply default options for any that are missing.
%
% Copyright 2008-2023.
%
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code 

if ~isfield(options,'format')
    options.format = 'html';
end
format = options.format;
if ~isfield(options,'stylesheet') || isempty(options.stylesheet)
    switch format
        case 'html'
            codepadDir = fileparts(which(mfilename));
            styleSheet = fullfile(codepadDir,'mxdom2simplehtml.xsl');
            options.stylesheet = styleSheet;
        case 'latex'
            codepadDir = fileparts(which(mfilename));
            styleSheet = fullfile(codepadDir,'private','mxdom2latex.xsl');
            options.stylesheet = styleSheet;
        case {'docbook','pdf'}
            codepadDir = fileparts(which(mfilename));
            styleSheet = fullfile(codepadDir,'private','mxdom2docbook.xsl');
            options.stylesheet = styleSheet;
        otherwise
            options.stylesheet = '';
    end
end
if ~isfield(options,'figureSnapMethod')
    options.figureSnapMethod = 'entireGUIWindow';
end
if ~isfield(options,'imageFormat') || isempty(options.imageFormat)
    options.imageFormat = '';
elseif strcmp(options.imageFormat,'jpg')
    options.imageFormat = 'jpeg';
elseif strcmp(options.imageFormat,'tif')
    options.imageFormat = 'tiff';
elseif strcmp(options.imageFormat,'gif')
        error('FSDA:supplyDefaultOptions:NoGIFs','no gif images');
end
if ~isfield(options,'useNewFigure')
    options.useNewFigure = true;
end
if ~isfield(options,'maxHeight')
    options.maxHeight = [];
end
if ~isfield(options,'maxWidth')
    options.maxWidth = [];
end
if ~isfield(options,'maxThumbnailHeight')
    options.maxThumbnailHeight = 64;
end
if ~isfield(options,'maxThumbnailWidth')
    options.maxThumbnailWidth = 85;
end
if ~isfield(options,'showCode')
    options.showCode = true;
end
if ~isfield(options,'evalCode')
    options.evalCode = true;
end
if ~isfield(options,'stopOnError')
    options.stopOnError = true;
end
if ~isfield(options,'catchError')
    options.catchError = true;
end
if ~isfield(options,'displayError')
    options.displayError = true;
end
if ~isfield(options,'createThumbnail')
    options.createThumbnail = true;
end
if ~isfield(options,'maxOutputLines')
    options.maxOutputLines = Inf;
end
if ~isfield(options,'codeToEvaluate')
    options.codeToEvaluate = '';
end
if ~isfield(options,'font')
    options.font = '';
end
if ~isfield(options,'titleFont')
    options.titleFont = options.font;
end
if ~isfield(options,'bodyFont')
    options.bodyFont = options.font;
end
if ~isfield(options,'monospaceFont')
    options.monospaceFont = options.font;
end
end
