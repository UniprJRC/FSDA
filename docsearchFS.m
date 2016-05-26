function docsearchFS(varargin)
%docsearchFS opens inside Help browser html documentation file.
%
%<a href="matlab: docsearchFS('docsearchFS')">Link to the help function</a>
%
%   docsearchFS opens the Help browser and displays the documentation home
%   page of FS.
%
%   docsearchFS('namehtmlhelpfile') searches FSDA documentation page
%   of namehtmlhelpfile
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('docsearchFS')">Link to the help page for this function</a>
% Last modified 06-Feb-2015

% Examples:

%
%{
    % Open main documentation page of FSDA (file index.html)
    docsearchFS
%}
%
%
%{
    % Open html documentation page of FSDA function named LXS
    docsearchFS('LXS')
%}

%% Beginning of code
if nargin <1
    namehtmlhelpfile = deblank(sprintf('%s ', varargin{:}));
elseif nargin == 1
    namehtmlhelpfile = varargin{1};
else
    namehtmlhelpfile = '';
end

a=ver('matlab');

if str2double(a.Version)>7.14
    
    if isempty(namehtmlhelpfile)
        web([docroot '/FSDA/index.html'])
    else
        
        [~,~,ext]=fileparts(namehtmlhelpfile);
        
        if isempty(ext)
            namehtmlhelpfile=[namehtmlhelpfile '.html'];
        elseif strcmp(ext,'html')==0
            error('FSDA:docsearchFS','Wrong file extension')
        end
        web([docroot '/FSDA/' namehtmlhelpfile])
    end
    
else
    % If installed version of MATLAB is 2012a or older, matlab function web
    % is called
    FileWithFullPath=which('docsearchFS.m');
    [pathFSDAstr]=fileparts(FileWithFullPath);
    fsep=filesep;
    
    if isempty(namehtmlhelpfile)
        outputOFHtmlHelpFile=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA\index.html'];
        web(outputOFHtmlHelpFile);
        
    else
        [~,~,ext]=fileparts(namehtmlhelpfile);
        if isempty(ext)
            outputOFHtmlHelpFile=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA' filesep namehtmlhelpfile '.html'];
        else
            outputOFHtmlHelpFile=[pathFSDAstr fsep 'helpfiles' fsep 'FSDA' filesep namehtmlhelpfile];
        end
        web(outputOFHtmlHelpFile);
    end
    
    
end

end