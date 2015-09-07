function docsearchFS(varargin)
%docsearchFS opens inside Help browser html documentation file.
%
%<a href="matlab: docsearchFS('docsearchFS')">Link to the help function</a>
%
%   docsearchFS opens the Help browser and displays the documentation home
%   page of FS. 
%
%   docsearchFS('namehtmlhelpfile') searches FSDA documentation for pages
%   with words that match the specified expression (version of matlab
%   <=2012a) or opens FSDA documentation associated to namehtmlhelpfile.
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
    % Open main documentation page of FSDA
    docsearchFS
%}
%
%
%{
    % Open html documentation page of FSDA function named LXS
    docsearchFS('lxs')
%}

%% Beginning of code
if nargin > 1
    namehtmlhelpfile = deblank(sprintf('%s ', varargin{:}));
elseif nargin == 1
    namehtmlhelpfile = varargin{1};
else
    namehtmlhelpfile = '';
end

% Find version of MATLAB which is installed 
a=ver('matlab');

if str2double(a.Version)>7.14
    % If installed version of MATLAB is newer than 2012a 
    % function helpbrowser is called
    
    % Find path of docsearchFS.m
    mname=mfilename('fullpath');
    
    if isempty(namehtmlhelpfile)
        % if namehtmlhelpfile is not specified main page of FSDA
        % documentation is opened
        rname = ['' mname(1:end-11) 'helpfiles' filesep 'FSDA' filesep 'fsda_product_page.html' ''];
        web(rname,'-new')
    else
        
        % rnames contains full path of associated html file
        rname = ['' mname(1:end-11) 'helpfiles' filesep 'FSDA' filesep namehtmlhelpfile '.html' ''];
        % Open html help file
        web(rname,'-new')
    end
else
    % If installed version of MATLAB is 2012a or older, matlab function docsearch
    % is called
    if isempty(namehtmlhelpfile)
        docsearch('FSDA')
    else
        docsearch(namehtmlhelpfile)
    end
end