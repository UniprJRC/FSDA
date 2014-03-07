function docsearchFS(namehtmlhelpfile)
% Find version of MATLAB which is installed
a=ver('matlab');

if str2double(a.Version)>7.14
    % If installed version of MATLAB is newer than 2012a command
    % function helpbrowser is called
    
    % Find path of docsearchFS.m
    mname=mfilename('fullpath');
    % Find full path of associated html file
    rname = ['' mname(1:end-11) 'helpfiles\FSDA\' namehtmlhelpfile '.html' ''];
    % Open html help file
    web(rname,'-helpbrowser')
else
    % If installed version of MATLAB is 2012a or older function docsearch
    % is called
    docsearch(namehtmlhelpfile)
end