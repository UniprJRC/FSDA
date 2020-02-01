function setuphelp(FSDApath)
%Load appropriate HTML help documentation files depending on the installed version of MATLAB
%
%  Required input arguments:
%
%    FSDA:      A string containing the path which contains the root folder of FSDA toolbox
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%
%{
%      If FSDA has been installed in D:\tmp13\FSDA 
%      in order to include the HTML help documentation system of FSDA use the following sintax
 
       setuphelp('D:\tmp13\FSDA')
%}

% The purpose of this file is to extract the appropriate help
% find the version of MATLAB currently  in use


%% Beginning of code
% identify FSDA root
%FSDAroot=which('logo.png');
% if ~isequal(FSDAroot(1:end-9),pwd)
%     disp('please go the root folder where FSDA is installed before running this file')
%     return
% else
%     str = ['Is ' FSDAroot ' the path where you installed FSDA? Y/N [Y]: '];
%     if isempty(str)
%         str = 'Y';
%     end
%     if strcmp(str','Y');
%     end
% end

a=version;

if str2double(a(1)) >= 8
    source   = [FSDApath filesep 'helpfiles' filesep 'FSDAR8'];
    toremove = [FSDApath filesep 'helpfiles' filesep 'FSDAR7'];
else
    source   = [FSDApath filesep 'helpfiles' filesep 'FSDAR7'];
    toremove = [FSDApath filesep 'helpfiles' filesep 'FSDAR8'];
end

% destination folder
destination= [FSDApath '\helpfiles\FSDA'];

% change name from FSDARX to FSDA
if exist(source,'file')
    movefile(source,destination)
else
    disp(['path ' source ' does not exist'])
    return
end

% Maybe this instruction in unnecessary
rmdir(toremove,'s')

end
