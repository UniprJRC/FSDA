function [FSDAroot, cleanup] = changeDirToRootWithCleanup
% This helper function changes to the root folder of the toolbox and
% creates a cleanup object so that after execution MATLAB reverts to the
% original folder it was in

% Find my folder location - assumed to be TOOLBOX/utilities_help/build by
% getting my mfilename and getting the path from it
thisPath = fileparts(mfilename('fullpath'));
% Extract the root by getting the path twice (equivalent to ../..)
FSDAroot = string(fileparts(fileparts(thisPath)));

cdir = pwd;
cleanup = onCleanup(@() cd(cdir));
cd(FSDAroot)

end