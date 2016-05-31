function [ss]=waitforbuttonpressFS
%waitforbuttonpressFS produces a warning message when user closes the figure he is interacting with 
%
% Copyright 2008-2016.
% Written by FSDA team
%
% Last modified ven  6 feb 2015 12:22:42
%
%

%% Beginning of code
try
    ss=1;
    ss=waitforbuttonpress;
catch ME
    if strcmp(ME.message,'waitforbuttonpress exit because all figures have been deleted')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        disp('The figure you are interacting with has been deleted: PROCESS ENDED');
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        ss=1;
    end
end
end