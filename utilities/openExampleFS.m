function openExampleFS(filename)
%openExampleFS Open an example for modification and execution.
% Copyright 2008-2021.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

    examp=which(filename);
    
    examp1=strrep(examp,'\','\\');
    opentoline(examp1,5)
end
