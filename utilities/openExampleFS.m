function openExampleFS(filename)
%openExampleFS Open an example for modification and execution.
% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

    examp=which(filename);
    
    examp1=strrep(examp,'\','\\');
    opentoline(examp1,5)
end
