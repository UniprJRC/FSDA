function chkoptions(options,UserOptions)
%chkoptions checks whether all supplied options exist
%
% Required input arguments:
%
% options           : a structure
% UserOptions       : cell array of strings
%
% Output:
%
%  The program checks if each string inside UserOptions is present in
%  structure option. If this condition is not fulfilled the execution
%  terminates and an error message is produced
%
% See also  chkinputR.m
%
% Copyright 2008-2014.
% Written by FSDA team
%
%
% Last modified 08-Dec-2013

%% Beginning of code
inpchk=isfield(options,UserOptions);
WrongOptions=UserOptions(inpchk==0);
if ~isempty(WrongOptions)
    disp('-------------------------');
    disp(['Non existent user option found-> ' cellstr(WrongOptions)])
    error('Error: in total %d non-existent user options found.', length(WrongOptions));
end

end
