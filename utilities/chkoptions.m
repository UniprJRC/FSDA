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
% Copyright 2008-2011.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti 
%            and Vytis Kopustinskas (2009-2010)
%
% Last modified 15-Nov-2011

    inpchk=isfield(options,UserOptions);
    WrongOptions=UserOptions(inpchk==0);
    if ~isempty(WrongOptions)
        disp('-------------------------');
        disp(['Non existent user option found-> ' cellstr(WrongOptions)])
        error('Error: in total %d non-existent user options found.', length(WrongOptions));
    end

end
