function chkoptions(options,UserOptions)
%chkoptions checks whether all supplied options exist
%
% Required input arguments:
%
% options           : Input structure. Structure
% UserOptions       : cell array of strings containing user options
%
% Output:
%
%  The program checks if each string inside UserOptions is present in
%  structure option. If this condition is not fulfilled the execution
%  terminates and an error message is produced. This function is called
%  inside every routine of FSDA to check whether input options exist
%
% See also  chkinputR.m
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
%% example_producing_error
    %To examplify the behaviour of chkoptions, we call function FSR with a
    %non existing optional parameter ('namez').
    n=200;
    p=3;
    state1=123498;
    randn('state', state1);
    X=randn(n,p);
    y=randn(n,1);
    kk=33;
    nameX={'age', 'salary', 'position'};
    namey='salary';
    namez='error';
    [out]=FSR(y,X,'nameX',nameX,'namey',namey,'namez',namez);
%}

%% Beginning of code
inpchk=isfield(options,UserOptions);
WrongOptions=UserOptions(inpchk==0);
if ~isempty(WrongOptions)
    disp('-------------------------');
    for i=1:length(WrongOptions)
        disp(['Non existent user option found-> ' cellstr(WrongOptions(i))])
    end
    error('FSDA:chkoptions:WrongInputOpt','In total %d non-existent user options found.', length(WrongOptions));
end

end
