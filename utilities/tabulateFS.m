function tbl = tabulateFS(x)
%Create frequency table of unique values of x, excluding possible 0 counts
%
%<a href="matlab: docsearchFS('tabulateFS')">Link to the help function</a>
%
%   tbl = tabulateFS(x) takes a vector x and returns a matrix, tbl. The
%   first column of table contains the unique values of x.  The second is
%   the number of instances of each value.  The last column contains the
%   percentage of each value.  This function differs from MATLAB function
%   tabulate because it excludes 0 counts.
%   Remark: tabulateFS with no output arguments returns a formatted table
%   in the command window.
%
%  Required input arguments:
%
%            x: vector for which frequency table has to be calculated.
%               vector of numeric data or  categorical variable, character
%               array, or cell array of strings of length n.
%               Data Types - double | single| categorical variable | character
%               array, | cell array of strings
%
%  Optional input arguments:
%
%  Output:
%
%    tbl :  frequency table of data in vector x.
%           Matrix of size unique(x)-by-3.
%           Information in tbl is arranged as follows:
%               1st column -- The unique values of x;
%               2nd column -- The number of instances of each value;
%               3rd column -- The percentage of each value.
%           If x is a categorical variable, character array, or cell array of
%           strings, tbl is a cell array.
%
%  Remark: tabulateFS with no output arguments returns a formatted table
%          in the command window.
%
% See also: tabulate
%
% References:
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('tabulateFS')">Link to the help page for this function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{
    % Tabulate Fisher iris data.
    load fisheriris;
    tbl=tabulateFS(species);
    disp(tbl)
%}

%{
    %% Explore the difference between tabulate and tabulateFS.
    % Run this code to see the output shown in the help file
    rng(100)
    x=randi([1 10],100,1);
    x(100)=30;
    % Output of tabulate
    disp('Output of MATLAB function tabulate')
    disp(tabulate(x));
    % Output of tabulateFS
    disp('Output of FSDA function tabulateFS')
    disp(tabulateFS(x));
%}

%% Beginning of code

tbl = tabulate(x); % frequency count

% exclude 0 counts
if isnumeric(tbl)
    tbl = tbl(tbl(:,2)~=0, :);
    if nargout == 0
        values=tbl(:,1);
        counts=tbl(:,2);
        percents=tbl(:,3);
        fprintf(1,'  %5d    %5d    %6.2f%%\n',[values counts percents]');
    end
else
    if nargout == 0
        disp(tbl)
    end
end