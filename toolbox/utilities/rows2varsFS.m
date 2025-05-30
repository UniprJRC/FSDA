function T2=rows2varsFS(T1, varargin)
%rows2varsFS calls rows2vars but the the names of the cols of the original table become rownames
%
%<a href="matlab: docsearchFS('rows2varsFS')">Link to the help function</a>
%
% MATLAB function rows2vars transforms T1 of size n-by-p into T2 of size
% p-by-(n+1). The first column of T2 is called OriginalVariableNames and
% contains the names of the table variables from T1 and the RowNames of T2
% are empty. FSDA function rows2varsFS transforms T1 of size n-by-p into T2 of
% size p-by-n. The names of the variables of T1 become the names of the
% rows of T2. For example, if T1 has two rows and five columns,
% after the call T2=rows2varsFS(T1) has five rows and 5 columns.
% The RowNames of T2 are the original variable names of T1 and viceversa.
%
%
% Required input arguments:
%
%    T1:         Origianal table. table or timetable. 
%                Input table, specified as a table or timetable of size n-by-p.
%
% Optional input arguments:
%
% VariableNamesSource. char, string, numeric or lofical.
%                 Variable in T1 that contains variable names
%                 See help of rows2vars for additional details
%               Example - 'VariableNamesSource',4
%               Data Types - character vector | string scalar | positive integer | logical vector
%
%
%  DataVariables   : Selected variables from T1. String array or character vector.
%               Indication of the variables that have to be reoriented.
%               See help of rows2vars for additional details
%               Example - 'bonflev',0.99
%               Data Types - string array | character vector | cell array of character vectors | pattern scalar | positive integer | vector of positive integers | logical vector
%
% VariableNamingRule :  Rule for naming variables in T2.  
% |             modify' (default) | 'preserve'.
%               See help of rows2vars for additional details
%               Example - 'VariableNamingRule','preserve'
%               Data Types - string  | char
%
%
% Output:
%
%         T2 :   table with p rows and n columns. Where p is the number of columns of T1
%                or length(DataVariables). n is the number of rows of the input table.
%                The RowNames of T2 are the VariableNames of T1 the VariableNames specified in DataVariables 
%
% See also: rows2vars
%
% References:
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('rows2varsFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% rows2varsFS with all default options.
    load patients
    T1 = table(LastName,Gender,Age,Height,Weight);
    disp('Original table (first 5 rows)')
    disp(head(T1,5))
    disp('Output of rows2vars (first 5 columns)')
    T2=rows2vars(T1)
    disp(head(T2(:,1:5)))
    disp('Output of rows2varsFS (first 5 columns)')
    T2fs=rows2varsFS(T1);
    disp(head(T2fs(:,1:5)))
%}

%{
    %% Example with RowNames in the original table.
    load patients
    T1 = table(Gender,Age,Height,Weight,'RowNames',LastName);
    disp('Original table (first 5 rows)')
    head(T1,5)
    T2=rows2vars(T1);
    disp('Output of rows2vars (first 5 columns)')
    disp(head(T2(:,1:5)))
    disp('Output of rows2varsFS (first 5 columns)')
    T2fs=rows2varsFS(T1);
    disp(head(T2fs(:,1:5)))
%}

%{
    %% Example of use of option "DataVariables".
    load patients
    T1 = table(Gender,Age,Height,Weight,'RowNames',LastName);
    disp('Original table (first 5 rows)')
    head(T1,5)
    disp('Output of rows2vars (first 5 columns)')
    T2=rows2vars(T1,'DataVariables',["Gender" "Age"]);
    disp(head(T2(:,1:5)))
    disp('Output of rows2varsFS ')
    T2fs=rows2varsFS(T1,'DataVariables',["Gender" "Age"]);
    disp(head(T2fs(:,1:5)))
%}

%{
    %% Example of use of option "VariableNamesSource".
    T1 = readtable('patients.xls');
    disp('Original table (first 5 rows)')
    head(T1,5)
    T2=rows2vars(T1,'VariableNamesSource','LastName');
    disp('Output of rows2vars (first 5 columns)')
    disp(head(T2(:,1:5)))
    disp('Output of rows2varsFS (first 5 columns)')
    T2fs=rows2varsFS(T1,'VariableNamesSource','LastName');
    disp(head(T2fs(:,1:5)))
%}

%{
    %% Example of use of option 'VariableNamingRule'.
    nam={'Temp' 'WindSpeed' 'Rain'};  
    BB=[59.5000    0.1000    0.0500
       63.0000    2.3000    0.0800
       61.7000    3.1000    0.1300
       55.4000    5.7000    0.1500
       62.3000    2.6000    0.8700
       58.8000    6.2000    0.3300];
    T1= array2timetable(BB,'VariableNames',nam,'RowTimes',datetime(2024,4,1:6));
    disp('Original table (first 5 rows)')
    head(T1,5)
    T2=rows2vars(T1,'VariableNamingRule','preserve');
    disp('Output of rows2vars (first 5 columns)')
    disp(head(T2(:,1:5)))
    disp('Output of rows2varsFS (first 5 columns)')
    T2fs=rows2varsFS(T1,'VariableNamingRule','preserve');
    disp(head(T2fs(:,1:5)))
%}

%% Beginning of code
if nargin>1
    Tmp=rows2vars(T1,varargin{:});
else
    Tmp=rows2vars(T1);
end
Tmp.Properties.RowNames=Tmp{:,1};
T2=Tmp(:,2:end);
end
%FScategory:UTIGEN