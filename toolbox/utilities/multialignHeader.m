function [T, WrongCountry, indexEmptyDate]= multialignHeader(Seqs, varargin)
%multialignHeader extracts date and country name from fasta sequence
%
%<a href="matlab: docsearchFS('multialignHeader')">Link to the help function</a>
%
% This function takes as input a vector of structures or a table (Seqs) with fields
% Sequence and Header. Inside field Header it extracts the name of the
% country and the submission date. The name of the country is associated
% with the corresponding continent (unless optional input argument
% addContinent is false). Inside field Sequence it counts the presence of
% particular symbols. The output is a table with height equal to
% height(Seqs). The number of columns in output table is equal to the number of
% fields in Seqs + 5 extra columns for the particular symbols, + 3 extra
% columns referred to Country_name, Continent_name and Submission_date.
%
%
% Required input arguments:
%
%
%        Seqs : Sequences which have to be analyzed. 
%               Vector of structures or table.
%               Vector of structures of length n or table with n rows with
%               the fields 
%               Seqs.Sequence = for the residues and
%               Seqs.Header = or (or 'Seqs.Name') for the labels.
%               Remark: note that Seqs can have other fields. These
%               additional fields will appear as columns in the output
%               table.
%                 Data Types - array of struct or table
%
% Optional input arguments:
%
%  addContinent  : add or not the name of the continent.
%                Boolean. If addContinent is true (default)
%                the categorical vector containing the corresponding
%                continent for each country is added.
%                 Example - 'addContinent','false'
%                 Data Types - boolean
%
% countPartSymb : add or not the 5 extra columns for the particular
%                   symbols. Boolean. 
%                   If countPartSymb is true (default)
%                   five additional columns named
%                   "n?"    "nX"    "n-"    "n*"    "n$"
%                  area added to output table which count for each sequence
%                  the number of  ?, X, -, *, $.
%                 Example - 'countPartSymb','false'
%                 Data Types - boolean
%
% Output:
%
% T : table of height n (same length of input vector of structures or input table).
%       Table containing detailed information (country name, continent and date) 
%       for each sequence. The input sequences can be aligned or not aligned. 
%       1st col = Header (in cell format)
%       2nd col = Sequence (in cell format)
%       3rd col = Date of submission (in datetime format)
%       4th col = Country of the lab (in string format)
%       5th col = Continent associated to the country of the lab (in categorical format).
%       6th-10th = Number of ?, X, -, *, $ in the sequence (numeric format).
%       Remaining columns are referred to the fields present in the input
%       structure Seqs. For example, if the input structure is obtained as
%       the output of function multialign2ref then two additional columns
%       named usedGap and usedDeletion are present (boolean format).
%
%    WrongCountry :  string array containing the strings referred to
%       the country names for which it was not possible to find the
%       reference Country name. A warning is given if this vector is not
%       empty.
%
%  indexEmptyDate    :  numeric vector containing the numbers referred to
%       the rows of T for which it was not possible to extract the date.
%       A warning is given if this vector is not empty. 
%
% See also: multialign, multialign2ref
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('multialignHeader')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example of multialignHeader without optional arguments.
    % Load fasta file containing original covid and other sequences
    Seqs = fastaread("X01sel.txt");
    
    % Call of multialignHeader with all default arguments
    T=multialignHeader(Seqs);
    
    disp('Output table (first 8 rows)')
    disp(head(T))
%}

%{
    %% Example of multialignHeader with 3 output arguments.
    % Load fasta file containing original covid and other sequences
    Seqs = fastaread("X01sel.txt");
    
    % Call of multialignHeader with all default arguments
    [T,indexWrongCountry, indexEmptyDate]=multialignHeader(Seqs);
    if ~isempty(indexWrongCountry)
        disp('Unable to extract some country names')
    end
    if ~isempty(indexEmptyDate)
        disp('Unable to extract some dates')
    end
%}

%{
    %% Example of multialignHeader with optional arguments.
    % Load fasta file containing original covid and other sequences
    Seqs = fastaread("X01sel.txt");
    T=multialignHeader(Seqs,'addContinent',false);
    disp('Output table (first 8 rows)')
    disp(head(T))
%}

%{
    % Example with wrong country name and missing date.
    % Load fasta file containing original covid and other sequences
    Seqs = fastaread("X01sel.txt");
    % Add a record where it is impossible to find both the date and 
    % the country
    Seqs(end+1).Header='Spike	hCoV-19/SSSS/NC-IBV-97020613/2021	wrong date  wrong date';
    Seqs(end).Sequence="XXXXXX-----XXX????";
    T=multialignHeader(Seqs,'addContinent',false);
%}

%{
    % Another example with wrong country names and dates.
    % Find the rows of input sequence for which it was not possible to find the
    % country
    Seqs = fastaread("X01sel.txt");
    % Add a record where it is impossible to find both the date and
    % the country
    Seqs(5).Header='Spike	hCoV-19/ItaIta/NC-IBV-97020613/2021	wrong date  wrong date';
    Seqs(end+1).Header='Spike	hCoV-19/SSSS/NC-IBV-97020613/2021	wrong date  wrong date';
    Seqs(end).Sequence="XXXXXX-----XXX????";
    [T, WrongCountry, indexEmptyDate]=multialignHeader(Seqs,'addContinent',false);
    
    s1n=1:length(Seqs);
    T1=T;
    indexnumber=(1:length(Seqs))';
    T1=addvars(T1,indexnumber,'before',1);
    for i=1:length(WrongCountry)
        boo=T1.Country==WrongCountry(i);
        disp(T1(boo,:))
    end
%}

%{
    % Example with input a table and option countPartSymb.
    Seqs = fastaread("X01sel.txt");
    Seqs=struct2table(Seqs);
    T=multialignHeader(Seqs,'addContinent',false,'countPartSymb',false);
%}

%% Beginning of code

addContinent=true;
countPartSymb=true;

options=struct('addContinent',addContinent,'countPartSymb',countPartSymb);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:multialignHeader:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    aux.chkoptions(options,UserOptions)
end

% Write in structure 'options' the options chosen by the user
if nargin > 2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end
addContinent=options.addContinent;
countPartSymb=options.countPartSymb;

% Check if Seqs is a table or a vector of structures
if istable(Seqs)
    T=Seqs;
    n=size(T,1);
else
    % n = number of sequences to analyze
    n=length(Seqs);
    T=struct2table(Seqs);
end

varNames=string(T.Properties.VariableNames);
T(:,varNames=="Country" | varNames=="Date")=[];
varNames=string(T.Properties.VariableNames);

if addContinent == true
    T(:,varNames=="Continent")=[];
    varNames=T.Properties.VariableNames;
end

if length(varNames)>2
       varNamesMoveToend=varNames(3:end);
else
    varNamesMoveToend='';
end


%% Country and data extraction loop
CountryIni=strings(n,1);
DateStr=strings(n,1);
if countPartSymb == true
    % Find symbols: ?, X, -, *, $ inside sequences
    NumberPartSymb=zeros(n,5);
end


%%
pw = aux.PoolWaitbar(n, 'Progress bar inside multialignHeader');

sequence=T.Sequence;
header=T.Header;
parfor i= 1:n

    % Add the waitbar
    increment(pw);

    if countPartSymb == true
    % strSeq=X{i,2};
    strSeq=sequence{i};

     NumberPartSymb(i,:)=[length(strfind(strSeq,'?')) length(strfind(strSeq,'X'))  length(strfind(strSeq,'-')) ...
         length(strfind(strSeq,'*'))    length(strfind(strSeq,'$')) ];
    end

    %% Extract the country and the date from header

    strDate=header{i};
    [matchSlash]=regexp(strDate,'/');

    strtest=string(strDate(matchSlash(1)+1:matchSlash(2)-1));
    if strtest=="env" || strtest==""
        strtest=string(strDate(matchSlash(2)+1:matchSlash(3)-1));
    end
    % Store the name of the country deleting extra spaces
    CountryIni(i)=strip(strtest);

    % Extract the date
    strspli=strsplit(strDate);
    tent=1;

    while tent<=min(length(strspli),7)

        strPdate=strspli(tent);

        [posPipe]=strfind(strPdate{:},"|"+digitsPattern(4)+"-"+digitsPattern(1,2) );
        if ~isempty(posPipe)
            strPdatechar=strPdate{:};
            strPdatechar=strPdatechar(posPipe(1)+1:posPipe(1)+10);

            posLastPipe=strfind(strPdatechar,'|');
            if ~isempty(posLastPipe)
                strPdatechar(posLastPipe:end)=[];
            end
            posHyphen=strfind(strPdatechar,'-');
            if isscalar(posHyphen)
                strPdatechar=[strPdatechar '-01'];
            end

            strPdate={strPdatechar};
        end

        strPdate=strrep(strPdate,'-00','-01');

        try
            DateStr(i)=datetime(strPdate,'InputFormat','uuuMM-dd');
            tent=10;
        catch
            try
                DateStr(i)=datetime(strPdate,'InputFormat','dd/MM/yy');
                tent=10;
            catch

                try
                    DateStr(i)=datetime(strPdate,'InputFormat','uuuu-MM-dd');
                    tent=10;
                catch
                    tent=tent+1;
                end

            end
        end

    end

end

%% For the country names named Country
% In this case we put as country name what is contained after the final pipe sign
seq=1:n;
badname=seq(CountryIni=="Country");
if ~isempty(badname)
    for j=1:length(badname)
        strDate=header{badname(j)};
        % strDate=strDate{:};
        [matchPipe]=strfind(strDate,'|');
        if ~isempty(matchPipe)
            CountryIni(badname(j))=strDate(matchPipe(end)+1:end);
        end
    end
end

%& Load reference country names  using file .xlsx
% which is inside subfolder multivariate_regression
FileName='continents-according-to-our-world-in-data.xlsx';

try
    Nam=readtable(FileName,"Sheet","WrongNames","TextType","string");
catch
    error('FSDA:multialignHeader:MissingFile','File "continents-according-to-our-world-in-data.xlsx" not found.');
end
NomiBadlyWritten=Nam.Wrong_Country_Name;
NomiCorrectlyWritten=Nam.Correct_Country_Name;
Country=CountryIni;

for i=1:length(NomiBadlyWritten)
    boo=NomiBadlyWritten(i)==CountryIni;
    if any(boo)
        Country(boo)=NomiCorrectlyWritten(i);
    end
end
T=addvars(T,Country,'After','Sequence');

Date=datetime(DateStr);
T=addvars(T,Date,'After','Sequence');

if  countPartSymb == true
    d1=array2table(NumberPartSymb,'VariableNames',["n?" "nX" "n-" "n*" "n$"]);

    T=[T d1];
end

%% Check whether there are no dates with missing values
indexEmptyDate=seq(isnat(T.Date));

if ~isempty(indexEmptyDate)
    disp('Records for which date was impossible to find date')
    Xempty=T(indexEmptyDate,:);
    Xempty.Properties.RowNames=string(indexEmptyDate);
    disp(Xempty)
    warning('FSDA:multialignHeader:MissDate','Missing dates in output table');
end

%% Check whether extracted countries have admissible names
% With the phrase "admissibile names" we mean all the country names that
% are contained in sheet
Wld=readtable(FileName,"Sheet","CountryContinent","TextType","string");

% Wld.Entity contains the list of all the country names
NomiAll=[Wld.Entity; "Animal"];
% wrongNameForCountry is the vector which contains all country names that
% have been found but are not present inside NomiAll
[WrongCountry,indCountry]=setdiff(T.Country,NomiAll);


if ~isempty(WrongCountry)
    disp('Records with wrong country name')
    XwrongNameForCountry=T(indCountry,:);
    XwrongNameForCountry.Properties.RowNames=string(WrongCountry);
    disp(XwrongNameForCountry)
    warning('FSDA:multialignHeader:WngCountry','Wrong country names in output table');
end

%% Add the continent
if addContinent == true
    Continent=strings(n,1);
    valueset=["Asia" "Europe" "Africa" "Oceania" "North America" "Antarctica" "South America"];
    Continent=categorical(Continent,valueset);

    for i=1:length(Wld.Entity)
        boo=Country==Wld.Entity(i);
        Continent(boo)=Wld.Continent(i);
    end
   T=addvars(T,Continent,'After','Country');
end

if ~isempty(varNamesMoveToend)
    T=movevars(T,varNamesMoveToend,'After',width(T));
end
%FScategory:UTIGEN

