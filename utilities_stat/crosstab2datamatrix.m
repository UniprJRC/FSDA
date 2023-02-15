function [X,T]=crosstab2datamatrix(N, varargin)
%crosstab2datamatrix recreates the original data matrix X from contingency table N
%
%
%<a href="matlab: docsearchFS('crosstab2datamatrix')">Link to the help function</a>
%
%  Required input arguments:
%
%       N    :    Contingency table (default).
%                 Matrix or Table which contains the input contingency
%                 table of size I-by-J. It contains the frequencies which
%                 have to be inflated. The output data matrix will have size
%                 sum(N(:))-by-2.
%
%  Optional input arguments:
%
%       Lr   :  Vector of row labels. Cell of length I.
%               Cell containing the labels of the rows of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table, because in this case  Lr=N.Properties.RowNames;
%               Example - 'Lr',{'a' 'b' 'c'}
%               Data Types - cell array of strings
%
%       Lc   :  Vector of column labels. Cell of lenght J.
%               Cell containing the labels of the columns of the input
%               contingency matrix N. This option is unnecessary if N is a
%               table, because in this case Lc=N.Properties.VariableNames;
%               Example - 'Lc',{'c1' c2' 'c3' 'c4'}
%               Data Types - cell array of strings
%
%
%  Output:
%
%      X   :  Original data matrix. cell or numeric matrix.
%             Object of class double or cell of size sum(N(:))-by-2
%             containing the original data matrix.
%             Original input which generated the contingency table.
%             Note that input N can be obtained using N=crosstab(X(:,1),X(:,2));
%
%      T   :  Original data matrix in table format.
%             Object of class table of size sum(N(:))-by-2
%             containing the original data matrix.
%             Original input which generated the contingency table.
%             Note that input N can be obtained using N=crosstab(T{:,1},X{:,2});
%             If the labels contained in Lr and Lr could be converted to
%             double the columns of T contain numeric values else their
%             class is categorical
%
% See also crosstab, rcontFS, CressieRead
%
% References:
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('crosstab2datamatrix')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:

%{
    % crosstab2datamatrix with all the default options.
    % In this case the input is a contingency table of class double.
    N=[26 26 23 18 9;
        6  7  9 14 23];
    % No labels for rows and columns are supplied
    X=crosstab2datamatrix(N);
%}

%{
    %% crosstab2datamatrix when input is a contingency table of class table.
    N = [24 23 30; 
         19 43 57;
         13 33 58];
    rownam={'Less_than_5000',  'Between_5000_and_25000' 'Greater_than_25000'};
    colnam= {'Dissatisfied' 'Moderately_satisfied' 'Very_satisfied'};
    if verLessThan('matlab','8.2.0') ==0
        Ntable=array2table(N,'RowNames',matlab.lang.makeValidName(rownam),'VariableNames',matlab.lang.makeValidName(colnam));
        %  Check relationship
        X=crosstab2datamatrix(Ntable);
    else
        X=crosstab2datamatrix(N,'Lr',rownam,'Lc',colnam);
    end
    disp('Compare original contingency table and the one passing through X')
    disp('Original contingency table')
    disp(N)
    disp('Contingency table obtained using crosstab applied to reconstructed original data matrix')
    disp(crosstab(X(:,1),X(:,2)))
%}

%{
    % Example of use of option Lc.
    % In this case just the column names are supplied
    % The default row labels 'r1' 'r2' 'r3' are used
    N = [24 23 30; 
         19 43 57;
         13 33 58];
    colnam= {'Dissatisfied' 'Moderately_satisfied' 'Very_satisfied'};
    X=crosstab2datamatrix(N,'Lc',colnam);
%}

%% Beginning of code

% Check MATLAB version. If it is not smaller than 2013b,  input can be a
% contingency table stored in table format
verMatlab=verLessThan('matlab','8.2.0');

Lr='';
Lc='';
options=struct('Lr',Lr,'Lc',Lc);

UserOptions=varargin(1:2:length(varargin));
if ~isempty(UserOptions)
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:crosstab2datamatrix:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end

    % Write in structure 'options' the options chosen by the user
    if nargin > 2
        for i=1:2:length(varargin)
            options.(varargin{i})=varargin{i+1};
        end
    end
    Lr  = options.Lr;
    Lc  = options.Lc;
end

% Extract labels for rows and columns
if verMatlab ==0 && istable(N)
    Lc=N.Properties.VariableNames;
    Lr=N.Properties.RowNames;
    N=table2array(N);
end

[I,J]=size(N);

n=sum(N(:));

if isempty(Lr) && isempty(Lc)
    Lr=1:I;
    Lc=1:J;
    X=zeros(n,2);
elseif isempty(Lr)
    % default labels for Lr are used
    % default labels for rows of contingency table
    Lr=cellstr(strcat('r',num2str((1:I)')));
    X=cell(n,2);
elseif isempty(Lc)
    % default labels for columns of contingency table
    Lc=cellstr(strcat('c',num2str((1:J)')));
    X=cell(n,2);
else
    X=cell(n,2);
end

% Check that the length of Lr is equal to I (number of rows of the
% contingency table)
if length(Lr)~=I
    error('FSDA:crosstab2datamatrix:WrongInputOpt',['length of cell containing row labels must be equal to ' num2str(I)]);
end
% Check that the length of Lc is equal to J (number of columns of the
% contingency table)
if length(Lc)~=J
    error('FSDA:crosstab2datamatrix:WrongInputOpt',['length of cell containing column labels must be equal to ' num2str(J)]);
end

vartype=["double" "double"];

if isnumeric(X)
    isnumericCol1=true;
    isnumericCol2=true;
    Lrnumeric=Lr;
    Lcnumeric=Lc;
else
    try
        isnumericCol1=true;
        Lrnumeric=cellfun(@str2num,Lr);
    catch
        isnumericCol1=false;
        vartype(1)="categorical";
    end

    try
        isnumericCol2=true;
        Lcnumeric=cellfun(@str2num,Lc);
    catch
        vartype(2)="categorical";
        isnumericCol2=false;
    end

end

T = table('Size',size(X),'VariableTypes',vartype);

% Reconstruct the original data matrix
nij=0;
for i=1:I
    for j=1:J
        %  X(nij+1:nij+N(i,j),:)= [repmat(Lr(i),N(i,j),1) repmat(Lc(j),N(i,j),1)];
        for ij=1:N(i,j)
            if isnumericCol1
                if iscell(X)
                    X{nij+ij,1}= Lrnumeric(i);
                else
                    X(nij+ij,1)= Lrnumeric(i);
                end
                T{nij+ij,1}= Lrnumeric(i);
            else
                if iscell(X)
                    X{nij+ij,1}=Lr{i};
                else
                    X(nij+ij,1)=Lr{i};
                end
                T{nij+ij,1}=string(Lr{i});
            end
            if isnumericCol2
                if iscell(X)
                    X{nij+ij,2}= Lcnumeric(j);
                else
                    X(nij+ij,2)= Lcnumeric(j);
                end
                T{nij+ij,2}= Lcnumeric(j);
            else
                X{nij+ij,2}=Lc{j};
                T{nij+ij,2}=string(Lc{j});
            end
        end
        nij=nij+N(i,j);
    end
end

% outputtable=true;
% if outputtable==true
% end
%
% for jj=1:2
% try
%     na=cellfun(@str2num,X(:,jj));
%     if ~any(isnan(na))
%         for i=1:size(X,1)
%             X{i,jj}=na(i);
%             T{i,jj}=na(i);
%         end
%     end
% catch
%     T{:,jj}=X(:,jj);
% end
% end

% try
%     na=cellfun(@str2num,X(:,2));
%     if ~any(isnan(na))
%         for i=1:size(X,1)
%             X{i,2}=na(i);
%         end
%     end
% catch
% end



end

%FScategory:MULT-Categorical