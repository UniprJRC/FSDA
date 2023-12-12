function ysmo=repDupValWithMean(x,y,varargin)
%repDupValWithMean replaces values of y which have non unique elements in vector x with local means
%
%<a href="matlab: docsearchFS('repDupValWithMean')">Link to the help function</a>
%
%  Required input arguments:
%
%    x:         vector with values to analyze. Vector.
%               A vector with n elements that may contain duplicated
%               values.
%               It can be either a row or a column vector.
%    y :        Vector on which the computations have to be made.
%               Vector.
%               It can be either a row or a column vector.
%
%  Optional input arguments:
%
%   consec  :  how to compute the local means. Boolean. If consec is true
%               the unique entries in vector x are defined as those values
%               which are equal and consecutive.
%               The default value of consec is false, therefore the unique
%               entries in vector x are defined as those values which are
%               equal and consecutive. When x is already sorted in order to
%               speed up calculations it is efficient to call the procedure
%               with the third argument consec set to true, because in this
%               case it avoids calling the MATLAB routine issorted to check
%               whether input vector x is sorted.
%               Example - false
%               Data Types - Boolean
%
%       w  : weights for the observations. Vector. Row or column vector of
%           length n containing the weights associated to each
%           observations. If w is not specified we assum $w=1$ for $i=1,
%           2, \ldots, n$.
%           Example - 'w',1:n
%           Data Types - double
%
%
% Output:
%
%   ysmo    :  Smoothed vector y with local means for non unique x values. Vector.
%               ysmo is a vector with the same dimension of y containing in
%               correspondence of the values of y which have non unique
%               entries in vector x, the arithmetic means of y for the
%               corresponding elements. The unique entries in vector x can
%               be defined as the values which are equal and consecutive
%               or simply equal but non necessarily consecutive (depending
%               on optional input argument consec).
%
% More About:
%
% This function does not use loops and is based just on built in MATLAB
% functions: diff, cumsum and accumarray. See also function accumulator
% from John D'Errico in the file exchange.
%
%
%
% See also: diff, accumarray
%
%
% References:
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('repDupValWithMean')">Link to the help function</a>
%
%$LastChangedDate:: 2019-02-06 18:12:59 #$: Date of the last commit


% Examples:

%{
    %% Case 1: x is already ordered.
    % Note that in this case x is ordered therefore the average between
    % consecutive values which are equal or the average of equal values is
    % the same.
    x=[ones(5,1); 6; 7; 8.2; 8.2; 10];
    % y is a vector containing any real number.
    y=(1:10)';
    ysmo=repDupValWithMean(x,y);
    disp(['      x   '   '      y   ' '      ysmo  '])
    disp([x y ysmo])
    % The first 5 elements of ysmo are equal to mean(y(1:5)) because the
    % corresponding elements of x share the same value.
    % The elements in position 8 and 9 of ysmo are equal to mean(y([8:9]))
    % because the corresponding elements of x share the same value.
    % All the other elements of vector ysmo are equal to y.
%}

%{
    %% Case 2: x is non ordered.
    % Note that in this case x is not ordered therefore if the function is
    % called with just two arguments it takes the average of the elements of
    % y which correspond to elements of x which are equal
    x=[8.2; 1.0; 1.0; 6.0; 7.0; 10.0; 1.0; 8.2; 1.0; 1.0];
    % y is a vector containing any real number.
    y=(11:20)';
    ysmo=repDupValWithMean(x,y);
    disp(['      x   '   '      y   ' '      ysmo  '])
    disp([x y ysmo])
    % Elements 1 and 8  share the same value of x therefore
    % ysmo(1)=smo(8) = mean(y([1 8])).
    % Elements 2, 3, 7, 9 and 10 share the same value of x therefore
    % ysmo(2)=ysmo(3)=ysmo(7)=ysmo(9)=ysmo(10)= mean(y([2 3 7 9 10])).
    % All the other elements of vector ysmo are equal to y.
%}

%{
    %% Case 3: x is non ordered and optional argument consec is true.
    % Note that in this case x is not ordered therefore if the function is
    % called with the third argument consec equal to true, this function
    % computes the average of the elements of y which correspond to
    % elements of x which are equal and consecutive.
    x=[8.2; 1.0; 1.0; 6.0; 7.0; 10.0; 1.0; 1.0; 1.0; 8.2];
    % y is a vector containing any real number.
    y=(11:20)';
    ysmo=repDupValWithMean(x,y,'consec',true);
    disp(['      x   '   '      y   ' '      ysmo  '])
    disp([x y ysmo])
    % Elements 2 and 3 of x are equal and consecutive therefore
    % ysmo(2)=ysmo(3) = mean(y([2 3])).
    % Elements 7, 8 and 9 of x are equal and consecutive therefore
    % ysmo(7)=ysmo(8)=ysmo(9) = mean(y([7 8 9])).
    % All the other elements of vector ysmo are equal to y.
%}

%{
    %% Simulation study to compare repDupValWithMean with and without loops.
    % Create function 'repDupValWithMeanLoop.m' and write it into a file.
    % At the end this file will be deleted.
    name='repDupValWithMeanLoop.m';
    filetmpID=fopen([pwd filesep name],'w');
    % % The implementation of this function using loops is given below.
    %     function  ysmoC=repDupValWithMeanLoop(x,y)
    %     n=length(x);
    %     ysmoC=y;
    %     j0=1;
    %     salta=0;
    %     for j=1:n-1
    %         if salta==0
    %             sm=ysmoC(j);
    %         end
    %         if x(j+1) <=x(j)
    %             sm=sm+ysmoC(j+1);
    %             salta=1;
    %             if j==n-1
    %                 sm=sm/(j-j0+2);
    %                 ysmoC(j0:j+1)=sm;
    %             end
    %         else
    %             salta=0;
    %             sm=sm/(j-j0+1);
    %             ysmoC(j0:j)=sm;
    %             j0=j+1;
    %         end
    %     end
    %     end
    %
    % In principle one should take the above function and save it in a file
    % In order to speed up thing we put it inside variable outstring and write
    % the content of outstring into a file.
    outstring=sprintf(['function  ysmoC=repDupValWithMeanLoop(x,y) \r' ...
        'n=length(x); \r' ...
        ' ysmoC=y; \r' ...
        'j0=1;     \r' ...
        'salta=0;      \r' ...
        'for j=1:n-1     \r' ...
        '    if salta==0  \r' ...
        '        sm=ysmoC(j); \r' ...
        '    end               \r ' ...
        '    if x(j+1) <=x(j)   \r' ...
        '       sm=sm+ysmoC(j+1);   \r' ...
        '        salta=1;          \r ' ...
        '        if j==n-1          \r' ...
        '            sm=sm/(j-j0+2); \r' ...
        '            ysmoC(j0:j+1)=sm; \r' ...
        '        end                    \r' ...
        '    else                       \r' ...
        '        salta=0;               \r' ...
        '        sm=sm/(j-j0+1);       \r ' ...
        '        ysmoC(j0:j)=sm;       \r ' ...
        '        j0=j+1;               \r ' ...
        '    end                        \r' ...
        'end                            \r' ...
        'end']);
    fprintf(filetmpID,'%s',outstring);
    fclose(filetmpID);

    % Simulation study to compare the two implementations.
    nsimul=10000;
    n=10000;
    imax=20;
    totimeOpt2=0;
    totimeOpt3=0;
    for j=1:nsimul
        x=randi(imax,n,1);
        y=randn(n,1);
        x=sort(x);
        tic
        ysmo2=repDupValWithMean(x,y);
        totimeOpt2=toc+totimeOpt2;
        tic
        ysmo3=repDupValWithMeanLoop(x,y);
        totimeOpt3=toc+totimeOpt3;

        % Check that the two implementations produce the same results.
        if max(abs(ysmo2-ysmo3))>1e-9
            error('The two implementations do not produCe equal results')
        end
    end
    disp('Comparison of times based on 10000 replicates')
    disp('Implementation without loops')
    disp(totimeOpt2)
    disp('Implementation using loops')
    disp(totimeOpt3)
    % Remove temporary file repDupValWithMeanLoop.m
    delete repDupValWithMeanLoop.m
%}


%{
   % Case 4: x is already ordered and weights are supplied.
    % Note that in this case x is ordered therefore the average between
    % consecutive values which are equal or the average of equal values is
    % the same.
    x=[ones(5,1); 6; 7; 8.2; 8.2; 10];
    % y is a vector containing any real number.
    y=(1:10)';
    w=(21:30)';
    ysmo=repDupValWithMean(x,y,'w',w);
    disp(['      x   '   '      y   ' '     w' '      ysmo  '])
    disp([x y w ysmo])
    % The first 5 elements of ysmo are equal to sum(y(1:5).*w(1:5))/sum(w(1:5))
    % because the corresponding elements of x share the same value.
    % The elements in position 8 and 9 of ysmo are equal to
    % sum(y(8:9).*w(8:9))/sum(w(8:9)) because the corresponding elements of x
    % share the same value.
    % All the other elements of vector ysmo are equal to y.
%}

%% Beginning of code

w=[];
consec=false;

if nargin > 2
options=struct('w',w,'consec',consec);

[varargin{:}] = convertStringsToChars(varargin{:});
UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSR:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        chkoptions(options,UserOptions)
    end
    
    % Write in structure 'options' the options chosen by the user
    
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end


consec=options.consec;
w=options.w(:);
end

if consec == true
    isSortedx = true;
else
    isSortedx = issorted(x);
end

n=length(x);
if isempty(w)
    w=ones(n,1);
end

% Determine if y is a row vector.
rowvec = isrow(y);

% Convert to column
x=x(:);
y=y(:);

% Using the original vector x or the sorted version
if isSortedx == true
    Sortx=x;
else
    [Sortx, indSortx] = sort(x);
    y=y(indSortx);
    w=w(indSortx);
end

% Implementation without loops using diff, cumsum and accumarray
dSortx = diff(Sortx);
groupsSortA = dSortx ~= 0;
groupsSortA = [true; groupsSortA];
idx = cumsum(groupsSortA);  % Lists position, starting at 1.
% The above lines can be replaced by the slower instruction
% [~,~,idx]=unique(x);
% Compute the (weighted) mean of the values of vector y which have duplicate values in
% vector x
% b=accumarray(idx,y,[],@mean);
b=accumarray(idx,1:n,[],@(x)sum(y(x) .* w(x)./ sum(w(x))));

ysmo=b(idx);

if isSortedx==false
    invIndSortX = indSortx;
    invIndSortX(invIndSortX) = 1:length(x);  % Find inverse permutation.
    ysmo=ysmo(invIndSortX);
end

% If y is row vector, return ysmo as row vector.
if rowvec
    ysmo = ysmo';
end
end
%FScategory:UTIGEN