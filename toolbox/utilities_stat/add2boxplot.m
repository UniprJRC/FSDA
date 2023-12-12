function add2boxplot(Y,g)
%add2boxplot adds labels to the boxplot figure
%
%<a href="matlab: docsearchFS('add2boxplot')">Link to the help function</a>
%
% add2boxplot assumes that one or more boxplots have been created
% using function boxplot. This function adds to the plot the labels
% associated with the outliers. Note that Y and g must have exactly the
% same dimensions to those called by boxplot. Note also that if Y is the
% associated table, the labels are referred to the rownames of Y.
%
%
%  Required input arguments:
%
% Y :           Input data. Vector or matrix or table.
%               This is the first input argument which has been passed to
%               function boxplot. Note that while boxplot requires an
%               array, Y can also be a table with the same rows and columns
%               to that passed to boxplot.
%                Data Types - single|double
%
% Optional input arguments:
%
%        g :     Grouping variable. This is exactly equal to g which has
%                been supplied to function boxplot.
%                 Example - [ones(10,1);ones(20,1)]
%                 Data Types - numeric vector | character array | string array | cell array | categorical array
% Output:
%
%
%
%
% See also: boxplot
%
% References:
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('add2boxplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

%{
    % Example 1: Y vector and no groups.
    % boxplot is called with input variable Y as a vector and there is no
    % grouping variable.
    rng(2)
    y=randn(1300,1);
    boxplot(y)
    add2boxplot(y);
    % Note that the input argument of add2boxplot could be a table
    % ytable=array2table(y);
    % add2boxplot(ytable);
%}

%{
    % Example 2: Y matrix and no groups.
    % boxplot is called with input variable Y as a matrix and there is no
    % grouping variable.
    rng(2)
    y=randn(300,5);
    boxplot(y)
    add2boxplot(y);
%}

%{
    % Example 3: Y vector and there are groups.
    % boxplot is called with input variable Y as a vector and there is a
    % grouping variable.
    load fisheriris
    boxplot(meas(:,1),species)
    add2boxplot(meas(:,1),species)
%}


%{
    %% Example 4: Y matrix and there are groups.
    % boxplot is called with input variable Y as a matrix and there is a
    % grouping variable
    Y=trnd(5,100,4);
    % define the grouping variable.
    g={'zz';'bb';'cc'; 'dd'};
    boxplot(Y,g)
    add2boxplot(Y,g)
%}

%{
    %% An example where Y is table.
    % load dataset referred to Italian cities about quality of life.
    load citiesItaly.mat
    Yst=citiesItaly;
    Yst{:,:}=zscore(citiesItaly{:,:});
    boxplot(Yst{:,:},'Labels',citiesItaly.Properties.VariableNames,'Jitter',0);
    add2boxplot(Yst)
%}

%{
    %% An example where Y is table and there are subplots.
    load citiesItaly.mat
    close all
    citiesItaly{:,:}=zscore(citiesItaly{:,:});
    
    subplot(1,2,1)
    Yst=citiesItaly{1:46,:};
    Ystt=citiesItaly(1:46,:);
    boxplot(Yst,'Labels',citiesItaly.Properties.VariableNames,'Jitter',0);
    add2boxplot(Ystt)
    title("Towns from north of Italy")

    subplot(1,2,2)
    Yst1=citiesItaly{47:end,:};
    Yst1t=citiesItaly(47:end,:);
    boxplot(Yst1,'Labels',citiesItaly.Properties.VariableNames,'Jitter',0);
    add2boxplot(Yst1t)
    title("Towns from center and south of Italy")

%}
%% Beginning of code
if nargin>1
    unig=unique(g,'stable');
else
    unig='';
end

% Find in current figure the objects which are tagged as outliers
objOutliers = findobj(gcf,'Tag','Outliers');
% Extract the x and y coordinates of the outliers
outliersYcoo = get(objOutliers,'YData');
outliersXcoo = get(objOutliers,'XData');

[n,p]=size(Y);
lunig=length(unig);
% Just in case the user has produced many subplots it is necessary
% to use just the information about the last length(unig) groups
% which have been created
if p>1
    outliersYcoo=outliersYcoo(1:p,:);
    outliersXcoo=outliersXcoo(1:p,:);
else
    if lunig>0
        outliersYcoo=outliersYcoo(1:lunig,:);
        outliersXcoo=outliersXcoo(1:lunig,:);
    end
end

if ~isempty(outliersYcoo)


    seq=(1:n)';
    % Define the row labels to add to the plot
    % If Y is a table we add the corresponding rownames else we just add the row
    if istable(Y)
        Yd=Y{:,:};
        rownam=Y.Properties.RowNames;
        if isempty(rownam)
            rownam=string(seq);
        end
    else
        Yd=Y;
        rownam=string(seq);
    end

    if isempty(unig)
        % This is the case when second argument is missing (there is no
        % classification variable)
        for j=1:p
            if iscell(outliersYcoo)
                outjYcoo=outliersYcoo{end+1-j};
                outjXcoo=outliersXcoo{end+1-j};
                outjXcoo=outjXcoo(:);
            else
                outjYcoo=outliersYcoo;
                outjXcoo=outliersXcoo;

            end
            minc=min(abs(Yd(:,j)-outjYcoo),[],2);
            rown=seq(minc<1e-10);
            if ~isempty(rown)
                [Ydjsorted,indsorj]=sort(Yd(rown,j));
                rownamj=rownam(rown);
                rownamjsor=rownamj(indsorj);

                % Alternate the labels left and right
                lrown=length(rown);
                even=2:2:lrown;
                odd=1:2:lrown;
                text(outjXcoo(odd)+0.1,Ydjsorted(odd),rownamjsor(odd),'HorizontalAlignment','left')
                text(outjXcoo(even)-0.1,Ydjsorted(even),rownamjsor(even),'HorizontalAlignment','right')
            end

        end
    else

        % This is the case when second argument is present (there is a grouping variable)
        % and Y has just one column.

        if size(Yd,2)==1


            for j=1:length(unig)
                booj=strcmp(g,unig(j));

                Ydj=Yd(booj);
                rownamj=rownam(booj);
                outjYcoo=outliersYcoo{end+1-j};
                outjXcoo=outliersXcoo{end+1-j};
                outjXcoo=outjXcoo(:);

                minc=min(abs(Ydj-outjYcoo),[],2);
                rown=seq(minc<1e-10);
                if ~isempty(rown)
                    [Ydjsorted,indsorj]=sort(Ydj(rown));
                    rownamjsel=rownamj(rown);
                    rownamjselsor=rownamjsel(indsorj);
                    % Alternate the labels left and right
                    lrown=length(rown);
                    even=2:2:lrown;
                    odd=1:2:lrown;
                    text(outjXcoo(odd)+0.1,Ydjsorted(odd),rownamjselsor(odd),'HorizontalAlignment','left')
                    text(outjXcoo(even)-0.1,Ydjsorted(even),rownamjselsor(even),'HorizontalAlignment','right')
                end
            end
        else
            % There is a grouping variable and Y has more than one column
            for j=1:length(unig)
                Ydj=Yd(:,j);
                outjYcoo=outliersYcoo{end+1-j};
                outjXcoo=outliersXcoo{end+1-j};
                outjXcoo=outjXcoo(:);

                minc=min(abs(Ydj-outjYcoo),[],2);
                rown=seq(minc<1e-10);

                [Ydjsorted,indsorj]=sort(Ydj(rown));
                rownamj=rownam(rown);
                rownamjsor=rownamj(indsorj);

                if ~isempty(rown)
                    % Alternate the labels left and right
                    lrown=length(rown);
                    even=2:2:lrown;
                    odd=1:2:lrown;
                    text(outjXcoo(odd)+0.1,Ydjsorted(odd),rownamjsor(odd),'HorizontalAlignment','left')
                    text(outjXcoo(even)-0.1,Ydjsorted(even),rownamjsor(even),'HorizontalAlignment','right')
                end
            end
        end
    end  % close if linked to the presence of classification variable

end
end

%FScategory:VIS-Mult